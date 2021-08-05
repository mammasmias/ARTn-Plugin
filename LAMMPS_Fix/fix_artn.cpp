/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_artn.h"

#include "atom.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "region.h"
#include "respa.h"
#include "update.h"
#include "variable.h"

#include "min_fire.h"

#include <cstring>

#include <fstream>
#include <iomanip>
#include <iostream>

using namespace std;

using namespace LAMMPS_NS;
using namespace FixConst;

extern "C"{
  void artn_( double *const f, double* etot, const int nat, const int *ityp, const char *elt, double *const tau, const int *order, const double *lat, const int *if_pos, int* disp, bool* lconv );
  void move_mode_( const int nat, double *const f, double *const vel, double* etot, int* nsteppos, double* dt_curr, double* alpha, const double* alpha_init, const double* dt_init, int* disp );
}

#ifdef INTEL 
extern int __artn_params_mp_iperp;
extern int __artn_params_mp_perp;
extern int __artn_params_mp_relx;
#else
extern int __artn_params_MOD_iperp;
extern int __artn_params_MOD_perp;
extern int __artn_params_MOD_relx;
#endif

/* ---------------------------------------------------------------------- */

FixARTn::FixARTn( LAMMPS *lmp, int narg, char **arg ): Fix( lmp, narg, arg ) 
{

  if (narg < 3) error->all(FLERR,"Illegal fix ARTn command");


  // ...Set mpi parameters
  MPI_Comm_rank( world, &me );
  MPI_Comm_size( world, &nproc );

  nloc = nullptr;
  ftot = nullptr;
  xtot = nullptr;
  vtot = nullptr;
  order_tot = nullptr;

  // ...Set to null/0 all the variable of the class
  istep = 0;
  nword = 0;
  word = nullptr;

  nextblank = 0;
  f_prev = nullptr;
  order = nullptr;
  elt = nullptr;
  if_pos = nullptr;

  alpha_init = 0.0;
  alphashrink = 0.0;
  dt_init = 0.0;
  dtsk = 0.0; 
  dtgrow = 0.0;
  tmax = 0.0;
  tmin = 0.0;
  dtmax = 0.0;
  dtmin = 0.0;

  fire_integrator = 0.0; 
  ntimestep_start = 0.0;

  delaystep_start_flag = 1;

  pe_compute = nullptr;

  etol = 0.0;
  ftol = 0.0;



  // ...Save the Fire Parameter - Init:
  alpha_init = 0.1;
  alphashrink = 0.99;

  //dt_init = 0.001; //update->dt;
  dtsk = 0.5;
  dtgrow = 1.1;

  dmax = 0.5;

  tmax = 20.;
  tmin = 0.02;

  fire_integrator = 0;

  //if( narg == 3 )return;

  int iarg(3);
  while( iarg < narg ){
    
    if( strcmp(arg[iarg], "alpha0") == 0 ){
      if (iarg+2 > narg) error->all(FLERR,"Illegal min_modify command");
      alpha_init = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"alphashrink") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal min_modify command");
      alphashrink = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"dtgrow") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal min_modify command");
      dtgrow = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"dtshrink") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal min_modify command");
      dtsk = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"tmax") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal min_modify command");
      tmax = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"tmin") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal min_modify command");
      tmin = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    //} else if (strcmp(arg[iarg],"halfstepback") == 0) {
    //  if (iarg+2 > narg) error->all(FLERR,"Illegal min_modify command");
    //  if (strcmp(arg[iarg+1],"yes") == 0) halfstepback_flag = 1;
    //  else if (strcmp(arg[iarg+1],"no") == 0) halfstepback_flag = 0;
    //  else error->all(FLERR,"Illegal min_modify command");
    //  iarg += 2;
    } else if (strcmp(arg[iarg],"initialdelay") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal min_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0) delaystep_start_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) delaystep_start_flag = 0;
      else error->all(FLERR,"Illegal min_modify command");
      iarg += 2;
    //} else if (strcmp(arg[iarg],"vdfmax") == 0) {
    //  if (iarg+2 > narg) error->all(FLERR,"Illegal min_modify command");
    //  max_vdotf_negatif = utils::numeric(FLERR,arg[iarg+1],false,lmp);
    //  iarg += 2;
    } else if (strcmp(arg[iarg],"integrator") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal min_modify command");
      if (strcmp(arg[iarg+1],"eulerimplicit") == 0) fire_integrator = 0;
      else error->all(FLERR,"Illegal min_modify command");
      iarg += 2;
    //} else if (strcmp(arg[iarg],"norm") == 0) {
    //  if (iarg+2 > narg) error->all(FLERR,"Illegal min_modify command");
    //  if (strcmp(arg[iarg+1],"two") == 0) normstyle = TWO;
    //  else if (strcmp(arg[iarg+1],"max") == 0) normstyle = MAX;
    //  else if (strcmp(arg[iarg+1],"inf") == 0) normstyle = INF;
    //  else error->all(FLERR,"Illegal min_modify command");
    //  iarg += 2;
    } else {
      int n = modify_param(narg-iarg,&arg[iarg]);
      if (n == 0) error->all(FLERR,"Illegal fix_modify command");
      iarg += n;
    }

  }
  

}


/* ---------------------------------------------------------------------- */

FixARTn::~FixARTn() {

  /* deallocate the array */
  memory->destroy( word );
  memory->destroy( order );
  memory->destroy( elt );
  memory->destroy( f_prev );
  memory->destroy( if_pos );

  memory->destroy( nloc );
  memory->destroy( ftot );
  memory->destroy( xtot );
  memory->destroy( vtot );
  memory->destroy( order_tot );
  
  pe_compute = nullptr;

}


/* ---------------------------------------------------------------------- */

int FixARTn::setmask()
{
  int mask = 0;
  //mask |= POST_FORCE;
  //mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixARTn::init() {


  // Check if FIRE minimization is well define

  if( strcmp(update->minimize_style, "fire") != 0 ) 
    error->all(FLERR,"Fix/ARTn must be used with the FIRE minimization"); 


  // compute for potential energy

  int id = modify->find_compute("thermo_pe");
  if (id < 0) error->all(FLERR,"Fix/ARTn could not find thermo_pe compute");
  pe_compute = modify->compute[id];


  // Initialize some variable

  istep = 0;
  nword = 6;

}





/* ---------------------------------------------------------------------- */

void FixARTn::setup( int /*vflag*/ ) {

  // Depending of the module environment we can change the move interface
  // and other thing

  //cout<<" * IN SETUP..." << endl;


}





/* ---------------------------------------------------------------------- */

void FixARTn::min_setup( int vflag ) {

  // Call here if it is needed - To confirm
  //cout<<" * IN MIN_SETUP..." << endl;
  //post_force( vflag );

  class Min *minimize = update-> minimize;

  cout<< " * FIX/ARTn::CHANGE PARAM..."<<endl;

  /*
  -- Change & Save the initial Fire Parameter 
     Exept: delaystep_start_flag = 1 (ALWAYS)
  */

  nword = 12;
  if( word )memory->destroy( word );
  memory->create( word, nword, 20, "fix:word" );

  strcpy( word[0], "alpha0" );
  strcpy( word[1], "0.0");   // at the begining
  strcpy( word[2], "alphashrink") ;
  strcpy( word[3], "0.99") ;
  strcpy( word[4], "delaystep") ;
  strcpy( word[5], "20") ;
  strcpy( word[6], "halfstepback") ;
  strcpy( word[7], "no") ;
  strcpy( word[8], "integrator") ;
  strcpy( word[9], "eulerimplicit") ;
  strcpy( word[10], "dmax") ;
  strcpy( word[11], "0.5") ;

  minimize-> modify_params( nword, word );


  dt_init = update->dt;

  dtmax = tmax*dt_init;
  dtmin = tmin*dt_init;

  //fire_integrator = 0;
  ntimestep_start = update->ntimestep;
  
  etol = update->etol;
  ftol = update->ftol;

  // ...We Change the convergence criterium to control it
  update->etol = 1e-18;
  update->ftol = 1e-18;

  // ...Print the new Parameters
  minimize-> setup_style();

  // ...Print the Initial Fire Parameters
  cout<< " * Alpha0->"<< alpha_init<< endl;
  cout<< " * dt0->"<< dt_init<< endl;
  cout<< " * dtmin->"<< dtmin<< endl;
  cout<< " * dtmax->"<< dtmax<< endl;


  // ...Copy the order of atom
  int nat = atom-> natoms;
  natoms = nat; // Save for the rescale step
  int nlocal = atom->nlocal;
  tagint *itag = atom->tag;
  memory->create( order, nat, "Fix/artn::order" );
  for( int i(0); i < nat; i++ ) order[i] = itag[i];


  // ...Allocate the previous force table :: IS LOCAL
  memory->create( f_prev, nlocal, 3, "fix/artn:f_prev");
  nextblank = 0;

  // Allocate memory for global force array
  memory->create( ftot, natoms, 3, "fix/artn:ftot");


  // ...Define the Element array for each type
  memory->create( elt, nat, "fix/artn:");
  const int *ityp = atom->type;
  for( int i(0); i < nat; i++ ) elt[i] = alphab[ ityp[i] ];


  // ...Define the constrains on the atomic movement
  memory->create(if_pos,nat,3,"fix/artn:if_pos");
  //memset( if_pos, 1, 3*nat ); 
  for( int i(0); i < nat; i++ ){
     if_pos[ i ][0] = 1;
     if_pos[ i ][1] = 1;
     if_pos[ i ][2] = 1;
  }


  // ...Parallelization
  memory->create( nloc, nproc, "fix/artn:nloc" );
  memory->create( ftot, nat, 3, "fix/artn:ftot" );
  memory->create( xtot, nat, 3, "fix/artn:xtot" );
  memory->create( vtot, nat, 3, "fix/artn:xtot" );
  memory->create( order_tot, nat, "fix/artn:order_tot" );


}










/* ---------------------------------------------------------------------- */

void FixARTn::min_post_force( int vflag ) {

  //cout<<" * IN MIN_POST_FORCES..." << endl;
  post_force( vflag );
  //cout<<" * OUT MIN_POST_FORCES..." << endl;

  //if( istep == 9 )exit(0);

}












/* ---------------------------------------------------------------------- */

void FixARTn::post_force( int /*vflag*/ ){

  // call pARTn library...


  // ...Link the minimizer
  class Min *minimize = update-> minimize; 
  if( !minimize )error->all(FLERR,"fix/ARTn::Min_vector is not linked");


  // ...We Change the convergence criterium to control it
  update->etol = 1e-18;
  update->ftol = 1e-18;
  //cout<< " * FIX_ARTn::MINIMIZE CONV:"<< update-> etol<< " | "<< update->ftol<< endl;


  // ...Basic Array to work
  int nlocal = atom->nlocal;
  int nat = atom->natoms;
  double **tau = atom->x;
  double **f = atom->f;
  double **vel = atom->v ;
  const int *ityp = atom->type;




  // ...Comput V.F to know in which part of alogrithm energy_force() is called
  double vdotf = 0.0, vdotfall;
  for( int i = 0; i < nlocal; i++ )
    vdotf += vel[i][0]*f[i][0] + vel[i][1]*f[i][1] + vel[i][2]*f[i][2];
  MPI_Allreduce( &vdotf, &vdotfall, 1, MPI_DOUBLE, MPI_SUM, world );
  //cout<<" * FIX_ARTN::Computed v.f = "<< vdotfall<<" | step "<< istep <<endl;






  // ---------------------------------------------------------------------  
  // ...If v.f is 0 or under min_fire call the force to ajust 
  // the integrator step parameter dtv
  if( !(vdotfall > 0) && update-> ntimestep > 1 && nextblank ){


    // ...Rescale the force if the dt has been change 
    double rscl = dt_curr / update->dt;

    // ...Reload the previous ARTn-force 
    for( int i(0); i < nlocal; i++){
      f[i][0] = f_prev[i][0] * rscl*rscl ;
      f[i][1] = f_prev[i][1] * rscl*rscl ;
      f[i][2] = f_prev[i][2] * rscl*rscl ;
    }

    //cout<< " ********* RETURN WITHOUT COMPUTE ARTn | v.f = "<< vdotfall<< " | Rescale "<< rescale<<" "<< dt_curr<< " " <<update->dt<<endl;


    // ...Next call should be after the integration
    nextblank = 0;
    return;

  } // ---------------------------------------------------------------------------------------- RETURN











  /*****************************************
   *  Now we enter in the ARTn Algorithm
   *****************************************/



  // ...Extract the energy in Ry for ARTn  
  double etot = pe_compute->compute_scalar();


  /* ...Convergence and displacement */
  bool lconv;
  int disp;



  // ...Build it with : domain-> boxlo[3], boxhi[3] and xy, xz, yz  
  double lat[3][3];
  double dx = domain->boxhi[0] - domain->boxlo[0], 
         dy = domain->boxhi[1] - domain->boxlo[1], 
         dz = domain->boxhi[2] - domain->boxlo[2];

  lat[0][0] = dx;     lat[0][1] = domain->xy    ; lat[0][2] = domain->xz ;
  lat[1][0] = 0.0 ;   lat[1][1] = dy            ; lat[1][2] = domain->yz ;
  lat[2][0] = 0.0 ;   lat[2][1] = 0.0           ; lat[2][2] = dz ;



  
  // PRE ARTn
  double *rmass = atom->rmass;
  double *mass = atom->mass;



  // ...Verification of the local size
  MPI_Allgather( &nlocal, 1, MPI_INT, nloc, nproc, MPI_INT, world );
  int ntot(0), lresize(0);
  for( int ipc(0); ipc < nproc; ipc++ ){
    ntot += nloc[ipc];
    //printf("[%d] nlocal %d \n", ipc, nloc[ipc] ); 
  }
  if( natoms != ntot )lresize = 1;

  if( lresize ){
    // Resize FTOT
    memory->destroy( ftot );
    memory->destroy( xtot );
    memory->destroy( vtot );
    memory->destroy( order_tot );
    natoms = atom->natoms;
    nat = natoms;
    memory->create( ftot, natoms, 3, "fix/artn:ftot");
    memory->create( xtot, natoms, 3, "fix/artn:xtot");
    memory->create( vtot, natoms, 3, "fix/artn:vtot");
    memory->create( order_tot, natoms, "fix/artn:order_tot");
    lresize = 0 ;
  }



  // ...Collect the position and force
  // Collect_Force( nloc, f, nat, ftot )
  double **fbuff;
  int *ibuff;
  const int iat = 251;
  if( !me ){

    int n0(0), nmax(0);
    for( int ipc(0); ipc < nproc; ipc++ ) nmax = max( nmax, nloc[ipc] );
    memory->create( fbuff, nmax, 3, "fix/artn:fbuff");
    memory->create( ibuff, nmax, "fix/artn:ibuff");
   


    // ...We fill with the locale force of proc = 0

    //cout<< " * Comm at the beginning"<<endl;
    for( int i(0); i < nloc[me]; i++ ){
         int idx = n0 + i;
         ftot[idx][0] = f[i][0];
         ftot[idx][1] = f[i][1];
         ftot[idx][2] = f[i][2];

         xtot[idx][0] = tau[i][0];
         xtot[idx][1] = tau[i][1];
         xtot[idx][2] = tau[i][2];

         vtot[idx][0] = vel[i][0];
         vtot[idx][1] = vel[i][1];
         vtot[idx][2] = vel[i][2];

         order_tot[idx] = order[i];

         //if( order[i]== iat )cout<< " * "<<order[i] <<" x: "<< xtot[idx][0] <<" v: "<< vtot[idx][0] <<" f: "<< ftot[idx][0]<<endl;
    }
    n0 += nloc[me];

  //if( istep == 1 )exit(0);

    // ...We fill woth the local force from proc > 0

    for( int iproc(1); iproc < nproc; iproc++ ){
      cout<<" * ENDTER IN IPROC LOOP 0"<<endl;
      int n = nloc[iproc];
      MPI_Recv( fbuff, 3*n, MPI_DOUBLE, iproc, 0, world, MPI_STATUS_IGNORE );
      MPI_Recv( ibuff, n, MPI_INT, iproc, 0, world, MPI_STATUS_IGNORE );

      // collect force in desorder
      for( int i(0); i < n; i++){
         int idx = n0 + i;
         ftot[idx][0] = fbuff[i][0];
         ftot[idx][1] = fbuff[i][1];
         ftot[idx][2] = fbuff[i][2];

         order_tot[idx] = ibuff[i];
      }

      // collect position in desorder
      MPI_Recv( fbuff, 3*n, MPI_DOUBLE, iproc, 10, world, MPI_STATUS_IGNORE );
      for( int i(0); i < n; i++){
         int idx = n0 + i;
         xtot[idx][0] = fbuff[i][0];
         xtot[idx][1] = fbuff[i][1];
         xtot[idx][2] = fbuff[i][2];
      }

      // collect velocity in desorder
      MPI_Recv( fbuff, 3*n, MPI_DOUBLE, iproc, 20, world, MPI_STATUS_IGNORE );
      for( int i(0); i < n; i++){
         int idx = n0 + i;
         vtot[idx][0] = fbuff[i][0];
         vtot[idx][1] = fbuff[i][1];
         vtot[idx][2] = fbuff[i][2];
      }
      n0 += n;

    }


    memory->destroy( fbuff );
    memory->destroy( ibuff );
      

  }else{

    // ...All the proc send own array to proc = 0
    MPI_Send( &f[0][0], 3*nloc[me], MPI_DOUBLE, 0, 0, world );
    MPI_Send( order, nloc[me], MPI_INT, 0, 0, world );

    MPI_Send( &tau[0][0], 3*nloc[me], MPI_DOUBLE, 0, 10, world );
    MPI_Send( &vel[0][0], 3*nloc[me], MPI_DOUBLE, 0, 20, world );
  }


  lconv = false;
  // ...ARTn
  //artn_( &f[0][0], &etot, nat, ityp, elt, &tau[0][0], order, &lat[0][0], &if_pos[0][0], &disp, &lconv );
  if( !me )artn_( &ftot[0][0], &etot, nat, ityp, elt, &xtot[0][0], order_tot, &lat[0][0], &if_pos[0][0], &disp, &lconv );




  // ----------------------------------------------------------------------COMVERGENCE 
  int iconv = int(lconv);
  MPI_Bcast( &iconv, 1, MPI_INT, 0, world );
  //MPI_Bcast( &lconv, 1, MPI_C_BOOL, 0, world );
  if( iconv ){

    // ...Reset the energy force tolerence
    update-> etol = etol;
    update-> ftol = ftol;

    // ...Spread the force 
    //Spread_Force( nloc, &f[0][0], &ftot[0][0] );
    
    if( !me ){

      // ...We fill with me = 0
      for( int i(0); i < nloc[me]; i++){
        f[i][0] = ftot[i][0];
        f[i][1] = ftot[i][1];
        f[i][2] = ftot[i][2];

        tau[i][0] = xtot[i][0];
        tau[i][1] = xtot[i][1];
        tau[i][2] = xtot[i][2];
      }

    
      // ...We send the part of ftot specific to iproc 
      int n0 = 0;
      MPI_Request send_request;
      for( int iproc(1); iproc < nproc; iproc++ ){
        cout<<" * ENDTER IN IPROC LOOP 1"<<endl;
        n0 += nloc[iproc-1];
        //MPI_Send( &ftot[n0][0], 3*nloc[iproc], MPI_DOUBLE, iproc, 0, world );
        MPI_Isend( &ftot[n0][0], 3*nloc[iproc], MPI_DOUBLE, iproc, 0, world, &send_request );
        MPI_Isend( &xtot[n0][0], 3*nloc[iproc], MPI_DOUBLE, iproc, 10, world, &send_request );
      }

      //memory->destroy( ftot );

    }else{

      int n = nloc[me];
      MPI_Recv( &f[0][0], 3*n, MPI_DOUBLE, 0, 0, world, MPI_STATUS_IGNORE );
      MPI_Recv( &tau[0][0], 3*n, MPI_DOUBLE, 0, 10, world, MPI_STATUS_IGNORE );

    }


    cout<< " ************************** ARTn CONVERGED"<<endl;
    return;
  } // --------------------------------------------------------------------------------



  // POST ARTn/PRE MOVE_MODE
  //cout<< " * POST_ARTn/PRE_MOVE_MODE::DISP::"<< disp<<endl;


  // ...Convert the movement to the force
  //move_mode_( nat, &f[0][0], &vel[0][0], &etot, &nsteppos, &dt_curr, &alpha, &alpha_init, &dt_init, &disp );
  if( !me )move_mode_( nat, &ftot[0][0], &vtot[0][0], &etot, &nsteppos, &dt_curr, &alpha, &alpha_init, &dt_init, &disp );


  // ...Spread the force 
  //Spread_Force( nloc, &f[0][0], &ftot[0][0] );

  if( !me ){

    // ...We fill with me = 0
    //cout<< " * Comm After move_mode"<<endl;
    for( int i(0); i < nloc[me]; i++){
      f[i][0] = ftot[i][0];
      f[i][1] = ftot[i][1];
      f[i][2] = ftot[i][2];
 
      tau[i][0] = xtot[i][0];
      tau[i][1] = xtot[i][1];
      tau[i][2] = xtot[i][2];
    
      vel[i][0] = vtot[i][0];
      vel[i][1] = vtot[i][1];
      vel[i][2] = vtot[i][2];

      //if( order[i] == iat )cout<< disp << " * "<<order[i] <<" x: "<< xtot[i][0] <<" v: "<< vtot[i][0] <<" f: "<< ftot[i][0]<<endl;
      //if( order[i] == iat )cout<< disp << " * "<<order[i] <<" x: "<< tau[i][0] <<" v: "<< vel[i][0] <<" f: "<< f[i][0]<<endl;
    }


    // ...We send the part of ftot specific to iproc 
    int n0 = 0;
    MPI_Request send_request;
    for( int iproc(1); iproc < nproc; iproc++ ){
       cout<<" * ENDTER IN IPROC LOOP 2"<<endl;
       n0 += nloc[iproc-1];
       //MPI_Send( &ftot[n0][0], 3*nloc[iproc], MPI_DOUBLE, iproc, 0, world );
       MPI_Isend( &ftot[n0][0], 3*nloc[iproc], MPI_DOUBLE, iproc, 0, world, &send_request );
       MPI_Isend( &xtot[n0][0], 3*nloc[iproc], MPI_DOUBLE, iproc, 10, world, &send_request );
       MPI_Isend( &vtot[n0][0], 3*nloc[iproc], MPI_DOUBLE, iproc, 20, world, &send_request );
    }

    //memory->destroy( ftot );

  }else{

    int n = nloc[me];
    MPI_Recv( &f[0][0], 3*n, MPI_DOUBLE, 0, 0, world, MPI_STATUS_IGNORE );
    MPI_Recv( &tau[0][0], 3*n, MPI_DOUBLE, 0, 10, world, MPI_STATUS_IGNORE );
    MPI_Recv( &vel[0][0], 3*n, MPI_DOUBLE, 0, 20, world, MPI_STATUS_IGNORE );
  
  }



  //cout<<" * POST_MOVE_MODE::"<<endl;
  // ...Convert to the LAMMPS units
  if( !(disp == __artn_params_MOD_perp || disp == __artn_params_MOD_relx) ){

    // Comvert the force Ry to LAMMPS units
    if( rmass ){
      for( int i(0); i<nat; i++){
        f[i][0] *= rmass[i]; 
        f[i][1] *= rmass[i]; 
        f[i][2] *= rmass[i]; 
      }
    }else{
      for( int i(0); i<nat; i++){
	f[i][0] *= mass[ityp[i]];      
	f[i][1] *= mass[ityp[i]];      
	f[i][2] *= mass[ityp[i]];
      }
    }

  }




  // CHANGE FIRE PARAMETER as function of the value of 

  if( (disp == __artn_params_MOD_perp && __artn_params_MOD_iperp == 1) 
    || (disp != __artn_params_MOD_perp && disp != __artn_params_MOD_relx) ){

    //cout<< " FIX_ARTN::Params Actuallization"<<endl;

    // ...Update the time
    update->dt = dt_curr;
    string str;


    // ...Allocate/Deallocate word
    nword = 6;
    if( word )memory->destroy( word );
    memory->create( word, nword, 20, "fix:word" );

    strcpy( word[0], "alpha0" );
    str = to_string(alpha); 
    strcpy( word[1], str.c_str() );

    //strcpy( word[2], "alphashrink") ;
    //strcpy( word[3], "1.0" ) ;

    strcpy( word[2], "delaystep" );
    str = to_string(nsteppos);
    strcpy( word[3], str.c_str() );

    // ...RELAX step -> halfstepback = yes
    if( disp == __artn_params_MOD_relx || disp == __artn_params_MOD_perp ){
      strcpy( word[4], "halfstepback" );
      strcpy( word[5], "yes" );
    } else {
      strcpy( word[4], "halfstepback" );
      strcpy( word[5], "no" );
    }

    // ...Print the command
    //cout<<" * -> Displacement: "<< disp <<endl;
    //cout<<" * -> dt_curr: "<< update->dt <<endl;
    //for( int i(0); i < nword; i = i+2 )
    //  cout<<" * -> "<< word[i] << ": "<< word[i+1] <<endl;


    // ...Send the new parameter to minmize
    minimize-> modify_params( nword, word );
    minimize-> init();
  }



  // ...Compute V.F: Allow to know if the next call come from the vdotf < 0 condition
  // or as normally after the integration step
  vdotf = 0.0;
  //for (int i = 0; i < nlocal; i++)
  for( int i(0); i < nloc[me]; i++ )
    vdotf += vel[i][0]*f[i][0] + vel[i][1]*f[i][1] + vel[i][2]*f[i][2];
  MPI_Allreduce( &vdotf, &vdotfall, 1, MPI_DOUBLE, MPI_SUM, world );

  if( !(vdotfall > 0) )nextblank = 1;
  if( istep == 0 )nextblank = 0;



  // ...Store the actual force
  //cout<< " * Before to leave "<< endl;
  //for( int i(0); i < nlocal; i++ ){
  for( int i(0); i < nloc[me]; i++ ){
    f_prev[i][0] = f[i][0];
    f_prev[i][1] = f[i][1];
    f_prev[i][2] = f[i][2];
    //if( order[i] == iat )cout<< " * "<<order[i] <<" x: "<< tau[i][0] <<" v: "<< vel[i][0] <<" f: "<< f[i][0]<<endl;
  }




  // ...Increment & return
  istep++;
  return;

}








