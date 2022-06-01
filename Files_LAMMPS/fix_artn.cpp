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
#include "artn.h"

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

  tab_comm = NULL;

  // ...Set to null/0 all the variable of the class
  istep = 0;
  nword = 0;
  nmax = 0;
  word = nullptr;

  nextblank = 0;
  f_prev = nullptr;
  v_prev = nullptr;
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

  // ...Define delaystep for the relaxation
  nsteppos0 = 5;

  //dt_init = 0.001; //update->dt;
  dtsk = 0.5;
  dtgrow = 1.1;

  dmax = 0.1;

  tmax = 20.;
  tmin = 0.02;

  fire_integrator = 0;

  //if( narg == 3 )return;

  int iarg(3);
  while( iarg < narg ){

    /* Here we change the min_fire parameter 
    */
    
    if (strcmp(arg[iarg],"dmax") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal min_modify command");
      dmax = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if( strcmp(arg[iarg], "alpha0") == 0 ){
      if (iarg+2 > narg) error->all(FLERR,"Illegal min_modify command");
      alpha_init = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"delaystep") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal min_modify command");
      nsteppos0 = utils::numeric(FLERR,arg[iarg+1],false,lmp);
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
  memory->destroy( v_prev );
  memory->destroy( if_pos );

  memory->destroy( nloc );
  memory->destroy( ftot );
  memory->destroy( xtot );
  memory->destroy( vtot );
  memory->destroy( order_tot );
  
  memory->destroy(tab_comm);

  pe_compute = nullptr;

}


/* ---------------------------------------------------------------------- */

int FixARTn::setmask()
{
  int mask = 0;
  //mask |= POST_FORCE;
  //mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  mask |= POST_RUN;
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


  // Communication
  nmax = atom->nmax;
  memory->create(tab_comm,nmax,"pair:tab_comm");

}







/* ---------------------------------------------------------------------- */

void FixARTn::min_setup( int vflag ) {

  // Call here if it is needed - To confirm
  //cout<<" * IN MIN_SETUP..." << endl;
  //post_force( vflag );

  class Min *minimize = update-> minimize;

  if( !me )cout<< " * FIX/ARTn::CHANGE PARAM..."<<endl;

  /*
  -- Change & Save the initial Fire Parameter 
     Exept: delaystep_start_flag = 1 (ALWAYS)
  */

  nword = 12;
  if( word )memory->destroy( word );
  memory->create( word, nword, 20, "fix:word" );

  string str;
  strcpy( word[0], "alpha0" );
  str = to_string(alpha_init);
  strcpy( word[1], str.c_str() );
  //strcpy( word[1], "0.0");   // at the begining
  strcpy( word[2], "alphashrink") ;
  str = to_string(alphashrink);
  strcpy( word[3], str.c_str() );
  //strcpy( word[3], "0.99") ;
  strcpy( word[4], "delaystep") ;
  str = to_string(nsteppos0);
  strcpy( word[5], str.c_str() );
  //strcpy( word[5], "5") ;
  strcpy( word[6], "halfstepback") ;
  strcpy( word[7], "no") ;
  strcpy( word[8], "integrator") ;
  strcpy( word[9], "eulerimplicit") ;
  strcpy( word[10], "dmax") ;
  str = to_string(dmax);
  strcpy( word[11], str.c_str() );
  //strcpy( word[11], "0.15") ;

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
  if( !me ){
    cout<< " * Alpha0->"<< alpha_init<< endl;
    cout<< " * dt0->"<< dt_init<< endl;
    cout<< " * dtmin->"<< dtmin<< endl;
    cout<< " * dtmax->"<< dtmax<< endl;
    cout<< " * ftm2v->"<< force->ftm2v << endl;
    cout<< " * dmax->"<< dmax << endl;
    cout<< " * delaystep->"<< nsteppos0 << endl;
  }


  // ...Copy the order of atom
  int nat = atom-> natoms;
  natoms = nat; // Save for the rescale step
  int nlocal = atom->nlocal;
  oldnloc = nlocal;

  tagint *itag = atom->tag;
  memory->create( order, nlocal, "Fix/artn::order" );
  for( int i(0); i < nlocal; i++ ) order[i] = itag[i];


  // ...Allocate the previous force table :: IS LOCAL
  memory->create( f_prev, nlocal, 3, "fix/artn:f_prev" );
  memory->create( v_prev, nlocal, 3, "fix/artn:v_prev" );
  nextblank = 0;



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
  memory->create( nlresize, nproc, "fix/artn:nlresize" );
  memory->create( ftot, nat, 3, "fix/artn:ftot" );
  memory->create( xtot, nat, 3, "fix/artn:xtot" );
  memory->create( vtot, nat, 3, "fix/artn:xtot" );
  memory->create( order_tot, nat, "fix/artn:order_tot" );


}










/* ---------------------------------------------------------------------- */

void FixARTn::min_post_force( int /*vflag*/ ){

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



  // ...Update and share the local system size

  MPI_Allgather( &nlocal, 1, MPI_INT, nloc, 1, MPI_INT, world );
  int ntot(0), lresize(0);
  for( int ipc(0); ipc < nproc; ipc++ )ntot += nloc[ipc];




  // ------------------------------------------------------------------- RESIZE SYSTEM SIZE


  // ...Resize total system 

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



  // ...Resize local system

  // verification of local size
  lresize = ( nloc[me] != oldnloc );
  for( int ipc(0); ipc < nproc; ipc++ )nlresize[ me ] = 0;
  MPI_Allgather( &lresize, 1, MPI_INT, nlresize, 1, MPI_INT, world );
  ntot = 0;
  for( int ipc(0); ipc < nproc; ipc++ )ntot += nlresize[ ipc ];


  // ...One of the local size change
  if( ntot > 0 ){


    // ...Array of old local size
    int *oldloc;
    memory->create( oldloc, nproc, "fix/artn:oldloc" );
    MPI_Allgather( &oldnloc, 1, MPI_INT, oldloc, 1, MPI_INT, world );


    // ...Create temporary Arrays
    int *inew;
    memory->create( inew, natoms, "fix/artn:inew" );
    memory->create( istart, natoms, "fix/artn:istart" );
    memory->create( length, natoms, "fix/artn:length" );


    // ---------------------------------- Use AllGatherv for f_prev to ftot
    //                                                       v_prev to vtot
    // ...Starting point:
    for( int ipc(0); ipc < nproc; ipc++ )
      istart[ ipc ] = (ipc > 0) ? istart[ ipc - 1 ] + 3*oldloc[ ipc - 1 ] : 0 ;

    // ...Length of receiv buffer
    for( int ipc(0); ipc < nproc; ipc++ ) length[ ipc ] = 3*oldloc[ ipc ];


    MPI_Gatherv( &f_prev[0][0], 3*oldnloc, MPI_DOUBLE, 
                 &ftot[0][0], length, istart, MPI_DOUBLE, 0, world );

    MPI_Gatherv( &v_prev[0][0], 3*oldnloc, MPI_DOUBLE, 
                 &vtot[0][0], length, istart, MPI_DOUBLE, 0, world );


    // --------------------------------- Use AllGatherv for order to order_tot
    // ...Starting point:
    for( int ipc(0); ipc < nproc; ipc++ )
      istart[ ipc ] = ( ipc > 0 ) ? istart[ ipc - 1 ] + oldloc[ ipc - 1 ] : 0 ;

    MPI_Gatherv( order, oldnloc, MPI_INT, 
                    order_tot, oldloc, istart, MPI_INT, 0, world );


  
    // ...Resize order
    tagint *itag = atom->tag;
    memory->destroy( order );
    memory->create( order, nlocal, "fix/artn:order" );

    // ...Fill now order
    for( int i(0); i < nlocal; i++ )order[ i ] = itag[ i ];


    // --------------------------------- Use AllGatherv for order to inew
    // ...Starting point:
    for( int ipc(0); ipc < nproc; ipc++ )
      istart[ ipc ] = ( ipc > 0 ) ? istart[ ipc - 1 ] + nloc[ ipc - 1 ] : 0 ;

    MPI_Gatherv( order, nlocal, MPI_INT,
                    inew, nloc, istart, MPI_INT, 0, world );


    // ...Change the order of force 
    if( !me ){
      if( !xtot )printf(" ERROR - *XTOT => NULL \n");
      if( !ftot )printf(" ERROR - *FTOT => NULL \n");
      for( int i(0); i < natoms; i++ ){

        int j;
        for( j = 0; j < natoms; j++ )
          if( order_tot[ j ] == inew[ i ] ){
            xtot[i][0] = ftot[j][0];
            xtot[i][1] = ftot[j][1];
            xtot[i][2] = ftot[j][2];
            break;
          }
      }  
    } // ::: ME = 0


    // ...Resize f_prev
    memory->destroy( f_prev );
    memory->create( f_prev, nlocal, 3, "fix/artn:f_prev" ); 


    // ...Starting point:
    for( int ipc(0); ipc < nproc; ipc++ )
      istart[ ipc ] = (ipc > 0) ? istart[ ipc - 1 ] + 3*nloc[ ipc - 1 ] : 0 ;

    // ...Length of receiv buffer
    for( int ipc(0); ipc < nproc; ipc++ )length[ ipc ] = 3*nloc[ ipc ];

    MPI_Scatterv( &xtot[0][0], length, istart, MPI_DOUBLE, 
                  &f_prev[0][0], 3*nlocal, MPI_DOUBLE, 0, world );


/*-------
    // ...Change the order of velocity
    if( !me ){
      if( !xtot )printf(" ERROR - *XTOT => NULL \n");
      if( !vtot )printf(" ERROR - *VTOT => NULL \n");
      for( int i(0); i < natoms; i++ ){

        int j;
        for( j = 0; j < natoms; j++ )
          if( order_tot[ j ] == inew[ i ] )break;

        xtot[i][0] = vtot[j][0];
        xtot[i][1] = vtot[j][1];
        xtot[i][2] = vtot[j][2];
      }
    } // ::: ME = 0


    // ...Resize f_prev
    memory->destroy( v_prev );
    memory->create( v_prev, nlocal, 3, "fix/artn:v_prev" );

    MPI_Scatterv( &xtot[0][0], length, istart, MPI_DOUBLE, 
                  &v_prev[0][0], 3*nlocal, MPI_DOUBLE, 0, world );
------*/





    // ...Save the new value
    oldnloc = nloc[ me ];

    // ...Destroy temporary arrays
    memory->destroy( oldloc );
    memory->destroy( inew );
    memory->destroy( length );
    memory->destroy( istart );

  } // -------------------------------------------------------------------------- END RESIZE SYSTEM






  // ...Comput V.F to know in which part of alogrithm energy_force() is called
  double vdotf = 0.0, vdotfall;
  for( int i = 0; i < nlocal; i++ )
    vdotf += vel[i][0]*f[i][0] + vel[i][1]*f[i][1] + vel[i][2]*f[i][2];
  MPI_Allreduce( &vdotf, &vdotfall, 1, MPI_DOUBLE, MPI_SUM, world );
  //if(istep == 48)cout<< me<<" * FIX_ARTN::Computed v.f = "<< vdotfall<<" | step "<< istep <<endl;






  // ---------------------------------------------------------------------  
  // ...If v.f is 0 or under min_fire call the force to ajust 
  // the integrator step parameter dtv
  //cout<< me<< " * CONDITION: "<< vdotfall<< "  " << update->ntimestep<<"  " << nextblank<<endl;
  if( !(vdotfall > 0) && update-> ntimestep > 1 && nextblank ){

    // ...Rescale the force if the dt has been change 
    double rscl = dt_curr / update->dt;

    // ...Reload the previous ARTn-force 
    for( int i(0); i < nlocal; i++){
      f[i][0] = f_prev[i][0] * rscl*rscl ;
      f[i][1] = f_prev[i][1] * rscl*rscl ;
      f[i][2] = f_prev[i][2] * rscl*rscl ;
      //vel[i][0] = 0.0; //v_prev[i][0] * rscl ;
      //vel[i][1] = 0. ;
      //vel[i][2] = 0. ;
      //vel[i][1] = v_prev[i][1] * rscl ;
      //vel[i][2] = v_prev[i][2] * rscl ;
      //cout<< me<<" * Rescale Force: "<< i<<" f:"<< f[i][0]<<endl;
    }


    //cout<< " ********* RETURN WITHOUT COMPUTE ARTn | v.f = "<< vdotfall<< " | Rescale "<< rscl<<" "<< dt_curr<< " " <<update->dt<<endl;


    // ...Next call should be after the integration
    nextblank = 0;
    return;

  } // ---------------------------------------------------------------------------------------- RETURN






  if( atom->nmax > nmax ){
    memory->destroy(tab_comm);
    nmax = atom->nmax;
    memory->create(tab_comm,nmax,"pair:tab_comm");
  }



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



  



  // ...Collect the position and force
  int *typ_tot;
  memory->create( typ_tot, natoms, "fix/artn:typ_tot");
  Collect_Arrays( nloc, tau, vel, f, nat, xtot, vtot, ftot, order_tot, typ_tot );

  //for( int i = 0; i < natoms-1; i++) cout<< i << " Order " << order_tot[i] <<endl;


  // ...ARTn
  lconv = false;
  double **disp_vec;
  if( !me ){
    memory->create( disp_vec, natoms, 3, "fix/artn:disp_vec");
    //artn_( &ftot[0][0], &etot, nat, ityp, elt, &xtot[0][0], order_tot, &lat[0][0], &if_pos[0][0], &disp, &disp_vec[0][0], &lconv );
    artn_( &ftot[0][0], &etot, nat, typ_tot, elt, &xtot[0][0], order_tot, &lat[0][0], &if_pos[0][0], &disp, &disp_vec[0][0], &lconv );
  }
  memory->destroy( typ_tot );


  // ...Spread the ARTn_Step (DISP) & Convergence
  int iconv = int(lconv);
  MPI_Bcast( &iconv, 1, MPI_INT, 0, world );
  MPI_Bcast( &disp, 1, MPI_INT, 0, world );


  // ...Convert the movement to the force
  if( !me ){
    move_mode_( nat, order_tot, &ftot[0][0], &vtot[0][0], &etot, &nsteppos, &dt_curr, &alpha, &alpha_init, &dt_init, &disp, &disp_vec[0][0] );
    memory->destroy( disp_vec );
  }



  // ---------------------------------------------------------------------- COMVERGENCE 
  if( iconv ){

    // ...Reset the energy force tolerence
    update-> etol = etol;
    update-> ftol = ftol;

    // ...Spread the force 
    Spread_Arrays( nloc, xtot, vtot, ftot, nat, tau, vel, f );


    MPI_Barrier( world );
    if( !me )cout<< "     ************************** ARTn CONVERGED"<<endl;
    return;
  } // --------------------------------------------------------------------------------



  // ...Convert the movement to the force
  //if( !me ){
  //  move_mode_( nat, order_tot, &ftot[0][0], &vtot[0][0], &etot, &nsteppos, &dt_curr, &alpha, &alpha_init, &dt_init, &disp, &disp_vec[0][0] );
  //  memory->destroy( disp_vec );
  //}


  // ...Spread the FIRE parameters
  MPI_Bcast( &dt_curr, 1, MPI_DOUBLE, 0, world );
  MPI_Bcast( &alpha, 1, MPI_DOUBLE, 0, world );



  // ...Spread the force 
  Spread_Arrays( nloc, xtot, vtot, ftot, nat, tau, vel, f );


  // ...Convert to the LAMMPS units
  if( !(disp == get_perp_() || disp == get_relx_()) ){

    double *rmass = atom->rmass;
    double *mass = atom->mass;

    // Comvert the force Ry to LAMMPS units
    if( rmass ){
      for( int i(0); i < nloc[me]; i++){
        f[i][0] *= rmass[i]; 
        f[i][1] *= rmass[i]; 
        f[i][2] *= rmass[i]; 
      }
    }else{
      for( int i(0); i < nloc[me]; i++){
	f[i][0] *= mass[ityp[i]];      
	f[i][1] *= mass[ityp[i]];      
	f[i][2] *= mass[ityp[i]];
      }
    }

  }




  // ...CHANGE FIRE PARAMETER as function of the value of 

  //if( (disp == get_perp_() && get_iperp_() == 1) 
  //  || (disp != get_perp_() && disp != get_relx_()) ){


    // ...Update the time
    //update->dt = dt_curr;
    string str;


    // ...Allocate/Deallocate word
    nword = 6;
    if( word )memory->destroy( word );
    memory->create( word, nword, 20, "fix:word" );

    strcpy( word[0], "alpha0" );
    str = to_string(alpha); 
    strcpy( word[1], str.c_str() );

    strcpy( word[2], "delaystep" );
    if( nsteppos != 0 )nsteppos = nsteppos0;
    str = to_string(nsteppos);
    strcpy( word[3], str.c_str() );

    // ...RELAX step -> halfstepback = yes
    //if( disp == __artn_params_MOD_relx || disp == __artn_params_MOD_perp ){
    if( disp == get_relx_() || disp == get_perp_() ){
      strcpy( word[4], "halfstepback" );
      strcpy( word[5], "yes" );
    } else {
      strcpy( word[4], "halfstepback" );
      strcpy( word[5], "no" );
    }

  // ...Launch modification of FIRE parameter
  if( (disp == get_perp_() && get_iperp_() == 1) 
    || (disp != get_perp_() && disp != get_relx_()) ){

    // ...Update the time
    update->dt = dt_curr;

    //printf(" CHange Fire param: %d at %ld\n", disp, update->ntimestep );
    //printf(" dt %f \n %s %s \n %s %s \n %s %s \n ", dt_curr, word[0], word[1], word[2], word[3], word[4], word[5]);

    // ...Send the new parameter to minmize
    minimize-> modify_params( nword, word );
    minimize-> init();
  }



  // ...Compute V.F: Allow to know if the next call come from the vdotf < 0 condition
  // or as normally after the integration step
  vdotf = 0.0;
  for( int i(0); i < nloc[me]; i++ )
    vdotf += vel[i][0]*f[i][0] + vel[i][1]*f[i][1] + vel[i][2]*f[i][2];
  MPI_Allreduce( &vdotf, &vdotfall, 1, MPI_DOUBLE, MPI_SUM, world );

  if( !(vdotfall > 0) )nextblank = 1;
  if( istep == 0 )nextblank = 0;



  // ...Store the actual force/velocity
  for( int i(0); i < nloc[me]; i++ ){
    f_prev[i][0] = f[i][0];
    f_prev[i][1] = f[i][1];
    f_prev[i][2] = f[i][2];
   // v_prev[i][0] = vel[i][0];
   // v_prev[i][1] = vel[i][1];
   // v_prev[i][2] = vel[i][2];
  }



  // ...Increment & return
  istep++;
  return;

}





/* ---------------------------------------------------------------------- */

void FixARTn::post_run(){

  // End of the ARTn research - we reset the ARTn counters & flag
  if( !me )clean_artn_(); // Only proc 0

}







/* ============================================================================ COMMUNICATION */


void FixARTn::Collect_Arrays( int* nloc, double **x, double **v, double **f, int nat, double **xtot, double **vtot, double **ftot, int *order_tot, int *typ_tot ){



  // ...Alloc temporary memory
  memory->create( istart, nproc, "fix/artn:istart");
  memory->create( length, nproc, "fix/artn:length");

  // ...Starting point:
  for( int ipc(0); ipc < nproc; ipc++ )
    istart[ ipc ] = (ipc > 0) ? istart[ ipc - 1 ] + 3*nloc[ ipc - 1 ] : 0 ;

  // ...Length of receiv buffer
  for( int ipc(0); ipc < nproc; ipc++ )length[ ipc ] = 3*nloc[ ipc ];


  // ...Gatherv ftot, vtot, xtot

  MPI_Gatherv( &f[0][0], 3*nloc[me], MPI_DOUBLE,
               &ftot[0][0], length, istart, MPI_DOUBLE, 0, world );

  MPI_Gatherv( &v[0][0], 3*nloc[me], MPI_DOUBLE,
               &vtot[0][0], length, istart, MPI_DOUBLE, 0, world );

  MPI_Gatherv( &x[0][0], 3*nloc[me], MPI_DOUBLE,
               &xtot[0][0], length, istart, MPI_DOUBLE, 0, world );



  // ...Starting point:
  for( int ipc(0); ipc < nproc; ipc++ )
    istart[ ipc ] = (ipc > 0) ? istart[ ipc - 1 ] + nloc[ ipc - 1 ] : 0 ;

  // ...Gatherv ftot, vtot, xtot

  MPI_Gatherv( order, nloc[me], MPI_INT,
               order_tot, nloc, istart, MPI_INT, 0, world );

  int *ityp = atom->type;
  MPI_Gatherv( ityp, nloc[me], MPI_INT,
               typ_tot, nloc, istart, MPI_INT, 0, world );


  memory->destroy( istart );
  memory->destroy( length );

}

/* --------------------------------------------------------------------------------------------------------------------------------- */

void FixARTn::Spread_Arrays( int *nloc, double **xtot, double **vtot, double **ftot, int nat, double **x, double **v, double **f ){

  // ...Alloc temporary memory
  memory->create( istart, nproc, "fix/artn:istart");
  memory->create( length, nproc, "fix/artn:length");

  // ...Starting point:
  for( int ipc(0); ipc < nproc; ipc++ )
    istart[ ipc ] = (ipc > 0) ? istart[ ipc - 1 ] + 3*nloc[ ipc - 1 ] : 0 ;

  // ...Length of receiv buffer
  for( int ipc(0); ipc < nproc; ipc++ )length[ ipc ] = 3*nloc[ ipc ];


  // ...Scatter ftot, vtot, xtot

  MPI_Scatterv( &ftot[0][0], length, istart, MPI_DOUBLE,
                &f[0][0], 3*nloc[me], MPI_DOUBLE, 0, world );

  MPI_Scatterv( &vtot[0][0], length, istart, MPI_DOUBLE,
                &v[0][0], 3*nloc[me], MPI_DOUBLE, 0, world );

  MPI_Scatterv( &xtot[0][0], length, istart, MPI_DOUBLE,
                &x[0][0], 3*nloc[me], MPI_DOUBLE, 0, world );

  memory->destroy( istart );
  memory->destroy( length );


}


/*
// ---------------------------------------------------------------------- 

int FixARTn::pack_reverse_comm(int n, int first, double *buf)
{ 
  int i,m,last;
  
  m = 0; 
  last = first + n; 
  for (i = first; i < last; i++) {
    buf[m++] = tab_comm[i];
    //    if (i<first+3) printf ("[%d] prc %d %d buf_send = %f \n",me,i,m-1,buf[m-1]);
  }
  return m;
}

// ---------------------------------------------------------------------- 

void FixARTn::unpack_reverse_comm(int n, int *list, double *buf)
{ 
  int i,j,m;
  
  m = 0; 
  for (i = 0; i < n; i++) {
    j = list[i];
    //  tab_comm[j] += buf[m++];
    tab_comm[j] = buf[m++];
    //    if (j<3) printf ("[%d] %d urc %d %d buf_recv = %f \n",me,n,i,m-1,buf[m-1]);
  }
}

// ----------------------------------------------------------------------

void FixARTn::forward(double *tab)
{ 
  int i;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  
  for (i=0; i<nlocal+nghost; i++) tab_comm[i] = tab[i];
  
  comm->forward_comm_fix(this);
  
  for (i=0; i<nlocal+nghost; i++) tab[i] = tab_comm[i];
}

*/


