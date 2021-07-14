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
extern int __artn_params_MOD_iperp;
extern int __artn_params_MOD_perp;
extern int __artn_params_MOD_relx;

/* ---------------------------------------------------------------------- */

FixARTn::FixARTn( LAMMPS *lmp, int narg, char **arg ): Fix( lmp, narg, arg ) 
{

  if (narg < 3) error->all(FLERR,"Illegal fix ARTn command");

  istep = 0;
  nword = 0;
  word = nullptr;

  nextblank = 0;
  f_prev = nullptr;

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

  pe_compute = nullptr;

}


/* ---------------------------------------------------------------------- */

FixARTn::~FixARTn() {

  /* deallocate the array */
  memory->destroy( word );
  memory->destroy( order );
  memory->destroy( f_prev );
  
  pe_compute = nullptr;

}


/* ---------------------------------------------------------------------- */

int FixARTn::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  //mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixARTn::init() {

  // Check if the other nodule/fix whatever are well define

  // compute for potential energy

  int id = modify->find_compute("thermo_pe");
  if (id < 0) error->all(FLERR,"Minimization could not find thermo_pe compute");
  pe_compute = modify->compute[id];

  istep = 0;

  nword = 6;
  //memory->create( word, nword, 20, "fix:word");

}


/* ---------------------------------------------------------------------- */

void FixARTn::setup( int /*vflag*/ ) {

  // Depending of the module environment we can change the move interface
  // and other thing

  cout<<" * IN SETUP..." << endl;


}

/* ---------------------------------------------------------------------- */

void FixARTn::min_setup( int vflag ) {

  // Call here if it is needed - To confirm
  cout<<" * IN MIN_SETUP..." << endl;
  //post_force( vflag );

  class Min *minimize = update-> minimize;

  cout<< " * CHANGE PARAM..."<<endl;

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


  // ...Save the Fire Parameter - Init:
  alpha_init = 0.1;
  alphashrink = 0.99;

  dt_init = 0.1; //update->dt;
  dtsk = 0.5;
  dtgrow = 1.1;

  dmax = 0.1;

  tmax = 20.;
  tmin = 0.02;
  dtmax = tmax*dt_init;
  dtmin = tmin*dt_init;

  fire_integrator = 0;
  ntimestep_start = update->ntimestep;
  
  etol = update->etol;
  ftol = update->ftol;

  // ...Print the new Parameters
  minimize-> setup_style();

  // ...Print the Initial Fire Parameters
  cout<< " * Alpha_init->"<< alpha_init<< endl;
  cout<< " * dt_init->"<< dt_init<< endl;
  cout<< " * dtmin->"<< dtmin<< endl;
  cout<< " * dtmax->"<< dtmax<< endl;


  // Copy the order of atom
  int nat = atom-> natoms;
  tagint *itag = atom->tag;
  memory->create( order, nat, "Fix_artn::order" );
  for( int i(0); i < nat; i++ ) order[i] = itag[i];


  // Allocate the previous force table
  memory->create( f_prev, nat, 3, "fix_artn:f_prev");
  nextblank = 0;


  

}


/* ---------------------------------------------------------------------- */

void FixARTn::min_post_force( int vflag ) {

  cout<<" * IN MIN_POST_FORCES..." << endl;
  post_force( vflag );

  cout<<" * OUT MIN_POST_FORCES..." << endl;

}



/* ---------------------------------------------------------------------- */

void FixARTn::post_force( int /*vflag*/ ){

  // call pARTn library...

  cout<<" * IN_POST_FORCES..." << endl;
  cout<<" * MIN::"<< update-> minimize_style<< " NTIMESTEP:: "<< update-> ntimestep << endl;

  class Min *minimize = update-> minimize; 

  if( !minimize )cout<<" * Min_vector is not linked"<<endl;

  //cout<<" * FIRE Param::"<< minimize->dt<<" | "<< minimize->alpha << endl;
  //cout<< " * FIX_ARTn::MINIMIZE CONV:"<< update-> etol<< " | "<< update->ftol<< endl;

  update->etol = 1e-12;
  update->etol = 1e-12;
  cout<< " * FIX_ARTn::MINIMIZE CONV:"<< update-> etol<< " | "<< update->ftol<< endl;

  //cout<< " * Param dt is accessible by update-> dt "<< update->dt << endl;
  //cout<< " * FIRE Params can be change through update->minimize->modify_params( int narg, char **args )"<<endl;


  double **tau = atom->x;
  double **f = atom->f;
  double **vel = atom->v ;
  const int *ityp = atom->type;


  // ...Comput V.F to know which call it is
  
  int nlocal = atom->nlocal;
  double vdotf = 0.0, vdotfall;

  for (int i = 0; i < nlocal; i++){
    //cout<< " * FIX_ARTn::position: "<< i<< "|"<< order[i]<< " | "<< ityp[i]<<" "<< tau[i][0]<< " "<< tau[i][1]<<" "<< tau[i][2]<<endl;
    vdotf += vel[i][0]*f[i][0] + vel[i][1]*f[i][1] + vel[i][2]*f[i][2];
  }
  MPI_Allreduce(&vdotf,&vdotfall,1,MPI_DOUBLE,MPI_SUM,world);
  cout<<" * FIX_ARTN::Computed v.f = "<< vdotfall<<endl;



  // ...If v.f is 0 or under min_fire call the force to ajust 
  // the integrator step parameter dtv
  if( !(vdotfall > 0) && update-> ntimestep > 1 && nextblank ){


    // ...Rescale the force if the dt has been change 
    double rescale = dt_curr / update->dt;

    // ...Reload the previous ARTn-force 
    for( int i(0); i < nlocal; i++){
      f[i][0] = f_prev[i][0] * rescale*rescale ;
      f[i][1] = f_prev[i][1] * rescale*rescale ;
      f[i][2] = f_prev[i][2] * rescale*rescale ;
    }

    cout<< " ********* RETURN WITHOUT COMPUTE ARTn | v.f = "<< vdotfall<< " | Rescale "<< rescale<<" "<< dt_curr<< " " <<update->dt<<endl;

    // ...Next call should be after the integration
    nextblank = 0;
    return;

  }

  /*...............................................................................*/





  double etot = pe_compute->compute_scalar();
  const int nat = atom-> natoms;

  // Conv. Criterium
  //double epse = update-> etol;
  //double epsf = update-> ftol;


  /*
     ...Array of element */
  char *elt ;
  elt = new char[nat];
  for( int i(0); i < nat; i++ ) elt[i] = alphab[ ityp[i] ];


  /*
    ...Build it with : domain-> boxlo[3], boxhi[3] and xy, xz, yz  */
  double lat[3][3];
  double dx = domain->boxhi[0] - domain->boxlo[0], 
         dy = domain->boxhi[1] - domain->boxlo[1], 
         dz = domain->boxhi[2] - domain->boxlo[2];

  lat[0][0] = dx;     lat[0][1] = domain->xy    ; lat[0][2] = domain->xz ;
  lat[1][0] = 0.0 ;   lat[1][1] = dy            ; lat[1][2] = domain->yz ;
  lat[2][0] = 0.0 ;   lat[2][1] = 0.0           ; lat[2][2] = dz ;


  /* 
    ...Array to fix the atom - could find a fix for that */
  int **if_pos;
  memory->create(if_pos,nat,3,"fix:if_pos");
  for( int i(0); i < nat; i++ ){
     if_pos[ i ][0] = 1;
     if_pos[ i ][1] = 1;
     if_pos[ i ][2] = 1;
     //cout<< " if_pos "<< if_pos[i][0]<< " | "<< if_pos[i][2]<< " | "<< if_pos[i][2] <<endl;
  }
  

  /* ...Convergence and displacement */
  bool lconv;
  int disp;

  
  // ...Store the actual force
  for( int i(0); i < nat; i++ ){
    f_prev[i][0] = f[i][0];
    f_prev[i][1] = f[i][1];
    f_prev[i][2] = f[i][2];
  }

  // PRE ARTn
  double *rmass = atom->rmass;
  double *mass = atom->mass;

  double EA2RB = eV2Ry / Ang2Bohr;
  for( int i(0); i<nat; i++){
    f[i][0] *= EA2RB;  
    f[i][1] *= EA2RB;  
    f[i][2] *= EA2RB;
    if( i == 243 )cout<< " Force("<<i<<"):: "<< f[i][0]<< " "<< f[i][1]<< " "<< f[i][2]<< endl;  
  }
/*  if( rmass ){
    for( int i(0); i<nat; i++){
      f[i][0] *= EA2RB / rmass[i]; 
      f[i][1] *= EA2RB / rmass[i]; 
      f[i][2] *= EA2RB / rmass[i]; 
    }
  }else{
    for( int i(0); i<nat; i++){
      f[i][0] *= EA2RB / mass[ityp[i]]; 
      f[i][1] *= EA2RB / mass[ityp[i]]; 
      f[i][2] *= EA2RB / mass[ityp[i]]; 
    }
  }
*/
  cout<< " * PRE_ARTn::Force convertion::"<< nat<< " | F *= "<< EA2RB <<endl;

  artn_( &f[0][0], &etot, nat, ityp, elt, &tau[0][0], order, &lat[0][0], &if_pos[0][0], &disp, &lconv );


  // Maybe should be a global variable
  memory->destroy(if_pos);

  // POST ARTn/PRE MOVE_MODE
  cout<< " * POST_ARTn/PRE_MOVE_MODE::DISP::"<< disp<<endl;

  double dt0 = dt_init*ps2aut ;
  dt_curr *= ps2aut;
  
  cout<< " * fixArtn:alpha_0 "<< alpha_init<< " | dt "<< dt_init<< endl;

  move_mode_( nat, &f[0][0], &vel[0][0], &etot, &nsteppos, &dt_curr, &alpha, &alpha_init, &dt0, &disp );

  dt_curr /= ps2aut ;

  cout<< " * fixArtn:alpha "<< alpha<< " | dt "<< dt_curr<< " | "<< disp << endl;


  double RB2EA = 1. / EA2RB ;// / force->ftm2v;
  if( disp == __artn_params_MOD_perp || disp == __artn_params_MOD_relx ){

    for( int i(0); i<nat; i++){
        if( i == 243 )cout<< " Force("<<i<<"):: "<< f[i][0]<< " "<< f[i][1]<< " "<< f[i][2]<< endl;  
        f[i][0] *= RB2EA;
        f[i][1] *= RB2EA;
        f[i][2] *= RB2EA;
        if( i == 243 )cout<< " Force("<<i<<"):: "<< f[i][0]<< " "<< f[i][1]<< " "<< f[i][2]<< endl;  
    }

  }else{

    // Comvert the force Ry to LAMMPS units
    if( rmass ){
      for( int i(0); i<nat; i++){
        f[i][0] *= RB2EA*rmass[i]; 
        f[i][1] *= RB2EA*rmass[i]; 
        f[i][2] *= RB2EA*rmass[i]; 
      }
    }else{
      for( int i(0); i<nat; i++){
	f[i][0] *= RB2EA*mass[ityp[i]];      
	f[i][1] *= RB2EA*mass[ityp[i]];      
	f[i][2] *= RB2EA*mass[ityp[i]];
      }
    }

  }



  cout<< " FIX_ARTN::Convert:: "<< EA2RB << " | "<< RB2EA <<endl;
  int iat = 243;
  cout<< " Force("<<iat<<"):: "<< f[iat][0]<< " "<< f[iat][1]<< " "<< f[iat][2]<< endl; 
  double dt2 = dt_curr*dt_curr/mass[ityp[iat]]*force->ftm2v;
  cout<< " v("<<iat<<"):: "<< vel[iat][0]<< " "<< vel[iat][1]<< " "<< vel[iat][2]<< std::endl; 
  cout<< " dx("<<iat<<"):: "<< f[iat][0]*dt2<< " "<< f[iat][1]*dt2<< " "<< f[iat][2]*dt2<< std::endl; 


  // ...Store the actual force
  for( int i(0); i < nat; i++ ){
    cout<< " ** DIFF_FORCE::"<<i<<" :: "<<f_prev[i][0] - f[i][0]<< " | "<<f_prev[i][1] - f[i][1]<< " | "<<f_prev[i][2] -f[i][2]<<endl;
  }



  // ...Compute V.F
/*  vdotf = 0.0;
  for (int i = 0; i < nlocal; i++)
    vdotf += vel[i][0]*force[i][0] + vel[i][1]*force[i][1] + vel[i][2]*force[i][2];
  MPI_Allreduce(&vdotf,&vdotfall,1,MPI_DOUBLE,MPI_SUM,world);
  cout<< " FIX_ARTN::VdotF:: "<< vdotfall<<endl;
*/



  // CHANGE FIRE PARAMETER as function of the value of 

  //if( (disp == 2 && __artn_params_MOD_iperp == 1) || disp != 2 ){
  if( (disp == __artn_params_MOD_perp && __artn_params_MOD_iperp == 1) 
    || (disp != __artn_params_MOD_perp && disp != __artn_params_MOD_relx) ){

    cout<< " FIX_ARTN::Params Actuallization"<<endl;

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
      //nword = nword + 2; 
      strcpy( word[4], "halfstepback" );
      strcpy( word[5], "yes" );
    } else {
      strcpy( word[4], "halfstepback" );
      strcpy( word[5], "no" );
    }

    // ...Print the command
    cout<<" * -> Displacement: "<< disp <<endl;
    cout<<" * -> dt_curr: "<< update->dt <<endl;
    for( int i(0); i < nword; i = i+2 )
      cout<<" * -> "<< word[i] << ": "<< word[i+1] <<endl;


    minimize-> modify_params( nword, word );
    minimize-> init();
  }

  // ...Print the new Parameters
  //minimize-> setup_style();



  // ...Compute V.F: Allow to know if the next call come from the vdotf < 0 condition
  // or as normally after the integration step
  vdotf = 0.0;
  for (int i = 0; i < nlocal; i++)
    vdotf += vel[i][0]*f[i][0] + vel[i][1]*f[i][1] + vel[i][2]*f[i][2];
  MPI_Allreduce(&vdotf,&vdotfall,1,MPI_DOUBLE,MPI_SUM,world);

  if( !(vdotfall > 0) )nextblank = 1;
  if( istep == 0 )nextblank = 0;

  cout<< " FIX_ARTN::END::VdotF:: "<< vdotfall<< " | NextBlank:: "<< nextblank <<endl;


  // ...Store the actual force
  for( int i(0); i < nat; i++ ){
    f_prev[i][0] = f[i][0];
    f_prev[i][1] = f[i][1];
    f_prev[i][2] = f[i][2];
  }

  // ...Store actual dt
  //dt_prev = 


  //if( istep == 10 )exit(0);

  istep++;
  return;

}








