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
  void artn_( double *const force, double* etot, const int nat, const int *ityp, const char *elt, double *const tau, const double *lat, const int *if_pos, int* disp, bool* lconv );
  void move_mode_( const int nat, double *const force, double *const vel, double* etot, int* nsteppos, double* dt_curr, double* alpha, const double* alpha_init, const double* dt_init, int* disp );
}

/* ---------------------------------------------------------------------- */

FixARTn::FixARTn( LAMMPS *lmp, int narg, char **arg ): Fix( lmp, narg, arg ) 
{

  if (narg < 3) error->all(FLERR,"Illegal fix ARTn command");

  istep = 0;


}


/* ---------------------------------------------------------------------- */

FixARTn::~FixARTn() {

  /* deallocate the array */

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

  int nword = 10;
  //char word[6][50];
  char **word;
  word = new char*[nword];

  word[0] = "alpha0";
  word[1] = "0.0";   // at the begining
  word[2]= "alphashrink" ;
  word[3]= "0.99" ;
  word[4]= "delaystep" ;
  word[5]= "20" ;
  word[6]= "halfstepback" ;
  word[7]= "no" ;
  word[8]= "integrator" ;
  word[9]= "eulerimplicit" ;

  minimize-> modify_params( nword, word );


  // ...Save the Fire Parameter - Init:
  alpha_init = 0.25;
  alphashrink = 0.99;

  dt_init = update->dt;
  dtsk = 0.5;
  dtgrow = 1.1;

  tmax = 20.;
  tmin = 0.02;
  dtmax = tmax*dt_init;
  dtmin = tmin*dt_init;

  fire_integrator = 0;
  ntimestep_start = update->ntimestep;
  

  // ...Print the new Parameters
  minimize-> setup_style();


  // ...Deallocate char**
  for( int i(0); i<nword; i++)
    word[i] = NULL;
  delete [] word;


  // Fuck! Protected : Accessible only by derived class 
  // We have to define our fire parameters.
  //dt_init = update-> minimize-> dt ;
  //alpha_init = update-> minimize-> alpha0;
  

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
  cout<<" * MIN::"<< update-> minimize_style<< endl;

  class Min *minimize = update-> minimize; 

  if( minimize )cout<<" * Min_vector is linked"<<endl;

  //cout<<" * FIRE Param::"<< minimize->dt<<" | "<< minimize->alpha << endl;
  cout<< " * MINIMIZE CONV:"<< update-> etol<< " | "<< update->max_eval<< endl;

  cout<< " * Param dt is accessible by update-> dt "<< update->dt << endl;
  cout<< " * FIRE Params can be change through update->minimize->modify_params( int narg, char **args )"<<endl;



  double **force = atom->f;
  double etot = pe_compute->compute_scalar();

  // Conv. Criterium
  double epse = update-> etol;
  double epsf = update-> ftol;

  const int nat = atom-> natoms;
  const int *ityp = atom->type;

  /*
     Array of element */
  char *elt ;
  elt = new char[nat];
  for( int i(0); i < nat; i++ ) elt[i] = alphab[ ityp[i] ];

  double **tau = atom->x;

  /*
     Build it with : domain-> boxlo[3], boxhi[3] and xy, xz, yz  */
  double lat[3][3];
  double dx = domain->boxhi[0] - domain->boxlo[0], 
         dy = domain->boxhi[1] - domain->boxlo[1], 
         dz = domain->boxhi[2] - domain->boxlo[2];

  lat[0][0] = dx;     lat[0][1] = domain->xy    ; lat[0][2] = domain->xz ;
  lat[1][0] = 0.0 ;   lat[1][1] = dy            ; lat[1][2] = domain->yz ;
  lat[2][0] = 0.0 ;   lat[2][1] = 0.0           ; lat[2][2] = dz ;


  /* 
     Array to fix the atom - could find a fix for that */
  int *if_pos;
  if_pos = new int[nat];
  for( int i(0); i < nat; i++ ) if_pos[ i ] = 1;
  double **vel = atom->v ;
  



  //bool lconv = boolian variable - take care of the Fortran/C++ interface
  bool lconv;
  int disp;


  // PRE ARTn
  cout<< " * PRE_ARTn:: "<< nat <<endl;
  double EA2RB = eV2Ry / Ang2Bohr;
  for( int i(0); i<nat; i++){
    force[i][0] *= EA2RB; 
    force[i][1] *= EA2RB; 
    force[i][2] *= EA2RB; 
  }

  artn_( &force[0][0], &etot, nat, ityp, elt, &tau[0][0], &lat[0][0], if_pos, &disp, &lconv );


  // Maybe should be a global variable
  delete [] if_pos;


  // POST ARTn/PRE MOVE_MODE
  cout<< " * POST_ARTn/PRE_MOVE_MODE::"<<endl;

  int nsteppos;
  double dt_curr = 1.;// = update->dt;
  double alpha = 1.;// = update->minimize->alpha0;
  double dt0 = dt_init*ps2aut ;
  
  cout<< " * alpa_0 "<< alpha_init<< " | dt "<< dt_init<< endl;

  //move_mode_( nat, &force[0][0], &vel[0][0], &etot, &nsteppos, &dt_curr, &alpha, &alpha_init, &dt_init, &disp );
  move_mode_( nat, &force[0][0], &vel[0][0], &etot, &nsteppos, &dt_curr, &alpha, &alpha_init, &dt0, &disp );

  double RB2EA = 1./EA2RB;
  for( int i(0); i<nat; i++){
    force[i][0] *= RB2EA; 
    force[i][1] *= RB2EA; 
    force[i][2] *= RB2EA; 
  }
  dt_curr /= ps2aut ;

  for( int i(0); i<30; i++)
    cout<< " * After move::f "<< i<< " | "<<force[i][0]<<" "<< force[i][1]<<" "<< force[i][2]<<endl;
  


  // CHANGE FIRE PARAMETER as function of the value of 
  char **word;
  int nword = 6, n ;
  word = new char*[nword];
  string str;

  update->dt = dt_curr;

  word[0] = "alpha0";
  //strncpy(word[0],  "alpha0", sizeof("alpha0"));
  str = to_string(alpha); n = str.length();
  cout<<" str: "<< str<< " | "<< alpha<< " | "<< n<<endl;
  //strncpy(word[1], str.c_str(), n*sizeof(char));
  //strncpy(word[1], str.c_str(), sizeof(str.c_str()));
  word[1] = "0.0";

  //exit(0);

  //word[2] = "alphashrink" ;
  //word[3] = "1.0" ;

  word[2]= "delaystep" ;
  str = to_string(nsteppos);
  n = str.length();
  //strncpy( word[3], str.c_str(), n*sizeof(char) );
  word[3]= "5" ;

  // ...RELAX step -> halfstepback = yes
  if( disp == 5 ){
    //nword = nword + 2; 
    word[4] = "halfstepback" ;
    word[5] = "yes" ;
  } else {
    word[4] = "halfstepback" ;
    word[5] = "no" ;
  }

  for( int i(0); i < nword; i = i+2 )
    cout<< word[i] << ": "<< word[i+1] <<endl;

  minimize-> modify_params( nword, word );

  minimize-> init();

  // ...Print the new Parameters
  minimize-> setup_style();

  for( int i(0); i<nword; i++) word[i] = NULL;
  delete [] word;

  exit(0);
  return;

}








