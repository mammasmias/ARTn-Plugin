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
  //void artn_( double **force, double etot, int nat, int *ityp, char *elt, double **tau, double lat[3][3], int *if_pos, char *move, bool lconv );
  //void move_mode_( int nat, double **force, double **vel, double etot, int nsteppos, double dt_curr, double alpha, double alpha_init, double dt_init, char *cmode );
  //void artn_( double **force, double* etot, int* nat, int *ityp, char *elt, double **tau, double lat[3][3], int *if_pos, int* disp, bool* lconv );
  void artn_( double **const force, double* etot, const int nat, const int *ityp, const char *elt, double **tau, const double lat[3][3], const int *if_pos, int* disp, bool* lconv );
  void move_mode_( const int nat, double **force, double **vel, double* etot, int* nsteppos, double* dt_curr, double* alpha, const double* alpha_init, const double* dt_init, int* disp );
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

  int nword = 8;
  //char word[6][50];
  char **word;
  word = new char*[nword];

  word[0] = "alpha0";
  word[1] = "0.0";
  word[2]= "alphashrink" ;
  word[3]= "1.0" ;
  word[4]= "delaystep" ;
  word[5]= "5" ;
  word[6]= "halfstepback" ;
  word[7]= "no" ;

  minimize-> modify_params( nword, word );
  //minimize-> modify_params( 0, word );

  minimize-> setup_style();

  for( int i(0); i<nword; i++)
    word[i] = NULL;
  delete [] word;


  /*
  -- Save the initial Fire Parameter 
  */
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

  double **x = atom->x;
  double **f = atom->f;
  int nat = atom->natoms;
  for( int i(0); i<10; i++){
    cout<< " * pos: "<< i<< " | "<< x[i][0]<< " "<<x[i][1]<< " "<<x[i][2]<< " "<<endl;
    //cout<< " * frc: "<< i<< " | "<<f[i][0]<< " "<<f[i][1]<< " "<<f[i][2]<< " "<<endl;
  }

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

  char *elt ;
  elt = new char[nat];
  for( int i(0); i < nat; i++ ) elt[i] = alphab[ ityp[i] ];

  double **tau = atom->x;

  //Build it with : domain-> boxlo[3], boxhi[3] and xy, xz, yz
  double lat[3][3];
  double dx = domain->boxhi[0] - domain->boxlo[0], 
         dy = domain->boxhi[1] - domain->boxlo[1], 
         dz = domain->boxhi[2] - domain->boxlo[2];

  lat[0][0] = dx;     lat[0][1] = domain->xy    ; lat[0][2] = domain->xz ;
  lat[1][0] = 0.0 ;   lat[1][1] = dy            ; lat[1][2] = domain->yz ;
  lat[2][0] = 0.0 ;   lat[2][1] = 0.0           ; lat[2][2] = dz ;

  //double alat = see the units of the simulation real metal SI : UPDATE->unit_style
  //int istep = maybe we don't need
  int *if_pos;
  if_pos = new int[nat];
  for( int i(0); i < nat; i++ ) if_pos[ i ] = 1;
  double **vel = atom->v ;
  

  int nsteppos;
  double dt_curr;// = update->dt;
  //double fire_alpha_init = private variable of Min_Fire;
  double alpha;// = update->minimize->alpha0;
  //double fire_alphashrink = update->minimize->alphashrink;


  //bool lconv = boolian variable - take care of the Fortran/C++ interface
  bool lconv;
  int disp;
  //char *move;
  //char* prefix = prefix for scratch files of engine
  //char* prefix;
  //char* tmp_dir = scratch directory of engine 
  //char* tmp_dir;


  // PRE ARTn
  cout<< " * PRE_ARTn:: "<< nat <<endl;

  //artn( force, etot, forc_conv_thr_qe, nat, ityp, atm, tau, at, alat, istep, if_pos, vel, dt, fire_alpha_init, lconv, prefix, tmp_dir )
  //artn_( force, etot, nat, ityp, elt, tau, lat, if_pos, move, lconv );
  //artn_( force, &etot, &nat, ityp, elt, tau, lat, if_pos, &disp, &lconv );
  artn_( force, &etot, nat, ityp, elt, tau, lat, if_pos, &disp, &lconv );

  delete [] if_pos;

  //return;

  // POST ARTn/PRE MOVE_MODE
  cout<< " * POST_ARTn/PRE_MOVE_MODE::"<<endl;
  for( int i(0); i<10; i++){
    cout<< " * pos: "<< i<< " | "<< tau[i][0]<< " "<<tau[i][1]<< " "<<tau[i][2]<< " "<<endl;
    cout<< " * frc: "<< i<< " | "<<force[i][0]<< " "<<force[i][1]<< " "<<force[i][2]<< " "<<endl;
  }

  //disp = 5; //RELX
  //move_mode_( nat, force, vel, fire_alpha_init, dt );
  //move_mode_( nat, force, vel, etot, nsteppos, dt_curr, alpha, alpha_init, dt_init, move );
  move_mode_( nat, force, vel, &etot, &nsteppos, &dt_curr, &alpha, &alpha_init, &dt_init, &disp );

  if( !force )cout<<" * FORCE loose the address"<<endl;

  // POST MOVE_MODE
  cout<< " * POST_MOVE_MODE::alpha "<< alpha <<" | "<< etot << " | "<< dt_curr<<" | "<<force[1][0]<<endl;
  for( int i(0); i<10; i++){
    cout<< " * pos: "<< i<< " | "<< tau[i][0]<< " "<<tau[i][1]<< " "<<tau[i][2]<< " "<<endl;
    //cout<< " * frc: "<< i<< " | "<<force[i][0]<< " "<<force[i][1]<< " "<<force[i][2]<< " "<<endl;
  }

  return;

  // CHANGE FIRE PARAMETER
  int nword = 4;
  char **word;
  word = new char*[nword];

  word[0] = "alpha0";
  word[1] = "0.0";
  word[2] = "alphashrink" ;
  word[3] = "1.0" ;
  //word[4]= "delaystep" ;
  //word[5]= "5" ;
  //word[6]= "halfstepback" ;
  //word[7]= "no" ;

  //minimize-> modify_params( nword, word );

  for( int i(0); i<nword; i++) word[i] = NULL;
  delete [] word;

  return;

}








