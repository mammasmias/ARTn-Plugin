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

}


/* ---------------------------------------------------------------------- */

void FixARTn::setup( int /*vflag*/ ) {

  // Depending of the module environment we can change the move interface
  // and other thing

  cout<<" * IN SETUP..." << endl;

  class Min *minimize = update-> minimize;

  cout<< " * CHANGE PARAM..."<<endl;

  int nword = 6;
  //char word[6][50];
  char **word;
  word = new char*[nword];

  word[0] = "alpha0";
  word[1] = "0.0";
  word[2]= "alphashrink" ;
  word[3]= "1.0" ;
  word[4]= "delaystep" ;
  word[5]= "5" ;

  //minimize-> modify_params( nword, word );
  minimize-> modify_params( 0, word );

  for( int i(0); i<nword; i++)
    word[i] = NULL;
  delete [] word;


}

/* ---------------------------------------------------------------------- */

void FixARTn::min_setup( int vflag ) {

  // Call here if it is needed - To confirm
  cout<<" * IN MIN_SETUP..." << endl;
  //post_force( vflag );

  class Min *minimize = update-> minimize;

  cout<< " * CHANGE PARAM..."<<endl;

  int nword = 6;
  //char word[6][50];
  char **word;
  word = new char*[nword];

  word[0] = "alpha0";
  word[1] = "0.0";
  word[2]= "alphashrink" ;
  word[3]= "1.0" ;
  word[4]= "delaystep" ;
  word[5]= "5" ;

  minimize-> modify_params( nword, word );
  //minimize-> modify_params( 0, word );

  minimize-> setup_style();

  for( int i(0); i<nword; i++)
    word[i] = NULL;
  delete [] word;


}


/* ---------------------------------------------------------------------- */

void FixARTn::min_post_force( int vflag ) {

  cout<<" * IN MIN_POST_FORCES..." << endl;
  post_force( vflag );



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

  cout<< " * Param dt is accessible by update-> dt"<< endl;
  cout<< " * FIRE Params can be change through update->minimize->modify_params( int narg, char **args )"

  /*
  double **force = atoms->f;
  double etot ?
  double forc_conv_thr_qe = update-> etol, ftol, max_eval;
  int nat = atom-> natoms;
  int *ityp = atom->type;
  char *atm ?
  double **tau = atom->x;
  double at  

  artn( force, etot, forc_conv_thr_qe, nat, ityp, atm, tau, at, alat, istep, if_pos, vel, dt, fire_alpha_init, lconv, prefix, tmp_dir )

  */

}













