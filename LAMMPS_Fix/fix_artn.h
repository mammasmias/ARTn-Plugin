/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(artn,FixARTn)

#else

#ifndef LMP_FIX_ARTN_H
#define LMP_FIX_ARTN_H

#include "fix.h"

namespace LAMMPS_NS {

class FixARTn : public Fix {
 public:
  FixARTn(class LAMMPS *, int, char **);
  virtual ~FixARTn();
  int setmask();
  virtual void init();
  void setup(int);
  void min_setup(int);
  virtual void post_force(int);
  //void post_force_respa(int, int, int);
  void min_post_force(int);

/*  extern "C"{
    void artn_( double **force, double etot, int nat, int *ityp, char *elt, double **tau, double lat[3][3], int *if_pos, char *move, bool lconv );
    void move_mode_( int nat, double **force, double **vel, double etot, int nsteppos, double dt_curr, double alpha, double alpha_init, double dt_init, char *cmode );
  }
*/

  //double compute_vector(int);
  //double memory_usage();

 protected:
  int istep;
  double alpha_init, dt_init;

  class Compute *pe_compute;        // compute for potential energy
  const char *alphab= "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

  constexpr unsigned int str2int(const char* str, int h = 0) {
    return !str[h] ? 5381 : (str2int(str, h+1) * 33) ^ str[h];
  }


};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Region ID for fix setforce does not exist

Self-explanatory.

E: Variable name for fix setforce does not exist

Self-explanatory.

E: Variable for fix setforce is invalid style

Only equal-style variables can be used.

E: Cannot use non-zero forces in an energy minimization

Fix setforce cannot be used in this manner.  Use fix addforce
instead.

*/
