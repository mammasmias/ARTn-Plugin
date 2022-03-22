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
#include "artn.h"

namespace LAMMPS_NS {

class FixARTn : public Fix {
 public:
  FixARTn(class LAMMPS *, int, char **);
  virtual ~FixARTn();
  int setmask();
  virtual void init();
  void min_setup(int);
  void min_post_force(int);
  void post_run();


  // Communication
  void Collect_Arrays( int*, double**, double**, double**, int, double**, double**, double**, int*, int* );
  void Spread_Arrays( int*, double**, double**, double**, int, double**, double**, double** );




 protected:

  // Following and interaction with lammps
  int istep, nword, natoms, nmax;
  char **word;

  // Engine atomic order
  int *order, *order_tot;

  // Constrains on atomic movement
  int **if_pos;

  // Element of type
  char *elt;

  // Store the previous force
  int nextblank;
  double **f_prev, **v_prev, **ftot, **xtot, **vtot;

  // energy/force tolerance
  double etol, ftol;

  // Fire parameters
  double alpha_init, alphashrink;
  double dt_init, dtsk, dtgrow;
  double tmax, tmin, dtmax, dtmin, dmax;
  int fire_integrator, ntimestep_start;
  int delaystep_start_flag;

  int nsteppos;
  double dt_curr;
  double alpha;

  class Compute *pe_compute;        // compute for potential energy
  const char *alphab = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";  // Should disappear

  // Parallelisation
  int me, nproc;
  int *nloc, oldnloc;
  int *istart, *length, *nlresize;
  double *tab_comm;

};

}

// Routines from the library pARTn
/*
extern "C"{
  void artn_( double *const f, double* etot, const int nat, const int *ityp, const char *elt, double *const tau, const int *order, const double *lat, const int *if_pos, int* disp, double *disp_vec, bool* lconv );
  void move_mode_( const int nat, const int* order, double *const f, double *const vel, double* etot, int* nsteppos, double* dt_curr, double* alpha, const double* alpha_init, const double* dt_init, int* disp, double *disp_vec );
  int get_iperp_();
  int get_perp_();
  int get_relx_();
}*/



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