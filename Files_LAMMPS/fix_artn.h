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

/** 
 * @authors 
 *   Matic Poberznic,
 *   Miha Gunde,
 *   Nicolas Salles
 * 
 * @brief Class Fix/ARTn for LAMMPS 
 * 
 * @par Purpose
 * ============
 * Interface for LAMMPS to apply the ARTn algorithm 
 * associate with the FIRE algorithm
 * 
 * @ingroup Interface 
*/
class FixARTn : public Fix {
 public:
  FixARTn(class LAMMPS *, int, char **);
  virtual ~FixARTn();
  int setmask();
  virtual void init();
  void min_setup(int);
  void min_post_force(int);
  void post_run();


 protected:

  // Communication
  void Collect_Arrays( int*, double**, double**, double**, int, double**, double**, double**, int*, int* );
  void Spread_Arrays( int*, double**, double**, double**, int, double**, double**, double** );

  // Resize routine
  void resize_total_system( int );
  void resize_local_system( int );

  // Following and interaction with lammps
  int istep,        //!< Number of time lammps call min_post_force() 
      nword,        //!< Number of arguments  for min->modify_params()
      natoms,       //!< Total number of atoms
      nmax;         //!< Number of atoms and ghost
  char **word;      //!< Array of string for min->modify_params()

  // Engine atomic order
  int *order,              //!< Array with the local order of atoms
      *order_tot;          //!< Array with global order of atoms following the ascending order of id proc

  // Constrains on atomic movement
  int **if_pos;            //!< array of 0/1 if the atoms can move or not

  // Element of type
  char *elt;               //!< Array of Atoms Name (Chemical name)

  // Store the previous force
  int nextblank;           //!< Flag telling if the nexcall is at the same step
  double **f_prev,         //!< 2D array store the previous forces field
         **v_prev,         //!< 2D array store the previous velocity field
         **ftot,           //!< 2D array of global forces over the N procs
         **xtot,           //!< 2D array of global positions over the N procs
         **vtot;           //!< 2D array of global velocities over the N procs

  // energy/force tolerance
  double etol,    //!< Tolerence on energies
         ftol;    //!< Tolerence on forces

  // Fire parameters
  double alpha_init,      //!<  FIRE: \f$ \alpha_0 \f$
         alphashrink;     //!<  FIRE: alpha reduction parameter
  double dt_init,         //!<  FIRE: dt\f$_0\f$
         dtsk,            //!<  FIRE: dt reduction parameter
         dtgrow;          //!<  FIRE: dt increase parameter
  double tmax,            //!<  FIRE: 
         tmin, 
         dtmax, 
         dtmin, 
         dmax;
  int fire_integrator,       //!< FIRE integrator Selector 
      ntimestep_start;       //!< Save the time step when it start
  int delaystep_start_flag;  //!< FIRE control the parameters initialization

  int nsteppos,         //!< Current delaystep (for the relax step)
      nsteppos0;        //!< Defined delaystep (for the relax step)
  double dt_curr;       //!< Current dt
  double alpha;         //!< Current alpha

  class Compute *pe_compute;                          //!< compute for potential energy
  const char *alphab = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";  //!< Should disappear

  // Parallelisation
  int me,                //!< Rank of the procs
      nproc;             //!< Number of Procs
  int *nloc,             //!< Local number of atoms
      oldnloc;           //!< Array to store the local number of atoms
  int *istart,           //!< Array for the resize routines
      *length,           //!< Array for the resize routines
      *nlresize;         //!< Array for the resize routines
  double *tab_comm;      //!< Array for the communication


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
