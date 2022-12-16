## How to use plugin-ARTn

### Philosophy

Plugin-ARTn is linked with Energy/Forces calculation engine through the minimization algorithm FIRE. The idea is to launch the engine for a FIRE minimization, and bias the minimization with the plugin-ARTn, to hijack FIRE and perform the ARTn method instead.

### How to use

ARTn can be used in different ways, with some optional variables.

- Research from the local minimum...
- Saddle refine...

### Input and Parameters

Once the Engine is compiled with the pARTn library, the ARTn input file `artn.in` is automatically read at the first step of the FIRE minimization. All parameters related to the ARTn calculation are defined in the file `artn.in`, which should be located in the working directory of the calculation. 

Depending of the engine, the working units change, and it is up to the user to be coherent between the units of input parameters, and the units of the engine. 
This warning is mainly for the **saddle point convergence**. The units are specified by the parameter `engine_units`. 
<!-- The list of ARTn's parameters are: !-->
**The values given in the input file should be in engine units**

###### General features:

- `lpush_final`: Values `.true./.false.`, default is `.true.`. 
  Flag to push to adjacent minimum along eigenvector. Flag to push to the second minimum.

- `lrestart`: Values `.true./.false.`, default is `.false.`.
  Flag for restarting a ARTn calculation.  initialize the ARTn parameters and configuration from the file `artn.restart`.

- `lmove_nextmin`: Value .true./.false., default is .false. 

  Flag to leave ARTn with the new minimum found. 

- `ninit`: Value integer, by default is `3`. Number of initial pushes before lanczos start.

- `neigen`: Value integer, by default is `1`. Number of steps made with eigenvector before perpendicular relax.

- `nperp`: Value integer, by default is `-1`. Maximum number of relaxation perpendicular to the move direction after an `init` or `eigen` push.

- `lnperp_limitation`: Values `.true./.false.`, default is `.true.`. 

  this option allows to constrain the number of perpendicular relaxation during the convergence to the saddle point, out of the basin. The limitation is incremental starting by 8, 12, 16, -1 (infinite). These values are stored in arrays `nperp_limitation(:)` where the first value is `nperp` in the basin. These list can be customized in input giving the values of the array as: `nperp_limitation = [...custom values]` 

- `lanc_mat_size`: Value integer, by default is `16`. Maximum number of Lanczos iterations

- `nsmooth`: Value integer, by default is `0`. Number of smoothing steps from push to eigenvector.

- `struc_format_out`: Value character, default is `"xsf"`. Output structure format. Value accepted `"xyz"` .
  Engine specific flag:

- `engine_units`: Value character, default is `qe`. For LAMMPS it is needed to specify `lammps/<units>` where `<units>` correspond to the units keywords in LAMMPS input: (`metal`, `charge`, ...)

###### The push mode:

- `push_mode`: Value character, by default is `all`. Type of initial push (`all` , `list` or `rad`)

  - `all`: The initial push is generated on all the atoms in the box with random direction.
  - `list`: The initial push is only on atoms defined in the list, defined by the array `push_ids`
  - `rad`: The initial push is on the atoms defined by the array `push_ids` and the atoms within the radial distance defined by `dist_thr` parameter.

- `push_ids`: Value integer, by default is empty. Defines the list of atom ids on which the initial push is generated, for `push_mode = list or rad`. The ids are separated by a coma, e.g.:

  `push_ids = 23, 201, 35`

- `dist_thr`: Value is real, by default is `0.0`. Radial distance threshold for `push_mode=rad`.

- `add_const`: Array of real, by default empty. Contains the constrain on the initial push.
The constrain contains 4 real values, 3 for the directions, and 1 for the solid angle around this direction. For example, to generate a push within solid angle of 30 degrees in the direction (-1, 0, 1), for the atom index 23: 

  `add_const(:,23)=-1.0, 0.0, 1.0, 30`

###### The saddle point convergence:

- `forc_thr`: Value is real, by default is `1e-3 Ry/bohr`. Force criteria for convergence to saddle.
- `fpara_thr`: Value is real, by default is `5e-3 Ry/bohr`. Initial force convergence criteria. Used for the parallel relaxation. UNUSED!! Remove from code.
- `eigval_thr`: Is a real value, by default is `-0.01 Ry/bohr^2` . Threshold for the Hessian eigenvalue obtained by Lanczos algorithm, which decides when to start following the corresponding eigenvector. The eigenvalue relative to the saddle point should be negative.
- `frelax_ene_thr`: Is a real value, by default is `-0.01 Ry`. Energy threshold at the saddle point to start relaxation to adjacent minima.
- `push_step_size`: Is a real value, by default is `0.3 bohr`. Step size of the inital push (note: the step size is limited by the engine) 
- `lanczos_disp`: Is a real value, by default is `1e-2 bohr`. Step size in the lanczos algorithm.
- `lanczos_eval_conv_thr`: Is a real value, default is `1e-2`. Threshold for convergence of eignevalue within each Lanczos evaluation, in relative units.
- `lanczos_max_size`: Integer, default is 16. Maximal number of iteration of each Lanczos evaluation.
- `lanczos_min_size`: Integer, default is 0. Minimal number of iteration of each Lanczos evaluation.
- `eigen_step_size`:  Is a real value, by default is `0.2 bohr`. Maximal step size for a step with the lanczos eigenvector (note: the step size is limited by the engine).

###### Output:

- `verbose`: Value is integer, by default is `0`. Level `0`  print in output file at each ARTn step without flag information, at `1`  it will add the information flag and at `2`  will print at each step: define push, push and perprelax.  

## The output

Various files can be found in output. 

- The details of ARTn research is given in file named `artn.out`. This file follows the engine units defined thank to the variable `engine_units`.  It gives the system parameters evolution used by ARTn to find the saddle point as well as the the balance energies between the locals minimum and saddle point.
- At the beginnning ARTn write the initial configuration in file `initp.*` where `*` is the format indicated by the variable `struc_format_out`.
- During the ARTn path toward the saddle point, the currently followed eigenvector is stored in the file `lastest_engenvec.*`.
- For each convergence reached, saddle point and local minima, the configurations are stored in files `sad####.*`, and `min####.*`
