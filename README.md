# plugin-ARTn

This is a working repo for the current version of the plugin-ARTn; currently it can be used with Quantum ESPRESSO and LAMMPS   

## Contains:

- `Modules-modified/`: Modification on QE's Modules (see Modules-modified/README.md for a description )  

- `PW-src-modified/`: Routines modified in PW/src/ of QE mainly to add  


- `examples/`: Contains `Al-vacancy.d` and `H2+H.d` 
- `qe-6.6/`: Quantum ESPRESSO software
- `src/`: ARTn plugin subroutines 
- `patch-ARTn.sh`: Patch QE with the ARTn plugin 
- `patch-FIRE.sh`: Patch QE with the FIRE minimization algorithm 
- `README.md`: The file you are reading
- `LAMMPS_Fix/`: Contains the fix of lammps to interface LAMMPS/ARTn

## TODO

- nsteppos in ARTn doesn't have the same meaning for QE and LAMMPS
- Work on the parallelization:
  - one ARTn run with many proc
  - many ARTn run with one proc
  - many ARTn run with many proc
- Generalize the unit conversion -> `MODULE units`

## QE/pARTn Interface 

### How to patch QE:

**First installation:**  Execute  the `patch-FIRE.sh` script (adds the FIRE minimization option to QE) 

**Install/update the ARTn-plugin**:
First configure QE, then put correct paths of QE and ART in the environment variables file, then run "make" to compile the libartn.a, and then run "make patch" to patch QE and make pw.

### How to run ARTn with QE:

For example calculations please see the `examples` directory. A brief
description of the minimum requirements is provided here.

In order to run an ARTn calculation you have to add the following
lines to the QE input file `pwscf.in`:

```fortran
&CONTROL
 calculation = 'relax' 
/

&SYSTEM
  nosym = .true. 
/

&IONS
  ion_dynamics = 'fire' 
/
```

Finally Quantum ESPRESSO must be launched with the flag -partn as follow:

```bash
./pw.x -partn -inp input_pw.txt
```

 



## LAMMPS/pARTn Interface

### Installation/Compilation

For the moment the pARTn library has to be compiled with **gfortran** only.

ARTn, in LAMMPS, is defined as a FIX, `fix_artn.h` and `fix_artn.cpp`. So  you copy and paste these two files in the `LAMMPS/src/`.

In the Makefile, i.e. `LAMMPS/src/MAKE/Makefile.serial`, you need to had the library PATH. For pARTn it needs the openblas library with pthread library and the gfortran library for the C++/fortran interface. Of course the ARTn library is built at ARTn compilation and placed in the `src/` folder.
An example:

```makefile
BLAS_LIB := /usr/lib/x86_64-linux-gnu/libopenblas.a -lpthread
FORT_LIB := /usr/lib/x86_64-linux-gnu/libgfortran.so.5

ART_PATH := /home/src/artn-plugin-qe
ART_LIB := $(ART_PATH)/src/libartn.a  $(FORT_LIB) $(BLAS_LIB)
```

And the VARIABLE `ART_LIB` should be added at the moment of the executable creation:

```makefile
$(EXE): main.o $(LMPLIB) $(EXTRA_LINK_DEPENDS)
	$(LINK) $(LINKFLAGS) main.o $(EXTRA_PATH) $(LMPLINK) $(EXTRA_LIB) $(LIB) -o $@ $(ART_LIB)
	$(SIZE) $@
```

Now you can compile LAMMPS using the normal command with the good name of the Makefile (serial/mpi)

```bash
make serial
```



## How to use plugin-ARTn

Plugin-ARTn enters in the FIRE algorithm minimization loop.
The parameters specific to ARTn algorithm must be specified in the `artn.in` file, located in the working directory of the calculation.
The list of parameters are:

###### General features:

- `lrelax`: Values `.true./.false`, default is `false.`.
Flag if we want to relax to adjacent minima from the saddle point.  
- `lpush_final`: Values `.true./.false.`, default is `.true.`.  
Flag to push to adjacent minimum along eigenvector.
- `lrestart`: Values `.true./.false.`, default is `.false.`.
Flag for restarting a ARTn calculation.  
- `npush`: Value integer, by default is `3`. Number of initial pushes before lanczos start.
- `neigen`: Value integer, by default is `1`. Number of steps made with eigenvector before perpendicular relax.
- `nlanc_init`: Value integer, by default is `16`. Maximum number of lanczos iterations
- `nsmooth`: Value integer, by default is `1`. Number of smoothing steps from push to eigenvector.
- `struc_format_out`: Value character, default is `"xsf"`. Output structure format. Value accepted `"xyz"` .
Engine specific flag:
- `engine_units`: Value character, default is `qe`. For LAMMPS it is needed to specify `lammps/<units>` where `<units>` correspond to the units keywords in LAMMPS input: (`metal`, `charge`, ...)
###### The push mode:

- `push_mode`: Value character, by default is `all`. Type of initial push (`all` , `list` or `rad`)

  - `all`: The initial push is on all the atom in the box with random direction.
  - `list`: The initial push is only on atoms define in the list define by the parameter `push_ids`
  - `rad`: The initial push is on the atoms define by the parameter `push_ids` and the atoms arounds separated by the distance define by `dist_thr` parameters.

- `push_ids`: Value integer, by default is empty. Define the list of atom's id on which we define an initial push for `push_mode = list or rad`. Each atom's id are separated by a coma:

  `push_ids = 23, 201, 35`

- `dist_thr`: Value is real, by default is `0`. The unit is in Angstrom. Distance Threshold between the atoms in `push_ids` parameters and the environment.

- `add_const`: Is an array of real, by default empty. Contains the constrain on the initial push if the user want to push in specific region of the space. The constrain contains 4 real value, 3 for the direction and 1 for the solid angle around this direction. An example for the atom 23 on which ask to go in directon (-x,0,z) with solid angle of 30 degrees: 

  `add_const(:,23)=-1.0, 0.0, 1.0, 30`

###### The saddle point convergence:

- `convcrit_init`: Value is real, by default is `1e-2`. The unit is eV/Angstrom. Initial force convergence criteria. Used for the perpendicular relax before the saddle point convergence
- `convcrit_final`: Value is real, by default is `1e-3`. The unit is eV/Angstrom. Initial force convergence criteria. Used for the perpendicular relax close to the saddle point.
- `fpara_convcrit`: Value is real, by default is `5e-3`. The unit is eV/Angstrom. Initial force convergence criteria. Used for the parallel relaxation.
- `eigval_thr`: Is a real value, by default is `-0.01` Ry/Angs^2. Threshold for the Hessian eigen value obtain by Lanczos algorithm to start to converge to. The eigen value relative to the saddle point should be negative.
- `relax_thr`: Is a real value, by default is `0.01`. Energy Threshold at the saddle point to start relaxation to adjacent minima.
- `push_step_size`: Is a real value, by default is `0.3`. The unit is in Angstrom. Step size of the inital push (note: the step size is limited by the engine) 
- `dlanc`: Is a real value, by default is `1e-2`. The unit is in Angstrom. Step size in the lanczos algorithm.
- `eigen_step_size`:  Is a real value, by default is `0.2`. The unit is in Angstrom. Step size for a step with the lanczos eigenvector (note: the step size is limited by the engine).



## The output

- explain the output file... More precisely explain the behavior by default tat pARTn does. 

  I mean the saddle point research first, after the minimization to the first minimum and after the minimization to the second minimum. 
