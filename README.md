# plugin-ARTn

This is a working repo for the current version of the plugin-ARTn; currently it can be used with Quantum ESPRESSO and LAMMPS   

## Contains:

- `Modules-modified/`: Modification on QE's Modules (see Modules-modified/README.md for a description )  ***(Has to disappear)***

- `PW-src-modified/`: Routines modified in PW/src/ of QE mainly to add  ***(Has to disappear)***


- `examples/`: Contains `Al-vacancy.d` and `H2+H.d` 
- `src/`: ARTn plugin subroutines 
- `patch-ARTn.sh`: Patch QE with the ARTn plugin 
- `patch-FIRE.sh`: Patch QE with the FIRE minimization algorithm ***(Has to disappear)***
- `README.md`: The file you are reading
- `LAMMPS_Fix/`: Contains the fix of lammps to interface LAMMPS/ARTn



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



### Use fix_artn

To activate fix ARTn is like all other fix in lammps:

```bash
fix ID group-ID style args value
```

with `style = artn`. For the moment we only test `group-ID = all`. It is possible to custom the FIRE parameters you want to use with the fix ARTn. For each parameter  you give the `name` following by the `value`. The parameters can be:

-  `alpha `
- `alphashrink` 
- `dtshrink`  
- `dmax `
- `tmax `
- `tmin` 

To see the meaning of these parameters refere to the min_fire web page of LAMMPS.



## How to use plugin-ARTn

### Philosophy

Plugin-ARTn is linked with Energy/Forces calculation Engine through the minimization algorithm FIRE. The engine should have this algorithm.  The idea is to launch the engine for a FIRE minimization and the activation of plugin-ARTn bias the minimization to apply the ARTn method.

### Input and Parameters

Once the Engine is compiled with the pARTn library the ARTn input is automatically red at the first moment of the engine minimization step. The ARTn input calls `artn.in`  which allows to change all the  ARTn's parameters. It should be located in the working directory of the calculation. 

Depending of the engine the works units changes and it is to the user to be coherent between the parameters he gives and the units of the engine he uses. This warning is mainly for the **the saddle point convergence**. The user has one keywords to specify the engine he uses (`engine_units`). The list of ARTn's parameters are:

###### General features:

**The values gives by the user through the input file should be in engine units**

- `lrelax`: Values `.true./.false`, default is `false.`.
Flag if we want to relax to adjacent minima from the saddle point.  
- `lpush_final`: Values `.true./.false.`, default is `.true.`.  
Flag to push to adjacent minimum along eigenvector. Flag to push to the second minimum.
- `lrestart`: Values `.true./.false.`, default is `.false.`.
Flag for restarting a ARTn calculation.  
- `ninit`: Value integer, by default is `3`. Number of initial pushes before lanczos start.
- `neigen`: Value integer, by default is `1`. Number of steps made with eigenvector before perpendicular relax.
- `npush`: Value integer, by default is `3`. Maximum number of relaxation perpendicular to the move direction 
- `lanc_mat_size`: Value integer, by default is `16`. Maximum number of Lanczos iterations
- `nsmooth`: Value integer, by default is `1`. Number of smoothing steps from push to eigenvector.
- `struc_format_out`: Value character, default is `"xsf"`. Output structure format. Value accepted `"xyz"` .
Engine specific flag:
- `engine_units`: Value character, default is `qe`. For LAMMPS it is needed to specify `lammps/<units>` where `<units>` correspond to the units keywords in LAMMPS input: (`metal`, `charge`, ...)
###### The push mode:

**The values gives by the user through the input file should be in engine units**

- `push_mode`: Value character, by default is `all`. Type of initial push (`all` , `list` or `rad`)

  - `all`: The initial push is on all the atom in the box with random direction.
  - `list`: The initial push is only on atoms define in the list define by the parameter `push_ids`
  - `rad`: The initial push is on the atoms define by the parameter `push_ids` and the atoms arounds separated by the distance define by `dist_thr` parameters.

- `push_ids`: Value integer, by default is empty. Define the list of atom's id on which we define an initial push for `push_mode = list or rad`. Each atom's id are separated by a coma:

  `push_ids = 23, 201, 35`

- `dist_thr`: Value is real, by default is `0 bohr`. Distance Threshold between the atoms in `push_ids` parameters and the environment.

- `add_const`: Is an array of real, by default empty. Contains the constrain on the initial push if the user want to push in specific region of the space. The constrain contains 4 real value, 3 for the direction and 1 for the solid angle around this direction. An example for the atom 23 on which ask to go in directon (-x,0,z) with solid angle of 30 degrees: 

  `add_const(:,23)=-1.0, 0.0, 1.0, 30`

###### The saddle point convergence:

**The values gives by the user through the input file should be in engine units**

- `init_forc_thr`: Value is real, by default is `1e-2 Ry/bohr` . Initial force convergence criteria. Used for the perpendicular relax before the saddle point convergence
- `final_forc_thr`: Value is real, by default is `1e-3 Ry/bohr`. Final force convergence criteria. Used for the perpendicular relax close to the saddle point.
- `fpara_thr`: Value is real, by default is `5e-3 Ry/bohr`. Initial force convergence criteria. Used for the parallel relaxation.
- `eigval_thr`: Is a real value, by default is `-0.01 Ry/bohr^2` . Threshold for the Hessian eigen value obtain by Lanczos algorithm to start to converge to. The eigen value relative to the saddle point should be negative.
- `frelax_ene_thr`: Is a real value, by default is `-0.01 Ry`. Energy threshold at the saddle point to start relaxation to adjacent minima.
- `push_step_size`: Is a real value, by default is `0.3 bohr`. Step size of the inital push (note: the step size is limited by the engine) 
- `dlanc`: Is a real value, by default is `1e-2 bohr`. Step size in the lanczos algorithm. Should not be (`lanc_step_size`)
- `eigen_step_size`:  Is a real value, by default is `0.2 bohr`. Step size for a step with the lanczos eigenvector (note: the step size is limited by the engine).

###### Ouput:

- `verbose`: Value is integer, by default is `0`. Level `0`  print in output file at each ARTn step without flag information, at `1`  it will add the information flag and at `2`  will print at each step: define push, push and perprelax.  







## The output

- explain the output file... More precisely explain the behavior by default tat pARTn does. 

  I mean the saddle point research first, after the minimization to the first minimum and after the minimization to the second minimum. 

## TODO

- nsteppos in ARTn doesn't have the same meaning for QE and LAMMPS in FIRE algo

- Work on verbose debug mode :ok:

- **Interface LAMMPS**: It's work but still have some noise coming from the velocity. Has to be solve but actually I (Nico) do a break...

- Do pARTn output as ARTn output :ok:

- Adapt the Units output with the Engine input :ok:

- Verify the parameter NAME :ok:

- Add the output filename custom :ok:

- `nperp` parameter it is deactivated when you it converge on the saddle point. Should be Activated when the system return in Basin. :ok:

- **RESTART**: Fast Restart procedure for lammps and binary - Write it only at the end of ARTn-step. :ok:

- **WARNING**: Create a error/warning log file to write all the step does not follow the normal behavior of ARTn.
	- **Transition INIT/PERP**: If the initial push is not enought, the perp-relax is not activated. So the `iinit` is incremented and can reach the lanczos step without never do perp-relax
	- **Transition Saddle/Relax**: If the `eigen_step_size` is too small ARTn can be blocked in PUSH_OVER mode.
	- **Kill the simulation**: We should be able to kill the run when the configuration goes banana!!
  
