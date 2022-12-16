## LAMMPS/pARTn Interface






### LAMMPS version after June 2022

#### Installation/Compilation

For the version after June 2022, LAMMPS include the a `Plugin` Class which allows to link LAMMPS with a dynamical library without to recompile at each time.

- **In LAMMPS folder**:
  So first step is to compile LAMMPS in *"shared library"* mode in mpi or serial.

```bash
$ make mode=shared mpi
```

- **In the plugin-ARTn repository**:
  Put the correct path in variable `LAMMPS_PATH` in file `environment_variables` as well as the fortran compiler use to compile library libartn.a in variable  `F90`.
  Therefore in the variable `CXX` put the sample compiler you used to compile LAMMPS.
  Then compile ARTn with the command:

```bash
$ make sharelib
```

At the end of the compilation the file `libartn.so` must appear in the folder `artn-plugin/`.


#### Use fix/artn

To be able to use the `Fix/ARTn` the plugin ARTn has to be loaded.
To load the library `libartn.so` use the command:

```bash
plugin  load  /Path-to-artn-plugin/libartn.so
```

Then you can activate the `Fix/ARTn` like all other fix in lammps:

```bash
fix ID group-ID style args value
```

with `style = artn`. For the moment we only test `group-ID = all`. It is possible to custom the FIRE parameters you want to use with the fix ARTn. For each parameter  you give the `name` following by the `value`. The parameters can be:

-  `alpha `
-  `alphashrink` 
-  `dtshrink`  
-  `dmax `
-  `tmax `
-  `tmin` 

To see the meaning of these parameters refere to the min_fire web page of LAMMPS.