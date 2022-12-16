## LAMMPS/pARTn Interface



### LAMMPS version before June 2022

We tested this interface only with gnu compiler. 

**In the plugin-ARTn repository**: First put the correct path in variable `LAMMPS_PATH` in file `environment_variables` as well as the fortran compiler use to compile library libartn.a in variable  `F90`. Afterwards compile the library `libartn.a`:

```bash
make lib
```

patch lammps, means copy the files in `Files_LAMMPS/` to `LAMMPS_PATH/src` doing:

```bash
make patch-lammps
```

**In LAMMPS folder**: The plugin-ARTn library has to be called during the LAMMPS compilation. In the Makefile you want to use to compile it, i.e. `LAMMPS_PATH/src/MAKE/` the following lines has to be added. One variable to inform on library OPENBLAS, one to inform on the fortran library because the C++/Fortran interface and one last  to imform on the path of ART repository.
An example:

```makefile
BLAS_LIB := /usr/lib/x86_64-linux-gnu/libopenblas.a -lpthread
FORT_LIB := /usr/lib/x86_64-linux-gnu/libgfortran.so.5

ART_PATH := /home/src/artn-plugin-qe
ART_LIB := $(ART_PATH)/src/libartn.a  $(FORT_LIB) $(BLAS_LIB)
```

And the VARIABLE `ART_LIB` should be added at the moment of the executable creation, for example at the end of the line as it is shown:

```makefile
$(EXE): main.o $(LMPLIB) $(EXTRA_LINK_DEPENDS)
	$(LINK) $(LINKFLAGS) main.o $(EXTRA_PATH) $(LMPLINK) $(EXTRA_LIB) $(LIB) -o $@ $(ART_LIB)
	$(SIZE) $@
```

Now you can compile LAMMPS using the normal command with the good name of the Makefile (serial/mpi)

```bash
make yourmakefile
```



#### Use fix/artn

To activate fix ARTn is like all other fix in lammps:

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