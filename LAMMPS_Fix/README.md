
# ARTn/LAMMPS interface

ARTn in LAMMPS is defined as a FIX, fix\_artn.h and fix\_artn.cpp.
In the Makefile, i.e. MAKE/Makefile.serial, you need to had the library PATH.
An example:

```
BLAS_LIB := /usr/lib/x86_64-linux-gnu/libopenblas.a -lpthread
FORT_LIB := /usr/lib/x86_64-linux-gnu/libgfortran.so.5

ART_PATH := /home/mammasmias-1/Nico/Project_pARTn/artn-plugin-qe
ART_LIB := $(ART_PATH)/src/libartn.a  $(FORT_LIB) $(BLAS_LIB)
```

And the VARIABLE `ART_LIB` should be added at the moment of the executable creation:
```
$(EXE): main.o $(LMPLIB) $(EXTRA_LINK_DEPENDS)
        $(LINK) $(LINKFLAGS) main.o $(EXTRA_PATH) $(LMPLINK) $(EXTRA_LIB) $(LIB) -o $@ $(ART_LIB)
        $(SIZE) $@
```


