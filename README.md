# ARTn-plugin-QE

This is a working repo for the current version of the ARTn plugin; currently it can only be used with Quantum ESPRESSO    

## Contains:

- `Modules-modified/`: Modification on QE's Modules (see Modules-modified/README.md for a description )  

- `PW-src-modified/`: Routines modified in PW/src/ of QE mainly to add  


- `examples/`: Contains `Al-vacancy.d` and `H2+H.d` 

- `qe-6.6/`: Quantum ESPRESSO software

- `src/`: ARTn plugin subroutines 

- `patch-ARTn.sh`: Patch QE with the ARTn plugin 

- `patch-FIRE.sh`: Patch QE with the FIRE minimization algorithm 

- `README.md`: The file you are reading 

## How to patch QE:

**First installation:**  Execute  the `patch-FIRE.sh` script (adds the FIRE minimization option to QE) 

**Install/update the ARTn-plugin**:
First configure QE, then put correct paths of QE and ART in the environment variables file, then run "make" to compile the libartn.a, and then run "make patch" to patch QE and make pw.

## How to run ARTn:

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



