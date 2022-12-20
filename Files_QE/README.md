# QE-7.0/pARTn Interface 

## How to patch QE:



**Install/update the ARTn-plugin**:
First **configure** QE, then put correct paths of QE and ART in the `environment variables` file, then compile the libartn.a:

```bash
make lib
```

and then run "make patch" to patch QE

```bash
make patch-qe
```

This command recompile automatically PW.

For the QE version older than 7.0 the ARTn patch need more files.

## How to run ARTn with QE:

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

