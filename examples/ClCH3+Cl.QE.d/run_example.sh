#!/bin/bash

#QE_PATH should be define in environment_variables
source ../../environment_variables

mpirun -np 4 ${QE_PATH}/bin/pw.x -partn < relax.ClCH3+Cl.in | tee relax.ClCH3+Cl.out
