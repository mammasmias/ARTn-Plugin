#!/bin/bash

#QE_PATH should be define in environment_variables
source ../../environment_variables

$PARA_PREFIX ${QE_PATH}/bin/pw.x -partn < relax.Al-vacancy.in | tee relax.Al-vacancy.out
