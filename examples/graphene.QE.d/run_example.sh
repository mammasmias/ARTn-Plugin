#!/bin/bash

#QE_PATH should be define in environment_variables
source ../../environment_variables

${PARA_PREFIX} ${QE_PATH}/bin/pw.x -partn < relax.graphene-3x2.C-vac.in | tee relax.graphene-3x2.C-vac.out
