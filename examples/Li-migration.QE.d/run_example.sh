#!/bin/bash

#QE_PATH should be define in environment_variables
source ../../environment_variables

${PARA_PREFIX} ${QE_PATH}/bin/pw.x -partn < relax.Li-migration.graphite-3x2x2.in | tee relax.Li-migration.graphite-3x2x2.out 
