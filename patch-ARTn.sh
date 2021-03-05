#!/bin/bash

#
# A shell script to patch QE with ARTn subroutines ;
# currently all the subroutines are copied into "PW/src/plugin_ext_forces.f90" 
#

# overwrite modified files
cat PW-src-modified/plugin_ext_forces.f90 \
    src/artn_params_mod.f90 \
    src/center.f90 \
    src/diag.f90 \
    src/lanczos.f90 \
    src/displacement_validation.f90 \
    src/pbc.f90  \
    src/invmat3x3.f90  \
    src/push_init.f90 \
    src/sum_force.f90 \
    src/perpforce.f90 \
    src/move_mode.f90 \
    src/write_report.f90 \
    src/write_struct.f90 \
    src/artn.f90 \
    > qe-6.6/PW/src/plugin_ext_forces.f90

# configure and compile q-e
(
    cd qe-6.6/
    if [ -f bin/pw.x ]; then
      echo "q-e already configured"
      make pw
    else
      #./configure
      make pw
    fi
)
