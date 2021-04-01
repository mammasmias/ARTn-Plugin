#!/bin/bash

#
# A shell script to patch QE with ARTn subroutines ;
# currently all the subroutines are copied into "PW/src/plugin_ext_forces.f90" 
#

# overwrite modified files
cp qe-6.6/PW/src/plugin_ext_forces.f90 .
cat PW-src-modified/plugin_ext_forces.f90 > qe-6.6/PW/src/plugin_ext_forces.f90

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
