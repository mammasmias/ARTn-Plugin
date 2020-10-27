#!/bin/bash

#
# A small shell script to patch QE with ARTn; warning this will eventually be done completely differently, only a working version ... 
#

# overwrite modified files 
cp  Modules-modified/*f90 qe-6.6/Modules/
cp  PW-src-modified/*f90  qe-6.6/PW/src/

# configure and compile q-e 
(
    cd qe-6.6/
    if [ -f bin/pw.x ]; then
      echo "q-e already configured" 
      make pw
    else
      ./configure
      make pw 
    fi 
)
