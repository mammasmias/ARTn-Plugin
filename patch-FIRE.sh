#!/bin/bash

#
# A shell script to add the FIRE minimization algorithm to QE; 
# should be run before attempting to add the ARTn plugin to QE 
# 

# overwrite modified files 
cp  Modules-modified/*f90 qe-6.6/Modules/

cp  PW-src-modified/dynamics_module.f90  qe-6.6/PW/src/
cp  PW-src-modified/move_ions.f90  qe-6.6/PW/src/
cp  PW-src-modified/input.f90  qe-6.6/PW/src/
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
