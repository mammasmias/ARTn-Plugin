#!/bin/bash

#
# A shell script to add the FIRE minimization algorithm to QE; 
# should be run before attempting to add the ARTn plugin to QE 
# 

if[ -z ${QE_PATH+x} ] ; 
then 
  echo " You have to declare the environmental_variable QE_PATH"
  stop
fi

# overwrite modified files 
#cp  Modules-modified/*f90 qe-6.6/Modules/
cp  Modules-modified/*f90 $QE_PATH/Modules/

#cp  PW-src-modified/dynamics_module.f90  qe-6.6/PW/src/
#cp  PW-src-modified/move_ions.f90  qe-6.6/PW/src/
#cp  PW-src-modified/input.f90  qe-6.6/PW/src/
cp  PW-src-modified/dynamics_module.f90  $QE_PATH/PW/src/
cp  PW-src-modified/move_ions.f90  $QE_PATH/PW/src/
cp  PW-src-modified/input.f90  $QE_PATH/PW/src/
# configure and compile q-e 
# configure and compile q-e 
(
    #cd qe-6.6/
    cd $QE_PATH
    if [ -f bin/pw.x ]; then
      echo "q-e already configured" 
      make pw
    else
      ./configure
      make pw 
    fi 
)
