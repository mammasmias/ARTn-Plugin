#!/bin/bash

#
# A small shell script to patch QE with ARTn; warning this will eventually be done completely differently, only a working version ... 
#

# overwrite modified files 
cp  Modules-modified/*f90 qe-6.6/Modules/

cat PW-src-modified/plugin_ext_forces.f90 \
    PW-src-modified/artn_params_mod.f90 \
    PW-src-modified/center.f90 \
    PW-src-modified/diag.f90 \
    PW-src-modified/lanczos.f90 \
    PW-src-modified/displacement_validation.f90 \
    PW-src-modified/push_init.f90 \
    PW-src-modified/sum_force.f90 \
    PW-src-modified/perpforce.f90 \
    PW-src-modified/move_mode.f90 \
    PW-src-modified/write_report.f90 \
    PW-src-modified/write_struct.f90 \
    PW-src-modified/artn.f90 \
    > qe-6.6/PW/src/plugin_ext_forces.f90

cp  PW-src-modified/dynamics_module.f90  qe-6.6/PW/src/
cp  PW-src-modified/forces.f90  qe-6.6/PW/src/
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
