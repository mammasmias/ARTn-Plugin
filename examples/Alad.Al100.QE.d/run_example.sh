#!/bin/bash
#QE_PATH should be defined in environment_variables
source ../../../environment_variables
cat > artn.in << EOF 
&ARTN_PARAMETERS
  ! parameters for the push
  push_mode = 'list'
  ! define which atoms are to be pushed and the constraints ...
  push_ids = 1, 14
  add_const(:,1) = 1.0, 1.0, -1.0, 0.
  add_const(:,14) = 1.0, 1.0,  1.0, 0.
  ! eigenpush parms
  eigval_thr = -0.005 
  lpush_final = .true.
/
EOF

mpirun -np 16 ${QE_PATH}/bin/pw.x -partn  < artn.Al-hollow.Al100-5x5-6l.in > artn.Al-hollow.Al100-5x5-6l.out  

cp artn.out artn-exchange.out

cat > artn.in << EOF 
&ARTN_PARAMETERS
  ! parameters for the push
  push_mode = 'list'
  ! define which atoms are to be pushed and the constraints ...
  push_ids = 1
  add_const(:,1) = 1.0, 1.0, -1.0, 0.
  ! eigenpush parms
  lpush_final = .true.
/
EOF

mpirun -np 16 ${QE_PATH}/bin/pw.x -partn  < artn.Al-hollow.Al100-5x5-6l.in > artn.Al-hollow.Al100-5x5-6l.out  

cp artn.out artn-hopping.out

