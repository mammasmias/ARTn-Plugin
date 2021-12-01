#!/bin/bash
QE_PATH=../../qe-6.6/bin

mpirun -np 2 ${QE_PATH}/pw.x -partn < relax.Al-vacancy.in | tee relax.Al-vacancy.out
