#=======================================
# WANNIER90 	
#=======================================

TOPDIR = /home/matic/programming/ARTn-pluggin/q-e.fire

F90=mpif90
FCOPTS=$(FFLAGS) -cpp 
LDOPTS=-g

LIBS = -L/home/matic/libs/lapac/lapack-3.7.0 -llapack  -lblas   -lblas  
