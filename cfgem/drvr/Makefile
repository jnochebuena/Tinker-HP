# If all three source definitions are blank, libdir.mak
# will automatically include all *.c, *.f, and *.f90 files.
#FC= g95
#FC= gfortran
#FOPTFLAGS= -O3
##FOPTFLAGS= -O0 -g
#LOAD= gfortran
#LOADLIB= 
CFGEM_DIR=../lib
FC= ifort
#FOPTFLAGS=  -w95 -ip -O3 
LOAD= ifort
#LOADLIB=  -lsvml -L/opt/intel/mkl721/lib/32 -lvml -lmkl_lapack -lmkl -lguide -lpthread
VECLIB= 
SHELL=/bin/sh

cfgem_drvr:  cfgem_drvr.o
	$(LOAD) -o cfgem_drvr cfgem_drvr.o -L$(CFGEM_DIR) -lcfgem

cfgem_drvr.o:  cfgem_drvr.f90
	$(FC) -c $(FOPTFLAGS) -module $(CFGEM_DIR) cfgem_drvr.f90

clean:
	-/bin/rm cfgem_drvr
	-/bin/rm cfgem_drvr.o
