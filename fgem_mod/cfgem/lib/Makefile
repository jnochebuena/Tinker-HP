# If all three source definitions are blank, libdir.mak
# will automatically include all *.c, *.f, and *.f90 files.
#FC= g95
#FC= gfortran
#FOPTFLAGS= -O3
#FOPTFLAGS= -O0 -g
#LOAD= gfortran
#LOADLIB= 
#VECLIB= veclib.o
FC= ifort
FOPTFLAGS=   -ip -O3 
LOAD= #ifort
LOADLIB=  #-L/share/apps/INTEL/psxe17/compilers_and_libraries_2017/linux/mkl/lib/intel64 -lmkl_vml_def -lmkl_lapack95_lp64 -lmkl_core -lmkl_gnu_thread
#LOADLIB=  -lsvml -L/share/apps/INTEL/psxe17/compilers_and_libraries_2017/linux/mkl/lib/intel64 -lvml -lmkl_lapack -lmkl -lguide -lpthread
VECLIB= 
SHELL=/bin/sh
AR= ar rv
RANLIB=ranlib

MOD_OBJ = \
	user.o auxiliary.o atoms.o sites.o \
	hermite.o site_site_list.o unit_cell.o \
	recip_sum.o recip_grids.o pme_recip.o \
	coulomb_site_site.o exchange_site_site.o ffp_recip.o \
	ewald_recip.o

OTHER_OBJ = \
	utility.o atoms_to_sites_local.o\
        xform_mpoles.o mcmurchie_davidson.o\
        fft_routines.o cfgem_pubfft.o

libcfgem.a:  cfgem.o $(MOD_OBJ) $(OTHER_OBJ) $(VECLIB)
	$(AR) libcfgem.a cfgem.o $(MOD_OBJ) $(OTHER_OBJ)\
		$(VECLIB) $(LOADLIB)
	$(RANLIB) libcfgem.a

cfgem.o:   cfgem.f90 $(MOD_OBJ) $(OTHER_OBJ) $(VECLIB)
	$(FC) -c $(FOPTFLAGS) cfgem.f90

user.o:   user.f90
	$(FC) -c $(FOPTFLAGS) user.f90

auxiliary.o:   user.o auxiliary.f90
	$(FC) -c $(FOPTFLAGS) auxiliary.f90

atoms.o:       atoms.f90
	$(FC) -c $(FOPTFLAGS) atoms.f90

sites.o:   user.o atoms.o auxiliary.o sites.f90
	$(FC) -c $(FOPTFLAGS) sites.f90

hermite.o:   user.o sites.o auxiliary.o hermite.f90
	$(FC) -c $(FOPTFLAGS) hermite.f90

unit_cell.o:   user.o auxiliary.o hermite.o unit_cell.f90
	$(FC) -c $(FOPTFLAGS) unit_cell.f90

site_site_list.o:   sites.o user.o hermite.o unit_cell.o \
		site_site_list.f90
	$(FC) -c $(FOPTFLAGS) site_site_list.f90

coulomb_site_site.o:  user.o sites.o unit_cell.o hermite.o\
		site_site_list.o coulomb_site_site.f90
	$(FC) -c $(FOPTFLAGS) coulomb_site_site.f90

exchange_site_site.o:  user.o sites.o unit_cell.o hermite.o\
		site_site_list.o exchange_site_site.f90
	$(FC) -c $(FOPTFLAGS) exchange_site_site.f90

recip_grids.o:   hermite.o auxiliary.o unit_cell.o \
		recip_grids.f90
	$(FC) -c $(FOPTFLAGS) recip_grids.f90

pme_recip.o:   hermite.o auxiliary.o unit_cell.o \
		recip_grids.o pme_recip.f90
	$(FC) -c $(FOPTFLAGS) pme_recip.f90

ffp_recip.o:   hermite.o auxiliary.o unit_cell.o \
		recip_grids.o ffp_recip.f90
	$(FC) -c $(FOPTFLAGS) ffp_recip.f90

ewald_recip.o:   hermite.o auxiliary.o unit_cell.o \
		recip_grids.o ewald_recip.f90
	$(FC) -c $(FOPTFLAGS) ewald_recip.f90

recip_sum.o:	user.o \
		recip_grids.o pme_recip.o\
		ffp_recip.o ewald_recip.o\
		recip_sum.f90
	$(FC) -c $(FOPTFLAGS) recip_sum.f90

utility.o:   utility.f90
	$(FC) -c $(FOPTFLAGS) utility.f90

atoms_to_sites_local.o:   atoms_to_sites_local.f90
	$(FC) -c $(FOPTFLAGS) atoms_to_sites_local.f90

xform_mpoles.o:   xform_mpoles.f90
	$(FC) -c $(FOPTFLAGS) xform_mpoles.f90

mcmurchie_davidson.o:   mcmurchie_davidson.f90
	$(FC) -c $(FOPTFLAGS) mcmurchie_davidson.f90

fft_routines.o:   fft_routines.f90
	$(FC) -c $(FOPTFLAGS) fft_routines.f90

cfgem_pubfft.o:   cfgem_pubfft.f
	$(FC) -c $(FOPTFLAGS) cfgem_pubfft.f

veclib.o:   veclib.f90
	$(FC) -c $(FOPTFLAGS) veclib.f90
clean:
	-/bin/rm libcfgem.a
	-/bin/rm *.o
	-/bin/rm *.mod
