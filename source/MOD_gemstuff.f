c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
!     University of North Texas
c
c     ####################################################################
c     ##                                                                ##
c     ##  module gemstuff  --  Data for GEM calculations                ##
c     ##                                                                ##
c     ####################################################################
c
c
c
c     coul_coefs  Hermite coefficients for Coulomb interactions
c     exch_coefs  Hermite coefficients for Exchange-Rep. interactions
c     coul_file   File with Hermite coefficients for Coulomb
c     exch_file   File with Hermite coefficients for Exchange-Rep. 
c     atom_crds   File for coords from TINKER
c     gem_input   File with GEM input parameters (nfft, grid size, etc)
c     gem_list    File with atom lists for GEM
c     list_file   File with frames, and molecule/atom mapping for GEM
c 
c
      module gemstuff
      implicit none
      integer first_time_gem
      integer crd_file_lun, free_lun ! fortran tapes for outputs
      integer atom_cnt
      double precision  pbc_box_size(3)
      double precision, allocatable  :: coul_coefs(:)
      double precision, allocatable  :: exch_coefs(:)
      double precision, allocatable  :: cfgem_coul_frc(:,:)
      double precision, allocatable  :: cfgem_exch_frc(:,:)
      double precision, allocatable  :: atom_crds(:,:)
      character*80 user_coul_file 
      character*80 user_exch_file 
      character*80 list_file 
      character*80 aux_file 
      character*80 aux_exch_file 
      character*80 crd_file 
      character*80 gem_input 
      !logical useGEM
c     S:JORGE
      logical usarGEM ! This variable has the same purpose as useGEM
      logical onlyHAL ! This variable allows the calculation of the modified Halgren term
      logical onlyXC  ! This variable allows the calculation of exchange and Coulomb terms
      !logical QMMM 
      real*8 ene_gem_coul
      real*8 ene_gem_exch
      real*8 ene_hal
c     E:JORGE
      save 
      end
