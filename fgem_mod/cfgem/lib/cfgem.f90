module cfgem_libmod

  implicit none

  private

  public        :: cfgem_ene_rec, null_cfgem_ene_rec, cfgem_coul, cfgem_exch

  ! This is the record structure used to return energies from the cfgem_*
  ! evaluation subroutines:

  ! mp - mpole
  ! ch - compact hermite
  ! dh - diffuse hermite

  type cfgem_ene_rec
    double precision    :: mp_mp
    double precision    :: mp_ch
    double precision    :: ch_mp
    double precision    :: ch_ch
    double precision    :: mp_dh
    double precision    :: dh_mp
    double precision    :: ch_dh
    double precision    :: dh_ch
    double precision    :: dh_dh
  end type cfgem_ene_rec

  type(cfgem_ene_rec), parameter      :: null_cfgem_ene_rec = &
    cfgem_ene_rec(0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0)

contains

! NOTE - The pbc_box_size argument is optional; if not provided, then pbc
!        box dimensions are taken from the user file.

!******************************************************************************
!
! Subroutine:   cfgem_coul
!
! Description:  Evaluates a described system for gem coulomb energies and
!               forces at a single point in time..
!
! Input Arguments:
!
!   aux_file, list_file, user_coul_file - gem input parameters file names
!   (see samples under test)
!
!   atom_cnt - count of atoms in system
!   atom_crds - coordinates for all atoms (Angstrom)
!   pbc_box_size - size of unit cell (Angstrom).  OPTIONAL, at end of arg list;
!                  if not present, values are taken from user_coul_file.
!   free_lun - a file logical unit number to be used by library temp file
!              operations (ie., you are not using it in the caller).
!   out_lun  - file logical unit number to be used by library for printed
!              output, if necessary.  The library will return actual data in
!              the output arguments, so this is used when something goes wrong.
!              Typical values are 6, or mdout in amber code.
!
! Output Arguments:
!
!   (ALL output arguments are zero'd on input)
!
!   cfgem_tot_coul_ene - total coulomb energy, Hartree
!   cfgem_coul_ene_fld - total coulomb field energy, Hartree
!   cfgem_coul_recip_ene - reciprocal space coulomb energy components, Hartree
!   cfgem_coul_direct_ene - direct space coulomb energy components, Hartree
!   cfgem_coul_adj_recip_ene - reciprocal space excluded interactions adjust
!                              coulomb energy components, Hartree
!   cfgem_coul_self_ene - self coulomb energy components, Hartree
!   cfgem_coul_wigner_ene - Wigner correction coulomb energy components, Hartree
!   cfgem_coul_frc - Array of total coulomb forces, Hartree/Bohr
!
!******************************************************************************

subroutine cfgem_coul(aux_file, list_file, user_coul_file, &
                      atom_cnt, atom_crds, free_lun, out_lun, &
                      cfgem_tot_coul_ene, &
                      cfgem_coul_ene_fld, &
                      cfgem_coul_recip_ene, &
                      cfgem_coul_direct_ene, &
                      cfgem_coul_adj_recip_ene, &
                      cfgem_coul_self_ene, &
                      cfgem_coul_wigner_ene, &
                      cfgem_coul_frc, &
                      pbc_box_size)

   use atoms,only : FH_ATOMS_init, &
                    FH_ATOMS_deallocate

   use user,only : FH_USER_read, &
                   FH_USER_first_check, &
                   FH_USER_second_check, &
                   dump_atom_forces, &
                   do_coulomb, &
                   do_overlap, &
                   calculation_method

   use auxiliary,only : FH_AUX_read, &
                        FH_AUX_get_rec_expon, &
                        FH_AUX_recip_auto_setup, &
                        FH_AUX_deallocate

   use sites,only : FH_SITES_readfile, &
                    FH_SITES_build, &
                    FH_SITES_deallocate

   use hermite, only : FH_hermite_load, &
                       FH_hermite_deallocate, &
                       FH_hermite_density_extent

   use site_site_list, only : FH_SITE_SITE_LIST_init, &
                              FH_SITE_SITE_LIST_deallocate

   use unit_cell, only : FH_UNIT_CELL_setup, &
                         FH_UNIT_CELL_deallocate

   use recip_sum, only : FH_RECIP_setup, &
                         FH_RECIP_deallocate

!! GAC: store GEM Coulomb and Exchange total Energy and Forces in 
!!      TINKER multipole vars
   !use energi,    only:  em
   !use deriv,     only:  dem
                        
   implicit none

! Formal arguments:

   ! atom_crds() in angstrom..
   ! free_lun is used by library as lun for reading input files.
   ! out_lun is lun for any printed output (should be error msgs only)
   ! pbc_box_size() in angstrom; if not present, use user_coul_file value.

   character(len=80), intent(in)        :: aux_file
   character(len=80), intent(in)        :: list_file
   character(len=80), intent(in)        :: user_coul_file
   integer,           intent(in)        :: atom_cnt
   double precision,  intent(in)        :: atom_crds(3, atom_cnt)
   integer,           intent(in)        :: free_lun
   integer,           intent(in)        :: out_lun

   double precision,    intent(out)     :: cfgem_tot_coul_ene ! ALL coulomb

   double precision,    intent(out)     :: cfgem_coul_ene_fld ! (for torque)

   type(cfgem_ene_rec), intent(out)     :: cfgem_coul_recip_ene
   type(cfgem_ene_rec), intent(out)     :: cfgem_coul_direct_ene
   type(cfgem_ene_rec), intent(out)     :: cfgem_coul_adj_recip_ene
   type(cfgem_ene_rec), intent(out)     :: cfgem_coul_self_ene
   type(cfgem_ene_rec), intent(out)     :: cfgem_coul_wigner_ene

   double precision,    intent(inout)     :: cfgem_coul_frc(3, atom_cnt)

   double precision, optional,  intent(in)        :: pbc_box_size(3)

! Local variables:
   
   include "interact_type.fh"

   double precision, allocatable, save  :: atom_crds_bohr(:,:)

   ! NOTE: While the subroutines used in this function return a virial, the
   !       recip space coulomb virial is known to be incorrect; thus we
   !       currently return no virial information..

   double precision                     :: cfgem_coul_virial(3,3)

   integer              :: i, j, k
   integer              :: alloc_failed

   ! Go ahead and convert coordinates from Anstrom to Bohr, in locally
   ! allocated array; don't want to modify input!

   call FH_ATOMS_init(atom_cnt, atom_crds, out_lun)

   call FH_SITE_SITE_LIST_init()

   if (present(pbc_box_size)) then
      call FH_USER_read(user_coul_file, free_lun, out_lun, pbc_box_size)
   else
      call FH_USER_read(user_coul_file, free_lun, out_lun)
   end if

   call FH_USER_first_check(do_overlap, out_lun)

   ! We only enforce that the wrong flag is not set..  Thus, user can skip
   ! inputting either if desired..

   if (do_overlap .ne. 0) then
     write(out_lun)'CFGEM_COUL: do_overlap set in user coulomb file!'
     stop
   end if

   ! Zero all output arguments:

   cfgem_tot_coul_ene = 0.d0                            ! SET
   cfgem_coul_ene_fld = 0.d0                            ! SET
   cfgem_coul_recip_ene = null_cfgem_ene_rec            ! SET
   cfgem_coul_direct_ene = null_cfgem_ene_rec           ! SET
   cfgem_coul_adj_recip_ene = null_cfgem_ene_rec        ! SET
   cfgem_coul_self_ene = null_cfgem_ene_rec             ! SET
   cfgem_coul_wigner_ene = null_cfgem_ene_rec           ! SET
   cfgem_coul_frc(:,:) = 0.d0                           ! SET
   cfgem_coul_virial(:,:) = 0.d0                        ! SET

   call FH_AUX_read(aux_file, free_lun, out_lun)

   call FH_AUX_recip_auto_setup(out_lun) ! setup based on auxiliary info

   call FH_USER_second_check(out_lun) ! check on the recip settings

   call FH_AUX_get_rec_expon(out_lun)

   call FH_SITES_readfile(list_file, free_lun, out_lun)

   call FH_hermite_load(out_lun)

   call FH_MCMUR_DAV_init()

   call FH_SITES_build()

   call FH_hermite_density_extent(out_lun)

   call FH_UNIT_CELL_setup(out_lun)

   call FH_RECIP_setup(out_lun)

   call nonbond_eval_coul()

   !em = em + cfgem_tot_coul_ene
   !dem = dem + cfgem_coul_frc

   ! Deallocate:

   call FH_ATOMS_deallocate()
   call FH_AUX_deallocate()
   call FH_SITES_deallocate()
   call FH_hermite_deallocate()
   call FH_UNIT_CELL_deallocate()
   call FH_RECIP_deallocate()
   call FH_SITE_SITE_LIST_deallocate()

   return

contains

subroutine nonbond_eval_coul()

   use user, only : verbose,do_coulomb,calculation_method, &
                    do_overlap

   use site_site_list, only : FH_SITE_SITE_LIST_eval

   use recip_sum, only : FH_RECIP_eval

   use sites, only : site_frc

   use hermite, only : FH_hermite_local_to_global, &
                       tot_num_mpole_coeffs, &
                       tot_num_sumCH_coeffs, &
                       tot_num_herm_Cprims, &
                       tot_num_herm_Dprims, &
                       Global_multipole_field, &
                       Global_sumCH_field, &
                       Global_Chermite_field, &
                       Global_Dhermite_field

   use coulomb_site_site, only : FH_COULOMB_SITE_SITE_CC, &
                                 FH_COULOMB_SITE_SITE_self_CC, &
                                 FH_COULOMB_SITE_SITE_CD_DD, &
                                 FH_COULOMB_SITE_SITE_self_CD_DD

   use site_site_list, only : &
                            max_unmasked_pairs_num_nghbrs, &
                            unmasked_pairs_list, &
                            unmasked_pairs_ucell_transind, &
                            unmasked_pairs_num_nghbrs, &
                            unmasked_pairs_list_offset, &
                            max_cutoff_direct_num_nghbrs, &
                            cutoff_direct_list, &
                            cutoff_direct_ucell_transind, &
                            cutoff_direct_num_nghbrs, &
                            cutoff_direct_list_offset, &
                            max_correct_recip_num_nghbrs, &
                            correct_recip_list, &
                            correct_recip_ucell_transind, &
                            correct_recip_num_nghbrs, &
                            correct_recip_list_offset

   implicit none

! Formal arguments:

! Local variables:

   double precision :: cutoff_mp_mp_ene, &
                       cutoff_mp_ch_ene, &
                       cutoff_ch_mp_ene, &
                       cutoff_ch_ch_ene, &
                       cutoff_mpole_mpole_ene, &
                       cutoff_mpole_herm_ene, &
                       cutoff_herm_herm_ene, &
                       ! self has all 9 components
                       self_mp_mp_ene, &
                       self_mp_ch_ene, &
                       self_mp_dh_ene, &
                       self_ch_mp_ene, &
                       self_ch_ch_ene, &
                       self_ch_dh_ene, &
                       self_dh_mp_ene, &
                       self_dh_ch_ene, &
                       self_dh_dh_ene, &
                       self_mpole_mpole_ene, &
                       self_mpole_herm_ene, &
                       self_herm_herm_ene, &
                       ! correct_recip has all 9 components
                       correct_recip_mp_mp_ene, &
                       correct_recip_mp_ch_ene, &
                       correct_recip_mp_dh_ene, &
                       correct_recip_ch_mp_ene, &
                       correct_recip_ch_ch_ene, &
                       correct_recip_ch_dh_ene, &
                       correct_recip_dh_mp_ene, &
                       correct_recip_dh_ch_ene, &
                       correct_recip_dh_dh_ene, &
                       correct_recip_mpole_mpole_ene, &
                       correct_recip_mpole_herm_ene, &
                       correct_recip_herm_herm_ene, &
                       ! recip has all 9 components
                       recip_mp_mp_ene, &
                       recip_mp_ch_ene, &
                       recip_mp_dh_ene, &
                       recip_ch_mp_ene, &
                       recip_ch_ch_ene, &
                       recip_ch_dh_ene, &
                       recip_dh_mp_ene, &
                       recip_dh_ch_ene, &
                       recip_dh_dh_ene, &
                       recip_mpole_mpole_ene, &
                       recip_mpole_herm_ene, &
                       recip_herm_herm_ene, &
                       ! wigner has all 9
                       wigner_mp_mp_ene, &
                       wigner_mp_ch_ene, &
                       wigner_mp_dh_ene, &
                       wigner_ch_mp_ene, &
                       wigner_ch_ch_ene, &
                       wigner_ch_dh_ene, &
                       wigner_dh_mp_ene, &
                       wigner_dh_ch_ene, &
                       wigner_dh_dh_ene, &
                       wigner_mpole_mpole_ene, &
                       wigner_mpole_herm_ene, &
                       wigner_herm_herm_ene, &
                       ! totals
                       tot_wigner_ene, &
                       tot_mpole_mpole_ene, &
                       tot_mpole_herm_ene, &
                       tot_herm_herm_ene, &
                       factor

   integer :: correcting_energies, &
              use_cutoff, &
              full_interactions

   include "interact_type.fh"

   cutoff_mp_mp_ene = 0.d0
   cutoff_mp_ch_ene = 0.d0
   cutoff_ch_mp_ene = 0.d0
   cutoff_ch_ch_ene = 0.d0

   recip_mp_mp_ene = 0.d0
   recip_mp_ch_ene = 0.d0
   recip_mp_dh_ene = 0.d0
   recip_ch_mp_ene = 0.d0
   recip_ch_ch_ene = 0.d0
   recip_ch_dh_ene = 0.d0
   recip_dh_mp_ene = 0.d0
   recip_dh_ch_ene = 0.d0
   recip_dh_dh_ene = 0.d0
   recip_mpole_mpole_ene = 0.d0
   recip_mpole_herm_ene = 0.d0
   recip_herm_herm_ene = 0.d0

   wigner_mp_mp_ene = 0.d0
   wigner_mp_ch_ene = 0.d0
   wigner_mp_dh_ene = 0.d0
   wigner_ch_mp_ene = 0.d0
   wigner_ch_ch_ene = 0.d0
   wigner_ch_dh_ene = 0.d0
   wigner_dh_mp_ene = 0.d0
   wigner_dh_ch_ene = 0.d0
   wigner_dh_dh_ene = 0.d0
   wigner_mpole_mpole_ene = 0.d0
   wigner_mpole_herm_ene = 0.d0
   wigner_herm_herm_ene = 0.d0

   self_mp_mp_ene = 0.d0
   self_mp_ch_ene = 0.d0
   self_mp_dh_ene = 0.d0
   self_ch_mp_ene = 0.d0
   self_ch_ch_ene = 0.d0
   self_ch_dh_ene = 0.d0
   self_dh_mp_ene = 0.d0
   self_dh_ch_ene = 0.d0
   self_dh_dh_ene = 0.d0
   self_mpole_mpole_ene = 0.d0
   self_mpole_herm_ene = 0.d0
   self_herm_herm_ene = 0.d0

   correct_recip_mp_mp_ene = 0.d0
   correct_recip_mp_ch_ene = 0.d0
   correct_recip_mp_dh_ene = 0.d0
   correct_recip_ch_mp_ene = 0.d0
   correct_recip_ch_ch_ene = 0.d0
   correct_recip_ch_dh_ene = 0.d0
   correct_recip_dh_mp_ene = 0.d0
   correct_recip_dh_ch_ene = 0.d0
   correct_recip_dh_dh_ene = 0.d0
   correct_recip_mpole_mpole_ene = 0.d0
   correct_recip_mpole_herm_ene = 0.d0
   correct_recip_herm_herm_ene = 0.d0

   tot_mpole_mpole_ene = 0.d0
   tot_mpole_herm_ene = 0.d0
   tot_herm_herm_ene = 0.d0

!  TDDEL factor = exchange_factor
!  TDDEL exch_cutoff = exchange_cutoff

   !!! FOR EXCHANGE CALCULATION THE FIELDS INCLUDE THE FACTOR 
!  TBDEL Global_Chermite_field_EXCH = 0.d0
!  TBDEL tot_exchange_ene = 0.d0

   if ( tot_num_mpole_coeffs > 0 )then
      Global_multipole_field = 0.d0
   endif
   if ( tot_num_sumCH_coeffs > 0 )then
      Global_sumCH_field = 0.d0
   endif
   if ( tot_num_herm_Cprims > 0 )then
      Global_Chermite_field = 0.d0
   endif
   if ( tot_num_herm_Dprims > 0 )then
      Global_Dhermite_field = 0.d0
   endif

   site_frc(:,:) = 0.d0

   call FH_hermite_local_to_global()
   call Map_site_coords()
   call FH_SITE_SITE_LIST_eval(out_lun)

   factor = 1.d0 ! factor for Coulomb interactions

   correcting_energies = 0
   use_cutoff = 1
   full_interactions = 0
   call FH_COULOMB_SITE_SITE_CC( &
                 correcting_energies, &
                 use_cutoff, &
                 full_interactions, &
                 max_cutoff_direct_num_nghbrs, &
                 cutoff_direct_num_nghbrs, &
                 cutoff_direct_list_offset, &
                 cutoff_direct_list, &
                 cutoff_direct_ucell_transind, &
                 cutoff_mp_mp_ene, &
                 cutoff_mp_ch_ene, &
                 cutoff_ch_mp_ene, &
                 cutoff_ch_ch_ene,factor,cfgem_coul_virial)
   cutoff_mpole_mpole_ene = cutoff_mp_mp_ene
   cutoff_mpole_herm_ene = cutoff_mp_ch_ene + &
                           cutoff_ch_mp_ene 
   cutoff_herm_herm_ene = cutoff_ch_ch_ene
   ! next the recip sum terms
   call FH_RECIP_eval( &
                      recip_mp_mp_ene, &
                      recip_mp_ch_ene, &
                      recip_ch_mp_ene, &
                      recip_ch_ch_ene, &
                      recip_mp_dh_ene, &
                      recip_dh_mp_ene, &
                      recip_ch_dh_ene, &
                      recip_dh_ch_ene, &
                      recip_dh_dh_ene,factor,cfgem_coul_virial,out_lun)

   recip_mpole_mpole_ene = recip_mp_mp_ene
   recip_mpole_herm_ene = recip_mp_ch_ene + &
                          recip_mp_dh_ene + &
                          recip_ch_mp_ene + &
                          recip_dh_mp_ene
   recip_herm_herm_ene = recip_ch_ch_ene + &
                         recip_ch_dh_ene + &
                         recip_dh_ch_ene + &
                         recip_dh_dh_ene
   ! next the ewald potential volume dependent wigner terms
   call  Wigner_correction( &
                      wigner_mp_mp_ene, &
                      wigner_mp_ch_ene, &
                      wigner_ch_mp_ene, &
                      wigner_ch_ch_ene,factor)
   wigner_mpole_mpole_ene = wigner_mp_mp_ene
   wigner_mpole_herm_ene = wigner_mp_ch_ene + &
                          wigner_mp_dh_ene + &
                          wigner_ch_mp_ene + &
                          wigner_dh_mp_ene
   wigner_herm_herm_ene = wigner_ch_ch_ene + &
                         wigner_ch_dh_ene + &
                         wigner_dh_ch_ene + &
                         wigner_dh_dh_ene
   tot_wigner_ene = wigner_mpole_mpole_ene + &
                    wigner_mpole_herm_ene + &
                    wigner_herm_herm_ene
   ! now correct the masked interactions
   ! first the compact compact
   correcting_energies = 1
   use_cutoff = 0
   full_interactions = 0
   call FH_COULOMB_SITE_SITE_CC( &
                 correcting_energies, &
                 use_cutoff, &
                 full_interactions, &
                 max_correct_recip_num_nghbrs, &
                 correct_recip_num_nghbrs, &
                 correct_recip_list_offset, &
                 correct_recip_list, &
                 correct_recip_ucell_transind, &
                 correct_recip_mp_mp_ene, &
                 correct_recip_mp_ch_ene, &
                 correct_recip_ch_mp_ene, &
                 correct_recip_ch_ch_ene,factor,cfgem_coul_virial)
   call FH_COULOMB_SITE_SITE_CD_DD( &
                 correcting_energies, &
                 max_correct_recip_num_nghbrs, &
                 correct_recip_num_nghbrs, &
                 correct_recip_list_offset, &
                 correct_recip_list, &
                 correct_recip_ucell_transind, &
                 correct_recip_mp_dh_ene, &
                 correct_recip_dh_mp_ene, &
                 correct_recip_ch_dh_ene, &
                 correct_recip_dh_ch_ene, &
                 correct_recip_dh_dh_ene,factor,cfgem_coul_virial)
   correct_recip_mpole_mpole_ene = correct_recip_mp_mp_ene
   correct_recip_mpole_herm_ene = correct_recip_mp_ch_ene + &
                         correct_recip_mp_dh_ene + &
                         correct_recip_ch_mp_ene + &
                         correct_recip_dh_mp_ene
   correct_recip_herm_herm_ene = correct_recip_ch_ch_ene + &
                        correct_recip_ch_dh_ene + &
                        correct_recip_dh_ch_ene + &
                        correct_recip_dh_dh_ene
   ! finally the self terms
   correcting_energies = 1
   call FH_COULOMB_SITE_SITE_self_CC( &
                               correcting_energies, &
                               self_mp_mp_ene, &
                               self_mp_ch_ene, &
                               self_ch_mp_ene, &
                               self_ch_ch_ene,factor,cfgem_coul_virial)
   call FH_COULOMB_SITE_SITE_self_CD_DD( &
                               correcting_energies, &
                               self_mp_dh_ene, &
                               self_dh_mp_ene, &
                               self_ch_dh_ene, &
                               self_dh_ch_ene, &
                               self_dh_dh_ene,factor,cfgem_coul_virial)
   self_mpole_mpole_ene = self_mp_mp_ene
   self_mpole_herm_ene = self_mp_ch_ene + &
                         self_mp_dh_ene + &
                         self_ch_mp_ene + &
                         self_dh_mp_ene
   self_herm_herm_ene = self_ch_ch_ene + &
                        self_ch_dh_ene + &
                        self_dh_ch_ene + &
                        self_dh_dh_ene
   tot_mpole_mpole_ene = cutoff_mpole_mpole_ene + &
                         recip_mpole_mpole_ene + &
                         wigner_mpole_mpole_ene + &
                         correct_recip_mpole_mpole_ene + &
                         self_mpole_mpole_ene
   tot_mpole_herm_ene = cutoff_mpole_herm_ene + &
                         recip_mpole_herm_ene + &
                         wigner_mpole_herm_ene + &
                         correct_recip_mpole_herm_ene + &
                         self_mpole_herm_ene
   tot_herm_herm_ene = cutoff_herm_herm_ene + &
                         recip_herm_herm_ene + &
                         wigner_herm_herm_ene + &
                         correct_recip_herm_herm_ene + &
                         self_herm_herm_ene

!  if ( verbose == 1 )then
!     write(6,*)'----------------------------------------'
!     write(6,*)'mpole-mpole components: '
!     write(6,*)'cutoff_mpole_mpole_ene = ',cutoff_mpole_mpole_ene
!     write(6,*)'recip_mpole_mpole_ene = ',recip_mpole_mpole_ene
!     write(6,*)'wigner_mpole_mpole_ene = ',wigner_mpole_mpole_ene
!     write(6,*)'correct_recip_mpole_mpole_ene = ', &
!                correct_recip_mpole_mpole_ene
!     write(6,*)'self_mpole_mpole_ene = ',self_mpole_mpole_ene
!     write(6,*)'----------------------------------------'
!     write(6,*)'mpole-herm components: '
!     write(6,*)'cutoff_mpole_herm_ene = ',cutoff_mpole_herm_ene
!     write(6,*)'recip_mpole_herm_ene = ',recip_mpole_herm_ene
!     write(6,*)'wigner_mpole_herm_ene = ',wigner_mpole_herm_ene
!     write(6,*)'correct_recip_mpole_herm_ene = ', &
!                correct_recip_mpole_herm_ene
!     write(6,*)'self_mpole_herm_ene = ',self_mpole_herm_ene
!     write(6,*)'----------------------------------------'
!     write(6,*)'herm-herm components: '
!     write(6,*)'cutoff_herm_herm_ene = ',cutoff_herm_herm_ene
!     write(6,*)'recip_herm_herm_ene = ',recip_herm_herm_ene
!     write(6,*)'wigner_herm_herm_ene = ',wigner_herm_herm_ene
!     write(6,*)'correct_recip_herm_herm_ene = ', &
!                correct_recip_herm_herm_ene
!     write(6,*)'self_herm_herm_ene = ',self_herm_herm_ene
!  endif

   ! write out Coulomb energy

   cfgem_tot_coul_ene = tot_mpole_mpole_ene + tot_mpole_herm_ene + &
                        tot_herm_herm_ene

   ! Finalize the nonbond forces

   call Torque_force(cfgem_coul_ene_fld,cfgem_coul_virial)

!--- Begin output variable loading here...

   cfgem_coul_ene_fld = cfgem_coul_ene_fld + tot_wigner_ene

!!!

   cfgem_coul_recip_ene%mp_mp = recip_mp_mp_ene
   cfgem_coul_recip_ene%mp_ch = recip_mp_ch_ene
   cfgem_coul_recip_ene%ch_mp = recip_ch_mp_ene
   cfgem_coul_recip_ene%ch_ch = recip_ch_ch_ene
   cfgem_coul_recip_ene%mp_dh = recip_mp_dh_ene
   cfgem_coul_recip_ene%dh_mp = recip_dh_mp_ene
   cfgem_coul_recip_ene%ch_dh = recip_ch_dh_ene
   cfgem_coul_recip_ene%dh_ch = recip_dh_ch_ene
   cfgem_coul_recip_ene%dh_dh = recip_dh_dh_ene

!!!

   cfgem_coul_direct_ene%mp_mp = cutoff_mp_mp_ene
   cfgem_coul_direct_ene%mp_ch = cutoff_mp_ch_ene
   cfgem_coul_direct_ene%ch_mp = cutoff_ch_mp_ene
   cfgem_coul_direct_ene%ch_ch = cutoff_ch_ch_ene
   cfgem_coul_direct_ene%mp_dh = 0.d0
   cfgem_coul_direct_ene%dh_mp = 0.d0
   cfgem_coul_direct_ene%ch_dh = 0.d0
   cfgem_coul_direct_ene%dh_ch = 0.d0
   cfgem_coul_direct_ene%dh_dh = 0.d0

!!!

   cfgem_coul_adj_recip_ene%mp_mp = correct_recip_mp_mp_ene
   cfgem_coul_adj_recip_ene%mp_ch = correct_recip_mp_ch_ene
   cfgem_coul_adj_recip_ene%ch_mp = correct_recip_ch_mp_ene
   cfgem_coul_adj_recip_ene%ch_ch = correct_recip_ch_ch_ene
   cfgem_coul_adj_recip_ene%mp_dh = correct_recip_mp_dh_ene
   cfgem_coul_adj_recip_ene%dh_mp = correct_recip_dh_mp_ene
   cfgem_coul_adj_recip_ene%ch_dh = correct_recip_ch_dh_ene
   cfgem_coul_adj_recip_ene%dh_ch = correct_recip_dh_ch_ene
   cfgem_coul_adj_recip_ene%dh_dh = correct_recip_dh_dh_ene

!!!

   cfgem_coul_self_ene%mp_mp = self_mp_mp_ene
   cfgem_coul_self_ene%mp_ch = self_mp_ch_ene
   cfgem_coul_self_ene%ch_mp = self_ch_mp_ene
   cfgem_coul_self_ene%ch_ch = self_ch_ch_ene
   cfgem_coul_self_ene%mp_dh = self_mp_dh_ene
   cfgem_coul_self_ene%dh_mp = self_dh_mp_ene
   cfgem_coul_self_ene%ch_dh = self_ch_dh_ene
   cfgem_coul_self_ene%dh_ch = self_dh_ch_ene
   cfgem_coul_self_ene%dh_dh = self_dh_dh_ene

!!!


   cfgem_coul_wigner_ene%mp_mp = wigner_mp_mp_ene
   cfgem_coul_wigner_ene%mp_ch = wigner_mp_ch_ene
   cfgem_coul_wigner_ene%ch_mp = wigner_ch_mp_ene
   cfgem_coul_wigner_ene%ch_ch = wigner_ch_ch_ene
   cfgem_coul_wigner_ene%mp_dh = wigner_mp_dh_ene
   cfgem_coul_wigner_ene%dh_mp = wigner_dh_mp_ene
   cfgem_coul_wigner_ene%ch_dh = wigner_ch_dh_ene
   cfgem_coul_wigner_ene%dh_ch = wigner_dh_ch_ene
   cfgem_coul_wigner_ene%dh_dh = wigner_dh_dh_ene

!!!

   cfgem_coul_frc(:,:) = site_frc(:,:)

!!!

! cfgem_coul_virial(:,:) was direct-loaded..

!!!

!--- End output variable loading...

   return

end subroutine nonbond_eval_coul

!******************************************************************************
!
! Subroutine:   cfgem_exch
!
! Description:  Evaluates a described system for gem exchange energies and
!               forces at a single point in time..  Exchange is also known as
!               "overlap" in the code.
!
! Input Arguments:
!
!   aux_file, list_file, user_exch_file - gem input parameters file names
!   (see samples under test)
!
!   atom_cnt - count of atoms in system
!   atom_crds - coordinates for all atoms (Angstrom)
!   pbc_box_size - size of unit cell (Angstrom).  OPTIONAL, at end of arg list;
!                  if not present, values are taken from user_exch_file.
!
!   free_lun - a file logical unit number to be used by library temp file
!              operations (ie., you are not using it in the caller).
!   out_lun  - file logical unit number to be used by library for printed
!              output, if necessary.  The library will return actual data in
!              the output arguments, so this is used when something goes wrong.
!              Typical values are 6, or mdout in amber code.
!
! Output Arguments:
!
!   (ALL output arguments are zero'd on input)
!
!   cfgem_exch_ch_ch_ene - total exchange energy, Hartree
!   cfgem_exch_ene_fld - total exchange field energy, Hartree
!   cfgem_exch_frc - Array of total exchange forces, Hartree/Bohr
!
!******************************************************************************

end subroutine cfgem_coul

subroutine cfgem_exch(aux_file, list_file, user_exch_file, &
                      atom_cnt, atom_crds, free_lun, out_lun, &
                      cfgem_exch_ch_ch_ene, &
                      cfgem_exch_ene_fld, &
                      cfgem_exch_frc, &
                      pbc_box_size)

   use atoms,only : FH_ATOMS_init, &
                    FH_ATOMS_deallocate

   use user,only : FH_USER_read, &
                   FH_USER_first_check, &
                   FH_USER_second_check, &
                   dump_atom_forces, &
                   do_coulomb, &
                   do_overlap, &
                   calculation_method

   use auxiliary,only : FH_AUX_read, &
                        FH_AUX_get_rec_expon, &
                        FH_AUX_recip_auto_setup, &
                        FH_AUX_deallocate

   use sites,only : FH_SITES_readfile, &
                    FH_SITES_build, &
                    FH_SITES_deallocate

   use hermite, only : FH_hermite_load, &
                       FH_hermite_deallocate, &
                       FH_hermite_density_extent

   use site_site_list, only : FH_SITE_SITE_LIST_init, &
                              FH_SITE_SITE_LIST_deallocate

   use unit_cell, only : FH_UNIT_CELL_setup, &
                         FH_UNIT_CELL_deallocate

   use recip_sum, only : FH_RECIP_setup, &
                         FH_RECIP_deallocate

!! GAC: store GEM Coulomb and Exchange total Energy and Forces in 
!!      TINKER multipole vars
   !use energi,    only:  em
   !use deriv,     only:  dem
                        
   implicit none

! Formal arguments:

   ! atom_crds() in angstrom..
   ! free_lun is used by library as lun for reading input files.
   ! out_lun is lun for any printed output (should be error msgs only)
   ! pbc_box_size() in angstrom; if not present, use user_exch_file value.

   character(len=80), intent(in)        :: aux_file
   character(len=80), intent(in)        :: list_file
   character(len=80), intent(in)        :: user_exch_file
   integer,           intent(in)        :: atom_cnt
   double precision,  intent(in)        :: atom_crds(3, atom_cnt)
   integer,           intent(in)        :: free_lun
   integer,           intent(in)        :: out_lun

   double precision,    intent(out)     :: cfgem_exch_ch_ch_ene

   double precision,    intent(out)     :: cfgem_exch_ene_fld ! (for torque)

   double precision,    intent(inout)     :: cfgem_exch_frc(3, atom_cnt)


   double precision, optional,  intent(in)        :: pbc_box_size(3)

! Local variables:
   
   include "interact_type.fh"

   double precision, allocatable, save  :: atom_crds_bohr(:,:)

   ! NOTE: While the subroutines used in this function return a virial, the
   !       recip space coulomb virial is known to be incorrect; thus we
   !       currently return no virial information..

   double precision                     :: cfgem_exch_virial(3,3)

   integer              :: i, j, k
   integer              :: alloc_failed

   ! Go ahead and convert coordinates from Anstrom to Bohr, in locally
   ! allocated array; don't want to modify input!

   call FH_ATOMS_init(atom_cnt, atom_crds, out_lun)

   call FH_SITE_SITE_LIST_init()

   if (present(pbc_box_size)) then
      call FH_USER_read(user_exch_file, free_lun, out_lun, pbc_box_size)
   else
      call FH_USER_read(user_exch_file, free_lun, out_lun)
   end if

   call FH_USER_first_check(do_overlap, out_lun)

   if (do_coulomb .ne. 0) then
     write(out_lun)'CFGEM_EXCH: do_coulomb set in user exchange file!'
     stop
   end if

   ! Zero all output arguments:

   cfgem_exch_ch_ch_ene = 0.d0
   cfgem_exch_ene_fld = 0.d0
   cfgem_exch_frc(:,:) = 0.d0
   cfgem_exch_virial(:,:) = 0.d0

   call FH_AUX_read(aux_file, free_lun, out_lun)

   call FH_AUX_recip_auto_setup(out_lun) ! setup based on auxiliary info

   call FH_USER_second_check(out_lun) ! check on the recip settings

   call FH_AUX_get_rec_expon(out_lun)

   call FH_SITES_readfile(list_file, free_lun, out_lun)

   call FH_hermite_load(out_lun)

   call FH_MCMUR_DAV_init()

   call FH_SITES_build()

   call FH_hermite_density_extent(out_lun)

   call FH_UNIT_CELL_setup(out_lun)

   call FH_RECIP_setup(out_lun)

   call nonbond_eval_exch()

   ! Add Exchange Energy/Force to Coulomb Contribution in TINKER vars
   !em = em + cfgem_exch_ch_ch_ene
   !dem = dem + cfgem_exch_frc

   ! Deallocate:

   call FH_ATOMS_deallocate()
   call FH_AUX_deallocate()
   call FH_SITES_deallocate()
   call FH_hermite_deallocate()
   call FH_UNIT_CELL_deallocate()
   call FH_RECIP_deallocate()
   call FH_SITE_SITE_LIST_deallocate()

   return

contains

subroutine nonbond_eval_exch()

   use user, only : verbose,do_coulomb,calculation_method,&
                    do_overlap,exchange_factor,exchange_cutoff

   use site_site_list, only : FH_SITE_SITE_LIST_eval

   use recip_sum, only : FH_RECIP_eval

   use sites, only : site_frc

   use hermite, only : FH_hermite_local_to_global, &
                       tot_num_mpole_coeffs, &
                       tot_num_sumCH_coeffs, &
                       tot_num_herm_Cprims, &
                       tot_num_herm_Dprims, &
                       Global_multipole_field, &
                       Global_sumCH_field, &
                       Global_Chermite_field, &
                       Global_Chermite_field_EXCH, &
                       Global_Dhermite_field

   use coulomb_site_site, only : FH_COULOMB_SITE_SITE_CC, &
                                 FH_COULOMB_SITE_SITE_self_CC, &
                                 FH_COULOMB_SITE_SITE_CD_DD, &
                                 FH_COULOMB_SITE_SITE_self_CD_DD

   use exchange_site_site, only : FH_EXCHANGE_SITE_SITE_CC

   use site_site_list, only : &
                            max_unmasked_pairs_num_nghbrs, &
                            unmasked_pairs_list, &
                            unmasked_pairs_ucell_transind, &
                            unmasked_pairs_num_nghbrs, &
                            unmasked_pairs_list_offset, &
                            max_cutoff_direct_num_nghbrs, &
                            cutoff_direct_list, &
                            cutoff_direct_ucell_transind, &
                            cutoff_direct_num_nghbrs, &
                            cutoff_direct_list_offset, &
                            max_emulate_recip_num_nghbrs, &
                            emulate_recip_list, &
                            emulate_recip_ucell_transind, &
                            emulate_recip_num_nghbrs, &
                            emulate_recip_list_offset, &
                            max_correct_recip_num_nghbrs, &
                            correct_recip_list, &
                            correct_recip_ucell_transind, &
                            correct_recip_num_nghbrs, &
                            correct_recip_list_offset

   implicit none

! Formal arguments:

! Local variables:

   integer :: correcting_energies, &
              use_cutoff, &
              full_interactions

   include "interact_type.fh"

   !!! FOR EXCHANGE CALCULATION THE FIELDS INCLUDE THE FACTOR 

   Global_Chermite_field_EXCH(:) = 0.d0

   if ( tot_num_mpole_coeffs > 0 )then
      Global_multipole_field = 0.d0
   endif
   if ( tot_num_sumCH_coeffs > 0 )then
      Global_sumCH_field = 0.d0
   endif
   if ( tot_num_herm_Cprims > 0 )then
      Global_Chermite_field = 0.d0
   endif
   if ( tot_num_herm_Dprims > 0 )then
      Global_Dhermite_field = 0.d0
   endif

   site_frc(:,:) = 0.d0

   call FH_hermite_local_to_global()
   call Map_site_coords()
   call FH_SITE_SITE_LIST_eval(out_lun)

   ! DO EXCHANGE CALCULATION

   ! DO ALL CALCULATIONS IN REAL SPACE
   ! in this case all hermites are compact since CD_split_expon is 0
   !
   ! NOTE: Input allows selection of cutoff to avoid calculation of
   !       integrals. For a box of 4096 waters, the energy error using
   !       a cutoff of 5 Angstrom is .216 kcal/mol (0.0004 %), 
   !       errors in forces are =< 10^-4 kcal/mol, with
   !       a speedup of 2 orders of magnitude (2 s vs. 125 s)

   correcting_energies = 0
   use_cutoff = 1
   full_interactions = 1

   call FH_EXCHANGE_SITE_SITE_CC( &
                    correcting_energies, &
                    use_cutoff, &
                    full_interactions, &
                    max_unmasked_pairs_num_nghbrs, &
                    unmasked_pairs_num_nghbrs, &
                    unmasked_pairs_list_offset, &
                    unmasked_pairs_list, &
                    unmasked_pairs_ucell_transind, &
                    cfgem_exch_ch_ch_ene, &
                    exchange_factor, &
                    exchange_cutoff, &
                    cfgem_exch_virial)
      
   ! Finalize the nonbond forces

   call Torque_force_EXCH(cfgem_exch_ene_fld,cfgem_exch_virial)

   cfgem_exch_frc(:,:) = site_frc(:,:)

   return

end subroutine nonbond_eval_exch

end subroutine cfgem_exch

!---------------------------------------------------------
subroutine Map_site_coords()

   use sites, only : num_sites,site_crd,map_site_crd,frac_site_crd
   use user, only : pbc
   use unit_cell, only : ucell,recip

   implicit none

   integer n
   double precision f1,f2,f3

   if ( pbc == 1 )then
      do n = 1,num_sites
         f1 = site_crd(1,n)*recip(1,1) + site_crd(2,n)*recip(2,1) + &
              site_crd(3,n)*recip(3,1)
         f2 = site_crd(1,n)*recip(1,2) + site_crd(2,n)*recip(2,2) + &
              site_crd(3,n)*recip(3,2)
         f3 = site_crd(1,n)*recip(1,3) + site_crd(2,n)*recip(2,3) + &
              site_crd(3,n)*recip(3,3)
         f1 = f1 - dnint(f1)
         f2 = f2 - dnint(f2)
         f3 = f3 - dnint(f3)
         frac_site_crd(1,n) = f1
         frac_site_crd(2,n) = f2
         frac_site_crd(3,n) = f3
         map_site_crd(1,n) = ucell(1,1)*f1 + ucell(1,2)*f2 + ucell(1,3)*f3
         map_site_crd(2,n) = ucell(2,1)*f1 + ucell(2,2)*f2 + ucell(2,3)*f3
         map_site_crd(3,n) = ucell(3,1)*f1 + ucell(3,2)*f2 + ucell(3,3)*f3
      enddo! n = 1,num_sites
   else
      do n = 1,num_sites
         map_site_crd(1,n) = site_crd(1,n)
         map_site_crd(2,n) = site_crd(2,n)
         map_site_crd(3,n) = site_crd(3,n)
      enddo
   endif! pbc == 1 
end subroutine Map_site_coords

!---------------------------------------------------------
subroutine Torque_force(ene_fld, virial)
   use sites, only : FH_SITES_de_drotsite_to_sitefrc, &
                     FH_SITES_sitefrc_to_frc,FH_SITES_add_extrapt_torque
   use hermite, only : FH_hermite_field_energy, &
                       FH_hermite_field_de_drot

   implicit none

   double precision,intent(out) :: ene_fld
   double precision,intent(inout) :: virial(3,3)

   integer :: pass
   ene_fld = 0.d0
   call FH_hermite_field_energy(ene_fld) ! THIS GIVES TORQUE ENERGY
   ! global hermite field dotted against coeffs gives energy
   ! de_drot given by hermite field dotted against derivs of coeffs wrt rot
   call FH_hermite_field_de_drot() ! THIS GIVES TORQUE FIELD
   pass = 1
   call FH_SITES_de_drotsite_to_sitefrc(pass,virial)
   ! next add the extra point site_frc torque contribution
   call FH_SITES_add_extrapt_torque()
   pass = 2
   call FH_SITES_de_drotsite_to_sitefrc(pass,virial)
   ! finally restrict forces to atoms
   call FH_SITES_sitefrc_to_frc()

end subroutine Torque_force

!---------------------------------------------------------
subroutine Torque_force_EXCH(ene_fld, virial)

   use sites, only : FH_SITES_de_drotsite_to_sitefrc, &
                     FH_SITES_sitefrc_to_frc,FH_SITES_add_extrapt_torque
   use hermite, only : FH_hermite_field_energy_EXCH, &
                       FH_hermite_field_de_drot_EXCH

   implicit none

   double precision,intent(out) :: ene_fld
   double precision,intent(inout) :: virial(3,3)

   integer :: pass
   ene_fld = 0.d0
   call FH_hermite_field_energy_EXCH(ene_fld) ! THIS GIVES TORQUE ENERGY
   ! global hermite field dotted against coeffs gives energy
   ! de_drot given by hermite field dotted against derivs of coeffs wrt rot
   call FH_hermite_field_de_drot_EXCH() ! THIS GIVES TORQUE FIELD
   pass = 1
   call FH_SITES_de_drotsite_to_sitefrc(pass,virial)
   ! next add the extra point site_frc torque contribution
   call FH_SITES_add_extrapt_torque()
   pass = 2
   call FH_SITES_de_drotsite_to_sitefrc(pass,virial)
   ! finally restrict forces to atoms
   call FH_SITES_sitefrc_to_frc()

end subroutine Torque_force_EXCH

!---------------------------------------------------------
subroutine Wigner_correction( &
                            mp_mp_ene, &
                            mp_ch_ene, &
                            ch_mp_ene, &
                            ch_ch_ene,factor)

   use user, only : pbc,CD_split_expon
   use sites, only : num_sites
   use hermite, only :  tot_num_mpole_coeffs, &
                        mpole_coeff_off_for_site, &
                        mpole_order_for_site, &
                        Global_multipole_coeff, &
                        tot_num_herm_Cprims, &
                        hermite_order_for_Cprim, &
                        herm_coeff_offset_of_Cprim, &
                        hermite_expon_for_Cprim, &
                        Global_Chermite_coeff
   use unit_cell, only : volume

   implicit none

   double precision,intent(in) :: factor
   double precision,intent(out) :: &
                                  mp_mp_ene, &
                                  mp_ch_ene, &
                                  ch_mp_ene, &
                                  ch_ch_ene

   double precision :: ewald_expon,pi, &
                       ch_expon1,ch_expon2,charge1,charge2
   integer :: n1,n2,np1,np2,off1,off2,offc1,offc2

   if ( pbc /= 1 )return !nothing to do here

   pi = 3.14159265358979323846d0
   ewald_expon = 0.5d0*CD_split_expon

   ! first mpole-mpole
   mp_mp_ene = 0.d0
   if ( tot_num_mpole_coeffs > 0 )then
      do n1 = 1,num_sites
         off1 = mpole_coeff_off_for_site(n1)
         if ( mpole_order_for_site(n1) > 0 )then
            charge1 = Global_multipole_coeff(off1+1)
            do n2 = 1,num_sites
               off2 = mpole_coeff_off_for_site(n2)
               if ( mpole_order_for_site(n2) > 0 )then
                  charge2 = Global_multipole_coeff(off2+1)
                  ! positive interactions e.g. nuclear-nuclear
                  mp_mp_ene = mp_mp_ene + charge1*charge2/ewald_expon
               endif !( mpole_order_for_site(n2) > 0 )then
            enddo !n2 = 1,num_sites
         endif !( mpole_order_for_site(n1) > 0 )then
      enddo !n1 = 1,num_sites
      ! negative sign since this is a compensating term
      mp_mp_ene = -0.5d0 * pi * mp_mp_ene / volume
   endif !( tot_num_mpole_coeffs > 0 )then

   ! next the mpole-compact hermite
   mp_ch_ene = 0.d0
   ch_mp_ene = 0.d0
   if ( (tot_num_herm_Cprims > 0) .and. (tot_num_mpole_coeffs > 0) )then
      do n1 = 1,num_sites
         off1 = mpole_coeff_off_for_site(n1)
         if ( mpole_order_for_site(n1) > 0 )then
            charge1 = Global_multipole_coeff(off1+1)
            do np2 = 1,tot_num_herm_Cprims
               offc2 = herm_coeff_offset_of_Cprim(np2)
               if ( hermite_order_for_Cprim(np2) > 0 )then
                  charge2 = Global_Chermite_coeff(offc2+1)
                  ch_expon2 = hermite_expon_for_Cprim(np2)
                  ! negative interactions e.g. nuclear-electron
                  mp_ch_ene = mp_ch_ene - &
                     charge1*charge2*(1.d0/ewald_expon - 1.d0/ch_expon2)
               endif !( hermite_order_for_Cprim(np1) > 0 )then
            enddo !np1 = 1,tot_num_herm_Cprims
         endif !( mpole_order_for_site(n1) > 0 )then
      enddo !n1 = 1,num_sites
      ! negative sign since this is a compensating term
      mp_ch_ene = -0.5d0 * pi * mp_ch_ene / volume
      ch_mp_ene = mp_ch_ene ! by symmetry (just reverse the above 2 loops)
   endif !( (tot_num_herm_Cprims > 0) .and. (tot_num_mpole_coeffs > 0) )then

   ! finally the compact hermite-compact hermite
   ch_ch_ene = 0.d0
   if ( (tot_num_herm_Cprims > 0) )then
      do np1 = 1,tot_num_herm_Cprims
         offc1 = herm_coeff_offset_of_Cprim(np1)
         if ( hermite_order_for_Cprim(np1) > 0 )then
            charge1 = Global_Chermite_coeff(offc1+1)
            ch_expon1 = hermite_expon_for_Cprim(np1)
            do np2 = 1,tot_num_herm_Cprims
               offc2 = herm_coeff_offset_of_Cprim(np2)
               if ( hermite_order_for_Cprim(np2) > 0 )then
                  charge2 = Global_Chermite_coeff(offc2+1)
                  ch_expon2 = hermite_expon_for_Cprim(np2)
                  ! positive interactions e.g. electron-electron
                  ch_ch_ene = ch_ch_ene + charge1*charge2* &
                      (1.d0/ewald_expon - 1.d0/ch_expon1 - 1.d0/ch_expon2)
               endif !( hermite_order_for_Cprim(np2) > 0 )then
            enddo !np2 = 1,tot_num_herm_Cprims
         endif !( hermite_order_for_Cprim(np1) > 0 )then
      enddo !np1 = 1,tot_num_herm_Cprims
      ! negative sign since this is a compensating term
      ch_ch_ene = -0.5d0 * pi * ch_ch_ene / volume
   endif !( (tot_num_herm_Cprims > 0) .and. (tot_num_mpole_coeffs > 0) )then

end subroutine Wigner_correction

end module cfgem_libmod
