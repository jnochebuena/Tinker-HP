module ewald_recip
   implicit none
   private

   integer, save :: ewald_mmax
   double precision,save, allocatable :: ewald_cosf(:,:,:), &
                                         ewald_sinf(:,:,:)
   public FH_EWALD_RECIP_setup, &
          FH_EWALD_RECIP_FT_density, &
          FH_EWALD_RECIP_ene_force_field, &
          FH_EWALD_RECIP_deallocate

contains
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!-----------------------------------------------------------------
! SETUP SUBROUTINES
!-----------------------------------------------------------------
!-----------------------------------------------------------------

!-----------------------------------------------------------------

subroutine FH_EWALD_RECIP_setup(out_lun)

   use user,only : num_ewald_grid_types, &
                   struc_fac_method_for_grid_type, &
                   num_prim_grid_types, &
                   nfft1_for_grid_type, &
                   nfft2_for_grid_type, &
                   nfft3_for_grid_type
   use sites, only : num_sites

   implicit none

   integer, intent(in)  :: out_lun

   include 'structure_factor_type.fh'

   integer :: gt,ier1,ier2,mmax

   if ( num_ewald_grid_types == 0 )return ! nothing further to do here

   ewald_mmax = -1
   do gt = 1,num_prim_grid_types
      if ( struc_fac_method_for_grid_type(gt) == SF_EWALD )then
         call FH_EWALD_RECIP_get_mmax(nfft1_for_grid_type(gt), &
                                      nfft2_for_grid_type(gt), &
                                      nfft3_for_grid_type(gt), &
                                      mmax)
         if ( mmax > ewald_mmax )ewald_mmax = mmax
      endif
   enddo
  
   if ( ewald_mmax > 0 )then
      allocate(ewald_cosf(3,num_sites,0:ewald_mmax),stat=ier1)
      allocate(ewald_sinf(3,num_sites,0:ewald_mmax),stat=ier2)
      if ( ier1 /= 0 .or. ier2 /= 0 )then
         write(out_lun,*)'FH_RECIP_ewald_grid_setup: allocate fails!'
         stop
      endif
   endif

end subroutine FH_EWALD_RECIP_setup
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!-----------------------------------------------------------------
! EVALUATION SUBROUTINES
!-----------------------------------------------------------------
!-----------------------------------------------------------------

subroutine FH_EWALD_RECIP_FT_density()

   use sites, only : num_sites,site_crd
   use unit_cell, only : recip,volume
   use hermite, only :  &
                       ! mpoles
                       tot_num_mpole_coeffs, &
                       Global_multipole_coeff, &
                       compact_rec_mpole_expon, &
                       diffuse_rec_mpole_expon, &
                       ! sumCH
                       tot_num_sumCH_coeffs, &
                       Global_sumCH_coeff, &
                       compact_rec_sumCH_expon, &
                       diffuse_rec_sumCH_expon, &
                       ! compact hermites
                       tot_num_herm_Cprims, &
                       Global_Chermite_coeff, &
                       ! diffuse hermites
                       tot_num_herm_Dprims, &
                       Global_Dhermite_coeff
   use user, only : struc_fac_method_for_grid_type, &
                    nfft1_for_grid_type, &
                    nfft2_for_grid_type, &
                    nfft3_for_grid_type, &
                    num_HC_prim_grids, &
                    grid_type_for_HC_prim_grid, &
                    num_HD_prim_grids, &
                    grid_type_for_HD_prim_grid, &
                    grid_type_for_MPOLES, &
                    grid_type_for_sumCH, &
                    num_ewald_grid_types
   use recip_grids, only : nfftdim1_for_gridtype, &
                           nfftdim2_for_gridtype, &
                           nfftdim3_for_gridtype, &
                           ! compact hermites
                           num_prims_for_HC_grid, &
                           prim_offset_for_HC_grid, &
                           site_that_owns_HC_grid_prim, &
                           Gcoeff_offset_for_HC_grid_prim, &
                           herm_order_for_HC_grid_prim, &
                           herm_expon_for_HC_grid_prim, &
                           HC_grid_offset_for_HC_grid,HC_grid, &
                           ! diffuse hermites
                           num_prims_for_HD_grid, &
                           prim_offset_for_HD_grid, &
                           site_that_owns_HD_grid_prim, &
                           herm_order_for_HD_grid_prim, &
                           herm_expon_for_HD_grid_prim, &
                           Gcoeff_offset_for_HD_grid_prim, &
                           HD_grid_offset_for_HD_grid,HD_grid, &
                           ! mpoles
                           num_mpoles_for_mpole_grid, &
                           site_that_owns_gridded_mpole, &
                           Mpole_order_for_gridded_mpole, &
                           Gcoeff_offset_for_gridded_mpole, &
                           MC_grid, MD_grid, &
                           ! sumCH
                           num_sumCH_for_sumCH_grid, &
                           site_that_owns_gridded_sumCH, &
                           sumCH_order_for_gridded_sumCH, &
                           Gcoeff_offset_for_gridded_sumCH, &
                           sumCH_C_grid, sumCH_D_grid

   implicit none

   include 'structure_factor_type.fh'

   integer g,gt,off,off1,n
   double precision :: mpole_exponent(num_mpoles_for_mpole_grid)
   double precision :: sumCH_exponent(num_sumCH_for_sumCH_grid)

   if ( num_ewald_grid_types == 0 )return ! nothing further to do here

   call FH_EWALD_RECIP_phase_factors(num_sites,ewald_mmax, &
                                     recip,site_crd, &
                                     ewald_cosf,ewald_sinf)

   ! first the mpole grids
   if ( tot_num_mpole_coeffs > 0 )then
      gt = grid_type_for_MPOLES
      if ( struc_fac_method_for_grid_type(gt) == SF_EWALD )then
         off = 0 ! only one grid-grid offset is zero
         off1 = 0 ! only one grid-prim_offset for grid is zero
         ! first the compact mpoles
         do n = 1,num_mpoles_for_mpole_grid
            mpole_exponent(n) = compact_rec_mpole_expon
         enddo
         call FH_EWALD_RECIP_struc_fac(  &
                               num_sites,ewald_mmax, &
                               num_mpoles_for_mpole_grid, &
                               off1, &
                               Global_multipole_coeff, &
                               Gcoeff_offset_for_gridded_mpole, &
                               site_that_owns_gridded_mpole, &
                               Mpole_order_for_gridded_mpole, &
                               mpole_exponent, &
                               recip,volume, &
                               ewald_cosf,ewald_sinf, &
                               nfft1_for_grid_type(gt), &
                               nfft2_for_grid_type(gt), &
                               nfft3_for_grid_type(gt), &
                               nfftdim1_for_gridtype(gt), &
                               nfftdim2_for_gridtype(gt), &
                               nfftdim2_for_gridtype(gt), &
                               MC_grid(off+1))
         ! next the diffuse mpoles
         do n = 1,num_mpoles_for_mpole_grid
            mpole_exponent(n) = diffuse_rec_mpole_expon
         enddo
         call FH_EWALD_RECIP_struc_fac(  &
                               num_sites,ewald_mmax, &
                               num_mpoles_for_mpole_grid, &
                               off1, &
                               Global_multipole_coeff, &
                               Gcoeff_offset_for_gridded_mpole, &
                               site_that_owns_gridded_mpole, &
                               Mpole_order_for_gridded_mpole, &
                               mpole_exponent, &
                               recip,volume, &
                               ewald_cosf,ewald_sinf, &
                               nfft1_for_grid_type(gt), &
                               nfft2_for_grid_type(gt), &
                               nfft3_for_grid_type(gt), &
                               nfftdim1_for_gridtype(gt), &
                               nfftdim2_for_gridtype(gt), &
                               nfftdim2_for_gridtype(gt), &
                               MD_grid(off+1))
      endif!struc_fac_method_for_grid_type(gt) == SF_EWALD
   endif !( tot_num_mpole_coeffs > 0 )then

   ! next the sumCH grids
   if ( tot_num_sumCH_coeffs > 0 )then
      gt = grid_type_for_sumCH
      if ( struc_fac_method_for_grid_type(gt) == SF_EWALD )then
         off = 0 ! only one grid-grid offset is zero
         off1 = 0 ! only one grid-prim_offset for grid is zero
         ! first the compact mpoles
         do n = 1,num_sumCH_for_sumCH_grid
            sumCH_exponent(n) = compact_rec_sumCH_expon
         enddo
         call FH_EWALD_RECIP_struc_fac(  &
                               num_sites,ewald_mmax, &
                               num_sumCH_for_sumCH_grid, &
                               off1, &
                               Global_sumCH_coeff, &
                               Gcoeff_offset_for_gridded_sumCH, &
                               site_that_owns_gridded_sumCH, &
                               sumCH_order_for_gridded_sumCH, &
                               sumCH_exponent, &
                               recip,volume, &
                               ewald_cosf,ewald_sinf, &
                               nfft1_for_grid_type(gt), &
                               nfft2_for_grid_type(gt), &
                               nfft3_for_grid_type(gt), &
                               nfftdim1_for_gridtype(gt), &
                               nfftdim2_for_gridtype(gt), &
                               nfftdim2_for_gridtype(gt), &
                               sumCH_C_grid(off+1))
         ! next the diffuse mpoles
         do n = 1,num_sumCH_for_sumCH_grid
            sumCH_exponent(n) = diffuse_rec_sumCH_expon
         enddo
         call FH_EWALD_RECIP_struc_fac(  &
                               num_sites,ewald_mmax, &
                               num_sumCH_for_sumCH_grid, &
                               off1, &
                               Global_sumCH_coeff, &
                               Gcoeff_offset_for_gridded_sumCH, &
                               site_that_owns_gridded_sumCH, &
                               sumCH_order_for_gridded_sumCH, &
                               sumCH_exponent, &
                               recip,volume, &
                               ewald_cosf,ewald_sinf, &
                               nfft1_for_grid_type(gt), &
                               nfft2_for_grid_type(gt), &
                               nfft3_for_grid_type(gt), &
                               nfftdim1_for_gridtype(gt), &
                               nfftdim2_for_gridtype(gt), &
                               nfftdim2_for_gridtype(gt), &
                               sumCH_D_grid(off+1))
      endif!struc_fac_method_for_grid_type(gt) == SF_EWALD
   endif !( tot_num_sumCH_coeffs > 0 )then

   if ( tot_num_herm_Cprims > 0 )then
      do g = 1,num_HC_prim_grids
         off = HC_grid_offset_for_HC_grid(g)
         gt = grid_type_for_HC_prim_grid(g)
         if ( struc_fac_method_for_grid_type(gt) == SF_EWALD )then
            call FH_EWALD_RECIP_struc_fac(  &
                               num_sites,ewald_mmax, &
                               num_prims_for_HC_grid(g), &
                               prim_offset_for_HC_grid(g), &
                               Global_CHermite_coeff, &
                               Gcoeff_offset_for_HC_grid_prim, &
                               site_that_owns_HC_grid_prim, &
                               herm_order_for_HC_grid_prim, &
                               herm_expon_for_HC_grid_prim, &
                               recip,volume, &
                               ewald_cosf,ewald_sinf, &
                               nfft1_for_grid_type(gt), &
                               nfft2_for_grid_type(gt), &
                               nfft3_for_grid_type(gt), &
                               nfftdim1_for_gridtype(gt), &
                               nfftdim2_for_gridtype(gt), &
                               nfftdim2_for_gridtype(gt), &
                               HC_grid(off+1))
         endif!struc_fac_method_for_grid_type(gt) == SF_EWALD
      enddo !g = 1,num_HC_prim_grids
   endif !( tot_num_herm_Cprims > 0 )then

   if ( tot_num_herm_Dprims > 0 )then
      do g = 1,num_HD_prim_grids
         off = HD_grid_offset_for_HD_grid(g)
         gt = grid_type_for_HD_prim_grid(g)
         if ( struc_fac_method_for_grid_type(gt) == SF_EWALD )then
            call FH_EWALD_RECIP_struc_fac(  &
                               num_sites,ewald_mmax, &
                               num_prims_for_HD_grid(g), &
                               prim_offset_for_HD_grid(g), &
                               Global_DHermite_coeff, &
                               Gcoeff_offset_for_HD_grid_prim, &
                               site_that_owns_HD_grid_prim, &
                               herm_order_for_HD_grid_prim, &
                               herm_expon_for_HD_grid_prim, &
                               recip,volume, &
                               ewald_cosf,ewald_sinf, &
                               nfft1_for_grid_type(gt), &
                               nfft2_for_grid_type(gt), &
                               nfft3_for_grid_type(gt), &
                               nfftdim1_for_gridtype(gt), &
                               nfftdim2_for_gridtype(gt), &
                               nfftdim2_for_gridtype(gt), &
                               HD_grid(off+1))
         endif!struc_fac_method_for_grid_type(gt) == SF_EWALD
      enddo !g = 1,num_HD_prim_grids
   endif !( tot_num_herm_Dprims > 0 )then

end subroutine FH_EWALD_RECIP_FT_density
!-----------------------------------------------------------------
subroutine FH_EWALD_RECIP_ene_force_field(energy)

   use sites, only : num_sites,site_crd,site_frc
   use unit_cell, only : recip,volume
   use hermite, only :  &
                       !mpoles
                       tot_num_mpole_coeffs, &
                       Global_multipole_coeff, &
                       Global_multipole_field, &
                       compact_rec_mpole_expon, &
                       !sumCH
                       tot_num_sumCH_coeffs, &
                       Global_sumCH_coeff, &
                       Global_sumCH_field, &
                       compact_rec_sumCH_expon, &
                       !compact hermites
                       tot_num_herm_Cprims, &
                       Global_CHermite_coeff, &
                       Global_CHermite_field, &
                       !diffuse hermites
                       tot_num_herm_Dprims, &
                       Global_DHermite_coeff, &
                       Global_DHermite_field
                       
                       
   use user, only : struc_fac_method_for_grid_type, &
                    nfft1_for_grid_type, &
                    nfft2_for_grid_type, &
                    nfft3_for_grid_type, &
                    num_HC_prim_grids, &
                    grid_type_for_HC_prim_grid, &
                    num_HD_prim_grids, &
                    grid_type_for_HD_prim_grid, &
                    grid_type_for_MPOLES, &
                    grid_type_for_sumCH, &
                    num_ewald_grid_types
   use recip_grids, only : nfftdim1_for_gridtype, &
                           nfftdim2_for_gridtype, &
                           nfftdim3_for_gridtype, &
                           ! mpoles
                           num_mpoles_for_mpole_grid, &
                           site_that_owns_gridded_mpole, &
                           Mpole_order_for_gridded_mpole, &
                           Gcoeff_offset_for_gridded_mpole, &
                           MC_grid,  &
                           MD_grid, &
                           ! sumCH
                           num_sumCH_for_sumCH_grid, &
                           site_that_owns_gridded_sumCH, &
                           sumCH_order_for_gridded_sumCH, &
                           Gcoeff_offset_for_gridded_sumCH, &
                           sumCH_C_grid,  &
                           sumCH_D_grid, &
                           ! compact hermites
                           num_prims_for_HC_grid, &
                           prim_offset_for_HC_grid, &
                           site_that_owns_HC_grid_prim, &
                           Gcoeff_offset_for_HC_grid_prim, &
                           herm_order_for_HC_grid_prim, &
                           herm_expon_for_HC_grid_prim, &
                           HC_grid_offset_for_HC_grid, &
                           HC_grid, &
                           ! diffuse hermites
                           num_prims_for_HD_grid, &
                           prim_offset_for_HD_grid, &
                           site_that_owns_HD_grid_prim, &
                           herm_order_for_HD_grid_prim, &
                           herm_expon_for_HD_grid_prim, &
                           Gcoeff_offset_for_HD_grid_prim, &
                           HD_grid_offset_for_HD_grid, &
                           HD_grid

   implicit none

   double precision,intent(inout) :: energy 

   include 'structure_factor_type.fh'

   integer :: n,g,gt,off,off1
   double precision :: mpole_exponent(num_mpoles_for_mpole_grid)
   double precision :: sumCH_exponent(num_sumCH_for_sumCH_grid)
   double precision :: factor
   !if ( energy_type == OVERLAP_ene_type )then
      !factor = exchange_factor
   !else
      factor = 1.d0
   !endif

   ! first the multipoles
   if ( tot_num_mpole_coeffs > 0 )then
      gt = grid_type_for_MPOLES
      if ( struc_fac_method_for_grid_type(gt) == SF_EWALD )then
         off = 0 ! just one grid
         off1 = 0 ! mpole offset is also 0, since one grid
         do n = 1,num_mpoles_for_mpole_grid
            mpole_exponent(n) = compact_rec_mpole_expon
         enddo
         call FH_EWALD_RECIP_one_ene_frc(  &
                               num_sites, &
                               ewald_mmax, &
                               num_mpoles_for_mpole_grid, &
                               off1, &
                               Global_multipole_coeff, &
                               Gcoeff_offset_for_gridded_mpole, &
                               site_that_owns_gridded_mpole, &
                               Mpole_order_for_gridded_mpole, &
                               mpole_exponent, &
                               recip, &
                               volume, &
                               ewald_cosf, &
                               ewald_sinf, &
                               nfft1_for_grid_type(gt), &
                               nfft2_for_grid_type(gt), &
                               nfft3_for_grid_type(gt), &
                               nfftdim1_for_gridtype(gt), &
                               nfftdim2_for_gridtype(gt), &
                               nfftdim2_for_gridtype(gt), &
                               MC_grid(off+1), &
                               factor, &
                               energy, &
                               site_frc, &
                               Global_multipole_field)
      endif !( struc_fac_method_for_grid_type(gt) == SF_EWALD )then
   endif !( tot_num_mpole_coeffs > 0 )then

   ! next the sumCH
   if ( tot_num_sumCH_coeffs > 0 )then
      gt = grid_type_for_sumCH
      if ( struc_fac_method_for_grid_type(gt) == SF_EWALD )then
         off = 0 ! just one grid
         off1 = 0 ! sumCH offset is also 0, since one grid
         do n = 1,num_sumCH_for_sumCH_grid
            sumCH_exponent(n) = compact_rec_sumCH_expon
         enddo
         call FH_EWALD_RECIP_one_ene_frc(  &
                               num_sites, &
                               ewald_mmax, &
                               num_sumCH_for_sumCH_grid, &
                               off1, &
                               Global_sumCH_coeff, &
                               Gcoeff_offset_for_gridded_sumCH, &
                               site_that_owns_gridded_sumCH, &
                               sumCH_order_for_gridded_sumCH, &
                               sumCH_exponent, &
                               recip, &
                               volume, &
                               ewald_cosf, &
                               ewald_sinf, &
                               nfft1_for_grid_type(gt), &
                               nfft2_for_grid_type(gt), &
                               nfft3_for_grid_type(gt), &
                               nfftdim1_for_gridtype(gt), &
                               nfftdim2_for_gridtype(gt), &
                               nfftdim2_for_gridtype(gt), &
                               sumCH_C_grid(off+1), &
                               factor, &
                               energy, &
                               site_frc, &
                               Global_sumCH_field)
      endif !( struc_fac_method_for_grid_type(gt) == SF_EWALD )then
   endif !( tot_num_sumCH_coeffs > 0 )then

   ! next the compact hermites
   if ( tot_num_herm_Cprims > 0 )then
      do g = 1,num_HC_prim_grids
         off = HC_grid_offset_for_HC_grid(g)
         gt = grid_type_for_HC_prim_grid(g)
         if ( struc_fac_method_for_grid_type(gt) == SF_EWALD )then
            call FH_EWALD_RECIP_one_ene_frc(  &
                               num_sites, &
                               ewald_mmax, &
                               num_prims_for_HC_grid(g), &
                               prim_offset_for_HC_grid(g), &
                               Global_CHermite_coeff, &
                               Gcoeff_offset_for_HC_grid_prim, &
                               site_that_owns_HC_grid_prim, &
                               herm_order_for_HC_grid_prim, &
                               herm_expon_for_HC_grid_prim, &
                               recip, &
                               volume, &
                               ewald_cosf, &
                               ewald_sinf, &
                               nfft1_for_grid_type(gt), &
                               nfft2_for_grid_type(gt), &
                               nfft3_for_grid_type(gt), &
                               nfftdim1_for_gridtype(gt), &
                               nfftdim2_for_gridtype(gt), &
                               nfftdim2_for_gridtype(gt), &
                               HC_grid(off+1), &
                               factor, &
                               energy, &
                               site_frc, &
                               Global_CHermite_field)
         endif !( struc_fac_method_for_grid_type(gt) == SF_EWALD )then
      enddo !g = 1,num_HC_prim_grids
   endif !( tot_num_herm_Cprims > 0 )then

   ! finally the diffuse hermites
   if ( tot_num_herm_Dprims > 0 )then
      do g = 1,num_HD_prim_grids
         off = HD_grid_offset_for_HD_grid(g)
         gt = grid_type_for_HD_prim_grid(g)
         if ( struc_fac_method_for_grid_type(gt) == SF_EWALD )then
            call FH_EWALD_RECIP_one_ene_frc(  &
                               num_sites, &
                               ewald_mmax, &
                               num_prims_for_HD_grid(g), &
                               prim_offset_for_HD_grid(g), &
                               Global_DHermite_coeff, &
                               Gcoeff_offset_for_HD_grid_prim, &
                               site_that_owns_HD_grid_prim, &
                               herm_order_for_HD_grid_prim, &
                               herm_expon_for_HD_grid_prim, &
                               recip, &
                               volume, &
                               ewald_cosf, &
                               ewald_sinf, &
                               nfft1_for_grid_type(gt), &
                               nfft2_for_grid_type(gt), &
                               nfft3_for_grid_type(gt), &
                               nfftdim1_for_gridtype(gt), &
                               nfftdim2_for_gridtype(gt), &
                               nfftdim2_for_gridtype(gt), &
                               HD_grid(off+1), &
                               factor, &
                               energy, &
                               site_frc, &
                               Global_DHermite_field)
         endif !( struc_fac_method_for_grid_type(gt) == SF_EWALD )then
      enddo !g = 1,num_HD_prim_grids
   endif !( tot_num_herm_Dprims > 0 )then

end subroutine FH_EWALD_RECIP_ene_force_field

!-----------------------------------------------------------------
!-----------------------------------------------------------------
! DEALLOCATE SUBROUTINES
!-----------------------------------------------------------------
!-----------------------------------------------------------------

!-----------------------------------------------------------------
subroutine FH_EWALD_RECIP_deallocate()

   implicit none

   if ( allocated(ewald_cosf) )deallocate(ewald_cosf)
   if ( allocated(ewald_sinf) )deallocate(ewald_sinf)
end subroutine FH_EWALD_RECIP_deallocate
!-----------------------------------------------------------------
end module ewald_recip
!-----------------------------------------------------------------
!-----------------------------------------------------------------
! NON-MODULE (LOCAL) SUBROUTINES
! THESE HANDLE ALL ARGUMENTS THROUGH ARGUMENT LIST
!-----------------------------------------------------------------
!-----------------------------------------------------------------
subroutine FH_EWALD_RECIP_get_mmax(nfft1,nfft2,nfft3,mmax)

   implicit none

   integer, intent(in) :: nfft1,nfft2,nfft3
   integer, intent(out) :: mmax

   integer :: m,nm,mm,nf1,nf2,nf3

   nf1 = nfft1/2
   if ( 2*nf1 < nfft1 )nf1 = nf1+1
   nf2 = nfft2/2
   if ( 2*nf2 < nfft2 )nf2 = nf2+1
   nf3 = nfft3/2
   if ( 2*nf3 < nfft3 )nf3 = nf3+1

   m = nf1 - 1
   nm = abs(nf1-nfft1)
   mm = max(m,nm)
   mmax = mm
   m = nf2 - 1
   nm = abs(nf2-nfft2)
   mm = max(m,nm)
   if ( mm > mmax )mmax = mm
   m = nf3 - 1
   nm = abs(nf3-nfft3)
   mm = max(m,nm)
   if ( mm > mmax )mmax = mm

end subroutine FH_EWALD_RECIP_get_mmax
!--------------------------------------------------------------------
subroutine FH_EWALD_RECIP_phase_factors(nsites,mmax,recip,site_crd,cosf,sinf)

   implicit none

   integer,intent(in) :: nsites,mmax
   double precision,intent(in) :: recip(3,3),site_crd(3,nsites)
   double precision,intent(out) :: cosf(3,nsites,0:mmax), &
                                   sinf(3,nsites,0:mmax)

   double precision :: frac(3,nsites)
   double precision :: pi,twopi,w
   integer :: j,n,m

   pi = 3.14159265358979323846d0
   twopi = 2.d0*3.14159265358979323846d0
   ! get frac coords
   do n = 1,nsites
      w = site_crd(1,n)*recip(1,1)+site_crd(2,n)*recip(2,1)+ &
          site_crd(3,n)*recip(3,1)
      frac(1,n) = w - dnint(w) + 0.5d0
      w = site_crd(1,n)*recip(1,2)+site_crd(2,n)*recip(2,2)+ &
          site_crd(3,n)*recip(3,2)
      frac(2,n) = w - dnint(w) + 0.5d0
      w = site_crd(1,n)*recip(1,3)+site_crd(2,n)*recip(2,3)+ &
          site_crd(3,n)*recip(3,3)
      frac(3,n) = w - dnint(w) + 0.5d0
   enddo
   ! next initial factors
   do n = 1,nsites
      do j = 1,3
         cosf(j,n,0) = 1.d0
         sinf(j,n,0) = 0.d0
      enddo
      do j = 1,3
         cosf(j,n,1) = cos(twopi*frac(j,n))
         sinf(j,n,1) = sin(twopi*frac(j,n))
      enddo
   enddo
   ! build the rest by recursion
   do m = 2,mmax
      do n = 1,nsites
         do j = 1,3
            cosf(j,n,m) = cosf(j,n,m-1)*cosf(j,n,1)-sinf(j,n,m-1)*sinf(j,n,1)
            sinf(j,n,m) = sinf(j,n,m-1)*cosf(j,n,1)+cosf(j,n,m-1)*sinf(j,n,1)
         enddo
      enddo
   enddo
end subroutine FH_EWALD_RECIP_phase_factors
!--------------------------------------------------------------------
subroutine FH_EWALD_RECIP_struc_fac( nsites,mmax, &
                               num_prims_for_grid, &
                               prim_offset_for_grid, &
                               Global_hermite_coeff, &
                               Gcoeff_offset_for_grid_prim, &
                               site_that_owns_grid_prim, &
                               hermite_order_for_grid_prim, &
                               hermite_expon_for_grid_prim, &
                               recip,volume,cosf,sinf, &
                               nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3, &
                               ewald_grid)

   implicit none

   integer,intent(in) :: nsites,mmax,num_prims_for_grid,prim_offset_for_grid
   double precision,intent(in) :: Global_hermite_coeff(*)
   integer,intent(in) :: Gcoeff_offset_for_grid_prim(*), &
                         site_that_owns_grid_prim(*), &
                         hermite_order_for_grid_prim(*)
   double precision,intent(in) :: hermite_expon_for_grid_prim(*), &
                                  recip(3,3),volume
   double precision,intent(in) :: cosf(3,nsites,0:mmax), &
                                  sinf(3,nsites,0:mmax)
   integer,intent(in) :: nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3
   double precision,intent(out) :: ewald_grid(2,nfft3,nfftdim1,nfft2)
                               
   include "direct_pointers.fh"
   integer :: kp,np,j,n,k1,k2,k3,m1,m2,m3,k10,nf1,nf2,nf3,t,u,v, &
              mp1,mp2,mp3,order,offset
   double precision :: pi,twopi,fac,mhat1,mhat2,mhat3,& 
                       msq,term,deriv_part(3,0:4),deriv(35), &
                       struc_fac(2),pre_fac(2),expon
   integer :: site_used(nsites)
   double precision :: c2(nsites),s2(nsites), &
                       c23(nsites),s23(nsites), &
                       c231(nsites),s231(nsites)

   pi = 3.14159265358979323846d0
   twopi = 2.d0*pi
   ! mark those sites covered by prims in this grid
   site_used = 0
   do kp = 1,num_prims_for_grid
      np = prim_offset_for_grid + kp
      n = site_that_owns_grid_prim(np)
      site_used(n) = 1
   enddo

   nf1 = nfft1/2
   if ( 2*nf1 < nfft1 )nf1 = nf1+1
   nf2 = nfft2/2
   if ( 2*nf2 < nfft2 )nf2 = nf2+1
   nf3 = nfft3/2
   if ( 2*nf3 < nfft3 )nf3 = nf3+1

   do k2 = 1, nfft2
      m2 = k2 - 1
      if ( k2 > nf2 )m2 = k2 - 1 - nfft2
      mp2 = abs(m2)
      if ( m2 < 0 )then ! we are using a back transform using exp(-2pim*u)
         do n = 1,nsites
            if ( site_used(n) == 1 )then
                c2(n) = cosf(2,n,mp2)
                s2(n) = sinf(2,n,mp2)! we are using exp(-2pim*u)
            endif
         enddo
      else
         do n = 1,nsites
            if ( site_used(n) == 1 )then
                c2(n) = cosf(2,n,mp2)
                s2(n) = -sinf(2,n,mp2)! we are using exp(-2pim*u)
            endif
         enddo
      endif
      do k3 = 1,nfft3
         m3 = k3 - 1
         if ( k3 > nf3 )m3 = k3 - 1 - nfft3
         mp3 = abs(m3)
         k10 = 1
         ! need (1,1,1) case also
         !if(k3+k2 == 2) k10 = 2
         ! get phase factors for marked sites
         if ( m3 < 0 )then! we are using a back transform using exp(-2pim*u)
            do n = 1,nsites
               if ( site_used(n) == 1 )then
                  c23(n) = c2(n)*cosf(3,n,mp3)-s2(n)*sinf(3,n,mp3)
                  s23(n) = s2(n)*cosf(3,n,mp3)+c2(n)*sinf(3,n,mp3)
               endif
            enddo
         else
            do n = 1,nsites
               if ( site_used(n) == 1 )then
                  c23(n) = c2(n)*cosf(3,n,mp3)+s2(n)*sinf(3,n,mp3)
                  s23(n) = s2(n)*cosf(3,n,mp3)-c2(n)*sinf(3,n,mp3)
               endif
            enddo
         endif
         do k1 = k10, nf1+1
            m1 = k1 - 1
            if ( k1 > nf1 )m1 = k1 - 1 - nfft1
            mp1 = abs(m1)
            if ( m1 < 0 )then ! we are using a back transform exp(-2pim*u)
               do n = 1,nsites
                  if ( site_used(n) == 1 )then
                     c231(n) = c23(n)*cosf(1,n,mp1)-s23(n)*sinf(1,n,mp1)
                     s231(n) = s23(n)*cosf(1,n,mp1)+c23(n)*sinf(1,n,mp1)
                  endif
               enddo
            else
               do n = 1,nsites
                  if ( site_used(n) == 1 )then
                     c231(n) = c23(n)*cosf(1,n,mp1)+s23(n)*sinf(1,n,mp1)
                     s231(n) = s23(n)*cosf(1,n,mp1)-c23(n)*sinf(1,n,mp1)
                  endif
               enddo
            endif
            ! we are using a back transform using exp(-2pim*u)
            mhat1 = -(recip(1,1)*m1+recip(1,2)*m2+recip(1,3)*m3)
            mhat2 = -(recip(2,1)*m1+recip(2,2)*m2+recip(2,3)*m3)
            mhat3 = -(recip(3,1)*m1+recip(3,2)*m2+recip(3,3)*m3)
            msq = mhat1*mhat1+mhat2*mhat2+mhat3*mhat3
            deriv_part(1,0) = 1.d0
            deriv_part(2,0) = 1.d0
            deriv_part(3,0) = 1.d0
            do j = 1,4
               deriv_part(1,j) = twopi*mhat1*deriv_part(1,j-1)
               deriv_part(2,j) = twopi*mhat2*deriv_part(2,j-1)
               deriv_part(3,j) = twopi*mhat3*deriv_part(3,j-1)
            enddo
            do j = 1,35
               t = t_part(j)
               u = u_part(j)
               v = v_part(j)
               deriv(j) = deriv_part(1,t)*deriv_part(2,u)*deriv_part(3,v)
            enddo
            struc_fac(1) = 0.d0
            struc_fac(2) = 0.d0
            do kp = 1,num_prims_for_grid
               np = prim_offset_for_grid + kp
               n = site_that_owns_grid_prim(np)
               order = hermite_order_for_grid_prim(np)
               offset = Gcoeff_offset_for_grid_prim(np)
               pre_fac(1) = 0.d0
               pre_fac(2) = 0.d0
               if ( order > 0 )then
                  pre_fac(1) = pre_fac(1) + Global_hermite_coeff(offset+1)
                  if ( order > 1 )then
                     do j = 2,4
                        pre_fac(2) = pre_fac(2) + &
                           deriv(j)*Global_hermite_coeff(offset+j)
                     enddo
                     if ( order > 4 )then
                        do j = 5,10
                           pre_fac(1) = pre_fac(1) - &
                              deriv(j)*Global_hermite_coeff(offset+j)
                        enddo
                        if ( order > 10 )then
                           do j = 11,20
                              pre_fac(2) = pre_fac(2) - &
                                 deriv(j)*Global_hermite_coeff(offset+j)
                           enddo
                           if ( order > 20 )then
                              do j = 21,35
                                 pre_fac(1) = pre_fac(1) + &
                                    deriv(j)*Global_hermite_coeff(offset+j)
                              enddo
                           endif !( order > 20 )
                        endif !( order > 10 )
                     endif !( order > 4 )
                  endif !( order > 1 )
               endif !( order > 0 )
               expon = hermite_expon_for_grid_prim(np)
               term = exp(-pi*pi*msq / expon) / volume
               struc_fac(1) = struc_fac(1) + &
                         term*(pre_fac(1)*c231(n) - pre_fac(2)*s231(n))
               struc_fac(2) = struc_fac(2) + &
                         term*(pre_fac(1)*s231(n) + pre_fac(2)*c231(n))
            enddo
            ewald_grid(1,k3,k1,k2) = struc_fac(1)
            ewald_grid(2,k3,k1,k2) = struc_fac(2)
         enddo
      enddo
   enddo
end subroutine FH_EWALD_RECIP_struc_fac
!--------------------------------------------------------------------
subroutine FH_EWALD_RECIP_one_ene_frc( nsites,mmax, &
                               num_prims_for_grid, &
                               prim_offset_for_grid, &
                               Global_hermite_coeff, &
                               Gcoeff_offset_for_grid_prim, &
                               site_that_owns_grid_prim, &
                               hermite_order_for_grid_prim, &
                               hermite_expon_for_grid_prim, &
                               recip,volume,cosf,sinf, &
                               nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3, &
                               ewald_grid,factor,energy,force, &
                               Global_hermite_field)

   implicit none

   integer,intent(in) :: nsites,mmax,num_prims_for_grid,prim_offset_for_grid
   double precision,intent(in) :: Global_hermite_coeff(*)
   integer,intent(in) :: Gcoeff_offset_for_grid_prim(*), &
                         site_that_owns_grid_prim(*), &
                         hermite_order_for_grid_prim(*)
   double precision,intent(in) :: hermite_expon_for_grid_prim(*), &
                                  recip(3,3),volume
   double precision,intent(in) :: cosf(3,nsites,0:mmax), &
                                  sinf(3,nsites,0:mmax)
   integer,intent(in) :: nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3
   double precision,intent(in) :: ewald_grid(2,nfft3,nfftdim1,nfft2), &
                                  factor
   double precision, intent(inout) :: energy,force(3,nsites)
   double precision,intent(inout) :: Global_hermite_field(*)
                               
   include "direct_pointers.fh"
   integer :: kp,np,j,j1,j2,j3,n,k1,k2,k3,m1,m2,m3,k10,nf1,nf2,nf3,t,u,v, &
              mp1,mp2,mp3,order,offset,field_order
   double precision :: pi,twopi,fac,mhat1,mhat2,mhat3,& 
                       msq,mult,term,deriv_part(3,0:5),deriv(56), &
                       expon,ene,postfac(2),field(56)
   integer :: site_used(nsites)
   double precision :: c2(nsites),s2(nsites), &
                       c23(nsites),s23(nsites), &
                       c231(nsites),s231(nsites)

   pi = 3.14159265358979323846d0
   twopi = 2.d0*3.14159265358979323846d0
   ene = 0.d0
   ! mark those sites covered by prims in this grid
   site_used = 0
   do kp = 1,num_prims_for_grid
      np = prim_offset_for_grid + kp
      n = site_that_owns_grid_prim(np)
      site_used(n) = 1
   enddo

   nf1 = nfft1/2
   if ( 2*nf1 < nfft1 )nf1 = nf1+1
   nf2 = nfft2/2
   if ( 2*nf2 < nfft2 )nf2 = nf2+1
   nf3 = nfft3/2
   if ( 2*nf3 < nfft3 )nf3 = nf3+1

   do k2 = 1, nfft2
      m2 = k2 - 1
      if ( k2 > nf2 )m2 = k2 - 1 - nfft2
      mp2 = abs(m2)
      if ( m2 < 0 )then ! we are using a forward transform using exp(2pim*u)
         do n = 1,nsites
            if ( site_used(n) == 1 )then
                c2(n) = cosf(2,n,mp2)
                s2(n) = -sinf(2,n,mp2)! we are using exp(2pim*u)
            endif
         enddo
      else
         do n = 1,nsites
            if ( site_used(n) == 1 )then
                c2(n) = cosf(2,n,mp2)
                s2(n) = sinf(2,n,mp2)! we are using exp(2pim*u)
            endif
         enddo
      endif
      do k3 = 1,nfft3
         m3 = k3 - 1
         if ( k3 > nf3 )m3 = k3 - 1 - nfft3
         mp3 = abs(m3)
         k10 = 1
         ! need (1,1,1) case also
         !if(k3+k2 == 2) k10 = 2
         ! get phase factors for marked sites
         if ( m3 < 0 )then! we are using a forward transform using exp(2pim*u)
            do n = 1,nsites
               if ( site_used(n) == 1 )then
                  c23(n) = c2(n)*cosf(3,n,mp3)+s2(n)*sinf(3,n,mp3)
                  s23(n) = s2(n)*cosf(3,n,mp3)-c2(n)*sinf(3,n,mp3)
               endif
            enddo
         else
            do n = 1,nsites
               if ( site_used(n) == 1 )then
                  c23(n) = c2(n)*cosf(3,n,mp3)-s2(n)*sinf(3,n,mp3)
                  s23(n) = s2(n)*cosf(3,n,mp3)+c2(n)*sinf(3,n,mp3)
               endif
            enddo
         endif
         do k1 = k10, nf1+1
            if ( k1 > 1 )then
               mult = 2.d0
            else
               mult = 1.d0
            endif
            m1 = k1 - 1
            if ( k1 > nf1 )m1 = k1 - 1 - nfft1
            mp1 = abs(m1)
            if ( m1 < 0 )then ! we are using a forward transform exp(2pim*u)
               do n = 1,nsites
                  if ( site_used(n) == 1 )then
                     c231(n) = c23(n)*cosf(1,n,mp1)+s23(n)*sinf(1,n,mp1)
                     s231(n) = s23(n)*cosf(1,n,mp1)-c23(n)*sinf(1,n,mp1)
                  endif
               enddo
            else
               do n = 1,nsites
                  if ( site_used(n) == 1 )then
                     c231(n) = c23(n)*cosf(1,n,mp1)-s23(n)*sinf(1,n,mp1)
                     s231(n) = s23(n)*cosf(1,n,mp1)+c23(n)*sinf(1,n,mp1)
                  endif
               enddo
            endif
            ! we are using a forward transform using exp(2pim*u)
            mhat1 = (recip(1,1)*m1+recip(1,2)*m2+recip(1,3)*m3)
            mhat2 = (recip(2,1)*m1+recip(2,2)*m2+recip(2,3)*m3)
            mhat3 = (recip(3,1)*m1+recip(3,2)*m2+recip(3,3)*m3)
            msq = mhat1*mhat1+mhat2*mhat2+mhat3*mhat3
            deriv_part(1,0) = 1.d0
            deriv_part(2,0) = 1.d0
            deriv_part(3,0) = 1.d0
            do j = 1,5
               deriv_part(1,j) = twopi*mhat1*deriv_part(1,j-1)
               deriv_part(2,j) = twopi*mhat2*deriv_part(2,j-1)
               deriv_part(3,j) = twopi*mhat3*deriv_part(3,j-1)
            enddo
            do j = 1,56
               t = t_part(j)
               u = u_part(j)
               v = v_part(j)
               deriv(j) = deriv_part(1,t)*deriv_part(2,u)*deriv_part(3,v)
            enddo
            do kp = 1,num_prims_for_grid
               np = prim_offset_for_grid + kp
               n = site_that_owns_grid_prim(np)
               expon = hermite_expon_for_grid_prim(np)
               term = mult*exp(-pi*pi*msq / expon) / volume
               postfac(1) = term*( c231(n)*ewald_grid(1,k3,k1,k2) - &
                                   s231(n)*ewald_grid(2,k3,k1,k2) )
               postfac(2) = term*( c231(n)*ewald_grid(2,k3,k1,k2) + &
                                   s231(n)*ewald_grid(1,k3,k1,k2) )
               order = hermite_order_for_grid_prim(np)
               field_order = 0
               if ( order == 1 )then
                  field_order = 4
               elseif ( order == 4 )then
                  field_order = 10
               elseif ( order == 10 )then
                  field_order = 20
               elseif ( order == 20 )then
                  field_order = 35
               elseif ( order == 35 )then
                  field_order = 56
               endif 
               offset = Gcoeff_offset_for_grid_prim(np)
               if ( field_order > 0 )then
                  field(1) = postfac(1)
                  if ( field_order > 1 )then
                     ! Real part of 2pi.i.m * (postfac(1)+i.postfac(2))
                     do j = 2,4 
                        field(j) = -deriv(j)*postfac(2)
                     enddo
                     if ( field_order > 4 )then
                        ! Real part of (2pi.i.m)**2 * (postfac(1)+i.postfac(2))
                        do j = 5,10
                           field(j) = -deriv(j)*postfac(1)
                        enddo
                        if ( field_order > 10 )then
                           ! Real part of (2pi.i.m)**3 * (postfac(1)+i.postfac(2))
                           do j = 11,20
                              field(j) = +deriv(j)*postfac(2)
                           enddo
                           if ( field_order > 20 )then
                              ! Real part of (2pi.i.m)**4 * (postfac(1)+i.postfac(2))
                              do j = 21,35
                                 field(j) = +deriv(j)*postfac(1)
                              enddo
                              if ( field_order > 35 )then
                                 do j = 36,56
                                    field(j) = -deriv(j)*postfac(2)
                                 enddo
                              endif !( field_order > 35 )
                           endif !( field_order > 20 )
                        endif !( field_order > 10 )
                     endif !( field_order > 4 )
                  endif !( field_order > 1 )
               endif !( field_order > 0 )
               do j = 1,order
                  ene = ene + factor*Global_hermite_coeff(offset+j)*field(j)
                  Global_hermite_field(offset+j) = &
                     Global_hermite_field(offset+j) + factor*field(j)
                  j1 = p_xder(j)
                  j2 = p_yder(j)
                  j3 = p_zder(j)
                  force(1,n) = force(1,n) - &
                               factor*Global_hermite_coeff(offset+j)*field(j1)
                  force(2,n) = force(2,n) - &
                               factor*Global_hermite_coeff(offset+j)*field(j2)
                  force(3,n) = force(3,n) - &
                               factor*Global_hermite_coeff(offset+j)*field(j3)
               enddo
            enddo !kp = 1,num_prims_for_grid
         enddo !do k1 = k10, nf1+1
      enddo !do k3 = 1,nfft3
   enddo !do k2 = 1, nfft2
   energy = energy + 0.5d0*ene
end subroutine FH_EWALD_RECIP_one_ene_frc
!--------------------------------------------------------------------
