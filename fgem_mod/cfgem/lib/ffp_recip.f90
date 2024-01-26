module ffp_recip
   implicit none
   private

!-----------------------------------------------------------------
   public         FH_FFP_RECIP_fill_grids, &
                  FH_FFP_RECIP_FT_density, &
                  FH_FFP_RECIP_fix_FT_phi, &
                  FH_FFP_RECIP_ene_force_field
contains
!-----------------------------------------------------------------
!-----------------------------------------------------------------
! EVALUATION ROUTINES---------------------------------------------
!-----------------------------------------------------------------
!-----------------------------------------------------------------
subroutine FH_FFP_RECIP_fill_grids(out_lun)

   use user,only : num_HC_prim_grids, &
                   num_HD_prim_grids, &
                   grid_type_for_HC_prim_grid, &
                   grid_type_for_HD_prim_grid, &
                   grid_type_for_MPOLES, &
                   grid_type_for_sumCH, &
                   struc_fac_method_for_grid_type, &
                   nfft1_for_grid_type, &
                   nfft2_for_grid_type, &
                   nfft3_for_grid_type, &
                   gaussian_recip_tol, &
                   num_ffp_grid_types
   use unit_cell, only : recip,ucell
   use sites, only : num_sites,site_crd
   use hermite, only : &
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
                       ! compact hermite
                       tot_num_herm_Cprims, &
                       Global_Chermite_coeff, &
                       ! diffuse hermite
                       tot_num_herm_Dprims, &
                       Global_Dhermite_coeff
   use recip_grids, only :  &
                           ! mpoles
                           num_mpoles_for_mpole_grid, &
                           site_that_owns_gridded_mpole, &
                           Mpole_order_for_gridded_mpole, &
                           Gcoeff_offset_for_gridded_mpole, &
                           MC_grid, &
                           MD_grid, &
                           ! sumCH
                           num_sumch_for_sumch_grid, &
                           site_that_owns_gridded_sumch, &
                           sumch_order_for_gridded_sumch, &
                           Gcoeff_offset_for_gridded_sumch, &
                           sumCH_C_grid, &
                           sumCH_D_grid, &
                           ! compact hermites
                           num_prims_for_HC_grid, &
                           prim_offset_for_HC_grid, &
                           site_that_owns_HC_grid_prim, &
                           herm_order_for_HC_grid_prim, &
                           herm_expon_for_HC_grid_prim, &
                           Gcoeff_offset_for_HC_grid_prim, &
                           HC_grid_offset_for_HC_grid,HC_grid, &
                           ! diffuse hermites
                           num_prims_for_HD_grid, &
                           prim_offset_for_HD_grid, &
                           site_that_owns_HD_grid_prim, &
                           herm_order_for_HD_grid_prim, &
                           herm_expon_for_HD_grid_prim, &
                           Gcoeff_offset_for_HD_grid_prim, &
                           HD_grid_offset_for_HD_grid,HD_grid, &
                           ! general
                           nfftdim1_for_gridtype, &
                           nfftdim2_for_gridtype, &
                           nfftdim3_for_gridtype

   implicit none

   integer, intent(in)  :: out_lun

   include 'structure_factor_type.fh'

   integer :: n,g,gt,off,off1
   integer :: ind_map(3,num_sites)
   double precision :: del_crd(3,num_sites),invtol,gauss_cut
   double precision :: mpole_exponent(num_mpoles_for_mpole_grid), &
                       sumCH_exponent(num_sumCH_for_sumCH_grid)

   if ( num_ffp_grid_types == 0 )return ! nothing further to do here
   invtol = 1.d0 / gaussian_recip_tol
   gauss_cut = sqrt(dlog(invtol))

   ! first fill the compact and diffuse multipole grid 
   if ( tot_num_mpole_coeffs > 0 )then
      gt = grid_type_for_MPOLES
      if ( struc_fac_method_for_grid_type(gt) == SF_FFP )then
         off = 0 ! just one grid
         off1 = 0 ! mpole offset is also 0, since one grid
         ! multipole exponents are diffuse adjusted to grid them
         do n = 1,num_mpoles_for_mpole_grid
            mpole_exponent(n) = compact_rec_mpole_expon
         enddo
         call FH_FFP_RECIP_map_coords( &
                                 num_sites, &
                                 nfft1_for_grid_type(gt), &
                                 nfft2_for_grid_type(gt), &
                                 nfft3_for_grid_type(gt), &
                                 recip, &
                                 ucell, &
                                 site_crd, &
                                 ind_map, &
                                 del_crd)
         call FH_FFP_RECIP_grid_GHerm_coeffs( &
                                 num_sites, &
                                 num_mpoles_for_mpole_grid, &
                                 off1, &
                                 Global_multipole_coeff, &
                                 site_that_owns_gridded_mpole, &
                                 ind_map, &
                                 del_crd, &
                                 gauss_cut, &
                                 ucell, &
                                 Gcoeff_offset_for_gridded_mpole, &
                                 Mpole_order_for_gridded_mpole, &
                                 mpole_exponent, &
                                 nfft1_for_grid_type(gt), &
                                 nfft2_for_grid_type(gt), &
                                 nfft3_for_grid_type(gt), &
                                 nfftdim1_for_gridtype(gt), &
                                 nfftdim2_for_gridtype(gt), &
                                 nfftdim3_for_gridtype(gt), &
                                 MC_grid(off+1),out_lun)
         ! next the diffuse
         do n = 1,num_mpoles_for_mpole_grid
            mpole_exponent(n) = diffuse_rec_mpole_expon
         enddo
         call FH_FFP_RECIP_grid_GHerm_coeffs( &
                                 num_sites, &
                                 num_mpoles_for_mpole_grid, &
                                 off1, &
                                 Global_multipole_coeff, &
                                 site_that_owns_gridded_mpole, &
                                 ind_map, &
                                 del_crd, &
                                 gauss_cut, &
                                 ucell, &
                                 Gcoeff_offset_for_gridded_mpole, &
                                 Mpole_order_for_gridded_mpole, &
                                 mpole_exponent, &
                                 nfft1_for_grid_type(gt), &
                                 nfft2_for_grid_type(gt), &
                                 nfft3_for_grid_type(gt), &
                                 nfftdim1_for_gridtype(gt), &
                                 nfftdim2_for_gridtype(gt), &
                                 nfftdim3_for_gridtype(gt), &
                                 MD_grid(off+1),out_lun)
      endif !( struc_fac_method_for_grid_type(gt) == SF_FFP )
   endif !( tot_num_mpole_coeffs > 0 )then

   ! first fill the compact and diffuse sumCH grid 
   if ( tot_num_sumCH_coeffs > 0 )then
      gt = grid_type_for_MPOLES
      if ( struc_fac_method_for_grid_type(gt) == SF_FFP )then
         off = 0 ! just one grid
         off1 = 0 ! sumCH offset is also 0, since one grid
         ! multipole exponents are diffuse adjusted to grid them
         do n = 1,num_sumCH_for_sumCH_grid
            sumCH_exponent(n) = compact_rec_sumCH_expon
         enddo
         call FH_FFP_RECIP_map_coords( &
                                 num_sites, &
                                 nfft1_for_grid_type(gt), &
                                 nfft2_for_grid_type(gt), &
                                 nfft3_for_grid_type(gt), &
                                 recip, &
                                 ucell, &
                                 site_crd, &
                                 ind_map, &
                                 del_crd)
         call FH_FFP_RECIP_grid_GHerm_coeffs( &
                                 num_sites, &
                                 num_sumCH_for_sumCH_grid, &
                                 off1, &
                                 Global_sumCH_coeff, &
                                 site_that_owns_gridded_sumCH, &
                                 ind_map, &
                                 del_crd, &
                                 gauss_cut, &
                                 ucell, &
                                 Gcoeff_offset_for_gridded_sumCH, &
                                 sumCH_order_for_gridded_sumCH, &
                                 sumCH_exponent, &
                                 nfft1_for_grid_type(gt), &
                                 nfft2_for_grid_type(gt), &
                                 nfft3_for_grid_type(gt), &
                                 nfftdim1_for_gridtype(gt), &
                                 nfftdim2_for_gridtype(gt), &
                                 nfftdim3_for_gridtype(gt), &
                                 sumCH_C_grid(off+1),out_lun)
         ! next the diffuse
         do n = 1,num_sumCH_for_sumCH_grid
            sumCH_exponent(n) = diffuse_rec_sumCH_expon
         enddo
         call FH_FFP_RECIP_grid_GHerm_coeffs( &
                                 num_sites, &
                                 num_sumCH_for_sumCH_grid, &
                                 off1, &
                                 Global_sumCH_coeff, &
                                 site_that_owns_gridded_sumCH, &
                                 ind_map, &
                                 del_crd, &
                                 gauss_cut, &
                                 ucell, &
                                 Gcoeff_offset_for_gridded_sumCH, &
                                 sumCH_order_for_gridded_sumCH, &
                                 sumCH_exponent, &
                                 nfft1_for_grid_type(gt), &
                                 nfft2_for_grid_type(gt), &
                                 nfft3_for_grid_type(gt), &
                                 nfftdim1_for_gridtype(gt), &
                                 nfftdim2_for_gridtype(gt), &
                                 nfftdim3_for_gridtype(gt), &
                                 sumCH_D_grid(off+1),out_lun)
      endif !( struc_fac_method_for_grid_type(gt) == SF_FFP )
   endif !( tot_num_sumCH_coeffs > 0 )then

   ! next the compact hermites
   if ( tot_num_herm_Cprims > 0 )then
      do g = 1,num_HC_prim_grids
         off = HC_grid_offset_for_HC_grid(g)
         gt = grid_type_for_HC_prim_grid(g)
         if ( struc_fac_method_for_grid_type(gt) == SF_FFP )then
            call FH_FFP_RECIP_map_coords( &
                                 num_sites, &
                                 nfft1_for_grid_type(gt), &
                                 nfft2_for_grid_type(gt), &
                                 nfft3_for_grid_type(gt), &
                                 recip,ucell,site_crd, &
                                 ind_map, &
                                 del_crd)
            call FH_FFP_RECIP_grid_GHerm_coeffs( &
                                 num_sites, &
                                 num_prims_for_HC_grid(g), &
                                 prim_offset_for_HC_grid(g), &
                                 Global_CHermite_coeff, &
                                 site_that_owns_HC_grid_prim, &
                                 ind_map, &
                                 del_crd, &
                                 gauss_cut, &
                                 ucell, &
                                 Gcoeff_offset_for_HC_grid_prim, &
                                 herm_order_for_HC_grid_prim, &
                                 herm_expon_for_HC_grid_prim, &
                                 nfft1_for_grid_type(gt), &
                                 nfft2_for_grid_type(gt), &
                                 nfft3_for_grid_type(gt), &
                                 nfftdim1_for_gridtype(gt), &
                                 nfftdim2_for_gridtype(gt), &
                                 nfftdim3_for_gridtype(gt), &
                                 HC_grid(off+1),out_lun)
         endif !( struc_fac_method_for_grid_type(gt) == SF_FFP )
      enddo !g = 1,num_HC_prim_grids
   endif !( tot_num_herm_Cprims > 0 )then

   ! finally the diffuse hermites
   if ( tot_num_herm_Dprims > 0 )then
      do g = 1,num_HD_prim_grids
         off = HD_grid_offset_for_HD_grid(g)
         gt = grid_type_for_HD_prim_grid(g)
         if ( struc_fac_method_for_grid_type(gt) == SF_FFP )then
            call FH_FFP_RECIP_map_coords( &
                                 num_sites, &
                                 nfft1_for_grid_type(gt), &
                                 nfft2_for_grid_type(gt), &
                                 nfft3_for_grid_type(gt), &
                                 recip,ucell,site_crd, &
                                 ind_map, &
                                 del_crd)
            call FH_FFP_RECIP_grid_GHerm_coeffs( &
                                 num_sites, &
                                 num_prims_for_HD_grid(g), &
                                 prim_offset_for_HD_grid(g), &
                                 Global_DHermite_coeff, &
                                 site_that_owns_HD_grid_prim, &
                                 ind_map, &
                                 del_crd, &
                                 gauss_cut, &
                                 ucell, &
                                 Gcoeff_offset_for_HD_grid_prim, &
                                 herm_order_for_HD_grid_prim, &
                                 herm_expon_for_HD_grid_prim, &
                                 nfft1_for_grid_type(gt), &
                                 nfft2_for_grid_type(gt), &
                                 nfft3_for_grid_type(gt), &
                                 nfftdim1_for_gridtype(gt), &
                                 nfftdim2_for_gridtype(gt), &
                                 nfftdim3_for_gridtype(gt), &
                                 HD_grid(off+1),out_lun)
         endif !( struc_fac_method_for_grid_type(gt) == SF_FFP )
      enddo !g = 1,num_HD_prim_grids
   endif !( tot_num_herm_Dprims > 0 )then

end subroutine FH_FFP_RECIP_fill_grids
!-----------------------------------------------------------------
subroutine FH_FFP_RECIP_FT_density()
   use user,only : num_HC_prim_grids, &
                   num_HD_prim_grids, &
                   grid_type_for_HC_prim_grid, &
                   grid_type_for_HD_prim_grid, &
                   grid_type_for_MPOLES, &
                   grid_type_for_sumCH, &
                   struc_fac_method_for_grid_type, &
                   nfft1_for_grid_type, &
                   nfft2_for_grid_type, &
                   nfft3_for_grid_type, &
                   num_ffp_grid_types
   use unit_cell, only : volume
   use recip_grids, only : HC_grid_offset_for_HC_grid, &
                           HC_grid, &
                           HD_grid_offset_for_HD_grid, &
                           HD_grid, &
                           nfftdim1_for_gridtype, &
                           MC_grid, &
                           MD_grid, &
                           sumCH_C_grid, &
                           sumCH_D_grid
   use hermite, only : tot_num_mpole_coeffs, &
                       tot_num_sumCH_coeffs, &
                       tot_num_herm_Cprims, &
                       tot_num_herm_Dprims

   implicit none

   include 'structure_factor_type.fh'

   integer :: g,gt,off

   if ( num_ffp_grid_types == 0 )return ! nothing further to do here

   ! first the compact and diffuse multipole grids
   if ( tot_num_mpole_coeffs > 0 )then
      gt = grid_type_for_MPOLES
      if ( struc_fac_method_for_grid_type(gt) == SF_FFP )then
         off = 0 ! only one grid
         call FH_FFP_RECIP_one_FT_density( &
            nfft1_for_grid_type(gt), &
            nfft2_for_grid_type(gt), &
            nfft3_for_grid_type(gt), &
            nfftdim1_for_gridtype(gt), &
            volume, &
            MC_grid(off+1))
         call FH_FFP_RECIP_one_FT_density( &
            nfft1_for_grid_type(gt), &
            nfft2_for_grid_type(gt), &
            nfft3_for_grid_type(gt), &
            nfftdim1_for_gridtype(gt), &
            volume, &
            MD_grid(off+1))
      endif !( struc_fac_method_for_grid_type(gt) == SF_FFP )then
   endif !( tot_num_mpole_coeffs > 0 )then

   ! next the compact and diffuse sumCH grids
   if ( tot_num_sumCH_coeffs > 0 )then
      gt = grid_type_for_sumCH
      if ( struc_fac_method_for_grid_type(gt) == SF_FFP )then
         off = 0 ! only one grid
         call FH_FFP_RECIP_one_FT_density( &
            nfft1_for_grid_type(gt), &
            nfft2_for_grid_type(gt), &
            nfft3_for_grid_type(gt), &
            nfftdim1_for_gridtype(gt), &
            volume, &
            sumCH_C_grid(off+1))
         call FH_FFP_RECIP_one_FT_density( &
            nfft1_for_grid_type(gt), &
            nfft2_for_grid_type(gt), &
            nfft3_for_grid_type(gt), &
            nfftdim1_for_gridtype(gt), &
            volume, &
            sumCH_D_grid(off+1))
      endif !( struc_fac_method_for_grid_type(gt) == SF_FFP )then
   endif !( tot_num_sumCH_coeffs > 0 )then

   ! next the compact hermite grids
   if ( tot_num_herm_Cprims > 0 )then
      do g = 1,num_HC_prim_grids
         gt = grid_type_for_HC_prim_grid(g)
         if ( struc_fac_method_for_grid_type(gt) == SF_FFP )then
            off = HC_grid_offset_for_HC_grid(g)
            call FH_FFP_RECIP_one_FT_density( &
               nfft1_for_grid_type(gt), &
               nfft2_for_grid_type(gt), &
               nfft3_for_grid_type(gt), &
               nfftdim1_for_gridtype(gt), &
               volume, &
               HC_grid(off+1))
         endif !( struc_fac_method_for_grid_type(gt) == SF_FFP )then
      enddo !g = 1,num_HC_prim_grids
   endif !( tot_num_herm_Cprims > 0 )then

   ! finally the diffuse hermite grids
   if ( tot_num_herm_Dprims > 0 )then
      do g = 1,num_HD_prim_grids
         gt = grid_type_for_HD_prim_grid(g)
         if ( struc_fac_method_for_grid_type(gt) == SF_FFP )then
            off = HD_grid_offset_for_HD_grid(g)
            call FH_FFP_RECIP_one_FT_density( &
               nfft1_for_grid_type(gt), &
               nfft2_for_grid_type(gt), &
               nfft3_for_grid_type(gt), &
               nfftdim1_for_gridtype(gt), &
               volume, &
               HD_grid(off+1))
         endif !( struc_fac_method_for_grid_type(gt) == SF_FFP )then
      enddo !g = 1,num_HD_prim_grids
   endif !( tot_num_herm_Dprims > 0 )then

end subroutine FH_FFP_RECIP_FT_density
!-----------------------------------------------------------------
subroutine FH_FFP_RECIP_fix_FT_phi()

   use user,only : num_HC_prim_grids, &
                   num_HD_prim_grids, &
                   grid_type_for_HC_prim_grid, &
                   grid_type_for_HD_prim_grid, &
                   grid_type_for_MPOLES, &
                   grid_type_for_sumCH, &
                   struc_fac_method_for_grid_type, &
                   nfft1_for_grid_type, &
                   nfft2_for_grid_type, &
                   nfft3_for_grid_type, &
                   num_ffp_grid_types
   use unit_cell, only : volume
   use recip_grids, only : HC_grid_offset_for_HC_grid, &
                           HC_grid, &
                           HD_grid_offset_for_HD_grid, &
                           HD_grid, &
                           nfftdim1_for_gridtype, &
                           MC_grid, &
                           MD_grid, &
                           sumCH_C_grid, &
                           sumCH_D_grid
   use hermite, only : tot_num_mpole_coeffs, &
                       tot_num_sumCH_coeffs, &
                       tot_num_herm_Cprims, &
                       tot_num_herm_Dprims

   implicit none

   include 'structure_factor_type.fh'

   integer :: g,gt,off

   if ( num_ffp_grid_types == 0 )return ! nothing further to do here

   ! first the compact and diffuse multipole grids
   if ( tot_num_mpole_coeffs > 0 )then
      gt = grid_type_for_MPOLES
      if ( struc_fac_method_for_grid_type(gt) == SF_FFP )then
         off = 0 ! only one grid
         call FH_FFP_RECIP_one_FT_phi( &
               nfft1_for_grid_type(gt), &
               nfft2_for_grid_type(gt), &
               nfft3_for_grid_type(gt), &
               nfftdim1_for_gridtype(gt), &
               volume, &
               MC_grid(off+1))
         call FH_FFP_RECIP_one_FT_phi( &
               nfft1_for_grid_type(gt), &
               nfft2_for_grid_type(gt), &
               nfft3_for_grid_type(gt), &
               nfftdim1_for_gridtype(gt), &
               volume, &
               MD_grid(off+1))
      endif !( struc_fac_method_for_grid_type(gt) == SF_FFP )then
   endif !( tot_num_mpole_coeffs > 0 )then

   ! next the compact and diffuse sumCH grids
   if ( tot_num_sumCH_coeffs > 0 )then
      gt = grid_type_for_sumCH
      if ( struc_fac_method_for_grid_type(gt) == SF_FFP )then
         off = 0 ! only one grid
         call FH_FFP_RECIP_one_FT_phi( &
               nfft1_for_grid_type(gt), &
               nfft2_for_grid_type(gt), &
               nfft3_for_grid_type(gt), &
               nfftdim1_for_gridtype(gt), &
               volume, &
               sumCH_C_grid(off+1))
         call FH_FFP_RECIP_one_FT_phi( &
               nfft1_for_grid_type(gt), &
               nfft2_for_grid_type(gt), &
               nfft3_for_grid_type(gt), &
               nfftdim1_for_gridtype(gt), &
               volume, &
               sumCH_D_grid(off+1))
      endif !( struc_fac_method_for_grid_type(gt) == SF_FFP )then
   endif !( tot_num_sumCH_coeffs > 0 )then

   ! next the compact hermite grids
   if ( tot_num_herm_Cprims > 0 )then
      do g = 1,num_HC_prim_grids
         gt = grid_type_for_HC_prim_grid(g)
         if ( struc_fac_method_for_grid_type(gt) == SF_FFP )then
            off = HC_grid_offset_for_HC_grid(g)
            call FH_FFP_RECIP_one_FT_phi( &
               nfft1_for_grid_type(gt), &
               nfft2_for_grid_type(gt), &
               nfft3_for_grid_type(gt), &
               nfftdim1_for_gridtype(gt), &
               volume, &
               HC_grid(off+1))
         endif !( struc_fac_method_for_grid_type(gt) == SF_FFP )then
      enddo !g = 1,num_HC_prim_grids
   endif !( tot_num_herm_Cprims > 0 )then

   ! finally the diffuse hermite grids
   if ( tot_num_herm_Dprims > 0 )then
      do g = 1,num_HD_prim_grids
         gt = grid_type_for_HD_prim_grid(g)
         if ( struc_fac_method_for_grid_type(gt) == SF_FFP )then
            off = HD_grid_offset_for_HD_grid(g)
            call FH_FFP_RECIP_one_FT_phi( &
               nfft1_for_grid_type(gt), &
               nfft2_for_grid_type(gt), &
               nfft3_for_grid_type(gt), &
               nfftdim1_for_gridtype(gt), &
               volume, &
               HD_grid(off+1))
         endif !( struc_fac_method_for_grid_type(gt) == SF_FFP )then
      enddo !g = 1,num_HD_prim_grids
   endif !( tot_num_herm_Dprims > 0 )then

end subroutine FH_FFP_RECIP_fix_FT_phi
!-----------------------------------------------------------------
subroutine FH_FFP_RECIP_ene_force_field(energy, out_lun)

   use unit_cell, only : recip,ucell
   use sites, only : num_sites,site_crd,site_frc
   use user,only : num_HC_prim_grids, &
                   num_HD_prim_grids, &
                   grid_type_for_HC_prim_grid, &
                   grid_type_for_HD_prim_grid, &
                   grid_type_for_MPOLES, &
                   grid_type_for_sumCH, &
                   struc_fac_method_for_grid_type, &
                   nfft1_for_grid_type, &
                   nfft2_for_grid_type, &
                   nfft3_for_grid_type, &
                   gaussian_recip_tol, &
                   num_ffp_grid_types
   use hermite, only : &
                       ! mpole
                       tot_num_mpole_coeffs, &
                       Global_multipole_field, &
                       Global_multipole_coeff, &
                       compact_rec_mpole_expon, &
                       diffuse_rec_mpole_expon, &
                       ! sumCH
                       tot_num_sumCH_coeffs, &
                       Global_sumCH_field, &
                       Global_sumCH_coeff, &
                       compact_rec_sumCH_expon, &
                       diffuse_rec_sumCH_expon, &
                       !compact hermites
                       tot_num_herm_Cprims, &
                       Global_CHermite_coeff, &
                       Global_CHermite_field, &
                       !diffuse hermites
                       tot_num_herm_Dprims, &
                       Global_DHermite_coeff, &
                       Global_DHermite_field

   use recip_grids, only :  &
                           ! compact hermites
                           num_prims_for_HC_grid, &
                           prim_offset_for_HC_grid, &
                           site_that_owns_HC_grid_prim, &
                           Ffield_offset_for_HC_grid_prim, &
                           Ffield_order_for_HC_grid_prim, &
                           herm_order_for_HC_grid_prim, &
                           herm_expon_for_HC_grid_prim, &
                           Gcoeff_offset_for_HC_grid_prim, &
                           HC_grid_offset_for_HC_grid, &
                           HC_grid, &
                           HC_grid_Ffield, &
                           ! diffuse hermites
                           num_prims_for_HD_grid, &
                           prim_offset_for_HD_grid, &
                           site_that_owns_HD_grid_prim, &
                           Ffield_offset_for_HD_grid_prim, &
                           Ffield_order_for_HD_grid_prim, &
                           herm_order_for_HD_grid_prim, &
                           herm_expon_for_HD_grid_prim, &
                           Gcoeff_offset_for_HD_grid_prim, &
                           HD_grid_offset_for_HD_grid, &
                           HD_grid, &
                           HD_grid_Ffield, &
                           ! mpoles
                           num_mpoles_for_mpole_grid, &
                           site_that_owns_gridded_mpole, &
                           Ffield_offset_for_gridded_mpole, &
                           Ffield_order_for_gridded_mpole, &
                           Mpole_order_for_gridded_mpole, &
                           Gcoeff_offset_for_gridded_mpole, &
                           MC_grid, &
                           Mpole_grid_Ffield, &
                           ! sumCH
                           num_sumCH_for_sumCH_grid, &
                           site_that_owns_gridded_sumCH, &
                           Ffield_offset_for_gridded_sumCH, &
                           Ffield_order_for_gridded_sumCH, &
                           sumCH_order_for_gridded_sumCH, &
                           Gcoeff_offset_for_gridded_sumCH, &
                           sumCH_C_grid, &
                           sumCH_grid_Ffield, &
                           ! general
                           nfftdim1_for_gridtype, &
                           nfftdim2_for_gridtype, &
                           nfftdim3_for_gridtype

   implicit none

   double precision,intent(inout) :: energy 
   integer, intent(in)  :: out_lun

   include 'structure_factor_type.fh'

   integer :: n,g,gt,off,off1
   integer :: ind_map(3,num_sites)
   double precision :: del_crd(3,num_sites),invtol,gauss_cut
   double precision :: mpole_exponent(num_mpoles_for_mpole_grid), &
                       sumCH_exponent(num_sumCH_for_sumCH_grid)
   double precision :: factor
   !if ( energy_type == OVERLAP_ene_type )then
      !factor = exchange_factor
   !else
      factor = 1.d0
   !endif


   if ( num_ffp_grid_types == 0 )return ! nothing further to do here
   invtol = 1.d0 / gaussian_recip_tol
   gauss_cut = sqrt(dlog(invtol))

   ! first the compact multipoles
   if ( tot_num_mpole_coeffs > 0 )then
      gt = grid_type_for_MPOLES
      if ( struc_fac_method_for_grid_type(gt) == SF_FFP )then
         off = 0 ! just one grid
         off1 = 0 ! mpole offset is also 0, since one grid
         ! multipole exponents are diffuse adjusted to grid them
         do n = 1,num_mpoles_for_mpole_grid
            mpole_exponent(n) = compact_rec_mpole_expon
         enddo
         call FH_FFP_RECIP_map_coords( &
                                 num_sites, &
                                 nfft1_for_grid_type(gt), &
                                 nfft2_for_grid_type(gt), &
                                 nfft3_for_grid_type(gt), &
                                 recip, &
                                 ucell, &
                                 site_crd, &
                                 ind_map, &
                                 del_crd)
         call FH_FFP_RECIP_get_GHerm_phi( &
                                 num_sites, &
                                 num_mpoles_for_mpole_grid, &
                                 off1, &
                                 site_that_owns_gridded_mpole, &
                                 ind_map, &
                                 del_crd, &
                                 gauss_cut, &
                                 ucell, &
                                 ! note use Ffield arrays here!!(same size)
                                 Ffield_offset_for_gridded_mpole, &
                                 Ffield_order_for_gridded_mpole, &
                                 mpole_exponent, &
                                 nfft1_for_grid_type(gt), &
                                 nfft2_for_grid_type(gt), &
                                 nfft3_for_grid_type(gt), &
                                 nfftdim1_for_gridtype(gt), &
                                 nfftdim2_for_gridtype(gt), &
                                 nfftdim3_for_gridtype(gt), &
                                 !note use Mpole_grid_Ffield here (same size)
                                 MC_grid(off+1), &
                                 Mpole_grid_Ffield,out_lun )
         call FH_FFP_RECIP_one_ene_frc_field( &
                                  num_sites, &
                                  num_mpoles_for_mpole_grid, &
                                  off1, &
                                  Global_multipole_coeff, &
                                  site_that_owns_gridded_mpole, &
                                  Gcoeff_offset_for_gridded_mpole, &
                                  Mpole_order_for_gridded_mpole, &
                                  ! note use Ffield arrays here!!(same size)
                                  Ffield_offset_for_gridded_mpole, &
                                  !note use Frac_hermite_field here (same size)
                                  Mpole_grid_Ffield, &
                                  factor, &
                                  energy, &
                                  site_frc, &
                                  Global_multipole_field)
      endif !( struc_fac_method_for_grid_type(gt) == SF_FFP )then
   endif !( tot_num_mpole_coeffs > 0 )then

   ! next the compact sumCH
   if ( tot_num_sumCH_coeffs > 0 )then
      gt = grid_type_for_sumCH
      if ( struc_fac_method_for_grid_type(gt) == SF_FFP )then
         off = 0 ! just one grid
         off1 = 0 ! sumCsumCHset is also 0, since one grid
         ! multipole exponents are diffuse adjusted to grid them
         do n = 1,num_sumCH_for_sumCH_grid
            sumCH_exponent(n) = compact_rec_sumCH_expon
         enddo
         call FH_FFP_RECIP_map_coords( &
                                 num_sites, &
                                 nfft1_for_grid_type(gt), &
                                 nfft2_for_grid_type(gt), &
                                 nfft3_for_grid_type(gt), &
                                 recip, &
                                 ucell, &
                                 site_crd, &
                                 ind_map, &
                                 del_crd)
         call FH_FFP_RECIP_get_GHerm_phi( &
                                 num_sites, &
                                 num_sumCH_for_sumCH_grid, &
                                 off1, &
                                 site_that_owns_gridded_sumCH, &
                                 ind_map, &
                                 del_crd, &
                                 gauss_cut, &
                                 ucell, &
                                 ! note use Ffield arrays here!!(same size)
                                 Ffield_offset_for_gridded_sumCH, &
                                 Ffield_order_for_gridded_sumCH, &
                                 sumCH_exponent, &
                                 nfft1_for_grid_type(gt), &
                                 nfft2_for_grid_type(gt), &
                                 nfft3_for_grid_type(gt), &
                                 nfftdim1_for_gridtype(gt), &
                                 nfftdim2_for_gridtype(gt), &
                                 nfftdim3_for_gridtype(gt), &
                                 !note use sumCH_grid_Ffield here (same size)
                                 sumCH_C_grid(off+1), &
                                 sumCH_grid_Ffield,out_lun )
         call FH_FFP_RECIP_one_ene_frc_field( &
                                  num_sites, &
                                  num_sumCH_for_sumCH_grid, &
                                  off1, &
                                  Global_sumCH_coeff, &
                                  site_that_owns_gridded_sumCH, &
                                  Gcoeff_offset_for_gridded_sumCH, &
                                  sumCH_order_for_gridded_sumCH, &
                                  ! note use Ffield arrays here!!(same size)
                                  Ffield_offset_for_gridded_sumCH, &
                                  !note use sumCH_grid_Ffield here (same size)
                                  sumCH_grid_Ffield, &
                                  factor, &
                                  energy, &
                                  site_frc, &
                                  Global_sumCH_field)
      endif !( struc_fac_method_for_grid_type(gt) == SF_FFP )then
   endif !( tot_num_sumCH_coeffs > 0 )then

   ! next the compact hermite
   if ( tot_num_herm_Cprims > 0 )then
      do g = 1,num_HC_prim_grids
         off = HC_grid_offset_for_HC_grid(g)
         gt = grid_type_for_HC_prim_grid(g)
         if ( struc_fac_method_for_grid_type(gt) == SF_FFP )then
            call FH_FFP_RECIP_map_coords( &
                                 num_sites, &
                                 nfft1_for_grid_type(gt), &
                                 nfft2_for_grid_type(gt), &
                                 nfft3_for_grid_type(gt), &
                                 recip, &
                                 ucell, &
                                 site_crd, &
                                 ind_map, &
                                 del_crd)
            call FH_FFP_RECIP_get_GHerm_phi( &
                                 num_sites, &
                                 num_prims_for_HC_grid(g), &
                                 prim_offset_for_HC_grid(g), &
                                 site_that_owns_HC_grid_prim, &
                                 ind_map, &
                                 del_crd, &
                                 gauss_cut, &
                                 ucell, &
                                 ! note use Ffield arrays here!!(same size)
                                 Ffield_offset_for_HC_grid_prim, &
                                 Ffield_order_for_HC_grid_prim, &
                                 herm_expon_for_HC_grid_prim, &
                                 nfft1_for_grid_type(gt), &
                                 nfft2_for_grid_type(gt), &
                                 nfft3_for_grid_type(gt), &
                                 nfftdim1_for_gridtype(gt), &
                                 nfftdim2_for_gridtype(gt), &
                                 nfftdim3_for_gridtype(gt), &
                                 !note use HC_grid_Ffield here (same size)
                                 HC_grid(off+1), &
                                 HC_grid_Ffield,out_lun )
            call FH_FFP_RECIP_one_ene_frc_field( &
                                  num_sites, &
                                  num_prims_for_HC_grid(g), &
                                  prim_offset_for_HC_grid(g), &
                                  Global_CHermite_coeff, &
                                  site_that_owns_HC_grid_prim, &
                                  Gcoeff_offset_for_HC_grid_prim, &
                                  herm_order_for_HC_grid_prim, &
                                  ! note use Ffield arrays here!!(same size)
                                  Ffield_offset_for_HC_grid_prim, &
                                  !note use Frac_hermite_field here (same size)
                                  HC_grid_Ffield, &
                                  factor, &
                                  energy, &
                                  site_frc, &
                                  Global_CHermite_field)
         endif !( struc_fac_method_for_grid_type(gt) == SF_FFP )then
      enddo !g = 1,num_HC_prim_grids
   endif !( tot_num_herm_Cprims > 0 )then

   ! finally the diffuse hermite
   if ( tot_num_herm_Dprims > 0 )then
      do g = 1,num_HD_prim_grids
         off = HD_grid_offset_for_HD_grid(g)
         gt = grid_type_for_HD_prim_grid(g)
         if ( struc_fac_method_for_grid_type(gt) == SF_FFP )then
            call FH_FFP_RECIP_map_coords( &
                                 num_sites, &
                                 nfft1_for_grid_type(gt), &
                                 nfft2_for_grid_type(gt), &
                                 nfft3_for_grid_type(gt), &
                                 recip, &
                                 ucell, &
                                 site_crd, &
                                 ind_map, &
                                 del_crd)
            call FH_FFP_RECIP_get_GHerm_phi( &
                                 num_sites, &
                                 num_prims_for_HD_grid(g), &
                                 prim_offset_for_HD_grid(g), &
                                 site_that_owns_HD_grid_prim, &
                                 ind_map, &
                                 del_crd, &
                                 gauss_cut, &
                                 ucell, &
                                 ! note use Ffield arrays here!!(same size)
                                 Ffield_offset_for_HD_grid_prim, &
                                 Ffield_order_for_HD_grid_prim, &
                                 herm_expon_for_HD_grid_prim, &
                                 nfft1_for_grid_type(gt), &
                                 nfft2_for_grid_type(gt), &
                                 nfft3_for_grid_type(gt), &
                                 nfftdim1_for_gridtype(gt), &
                                 nfftdim2_for_gridtype(gt), &
                                 nfftdim3_for_gridtype(gt), &
                                 !note use HD_grid_Ffield here (same size)
                                 HD_grid(off+1), &
                                 HD_grid_Ffield,out_lun )
            call FH_FFP_RECIP_one_ene_frc_field( &
                                  num_sites, &
                                  num_prims_for_HD_grid(g), &
                                  prim_offset_for_HD_grid(g), &
                                  Global_DHermite_coeff, &
                                  site_that_owns_HD_grid_prim, &
                                  Gcoeff_offset_for_HD_grid_prim, &
                                  herm_order_for_HD_grid_prim, &
                                  ! note use Ffield arrays here!!(same size)
                                  Ffield_offset_for_HD_grid_prim, &
                                  !note use Frac_hermite_field here (same size)
                                  HD_grid_Ffield, &
                                  factor, &
                                  energy, &
                                  site_frc, &
                                  Global_DHermite_field)
         endif !( struc_fac_method_for_grid_type(gt) == SF_FFP )then
      enddo !g = 1,num_HD_prim_grids
   endif !( tot_num_herm_Dprims > 0 )then

end subroutine FH_FFP_RECIP_ene_force_field
!-----------------------------------------------------------------
end module ffp_recip
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!-----------------------------------------------------------------
! NON-MODULE (LOCAL) SUBROUTINES
! THESE HANDLE ALL ARGUMENTS THROUGH ARGUMENT LIST
!-----------------------------------------------------------------
!-----------------------------------------------------------------

!-----------------------------------------------------------------
subroutine FH_FFP_RECIP_map_coords( &
                             num_sites, &
                             nfft1,nfft2,nfft3, &
                             recip,ucell,site_crd, &
                             ind_map,del_crd)

   implicit none

   integer, intent(in) :: num_sites,nfft1,nfft2,nfft3
   double precision,intent(in) :: recip(3,3),ucell(3,3), &
                                  site_crd(3,num_sites)
   integer,intent(out) :: ind_map(3,num_sites)
   double precision, intent(out) :: del_crd(3,num_sites)

   integer :: n,ix,jx,kx
   double precision :: w,frac(3),map(3),gmap(3)
   
   do n = 1,num_sites
      w = site_crd(1,n)*recip(1,1)+site_crd(2,n)*recip(2,1)+ &
          site_crd(3,n)*recip(3,1)
      frac(1) = w - dnint(w) + 0.5d0
      w = site_crd(1,n)*recip(1,2)+site_crd(2,n)*recip(2,2)+ &
          site_crd(3,n)*recip(3,2)
      frac(2) = w - dnint(w) + 0.5d0
      w = site_crd(1,n)*recip(1,3)+site_crd(2,n)*recip(2,3)+ &
          site_crd(3,n)*recip(3,3)
      frac(3) = w - dnint(w) + 0.5d0
      map(1) = ucell(1,1)*frac(1) + ucell(1,2)*frac(2) + ucell(1,3)*frac(3)
      map(2) = ucell(2,1)*frac(1) + ucell(2,2)*frac(2) + ucell(2,3)*frac(3)
      map(3) = ucell(3,1)*frac(1) + ucell(3,2)*frac(2) + ucell(3,3)*frac(3)
      ix = int(nfft1*frac(1))
      jx = int(nfft2*frac(2))
      kx = int(nfft3*frac(3))
      frac(1) = dble(ix)/nfft1
      frac(2) = dble(jx)/nfft2
      frac(3) = dble(kx)/nfft3
      gmap(1) = ucell(1,1)*frac(1) + ucell(1,2)*frac(2) + ucell(1,3)*frac(3)
      gmap(2) = ucell(2,1)*frac(1) + ucell(2,2)*frac(2) + ucell(2,3)*frac(3)
      gmap(3) = ucell(3,1)*frac(1) + ucell(3,2)*frac(2) + ucell(3,3)*frac(3)
      del_crd(1,n) = map(1) - gmap(1)
      del_crd(2,n) = map(2) - gmap(2)
      del_crd(3,n) = map(3) - gmap(3)
      ind_map(1,n) = ix
      ind_map(2,n) = jx
      ind_map(3,n) = kx
   enddo

end subroutine FH_FFP_RECIP_map_coords
!--------------------------------------------------------------
subroutine FH_FFP_RECIP_grid_GHerm_coeffs( &
                                 nsites, &
                                 num_prims_for_grid, &
                                 prim_offset_for_grid, &
                                 Glob_hermite_coeff, &
                                 site_that_owns_grid_prim, &
                                 ind_map, &
                                 del_crd,gauss_cut,unit_cell, &
                                 Gcoeff_offset_for_grid_prim, &
                                 hermite_order_for_grid_prim, &
                                 hermite_expon_for_grid_prim, &
                                 nfft1,nfft2,nfft3, &
                                 nfftdim1,nfftdim2,nfftdim3,ffp_grid,out_lun)

   implicit none

   integer, intent(in) :: nsites,num_prims_for_grid,prim_offset_for_grid
   double precision,intent(in) :: Glob_hermite_coeff(*)
   integer,intent(in) :: site_that_owns_grid_prim(*)
   integer,intent(in) :: ind_map(3,nsites)
   double precision,intent(in) :: del_crd(3,nsites),gauss_cut, &
                                  unit_cell(3,3)
   integer,intent(in) :: Gcoeff_offset_for_grid_prim(*), &
                         hermite_order_for_grid_prim(*)
   double precision,intent(in) :: hermite_expon_for_grid_prim(*)
   integer,intent(in) :: nfft1,nfft2,nfft3, &
                         nfftdim1,nfftdim2,nfftdim3
   double precision,intent(out) :: ffp_grid(2*nfftdim1,nfftdim2,nfftdim3)
   integer, intent(in)  :: out_lun

   include "mpole_index.fh"

   integer :: ntot,kp,np,n,hc_order,hc_off,imax,jmax,kmax,i,j,k,i0,j0,k0
              
   double precision :: hx,hy,hz,expon,sq_expon,rcut,dis_add, &
                       rcut2,rcut23,weight,pi
   double precision :: term0,term1,term2,term3,term4
   double precision :: t0,t1,t2,t3,t4
   double precision :: u0,u1,u2,u3,u4
   double precision :: v0,v1,v2,v3,v4
   double precision :: dx(-nfft1:nfft1),dx2(-nfft1:nfft1),polx1(-nfft1:nfft1),&
                       polx2(-nfft1:nfft1),polx3(-nfft1:nfft1), &
                       polx4(-nfft1:nfft1),gx(-nfft1:nfft1)
   double precision :: dy(-nfft2:nfft2),dy2(-nfft2:nfft2),poly1(-nfft1:nfft2),&
                       poly2(-nfft2:nfft2),poly3(-nfft2:nfft2), &
                       poly4(-nfft2:nfft2),gy(-nfft2:nfft2)
   double precision :: dz(-nfft3:nfft3),dz2(-nfft3:nfft3),polz1(-nfft1:nfft3),&
                       polz2(-nfft3:nfft3),polz3(-nfft3:nfft3), &
                       polz4(-nfft3:nfft3),gz(-nfft3:nfft3)

   ntot = 2*nfftdim1*nfftdim2*nfftdim3
   call UTIL_zero_real_array(ffp_grid,ntot)
   pi = 3.14159265358979323846d0

   ! assuming orthog unit cell
   hx = unit_cell(1,1) / nfft1
   hy = unit_cell(2,2) / nfft2
   hz = unit_cell(3,3) / nfft3
   do kp = 1,num_prims_for_grid
      np = prim_offset_for_grid + kp
      n = site_that_owns_grid_prim(np)
      hc_order = hermite_order_for_grid_prim(np)
      hc_off = Gcoeff_offset_for_grid_prim(np)
      expon = hermite_expon_for_grid_prim(np)
      sq_expon = sqrt(expon)
      rcut = gauss_cut / sq_expon
      dis_add = sqrt(hx**2+hy**2+hz**2)/2.d0
      rcut = rcut + dis_add
      rcut2 = rcut*rcut
      imax = rcut / hx
      jmax = rcut / hy
      kmax = rcut / hz
      if ( imax > nfft1 .or. jmax > nfft2 .or. kmax > nfft3 )then
         write(out_lun,*)'gaussian too wide for tmp arrays'
         write(out_lun,*)'rcut,imax,jmax,kmax,nfft1,nfft2,nfft3 = ', &
              rcut,imax,jmax,kmax,nfft1,nfft2,nfft3
         stop
      endif
      weight = (expon/pi)*sqrt(expon/pi)*hx*hy*hz !norm of gaussian times dvol
      ! get the one dim factors
      do i = -imax,imax
         dx(i) = sq_expon*(i*hx - del_crd(1,n))
         dx2(i) = dx(i)*dx(i)
      enddo
      do i = -imax,imax
         gx(i) = exp(-dx2(i))
      enddo
      if ( hc_order > 1 )then
         do i = -imax,imax
            polx1(i) = 2.d0*sq_expon*dx(i)
         enddo
         if ( hc_order > 4 )then
            do i = -imax,imax
               polx2(i) = expon*(4.d0*dx2(i) - 2.d0)
            enddo
            if ( hc_order > 10 )then
               do i = -imax,imax
                  polx3(i) = sq_expon*expon*(8.d0*dx(i)*dx2(i) - 12.d0*dx(i))
               enddo
               if ( hc_order > 20 )then
                  do i = -imax,imax
                     polx4(i) = expon*expon* &
                           (16.d0*dx2(i)*dx2(i) - 48.d0*dx2(i) + 12.d0)
                  enddo
               endif
            endif
         endif
      endif
      do j = -jmax,jmax
         dy(j) = sq_expon*(j*hy - del_crd(2,n))
         dy2(j) = dy(j)*dy(j)
      enddo
      do j = -jmax,jmax
         gy(j) = exp(-dy2(j))
      enddo
      if ( hc_order > 1 )then
         do j = -jmax,jmax
            poly1(j) = 2.d0*sq_expon*dy(j)
         enddo
         if ( hc_order > 4 )then
            do j = -jmax,jmax
               poly2(j) = expon*(4.d0*dy2(j) - 2.d0)
            enddo
            if ( hc_order > 10 )then
               do j = -jmax,jmax
                  poly3(j) = sq_expon*expon*(8.d0*dy(j)*dy2(j) - 12.d0*dy(j))
               enddo
               if ( hc_order > 20 )then
                  do j = -jmax,jmax
                     poly4(j) = expon*expon* &
                           (16.d0*dy2(j)*dy2(j) - 48.d0*dy2(j) + 12.d0)
                  enddo
               endif
            endif
         endif
      endif

      do k = -kmax,kmax
         dz(k) = sq_expon*(k*hz - del_crd(3,n))
         dz2(k) = dz(k)*dz(k)
      enddo
      do k = -kmax,kmax
         gz(k) = exp(-dz2(k))
      enddo
      if ( hc_order > 1 )then
         do k = -kmax,kmax
            polz1(k) = 2.d0*sq_expon*dz(k)
         enddo
         if ( hc_order > 4 )then
            do k = -kmax,kmax
               polz2(k) = expon*(4.d0*dz2(k) - 2.d0)
            enddo
            if ( hc_order > 10 )then
               do k = -kmax,kmax
                  polz3(k) = sq_expon*expon*(8.d0*dz(k)*dz2(k) - 12.d0*dz(k))
               enddo
               if ( hc_order > 20 )then
                  do k = -kmax,kmax
                     polz4(k) = expon*expon* &
                           (16.d0*dz2(k)*dz2(k) - 48.d0*dz2(k) + 12.d0)
                  enddo
               endif
            endif
         endif
      endif
      if ( hc_order == 0 )then
         ! do nothing. this prim has no hermite coefficients
      elseif ( hc_order == 1 )then
         do k0 = -kmax,kmax
            v0 = weight*gz(k0)
            k = k0 + ind_map(3,n)
            if ( k < 0 )k = k + nfft3
            if ( k >= nfft3 )k = k - nfft3
            k = k + 1 !fortran indexing in fft
            do j0 = -jmax,jmax
               rcut23 = rcut2 - dy2(j0) - dz2(k0)
               if ( rcut23 > 0 )then
                  u0 = gy(j0)
                  term0 = Glob_hermite_coeff(hc_off+Ind_000)*u0*v0
                  j = j0 + ind_map(2,n)
                  if ( j < 0 )j = j + nfft2
                  if ( j >= nfft2 )j = j - nfft2
                  j = j + 1!fortran indexing in fft
                  do i0 = -imax,imax
                     if ( dx2(i0) < rcut23 ) then
                        t0 = gx(i0)
                        i = i0 + ind_map(1,n)
                        if ( i < 0 )i = i + nfft1
                        if ( i >= nfft1 )i = i - nfft1
                        i = i + 1 !fortran indexing in fft
                        ffp_grid(i,j,k) = ffp_grid(i,j,k) + term0*t0
                     endif
                  enddo
               endif
            enddo
         enddo
      elseif ( hc_order == 4 )then
         do k0 = -kmax,kmax
            v0 = weight*gz(k0)
            v1 = polz1(k0)*v0
            k = k0 + ind_map(3,n)
            if ( k < 0 )k = k + nfft3
            if ( k >= nfft3 )k = k - nfft3
            k = k + 1 !fortran indexing in fft
            do j0 = -jmax,jmax
               rcut23 = rcut2 - dy2(j0) - dz2(k0)
               if ( rcut23 > 0 )then
                  u0 = gy(j0)
                  u1 = poly1(j0)*u0
                  term0 = Glob_hermite_coeff(hc_off+Ind_000)*u0*v0 + &
                          Glob_hermite_coeff(hc_off+Ind_010)*u1*v0 + &
                          Glob_hermite_coeff(hc_off+Ind_001)*u0*v1
                  term1 = Glob_hermite_coeff(hc_off+Ind_100)*u0*v0
                  j = j0 + ind_map(2,n)
                  if ( j < 0 )j = j + nfft2
                  if ( j >= nfft2 )j = j - nfft2
                  j = j + 1!fortran indexing in fft
                  do i0 = -imax,imax
                     if ( dx2(i0) < rcut23 ) then
                        t0 = gx(i0)
                        t1 = polx1(i0)*t0
                        i = i0 + ind_map(1,n)
                        if ( i < 0 )i = i + nfft1
                        if ( i >= nfft1 )i = i - nfft1
                        i = i + 1 !fortran indexing in fft
                        ffp_grid(i,j,k) = ffp_grid(i,j,k) + term0*t0 + &
                                          term1*t1
                     endif
                  enddo
               endif
            enddo
         enddo
      elseif ( hc_order == 10 )then
         do k0 = -kmax,kmax
            v0 = weight*gz(k0)
            v1 = polz1(k0)*v0
            v2 = polz2(k0)*v0
            k = k0 + ind_map(3,n)
            if ( k < 0 )k = k + nfft3
            if ( k >= nfft3 )k = k - nfft3
            k = k + 1 !fortran indexing in fft
            do j0 = -jmax,jmax
               rcut23 = rcut2 - dy2(j0) - dz2(k0)
               if ( rcut23 > 0 )then
                  u0 = gy(j0)
                  u1 = poly1(j0)*u0
                  u2 = poly2(j0)*u0
                  term0 = Glob_hermite_coeff(hc_off+Ind_000)*u0*v0 + &
                          Glob_hermite_coeff(hc_off+Ind_010)*u1*v0 + &
                          Glob_hermite_coeff(hc_off+Ind_001)*u0*v1 + &
                          Glob_hermite_coeff(hc_off+Ind_020)*u2*v0 + &
                          Glob_hermite_coeff(hc_off+Ind_002)*u0*v2 + &
                          Glob_hermite_coeff(hc_off+Ind_011)*u1*v1
                  term1 = Glob_hermite_coeff(hc_off+Ind_100)*u0*v0 + &
                          Glob_hermite_coeff(hc_off+Ind_110)*u1*v0 + &
                          Glob_hermite_coeff(hc_off+Ind_101)*u0*v1
                  term2 = Glob_hermite_coeff(hc_off+Ind_200)*u0*v0
                  j = j0 + ind_map(2,n)
                  if ( j < 0 )j = j + nfft2
                  if ( j >= nfft2 )j = j - nfft2
                  j = j + 1!fortran indexing in fft
                  do i0 = -imax,imax
                     if ( dx2(i0) < rcut23 ) then
                        t0 = gx(i0)
                        t1 = polx1(i0)*t0
                        t2 = polx2(i0)*t0
                        i = i0 + ind_map(1,n)
                        if ( i < 0 )i = i + nfft1
                        if ( i >= nfft1 )i = i - nfft1
                        i = i + 1 !fortran indexing in fft
                        ffp_grid(i,j,k) = ffp_grid(i,j,k) + term0*t0 + &
                                          term1*t1 + term2*t2
                     endif
                  enddo
               endif
            enddo
         enddo
      elseif ( hc_order == 20 )then
         do k0 = -kmax,kmax
            v0 = weight*gz(k0)
            v1 = polz1(k0)*v0
            v2 = polz2(k0)*v0
            v3 = polz3(k0)*v0
            k = k0 + ind_map(3,n)
            if ( k < 0 )k = k + nfft3
            if ( k >= nfft3 )k = k - nfft3
            k = k + 1 !fortran indexing in fft
            do j0 = -jmax,jmax
               rcut23 = rcut2 - dy2(j0) - dz2(k0)
               if ( rcut23 > 0 )then
                  u0 = gy(j0)
                  u1 = poly1(j0)*u0
                  u2 = poly2(j0)*u0
                  u3 = poly3(j0)*u0
                  term0 = Glob_hermite_coeff(hc_off+Ind_000)*u0*v0 + &
                          Glob_hermite_coeff(hc_off+Ind_010)*u1*v0 + &
                          Glob_hermite_coeff(hc_off+Ind_001)*u0*v1 + &
                          Glob_hermite_coeff(hc_off+Ind_020)*u2*v0 + &
                          Glob_hermite_coeff(hc_off+Ind_002)*u0*v2 + &
                          Glob_hermite_coeff(hc_off+Ind_011)*u1*v1 + &
                          Glob_hermite_coeff(hc_off+Ind_030)*u3*v0 + &
                          Glob_hermite_coeff(hc_off+Ind_003)*u0*v3 + &
                          Glob_hermite_coeff(hc_off+Ind_021)*u2*v1 + &
                          Glob_hermite_coeff(hc_off+Ind_012)*u1*v2
                  term1 = Glob_hermite_coeff(hc_off+Ind_100)*u0*v0 + &
                          Glob_hermite_coeff(hc_off+Ind_110)*u1*v0 + &
                          Glob_hermite_coeff(hc_off+Ind_101)*u0*v1 + &
                          Glob_hermite_coeff(hc_off+Ind_120)*u2*v0 + &
                          Glob_hermite_coeff(hc_off+Ind_102)*u0*v2 + &
                          Glob_hermite_coeff(hc_off+Ind_111)*u1*v1
                  term2 = Glob_hermite_coeff(hc_off+Ind_200)*u0*v0 + &
                          Glob_hermite_coeff(hc_off+Ind_210)*u1*v0 + &
                          Glob_hermite_coeff(hc_off+Ind_201)*u0*v1
                  term3 = Glob_hermite_coeff(hc_off+Ind_300)*u0*v0
                  j = j0 + ind_map(2,n)
                  if ( j < 0 )j = j + nfft2
                  if ( j >= nfft2 )j = j - nfft2
                  j = j + 1!fortran indexing in fft
                  do i0 = -imax,imax
                     if ( dx2(i0) < rcut23 ) then
                        t0 = gx(i0)
                        t1 = polx1(i0)*t0
                        t2 = polx2(i0)*t0
                        t3 = polx3(i0)*t0
                        i = i0 + ind_map(1,n)
                        if ( i < 0 )i = i + nfft1
                        if ( i >= nfft1 )i = i - nfft1
                        i = i + 1 !fortran indexing in fft
                        ffp_grid(i,j,k) = ffp_grid(i,j,k) + term0*t0 + &
                                          term1*t1 + term2*t2 + &
                                          term3*t3
                     endif
                  enddo
               endif
            enddo
         enddo
      elseif ( hc_order == 35 )then
         do k0 = -kmax,kmax
            v0 = weight*gz(k0)
            v1 = polz1(k0)*v0
            v2 = polz2(k0)*v0
            v3 = polz3(k0)*v0
            v4 = polz4(k0)*v0
            k = k0 + ind_map(3,n)
            if ( k < 0 )k = k + nfft3
            if ( k >= nfft3 )k = k - nfft3
            k = k + 1 !fortran indexing in fft
            do j0 = -jmax,jmax
               rcut23 = rcut2 - dy2(j0) - dz2(k0)
               if ( rcut23 > 0 )then
                  u0 = gy(j0)
                  u1 = poly1(j0)*u0
                  u2 = poly2(j0)*u0
                  u3 = poly3(j0)*u0
                  u4 = poly4(j0)*u0
! hardwire our knowledge of layout of theta1,2,3 to pre-assemble factors
                  term0 = Glob_hermite_coeff(hc_off+Ind_000)*u0*v0 + &
                          Glob_hermite_coeff(hc_off+Ind_010)*u1*v0 + &
                          Glob_hermite_coeff(hc_off+Ind_001)*u0*v1 + &
                          Glob_hermite_coeff(hc_off+Ind_020)*u2*v0 + &
                          Glob_hermite_coeff(hc_off+Ind_002)*u0*v2 + &
                          Glob_hermite_coeff(hc_off+Ind_011)*u1*v1 + &
                          Glob_hermite_coeff(hc_off+Ind_030)*u3*v0 + &
                          Glob_hermite_coeff(hc_off+Ind_003)*u0*v3 + &
                          Glob_hermite_coeff(hc_off+Ind_021)*u2*v1 + &
                          Glob_hermite_coeff(hc_off+Ind_012)*u1*v2 + &
                          Glob_hermite_coeff(hc_off+Ind_040)*u4*v0 + &
                          Glob_hermite_coeff(hc_off+Ind_004)*u0*v4 + &
                          Glob_hermite_coeff(hc_off+Ind_031)*u3*v1 + &
                          Glob_hermite_coeff(hc_off+Ind_013)*u1*v3 + &
                          Glob_hermite_coeff(hc_off+Ind_022)*u2*v2
                  term1 = Glob_hermite_coeff(hc_off+Ind_100)*u0*v0 + &
                          Glob_hermite_coeff(hc_off+Ind_110)*u1*v0 + &
                          Glob_hermite_coeff(hc_off+Ind_101)*u0*v1 + &
                          Glob_hermite_coeff(hc_off+Ind_120)*u2*v0 + &
                          Glob_hermite_coeff(hc_off+Ind_102)*u0*v2 + &
                          Glob_hermite_coeff(hc_off+Ind_111)*u1*v1 + &
                          Glob_hermite_coeff(hc_off+Ind_130)*u3*v0 + &
                          Glob_hermite_coeff(hc_off+Ind_103)*u0*v3 + &
                          Glob_hermite_coeff(hc_off+Ind_121)*u2*v1 + &
                          Glob_hermite_coeff(hc_off+Ind_112)*u1*v2
                  term2 = Glob_hermite_coeff(hc_off+Ind_200)*u0*v0 + &
                          Glob_hermite_coeff(hc_off+Ind_210)*u1*v0 + &
                          Glob_hermite_coeff(hc_off+Ind_201)*u0*v1 + &
                          Glob_hermite_coeff(hc_off+Ind_220)*u2*v0 + &
                          Glob_hermite_coeff(hc_off+Ind_202)*u0*v2 + &
                          Glob_hermite_coeff(hc_off+Ind_211)*u1*v1
                  term3 = Glob_hermite_coeff(hc_off+Ind_300)*u0*v0 + &
                          Glob_hermite_coeff(hc_off+Ind_310)*u1*v0 + &
                          Glob_hermite_coeff(hc_off+Ind_301)*u0*v1
                  term4 = Glob_hermite_coeff(hc_off+Ind_400)*u0*v0
                  j = j0 + ind_map(2,n)
                  if ( j < 0 )j = j + nfft2
                  if ( j >= nfft2 )j = j - nfft2
                  j = j + 1!fortran indexing in fft
                  do i0 = -imax,imax
                     if ( dx2(i0) < rcut23 ) then
                        t0 = gx(i0)
                        t1 = polx1(i0)*t0
                        t2 = polx2(i0)*t0
                        t3 = polx3(i0)*t0
                        t4 = polx4(i0)*t0
                        i = i0 + ind_map(1,n)
                        if ( i < 0 )i = i + nfft1
                        if ( i >= nfft1 )i = i - nfft1
                        i = i + 1 !fortran indexing in fft
                        ffp_grid(i,j,k) = ffp_grid(i,j,k) + term0*t0 + &
                                          term1*t1 + term2*t2 + &
                                          term3*t3 + term4*t4
                     endif
                  enddo
               endif
            enddo
         enddo
      endif
   enddo

end subroutine FH_FFP_RECIP_grid_GHerm_coeffs
!------------------------------------------------------------------
subroutine FH_FFP_RECIP_one_FT_density( &
                                 nfft1,nfft2,nfft3,nfftdim1, &
                                 volume,ffp_grid)

   implicit none

   integer,intent(in) :: nfft1,nfft2,nfft3,nfftdim1
   double precision,intent(in) :: volume
   double precision, intent(inout) :: ffp_grid(2,nfft3,nfftdim1,nfft2)

   integer :: k1,k2,k3,k10,nf1
   double precision :: volinv

   volinv = 1.d0 / volume

   nf1 = nfft1/2
   if ( 2*nf1 < nfft1 )nf1 = nf1+1
   do k2 = 1, nfft2
      do k3 = 1,nfft3
         k10 = 1
         ! need (1,1,1) case also
         !if ( k3+k2 == 2 )k10 = 2
         do k1 = k10, nf1+1
            ffp_grid(1,k3,k1,k2) = ffp_grid(1,k3,k1,k2) * volinv
            ffp_grid(2,k3,k1,k2) = ffp_grid(2,k3,k1,k2) * volinv
         enddo
      enddo
   enddo
end subroutine FH_FFP_RECIP_one_FT_density
!------------------------------------------------------------------
subroutine FH_FFP_RECIP_one_FT_phi( &
                                 nfft1,nfft2,nfft3,nfftdim1, &
                                 volume,ffp_grid)

   implicit none

   integer,intent(in) :: nfft1,nfft2,nfft3,nfftdim1
   double precision,intent(in) :: volume
   double precision, intent(inout) :: ffp_grid(2,nfft3,nfftdim1,nfft2)

   integer :: k1,k2,k3,k10,nf1
   double precision :: volinv

   volinv = 1.d0 / volume

   nf1 = nfft1/2
   if ( 2*nf1 < nfft1 )nf1 = nf1+1
   do k2 = 1, nfft2
      do k3 = 1,nfft3
         k10 = 1
         ! need (1,1,1) case also
         !if ( k3+k2 == 2 )k10 = 2
         do k1 = k10, nf1+1
            ffp_grid(1,k3,k1,k2) = ffp_grid(1,k3,k1,k2) * volinv
            ffp_grid(2,k3,k1,k2) = ffp_grid(2,k3,k1,k2) * volinv
         enddo
      enddo
   enddo
end subroutine FH_FFP_RECIP_one_FT_phi
!------------------------------------------------------------------
subroutine FH_FFP_RECIP_get_GHerm_phi( &
                                 nsites, &
                                 num_prims_for_grid, &
                                 prim_offset_for_grid, &
                                 site_that_owns_grid_prim, &
                                 ind_map, &
                                 del_crd,gauss_cut,unit_cell, &
                                 Field_offset_for_grid_prim, &
                                 Field_order_for_grid_prim, &
                                 hermite_expon_for_grid_prim, &
                                 nfft1,nfft2,nfft3, &
                                 nfftdim1,nfftdim2,nfftdim3,ffp_grid, &
                                 Gphi,out_lun)

   implicit none

   integer, intent(in) :: nsites,num_prims_for_grid,prim_offset_for_grid
   integer,intent(in) :: site_that_owns_grid_prim(*)
   integer,intent(in) :: ind_map(3,nsites)
   double precision,intent(in) :: del_crd(3,nsites),gauss_cut, &
                                  unit_cell(3,3)
   integer,intent(in) :: Field_offset_for_grid_prim(*), &
                         Field_order_for_grid_prim(*)
   double precision,intent(in) :: hermite_expon_for_grid_prim(*)
   integer,intent(in) :: nfft1,nfft2,nfft3, &
                         nfftdim1,nfftdim2,nfftdim3
   double precision,intent(in) :: ffp_grid(2*nfftdim1,nfftdim2,nfftdim3)
   double precision, intent(out) :: Gphi(*)
   integer, intent(in)          :: out_lun

   include "mpole_index.fh"

   integer :: kp,np,n,field_order,field_off,imax,jmax,kmax, &
              i,j,k,i0,j0,k0
              
   double precision :: hx,hy,hz,expon,sq_expon,rcut,dis_add, &
                       rcut2,rcut23,weight,pi

   double precision :: t0,t1,t2,t3,t4,t5,tt
   double precision :: u0,u1,u2,u3,u4,u5
   double precision :: v0,v1,v2,v3,v4,v5
   double precision :: tu00,tu01,tu10,tu20,tu02,tu11,tu30, &
          tu03,tu21,tu12,tu40,tu04,tu31,tu13,tu22,tu50,tu05, &
          tu41,tu14,tu32,tu23
   double precision :: tuv000,tuv100,tuv010,tuv001
   double precision :: tuv200,tuv020,tuv002,tuv110,tuv101,tuv011
   double precision :: tuv300,tuv030,tuv003,tuv210,tuv201,tuv120, &
          tuv021,tuv102,tuv012,tuv111
   double precision :: tuv400,tuv040,tuv004,tuv310,tuv301, &
          tuv130,tuv031,tuv103,tuv013,tuv220,tuv202,tuv022, &
          tuv211, tuv121,tuv112
   double precision :: tuv500,tuv050,tuv005,tuv410,tuv401, &
          tuv140,tuv041,tuv104,tuv014,tuv320,tuv302,tuv230, &
          tuv032,tuv203,tuv023,tuv311,tuv131,tuv113,tuv221, &
          tuv212,tuv122
   double precision :: tq
   double precision :: dx(-nfft1:nfft1),dx2(-nfft1:nfft1),polx1(-nfft1:nfft1),&
                       polx2(-nfft1:nfft1),polx3(-nfft1:nfft1), &
                       polx4(-nfft1:nfft1),polx5(-nfft1:nfft1), &
                       gx(-nfft1:nfft1)
   double precision :: dy(-nfft2:nfft2),dy2(-nfft2:nfft2),poly1(-nfft1:nfft2),&
                       poly2(-nfft2:nfft2),poly3(-nfft2:nfft2), &
                       poly4(-nfft2:nfft2),poly5(-nfft2:nfft2), &
                       gy(-nfft2:nfft2)
   double precision :: dz(-nfft3:nfft3),dz2(-nfft3:nfft3),polz1(-nfft1:nfft3),&
                       polz2(-nfft3:nfft3),polz3(-nfft3:nfft3), &
                       polz4(-nfft3:nfft3),polz5(-nfft3:nfft3), &
                       gz(-nfft3:nfft3)

   ! assuming orthog unit cell
   hx = unit_cell(1,1) / nfft1
   hy = unit_cell(2,2) / nfft2
   hz = unit_cell(3,3) / nfft3
   pi = 3.14159265358979323846d0
   do kp = 1,num_prims_for_grid
      np = prim_offset_for_grid + kp
      n = site_that_owns_grid_prim(np)
      field_order = Field_order_for_grid_prim(np)
      field_off = Field_offset_for_grid_prim(np)
      expon = hermite_expon_for_grid_prim(np)
      sq_expon = sqrt(expon)
      rcut = gauss_cut / sq_expon
      dis_add = sqrt(hx**2+hy**2+hz**2)/2.d0
      rcut = rcut + dis_add
      rcut2 = rcut*rcut
      imax = rcut / hx
      jmax = rcut / hy
      kmax = rcut / hz
      if ( imax > nfft1 .or. jmax > nfft2 .or. kmax > nfft3 )then
         write(out_lun,*)'gaussian too wide for tmp arrays'
         write(out_lun,*)'rcut,imax,jmax,kmax,nfft1,nfft2,nfft3 = ', &
              rcut,imax,jmax,kmax,nfft1,nfft2,nfft3
         stop
      endif
      weight = (expon/pi)*sqrt(expon/pi)*hx*hy*hz !norm of gaussian times dvol
      ! get the one dim factors
      do i = -imax,imax
         dx(i) = sq_expon*(i*hx - del_crd(1,n))
         dx2(i) = dx(i)*dx(i)
      enddo
      do i = -imax,imax
         gx(i) = exp(-dx2(i))
      enddo
      if ( field_order > 1 )then
         do i = -imax,imax
            polx1(i) = 2.d0*sq_expon*dx(i)
         enddo
         if ( field_order > 4 )then
            do i = -imax,imax
               polx2(i) = expon*(4.d0*dx2(i) - 2.d0)
            enddo
            if ( field_order > 10 )then
               do i = -imax,imax
                  polx3(i) = sq_expon*expon*(8.d0*dx(i)*dx2(i) - 12.d0*dx(i))
               enddo
               if ( field_order > 20 )then
                  do i = -imax,imax
                     polx4(i) = expon*expon* &
                           (16.d0*dx2(i)*dx2(i) - 48.d0*dx2(i) + 12.d0)
                  enddo
                  if ( field_order > 35 )then
                     do i = -imax,imax
                        polx5(i) = expon*expon*sq_expon* &
                           (32.d0*dx2(i)*dx2(i)*dx(i) - 160.d0*dx2(i)*dx(i) + &
                            120.d0*dx(i))
                     enddo
                  endif
               endif
            endif
         endif
      endif
      do j = -jmax,jmax
         dy(j) = sq_expon*(j*hy - del_crd(2,n))
         dy2(j) = dy(j)*dy(j)
      enddo
      do j = -jmax,jmax
         gy(j) = exp(-dy2(j))
      enddo
      if ( field_order > 1 )then
         do j = -jmax,jmax
            poly1(j) = 2.d0*sq_expon*dy(j)
         enddo
         if ( field_order > 4 )then
            do j = -jmax,jmax
               poly2(j) = expon*(4.d0*dy2(j) - 2.d0)
            enddo
            if ( field_order > 10 )then
               do j = -jmax,jmax
                  poly3(j) = sq_expon*expon*(8.d0*dy(j)*dy2(j) - 12.d0*dy(j))
               enddo
               if ( field_order > 20 )then
                  do j = -jmax,jmax
                     poly4(j) = expon*expon* &
                           (16.d0*dy2(j)*dy2(j) - 48.d0*dy2(j) + 12.d0)
                  enddo
                  if ( field_order > 35 )then
                     do j = -jmax,jmax
                        poly5(j) = expon*expon*sq_expon* &
                           (32.d0*dy2(j)*dy2(j)*dy(j) - 160.d0*dy2(j)*dy(j) + &
                            120.d0*dy(j))
                     enddo
                  endif
               endif
            endif
         endif
      endif

      do k = -kmax,kmax
         dz(k) = sq_expon*(k*hz - del_crd(3,n))
         dz2(k) = dz(k)*dz(k)
      enddo
      do k = -kmax,kmax
         gz(k) = exp(-dz2(k))
      enddo
      if ( field_order > 1 )then
         do k = -kmax,kmax
            polz1(k) = 2.d0*sq_expon*dz(k)
         enddo
         if ( field_order > 4 )then
            do k = -kmax,kmax
               polz2(k) = expon*(4.d0*dz2(k) - 2.d0)
            enddo
            if ( field_order > 10 )then
               do k = -kmax,kmax
                  polz3(k) = sq_expon*expon*(8.d0*dz(k)*dz2(k) - 12.d0*dz(k))
               enddo
               if ( field_order > 20 )then
                  do k = -kmax,kmax
                     polz4(k) = expon*expon* &
                           (16.d0*dz2(k)*dz2(k) - 48.d0*dz2(k) + 12.d0)
                  enddo
                  if ( field_order > 35 )then
                     do k = -kmax,kmax
                        polz5(k) = expon*expon*sq_expon* &
                           (32.d0*dz2(k)*dz2(k)*dz(k) - 160.d0*dz2(k)*dz(k) + &
                            120.d0*dz(k))
                     enddo
                  endif
               endif
            endif
         endif
      endif
      if ( field_order == 0 )then
         ! do nothing. this site has no coefficients for this hermite set
      elseif ( field_order > 56 )then ! later maybe fill this in with general
         write(out_lun,*)'field order > 56!!!', field_order
         stop 
      else !field_ord > 0
         tuv000 = 0.d0
         if ( field_order > 1 )then
            tuv001 = 0.d0
            tuv010 = 0.d0
            tuv100 = 0.d0
            if ( field_order > 4 )then
               tuv200 = 0.d0
               tuv020 = 0.d0
               tuv002 = 0.d0
               tuv110 = 0.d0
               tuv101 = 0.d0
               tuv011 = 0.d0
               if ( field_order > 10 )then
                  tuv300 = 0.d0
                  tuv030 = 0.d0
                  tuv003 = 0.d0
                  tuv210 = 0.d0
                  tuv201 = 0.d0
                  tuv120 = 0.d0
                  tuv021 = 0.d0
                  tuv102 = 0.d0
                  tuv012 = 0.d0
                  tuv111 = 0.d0
                  if ( field_order > 20 )then
                     tuv400 = 0.d0
                     tuv040 = 0.d0
                     tuv004 = 0.d0
                     tuv310 = 0.d0
                     tuv301 = 0.d0
                     tuv130 = 0.d0
                     tuv031 = 0.d0
                     tuv103 = 0.d0
                     tuv013 = 0.d0
                     tuv220 = 0.d0
                     tuv202 = 0.d0
                     tuv022 = 0.d0
                     tuv211 = 0.d0
                     tuv121 = 0.d0
                     tuv112 = 0.d0
                     if ( field_order == 56 )then
                        tuv500 = 0.d0
                        tuv050 = 0.d0
                        tuv005 = 0.d0
                        tuv410 = 0.d0
                        tuv401 = 0.d0
                        tuv140 = 0.d0
                        tuv041 = 0.d0
                        tuv104 = 0.d0
                        tuv014 = 0.d0
                        tuv320 = 0.d0
                        tuv302 = 0.d0
                        tuv230 = 0.d0
                        tuv032 = 0.d0
                        tuv203 = 0.d0
                        tuv023 = 0.d0
                        tuv311 = 0.d0
                        tuv131 = 0.d0
                        tuv113 = 0.d0
                        tuv221 = 0.d0
                        tuv212 = 0.d0
                        tuv122 = 0.d0
                     endif !( field_ord == 56 )
                  endif !( field_ord > 20 )
               endif !( field_ord > 10 )
            endif !( field_ord > 4 )
         endif !(field_ord > 1 )
         do k0 = -kmax,kmax
            v0 = gz(k0) 
            tu00 = 0.d0
            if ( field_order > 1 )then
               v1 = polz1(k0)*v0 !1st deriv of gz
               tu10 = 0.d0
               tu01 = 0.d0
               if ( field_order > 4 )then
                  v2 = polz2(k0)*v0 !2nd deriv of gz
                  tu20 = 0.d0
                  tu11 = 0.d0
                  tu02 = 0.d0
                  if( field_order > 10 )then
                     v3 = polz3(k0)*v0 !3rd deriv of gz
                     tu30 = 0.d0
                     tu21 = 0.d0
                     tu12 = 0.d0
                     tu03 = 0.d0
                     if( field_order > 20 )then
                        v4 = polz4(k0)*v0 !4th deriv of gz
                        tu40 = 0.d0
                        tu31 = 0.d0
                        tu22 = 0.d0
                        tu13 = 0.d0
                        tu04 = 0.d0
                        if ( field_order == 56 )then
                           v5 = polz5(k0)*v0 !5th deriv of gz
                           tu50 = 0.d0
                           tu41 = 0.d0
                           tu32 = 0.d0
                           tu23 = 0.d0
                           tu14 = 0.d0
                           tu05 = 0.d0
                        endif !( field_order == 56 )
                     endif !( field_order > 20 )
                  endif !( field_order > 10 )
               endif !( field_order > 4 )
            endif !( field_order > 1 )
            k = k0 + ind_map(3,n)
            if ( k < 0 )k = k + nfft3
            if ( k >= nfft3 )k = k - nfft3
            k = k + 1 !fortran indexing in fft
            do j0 = -jmax,jmax
               j = j0 + ind_map(2,n)
               if ( j < 0 )j = j + nfft2
               if ( j >= nfft2 )j = j - nfft2
               j = j + 1 !fortran indexing in fft
               if ( field_order == 1 )then
                  u0 = gy(j0) 
                  t0 = 0.d0
                  do i0 = -imax,imax
                     i = i0 + ind_map(1,n)
                     if ( i < 0 )i = i + nfft1
                     if ( i >= nfft1 )i = i - nfft1
                     i = i + 1 !fortran indexing in fft
                     tq = ffp_grid(i,j,k)*gx(i0)
                     t0 = t0 + tq
                  enddo !i0 = -imax,imax
                  tu00 = tu00 + t0*u0
               elseif ( field_order == 4 )then
                  u0 = gy(j0) 
                  u1 = poly1(j0)*u0
                  t0 = 0.d0
                  t1 = 0.d0
                  do i0 = -imax,imax
                     i = i0 + ind_map(1,n)
                     if ( i < 0 )i = i + nfft1
                     if ( i >= nfft1 )i = i - nfft1
                     i = i + 1 !fortran indexing in fft
                     tq = ffp_grid(i,j,k)*gx(i0)
                     t0 = t0 + tq
                     t1 = t1 + tq*polx1(i0)
                  enddo !i0 = -imax,imax
                  tu00 = tu00 + t0*u0

                  tu10 = tu10 + t1*u0
                  tu01 = tu01 + t0*u1
               elseif ( field_order == 10 )then
                  u0 = gy(j0) 
                  u1 = poly1(j0)*u0
                  u2 = poly2(j0)*u0
                  t0 = 0.d0
                  t1 = 0.d0
                  t2 = 0.d0
                  do i0 = -imax,imax
                     i = i0 + ind_map(1,n)
                     if ( i < 0 )i = i + nfft1
                     if ( i >= nfft1 )i = i - nfft1
                     i = i + 1 !fortran indexing in fft
                     tq = ffp_grid(i,j,k)*gx(i0)
                     t0 = t0 + tq
                     t1 = t1 + tq*polx1(i0)
                     t2 = t2 + tq*polx2(i0)
                  enddo !i0 = -imax,imax
                  tu00 = tu00 + t0*u0

                  tu10 = tu10 + t1*u0
                  tu01 = tu01 + t0*u1

                  tu20 = tu20 + t2*u0
                  tu11 = tu11 + t1*u1
                  tu02 = tu02 + t0*u2
               elseif ( field_order == 20 )then
                  u0 = gy(j0) 
                  u1 = poly1(j0)*u0
                  u2 = poly2(j0)*u0
                  u3 = poly3(j0)*u0
                  t0 = 0.d0
                  t1 = 0.d0
                  t2 = 0.d0
                  t3 = 0.d0
                  do i0 = -imax,imax
                     i = i0 + ind_map(1,n)
                     if ( i < 0 )i = i + nfft1
                     if ( i >= nfft1 )i = i - nfft1
                     i = i + 1 !fortran indexing in fft
                     tq = ffp_grid(i,j,k)*gx(i0)
                     t0 = t0 + tq
                     t1 = t1 + tq*polx1(i0)
                     t2 = t2 + tq*polx2(i0)
                     t3 = t3 + tq*polx3(i0)
                  enddo !i0 = -imax,imax
                  tu00 = tu00 + t0*u0

                  tu10 = tu10 + t1*u0
                  tu01 = tu01 + t0*u1

                  tu20 = tu20 + t2*u0
                  tu11 = tu11 + t1*u1
                  tu02 = tu02 + t0*u2

                  tu30 = tu30 + t3*u0 
                  tu21 = tu21 + t2*u1 
                  tu12 = tu12 + t1*u2 
                  tu03 = tu03 + t0*u3
               elseif ( field_order == 35 )then
                  u0 = gy(j0) 
                  u1 = poly1(j0)*u0
                  u2 = poly2(j0)*u0
                  u3 = poly3(j0)*u0
                  u4 = poly4(j0)*u0
                  t0 = 0.d0
                  t1 = 0.d0
                  t2 = 0.d0
                  t3 = 0.d0
                  t4 = 0.d0
                  do i0 = -imax,imax
                     i = i0 + ind_map(1,n)
                     if ( i < 0 )i = i + nfft1
                     if ( i >= nfft1 )i = i - nfft1
                     i = i + 1 !fortran indexing in fft
                     tq = ffp_grid(i,j,k)*gx(i0)
                     t0 = t0 + tq
                     t1 = t1 + tq*polx1(i0)
                     t2 = t2 + tq*polx2(i0)
                     t3 = t3 + tq*polx3(i0)
                     t4 = t4 + tq*polx4(i0)
                  enddo !i0 = -imax,imax
                  tu00 = tu00 + t0*u0

                  tu10 = tu10 + t1*u0
                  tu01 = tu01 + t0*u1

                  tu20 = tu20 + t2*u0
                  tu11 = tu11 + t1*u1
                  tu02 = tu02 + t0*u2

                  tu30 = tu30 + t3*u0 
                  tu21 = tu21 + t2*u1 
                  tu12 = tu12 + t1*u2 
                  tu03 = tu03 + t0*u3

                  tu40 = tu40 + t4*u0
                  tu31 = tu31 + t3*u1
                  tu22 = tu22 + t2*u2
                  tu13 = tu13 + t1*u3
                  tu04 = tu04 + t0*u4
               elseif ( field_order == 56 )then
                  u0 = gy(j0) 
                  u1 = poly1(j0)*u0
                  u2 = poly2(j0)*u0
                  u3 = poly3(j0)*u0
                  u4 = poly4(j0)*u0
                  u5 = poly5(j0)*u0
                  t0 = 0.d0
                  t1 = 0.d0
                  t2 = 0.d0
                  t3 = 0.d0
                  t4 = 0.d0
                  t5 = 0.d0
                  do i0 = -imax,imax
                     i = i0 + ind_map(1,n)
                     if ( i < 0 )i = i + nfft1
                     if ( i >= nfft1 )i = i - nfft1
                     i = i + 1 !fortran indexing in fft
                     tq = ffp_grid(i,j,k)*gx(i0)
                     t0 = t0 + tq
                     t1 = t1 + tq*polx1(i0)
                     t2 = t2 + tq*polx2(i0)
                     t3 = t3 + tq*polx3(i0)
                     t4 = t4 + tq*polx4(i0)
                     t5 = t5 + tq*polx5(i0)
                  enddo !i0 = -imax,imax
                  tu00 = tu00 + t0*u0

                  tu10 = tu10 + t1*u0
                  tu01 = tu01 + t0*u1

                  tu20 = tu20 + t2*u0
                  tu11 = tu11 + t1*u1
                  tu02 = tu02 + t0*u2

                  tu30 = tu30 + t3*u0 
                  tu21 = tu21 + t2*u1 
                  tu12 = tu12 + t1*u2 
                  tu03 = tu03 + t0*u3

                  tu40 = tu40 + t4*u0
                  tu31 = tu31 + t3*u1
                  tu22 = tu22 + t2*u2
                  tu13 = tu13 + t1*u3
                  tu04 = tu04 + t0*u4

                  tu50 = tu50 + t5*u0
                  tu41 = tu41 + t4*u1
                  tu32 = tu32 + t3*u2
                  tu23 = tu23 + t2*u3
                  tu14 = tu14 + t1*u4
                  tu05 = tu05 + t0*u5
               endif !( field_ord == 1 )
            enddo !j0 = -jmax,jmax
            tuv000 = tuv000 + tu00*v0
            if ( field_order > 1 )then 
               tuv100 = tuv100 + tu10*v0
               tuv010 = tuv010 + tu01*v0
               tuv001 = tuv001 + tu00*v1
               if ( field_order > 4 )then
                  tuv200 = tuv200 + tu20*v0
                  tuv020 = tuv020 + tu02*v0
                  tuv002 = tuv002 + tu00*v2
                  tuv110 = tuv110 + tu11*v0
                  tuv101 = tuv101 + tu10*v1
                  tuv011 = tuv011 + tu01*v1
                  if ( field_order > 10 )then
                     tuv300 = tuv300 + tu30*v0
                     tuv030 = tuv030 + tu03*v0
                     tuv003 = tuv003 + tu00*v3
                     tuv210 = tuv210 + tu21*v0
                     tuv201 = tuv201 + tu20*v1
                     tuv120 = tuv120 + tu12*v0
                     tuv021 = tuv021 + tu02*v1
                     tuv102 = tuv102 + tu10*v2
                     tuv012 = tuv012 + tu01*v2
                     tuv111 = tuv111 + tu11*v1
                     if ( field_order > 20 )then
                        tuv400 = tuv400 + tu40*v0
                        tuv040 = tuv040 + tu04*v0
                        tuv004 = tuv004 + tu00*v4
                        tuv310 = tuv310 + tu31*v0
                        tuv301 = tuv301 + tu30*v1
                        tuv130 = tuv130 + tu13*v0
                        tuv031 = tuv031 + tu03*v1
                        tuv103 = tuv103 + tu10*v3
                        tuv013 = tuv013 + tu01*v3
                        tuv220 = tuv220 + tu22*v0
                        tuv202 = tuv202 + tu20*v2
                        tuv022 = tuv022 + tu02*v2
                        tuv211 = tuv211 + tu21*v1
                        tuv121 = tuv121 + tu12*v1
                        tuv112 = tuv112 + tu11*v2
                        if ( field_order == 56 )then
                           tuv500 = tuv500 + tu50*v0
                           tuv050 = tuv050 + tu05*v0
                           tuv005 = tuv005 + tu00*v5
                           tuv410 = tuv410 + tu41*v0
                           tuv401 = tuv401 + tu40*v1
                           tuv140 = tuv140 + tu14*v0
                           tuv041 = tuv041 + tu04*v1
                           tuv104 = tuv104 + tu10*v4
                           tuv014 = tuv014 + tu01*v4
                           tuv320 = tuv320 + tu32*v0
                           tuv302 = tuv302 + tu30*v2
                           tuv230 = tuv230 + tu23*v0
                           tuv032 = tuv032 + tu03*v2
                           tuv203 = tuv203 + tu20*v3
                           tuv023 = tuv023 + tu02*v3
                           tuv311 = tuv311 + tu31*v1
                           tuv131 = tuv131 + tu13*v1
                           tuv113 = tuv113 + tu11*v3
                           tuv221 = tuv221 + tu22*v1
                           tuv212 = tuv212 + tu21*v2
                           tuv122 = tuv122 + tu12*v2
                        endif !( field_order == 56 )
                     endif !( field_order > 20 )
                  endif !( field_order > 10 )
               endif !( field_order > 4 )
            endif !( field_order > 1 )
         enddo !k0 = -kmax,kmax
         ! now fill fields; field_off starts this atoms fields
         Gphi(field_off+Ind_000) = tuv000
         if ( field_order > 1 )then
            Gphi(field_off+Ind_100) = tuv100
            Gphi(field_off+Ind_010) = tuv010
            Gphi(field_off+Ind_001) = tuv001
            if ( field_order > 4 )then
               Gphi(field_off+Ind_200) = tuv200
               Gphi(field_off+Ind_020) = tuv020
               Gphi(field_off+Ind_002) = tuv002
               Gphi(field_off+Ind_110) = tuv110
               Gphi(field_off+Ind_101) = tuv101
               Gphi(field_off+Ind_011) = tuv011
               if ( field_order > 10 )then
                  Gphi(field_off+Ind_300) = tuv300
                  Gphi(field_off+Ind_030) = tuv030
                  Gphi(field_off+Ind_003) = tuv003
                  Gphi(field_off+Ind_210) = tuv210
                  Gphi(field_off+Ind_201) = tuv201
                  Gphi(field_off+Ind_120) = tuv120
                  Gphi(field_off+Ind_021) = tuv021
                  Gphi(field_off+Ind_102) = tuv102
                  Gphi(field_off+Ind_012) = tuv012
                  Gphi(field_off+Ind_111) = tuv111
                  if ( field_order > 20 )then
                     Gphi(field_off+Ind_400) = tuv400
                     Gphi(field_off+Ind_040) = tuv040
                     Gphi(field_off+Ind_004) = tuv004
                     Gphi(field_off+Ind_310) = tuv310
                     Gphi(field_off+Ind_301) = tuv301
                     Gphi(field_off+Ind_130) = tuv130
                     Gphi(field_off+Ind_031) = tuv031
                     Gphi(field_off+Ind_103) = tuv103
                     Gphi(field_off+Ind_013) = tuv013
                     Gphi(field_off+Ind_220) = tuv220
                     Gphi(field_off+Ind_202) = tuv202
                     Gphi(field_off+Ind_022) = tuv022
                     Gphi(field_off+Ind_211) = tuv211
                     Gphi(field_off+Ind_121) = tuv121
                     Gphi(field_off+Ind_112) = tuv112
                     if ( field_order == 56 )then
                        Gphi(field_off+Ind_500) = tuv500
                        Gphi(field_off+Ind_050) = tuv050
                        Gphi(field_off+Ind_005) = tuv005
                        Gphi(field_off+Ind_410) = tuv410
                        Gphi(field_off+Ind_401) = tuv401
                        Gphi(field_off+Ind_140) = tuv140
                        Gphi(field_off+Ind_041) = tuv041
                        Gphi(field_off+Ind_104) = tuv104
                        Gphi(field_off+Ind_014) = tuv014
                        Gphi(field_off+Ind_320) = tuv320
                        Gphi(field_off+Ind_302) = tuv302
                        Gphi(field_off+Ind_230) = tuv230
                        Gphi(field_off+Ind_032) = tuv032
                        Gphi(field_off+Ind_203) = tuv203
                        Gphi(field_off+Ind_023) = tuv023
                        Gphi(field_off+Ind_311) = tuv311
                        Gphi(field_off+Ind_131) = tuv131
                        Gphi(field_off+Ind_113) = tuv113
                        Gphi(field_off+Ind_221) = tuv221
                        Gphi(field_off+Ind_212) = tuv212
                        Gphi(field_off+Ind_122) = tuv122
                     endif !( field_order == 56 )
                  endif !( field_order > 20 )
               endif !( field_order > 10 )
            endif !( field_order > 4 )
         endif !( field_order > 1 ) 
         ! scale by weight
         do j = 1,field_order
            Gphi(field_off+j) = weight*Gphi(field_off+j)
         enddo
      endif !( field_order == 0 ) outermost branch on field order
   enddo !kp = 1,num_prims_for_grid
end subroutine FH_FFP_RECIP_get_GHerm_phi
!------------------------------------------------------------------
subroutine FH_FFP_RECIP_one_ene_frc_field( &
                                 nsites, &
                                 num_prims_for_grid, &
                                 prim_offset_for_grid, &
                                 Glob_hermite_coeff, &
                                 site_that_owns_grid_prim, &
                                 Gcoeff_offset_for_grid_prim, &
                                 hermite_order_for_grid_prim, &
                                 Field_offset_for_grid_prim, &
                                 Gphi,factor,&
                                 energy,force,Global_hermite_field)

   implicit none

   integer, intent(in) :: nsites,num_prims_for_grid,prim_offset_for_grid
   double precision,intent(in) :: Glob_hermite_coeff(*)
   integer,intent(in) :: site_that_owns_grid_prim(*)
   integer,intent(in) :: Gcoeff_offset_for_grid_prim(*), &
                         hermite_order_for_grid_prim(*)
   integer,intent(in) :: Field_offset_for_grid_prim(*)
   double precision, intent(in) :: Gphi(*),factor
   double precision, intent(inout) :: energy,force(3,*)
   double precision, intent(inout) :: Global_hermite_field(*)

   include "direct_pointers.fh"
   double precision :: ene,dfx,dfy,dfz
   integer :: j1,j2,j3,j,kp,n,np,order,off1,off2

   ene = 0.d0

   do kp = 1,num_prims_for_grid
      np = prim_offset_for_grid + kp
      n = site_that_owns_grid_prim(np)
      order = hermite_order_for_grid_prim(np)
      off1 = Gcoeff_offset_for_grid_prim(np)
      off2 = Field_offset_for_grid_prim(np)
      dfx = 0.d0
      dfy = 0.d0
      dfz = 0.d0
      do j = 1,order
         ene = ene + Glob_hermite_coeff(off1+j)*Gphi(off2+j)
         ! update the global hermite field
         Global_hermite_field(off1+j) = Global_hermite_field(off1+j) + &
                                        factor*Gphi(off2+j)
         j1 = p_xder(j) 
         j2 = p_yder(j) 
         j3 = p_zder(j)
         ! update cartesian gradients
         dfx = dfx + Glob_hermite_coeff(off1+j)*Gphi(off2+j1)
         dfy = dfy + Glob_hermite_coeff(off1+j)*Gphi(off2+j2)
         dfz = dfz + Glob_hermite_coeff(off1+j)*Gphi(off2+j3)
      enddo
      ! force is negative of gradient
      force(1,n) = force(1,n) - factor*dfx
      force(2,n) = force(2,n) - factor*dfy
      force(3,n) = force(3,n) - factor*dfz
   enddo
   energy = energy + 0.5d0*factor*ene
end subroutine FH_FFP_RECIP_one_ene_frc_field
!------------------------------------------------------------------
