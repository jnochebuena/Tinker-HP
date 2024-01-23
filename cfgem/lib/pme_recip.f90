module pme_recip
   implicit none
   private

   ! multipoles
   double precision,save, allocatable :: MC_pme_multiplier(:)
   double precision,save, allocatable :: MD_pme_multiplier(:)
   double precision, save, allocatable :: MPOLE_prefac1(:), &
                                          MPOLE_prefac2(:), &
                                          MPOLE_prefac3(:)
   integer, save               :: tot_siz_MPOLE_prefac, &
                                  tot_siz_MPOLE_pme_mult
   integer, save              :: num_MC_grid_expon, &
                                 num_MD_grid_expon
   double precision, save, dimension(1) :: MC_weight_grid_expon, &
                                           MC_grid_expon, &
                                           MD_weight_grid_expon, &
                                           MD_grid_expon

   ! sum of compact hermites
   double precision,save, allocatable :: sumCH_C_pme_multiplier(:)
   double precision,save, allocatable :: sumCH_D_pme_multiplier(:)
   double precision, save, allocatable :: sumCH_prefac1(:), &
                                          sumCH_prefac2(:), &
                                          sumCH_prefac3(:)
   integer, save               :: tot_siz_sumCH_prefac, &
                                  tot_siz_sumCH_pme_mult
   integer, save              :: num_sumCH_C_grid_expon, &
                                 num_sumCH_D_grid_expon
   double precision, save, dimension(1) :: sumCH_C_weight_grid_expon, &
                                           sumCH_C_grid_expon, &
                                           sumCH_D_weight_grid_expon, &
                                           sumCH_D_grid_expon

   ! compact hermites
   integer, save, allocatable :: HC_prefac_off_for_HC_grid(:), &
                                 HC_pme_mult_off_for_HC_grid(:)
   integer, save               :: tot_siz_HC_prefac, &
                                  tot_siz_HC_pme_mult
   double precision,save, allocatable :: HC_pme_multiplier(:)
   double precision, save, allocatable :: HC_prefac1(:), &
                                          HC_prefac2(:), &
                                          HC_prefac3(:)
   integer, save              :: siz_HC_grid_expon
   integer, save, allocatable :: num_expon_for_HC_grid(:),  &
                                 expon_off_for_HC_grid(:)
   double precision, save, allocatable :: HC_grid_expon(:), &
                                          HC_weight_grid_expon(:)

   ! diffuse hermites
   integer, save, allocatable :: HD_prefac_off_for_HD_grid(:), &
                                 HD_pme_mult_off_for_HD_grid(:)
   integer, save               :: tot_siz_HD_prefac, &
                                  tot_siz_HD_pme_mult
   double precision,save, allocatable :: HD_pme_multiplier(:)
   double precision, save, allocatable :: HD_prefac1(:), &
                                          HD_prefac2(:), &
                                          HD_prefac3(:)
   integer, save              :: siz_HD_grid_expon
   integer, save, allocatable :: num_expon_for_HD_grid(:), &
                                 expon_off_for_HD_grid(:)
   double precision, save, allocatable :: HD_grid_expon(:), &
                                          HD_weight_grid_expon(:)


   integer, save, allocatable :: Max_deriv_per_site_gridtype(:,:)
   integer, save, allocatable :: Spline_offset_for_site_gridtype(:,:)
   integer, save :: Bspline_table_size
   double precision, save, allocatable :: theta1(:),&
                                          theta2(:), &
                                          theta3(:)
   integer, save, allocatable :: init_grid_ind_for_site_gridtype(:,:,:)
!-----------------------------------------------------------------
   public FH_PME_RECIP_setup, &
          FH_PME_RECIP_deallocate, &
          FH_PME_RECIP_global_to_frac, &
          FH_PME_RECIP_Bspline_fill, &
          FH_PME_RECIP_fill_grids, &
          FH_PME_RECIP_grid_mult, &
          FH_PME_RECIP_FT_density, &
          FH_PME_RECIP_fix_FT_phi, &
          FH_PME_recip_ene_force_field
contains
!-----------------------------------------------------------------
!-----------------------------------------------------------------
! SETUP SUBROUTINES
!-----------------------------------------------------------------
!-----------------------------------------------------------------

subroutine FH_PME_RECIP_setup(out_lun)

   implicit none

   integer, intent(in)  :: out_lun

   call FH_PME_RECIP_mult_expon_grids(out_lun)
   call FH_PME_RECIP_setup_Bsplines(out_lun)
   call FH_PME_RECIP_grid_setup(out_lun)

   return

end subroutine FH_PME_RECIP_setup

!-----------------------------------------------------------------
subroutine FH_PME_RECIP_mult_expon_grids(out_lun)

   use user, only : num_HC_prim_grids, &
                    num_HD_prim_grids, &
                    grid_type_for_HC_prim_grid, &
                    grid_type_for_HD_prim_grid, &
                    struc_fac_method_for_grid_type
   use hermite, only : Local_Chermite_coeff, &
                       Local_Dhermite_coeff, &
                       compact_rec_mpole_expon, &
                       diffuse_rec_mpole_expon, &
                       compact_rec_sumCH_expon, &
                       diffuse_rec_sumCH_expon, &
                       tot_num_herm_Cprims, &
                       tot_num_herm_Dprims
                       
   use recip_grids, only : prim_offset_for_HC_grid, &
                           num_prims_for_HC_grid, &
                           herm_expon_for_HC_grid_prim, &
                           herm_order_for_HC_grid_prim, &
                           prim_offset_for_HD_grid, &
                           num_prims_for_HD_grid, &
                           herm_expon_for_HD_grid_prim, &
                           herm_order_for_HD_grid_prim, &
                           Gcoeff_offset_for_HC_grid_prim, &
                           Gcoeff_offset_for_HD_grid_prim

   implicit none

   integer, intent(in)  :: out_lun

   integer ier1,ier2,ier3,ier4
   integer :: off,gc,gd,gt,just_s,prob
   include "structure_factor_type.fh"
   include "mpole_exponent.fh"

   ! for MPOLE grids
   num_MC_grid_expon = 1
   MC_weight_grid_expon = 1.d0
   MC_grid_expon = compact_rec_mpole_expon
   
   num_MD_grid_expon = 1
   MD_weight_grid_expon = 1.d0
   MD_grid_expon = diffuse_rec_mpole_expon

   ! for sumCH grids
   num_sumCH_C_grid_expon = 1
   sumCH_C_weight_grid_expon = 1.d0
   sumCH_C_grid_expon = compact_rec_sumCH_expon
   
   num_sumCH_D_grid_expon = 1
   sumCH_D_weight_grid_expon = 1.d0
   sumCH_D_grid_expon = diffuse_rec_sumCH_expon

   ! compact hermites
   if ( tot_num_herm_Cprims > 0 )then
      if ( num_HC_prim_grids == 0 )then
         write(out_lun,*)'FH_PME_RECIP_mult_expon_grids:', &
                   ' tot_num_herm_Cprims > 0 but num_HC_prim_grids is 0' 
         stop
      endif
      if (allocated(num_expon_for_HC_grid))deallocate(num_expon_for_HC_grid)
      allocate(num_expon_for_HC_grid(num_HC_prim_grids),stat=ier1)
      if (allocated(expon_off_for_HC_grid))deallocate(expon_off_for_HC_grid)
      allocate(expon_off_for_HC_grid(num_HC_prim_grids),stat=ier2)
      if ( ier1 /= 0 .or. ier2 /= 0 )then
         write(out_lun,*)'FH_RECIP_pme_mult_expon_grids: allocate fails!'
         stop
      endif
      ! get the HC grid exponents
      off = 0
      do gc = 1,num_HC_prim_grids
         expon_off_for_HC_grid(gc) = off
         num_expon_for_HC_grid(gc) = 0
         gt = grid_type_for_HC_prim_grid(gc)
         if ( struc_fac_method_for_grid_type(gt) == SF_PME )then
            call FH_PME_RECIP_count_grid_expon( &
                        prim_offset_for_HC_grid(gc), &
                        num_prims_for_HC_grid(gc), &
                        herm_expon_for_HC_grid_prim, &
                        herm_order_for_HC_grid_prim, &
                        num_expon_for_HC_grid(gc), &
                        just_s )
            if ( (num_expon_for_HC_grid(gc) > 1) .and. (just_s == 0) )then
               write(out_lun,*)'FH_RECIP_setup_gridded_prims: ', &
                     'Bad pme grid: gc,num_expon_for_HC_grid(gc),just_s = ', &
                     gc,num_expon_for_HC_grid(gc),just_s
               stop
            endif
         endif
         off = off + num_expon_for_HC_grid(gc)
      enddo
      siz_HC_grid_expon = off
      ! allocate and fill
      if ( siz_HC_grid_expon > 0 )then
         allocate( HC_grid_expon(siz_HC_grid_expon),stat=ier1 )
         allocate( HC_weight_grid_expon(siz_HC_grid_expon),stat=ier2 )
         if ( ier1 /= 0 .or. ier2 /= 0 )then
            write(out_lun,*)'FH_RECIP_pme_mult_expon_grids: allocate fails!'
            stop
         endif
         HC_grid_expon = 0.d0
         HC_weight_grid_expon = 0.d0
      endif
      do gc = 1,num_HC_prim_grids
         gt = grid_type_for_HC_prim_grid(gc)
         if ( struc_fac_method_for_grid_type(gt) == SF_PME )then
            call FH_PME_RECIP_grid_expon_weight( &
                        prim_offset_for_HC_grid(gc), &
                        num_prims_for_HC_grid(gc), &
                        expon_off_for_HC_grid(gc), &
                        herm_expon_for_HC_grid_prim, &
                        Gcoeff_offset_for_HC_grid_prim, &
                        Local_Chermite_coeff, &
                        HC_grid_expon, &
                        HC_weight_grid_expon,prob)
            if ( prob == 1 )then
               write(out_lun,*)'FH_RECIP_pme_mult_expon_grids: gc = ',gc, &
                    'problem--,coeffs sum to near zero'
               stop
            endif
         endif !( struc_fac_method_for_grid_type(gt) == SF_PME )
      enddo !gc = 1,num_HC_prim_grids
   else
      siz_HC_grid_expon = 0
   endif !( tot_num_herm_Cprims > 0 )

   ! diffuse hermites
   if ( tot_num_herm_Dprims > 0 )then
      if ( num_HD_prim_grids == 0 )then
         write(out_lun,*)'FH_PME_RECIP_mult_expon_grids:', &
                   ' tot_num_herm_Dprims > 0 but num_HD_prim_grids is 0' 
         stop
      endif
      allocate(num_expon_for_HD_grid(num_HD_prim_grids),stat=ier1)
      allocate(expon_off_for_HD_grid(num_HD_prim_grids),stat=ier2)
      if ( ier1 /= 0 .or. ier2 /= 0 )then
         write(out_lun,*)'FH_RECIP_pme_mult_expon_grids: allocate fails!'
         stop
      endif
      ! get the HD grid exponents
      off = 0
      do gd = 1,num_HD_prim_grids
         expon_off_for_HD_grid(gd) = off
         num_expon_for_HD_grid(gd) = 0
         gt = grid_type_for_HD_prim_grid(gd)
         if ( struc_fac_method_for_grid_type(gt) == SF_PME )then
            call FH_PME_RECIP_count_grid_expon( &
                        prim_offset_for_HD_grid(gd), &
                        num_prims_for_HD_grid(gd), &
                        herm_expon_for_HD_grid_prim, &
                        herm_order_for_HD_grid_prim, &
                        num_expon_for_HD_grid(gd), &
                        just_s )
            if ( (num_expon_for_HD_grid(gd) > 1) .and. (just_s == 0) )then
               write(out_lun,*)'FH_RECIP_setup_gridded_prims: ', &
                     'Bad pme grid: gd,num_expon_for_HD_grid(gd),just_s = ', &
                     gd,num_expon_for_HD_grid(gd),just_s
               stop
            endif
         endif
         off = off + num_expon_for_HD_grid(gd)
      enddo
      siz_HD_grid_expon = off
      ! allocate and fill
      if ( siz_HD_grid_expon > 0 )then
         allocate( HD_grid_expon(siz_HD_grid_expon),stat=ier1 )
         allocate( HD_weight_grid_expon(siz_HD_grid_expon),stat=ier2 )
         if ( ier1 /= 0 .or. ier2 /= 0 )then
            write(out_lun,*)'FH_RECIP_pme_mult_expon_grids: allocate fails!'
            stop
         endif
         HD_grid_expon = 0.d0
         HD_weight_grid_expon = 0.d0
      endif
      do gd = 1,num_HD_prim_grids
         gt = grid_type_for_HD_prim_grid(gd)
         if ( struc_fac_method_for_grid_type(gt) == SF_PME )then
            call FH_PME_RECIP_grid_expon_weight( &
                        prim_offset_for_HD_grid(gd), &
                        num_prims_for_HD_grid(gd), &
                        expon_off_for_HD_grid(gd), &
                        herm_expon_for_HD_grid_prim, &
                        Gcoeff_offset_for_HD_grid_prim, &
                        Local_Dhermite_coeff, &
                        HD_grid_expon, &
                        HD_weight_grid_expon,prob)
            if ( prob == 1 )then
               write(out_lun,*)'FH_RECIP_pme_mult_expon_grids: gd = ',gd, &
                    'problem--,coeffs sum to near zero'
               stop
            endif
         endif !( struc_fac_method_for_grid_type(gt) == SF_PME )
      enddo !gd = 1,num_HD_prim_grids
   else
      siz_HD_grid_expon = 0
   endif !( tot_num_herm_Dprims > 0 )

end subroutine FH_PME_RECIP_mult_expon_grids
!-----------------------------------------------------------------
subroutine FH_PME_RECIP_setup_Bsplines(out_lun)

   use sites, only : num_sites
   use user,only : num_HC_prim_grids, &
                   num_HD_prim_grids, &
                   grid_type_for_HC_prim_grid, &
                   grid_type_for_HD_prim_grid, &
                   grid_type_for_MPOLES, &
                   grid_type_for_sumCH, &
                   num_prim_grid_types, &
                   struc_fac_method_for_grid_type, &
                   Bspline_order_for_grid_type, &
                   num_pme_grid_types
   use hermite,only : tot_num_mpole_coeffs, &
                      mpole_level_for_site, &
                      mpole_order_for_site, &
                      tot_num_sumCH_coeffs, &
                      sumCH_level_for_site, &
                      sumCH_order_for_site, &
                      tot_num_herm_Cprims, &
                      hermite_level_for_Cprim, &
                      site_that_owns_Cprim, &
                      hermite_order_for_Cprim, &
                      tot_num_herm_Dprims, &
                      hermite_level_for_Dprim, &
                      site_that_owns_Dprim, &
                      hermite_order_for_Dprim
   use recip_grids,only : prim_offset_for_HC_grid, &
                          prim_offset_for_HD_grid, &
                          num_prims_for_HC_grid, &
                          num_prims_for_HD_grid, &
                          prim_num_for_HC_grid_prim, &
                          prim_num_for_HD_grid_prim

   implicit none

   integer, intent(in)  :: out_lun

   include 'structure_factor_type.fh'

   integer :: ier1,ier2,ier3,ier4
   integer g,gt,off,num,n,np,ngp,deriv

   if ( num_pme_grid_types == 0 )return ! nothing further to do here

   allocate( Max_deriv_per_site_gridtype(num_sites,num_prim_grid_types), &
             stat=ier1)
   allocate( Spline_offset_for_site_gridtype(num_sites,num_prim_grid_types), &
             stat=ier2)
   if ( ier1 /= 0 .or. ier2 /= 0 )then
      write(out_lun,*)'FH_RECIP_setup_Bsplines: allocate fails!'
      stop
   endif

   ! setup the spline offsets and deriv val for 1st dimension
   Max_deriv_per_site_gridtype = 0
   
   ! multipoles
   if ( tot_num_mpole_coeffs > 0 )then
      gt = grid_type_for_MPOLES
      do n = 1,num_sites
         if ( mpole_order_for_site(n) > 0 )then
            deriv = mpole_level_for_site(n) + 1! add one for gradient
         else
            deriv = 0
         endif
         if ( deriv > Max_deriv_per_site_gridtype(n,gt) ) &
                      Max_deriv_per_site_gridtype(n,gt) = deriv
      enddo
   endif

   ! sumCH
   if ( tot_num_sumCH_coeffs > 0 )then
      gt = grid_type_for_sumCH
      do n = 1,num_sites
         if ( sumCH_order_for_site(n) > 0 )then
            deriv = sumCH_level_for_site(n) + 1! add one for gradient
         else
            deriv = 0
         endif
         if ( deriv > Max_deriv_per_site_gridtype(n,gt) ) &
                      Max_deriv_per_site_gridtype(n,gt) = deriv
      enddo
   endif

   ! compact hermites
   if ( tot_num_herm_Cprims > 0 )then
      do g = 1,num_HC_prim_grids
         gt = grid_type_for_HC_prim_grid(g)
         off = prim_offset_for_HC_grid(g)
         num = num_prims_for_HC_grid(g)
         do ngp = 1,num
            np = prim_num_for_HC_grid_prim(off+ngp)
            n = site_that_owns_Cprim(np)
            if ( hermite_order_for_Cprim(np) > 0 )then
               deriv = hermite_level_for_Cprim(np) + 1! add one for gradient
            else
               deriv = 0
            endif
            if ( deriv > Max_deriv_per_site_gridtype(n,gt) ) &
                         Max_deriv_per_site_gridtype(n,gt) = deriv
         enddo
      enddo
   endif !( tot_num_herm_Cprims > 0 )

   ! diffuse hermites
   if ( tot_num_herm_Dprims > 0 )then
      do g = 1,num_HD_prim_grids
         gt = grid_type_for_HD_prim_grid(g)
         off = prim_offset_for_HD_grid(g)
         num = num_prims_for_HD_grid(g)
         do ngp = 1,num
            np = prim_num_for_HD_grid_prim(off+ngp)
            n = site_that_owns_Dprim(np)
            if ( hermite_order_for_Dprim(np) > 0 )then
               deriv = hermite_level_for_Dprim(np) + 1! add one for gradient
            else
               deriv = 0
            endif
            if ( deriv > Max_deriv_per_site_gridtype(n,gt) ) &
                         Max_deriv_per_site_gridtype(n,gt) = deriv
         enddo
      enddo
   endif !( tot_num_herm_Dprims > 0 )

   ! get the spline offsets
   Spline_offset_for_site_gridtype = 0
   off = 0
   do gt = 1,num_prim_grid_types
      if ( struc_fac_method_for_grid_type(gt) == SF_PME )then !Bspline used
         do n = 1,num_sites
            Spline_offset_for_site_gridtype(n,gt) = off
            off = off + Bspline_order_for_grid_type(gt)* &
                (Max_deriv_per_site_gridtype(n,gt)+1) ! 1st dimension (0:deriv)
         enddo
      endif 
   enddo

   Bspline_table_size = off
   allocate(theta1(Bspline_table_size),stat=ier1)
   allocate(theta2(Bspline_table_size),stat=ier2)
   allocate(theta3(Bspline_table_size),stat=ier3)
   allocate(init_grid_ind_for_site_gridtype(3,num_sites,num_prim_grid_types),&
            stat=ier4)
   if ( ier1 /= 0 .or. ier2 /= 0 .or. ier3 /= 0 .or. ier4 /= 0 )then
      write(out_lun,*)'FH_RECIP_setup_Bsplines: allocate fails!'
      stop
   endif

end subroutine FH_PME_RECIP_setup_Bsplines
!-----------------------------------------------------------------
subroutine FH_PME_RECIP_grid_setup(out_lun)

   use user, only : grid_type_for_MPOLES, &
                    grid_type_for_sumCH, &
                    num_HC_prim_grids, &
                    grid_type_for_HC_prim_grid, &
                    num_HD_prim_grids, &
                    grid_type_for_HD_prim_grid, &
                    struc_fac_method_for_grid_type, &
                    nfft1_for_grid_type, &
                    nfft2_for_grid_type, &
                    nfft3_for_grid_type, &
                    Bspline_order_for_grid_type, &
                    num_pme_grid_types
   use recip_grids, only : nfftdim1_for_gridtype, &
                           nfftdim2_for_gridtype, &
                           nfftdim3_for_gridtype
   use hermite, only : tot_num_mpole_coeffs, &
                       tot_num_sumCH_coeffs, &
                       tot_num_herm_Cprims, &
                       tot_num_herm_Dprims

   implicit none

   integer, intent(in)  :: out_lun

   include 'structure_factor_type.fh'

   integer ier1,ier2,ier3,ier4,ier5
   integer off1,off2,nfft,g,gt

   if ( num_pme_grid_types == 0 )return ! nothing further to do here

   !multipoles
   if ( tot_num_mpole_coeffs > 0 )then
      gt = grid_type_for_MPOLES
      if ( struc_fac_method_for_grid_type(gt) == SF_PME )then!transform
         tot_siz_MPOLE_pme_mult = 2*nfftdim1_for_gridtype(gt)* &
                                    nfftdim2_for_gridtype(gt)* &
                                    nfftdim3_for_gridtype(gt)
         nfft = max(nfft1_for_grid_type(gt), &
                    nfft2_for_grid_type(gt), &
                    nfft3_for_grid_type(gt))
         tot_siz_MPOLE_prefac = 2*nfft ! prefac are complex not real
      else
         tot_siz_MPOLE_pme_mult = 0
         tot_siz_MPOLE_prefac = 0
      endif !( struc_fac_method_for_grid_type(gt) == SF_PME )
      if ( tot_siz_MPOLE_pme_mult > 0 )then
         allocate(MC_pme_multiplier(tot_siz_MPOLE_pme_mult),stat=ier1)
         allocate(MD_pme_multiplier(tot_siz_MPOLE_pme_mult),stat=ier2)
         allocate(MPOLE_prefac1(tot_siz_MPOLE_prefac),stat=ier3)
         allocate(MPOLE_prefac2(tot_siz_MPOLE_prefac),stat=ier4)
         allocate(MPOLE_prefac3(tot_siz_MPOLE_prefac),stat=ier5)
         if ( ier1 /= 0 .or. ier2 /= 0 .or. ier3 /= 0 .or. ier4 /= 0 &
             .or. ier5 /= 0 )then
            write(out_lun,*)'FH_RECIP_pme_grid_setup: allocate fails!'
            stop
         endif
         call FH_PME_RECIP_load_prefacs( &
                  nfft1_for_grid_type(gt), &
                  nfft2_for_grid_type(gt), &
                  nfft3_for_grid_type(gt), &
                  Bspline_order_for_grid_type(gt), &
                  MPOLE_prefac1,MPOLE_prefac2,MPOLE_prefac3)
      endif !( tot_siz_MPOLE_pme_mult > 0 )
   else
      tot_siz_MPOLE_pme_mult = 0
      tot_siz_MPOLE_prefac = 0
   endif !( tot_num_mpole_coeffs > 0 )

   !sumCH
   if ( tot_num_sumCH_coeffs > 0 )then
      gt = grid_type_for_sumCH
      if ( struc_fac_method_for_grid_type(gt) == SF_PME )then!transform
         tot_siz_sumCH_pme_mult = 2*nfftdim1_for_gridtype(gt)* &
                                    nfftdim2_for_gridtype(gt)* &
                                    nfftdim3_for_gridtype(gt)
         nfft = max(nfft1_for_grid_type(gt), &
                    nfft2_for_grid_type(gt), &
                    nfft3_for_grid_type(gt))
         tot_siz_sumCH_prefac = 2*nfft ! prefac are complex not real
      else
         tot_siz_sumCH_pme_mult = 0
         tot_siz_sumCH_prefac = 0
      endif !( struc_fac_method_for_grid_type(gt) == SF_PME )
      if ( tot_siz_sumCH_pme_mult > 0 )then
         allocate(sumCH_C_pme_multiplier(tot_siz_sumCH_pme_mult),stat=ier1)
         allocate(sumCH_D_pme_multiplier(tot_siz_sumCH_pme_mult),stat=ier2)
         allocate(sumCH_prefac1(tot_siz_sumCH_prefac),stat=ier3)
         allocate(sumCH_prefac2(tot_siz_sumCH_prefac),stat=ier4)
         allocate(sumCH_prefac3(tot_siz_sumCH_prefac),stat=ier5)
         if ( ier1 /= 0 .or. ier2 /= 0 .or. ier3 /= 0 .or. ier4 /= 0 &
             .or. ier5 /= 0 )then
            write(out_lun,*)'FH_RECIP_pme_grid_setup: allocate fails!'
            stop
         endif
         call FH_PME_RECIP_load_prefacs( &
                  nfft1_for_grid_type(gt), &
                  nfft2_for_grid_type(gt), &
                  nfft3_for_grid_type(gt), &
                  Bspline_order_for_grid_type(gt), &
                  sumCH_prefac1,sumCH_prefac2,sumCH_prefac3)
      endif !( tot_siz_sumCH_pme_mult > 0 )
   else
      tot_siz_sumCH_pme_mult = 0
      tot_siz_sumCH_prefac = 0
   endif !( tot_num_sumCH_coeffs > 0 )

   !compact hermites
   if ( tot_num_herm_Cprims > 0 )then
      allocate(HC_pme_mult_off_for_HC_grid(num_HC_prim_grids),stat=ier1)
      allocate(HC_prefac_off_for_HC_grid(num_HC_prim_grids),stat=ier2)
      if ( ier1 /= 0 .or. ier2 /= 0 )then
         write(out_lun,*)'FH_RECIP_pme_grid_setup: allocate fails!'
         stop
      endif
      off1 = 0
      off2 = 0
      do g = 1,num_HC_prim_grids
         HC_pme_mult_off_for_HC_grid(g) = off1
         HC_prefac_off_for_HC_grid(g) = off2
         gt = grid_type_for_HC_prim_grid(g)
         if ( struc_fac_method_for_grid_type(gt) == SF_PME )then!transform
            off1 = off1 + 2*nfftdim1_for_gridtype(gt)* &
                            nfftdim2_for_gridtype(gt)* &
                            nfftdim3_for_gridtype(gt)
            nfft = max(nfft1_for_grid_type(gt), &
                       nfft2_for_grid_type(gt), &
                       nfft3_for_grid_type(gt))
            off2 = off2 + 2*nfft ! prefac are complex not real
         endif !( struc_fac_method_for_grid_type(gt) == SF_PME )
      enddo !g = 1,num_HC_prim_grids
      tot_siz_HC_pme_mult = off1
      tot_siz_HC_prefac = off2
      if ( tot_siz_HC_pme_mult > 0 )then
         allocate(HC_pme_multiplier(tot_siz_HC_pme_mult),stat=ier1)
         allocate(HC_prefac1(tot_siz_HC_prefac),stat=ier2)
         allocate(HC_prefac2(tot_siz_HC_prefac),stat=ier3)
         allocate(HC_prefac3(tot_siz_HC_prefac),stat=ier4)
         if ( ier1 /= 0 .or. ier2 /= 0 .or. ier3 /= 0 .or. ier4 /= 0 )then
            write(out_lun,*)'FH_RECIP_pme_grid_setup: allocate fails!'
            stop
         endif
         do g = 1,num_HC_prim_grids
            off1 = HC_prefac_off_for_HC_grid(g)
            gt = grid_type_for_HC_prim_grid(g)
            if ( struc_fac_method_for_grid_type(gt) == SF_PME )then!transform
               call FH_PME_RECIP_load_prefacs( &
                        nfft1_for_grid_type(gt), &
                        nfft2_for_grid_type(gt), &
                        nfft3_for_grid_type(gt), &
                        Bspline_order_for_grid_type(gt), &
                        HC_prefac1(off1+1), &
                        HC_prefac2(off1+1), &
                        HC_prefac3(off1+1))
            endif !( struc_fac_method_for_grid_type(gt) == SF_PME )
         enddo !g = 1,num_HC_prim_grids
      endif !( tot_siz_HC_pme_mult > 0 )
   else
      tot_siz_HC_pme_mult = 0
      tot_siz_HC_prefac = 0
   endif !( tot_num_herm_Cprims > 0 )then

   !diffuse hermites
   if ( tot_num_herm_Dprims > 0 )then
      allocate(HD_pme_mult_off_for_HD_grid(num_HD_prim_grids),stat=ier1)
      allocate(HD_prefac_off_for_HD_grid(num_HD_prim_grids),stat=ier2)
      if ( ier1 /= 0 .or. ier2 /= 0 )then
         write(out_lun,*)'FH_RECIP_pme_grid_setup: allocate fails!'
         stop
      endif
      off1 = 0
      off2 = 0
      do g = 1,num_HD_prim_grids
         HD_pme_mult_off_for_HD_grid(g) = off1
         HD_prefac_off_for_HD_grid(g) = off2
         gt = grid_type_for_HD_prim_grid(g)
         if ( struc_fac_method_for_grid_type(gt) == SF_PME )then!transform
            off1 = off1 + 2*nfftdim1_for_gridtype(gt)* &
                            nfftdim2_for_gridtype(gt)* &
                            nfftdim3_for_gridtype(gt)
            nfft = max(nfft1_for_grid_type(gt), &
                       nfft2_for_grid_type(gt), &
                       nfft3_for_grid_type(gt))
            off2 = off2 + 2*nfft ! prefac are complex not real
         endif !( struc_fac_method_for_grid_type(gt) == SF_PME )
      enddo !g = 1,num_HD_prim_grids
      tot_siz_HD_pme_mult = off1
      tot_siz_HD_prefac = off2
      if ( tot_siz_HD_pme_mult > 0 )then
         allocate(HD_pme_multiplier(tot_siz_HD_pme_mult),stat=ier1)
         allocate(HD_prefac1(tot_siz_HD_prefac),stat=ier2)
         allocate(HD_prefac2(tot_siz_HD_prefac),stat=ier3)
         allocate(HD_prefac3(tot_siz_HD_prefac),stat=ier4)
         if ( ier1 /= 0 .or. ier2 /= 0 .or. ier3 /= 0 .or. ier4 /= 0 )then
            write(out_lun,*)'FH_RECIP_pme_grid_setup: allocate fails!'
            stop
         endif
         do g = 1,num_HD_prim_grids
            off1 = HD_prefac_off_for_HD_grid(g)
            gt = grid_type_for_HD_prim_grid(g)
            if ( struc_fac_method_for_grid_type(gt) == SF_PME )then!transform
               call FH_PME_RECIP_load_prefacs( &
                        nfft1_for_grid_type(gt), &
                        nfft2_for_grid_type(gt), &
                        nfft3_for_grid_type(gt), &
                        Bspline_order_for_grid_type(gt), &
                        HD_prefac1(off1+1), &
                        HD_prefac2(off1+1), &
                        HD_prefac3(off1+1))
            endif !( struc_fac_method_for_grid_type(gt) == SF_PME )
         enddo !g = 1,num_HD_prim_grids
      endif !( tot_siz_HD_pme_mult > 0 )
   else
      tot_siz_HD_pme_mult = 0
      tot_siz_HD_prefac = 0
   endif !( tot_num_herm_Dprims > 0 )then

end subroutine FH_PME_RECIP_grid_setup
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!-----------------------------------------------------------------
! EVALUATION SUBROUTINES
!-----------------------------------------------------------------
!-----------------------------------------------------------------

!-----------------------------------------------------------------
subroutine FH_PME_RECIP_global_to_frac()

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
                   num_pme_grid_types
   use unit_cell, only : recip
   use sites, only : num_sites
   use hermite, only : tot_num_mpole_coeffs, &
                       Global_multipole_coeff, &
                       tot_num_sumCH_coeffs, &
                       Global_sumCH_coeff, &
                       tot_num_herm_Cprims, &
                       Global_Chermite_coeff, &
                       tot_num_herm_Dprims, &
                       Global_Dhermite_coeff
   use recip_grids, only : num_mpoles_for_mpole_grid, &
                           Fcoeff_offset_for_gridded_mpole, &
                           Gcoeff_offset_for_gridded_mpole, &
                           Mpole_order_for_gridded_mpole, &
                           Mpole_grid_Fcoeff, &
                           ! sumch
                           num_sumCH_for_sumCH_grid, &
                           Fcoeff_offset_for_gridded_sumCH, &
                           Gcoeff_offset_for_gridded_sumCH, &
                           sumCH_order_for_gridded_sumCH, &
                           sumCH_grid_Fcoeff, &
                           ! compact hermites
                           num_prims_for_HC_grid, &
                           prim_offset_for_HC_grid, &
                           Fcoeff_offset_for_HC_grid_prim, &
                           Gcoeff_offset_for_HC_grid_prim, &
                           herm_order_for_HC_grid_prim, &
                           HC_grid_Fcoeff, &
                           ! diffuse hermites
                           num_prims_for_HD_grid, &
                           prim_offset_for_HD_grid, &
                           Fcoeff_offset_for_HD_grid_prim, &
                           Gcoeff_offset_for_HD_grid_prim, &
                           herm_order_for_HD_grid_prim, &
                           HD_grid_Fcoeff

   implicit none
                           
   include 'structure_factor_type.fh'

   integer :: g,gt,off

   if ( num_pme_grid_types == 0 )return ! nothing further to do here

   ! xform the multipoles--
   if ( tot_num_mpole_coeffs > 0 )then
      gt = grid_type_for_MPOLES
      if ( struc_fac_method_for_grid_type(gt) == SF_PME )then!transform
         off = 0 ! only one grid
         call FH_PME_RECIP_one_glob_to_frac(  &
                   nfft1_for_grid_type(gt), &
                   nfft2_for_grid_type(gt), &
                   nfft3_for_grid_type(gt), &
                   recip, &
                   num_mpoles_for_mpole_grid, &
                   off, &
                   Fcoeff_offset_for_gridded_mpole, &
                   Gcoeff_offset_for_gridded_mpole, &
                   Mpole_order_for_gridded_mpole, &
                   Global_multipole_coeff, &
                   Mpole_grid_Fcoeff)
      endif!( struc_fac_method_for_grid_type(gt) == SF_PME )
   endif !( tot_num_mpole_coeffs > 0 )then
   ! next xform the sumch
   if ( tot_num_sumCH_coeffs > 0 )then
      gt = grid_type_for_sumCH
      if ( struc_fac_method_for_grid_type(gt) == SF_PME )then!transform
         off = 0 ! only one grid
         call FH_PME_RECIP_one_glob_to_frac(  &
                   nfft1_for_grid_type(gt), &
                   nfft2_for_grid_type(gt), &
                   nfft3_for_grid_type(gt), &
                   recip, &
                   num_sumCH_for_sumCH_grid, &
                   off, &
                   Fcoeff_offset_for_gridded_sumCH, &
                   Gcoeff_offset_for_gridded_sumCH, &
                   sumCH_order_for_gridded_sumCH, &
                   Global_sumCH_coeff, &
                   sumCH_grid_Fcoeff)
      endif!( struc_fac_method_for_grid_type(gt) == SF_PME )
   endif !( tot_num_sumCH_coeffs > 0 )then

   ! next the compact herms
   if ( tot_num_herm_Cprims > 0 )then
      do g = 1,num_HC_prim_grids
         gt = grid_type_for_HC_prim_grid(g)
         if ( struc_fac_method_for_grid_type(gt) == SF_PME )then!transform
            call FH_PME_RECIP_one_glob_to_frac(  &
                   nfft1_for_grid_type(gt), &
                   nfft2_for_grid_type(gt), &
                   nfft3_for_grid_type(gt), &
                   recip, &
                   num_prims_for_HC_grid(g), &
                   prim_offset_for_HC_grid(g), &
                   Fcoeff_offset_for_HC_grid_prim, &
                   Gcoeff_offset_for_HC_grid_prim, &
                   herm_order_for_HC_grid_prim, &
                   Global_Chermite_coeff, &
                   HC_grid_Fcoeff)
         endif !( struc_fac_method_for_grid_type(gt) == SF_PME )then!transform
      enddo !g = 1,num_HC_prim_grids
   endif !( tot_num_herm_Cprims > 0 )then

   ! finally the diffuse herms
   if ( tot_num_herm_Dprims > 0 )then
      do g = 1,num_HD_prim_grids
         gt = grid_type_for_HD_prim_grid(g)
         if ( struc_fac_method_for_grid_type(gt) == SF_PME )then!transform
            call FH_PME_RECIP_one_glob_to_frac(  &
                   nfft1_for_grid_type(gt), &
                   nfft2_for_grid_type(gt), &
                   nfft3_for_grid_type(gt), &
                   recip, &
                   num_prims_for_HD_grid(g), &
                   prim_offset_for_HD_grid(g), &
                   Fcoeff_offset_for_HD_grid_prim, &
                   Gcoeff_offset_for_HD_grid_prim, &
                   herm_order_for_HD_grid_prim, &
                   Global_Dhermite_coeff, &
                   HD_grid_Fcoeff)
         endif !( struc_fac_method_for_grid_type(gt) == SF_PME )then!transform
      enddo !g = 1,num_HD_prim_grids
   endif !( tot_num_herm_Dprims > 0 )then

end subroutine FH_PME_RECIP_global_to_frac
!-----------------------------------------------------------------
subroutine FH_PME_RECIP_Bspline_fill(out_lun)

   use user,only : num_prim_grid_types, &
                   struc_fac_method_for_grid_type, &
                   Bspline_order_for_grid_type, &
                   nfft1_for_grid_type, &
                   nfft2_for_grid_type, &
                   nfft3_for_grid_type, &
                   num_pme_grid_types
   use unit_cell, only : recip
   use sites, only : num_sites,site_crd

   implicit none

   integer, intent(in) :: out_lun

   include 'structure_factor_type.fh'
   include 'mpole_sizes.fh'

   integer :: gt
   double precision scratch(MAXSplineOrder*MAXSplineOrder)

   if ( num_pme_grid_types == 0 )return ! nothing further to do here

   do gt = 1,num_prim_grid_types
      if ( struc_fac_method_for_grid_type(gt) == SF_PME )then !Bspline used
         call FH_PME_RECIP_one_Bspline_fill( &
                              site_crd, &
                              Max_deriv_per_site_gridtype(:,gt), &
                              Spline_offset_for_site_gridtype(:,gt), &
                              scratch, &
                              num_sites, &
                              nfft1_for_grid_type(gt), &
                              nfft2_for_grid_type(gt), &
                              nfft3_for_grid_type(gt), &
                              Bspline_order_for_grid_type(gt), &
                              recip, &
                              init_grid_ind_for_site_gridtype(:,:,gt), &
                              theta1, &
                              theta2, &
                              theta3, &
                              out_lun)
      endif
   enddo

end subroutine FH_PME_RECIP_Bspline_fill
!-----------------------------------------------------------------
subroutine FH_PME_RECIP_fill_grids()

   use user,only : num_HC_prim_grids, &
                   num_HD_prim_grids, &
                   grid_type_for_HC_prim_grid, &
                   grid_type_for_HD_prim_grid, &
                   grid_type_for_MPOLES, &
                   grid_type_for_sumCH, &
                   struc_fac_method_for_grid_type, &
                   Bspline_order_for_grid_type, &
                   nfft1_for_grid_type, &
                   nfft2_for_grid_type, &
                   nfft3_for_grid_type, &
                   num_pme_grid_types
   use sites, only : num_sites
   use hermite, only : tot_num_mpole_coeffs, &
                       tot_num_sumCH_coeffs, &
                       tot_num_herm_Cprims, &
                       tot_num_herm_Dprims
   use recip_grids, only : num_prims_for_HC_grid, &
                           prim_offset_for_HC_grid, &
                           HC_grid_Fcoeff, &
                           site_that_owns_HC_grid_prim, &
                           Fcoeff_offset_for_HC_grid_prim, &
                           herm_order_for_HC_grid_prim, &
                           HC_grid_offset_for_HC_grid,HC_grid, &
                           ! diffuse hermites
                           num_prims_for_HD_grid, &
                           prim_offset_for_HD_grid, &
                           HD_grid_Fcoeff, &
                           site_that_owns_HD_grid_prim, &
                           Fcoeff_offset_for_HD_grid_prim, &
                           herm_order_for_HD_grid_prim, &
                           HD_grid_offset_for_HD_grid,HD_grid, &
                           ! multipoles
                           num_mpoles_for_mpole_grid, &
                           Mpole_grid_Fcoeff,&
                           site_that_owns_gridded_mpole, &
                           Fcoeff_offset_for_gridded_mpole, &
                           Mpole_order_for_gridded_mpole, &
                           MC_grid, &
                           ! sum of compact hermites
                           num_sumCH_for_sumCH_grid, &
                           sumCH_grid_Fcoeff, &
                           site_that_owns_gridded_sumCH, &
                           Fcoeff_offset_for_gridded_sumCH, &
                           sumCH_order_for_gridded_sumCH, &
                           sumCH_C_grid, &
                           ! for general pme grids
                           nfftdim1_for_gridtype,nfftdim2_for_gridtype, &
                           nfftdim3_for_gridtype

   include 'structure_factor_type.fh'

   integer :: g,gt,off

   if ( num_pme_grid_types == 0 )return ! nothing further to do here

   ! first fill the compact multipole grid 
   ! copy to MD_grid in recip space and modify gaussian
   if ( tot_num_mpole_coeffs > 0 )then
      gt = grid_type_for_MPOLES
      if ( struc_fac_method_for_grid_type(gt) == SF_PME )then
         off = 0
         call FH_PME_RECIP_grid_Fcoeffs( &
                               num_sites, &
                               num_mpoles_for_mpole_grid, &
                               off, & !only one grid, offset = 0
                               Mpole_grid_Fcoeff,&
                               theta1,theta2,theta3,&
                               site_that_owns_gridded_mpole,&
                               init_grid_ind_for_site_gridtype(:,:,gt), &
                               Spline_offset_for_site_gridtype(:,gt), &
                               Fcoeff_offset_for_gridded_mpole, &
                               Mpole_order_for_gridded_mpole, &
                               Max_deriv_per_site_gridtype(:,gt), &
                               nfft1_for_grid_type(gt), &
                               nfft2_for_grid_type(gt), &
                               nfft3_for_grid_type(gt), &
                               Bspline_order_for_grid_type(gt), &
                               nfftdim1_for_gridtype(gt),&
                               nfftdim2_for_gridtype(gt), &
                               nfftdim2_for_gridtype(gt), &
                               MC_grid(off+1))
      endif !( struc_fac_method_for_grid_type(gt) == SF_PME )then
   endif !( tot_num_mpole_coeffs > 0 )then
   ! next fill the compact sumCH grid 
   ! copy to diffuse sumCH in recip space and modify gaussian
   if ( tot_num_sumCH_coeffs > 0 )then
      gt = grid_type_for_sumCH
      if ( struc_fac_method_for_grid_type(gt) == SF_PME )then
         off = 0
         call FH_PME_RECIP_grid_Fcoeffs( &
                               num_sites, &
                               num_sumCH_for_sumCH_grid, &
                               off, & !only one grid, offset = 0
                               sumCH_grid_Fcoeff,&
                               theta1,theta2,theta3,&
                               site_that_owns_gridded_sumCH,&
                               init_grid_ind_for_site_gridtype(:,:,gt), &
                               Spline_offset_for_site_gridtype(:,gt), &
                               Fcoeff_offset_for_gridded_sumCH, &
                               sumCH_order_for_gridded_sumCH, &
                               Max_deriv_per_site_gridtype(:,gt), &
                               nfft1_for_grid_type(gt), &
                               nfft2_for_grid_type(gt), &
                               nfft3_for_grid_type(gt), &
                               Bspline_order_for_grid_type(gt), &
                               nfftdim1_for_gridtype(gt),&
                               nfftdim2_for_gridtype(gt), &
                               nfftdim2_for_gridtype(gt), &
                               sumCH_C_grid(off+1))
      endif !( struc_fac_method_for_grid_type(gt) == SF_PME )then
   endif !( tot_num_sumCH_coeffs > 0 )then

   ! next fill the compact herm grid 
   if ( tot_num_herm_Cprims > 0 )then
      do g = 1,num_HC_prim_grids
         off = HC_grid_offset_for_HC_grid(g)
         gt = grid_type_for_HC_prim_grid(g)
         if ( struc_fac_method_for_grid_type(gt) == SF_PME )then
            call FH_PME_RECIP_grid_Fcoeffs( &
                               num_sites, &
                               num_prims_for_HC_grid(g), &
                               prim_offset_for_HC_grid(g), &
                               HC_grid_Fcoeff,&
                               theta1,theta2,theta3,&
                               site_that_owns_HC_grid_prim,&
                               init_grid_ind_for_site_gridtype(:,:,gt), &
                               Spline_offset_for_site_gridtype(:,gt), &
                               Fcoeff_offset_for_HC_grid_prim, &
                               herm_order_for_HC_grid_prim, &
                               Max_deriv_per_site_gridtype(:,gt), &
                               nfft1_for_grid_type(gt), &
                               nfft2_for_grid_type(gt), &
                               nfft3_for_grid_type(gt), &
                               Bspline_order_for_grid_type(gt), &
                               nfftdim1_for_gridtype(gt),&
                               nfftdim2_for_gridtype(gt), &
                               nfftdim2_for_gridtype(gt), &
                               HC_grid(off+1))
         endif !( struc_fac_method_for_grid_type(gt) == SF_PME )then
      enddo !g = 1,num_HC_prim_grids
   endif !( tot_num_herm_Cprims > 0 )then

   ! finally fill the diffuse herm grid 
   if ( tot_num_herm_Dprims > 0 )then
      do g = 1,num_HD_prim_grids
         off = HD_grid_offset_for_HD_grid(g)
         gt = grid_type_for_HD_prim_grid(g)
         if ( struc_fac_method_for_grid_type(gt) == SF_PME )then
            call FH_PME_RECIP_grid_Fcoeffs( &
                               num_sites, &
                               num_prims_for_HD_grid(g), &
                               prim_offset_for_HD_grid(g), &
                               HD_grid_Fcoeff,&
                               theta1,theta2,theta3,&
                               site_that_owns_HD_grid_prim,&
                               init_grid_ind_for_site_gridtype(:,:,gt), &
                               Spline_offset_for_site_gridtype(:,gt), &
                               Fcoeff_offset_for_HD_grid_prim, &
                               herm_order_for_HD_grid_prim, &
                               Max_deriv_per_site_gridtype(:,gt), &
                               nfft1_for_grid_type(gt), &
                               nfft2_for_grid_type(gt), &
                               nfft3_for_grid_type(gt), &
                               Bspline_order_for_grid_type(gt), &
                               nfftdim1_for_gridtype(gt),&
                               nfftdim2_for_gridtype(gt), &
                               nfftdim2_for_gridtype(gt), &
                               HD_grid(off+1))
         endif !( struc_fac_method_for_grid_type(gt) == SF_PME )then
      enddo !g = 1,num_HD_prim_grids
   endif !( tot_num_herm_Dprims > 0 )then

end subroutine FH_PME_RECIP_fill_grids
!-----------------------------------------------------------------
subroutine FH_PME_RECIP_grid_mult(out_lun)

   use user,only : num_HC_prim_grids, &
                   grid_type_for_HC_prim_grid, &
                   num_HD_prim_grids, &
                   grid_type_for_HD_prim_grid, &
                   grid_type_for_MPOLES, &
                   grid_type_for_sumCH, &
                   struc_fac_method_for_grid_type, &
                   num_pme_grid_types, &
                   nfft1_for_grid_type, &
                   nfft2_for_grid_type, &
                   nfft3_for_grid_type, &
                   orthogonal_ucell
   use recip_grids, only : nfftdim1_for_gridtype
   use unit_cell, only : recip,volume
   use hermite, only : tot_num_mpole_coeffs, &
                       tot_num_sumCH_coeffs, &
                       tot_num_herm_Cprims, &
                       tot_num_herm_Dprims

   implicit none

   integer, intent(in)  :: out_lun

   include 'structure_factor_type.fh'

   integer :: g,gt,off,off1,off2,num2

   if ( num_pme_grid_types == 0 )return ! nothing further to do here

   ! first the multipole grid multipliers
   if ( tot_num_mpole_coeffs > 0 )then
      gt = grid_type_for_MPOLES
      if ( struc_fac_method_for_grid_type(gt) == SF_PME )then
         ! first MC_pme_multiplier
         off = 0 !off_MC_pme_mult
         off1 = 0 ! off_MPOLE_prefac
         off2 = 0 !off_MC_grid_expon
         num2 = num_MC_grid_expon
         if ( orthogonal_ucell == 1 )then
            call FH_PME_RECIP_grid_mult_orthog( &
                     nfft1_for_grid_type(gt), &
                     nfft2_for_grid_type(gt), &
                     nfft3_for_grid_type(gt), &
                     nfftdim1_for_gridtype(gt), &
                     MC_grid_expon(off2+1), &
                     MC_weight_grid_expon(off2+1), &
                     num2, &
                     recip, &
                     volume, &
                     MPOLE_prefac1(off1+1), &
                     MPOLE_prefac2(off1+1), &
                     MPOLE_prefac3(off1+1), &
                     MC_pme_multiplier(off+1))
         else
            write(out_lun,*)&
              'FH_PME_RECIP_grid_mult:non-orthog ucells not ok yet'
            stop
         endif !( orthogonal_ucell == 1 )then
         ! next MD_pme_multiplier
         off = 0 !off_MD_pme_mult
         off1 = 0 ! off_MPOLE_prefac
         off2 = 0 !off_MD_grid_expon
         num2 = num_MD_grid_expon
         if ( orthogonal_ucell == 1 )then
            call FH_PME_RECIP_grid_mult_orthog( &
                     nfft1_for_grid_type(gt), &
                     nfft2_for_grid_type(gt), &
                     nfft3_for_grid_type(gt), &
                     nfftdim1_for_gridtype(gt), &
                     MD_grid_expon(off2+1), &
                     MD_weight_grid_expon(off2+1), &
                     num2, &
                     recip, &
                     volume, &
                     MPOLE_prefac1(off1+1), &
                     MPOLE_prefac2(off1+1), &
                     MPOLE_prefac3(off1+1), &
                     MD_pme_multiplier(off+1))
         else
            write(out_lun,*)&
              'FH_PME_RECIP_grid_mult:non-orthog ucells not ok yet'
            stop
         endif !( orthogonal_ucell == 1 )then
      endif !( struc_fac_method_for_grid_type(gt) == SF_PME )then
   endif !( tot_num_mpole_coeffs > 0 )then

   ! next the sumCH grid multipliers
   if ( tot_num_sumCH_coeffs > 0 )then
      gt = grid_type_for_sumCH
      if ( struc_fac_method_for_grid_type(gt) == SF_PME )then
         ! first sumCH_C_pme_multiplier
         off = 0 !off_sumCH_C_pme_mult
         off1 = 0 ! off_sumCH_prefac
         off2 = 0 !off_sumCH_C_grid_expon
         num2 = num_sumCH_C_grid_expon
         if ( orthogonal_ucell == 1 )then
            call FH_PME_RECIP_grid_mult_orthog( &
                     nfft1_for_grid_type(gt), &
                     nfft2_for_grid_type(gt), &
                     nfft3_for_grid_type(gt), &
                     nfftdim1_for_gridtype(gt), &
                     sumCH_C_grid_expon(off2+1), &
                     sumCH_C_weight_grid_expon(off2+1), &
                     num2, &
                     recip, &
                     volume, &
                     sumCH_prefac1(off1+1), &
                     sumCH_prefac2(off1+1), &
                     sumCH_prefac3(off1+1), &
                     sumCH_C_pme_multiplier(off+1))
         else
            write(out_lun,*)&
              'FH_PME_RECIP_grid_mult:non-orthog ucells not ok yet'
            stop
         endif !( orthogonal_ucell == 1 )then
         ! next MD_pme_multiplier
         off = 0 !off_sumCH_D_pme_mult
         off1 = 0 ! off_sumCH_prefac
         off2 = 0 !off_sumCH_D_grid_expon
         num2 = num_sumCH_D_grid_expon
         if ( orthogonal_ucell == 1 )then
            call FH_PME_RECIP_grid_mult_orthog( &
                     nfft1_for_grid_type(gt), &
                     nfft2_for_grid_type(gt), &
                     nfft3_for_grid_type(gt), &
                     nfftdim1_for_gridtype(gt), &
                     sumCH_D_grid_expon(off2+1), &
                     sumCH_D_weight_grid_expon(off2+1), &
                     num2, &
                     recip, &
                     volume, &
                     sumCH_prefac1(off1+1), &
                     sumCH_prefac2(off1+1), &
                     sumCH_prefac3(off1+1), &
                     sumCH_D_pme_multiplier(off+1))
         else
            write(out_lun,*)&
              'FH_PME_RECIP_grid_mult:non-orthog ucells not ok yet'
            stop
         endif !( orthogonal_ucell == 1 )then
      endif !( struc_fac_method_for_grid_type(gt) == SF_PME )then
   endif !( tot_num_sumCH_coeffs > 0 )then

   ! next the compact hermite grid multipliers
   if ( tot_num_herm_Cprims > 0 )then
      do g = 1,num_HC_prim_grids
         gt = grid_type_for_HC_prim_grid(g)
         if ( struc_fac_method_for_grid_type(gt) == SF_PME )then
            off = HC_pme_mult_off_for_HC_grid(g)
            off1 = HC_prefac_off_for_HC_grid(g)
            off2 = expon_off_for_HC_grid(g)
            num2 = num_expon_for_HC_grid(g)
            if ( orthogonal_ucell == 1 )then
               call FH_PME_RECIP_grid_mult_orthog( &
                        nfft1_for_grid_type(gt), &
                        nfft2_for_grid_type(gt), &
                        nfft3_for_grid_type(gt), &
                        nfftdim1_for_gridtype(gt), &
                        HC_grid_expon(off2+1), &
                        HC_weight_grid_expon(off2+1),num2, &
                        recip,volume, &
                        HC_prefac1(off1+1),HC_prefac2(off1+1), &
                        HC_prefac3(off1+1), &
                        HC_pme_multiplier(off+1))
            else
               write(out_lun,*)&
                 'FH_PME_RECIP_grid_mult:non-orthog ucells not ok yet'
               stop
            endif !( orthogonal_ucell == 1 )then
         endif !( struc_fac_method_for_grid_type(gt) == SF_PME )then
      enddo !g = 1,num_HC_prim_grids
   endif !( tot_num_herm_Cprims > 0 )then

   ! finally the diffuse hermite grid multipliers
   if ( tot_num_herm_Dprims > 0 )then
      do g = 1,num_HD_prim_grids
         gt = grid_type_for_HD_prim_grid(g)
         if ( struc_fac_method_for_grid_type(gt) == SF_PME )then
            off = HD_pme_mult_off_for_HD_grid(g)
            off1 = HD_prefac_off_for_HD_grid(g)
            off2 = expon_off_for_HD_grid(g)
            num2 = num_expon_for_HD_grid(g)
            if ( orthogonal_ucell == 1 )then
               call FH_PME_RECIP_grid_mult_orthog( &
                        nfft1_for_grid_type(gt), &
                        nfft2_for_grid_type(gt), &
                        nfft3_for_grid_type(gt), &
                        nfftdim1_for_gridtype(gt), &
                        HD_grid_expon(off2+1), &
                        HD_weight_grid_expon(off2+1),num2, &
                        recip,volume, &
                        HD_prefac1(off1+1),HD_prefac2(off1+1), &
                        HD_prefac3(off1+1), &
                        HD_pme_multiplier(off+1))
            else
               write(out_lun,*)&
                 'FH_PME_RECIP_grid_mult:non-orthog ucells not ok yet'
               stop
            endif !( orthogonal_ucell == 1 )then
         endif !( struc_fac_method_for_grid_type(gt) == SF_PME )then
      enddo !g = 1,num_HD_prim_grids
   endif !( tot_num_herm_Dprims > 0 )then

end subroutine FH_PME_RECIP_grid_mult
!-----------------------------------------------------------------
subroutine FH_PME_RECIP_FT_density()

   use user,only : num_HC_prim_grids, &
                   grid_type_for_HC_prim_grid, &
                   num_HD_prim_grids, &
                   grid_type_for_HD_prim_grid, &
                   grid_type_for_MPOLES, &
                   grid_type_for_sumCH, &
                   struc_fac_method_for_grid_type, &
                   nfft1_for_grid_type, &
                   nfft2_for_grid_type, &
                   nfft3_for_grid_type, &
                   num_pme_grid_types
   use recip_grids, only : HC_grid_offset_for_HC_grid, &
                           HD_grid_offset_for_HD_grid, &
                           nfftdim1_for_gridtype, &
                           tot_siz_MPOLE_grid, &
                           tot_siz_sumCH_grid, &
                           HC_grid, &
                           HD_grid, &
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

   integer :: g,gt,off1,off2

   if ( num_pme_grid_types == 0 )return ! nothing further to do here

   ! first the mpole grids  
   if ( tot_num_mpole_coeffs > 0 )then
      gt = grid_type_for_MPOLES
      if ( struc_fac_method_for_grid_type(gt) == SF_PME )then
         off1 = 0
         off2 = 0
         !First copy MC_grid to MD_grid
         call UTIL_copy_real_array(MC_grid,MD_grid,tot_siz_MPOLE_grid)
         ! next multiply by prefac and expon
         call FH_PME_RECIP_fix_one_FT_density( &
                     nfft1_for_grid_type(gt), &
                     nfft2_for_grid_type(gt), &
                     nfft3_for_grid_type(gt), &
                     nfftdim1_for_gridtype(gt), &
                     MC_pme_multiplier(off1+1), &
                     MC_grid(off2+1))
         call FH_PME_RECIP_fix_one_FT_density( &
                     nfft1_for_grid_type(gt), &
                     nfft2_for_grid_type(gt), &
                     nfft3_for_grid_type(gt), &
                     nfftdim1_for_gridtype(gt), &
                     MD_pme_multiplier(off1+1), &
                     MD_grid(off2+1))
      endif !( struc_fac_method_for_grid_type(gt) == SF_PME )then
   endif !( tot_num_mpole_coeffs > 0 )then

   ! next the sumCH grids
   if ( tot_num_sumCH_coeffs > 0 )then
      gt = grid_type_for_sumCH
      if ( struc_fac_method_for_grid_type(gt) == SF_PME )then
         off1 = 0
         off2 = 0
         !First copy sumCH_C_grid to sumCH_D_grid
         call UTIL_copy_real_array(sumCH_C_grid,sumCH_D_grid, &
                                   tot_siz_sumCH_grid)
         ! next multiply by prefac and expon
         call FH_PME_RECIP_fix_one_FT_density( &
                     nfft1_for_grid_type(gt), &
                     nfft2_for_grid_type(gt), &
                     nfft3_for_grid_type(gt), &
                     nfftdim1_for_gridtype(gt), &
                     sumCH_C_pme_multiplier(off1+1), &
                     sumCH_C_grid(off2+1))
         call FH_PME_RECIP_fix_one_FT_density( &
                     nfft1_for_grid_type(gt), &
                     nfft2_for_grid_type(gt), &
                     nfft3_for_grid_type(gt), &
                     nfftdim1_for_gridtype(gt), &
                     sumCH_D_pme_multiplier(off1+1), &
                     sumCH_D_grid(off2+1))
      endif !( struc_fac_method_for_grid_type(gt) == SF_PME )then
   endif !( tot_num_sumCH_coeffs > 0 )then

   ! next the compact hermites
   if ( tot_num_herm_Cprims > 0 )then
      do g = 1,num_HC_prim_grids
         gt = grid_type_for_HC_prim_grid(g)
         if ( struc_fac_method_for_grid_type(gt) == SF_PME )then
            off1 = HC_pme_mult_off_for_HC_grid(g)
            off2 = HC_grid_offset_for_HC_grid(g)
            call FH_PME_RECIP_fix_one_FT_density( &
                        nfft1_for_grid_type(gt), &
                        nfft2_for_grid_type(gt), &
                        nfft3_for_grid_type(gt), &
                        nfftdim1_for_gridtype(gt), &
                        HC_pme_multiplier(off1+1), &
                        HC_grid(off2+1))
         endif !( struc_fac_method_for_grid_type(gt) == SF_PME )then
      enddo !g = 1,num_HC_prim_grids
   endif !( tot_num_herm_Cprims > 0 )then

   ! finally the diffuse hermites
   if ( tot_num_herm_Dprims > 0 )then
      do g = 1,num_HD_prim_grids
         gt = grid_type_for_HD_prim_grid(g)
         if ( struc_fac_method_for_grid_type(gt) == SF_PME )then
            off1 = HD_pme_mult_off_for_HD_grid(g)
            off2 = HD_grid_offset_for_HD_grid(g)
            call FH_PME_RECIP_fix_one_FT_density( &
                        nfft1_for_grid_type(gt), &
                        nfft2_for_grid_type(gt), &
                        nfft3_for_grid_type(gt), &
                        nfftdim1_for_gridtype(gt), &
                        HD_pme_multiplier(off1+1), &
                        HD_grid(off2+1))
         endif !( struc_fac_method_for_grid_type(gt) == SF_PME )then
      enddo !g = 1,num_HD_prim_grids
   endif !( tot_num_herm_Dprims > 0 )then

end subroutine FH_PME_RECIP_FT_density
!-----------------------------------------------------------------
subroutine FH_PME_RECIP_fix_FT_phi()

   use user,only : num_HC_prim_grids, &
                   grid_type_for_HC_prim_grid, &
                   num_HD_prim_grids, &
                   grid_type_for_HD_prim_grid, &
                   grid_type_for_MPOLES, &
                   grid_type_for_sumCH, &
                   struc_fac_method_for_grid_type, &
                   nfft1_for_grid_type, &
                   nfft2_for_grid_type, &
                   nfft3_for_grid_type, &
                   num_pme_grid_types
   use unit_cell, only : volume
   use recip_grids, only : HC_grid_offset_for_HC_grid, &
                           HD_grid_offset_for_HD_grid, &
                           nfftdim1_for_gridtype, &
                           HC_grid, &
                           HD_grid, &
                           tot_siz_MPOLE_grid, &
                           MC_grid, &
                           MD_grid, &
                           tot_siz_sumCH_grid, &
                           sumCH_C_grid, &
                           sumCH_D_grid
   use hermite, only : tot_num_mpole_coeffs, &
                       tot_num_sumCH_coeffs, &
                       tot_num_herm_Cprims, &
                       tot_num_herm_Dprims

   implicit none

   include 'structure_factor_type.fh'

   integer :: g,gt,off1,off2

   if ( num_pme_grid_types == 0 )return ! nothing further to do here

   ! first mpole grids  
   if ( tot_num_mpole_coeffs > 0 )then
      gt = grid_type_for_MPOLES
      if ( struc_fac_method_for_grid_type(gt) == SF_PME )then
         off1 = 0
         off2 = 0
         !First copy MC_grid to MD_grid
         call UTIL_copy_real_array(MC_grid,MD_grid,tot_siz_MPOLE_grid)
         ! next multiply by prefac and expon
         call FH_PME_RECIP_fix_one_FT_phi( &
                     nfft1_for_grid_type(gt), &
                     nfft2_for_grid_type(gt), &
                     nfft3_for_grid_type(gt), &
                     nfftdim1_for_gridtype(gt), &
                     volume, &
                     MC_pme_multiplier(off1+1), &
                     MC_grid(off2+1))
         call FH_PME_RECIP_fix_one_FT_phi( &
                     nfft1_for_grid_type(gt), &
                     nfft2_for_grid_type(gt), &
                     nfft3_for_grid_type(gt), &
                     nfftdim1_for_gridtype(gt), &
                     volume, &
                     MD_pme_multiplier(off1+1), &
                     MD_grid(off2+1))
      endif !( struc_fac_method_for_grid_type(gt) == SF_PME )then
   endif !( tot_num_mpole_coeffs > 0 )then

   ! next sumCH grids  
   if ( tot_num_sumCH_coeffs > 0 )then
      gt = grid_type_for_sumCH
      if ( struc_fac_method_for_grid_type(gt) == SF_PME )then
         off1 = 0
         off2 = 0
         !First copy MC_grid to MD_grid
         call UTIL_copy_real_array(sumCH_C_grid,sumCH_D_grid, &
                                   tot_siz_sumCH_grid)
         ! next multiply by prefac and expon
         call FH_PME_RECIP_fix_one_FT_phi( &
                     nfft1_for_grid_type(gt), &
                     nfft2_for_grid_type(gt), &
                     nfft3_for_grid_type(gt), &
                     nfftdim1_for_gridtype(gt), &
                     volume, &
                     sumCH_C_pme_multiplier(off1+1), &
                     sumCH_C_grid(off2+1))
         call FH_PME_RECIP_fix_one_FT_phi( &
                     nfft1_for_grid_type(gt), &
                     nfft2_for_grid_type(gt), &
                     nfft3_for_grid_type(gt), &
                     nfftdim1_for_gridtype(gt), &
                     volume, &
                     sumCH_D_pme_multiplier(off1+1), &
                     sumCH_D_grid(off2+1))
      endif !( struc_fac_method_for_grid_type(gt) == SF_PME )then
   endif !( tot_num_mpole_coeffs > 0 )then

   ! next the compact hermites
   if ( tot_num_herm_Cprims > 0 )then
      do g = 1,num_HC_prim_grids
         gt = grid_type_for_HC_prim_grid(g)
         if ( struc_fac_method_for_grid_type(gt) == SF_PME )then
            off1 = HC_pme_mult_off_for_HC_grid(g)
            off2 = HC_grid_offset_for_HC_grid(g)
            call FH_PME_RECIP_fix_one_FT_phi( &
                        nfft1_for_grid_type(gt), &
                        nfft2_for_grid_type(gt), &
                        nfft3_for_grid_type(gt), &
                        nfftdim1_for_gridtype(gt), &
                        volume, &
                        HC_pme_multiplier(off1+1), &
                        HC_grid(off2+1))
         endif !( struc_fac_method_for_grid_type(gt) == SF_PME )then
      enddo !g = 1,num_HC_prim_grids
   endif !( tot_num_herm_Cprims > 0 )then

   ! finally the diffuse hermites
   if ( tot_num_herm_Dprims > 0 )then
      do g = 1,num_HD_prim_grids
         gt = grid_type_for_HD_prim_grid(g)
         if ( struc_fac_method_for_grid_type(gt) == SF_PME )then
            off1 = HD_pme_mult_off_for_HD_grid(g)
            off2 = HD_grid_offset_for_HD_grid(g)
            call FH_PME_RECIP_fix_one_FT_phi( &
                        nfft1_for_grid_type(gt), &
                        nfft2_for_grid_type(gt), &
                        nfft3_for_grid_type(gt), &
                        nfftdim1_for_gridtype(gt), &
                        volume, &
                        HD_pme_multiplier(off1+1), &
                        HD_grid(off2+1))
         endif !( struc_fac_method_for_grid_type(gt) == SF_PME )then
      enddo !g = 1,num_HD_prim_grids
   endif !( tot_num_herm_Dprims > 0 )then

end subroutine FH_PME_RECIP_fix_FT_phi
!-----------------------------------------------------------------
subroutine FH_PME_RECIP_ene_force_field(energy, out_lun)

   use user,only : num_HC_prim_grids, &
                   num_HD_prim_grids, &
                   grid_type_for_HC_prim_grid, &
                   grid_type_for_HD_prim_grid, &
                   grid_type_for_MPOLES, &
                   grid_type_for_sumCH, &
                   struc_fac_method_for_grid_type, &
                   num_pme_grid_types, &
                   nfft1_for_grid_type, &
                   nfft2_for_grid_type, &
                   nfft3_for_grid_type, &
                   Bspline_order_for_grid_type
   use unit_cell, only : recip
   use sites, only : num_sites,site_crd,site_frc
   use hermite, only : tot_num_mpole_coeffs, &
                       Global_multipole_field, &
                       tot_num_sumCH_coeffs, &
                       Global_sumCH_field, &
                       tot_num_herm_Cprims, &
                       Global_CHermite_field, &
                       tot_num_herm_Dprims, &
                       Global_DHermite_field
   use recip_grids, only :  &
                           ! mpoles
                           num_mpoles_for_mpole_grid, &
                           site_that_owns_gridded_mpole, &
                           Fcoeff_offset_for_gridded_mpole, &
                           Gcoeff_offset_for_gridded_mpole, &
                           Ffield_offset_for_gridded_mpole, &
                           Mpole_order_for_gridded_mpole, &
                           Ffield_order_for_gridded_mpole, &
                           Mpole_grid_Fcoeff, &
                           Mpole_grid_Ffield, &
                           MC_grid, &
                           ! sumCH
                           num_sumCH_for_sumCH_grid, &
                           site_that_owns_gridded_sumCH, &
                           Fcoeff_offset_for_gridded_sumCH, &
                           Gcoeff_offset_for_gridded_sumCH, &
                           Ffield_offset_for_gridded_sumCH, &
                           sumCH_order_for_gridded_sumCH, &
                           Ffield_order_for_gridded_sumCH, &
                           sumCH_grid_Fcoeff, &
                           sumCH_grid_Ffield, &
                           sumCH_C_grid, &
                           ! compact hermites
                           num_prims_for_HC_grid, &
                           site_that_owns_HC_grid_prim, &
                           HC_grid_offset_for_HC_grid, &
                           prim_offset_for_HC_grid, &
                           Fcoeff_offset_for_HC_grid_prim, &
                           Gcoeff_offset_for_HC_grid_prim, &
                           Ffield_offset_for_HC_grid_prim, &
                           herm_order_for_HC_grid_prim, &
                           Ffield_order_for_HC_grid_prim, &
                           HC_grid_Fcoeff, &
                           HC_grid_Ffield, &
                           HC_grid, &
                           ! diffuse hermites
                           num_prims_for_HD_grid, &
                           site_that_owns_HD_grid_prim, &
                           HD_grid_offset_for_HD_grid, &
                           prim_offset_for_HD_grid, &
                           Fcoeff_offset_for_HD_grid_prim, &
                           Gcoeff_offset_for_HD_grid_prim, &
                           Ffield_offset_for_HD_grid_prim, &
                           herm_order_for_HD_grid_prim, &
                           Ffield_order_for_HD_grid_prim, &
                           HD_grid_Fcoeff, &
                           HD_grid_Ffield, &
                           HD_grid, &
                           ! general
                           nfftdim1_for_gridtype, &
                           nfftdim2_for_gridtype, &
                           nfftdim3_for_gridtype

   implicit none

   double precision,intent(inout)       :: energy
   integer, intent(in)                  :: out_lun

   include 'structure_factor_type.fh'
   include 'interact_type.fh'

   integer :: g,gt,off,prim_offset_for_MC,prim_offset_for_sumCH_C
   double precision :: factor

   if ( num_pme_grid_types == 0 )return ! nothing further to do here

   !if ( energy_type == OVERLAP_ene_type )then
      !factor = exchange_factor
   !else
      factor = 1.d0
   !endif

   ! first the mpoles
   if ( tot_num_mpole_coeffs > 0 )then
      gt = grid_type_for_MPOLES
      if ( struc_fac_method_for_grid_type(gt) == SF_PME )then
         off = 0 ! only one grid
         prim_offset_for_MC = 0
         call FH_PME_RECIP_get_one_FHerm_phi( &
                             num_sites, &
                             num_mpoles_for_mpole_grid, &
                             prim_offset_for_MC, &
                             theta1,theta2,theta3, &
                             site_that_owns_gridded_mpole, &
                             init_grid_ind_for_site_gridtype(:,:,gt), &
                             Spline_offset_for_site_gridtype(:,gt), &
                             Ffield_offset_for_gridded_mpole, &
                             Ffield_order_for_gridded_mpole, &
                             Max_deriv_per_site_gridtype(:,gt), &
                             nfft1_for_grid_type(gt), &
                             nfft2_for_grid_type(gt), &
                             nfft3_for_grid_type(gt), &
                             Bspline_order_for_grid_type(gt), &
                             nfftdim1_for_gridtype(gt),&
                             nfftdim2_for_gridtype(gt), &
                             nfftdim2_for_gridtype(gt), &
                             MC_grid(off+1), &
                             Mpole_grid_Ffield,out_lun )
         call FH_PME_RECIP_one_ene_force( &
                             nfft1_for_grid_type(gt), &
                             nfft2_for_grid_type(gt), &
                             nfft3_for_grid_type(gt), &
                             recip, &
                             num_mpoles_for_mpole_grid, &
                             prim_offset_for_MC, &
                             site_that_owns_gridded_mpole, &
                             Fcoeff_offset_for_gridded_mpole, &
                             Mpole_order_for_gridded_mpole, &
                             Ffield_offset_for_gridded_mpole, &
                             Mpole_grid_Fcoeff, &
                             Mpole_grid_Ffield, &
                             factor, &
                             site_frc, &
                             energy)
         call FH_PME_RECIP_one_torque( &
                             nfft1_for_grid_type(gt), &
                             nfft2_for_grid_type(gt), &
                             nfft3_for_grid_type(gt), &
                             recip, &
                             num_mpoles_for_mpole_grid, &
                             prim_offset_for_MC, &
                             Mpole_order_for_gridded_mpole, &
                             Ffield_offset_for_gridded_mpole, &
                             Gcoeff_offset_for_gridded_mpole, &
                             factor, &
                             Mpole_grid_Ffield, &
                             Global_multipole_field )
      endif !( struc_fac_method_for_grid_type(gt) == SF_PME )
   endif !( tot_num_mpole_coeffs > 0 )then

   ! next the sumCH
   if ( tot_num_sumCH_coeffs > 0 )then
      gt = grid_type_for_sumCH
      if ( struc_fac_method_for_grid_type(gt) == SF_PME )then
         off = 0 ! only one grid
         prim_offset_for_sumCH_C = 0
         call FH_PME_RECIP_get_one_FHerm_phi( &
                             num_sites, &
                             num_sumCH_for_sumCH_grid, &
                             prim_offset_for_sumCH_C, &
                             theta1,theta2,theta3, &
                             site_that_owns_gridded_sumCH, &
                             init_grid_ind_for_site_gridtype(:,:,gt), &
                             Spline_offset_for_site_gridtype(:,gt), &
                             Ffield_offset_for_gridded_sumCH, &
                             Ffield_order_for_gridded_sumCH, &
                             Max_deriv_per_site_gridtype(:,gt), &
                             nfft1_for_grid_type(gt), &
                             nfft2_for_grid_type(gt), &
                             nfft3_for_grid_type(gt), &
                             Bspline_order_for_grid_type(gt), &
                             nfftdim1_for_gridtype(gt),&
                             nfftdim2_for_gridtype(gt), &
                             nfftdim2_for_gridtype(gt), &
                             sumCH_C_grid(off+1), &
                             sumCH_grid_Ffield,out_lun )
         call FH_PME_RECIP_one_ene_force( &
                             nfft1_for_grid_type(gt), &
                             nfft2_for_grid_type(gt), &
                             nfft3_for_grid_type(gt), &
                             recip, &
                             num_sumCH_for_sumCH_grid, &
                             prim_offset_for_sumCH_C, &
                             site_that_owns_gridded_sumCH, &
                             Fcoeff_offset_for_gridded_sumCH, &
                             sumCH_order_for_gridded_sumCH, &
                             Ffield_offset_for_gridded_sumCH, &
                             sumCH_grid_Fcoeff, &
                             sumCH_grid_Ffield, &
                             factor, &
                             site_frc, &
                             energy)
         call FH_PME_RECIP_one_torque( &
                             nfft1_for_grid_type(gt), &
                             nfft2_for_grid_type(gt), &
                             nfft3_for_grid_type(gt), &
                             recip, &
                             num_sumCH_for_sumCH_grid, &
                             prim_offset_for_MC, &
                             sumCH_order_for_gridded_sumCH, &
                             Ffield_offset_for_gridded_sumCH, &
                             Gcoeff_offset_for_gridded_sumCH, &
                             factor, &
                             sumCH_grid_Ffield, &
                             Global_sumCH_field )
      endif !( struc_fac_method_for_grid_type(gt) == SF_PME )
   endif !( tot_num_sumCH_coeffs > 0 )then

   ! next the hermite compacts
   if ( tot_num_herm_Cprims > 0 )then
      do g = 1,num_HC_prim_grids
         gt = grid_type_for_HC_prim_grid(g)
         if ( struc_fac_method_for_grid_type(gt) == SF_PME )then
            off = HC_grid_offset_for_HC_grid(g)
            call FH_PME_RECIP_get_one_FHerm_phi( &
                             num_sites, &
                             num_prims_for_HC_grid(g), &
                             prim_offset_for_HC_grid(g), &
                             theta1,theta2,theta3, &
                             site_that_owns_HC_grid_prim, &
                             init_grid_ind_for_site_gridtype(:,:,gt), &
                             Spline_offset_for_site_gridtype(:,gt), &
                             Ffield_offset_for_HC_grid_prim, &
                             Ffield_order_for_HC_grid_prim, &
                             Max_deriv_per_site_gridtype(:,gt), &
                             nfft1_for_grid_type(gt), &
                             nfft2_for_grid_type(gt), &
                             nfft3_for_grid_type(gt), &
                             Bspline_order_for_grid_type(gt), &
                             nfftdim1_for_gridtype(gt),&
                             nfftdim2_for_gridtype(gt), &
                             nfftdim2_for_gridtype(gt), &
                             HC_grid(off+1), &
                             HC_grid_Ffield,out_lun )
            call FH_PME_RECIP_one_ene_force( &
                             nfft1_for_grid_type(gt), &
                             nfft2_for_grid_type(gt), &
                             nfft3_for_grid_type(gt), &
                             recip, &
                             num_prims_for_HC_grid(g), &
                             prim_offset_for_HC_grid(g), &
                             site_that_owns_HC_grid_prim, &
                             Fcoeff_offset_for_HC_grid_prim, &
                             herm_order_for_HC_grid_prim, &
                             Ffield_offset_for_HC_grid_prim, &
                             HC_grid_Fcoeff, &
                             HC_grid_Ffield, &
                             factor, &
                             site_frc, & 
                             energy)
         call FH_PME_RECIP_one_torque( &
                             nfft1_for_grid_type(gt), &
                             nfft2_for_grid_type(gt), &
                             nfft3_for_grid_type(gt), &
                             recip, &
                             num_prims_for_HC_grid(g), &
                             prim_offset_for_HC_grid(g), &
                             herm_order_for_HC_grid_prim, &
                             Ffield_offset_for_HC_grid_prim, &
                             Gcoeff_offset_for_HC_grid_prim, &
                             factor, &
                             HC_grid_Ffield, &
                             Global_CHermite_field )
         endif !( struc_fac_method_for_grid_type(gt) == SF_PME )
      enddo !g = 1,num_HC_prim_grids
   endif !( tot_num_herm_Cprims > 0 )then

   ! finally the hermite diffuse
   if ( tot_num_herm_Dprims > 0 )then
      do g = 1,num_HD_prim_grids
         gt = grid_type_for_HD_prim_grid(g)
         if ( struc_fac_method_for_grid_type(gt) == SF_PME )then
            off = HD_grid_offset_for_HD_grid(g)
            call FH_PME_RECIP_get_one_FHerm_phi( &
                             num_sites, &
                             num_prims_for_HD_grid(g), &
                             prim_offset_for_HD_grid(g), &
                             theta1,theta2,theta3, &
                             site_that_owns_HD_grid_prim, &
                             init_grid_ind_for_site_gridtype(:,:,gt), &
                             Spline_offset_for_site_gridtype(:,gt), &
                             Ffield_offset_for_HD_grid_prim, &
                             Ffield_order_for_HD_grid_prim, &
                             Max_deriv_per_site_gridtype(:,gt), &
                             nfft1_for_grid_type(gt), &
                             nfft2_for_grid_type(gt), &
                             nfft3_for_grid_type(gt), &
                             Bspline_order_for_grid_type(gt), &
                             nfftdim1_for_gridtype(gt),&
                             nfftdim2_for_gridtype(gt), &
                             nfftdim2_for_gridtype(gt), &
                             HD_grid(off+1), &
                             HD_grid_Ffield,out_lun )
            call FH_PME_RECIP_one_ene_force( &
                             nfft1_for_grid_type(gt), &
                             nfft2_for_grid_type(gt), &
                             nfft3_for_grid_type(gt), &
                             recip, &
                             num_prims_for_HD_grid(g), &
                             prim_offset_for_HD_grid(g), &
                             site_that_owns_HD_grid_prim, &
                             Fcoeff_offset_for_HD_grid_prim, &
                             herm_order_for_HD_grid_prim, &
                             Ffield_offset_for_HD_grid_prim, &
                             HD_grid_Fcoeff, &
                             HD_grid_Ffield, &
                             factor, &
                             site_frc, & 
                             energy)
         call FH_PME_RECIP_one_torque( &
                             nfft1_for_grid_type(gt), &
                             nfft2_for_grid_type(gt), &
                             nfft3_for_grid_type(gt), &
                             recip, &
                             num_prims_for_HD_grid(g), &
                             prim_offset_for_HD_grid(g), &
                             herm_order_for_HD_grid_prim, &
                             Ffield_offset_for_HD_grid_prim, &
                             Gcoeff_offset_for_HD_grid_prim, &
                             factor, &
                             HD_grid_Ffield, &
                             Global_DHermite_field )
         endif !( struc_fac_method_for_grid_type(gt) == SF_PME )
      enddo !g = 1,num_HD_prim_grids
   endif !( tot_num_herm_Dprims > 0 )then

end subroutine FH_PME_RECIP_ene_force_field
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!-----------------------------------------------------------------
! DEALLOCATE SUBROUTINES
!-----------------------------------------------------------------
!-----------------------------------------------------------------

!-----------------------------------------------------------------
subroutine FH_PME_RECIP_deallocate()

   implicit none

   call FH_PME_RECIP_mult_expon_deall()
   call FH_PME_RECIP_bspline_deall()
   call FH_PME_RECIP_grid_deall()

end subroutine FH_PME_RECIP_deallocate
!-----------------------------------------------------------------
subroutine FH_PME_RECIP_mult_expon_deall()

   implicit none

   ! mpole
       ! nothing
   ! sumCH
       ! nothing
   ! compact hermites
   if ( allocated(num_expon_for_HC_grid) )deallocate(num_expon_for_HC_grid)
   if ( allocated(expon_off_for_HC_grid) )deallocate(expon_off_for_HC_grid)
   if ( allocated(HC_grid_expon) )deallocate(HC_grid_expon)
   if ( allocated(HC_weight_grid_expon) )deallocate(HC_weight_grid_expon)
   ! diffuse hermites
   if ( allocated(num_expon_for_HD_grid) )deallocate(num_expon_for_HD_grid)
   if ( allocated(expon_off_for_HD_grid) )deallocate(expon_off_for_HD_grid)
   if ( allocated(HD_grid_expon) )deallocate(HD_grid_expon)
   if ( allocated(HD_weight_grid_expon) )deallocate(HD_weight_grid_expon)
end subroutine FH_PME_RECIP_mult_expon_deall
!-----------------------------------------------------------------
subroutine FH_PME_RECIP_bspline_deall()

   implicit none

   if ( allocated(Max_deriv_per_site_gridtype) ) &
          deallocate(Max_deriv_per_site_gridtype)
   if ( allocated(Spline_offset_for_site_gridtype) ) &
          deallocate(Spline_offset_for_site_gridtype)
   if ( allocated(theta1) )deallocate(theta1)
   if ( allocated(theta2) )deallocate(theta2)
   if ( allocated(theta3) )deallocate(theta3)
   if ( allocated(init_grid_ind_for_site_gridtype) ) &
          deallocate(init_grid_ind_for_site_gridtype)

end subroutine FH_PME_RECIP_bspline_deall
!-----------------------------------------------------------------
subroutine FH_PME_RECIP_grid_deall()

   implicit none

   ! mpoles
   if ( allocated(MC_pme_multiplier) ) &
                    deallocate(MC_pme_multiplier)
   if ( allocated(MD_pme_multiplier) ) &
                    deallocate(MD_pme_multiplier)
   if ( allocated(MPOLE_prefac1) )deallocate(MPOLE_prefac1)
   if ( allocated(MPOLE_prefac2) )deallocate(MPOLE_prefac2)
   if ( allocated(MPOLE_prefac3) )deallocate(MPOLE_prefac3)

   ! sumCH
   if ( allocated(sumCH_C_pme_multiplier) ) &
                    deallocate(sumCH_C_pme_multiplier)
   if ( allocated(sumCH_D_pme_multiplier) ) &
                    deallocate(sumCH_D_pme_multiplier)
   if ( allocated(sumCH_prefac1) )deallocate(sumCH_prefac1)
   if ( allocated(sumCH_prefac2) )deallocate(sumCH_prefac2)
   if ( allocated(sumCH_prefac3) )deallocate(sumCH_prefac3)

   ! compact hermites
   if ( allocated(HC_pme_mult_off_for_HC_grid) ) &
                    deallocate(HC_pme_mult_off_for_HC_grid)
   if ( allocated(HC_prefac_off_for_HC_grid) ) &
                   deallocate(HC_prefac_off_for_HC_grid)
   if ( allocated(HC_pme_multiplier) ) &
                    deallocate(HC_pme_multiplier)
   if ( allocated(HC_prefac1) )deallocate(HC_prefac1)
   if ( allocated(HC_prefac2) )deallocate(HC_prefac2)
   if ( allocated(HC_prefac3) )deallocate(HC_prefac3)

   ! diffuse hermites
   if ( allocated(HD_pme_mult_off_for_HD_grid) ) &
                   deallocate(HD_pme_mult_off_for_HD_grid)
   if ( allocated(HD_prefac_off_for_HD_grid) ) &
                   deallocate(HD_prefac_off_for_HD_grid)
   if ( allocated(HD_pme_multiplier) ) &
                    deallocate(HD_pme_multiplier)
   if ( allocated(HD_prefac1) )deallocate(HD_prefac1)
   if ( allocated(HD_prefac2) )deallocate(HD_prefac2)
   if ( allocated(HD_prefac3) )deallocate(HD_prefac3)

end subroutine FH_PME_RECIP_grid_deall
!-----------------------------------------------------------------
!-----------------------------------------------------------------
end module pme_recip
!-----------------------------------------------------------------
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!-----------------------------------------------------------------
! NON-MODULE (LOCAL) SUBROUTINES
! THESE HANDLE ALL ARGUMENTS THROUGH ARGUMENT LIST
!-----------------------------------------------------------------
!-----------------------------------------------------------------

!-----------------------------------------------------------------
subroutine FH_PME_RECIP_count_grid_expon(offG,numG,exponGH,orderGH, &
                                     num_diff_expon,just_s)

   implicit none

   integer,intent(in) :: offG,numG
   double precision, intent(in) :: exponGH(*)
   integer,intent(in) :: orderGH(*)
   integer,intent(out) :: num_diff_expon,just_s

   double precision :: buf(numG),tol,expo
   integer :: kp,np,nw,numnew,new
   buf = 0.d0
   tol = 1.d-10
   numnew = 0
   just_s = 1
   do kp = 1,numG
      np = offG + kp
      expo = exponGH(np)
      if ( orderGH(np) > 1 )just_s = 0 
      new = 1
      do nw = 1,numnew
         if ( dabs(expo - buf(nw)) < tol )then
            new = 0
            exit
         endif 
      enddo
      if ( new == 1 )then
         numnew = numnew + 1
         buf(numnew) = expo
      endif
   enddo
   num_diff_expon = numnew
end subroutine FH_PME_RECIP_count_grid_expon
!-----------------------------------------------------------------
subroutine FH_PME_RECIP_grid_expon_weight(offG,numG,off_expon_grid, &
                                       exponGH,coeff_offGH,LHerm_coeff, &
                                       expon_grid,weight_expon_grid,prob)

   implicit none

   integer,intent(in) :: offG,numG,off_expon_grid
   double precision, intent(in) :: exponGH(*)
   integer,intent(in) :: coeff_offGH(*)
   double precision, intent(in) :: LHerm_coeff(*)
   double precision, intent(out) :: expon_grid(*),weight_expon_grid(*)
   integer,intent(out) :: prob
   
   double precision :: expbuf(numG),coefbuf(numG),tol,expo,coef,sum_coef
   integer :: kp,np,nw,numnew,new,off_coef

   expbuf = 0.d0
   coefbuf = 0.d0
   tol = 1.d-10
   numnew = 0
   do kp = 1,numG
      np = offG + kp
      expo = exponGH(np)
      off_coef = coeff_offGH(np)
      coef = LHerm_coeff(off_coef + 1)
      new = 1
      do nw = 1,numnew
         if ( dabs(expo - expbuf(nw)) < tol )then
            new = 0
            exit
         endif 
      enddo
      if ( new == 1 )then
         numnew = numnew + 1
         expbuf(numnew) = expo
         coefbuf(numnew) = coef
      endif
   enddo
   sum_coef = 0.d0
   do nw = 1,numnew
      expon_grid(off_expon_grid+nw) = expbuf(nw)
   enddo
   do nw = 1,numnew
      sum_coef = sum_coef + coefbuf(nw)
   enddo
   if ( dabs(sum_coef) > tol )then
      prob = 0
      do nw = 1,numnew
         weight_expon_grid(off_expon_grid+nw) = coefbuf(nw) / sum_coef
      enddo
   else
      prob = 1
   endif
end subroutine FH_PME_RECIP_grid_expon_weight
!-----------------------------------------------------------------
subroutine FH_PME_RECIP_load_prefacs(nfft1,nfft2,nfft3,order, &
                     prefac1,prefac2,prefac3)

   implicit none

   integer,intent(in) :: nfft1,nfft2,nfft3,order
   double precision,intent(out) :: prefac1(2,nfft1),prefac2(2,nfft2), &
                                   prefac3(2,nfft3)

   call FH_PME_RECIP_bsp_interp_coeff(order,nfft1,prefac1)
   call FH_PME_RECIP_mult_by_lambda(order,nfft1,prefac1)

   call FH_PME_RECIP_bsp_interp_coeff(order,nfft2,prefac2)
   call FH_PME_RECIP_mult_by_lambda(order,nfft2,prefac2)

   call FH_PME_RECIP_bsp_interp_coeff(order,nfft3,prefac3)
   call FH_PME_RECIP_mult_by_lambda(order,nfft3,prefac3)

end subroutine FH_PME_RECIP_load_prefacs
!----------------------------------------------------------------------
subroutine FH_PME_RECIP_bsp_interp_coeff(order,nfft,coeff_array)

   implicit none

   integer, intent(in) :: order,nfft
   double precision,intent(out) :: coeff_array(2,nfft)

   integer j,k
   double precision :: bsp_array(0:order)
   double precision :: sumr,sumi,pi,twopi,arg,denom,eps
   include 'scale.fh'
   
   pi = 3.14159265358979323846d0
   twopi = 2.d0*pi
   ! fill bspline at integer values
   call FH_PME_RECIP_bspline_integers(order,bsp_array)
   do k = 1,nfft
      sumr = 0.d0
      sumi = 0.d0
      do j = 1,order-1
         ! modify since we are taking negative of mbar,where mbar = (k-1)/nfft
         ! due to nature of fft routines being backwards
         !arg = -twopi*dble(k-1)*dble(j) / nfft 
         arg = twopi*dble(k-1)*dble(j) / nfft 
         sumr = sumr + bsp_array(j)*cos(arg)
         sumi = sumi + bsp_array(j)*sin(arg)
      enddo
      ! invert sum to get coeff to force interpolation
      denom = sumr**2 + sumi**2
      if ( denom < small_value )then
         coeff_array(1,k) = 0.d0
         coeff_array(2,k) = 0.d0
      else
         coeff_array(1,k) = sumr / denom
         coeff_array(2,k) = -sumi / denom
      endif
   enddo 
end subroutine FH_PME_RECIP_bsp_interp_coeff
!----------------------------------------------------------------------
subroutine FH_PME_RECIP_bspline_integers(order,array)

   implicit none

   integer, intent(in) :: order
   double precision, intent(out) :: array(0:order)

   integer :: k
   ! order 2
   array(0) = 0.d0
   array(1) = 1.d0
   array(2) = 0.d0
   do k = 3,order
      call FH_PME_RECIP_bspline_int_rec(k,array)
   enddo
end subroutine FH_PME_RECIP_bspline_integers
!----------------------------------------------------------------------
subroutine FH_PME_RECIP_bspline_int_rec(order,array)

   implicit none

   integer, intent(in) :: order
   double precision, intent(inout) :: array(0:order)

   double precision :: temp(0:order-1)
   integer k

   do k = 0,order-1
      temp(k) = array(k)
   enddo
   array(0) = 0.d0
   array(order) = 0.d0
   do k = 1,order-1
      array(k) = (dble(k)/(order-1))*temp(k) + &
                 (dble(order-k)/(order-1))*temp(k-1)
   enddo
end subroutine FH_PME_RECIP_bspline_int_rec
!----------------------------------------------------------------------
subroutine FH_PME_RECIP_mult_by_lambda(order,nfft,prefac)

   implicit none

   integer,intent(in) :: order,nfft
   double precision, intent(inout) :: prefac(2,nfft)

   integer :: k,m,nf,order2
   double precision :: tol,lambda,gsum,gsum2
   tol = 1.d-8
   nf = nfft / 2
   if ( 2*nf < nfft )nf = nf+1
   order2 = 2*order

   do k = 1,nfft
      m = k - 1
      if ( k > nf )m = k - 1 - nfft
      if ( m == 0 )then
         lambda = 1.d0
      else
         call FH_PME_RECIP_gamma_sum(m,nfft,order,tol,gsum)
         call FH_PME_RECIP_gamma_sum(m,nfft,order2,tol,gsum2)
         lambda = gsum / gsum2
      endif
      prefac(1,k) = lambda*prefac(1,k)
      prefac(2,k) = lambda*prefac(2,k)
   enddo

end subroutine FH_PME_RECIP_mult_by_lambda
!----------------------------------------------------------------------
subroutine FH_PME_RECIP_gamma_sum(m,nfft,order,tol,gsum)

   implicit none

   integer,intent(in) :: m,nfft,order
   double precision, intent(in) :: tol
   double precision, intent(out) :: gsum

   double precision :: term,mbar
   integer :: k

   if ( m == 0 )then
      gsum = 1.d0
      return
   endif
   ! remember we are using back transform---exp(-2pi*i*mbar*u)
   mbar = -dble(m)/nfft
   ! k = 0 term
   gsum = 0.d0
   term = 1.d0
   gsum = gsum + term
   do k = 1,99999
      term = (mbar/(mbar+k))**order
      gsum = gsum + term
      if ( abs(term) < tol )exit
   enddo
   do k = 1,99999
      term = (mbar/(mbar-k))**order
      gsum = gsum + term
      if ( abs(term) < tol )exit
   enddo

end subroutine FH_PME_RECIP_gamma_sum
!----------------------------------------------------------------------
subroutine FH_PME_RECIP_one_glob_to_frac( &
                   nfft1, &
                   nfft2, &
                   nfft3, &
                   recip, &
                   num_prims_for_grid, &
                   prim_offset_for_grid, &
                   Fcoeff_offset_for_grid_prim, &
                   Gcoeff_offset_for_grid_prim, &
                   hermite_order_for_grid_prim, &
                   Global_hermite_coeffs, &
                   Frac_hermite_coeffs)

   implicit none

   integer, intent(in) :: nfft1,nfft2,nfft3
   double precision, intent(in) :: recip(3,3)
   integer, intent(in) :: num_prims_for_grid,prim_offset_for_grid
   integer, intent(in) :: Fcoeff_offset_for_grid_prim(*), &
                          Gcoeff_offset_for_grid_prim(*), &
                          hermite_order_for_grid_prim(*)
   double precision, intent(in) :: Global_hermite_coeffs(*)
   double precision, intent(out) :: Frac_hermite_coeffs(*)
   include "mpole_sizes.fh"

   double precision mpole_xform_3x3(3,3),field_xform_3x3(3,3)
   double precision Matrix_xyz_F(MAXMP,MAXMP)
   integer j,n,ngp,dimxyz_F,order,off_F,off_G

!  first get mpole_xform_3x3
   call FH_PME_RECIP_jacobian(nfft1,nfft2,nfft3,recip, &
                             mpole_xform_3x3,field_xform_3x3)
   call XFORM_MPOLE_matrix(mpole_xform_3x3,Matrix_xyz_F,MAXMP)

   dimxyz_F = MAXMP
   do n = 1,num_prims_for_grid
      ngp = prim_offset_for_grid + n
      order = hermite_order_for_grid_prim(ngp)
      if ( order > 0 )then
         off_F = Fcoeff_offset_for_grid_prim(ngp)
         off_G = Gcoeff_offset_for_grid_prim(ngp)
         call XFORM_MPOLE(Matrix_xyz_F,dimxyz_F, &
                             Global_hermite_coeffs(off_G+1), &
                             Frac_hermite_coeffs(off_F+1), order)
      endif
   enddo
   
end subroutine FH_PME_RECIP_one_glob_to_frac
!---------------------------------------------------------------
subroutine FH_PME_RECIP_jacobian(nfft1,nfft2,nfft3,recip, &
                      mpole_xform_3x3,field_xform_3x3)

  implicit none

  integer nfft1,nfft2,nfft3
  double precision recip(3,3)
  double precision mpole_xform_3x3(3,3),field_xform_3x3(3,3)

  double precision du1_dx,du1_dy,du1_dz,du2_dx,du2_dy,du2_dz, &
                   du3_dx,du3_dy,du3_dz

  ! get the matrix for change of variables  scaled frac in terms of cartesian
  du1_dx = nfft1*recip(1,1)
  du1_dy = nfft1*recip(2,1)
  du1_dz = nfft1*recip(3,1)
  du2_dx = nfft2*recip(1,2)
  du2_dy = nfft2*recip(2,2)
  du2_dz = nfft2*recip(3,2)
  du3_dx = nfft3*recip(1,3)
  du3_dy = nfft3*recip(2,3)
  du3_dz = nfft3*recip(3,3)

  field_xform_3x3(1,1) = du1_dx
  field_xform_3x3(1,2) = du2_dx
  field_xform_3x3(1,3) = du3_dx
  field_xform_3x3(2,1) = du1_dy
  field_xform_3x3(2,2) = du2_dy
  field_xform_3x3(2,3) = du3_dy
  field_xform_3x3(3,1) = du1_dz
  field_xform_3x3(3,2) = du2_dz
  field_xform_3x3(3,3) = du3_dz

  mpole_xform_3x3(1,1) = du1_dx
  mpole_xform_3x3(1,2) = du1_dy
  mpole_xform_3x3(1,3) = du1_dz
  mpole_xform_3x3(2,1) = du2_dx
  mpole_xform_3x3(2,2) = du2_dy
  mpole_xform_3x3(2,3) = du2_dz
  mpole_xform_3x3(3,1) = du3_dx
  mpole_xform_3x3(3,2) = du3_dy
  mpole_xform_3x3(3,3) = du3_dz
end subroutine FH_PME_RECIP_jacobian
!--------------------------------------------------------------
subroutine FH_PME_RECIP_one_Bspline_fill(  &
                              sitecrd, &
                              Deriv_order, &
                              Spline_offset, &
                              scratch,  &
                              nsites, &
                              nfft1, &
                              nfft2, &
                              nfft3, &
                              Spline_order, &
                              recip, &
                              init_grid_ind, &
                              theta1, &
                              theta2, &
                              theta3, &
                              out_lun)

   implicit none

   double precision,intent(in) :: sitecrd(3,*)
   integer, intent(in) :: Deriv_order(*),Spline_offset(*)
   integer,intent(in) :: nsites,nfft1,nfft2,nfft3,Spline_order
   double precision,intent(in) :: recip(3,3)
   integer,intent(out) :: init_grid_ind(3,nsites)
   double precision,intent(out) :: scratch(*),theta1(*),theta2(*),theta3(*)
   integer, intent(in)  :: out_lun

  integer n,off,dr_ord,ifr
  double precision w,fr
  do n = 1,nsites
    dr_ord = Deriv_order(n)
    if ( dr_ord > 0 )then
      off = Spline_offset(n) + 1
      w = sitecrd(1,n)*recip(1,1)+sitecrd(2,n)*recip(2,1)+ &
          sitecrd(3,n)*recip(3,1)
      fr = nfft1*(w - dnint(w) + 0.5d0)
      ifr = int(fr)
      w = fr - ifr
      init_grid_ind(1,n) = ifr - Spline_order
      call FH_PME_RECIP_bspline_fill_gen(w,Spline_order,scratch,dr_ord, &
                                         theta1(off),out_lun)
      w = sitecrd(1,n)*recip(1,2)+sitecrd(2,n)*recip(2,2)+ &
          sitecrd(3,n)*recip(3,2)
      fr = nfft2*(w - dnint(w) + 0.5d0)
      ifr = int(fr)
      w = fr - ifr
      init_grid_ind(2,n) = ifr - Spline_order
      call FH_PME_RECIP_bspline_fill_gen(w,Spline_order,scratch,dr_ord, &
                                         theta2(off),out_lun)
      w = sitecrd(1,n)*recip(1,3)+sitecrd(2,n)*recip(2,3)+ &
          sitecrd(3,n)*recip(3,3)
      fr = nfft3*(w - dnint(w) + 0.5d0)
      ifr = int(fr)
      w = fr - ifr
      init_grid_ind(3,n) = ifr - Spline_order
      call FH_PME_RECIP_bspline_fill_gen(w,Spline_order,scratch,dr_ord, &
                                         theta3(off),out_lun)
    endif
  enddo
  return
end subroutine FH_PME_RECIP_one_Bspline_fill
!--------------------------------------------------------------------
subroutine FH_PME_RECIP_bspline_fill_gen(w, order, array, dr_order, new_array, &
                                         out_lun)

  implicit none

  integer, intent(in)   :: out_lun

  integer order,dr_order
  double precision w,array(order,order),new_array(dr_order+1,order)

  integer k,j

! init order 2
  array(2,2) = w      
  array(1,2) = 1.d0 - w      
! one pass to order 3
  array(3,3) = 0.5d0*w*array(2,2)
  array(2,3) = 0.5d0*((w + 1.d0)*array(1,2)+(2.d0-w)*array(2,2))
  array(1,3) = 0.5d0*(1.d0-w)*array(1,2)
! compute standard b-spline recursion
  do k = 4,order
    call FH_PME_RECIP_bspline_recur(w,k,array(1,k-1),array(1,k))
  enddo
! do derivatives
  if ( dr_order > 0 )then
    call FH_PME_RECIP_bspline_diff(array(1,order-1),order)
    if ( dr_order > 1 )then
      call FH_PME_RECIP_bspline_diff(array(1,order-2),order-1)
      call FH_PME_RECIP_bspline_diff(array(1,order-2),order)
      if ( dr_order > 2 )then
        call FH_PME_RECIP_bspline_diff(array(1,order-3),order-2)
        call FH_PME_RECIP_bspline_diff(array(1,order-3),order-1)
        call FH_PME_RECIP_bspline_diff(array(1,order-3),order)
        if ( dr_order > 3 )then
          call FH_PME_RECIP_bspline_diff(array(1,order-4),order-3)
          call FH_PME_RECIP_bspline_diff(array(1,order-4),order-2)
          call FH_PME_RECIP_bspline_diff(array(1,order-4),order-1)
          call FH_PME_RECIP_bspline_diff(array(1,order-4),order)
          if ( dr_order > 4 )then
            call FH_PME_RECIP_bspline_diff(array(1,order-5),order-4)
            call FH_PME_RECIP_bspline_diff(array(1,order-5),order-3)
            call FH_PME_RECIP_bspline_diff(array(1,order-5),order-2)
            call FH_PME_RECIP_bspline_diff(array(1,order-5),order-1)
            call FH_PME_RECIP_bspline_diff(array(1,order-5),order)
            if ( dr_order > 5 )then
              write(out_lun,*)'derivs of order > 5 not implemented!'
              stop
            endif !( dr_order > 5 )then
          endif !( dr_order > 4 )then
        endif !( dr_order > 3 )then
      endif !( dr_order > 2 )then
    endif !( dr_order > 1 )then
  endif !( dr_order > 0 )then
! re-arrange array
  do k = 1,order
    do j = 1,dr_order+1
      new_array(j,k) = array(k,order-j+1)
    enddo
  enddo
  return
end subroutine FH_PME_RECIP_bspline_fill_gen
!---------------------------------------------------
subroutine FH_PME_RECIP_bspline_recur(w,n,old,new)

  implicit none

  double precision old(*),new(*),w
  integer n
! Using notation from Essmann et al; w = u-[u] and
! array(j) = M_n(w + order - j)  where n is order
! RECURSION:  M_n(w) = (w/(n-1))*M_n-1(w)+((n-w)/(n-1))*M_n-1(w-1)
! i.e.   M_n(w+n-j) = ((w+n-j)/(n-1))*M_n-1(w+n-j)+((j-w)/(n-1))*M_n-1(w+n-j-1)
! i.e.   new(j) = ((w+n-j)/(n-1))*old(j-1) + ((j-w)/(n-1))*old(j)
! where old is array before one_pass (thus n->n-1) and new is array afterwards

  double precision div
  integer j

  div = 1.d0 / (n-1)
  new(n) = div*w*old(n-1)
  do j = 1,n-2
    new(n-j) = div*((w+j)*old(n-j-1) + (n-j-w)*old(n-j))
  enddo
  new(1) = div*(1-w)*old(1)
  return
end subroutine FH_PME_RECIP_bspline_recur
!-------------------------------------------------------------
subroutine FH_PME_RECIP_bspline_diff(c,n)

  implicit none

  double precision c(*)
  integer n
! Using notation from Essmann et al; w = u-[u] and
! array(j) = M_n(w + order - j)  where n is order
! DERIVATIVE:    d/dw M_n(w) = M_n-1(w) - M_n-1(w-1)
! i.e.   d/dw M_n(w+n-j) = M_n-1(w+n-j) - M_n-1(w+n-j-1)
! i.e.   new(j) = old(j-1) - old(j)
! where old is array before one_pass (thus n->n-1) and new is array afterwards
! do backwards to do in place

  integer j
  c(n) = c(n-1)
  do j = n-1,2,-1
    c(j) = c(j-1) - c(j)
  enddo
  c(1) = -c(1)
  return
end subroutine FH_PME_RECIP_bspline_diff
!----------------------------------------------------------------------
subroutine FH_PME_RECIP_grid_Fcoeffs(nsites, &
                                   num_prims_for_grid, &
                                   prim_offset_for_grid, &
                                   Frac_hermite_coeff,&
                                   theta1,theta2,theta3, &
                                   site_that_owns_grid_prim, &
                                   init_grid_ind_for_site, &
                                   Bspline_offset_for_site, &
                                   Fcoeff_offset_for_grid_prim, &
                                   hermite_order_for_grid_prim, &
                                   Max_deriv_per_site, &
                                   nfft1,nfft2,nfft3,Spline_order, &
                                   nfftdim1,nfftdim2,nfftdim3,pme_grid)

  implicit none

  integer, intent(in) :: nsites,num_prims_for_grid,prim_offset_for_grid
  double precision,intent(in) :: Frac_hermite_coeff(*)
  double precision,intent(in) :: theta1(*),theta2(*),theta3(*)
  integer,intent(in) :: site_that_owns_grid_prim(*)
  integer,intent(in) :: init_grid_ind_for_site(3,nsites), & 
                        Bspline_offset_for_site(nsites)
  integer,intent(in) :: Fcoeff_offset_for_grid_prim(*), &
                        hermite_order_for_grid_prim(*)
  integer,intent(in) :: Max_deriv_per_site(nsites)
  integer,intent(in) :: nfft1,nfft2,nfft3,Spline_order, &
                        nfftdim1,nfftdim2,nfftdim3
  double precision,intent(out) :: pme_grid(2*nfftdim1,nfftdim2,nfftdim3)
  include "mpole_index.fh"

  integer kp,np,n,ntot,ithoff,hc_off,dr_order,hc_order
  integer i0,i,j0,j,k0,k,ind1,ind2,ind3
  integer igrd0,jgrd0,kgrd0
  integer ith1,ith2,ith3
  double precision term0,term1,term2,term3,term4
  double precision t0,t1,t2,t3,t4
  double precision u0,u1,u2,u3,u4
  double precision v0,v1,v2,v3,v4

  ntot = 2*nfftdim1*nfftdim2*nfftdim3
  call UTIL_zero_real_array(pme_grid,ntot)

  do kp = 1,num_prims_for_grid
    np = prim_offset_for_grid + kp
    n = site_that_owns_grid_prim(np)
    ithoff = Bspline_offset_for_site(n)
    dr_order = Max_deriv_per_site(n) + 1 
            ! note add one for deriv orders 0,1,2,..,Max_deriv_per_site(n)
    hc_off = Fcoeff_offset_for_grid_prim(np)
    hc_order = hermite_order_for_grid_prim(np) !multipole order of prim np
        ! note add one for deriv orders 0,1,2,..,deriv_order(np) 
    igrd0 = init_grid_ind_for_site(1,n) !begin index in 1st direction
    jgrd0 = init_grid_ind_for_site(2,n) !begin index in 2nd direction
    kgrd0 = init_grid_ind_for_site(3,n) !begin index in 3rd direction
    k0 = kgrd0
    if ( hc_order == 0 )then
    ! do nothing. this prim has no hermite coefficients
    elseif ( hc_order == 1 )then
      do ith3 = 1,Spline_order
        k0 = k0 + 1
        k = k0 + 1 + (nfft3 - isign(nfft3,k0))/2
        j0 = jgrd0
        ind3 = ithoff + (ith3-1)*dr_order + 1
        v0 = theta3(ind3)
        do ith2 = 1,Spline_order
          j0 = j0 + 1
          j = j0 + 1 + (nfft2 - isign(nfft2,j0))/2
          i0 = igrd0
          ind2 = ithoff + (ith2-1)*dr_order + 1
          u0 = theta2(ind2) 
          term0 = Frac_hermite_coeff(hc_off+Ind_000)*u0*v0
          do ith1 = 1,Spline_order
            i0 = i0 + 1
            i = i0 + 1 + (nfft1 - isign(nfft1,i0))/2
            ind1 = ithoff + (ith1-1)*dr_order + 1
            t0 = theta1(ind1)
            pme_grid(i,j,k) = pme_grid(i,j,k) + term0*t0
          enddo
        enddo !ith2 = 1,Spline_order
      enddo !ith3 = 1,Spline_order
    else if ( hc_order == 4 )then
      do ith3 = 1,Spline_order
        k0 = k0 + 1
        k = k0 + 1 + (nfft3 - isign(nfft3,k0))/2
        ind3 = ithoff + (ith3-1)*dr_order + 1
        v0 = theta3(ind3)  !theta3
        v1 = theta3(ind3+1) !1st deriv of theta3
        j0 = jgrd0
        do ith2 = 1,Spline_order
          j0 = j0 + 1
          j = j0 + 1 + (nfft2 - isign(nfft2,j0))/2
          ind2 = ithoff + (ith2-1)*dr_order + 1
          u0 = theta2(ind2)  !theta2
          u1 = theta2(ind2+1) !1st deriv of theta2
! hardwire our knowledge of layout of theta1,2,3 to pre-assemble factors
          term0 = Frac_hermite_coeff(hc_off+Ind_000)*u0*v0 + &
                  Frac_hermite_coeff(hc_off+Ind_010)*u1*v0 + &
                  Frac_hermite_coeff(hc_off+Ind_001)*u0*v1
          term1 = Frac_hermite_coeff(hc_off+Ind_100)*u0*v0
          i0 = igrd0
          do ith1 = 1,Spline_order
            i0 = i0 + 1
            i = i0 + 1 + (nfft1 - isign(nfft1,i0))/2
            ind1 = ithoff + (ith1-1)*dr_order + 1
            t0 = theta1(ind1) !theta1
            t1 = theta1(ind1+1) !1st deriv of theta1
            pme_grid(i,j,k) = pme_grid(i,j,k) + term0*t0 + term1*t1
          enddo
        enddo !ith2 = 1,Spline_order
      enddo !ith3 = 1,Spline_order
    else if ( hc_order == 10 )then
      do ith3 = 1,Spline_order
        k0 = k0 + 1
        k = k0 + 1 + (nfft3 - isign(nfft3,k0))/2
        j0 = jgrd0
        ind3 = ithoff + (ith3-1)*dr_order + 1
        v0 = theta3(ind3)  !theta3
        v1 = theta3(ind3+1) !1st deriv of theta3
        v2 = theta3(ind3+2) !2nd deriv of theta3
        do ith2 = 1,Spline_order
          j0 = j0 + 1
          j = j0 + 1 + (nfft2 - isign(nfft2,j0))/2
          ind2 = ithoff + (ith2-1)*dr_order + 1
          u0 = theta2(ind2)  !theta2
          u1 = theta2(ind2+1) !1st deriv of theta2
          u2 = theta2(ind2+2) !2nd deriv of theta2
          i0 = igrd0
! hardwire our knowledge of layout of theta1,2,3 to pre-assemble factors
          term0 = Frac_hermite_coeff(hc_off+Ind_000)*u0*v0 + &
                  Frac_hermite_coeff(hc_off+Ind_010)*u1*v0 + &
                  Frac_hermite_coeff(hc_off+Ind_001)*u0*v1 + &
                  Frac_hermite_coeff(hc_off+Ind_020)*u2*v0 + &
                  Frac_hermite_coeff(hc_off+Ind_002)*u0*v2 + &
                  Frac_hermite_coeff(hc_off+Ind_011)*u1*v1
          term1 = Frac_hermite_coeff(hc_off+Ind_100)*u0*v0 + &
                  Frac_hermite_coeff(hc_off+Ind_110)*u1*v0 + &
                  Frac_hermite_coeff(hc_off+Ind_101)*u0*v1
          term2 = Frac_hermite_coeff(hc_off+Ind_200)*u0*v0
          do ith1 = 1,Spline_order
            i0 = i0 + 1
            i = i0 + 1 + (nfft1 - isign(nfft1,i0))/2
            ind1 = ithoff + (ith1-1)*dr_order + 1
            t0 = theta1(ind1) !theta1
            t1 = theta1(ind1+1) !1st deriv of theta1
            t2 = theta1(ind1+2) !2nd deriv of theta1
            pme_grid(i,j,k) = pme_grid(i,j,k) + term0*t0 + term1*t1 + &
                                                    term2*t2
          enddo
        enddo !ith2 = 1,Spline_order
      enddo !ith3 = 1,Spline_order
    elseif ( hc_order == 20 )then
      do ith3 = 1,Spline_order
        k0 = k0 + 1
        k = k0 + 1 + (nfft3 - isign(nfft3,k0))/2
        ind3 = ithoff + (ith3-1)*dr_order + 1
        v0 = theta3(ind3)  !theta3
        v1 = theta3(ind3+1) !1st deriv of theta3
        v2 = theta3(ind3+2) !2nd deriv of theta3
        v3 = theta3(ind3+3) !3rd deriv of theta3
        j0 = jgrd0
        do ith2 = 1,Spline_order
          j0 = j0 + 1
          j = j0 + 1 + (nfft2 - isign(nfft2,j0))/2
          ind2 = ithoff + (ith2-1)*dr_order + 1
          u0 = theta2(ind2)  !theta2
          u1 = theta2(ind2+1) !1st deriv of theta2
          u2 = theta2(ind2+2) !2nd deriv of theta2
          u3 = theta2(ind2+3) !3rd deriv of theta2
! hardwire our knowledge of layout of theta1,2,3 to pre-assemble factors
          term0 = Frac_hermite_coeff(hc_off+Ind_000)*u0*v0 + &
                  Frac_hermite_coeff(hc_off+Ind_010)*u1*v0 + &
                  Frac_hermite_coeff(hc_off+Ind_001)*u0*v1 + &
                  Frac_hermite_coeff(hc_off+Ind_020)*u2*v0 + &
                  Frac_hermite_coeff(hc_off+Ind_002)*u0*v2 + &
                  Frac_hermite_coeff(hc_off+Ind_011)*u1*v1 + &
                  Frac_hermite_coeff(hc_off+Ind_030)*u3*v0 + &
                  Frac_hermite_coeff(hc_off+Ind_003)*u0*v3 + &
                  Frac_hermite_coeff(hc_off+Ind_021)*u2*v1 + &
                  Frac_hermite_coeff(hc_off+Ind_012)*u1*v2
          term1 = Frac_hermite_coeff(hc_off+Ind_100)*u0*v0 + &
                  Frac_hermite_coeff(hc_off+Ind_110)*u1*v0 + &
                  Frac_hermite_coeff(hc_off+Ind_101)*u0*v1 + &
                  Frac_hermite_coeff(hc_off+Ind_120)*u2*v0 + &
                  Frac_hermite_coeff(hc_off+Ind_102)*u0*v2 + &
                  Frac_hermite_coeff(hc_off+Ind_111)*u1*v1
          term2 = Frac_hermite_coeff(hc_off+Ind_200)*u0*v0 + &
                  Frac_hermite_coeff(hc_off+Ind_210)*u1*v0 + &
                  Frac_hermite_coeff(hc_off+Ind_201)*u0*v1
          term3 = Frac_hermite_coeff(hc_off+Ind_300)*u0*v0
          i0 = igrd0
          do ith1 = 1,Spline_order
            i0 = i0 + 1
            i = i0 + 1 + (nfft1 - isign(nfft1,i0))/2
            ind1 = ithoff + (ith1-1)*dr_order + 1
            t0 = theta1(ind1) !theta1
            t1 = theta1(ind1+1) !1st deriv of theta1
            t2 = theta1(ind1+2) !2nd deriv of theta1
            t3 = theta1(ind1+3) !3rd deriv of theta1
            pme_grid(i,j,k) = pme_grid(i,j,k) + term0*t0 + term1*t1 &
                                + term2*t2 + term3*t3
          enddo
        enddo !ith2 = 1,Spline_order
      enddo !ith3 = 1,Spline_order
    else if ( hc_order == 35 )then
      do ith3 = 1,Spline_order
        k0 = k0 + 1
        k = k0 + 1 + (nfft3 - isign(nfft3,k0))/2
        j0 = jgrd0
        ind3 = ithoff + (ith3-1)*dr_order + 1
        v0 = theta3(ind3)  !theta3
        v1 = theta3(ind3+1) !1st deriv of theta3
        v2 = theta3(ind3+2) !2nd deriv of theta3
        v3 = theta3(ind3+3) !3rd deriv of theta3
        v4 = theta3(ind3+4) !4th deriv of theta3
        do ith2 = 1,Spline_order
          j0 = j0 + 1
          j = j0 + 1 + (nfft2 - isign(nfft2,j0))/2
          ind2 = ithoff + (ith2-1)*dr_order + 1
          u0 = theta2(ind2)  !theta2
          u1 = theta2(ind2+1) !1st deriv of theta2
          u2 = theta2(ind2+2) !2nd deriv of theta2
          u3 = theta2(ind2+3) !3rd deriv of theta2
          u4 = theta2(ind2+4) !4th deriv of theta2
          i0 = igrd0
! hardwire our knowledge of layout of theta1,2,3 to pre-assemble factors
          term0 = Frac_hermite_coeff(hc_off+Ind_000)*u0*v0 + &
                  Frac_hermite_coeff(hc_off+Ind_010)*u1*v0 + &
                  Frac_hermite_coeff(hc_off+Ind_001)*u0*v1 + &
                  Frac_hermite_coeff(hc_off+Ind_020)*u2*v0 + &
                  Frac_hermite_coeff(hc_off+Ind_002)*u0*v2 + &
                  Frac_hermite_coeff(hc_off+Ind_011)*u1*v1 + &
                  Frac_hermite_coeff(hc_off+Ind_030)*u3*v0 + &
                  Frac_hermite_coeff(hc_off+Ind_003)*u0*v3 + &
                  Frac_hermite_coeff(hc_off+Ind_021)*u2*v1 + &
                  Frac_hermite_coeff(hc_off+Ind_012)*u1*v2 + &
                  Frac_hermite_coeff(hc_off+Ind_040)*u4*v0 + &
                  Frac_hermite_coeff(hc_off+Ind_004)*u0*v4 + &
                  Frac_hermite_coeff(hc_off+Ind_031)*u3*v1 + &
                  Frac_hermite_coeff(hc_off+Ind_013)*u1*v3 + &
                  Frac_hermite_coeff(hc_off+Ind_022)*u2*v2
          term1 = Frac_hermite_coeff(hc_off+Ind_100)*u0*v0 + &
                  Frac_hermite_coeff(hc_off+Ind_110)*u1*v0 + &
                  Frac_hermite_coeff(hc_off+Ind_101)*u0*v1 + &
                  Frac_hermite_coeff(hc_off+Ind_120)*u2*v0 + &
                  Frac_hermite_coeff(hc_off+Ind_102)*u0*v2 + &
                  Frac_hermite_coeff(hc_off+Ind_111)*u1*v1 + &
                  Frac_hermite_coeff(hc_off+Ind_130)*u3*v0 + &
                  Frac_hermite_coeff(hc_off+Ind_103)*u0*v3 + &
                  Frac_hermite_coeff(hc_off+Ind_121)*u2*v1 + &
                  Frac_hermite_coeff(hc_off+Ind_112)*u1*v2
          term2 = Frac_hermite_coeff(hc_off+Ind_200)*u0*v0 + &
                  Frac_hermite_coeff(hc_off+Ind_210)*u1*v0 + &
                  Frac_hermite_coeff(hc_off+Ind_201)*u0*v1 + &
                  Frac_hermite_coeff(hc_off+Ind_220)*u2*v0 + &
                  Frac_hermite_coeff(hc_off+Ind_202)*u0*v2 + &
                  Frac_hermite_coeff(hc_off+Ind_211)*u1*v1
          term3 = Frac_hermite_coeff(hc_off+Ind_300)*u0*v0 + &
                  Frac_hermite_coeff(hc_off+Ind_310)*u1*v0 + &
                  Frac_hermite_coeff(hc_off+Ind_301)*u0*v1
          term4 = Frac_hermite_coeff(hc_off+Ind_400)*u0*v0
          do ith1 = 1,Spline_order
            i0 = i0 + 1
            i = i0 + 1 + (nfft1 - isign(nfft1,i0))/2
            ind1 = ithoff + (ith1-1)*dr_order + 1
            t0 = theta1(ind1) !theta1
            t1 = theta1(ind1+1) !1st deriv of theta1
            t2 = theta1(ind1+2) !2nd deriv of theta1
            t3 = theta1(ind1+3) !3rd deriv of theta1
            t4 = theta1(ind1+4) !4th deriv of theta1
            pme_grid(i,j,k) = pme_grid(i,j,k) + term0*t0 + term1*t1 &
                                + term2*t2 + term3*t3 + term4*t4
          enddo
        enddo !ith2 = 1,Spline_order
      enddo !ith3 = 1,Spline_order
    endif   ! if branch over hc_order
  enddo !np = 1,num_prims_to_grid

end subroutine FH_PME_RECIP_grid_Fcoeffs
!-----------------------------------------------------------
subroutine FH_PME_RECIP_grid_mult_orthog( &
                           nfft1,nfft2,nfft3,nfftdim1, &
                           grid_expon,weight_grid_expon,num_exp, &
                           recip,volume,prefac1,prefac2,prefac3, &
                           pme_multiplier)

   implicit none

   integer, intent(in) :: nfft1,nfft2,nfft3,nfftdim1
   double precision, intent(in) :: grid_expon(*),weight_grid_expon(*)
   integer, intent(in) :: num_exp
   double precision, intent(in) :: recip(3,3),volume, &
                                   prefac1(2,*), prefac2(2,*),prefac3(2,*)
   double precision, intent(out) :: pme_multiplier(2,nfft3,nfftdim1,nfft2)
   include 'mpole_exponent.fh'

   integer :: k1,k2,k3,m1,m2,m3,k10,nf1,nf2,nf3,itot,kexp,siz
   double precision :: pi,fac,term,mhat1,mhat2,mhat3,gfunc23, &
                       gfunc123,tmp23r,tmp23i,tmp123r,tmp123i,vol_inv

   double precision :: m1_tbl(-(nfft1/2 + 1) : nfft1/2 + 1)
   double precision :: m2_tbl(-(nfft2/2 + 1) : nfft2/2 + 1)
   double precision :: m3_tbl(-(nfft3/2 + 1) : nfft3/2 + 1)

   pi = 3.14159265358979323846d0
   vol_inv = 1.d0 / volume
   siz = 2*nfft3*nfftdim1*nfft2
   call UTIL_zero_real_array(pme_multiplier,siz)

   nf1 = nfft1/2
   if ( 2*nf1 < nfft1 )nf1 = nf1+1
   nf2 = nfft2/2
   if ( 2*nf2 < nfft2 )nf2 = nf2+1
   nf3 = nfft3/2
   if ( 2*nf3 < nfft3 )nf3 = nf3+1

   do kexp = 1,num_exp
      term = weight_grid_expon(kexp)*vol_inv
      fac = pi*pi / grid_expon(kexp)
      if ( grid_expon(kexp) < MPOLE_exponent )then
         itot = 0
         do m1 = -(nfft1/2 + 1), nfft1/2 + 1
            itot = itot + 1
            mhat1 = recip(1,1)*m1
            m1_tbl(m1) = -fac*mhat1*mhat1
         enddo
         call vdexp(itot,m1_tbl(-nfft1/2-1),m1_tbl(-nfft1/2-1))

         itot = 0
         do m2 = -(nfft2/2 + 1), nfft2/2 + 1
            itot = itot + 1
            mhat2 = recip(2,2)*m2
            m2_tbl(m2) = -fac*mhat2*mhat2
         enddo
         call vdexp(itot,m2_tbl(-nfft2/2-1),m2_tbl(-nfft2/2-1))

         itot = 0
         fac = pi*pi / grid_expon(kexp)
         do m3 = -(nfft3/2 + 1), nfft3/2 + 1
            itot = itot + 1
            mhat3 = recip(3,3)*m3
            m3_tbl(m3) = -fac*mhat3*mhat3
         enddo
         call vdexp(itot,m3_tbl(-nfft3/2-1),m3_tbl(-nfft3/2-1))
      endif !( grid_expon(kexp) < MPOLE_exponent )

      do k2 = 1, nfft2
         m2 = k2 - 1
         if ( k2 > nf2 )m2 = k2 - 1 - nfft2
         do k3 = 1,nfft3
            m3 = k3 - 1
            if ( k3 > nf3 )m3 = k3 - 1 - nfft3
            k10 = 1
            ! need (1,1,1) case also
            !if(k3+k2 == 2) k10 = 2
         
            gfunc23 = term*m2_tbl(m2)*m3_tbl(m3)
            tmp23r = prefac2(1,k2)*prefac3(1,k3) - prefac2(2,k2)*prefac3(2,k3)
            tmp23i = prefac2(1,k2)*prefac3(2,k3) + prefac2(2,k2)*prefac3(1,k3)
         
            do k1 = k10, nf1+1
               m1 = k1 - 1
               if ( k1 > nf1 )m1 = k1 - 1 - nfft1
               if ( grid_expon(kexp) < MPOLE_exponent )then
                  gfunc123 = gfunc23*m1_tbl(m1)
               else
                  gfunc123 = term
               endif
               tmp123r = prefac1(1,k1)*tmp23r - prefac1(2,k1)*tmp23i
               tmp123i = prefac1(1,k1)*tmp23i + prefac1(2,k1)*tmp23r
               pme_multiplier(1,k3,k1,k2) =  &
               pme_multiplier(1,k3,k1,k2) + gfunc123*tmp123r
               pme_multiplier(2,k3,k1,k2) =  &
               pme_multiplier(2,k3,k1,k2) + gfunc123*tmp123i
            enddo ! k1 = k10, nf1+1
         enddo !k3 = 1,nfft3
      enddo !k2 = 1, nfft2
   enddo!kexp = 1,num_exp
end subroutine FH_PME_RECIP_grid_mult_orthog
!-----------------------------------------------------------
subroutine FH_PME_RECIP_fix_one_FT_density( &
                           nfft1,nfft2,nfft3,nfftdim1, &
                           pme_multiplier, &
                           pme_grid)

   implicit none

   integer, intent(in) :: nfft1,nfft2,nfft3,nfftdim1
   double precision, intent(in) :: pme_multiplier(2,nfft3,nfftdim1,nfft2)
   double precision, intent(inout) :: pme_grid(2,nfft3,nfftdim1,nfft2)

   integer :: k1,k2,k3,k10,nf1
   double precision :: tmpr,tmpi

   nf1 = nfft1/2
   if ( 2*nf1 < nfft1 )nf1 = nf1+1
   do k2 = 1, nfft2
      do k3 = 1,nfft3
         k10 = 1
         ! need (1,1,1) case also
         !if ( k3+k2 == 2 )k10 = 2
         do k1 = k10, nf1+1
            tmpr = pme_multiplier(1,k3,k1,k2)*pme_grid(1,k3,k1,k2) - &
                   pme_multiplier(2,k3,k1,k2)*pme_grid(2,k3,k1,k2)
            tmpi = pme_multiplier(1,k3,k1,k2)*pme_grid(2,k3,k1,k2) + &
                   pme_multiplier(2,k3,k1,k2)*pme_grid(1,k3,k1,k2)
            pme_grid(1,k3,k1,k2) = tmpr
            pme_grid(2,k3,k1,k2) = tmpi
         enddo
      enddo
   enddo
end subroutine FH_PME_RECIP_fix_one_FT_density
!-----------------------------------------------------------
subroutine FH_PME_RECIP_fix_one_FT_phi( &
                           nfft1_1,nfft2_1,nfft3_1,nfftdim1_1, &
                           volume, &
                           pme_grid_multiplier,pme_grid)

   implicit none

   integer, intent(in) :: nfft1_1,nfft2_1,nfft3_1,nfftdim1_1
   double precision, intent(in) :: volume, & 
                           pme_grid_multiplier(2,nfft3_1,nfftdim1_1,nfft2_1)
   double precision, intent(out) :: pme_grid(2,nfft3_1,nfftdim1_1,nfft2_1) 

   integer :: k1,k2,k3,k10,nf1
   double precision :: tmpr,tmpi

   nf1 = nfft1_1/2
   if ( 2*nf1 < nfft1_1 )nf1 = nf1+1
   do k2 = 1, nfft2_1
      do k3 = 1,nfft3_1
         k10 = 1
         ! need (1,1,1) case also
         !if ( k3+k2 == 2 )k10 = 2
         do k1 = k10, nf1+1
            ! use complex conjugate of pme_grid_multiplier
            tmpr = pme_grid_multiplier(1,k3,k1,k2)*pme_grid(1,k3,k1,k2) + &
                   pme_grid_multiplier(2,k3,k1,k2)*pme_grid(2,k3,k1,k2)
            tmpi = pme_grid_multiplier(1,k3,k1,k2)*pme_grid(2,k3,k1,k2) - &
                   pme_grid_multiplier(2,k3,k1,k2)*pme_grid(1,k3,k1,k2)
            pme_grid(1,k3,k1,k2) = tmpr
            pme_grid(2,k3,k1,k2) = tmpi
            !pme_grid(1,k3,k1,k2) = tmpr*volume !mult by volume for plancherel
            !pme_grid(2,k3,k1,k2) = tmpi*volume
         enddo
      enddo
   enddo

end subroutine FH_PME_RECIP_fix_one_FT_phi
!-----------------------------------------------------------
subroutine FH_PME_RECIP_get_one_FHerm_phi( &
                                   nsites, &
                                   num_prims_for_grid, &
                                   prim_offset_for_grid, &
                                   theta1,theta2,theta3, &
                                   site_that_owns_grid_prim, &
                                   init_grid_ind_for_site, &
                                   Bspline_offset_for_site, &
                                   Field_offset_for_grid_prim, &
                                   Field_order_for_grid_prim, &
                                   Max_deriv_per_site, &
                                   nfft1,nfft2,nfft3,Spline_order, &
                                   nfftdim1,nfftdim2,nfftdim3,pme_grid, &
                                   Fphi, out_lun)

  implicit none

  integer, intent(in) :: nsites,num_prims_for_grid,prim_offset_for_grid
  double precision,intent(in) :: theta1(*),theta2(*),theta3(*)
  integer,intent(in) :: site_that_owns_grid_prim(*)
  integer,intent(in) :: init_grid_ind_for_site(3,nsites), & 
                        Bspline_offset_for_site(nsites)
  integer,intent(in) :: Field_offset_for_grid_prim(*), &
                        Field_order_for_grid_prim(*)
  integer,intent(in) :: Max_deriv_per_site(nsites)
  integer,intent(in) :: nfft1,nfft2,nfft3,Spline_order, &
                        nfftdim1,nfftdim2,nfftdim3
  double precision,intent(in) :: pme_grid(2*nfftdim1,nfftdim2,nfftdim3)
  double precision, intent(out) :: Fphi(*)
  integer, intent(in)   :: out_lun

  include "mpole_index.fh"
  integer field_off,dr_order,field_ord
  integer ithoff
  integer kp,np,n,ith1,ith2,ith3,i0,j0,k0,i,j,k
  integer ind1,ind2,ind3
  integer igrd0,jgrd0,kgrd0
 
  double precision t0,t1,t2,t3,t4,t5
  double precision u0,u1,u2,u3,u4,u5
  double precision v0,v1,v2,v3,v4,v5
  double precision tu00,tu01,tu10,tu20,tu02,tu11,tu30, &
         tu03,tu21,tu12,tu40,tu04,tu31,tu13,tu22,tu50,tu05, &
         tu41,tu14,tu32,tu23
  double precision tuv000,tuv100,tuv010,tuv001
  double precision tuv200,tuv020,tuv002,tuv110,tuv101,tuv011
  double precision tuv300,tuv030,tuv003,tuv210,tuv201,tuv120, &
         tuv021,tuv102,tuv012,tuv111
  double precision tuv400,tuv040,tuv004,tuv310,tuv301, &
         tuv130,tuv031,tuv103,tuv013,tuv220,tuv202,tuv022, &
         tuv211, tuv121,tuv112
  double precision tuv500,tuv050,tuv005,tuv410,tuv401, &
         tuv140,tuv041,tuv104,tuv014,tuv320,tuv302,tuv230, &
         tuv032,tuv203,tuv023,tuv311,tuv131,tuv113,tuv221, &
         tuv212,tuv122
  double precision tq

  do kp = 1,num_prims_for_grid
    np = prim_offset_for_grid + kp
    field_ord = Field_order_for_grid_prim(np) !field order of prim np
    if ( field_ord == 0 )then
        ! do nothing. this site has no coefficients for this hermite set
    elseif ( field_ord > 56 )then ! later maybe fill this in with general
        write(out_lun,*)'field order > 56!!!', field_ord
        stop 
    else !field_ord > 0
      field_off = Field_offset_for_grid_prim(np)
      n = site_that_owns_grid_prim(np)
      ithoff = Bspline_offset_for_site(n)
      dr_order = Max_deriv_per_site(n) + 1 
            ! note add one for deriv orders 0,1,2,..,Max_deriv_per_site(n)
      igrd0 = init_grid_ind_for_site(1,n) !begin index in 1st direction
      jgrd0 = init_grid_ind_for_site(2,n) !begin index in 2nd direction
      kgrd0 = init_grid_ind_for_site(3,n) !begin index in 3rd direction

      tuv000 = 0.d0
      if ( field_ord > 1 )then
        tuv001 = 0.d0
        tuv010 = 0.d0
        tuv100 = 0.d0
        if ( field_ord > 4 )then
          tuv200 = 0.d0
          tuv020 = 0.d0
          tuv002 = 0.d0
          tuv110 = 0.d0
          tuv101 = 0.d0
          tuv011 = 0.d0
          if ( field_ord > 10 )then
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
            if ( field_ord > 20 )then
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
              if ( field_ord == 56 )then
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
      k0 = kgrd0
      do ith3 = 1,Spline_order
        ind3 = ithoff + (ith3 - 1)*dr_order + 1
        v0 = theta3(ind3) !theta3
        tu00 = 0.d0
        if ( field_ord > 1 )then
          v1 = theta3(ind3 + 1) !1st deriv of theta3
          tu10 = 0.d0
          tu01 = 0.d0
          if ( field_ord > 4 )then
            v2 = theta3(ind3 + 2) !2nd deriv of theta3
            tu20 = 0.d0
            tu11 = 0.d0
            tu02 = 0.d0
            if( field_ord > 10 )then
              v3 = theta3(ind3 + 3) !3rd deriv of theta3
              tu30 = 0.d0
              tu21 = 0.d0
              tu12 = 0.d0
              tu03 = 0.d0
              if( field_ord > 20 )then
                v4 = theta3(ind3 + 4) !4th deriv of theta3
                tu40 = 0.d0
                tu31 = 0.d0
                tu22 = 0.d0
                tu13 = 0.d0
                tu04 = 0.d0
                if ( field_ord == 56 )then
                  v5 = theta3(ind3 + 5) !5th deriv of theta3
                  tu50 = 0.d0
                  tu41 = 0.d0
                  tu32 = 0.d0
                  tu23 = 0.d0
                  tu14 = 0.d0
                  tu05 = 0.d0
                endif !( field_ord == 56 )
              endif !( field_ord > 20 )
            endif !( field_ord > 10 )
          endif !( field_ord > 4 )
        endif !( field_ord > 1 )
        k0 = k0 + 1
        k = k0 + 1 + (nfft3 - isign(nfft3,k0))/2
        j0 = jgrd0
        do ith2 = 1,Spline_order
          j0 = j0 + 1
          j = j0 + 1 + (nfft2 - isign(nfft2,j0))/2
          i0 = igrd0
          ind2 = ithoff + (ith2 - 1)*dr_order + 1
          if ( field_ord == 1 )then
            u0 = theta2(ind2) !theta2
            t0 = 0.d0
            do ith1 = 1,Spline_order
              i0 = i0 + 1
              i = i0 + 1 + (nfft1 - isign(nfft1,i0))/2
              ind1 = ithoff + (ith1-1)*dr_order + 1
              tq = pme_grid(i,j,k)
              t0 = t0 + tq*theta1(ind1)
            enddo !ith1 = 1,Spline_order
            tu00 = tu00 + t0*u0
          elseif ( field_ord == 4 )then
            u0 = theta2(ind2) !theta2
            u1 = theta2(ind2 + 1) !1st deriv of theta2
            t0 = 0.d0
            t1 = 0.d0
            do ith1 = 1,Spline_order
              i0 = i0 + 1
              i = i0 + 1 + (nfft1 - isign(nfft1,i0))/2
              ind1 = ithoff + (ith1-1)*dr_order + 1
              tq = pme_grid(i,j,k)
              t0 = t0 + tq*theta1(ind1)
              t1 = t1 + tq*theta1(ind1+1)
            enddo !ith1 = 1,Spline_order
            tu00 = tu00 + t0*u0

            tu10 = tu10 + t1*u0
            tu01 = tu01 + t0*u1
          elseif ( field_ord == 10 )then
            u0 = theta2(ind2) !theta2
            u1 = theta2(ind2 + 1) !1st deriv of theta2
            u2 = theta2(ind2 + 2) !2nd deriv of theta2
            t0 = 0.d0
            t1 = 0.d0
            t2 = 0.d0
            do ith1 = 1,Spline_order
              i0 = i0 + 1
              i = i0 + 1 + (nfft1 - isign(nfft1,i0))/2
              ind1 = ithoff + (ith1-1)*dr_order + 1
              tq = pme_grid(i,j,k)
              t0 = t0 + tq*theta1(ind1)
              t1 = t1 + tq*theta1(ind1+1)
              t2 = t2 + tq*theta1(ind1+2)
            enddo !ith1 = 1,Spline_order
            tu00 = tu00 + t0*u0

            tu10 = tu10 + t1*u0
            tu01 = tu01 + t0*u1

            tu20 = tu20 + t2*u0
            tu11 = tu11 + t1*u1
            tu02 = tu02 + t0*u2
          elseif ( field_ord == 20 )then
            u0 = theta2(ind2) !theta2
            u1 = theta2(ind2 + 1) !1st deriv of theta2
            u2 = theta2(ind2 + 2) !2nd deriv of theta2
            u3 = theta2(ind2 + 3) !3rd deriv of theta2
            t0 = 0.d0
            t1 = 0.d0
            t2 = 0.d0
            t3 = 0.d0
            do ith1 = 1,Spline_order
              i0 = i0 + 1
              i = i0 + 1 + (nfft1 - isign(nfft1,i0))/2
              ind1 = ithoff + (ith1-1)*dr_order + 1
              tq = pme_grid(i,j,k)
              t0 = t0 + tq*theta1(ind1)
              t1 = t1 + tq*theta1(ind1+1)
              t2 = t2 + tq*theta1(ind1+2)
              t3 = t3 + tq*theta1(ind1+3)
            enddo !ith1 = 1,Spline_order
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
          elseif ( field_ord == 35 )then
            u0 = theta2(ind2) !theta2
            u1 = theta2(ind2 + 1) !1st deriv of theta2
            u2 = theta2(ind2 + 2) !2nd deriv of theta2
            u3 = theta2(ind2 + 3) !3rd deriv of theta2
            u4 = theta2(ind2 + 4) !4th deriv of theta2
            t0 = 0.d0
            t1 = 0.d0
            t2 = 0.d0
            t3 = 0.d0
            t4 = 0.d0
            do ith1 = 1,Spline_order
              i0 = i0 + 1
              i = i0 + 1 + (nfft1 - isign(nfft1,i0))/2
              ind1 = ithoff + (ith1-1)*dr_order + 1
              tq = pme_grid(i,j,k)
              t0 = t1 + tq*theta1(ind1)
              t1 = t1 + tq*theta1(ind1+1)
              t2 = t2 + tq*theta1(ind1+2)
              t3 = t3 + tq*theta1(ind1+3)
              t4 = t4 + tq*theta1(ind1+4)
            enddo !ith1 = 1,Spline_order
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
          elseif ( field_ord == 56 )then
            u0 = theta2(ind2) !theta2
            u1 = theta2(ind2 + 1) !1st deriv of theta2
            u2 = theta2(ind2 + 2) !2nd deriv of theta2
            u3 = theta2(ind2 + 3) !3rd deriv of theta2
            u4 = theta2(ind2 + 4) !4th deriv of theta2
            u5 = theta2(ind2 + 5) !5th deriv of theta2
            t0 = 0.d0
            t1 = 0.d0
            t2 = 0.d0
            t3 = 0.d0
            t4 = 0.d0
            t5 = 0.d0
            do ith1 = 1,Spline_order
              i0 = i0 + 1
              i = i0 + 1 + (nfft1 - isign(nfft1,i0))/2
              ind1 = ithoff + (ith1-1)*dr_order + 1
              tq = pme_grid(i,j,k)
              t0 = t1 + tq*theta1(ind1)
              t1 = t1 + tq*theta1(ind1+1)
              t2 = t2 + tq*theta1(ind1+2)
              t3 = t3 + tq*theta1(ind1+3)
              t4 = t4 + tq*theta1(ind1+4)
              t5 = t5 + tq*theta1(ind1+5)
            enddo !ith1 = 1,Spline_order
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
        enddo !ith2 = 1,Spline_order
        tuv000 = tuv000 + tu00*v0
        if ( field_ord > 1 )then 
          tuv100 = tuv100 + tu10*v0
          tuv010 = tuv010 + tu01*v0
          tuv001 = tuv001 + tu00*v1
          if ( field_ord > 4 )then
            tuv200 = tuv200 + tu20*v0
            tuv020 = tuv020 + tu02*v0
            tuv002 = tuv002 + tu00*v2
            tuv110 = tuv110 + tu11*v0
            tuv101 = tuv101 + tu10*v1
            tuv011 = tuv011 + tu01*v1
            if ( field_ord > 10 )then
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
              if ( field_ord > 20 )then
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
                if ( field_ord == 56 )then
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
                endif !( field_ord == 56 )
              endif !( field_ord > 20 )
            endif !( field_ord > 10 )
          endif !( field_ord > 4 )
        endif !( field_ord > 1 )
      enddo !ith3 = 1,Spline_order
      ! now fill fields; field_off starts this atoms fields
      Fphi(field_off+Ind_000) = tuv000
      if ( field_ord > 1 )then
        Fphi(field_off+Ind_100) = tuv100
        Fphi(field_off+Ind_010) = tuv010
        Fphi(field_off+Ind_001) = tuv001
        if ( field_ord > 4 )then
          Fphi(field_off+Ind_200) = tuv200
          Fphi(field_off+Ind_020) = tuv020
          Fphi(field_off+Ind_002) = tuv002
          Fphi(field_off+Ind_110) = tuv110
          Fphi(field_off+Ind_101) = tuv101
          Fphi(field_off+Ind_011) = tuv011
          if ( field_ord > 10 )then
            Fphi(field_off+Ind_300) = tuv300
            Fphi(field_off+Ind_030) = tuv030
            Fphi(field_off+Ind_003) = tuv003
            Fphi(field_off+Ind_210) = tuv210
            Fphi(field_off+Ind_201) = tuv201
            Fphi(field_off+Ind_120) = tuv120
            Fphi(field_off+Ind_021) = tuv021
            Fphi(field_off+Ind_102) = tuv102
            Fphi(field_off+Ind_012) = tuv012
            Fphi(field_off+Ind_111) = tuv111
            if ( field_ord > 20 )then
              Fphi(field_off+Ind_400) = tuv400
              Fphi(field_off+Ind_040) = tuv040
              Fphi(field_off+Ind_004) = tuv004
              Fphi(field_off+Ind_310) = tuv310
              Fphi(field_off+Ind_301) = tuv301
              Fphi(field_off+Ind_130) = tuv130
              Fphi(field_off+Ind_031) = tuv031
              Fphi(field_off+Ind_103) = tuv103
              Fphi(field_off+Ind_013) = tuv013
              Fphi(field_off+Ind_220) = tuv220
              Fphi(field_off+Ind_202) = tuv202
              Fphi(field_off+Ind_022) = tuv022
              Fphi(field_off+Ind_211) = tuv211
              Fphi(field_off+Ind_121) = tuv121
              Fphi(field_off+Ind_112) = tuv112
              if ( field_ord == 56 )then
                Fphi(field_off+Ind_500) = tuv500
                Fphi(field_off+Ind_050) = tuv050
                Fphi(field_off+Ind_005) = tuv005
                Fphi(field_off+Ind_410) = tuv410
                Fphi(field_off+Ind_401) = tuv401
                Fphi(field_off+Ind_140) = tuv140
                Fphi(field_off+Ind_041) = tuv041
                Fphi(field_off+Ind_104) = tuv104
                Fphi(field_off+Ind_014) = tuv014
                Fphi(field_off+Ind_320) = tuv320
                Fphi(field_off+Ind_302) = tuv302
                Fphi(field_off+Ind_230) = tuv230
                Fphi(field_off+Ind_032) = tuv032
                Fphi(field_off+Ind_203) = tuv203
                Fphi(field_off+Ind_023) = tuv023
                Fphi(field_off+Ind_311) = tuv311
                Fphi(field_off+Ind_131) = tuv131
                Fphi(field_off+Ind_113) = tuv113
                Fphi(field_off+Ind_221) = tuv221
                Fphi(field_off+Ind_212) = tuv212
                Fphi(field_off+Ind_122) = tuv122
              endif !( field_ord == 56 )
            endif !( field_ord > 20 )
          endif !( field_ord > 10 )
        endif !( field_ord > 4 )
      endif !( field_ord > 1 ) 
    endif !( field_ord == 0 ) outermost branch on field order
  enddo !kp = 1,num_prims_for_grid
  return
end subroutine FH_PME_RECIP_get_one_FHerm_phi
!------------------------------------------------------------
subroutine FH_PME_RECIP_one_ene_force(nfft1,nfft2,nfft3,recip, &
                   num_prims_for_grid,prim_offset_for_grid, &
                   site_that_owns_grid_prim, &
                   Fcoeff_offset_for_grid_prim, &
                   hermite_order_for_grid_prim, &
                   Ffield_offset_for_grid_prim, &
                   Frac_hermite_coeffs,Frac_hermite_field, &
                   factor,force,energy)

   implicit none

   integer, intent(in) :: nfft1,nfft2,nfft3
   double precision, intent(in) :: recip(3,3)
   integer, intent(in) :: num_prims_for_grid,prim_offset_for_grid
   integer, intent(in) :: site_that_owns_grid_prim(*), &
                          Fcoeff_offset_for_grid_prim(*), &
                          hermite_order_for_grid_prim(*), &
                          Ffield_offset_for_grid_prim(*)
   double precision, intent(in) :: Frac_hermite_coeffs(*)
   double precision, intent(in) :: Frac_hermite_field(*),factor
   double precision, intent(inout) :: energy,force(3,*)

   include "direct_pointers.fh"
   integer :: j1,j2,j3,j,kp,n,np,order,off1,off2
   double precision :: f1,f2,f3,dfx,dfy,dfz,ene
   double precision mpole_xform_3x3(3,3),field_xform_3x3(3,3)

   ene = 0.d0
!  first get field_xform_3x3
   call FH_PME_RECIP_jacobian(nfft1,nfft2,nfft3,recip, &
                             mpole_xform_3x3,field_xform_3x3)
   do kp = 1,num_prims_for_grid
      np = prim_offset_for_grid + kp
      n = site_that_owns_grid_prim(np)
      order = hermite_order_for_grid_prim(np)
      off1 = Fcoeff_offset_for_grid_prim(np)
      off2 = Ffield_offset_for_grid_prim(np)
      f1 = 0.d0
      f2 = 0.d0
      f3 = 0.d0
      do j = 1,order
         ene = ene + Frac_hermite_coeffs(off1+j)*Frac_hermite_field(off2+j)
         j1 = p_xder(j) 
         j2 = p_yder(j) 
         j3 = p_zder(j) 
         f1 = f1 + Frac_hermite_coeffs(off1+j)*Frac_hermite_field(off2+j1)
         f2 = f2 + Frac_hermite_coeffs(off1+j)*Frac_hermite_field(off2+j2)
         f3 = f3 + Frac_hermite_coeffs(off1+j)*Frac_hermite_field(off2+j3)
      enddo
      ! force is negative of gradient
      ! transform from scaled fractional to cartesian--same as fields
      dfx = field_xform_3x3(1,1)*f1+field_xform_3x3(1,2)*f2+ &
            field_xform_3x3(1,3)*f3
      dfy = field_xform_3x3(2,1)*f1+field_xform_3x3(2,2)*f2+ &
            field_xform_3x3(2,3)*f3
      dfz = field_xform_3x3(3,1)*f1+field_xform_3x3(3,2)*f2+ &
            field_xform_3x3(3,3)*f3
      force(1,n) = force(1,n) - factor*dfx
      force(2,n) = force(2,n) - factor*dfy
      force(3,n) = force(3,n) - factor*dfz
   enddo
   energy = energy + 0.5d0*factor*ene
end subroutine FH_PME_RECIP_one_ene_force
!------------------------------------------------------------
subroutine FH_PME_RECIP_one_torque( &
                               nfft1,nfft2,nfft3,recip, &
                               num_prims_for_grid,prim_offset_for_grid, &
                               hermite_order_for_grid_prim, &
                               Ffield_offset_for_grid_prim, &
                               Gcoeff_offset_for_grid_prim, &
                               factor, &
                               Frac_hermite_field, &
                               Global_hermite_field )

   implicit none

   integer,intent(in) :: num_prims_for_grid,prim_offset_for_grid, &
                         nfft1,nfft2,nfft3
   double precision, intent(in) :: recip(3,3)
   integer,intent(in) :: hermite_order_for_grid_prim(*), &
                         Ffield_offset_for_grid_prim(*), &
                         Gcoeff_offset_for_grid_prim(*)
   double precision, intent(in) :: factor,Frac_hermite_field(*)
   double precision, intent(inout) :: Global_hermite_field(*)

   include "mpole_sizes.fh"
   double precision :: mpole_xform_3x3(3,3),field_xform_3x3(3,3)
   double precision :: Field_xy(MAXMP*MAXMP),Cfield(MAXMP)
   integer :: dimxy,order,off,off_G,j,kp,np

   Cfield = 0.d0
!  first get field_xform_3x3
   call FH_PME_RECIP_jacobian(nfft1,nfft2,nfft3,recip, &
                             mpole_xform_3x3,field_xform_3x3)
  ! get the higher order terms
   call XFORM_MPOLE_field_matrix(field_xform_3x3,Field_xy,MAXMP)
   dimxy = MAXMP
   do kp = 1,num_prims_for_grid
      np = prim_offset_for_grid + kp
      order = hermite_order_for_grid_prim(np)
      off = Ffield_offset_for_grid_prim(np)
      off_G = Gcoeff_offset_for_grid_prim(np)
      call XFORM_FIELD(Field_xy,dimxy,Frac_hermite_field(off+1),Cfield,order)
      do j = 1,order
         Global_hermite_field(off_G+j) = Global_hermite_field(off_G+j) + &
                        factor*CField(j)
      enddo
   enddo
end subroutine FH_PME_RECIP_one_torque
!------------------------------------------------------------
