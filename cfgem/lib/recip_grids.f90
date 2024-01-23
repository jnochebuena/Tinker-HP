!---------------------------------------------------------------------
! There are 4 types of main grids handled here:
! 1) mpole grids and supporting data structures
! 2) sumCH grids and supporting data structures
! 3) compact hermite grids and supporting data structures
! 4) diffuse hermite grids and supporting data structures
! 
! In direct space the following interactions are calculated:
! multipole <> multipole 
!        (damped coulomb i.e. erfc(beta*rij) / rij in place of 1 / rij  
!         where beta is ewald coefficient)    for site1,site2 pairs
! multipole <> compact hermites (damped coulomb) for site1,site2 pairs
! compact hermites <> multipoles (damped hermites) for site1,site2 pairs
! compact hermites <> compact hermites (damped coulomb) for site1,site2 pairs
!
!  The rest of the coulomb for these interactions are calculated in reciprocal
!  space:
 !  e.g. multipole <> multipole with erf(beta*rij) / rij in place of 1/rij
!  where beta is ewald coefficient
!  To calculate  multipole <> multipole in reciprocal space we fourier 
!  transform the multipole arrays and then for each reciprocal vector m
!  we multiply  exp(-pi^2 m^2 / beta) times 1/m^2 times mpole^(m) * mpole^(-m)
!  where  mpole^(m) is the fourier transformed multipole array, and sum
!  over m (1/m^2 is replaced by (1-cos(2pi |m|D)/m^2 for finite clusters
!  with diameter D)
!
!  To make other calculations easier we reformulate this by letting
!  MC_grid be the fourier transform of multipoles  and
!  MD_grid be exp(-pi^2 m^2 / beta) times MC_grid(m), which is the same as
!  the fourier transform of a corresponding hermite gaussian with exponent
!  coefficient beta
!
!  Then the multipole <> multipole interaction in reciprocal space is
!  calculated by summing MC_grid(m)*MD_grid(-m)*(1/m^2) over m
!  again replace  1/m^2 by (1-cos(2pi |m|D)/m^2 for clusters)
!  Finally we simplify one more time:  MD_grid(m)*(1/m^2) is the
!  Fourier transformed electrostatic potential due to the corresponding 
!  hermite gaussians - we store it in MD_FT_potential(m)
!  similarly MC_FT_potential holds the Fourier transformed electrostatic 
!  potential due to the multipoles---then the reciprocal space interaction
!  is obtained by summing  MD_FT_potential(m)*MC_grid(-m) over m or
!  summing MC_FT_potential(m)*MD_grid(m)

!  Similarly for any 2 compact hermites we calculate their remaining
!  interactions in reciprocal space as if they were multipoles interacting
!  through erf(beta*rij) / rij in place of 1/rij  ---see X-ray encyclopedia
!   book chapter for details --since they all act like multipoles, we can
!  accumulate their coefficients into one array callled the sumCH array that
!  acts like a multipole array
!  Therefore we have  sumCH_C_grid and sumCH_D_grid and sumCH_C_FT_pot_grid
!  and sumCH_D_FT_pot_grid as above and we can interact
!  sumCH_D_FT_pot_grid with sumCH_C_grid to get the reciprocal space
!  compact hermite <> compact hermite interactions
!  and also sumCH_D_FT_pot_grid with MC_grid, MD_FT_potential with sumCH_C_grid
!  to get the multipole <> compact hermite and compact hermite <> multipole
!  interactions in reciprocal space

!  Finally the multipole <> diffuse hermite, diffuse hermite <> multipole
!  compact hermite <> diffuse hermite, diffuse hermite <> compact hermite
!  and  diffuse hermite <> diffuse hermite interactions are calculated
!  entirely in reciprocal space
!  For these neither the compact nor diffuse hermites can be treated as 
!  multipoles since there is no ewald adjustment of the interactions
!  Thus all the hermites have to be separately transformed
!  For pme this entails a number of grids stored in HC_grid and HD_grid
!  For ffp or regular ewald just one HC_grid and one HD_grid are needed
!  The gridded Fourier transformed values at m are summed into 
!  collect_HC_grid and collect_HD_grid and then multiplied by 1/m^2
!  (or (1-cos(2pi |m|D))/m^2 for clusters) to give HC_FT_pot_grid or
!  HD_FT_pot_grid, the electrostatic potential due to the compact or
!  diffuse hermite arrays. The reciprocal space interaction energies are
!  obtained as above by multiplying by the appropriate "charge" grid and
!  summing over m 
!---------------------------------------------------------------------
module recip_grids
   implicit none
   private

   
   ! grid type and fft stuff
   integer, save, allocatable :: nfftdim1_for_gridtype(:), &
                                 nfftdim2_for_gridtype(:), &
                                 nfftdim3_for_gridtype(:)
   integer, save, allocatable :: nfftable_for_gridtype(:), &
                                 nffwork_for_gridtype(:), &
                                 sizfftab_for_gridtype(:), &
                                 off_fftab_for_gridtype(:), &
                                 sizffwrk_for_gridtype(:)
   double precision,save, allocatable :: fftable(:),ffwork(:),fft_tmpy(:), &
                                         fft_alpha(:),fft_beta(:)

   ! mpole stuff
   integer, save              :: num_mpoles_for_mpole_grid, &
                                 tot_siz_Mpole_grid_Ffield, &
                                 tot_siz_Mpole_grid_Fcoeff
   integer, save, allocatable :: Fcoeff_offset_for_gridded_mpole(:), &
                                 Gcoeff_offset_for_gridded_mpole(:), &
                                 Mpole_order_for_gridded_mpole(:), &
                                 Ffield_offset_for_gridded_mpole(:), &
                                 Ffield_order_for_gridded_mpole(:)
   integer, save, allocatable :: site_that_owns_gridded_mpole(:)
   double precision, save, allocatable :: Mpole_grid_Fcoeff(:), &
                                          Mpole_grid_Ffield(:)
   integer, save               :: tot_siz_MPOLE_grid
   double precision,save, allocatable :: MC_grid(:), & 
                                         MD_grid(:) 
   integer, save              :: nfft1_for_MPOLE_collect, &
                                 nfft2_for_MPOLE_collect, &
                                 nfft3_for_MPOLE_collect, &
                                 nfftdim1_for_MPOLE_collect, &
                                 nfftdim2_for_MPOLE_collect, &
                                 nfftdim3_for_MPOLE_collect
   integer, save              :: nfft1_for_MPOLE_FT_pot, &
                                 nfft2_for_MPOLE_FT_pot, &
                                 nfft3_for_MPOLE_FT_pot, &
                                 nfftdim1_for_MPOLE_FT_pot, &
                                 nfftdim2_for_MPOLE_FT_pot, &
                                 nfftdim3_for_MPOLE_FT_pot
   integer, save              :: tot_siz_MPOLE_collect_grid
   integer, save              :: tot_siz_MPOLE_FT_potgrid
   double precision, save, allocatable :: collect_MC_grid(:), &
                                          collect_MD_grid(:), &
                                          MC_FT_pot_grid(:), &
                                          MD_FT_pot_grid(:)
   ! sumCH stuff
   integer, save              :: num_sumCH_for_sumCH_grid, &
                                 tot_siz_sumCH_grid_Ffield, &
                                 tot_siz_sumCH_grid_Fcoeff
   integer, save, allocatable :: Fcoeff_offset_for_gridded_sumCH(:), &
                                 Gcoeff_offset_for_gridded_sumCH(:), &
                                 sumCH_order_for_gridded_sumCH(:), &
                                 Ffield_offset_for_gridded_sumCH(:), &
                                 Ffield_order_for_gridded_sumCH(:)
   integer, save, allocatable :: site_that_owns_gridded_sumCH(:)
   double precision, save, allocatable :: sumCH_grid_Fcoeff(:), &
                                          sumCH_grid_Ffield(:)
   integer, save               :: tot_siz_sumCH_grid
   double precision,save, allocatable :: sumCH_C_grid(:), & 
                                         sumCH_D_grid(:) 
   integer, save              :: nfft1_for_sumCH_collect, &
                                 nfft2_for_sumCH_collect, &
                                 nfft3_for_sumCH_collect, &
                                 nfftdim1_for_sumCH_collect, &
                                 nfftdim2_for_sumCH_collect, &
                                 nfftdim3_for_sumCH_collect
   integer, save              :: nfft1_for_sumCH_FT_pot, &
                                 nfft2_for_sumCH_FT_pot, &
                                 nfft3_for_sumCH_FT_pot, &
                                 nfftdim1_for_sumCH_FT_pot, &
                                 nfftdim2_for_sumCH_FT_pot, &
                                 nfftdim3_for_sumCH_FT_pot
   integer, save              :: tot_siz_sumCH_collect_grid
   integer, save              :: tot_siz_sumCH_FT_potgrid
   double precision, save, allocatable :: collect_sumCH_C_grid(:), &
                                          collect_sumCH_D_grid(:), &
                                          sumCH_C_FT_pot_grid(:), &
                                          sumCH_D_FT_pot_grid(:)
   ! compact hermite stuff
   integer, save              :: tot_num_gridded_HC_prims
   integer, save, allocatable :: num_prims_for_HC_grid(:), &
                                 prim_offset_for_HC_grid(:)
   integer, save, allocatable :: prim_num_for_HC_grid_prim(:)
   integer, save, allocatable :: herm_order_for_HC_grid_prim(:)
   double precision, save, allocatable :: herm_expon_for_HC_grid_prim(:)
   integer, save, allocatable :: Ffield_order_for_HC_grid_prim(:)
   integer, save, allocatable :: Fcoeff_offset_for_HC_grid_prim(:)
   integer, save, allocatable :: Ffield_offset_for_HC_grid_prim(:)
   integer, save, allocatable :: Gcoeff_offset_for_HC_grid_prim(:)
   integer, save, allocatable :: site_that_owns_HC_grid_prim(:)
   integer, save              :: tot_siz_HC_grid_Fcoeff, &
                                 tot_siz_HC_grid_Ffield
   double precision, save, allocatable :: HC_grid_Fcoeff(:), &
                                          HC_grid_Ffield(:)
   integer, save, allocatable :: HC_grid_offset_for_HC_grid(:)
   integer, save               :: tot_siz_HC_grid
   double precision,save, allocatable :: HC_grid(:)
   integer, save              :: nfft1_for_HC_collect, &
                                 nfft2_for_HC_collect, &
                                 nfft3_for_HC_collect, &
                                 nfftdim1_for_HC_collect, &
                                 nfftdim2_for_HC_collect, &
                                 nfftdim3_for_HC_collect
   integer, save              :: nfft1_for_HC_FT_pot, &
                                 nfft2_for_HC_FT_pot, &
                                 nfft3_for_HC_FT_pot, &
                                 nfftdim1_for_HC_FT_pot, &
                                 nfftdim2_for_HC_FT_pot, &
                                 nfftdim3_for_HC_FT_pot
   integer, save              :: tot_siz_collect_HC_grid
   integer, save              :: tot_siz_HC_FT_pot_grid
   double precision, save, allocatable :: collect_HC_grid(:), &
                                          HC_FT_pot_grid(:)
   ! diffuse hermite stuff
   integer, save              :: tot_num_gridded_HD_prims
   integer, save, allocatable :: num_prims_for_HD_grid(:), &
                                 prim_offset_for_HD_grid(:)
   integer, save, allocatable :: prim_num_for_HD_grid_prim(:)
   integer, save, allocatable :: herm_order_for_HD_grid_prim(:)
   double precision, save, allocatable :: herm_expon_for_HD_grid_prim(:)
   integer, save, allocatable :: Ffield_order_for_HD_grid_prim(:)
   integer, save, allocatable :: Fcoeff_offset_for_HD_grid_prim(:)
   integer, save, allocatable :: Ffield_offset_for_HD_grid_prim(:)
   integer, save, allocatable :: Gcoeff_offset_for_HD_grid_prim(:)
   integer, save, allocatable :: site_that_owns_HD_grid_prim(:)
   integer, save              :: tot_siz_HD_grid_Fcoeff, &
                                 tot_siz_HD_grid_Ffield
   double precision, save, allocatable :: HD_grid_Fcoeff(:), &
                                          HD_grid_Ffield(:)
   integer, save, allocatable :: HD_grid_offset_for_HD_grid(:)
   integer, save               :: tot_siz_HD_grid
   double precision,save, allocatable :: HD_grid(:)
   integer, save              :: nfft1_for_HD_collect, &
                                 nfft2_for_HD_collect, &
                                 nfft3_for_HD_collect, &
                                 nfftdim1_for_HD_collect, &
                                 nfftdim2_for_HD_collect, &
                                 nfftdim3_for_HD_collect
   integer, save              :: nfft1_for_HD_FT_pot, &
                                 nfft2_for_HD_FT_pot, &
                                 nfft3_for_HD_FT_pot, &
                                 nfftdim1_for_HD_FT_pot, &
                                 nfftdim2_for_HD_FT_pot, &
                                 nfftdim3_for_HD_FT_pot
   integer, save              :: tot_siz_collect_HD_grid
   integer, save              :: tot_siz_HD_FT_pot_grid
   double precision, save, allocatable :: collect_HD_grid(:), &
                                          HD_FT_pot_grid(:)

   ! Coulomb grid
   integer, save :: nfft1_for_coulomb, &
                    nfft2_for_coulomb, &
                    nfft3_for_coulomb
   integer, save :: nfftdim1_for_coulomb, &
                    nfftdim2_for_coulomb, &
                    nfftdim3_for_coulomb, &
                    tot_siz_coulomb
   double precision, save, allocatable :: coulomb_grid(:)
   ! diffuse adjust grid
   integer, save :: nfft1_for_diff_adj, &
                    nfft2_for_diff_adj, &
                    nfft3_for_diff_adj
   integer, save :: nfftdim1_for_diff_adj, &
                    nfftdim2_for_diff_adj, &
                    nfftdim3_for_diff_adj, &
                    tot_siz_diff_adj
   double precision, save, allocatable :: diff_adj_grid(:)
   ! RECIPROCAL VIRIAL FOR MD
   double precision, save :: rec_virial(3,3)
!-----------------------------------------------------------------
public    FH_RECIP_GRIDS_setup, &
          FH_RECIP_GRIDS_deallocate, &
          FH_RECIP_GRIDS_fft_R_to_K, &
          FH_RECIP_GRIDS_collective, &
          FH_RECIP_GRIDS_mult_coulomb, &
          FH_RECIP_GRIDS_virial, &
          FH_RECIP_GRIDS_FT_ene, &
          RECIP_GRIDS_load_FT_pot, &
          FH_RECIP_GRIDS_fft_K_to_R, &
          ! grid type stuff
          nfftdim1_for_gridtype, &
          nfftdim2_for_gridtype, &
          nfftdim3_for_gridtype, &
          ! mpole stuff
          num_mpoles_for_mpole_grid, &
          tot_siz_MPOLE_grid, &
          mpole_order_for_gridded_mpole, &
          site_that_owns_gridded_mpole, &
          Fcoeff_offset_for_gridded_mpole, &
          Gcoeff_offset_for_gridded_mpole, &
          Ffield_offset_for_gridded_mpole, &
          Ffield_order_for_gridded_mpole, &
          Mpole_grid_Fcoeff, &
          Mpole_grid_Ffield, &
          MC_grid, &
          MD_grid, &
          !sumCH stuff
          num_sumCH_for_sumCH_grid, &
          tot_siz_sumCH_grid, &
          sumCH_order_for_gridded_sumCH, &
          site_that_owns_gridded_sumCH, &
          Fcoeff_offset_for_gridded_sumCH, &
          Gcoeff_offset_for_gridded_sumCH, &
          Ffield_offset_for_gridded_sumCH, &
          Ffield_order_for_gridded_sumCH, &
          sumCH_grid_Fcoeff, &
          sumCH_grid_Ffield, &
          sumCH_C_grid, &
          sumCH_D_grid, &
          ! compact hermite stuff
          tot_num_gridded_HC_prims, &
          num_prims_for_HC_grid, &
          prim_offset_for_HC_grid, &
          prim_num_for_HC_grid_prim, &
          herm_expon_for_HC_grid_prim, &
          herm_order_for_HC_grid_prim, &
          Gcoeff_offset_for_HC_grid_prim, &
          Fcoeff_offset_for_HC_grid_prim, &
          site_that_owns_HC_grid_prim, &
          HC_grid_offset_for_HC_grid, &
          HC_grid_Fcoeff, &
          HC_grid_Ffield, &
          Ffield_offset_for_HC_grid_prim, &
          Ffield_order_for_HC_grid_prim, &
          HC_grid, &
          ! diffuse hermite stuff
          tot_num_gridded_HD_prims, &
          num_prims_for_HD_grid, &
          prim_offset_for_HD_grid, &
          prim_num_for_HD_grid_prim, &
          herm_expon_for_HD_grid_prim, &
          herm_order_for_HD_grid_prim, &
          Gcoeff_offset_for_HD_grid_prim, &
          Fcoeff_offset_for_HD_grid_prim, &
          site_that_owns_HD_grid_prim, &
          HD_grid_offset_for_HD_grid, &
          HD_grid_Fcoeff, &
          HD_grid_Ffield, &
          Ffield_offset_for_HD_grid_prim, &
          Ffield_order_for_HD_grid_prim, &
          HD_grid,&
       ! VIRIAL
          rec_virial
          

contains
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!-----------------------------------------------------------------
! SETUP SUBROUTINES
!-----------------------------------------------------------------
!-----------------------------------------------------------------

!-----------------------------------------------------------------
subroutine FH_RECIP_GRIDS_setup(out_lun)

   implicit none

   integer, intent(in)  :: out_lun

   call FH_RECIP_GRIDS_setup_fft(out_lun)
   call FH_RECIP_GRIDS_setup_ptrs(out_lun)
   call FH_RECIP_GRIDS_setup_CD_grids(out_lun)
   call FH_RECIP_GRIDS_setup_collects(out_lun)
   call FH_RECIP_GRIDS_setup_coulomb(out_lun)
   call FH_RECIP_GRIDS_setup_diff_adj(out_lun)

end subroutine FH_RECIP_GRIDS_setup
!-----------------------------------------------------------------
subroutine FH_RECIP_GRIDS_setup_fft(out_lun)

   use user,only : num_prim_grid_types, &
                   nfft1_for_grid_type, &
                   nfft2_for_grid_type, &
                   nfft3_for_grid_type

   implicit none

   integer, intent(in)  :: out_lun

   integer :: ier1,ier2,ier3,ier4,ier5,ier6,ier7,ier8
   integer :: gt,off,maxsizffw,maxnfft1,maxnfftdim1,totsizfftab
   double precision :: dummy

   if(allocated(nfftdim1_for_gridtype)) deallocate(nfftdim1_for_gridtype)
   allocate(nfftdim1_for_gridtype(num_prim_grid_types),stat=ier1)

   if(allocated(nfftdim2_for_gridtype)) deallocate(nfftdim2_for_gridtype)
   allocate(nfftdim2_for_gridtype(num_prim_grid_types),stat=ier2)

   if(allocated(nfftdim3_for_gridtype)) deallocate(nfftdim3_for_gridtype)
   allocate(nfftdim3_for_gridtype(num_prim_grid_types),stat=ier3)

   if(allocated(nfftable_for_gridtype)) deallocate(nfftable_for_gridtype)
   allocate(nfftable_for_gridtype(num_prim_grid_types),stat=ier4)

   if(allocated(nffwork_for_gridtype)) deallocate(nffwork_for_gridtype)
   allocate(nffwork_for_gridtype(num_prim_grid_types),stat=ier5)

   if(allocated(sizfftab_for_gridtype)) deallocate(sizfftab_for_gridtype)
   allocate(sizfftab_for_gridtype(num_prim_grid_types),stat=ier6)

   if(allocated(off_fftab_for_gridtype)) deallocate(off_fftab_for_gridtype)
   allocate(off_fftab_for_gridtype(num_prim_grid_types),stat=ier7)

   if(allocated(sizffwrk_for_gridtype)) deallocate(sizffwrk_for_gridtype)
   allocate(sizffwrk_for_gridtype(num_prim_grid_types),stat=ier8)
   if ( ier1 /= 0 .or. ier2 /= 0 .or. ier3 /= 0 .or. ier4 /= 0 .or. &
        ier5 /= 0 .or. ier6 /= 0 .or. ier7 /= 0 .or. ier8 /= 0 )then
      write(out_lun,*)'FH_RECIP_fft_setup: allocate fails!'
      stop
   endif

   maxsizffw = -1
   maxnfft1 = -1
   maxnfftdim1 = -1
   do gt = 1,num_prim_grid_types
      call get_fftdims( &
                   nfft1_for_grid_type(gt), &
                   nfft2_for_grid_type(gt), &
                   nfft3_for_grid_type(gt),     &
                   nfftdim1_for_gridtype(gt), &
                   nfftdim2_for_gridtype(gt), &
                   nfftdim3_for_gridtype(gt), &
                   nfftable_for_gridtype(gt), &
                   nffwork_for_gridtype(gt), &
                   sizfftab_for_gridtype(gt), &
                   sizffwrk_for_gridtype(gt), &
                   out_lun)
      if (  sizffwrk_for_gridtype(gt) > maxsizffw ) &
             maxsizffw = sizffwrk_for_gridtype(gt)
      if (  nfftdim1_for_gridtype(gt) > maxnfftdim1 ) &
             maxnfftdim1 = nfftdim1_for_gridtype(gt)
      if (  nfft1_for_grid_type(gt) > maxnfft1 ) &
          maxnfft1 = nfft1_for_grid_type(gt)
   enddo

   off = 0
   do gt = 1,num_prim_grid_types
      off_fftab_for_gridtype(gt) = off
      off = off + sizfftab_for_gridtype(gt)
   enddo
   totsizfftab = off

   allocate(fftable(totsizfftab),stat=ier1)
   allocate(ffwork(maxsizffw),stat=ier2)
   allocate(fft_tmpy(2*maxnfftdim1),stat=ier3)
   allocate(fft_alpha(maxnfft1),stat=ier4)
   allocate(fft_beta(maxnfft1),stat=ier5)
   if ( ier1 /= 0 .or. ier2 /= 0 .or. ier3 /= 0 .or. ier4 /= 0 .or. &
        ier5 /= 0 )then
      write(out_lun,*)'FH_RECIP_fft_setup: allocate fails!'
      stop
   endif
   ! fill fftable
   do gt = 1,num_prim_grid_types
      off = off_fftab_for_gridtype(gt)
      call fft_setup_gem(dummy,fftable(off+1),ffwork,  &
                     nfft1_for_grid_type(gt), &
                     nfft2_for_grid_type(gt), &
                     nfft3_for_grid_type(gt),   &
                     nfftdim1_for_gridtype(gt), &
                     nfftdim2_for_gridtype(gt), &
                     nfftdim3_for_gridtype(gt),    &
                     nfftable_for_gridtype(gt), &
                     nffwork_for_gridtype(gt), &
                     out_lun)
   enddo
end subroutine FH_RECIP_GRIDS_setup_fft
!-----------------------------------------------------------------
subroutine FH_RECIP_GRIDS_setup_ptrs(out_lun)

   use user, only : num_HC_prim_grids,num_HD_prim_grids,verbose
   use sites, only : num_sites
   use hermite, only : tot_num_mpole_coeffs, &
                       mpole_level_for_site, &
                       mpole_coeff_off_for_site, &
                       tot_num_sumCH_coeffs, &
                       sumCH_level_for_site, &
                       sumCH_coeff_off_for_site, &
                       tot_num_herm_Cprims, &
                       hermite_level_for_Cprim, &
                       hermite_rec_expon_for_Cprim, &
                       herm_coeff_offset_of_Cprim, &
                       site_that_owns_Cprim, &
                       HC_grid_number_for_Cprim, &
                       tot_num_herm_Dprims, &
                       hermite_level_for_Dprim, &
                       hermite_rec_expon_for_Dprim, &
                       herm_coeff_offset_of_Dprim, &
                       site_that_owns_Dprim, &
                       HD_grid_number_for_Dprim

   implicit none

   integer, intent(in)  :: out_lun

   include "interact_type.fh"
   integer :: ier1,ier2,ier3,ier4,ier5,ier6,ier7,ier8
   integer :: n,ngp,np,gc,gd,gt,offc,offd,lev_prim,lev_mp
   integer :: mp_order_for_level(-1:5)

   mp_order_for_level(-1) = 0
   mp_order_for_level(0) = 1
   mp_order_for_level(1) = 4
   mp_order_for_level(2) = 10
   mp_order_for_level(3) = 20
   mp_order_for_level(4) = 35
   mp_order_for_level(5) = 56

   ! first mpole stuff
   if ( tot_num_mpole_coeffs > 0 )then
      ngp = 0
      do n = 1, num_sites
         lev_mp = mpole_level_for_site(n)
         if ( lev_mp >= 0 )then
            ngp = ngp + 1
         endif
      enddo
      num_mpoles_for_mpole_grid = ngp
      if (allocated(Ffield_order_for_gridded_mpole)) &
              deallocate(Ffield_order_for_gridded_mpole)
      allocate( Ffield_order_for_gridded_mpole(num_mpoles_for_mpole_grid), &
                    stat=ier1)

      if (allocated(Ffield_offset_for_gridded_mpole)) &
              deallocate(Ffield_offset_for_gridded_mpole)
      allocate( Ffield_offset_for_gridded_mpole(num_mpoles_for_mpole_grid), &
                    stat=ier2)

      if (allocated(site_that_owns_gridded_mpole)) &
              deallocate(site_that_owns_gridded_mpole)
      allocate( site_that_owns_gridded_mpole(num_mpoles_for_mpole_grid), &
                    stat=ier3)

      if (allocated(Mpole_order_for_gridded_mpole)) &
              deallocate(Mpole_order_for_gridded_mpole)
      allocate( Mpole_order_for_gridded_mpole(num_mpoles_for_mpole_grid), &
                    stat=ier4)

      if (allocated(Fcoeff_offset_for_gridded_mpole)) &
              deallocate(Fcoeff_offset_for_gridded_mpole)
      allocate( Fcoeff_offset_for_gridded_mpole(num_mpoles_for_mpole_grid), &
                    stat=ier5)

      if (allocated(Gcoeff_offset_for_gridded_mpole)) &
              deallocate(Gcoeff_offset_for_gridded_mpole)
      allocate( Gcoeff_offset_for_gridded_mpole(num_mpoles_for_mpole_grid), &
                    stat=ier6)
      if ( ier1 /= 0 .or. ier2 /= 0 .or. ier3 /= 0 )then
          write(out_lun,*)'FH_RECIP_GRIDS_setup_ptrs: allocate fails!'
          stop
      endif
      ! next the site_that_owns_gridded_mpole
      ngp = 0
      do n = 1, num_sites
         lev_mp = mpole_level_for_site(n)
         if ( lev_mp >= 0 )then
            ngp = ngp + 1
            site_that_owns_gridded_mpole(ngp) = n
            Mpole_order_for_gridded_mpole(ngp) = mp_order_for_level(lev_mp)
            Ffield_order_for_gridded_mpole(ngp) = mp_order_for_level(lev_mp+1)
            Gcoeff_offset_for_gridded_mpole(ngp) = mpole_coeff_off_for_site(n)
         endif
      enddo
      ! get the Fcoeff, Ffield offsets for mpoles
      Fcoeff_offset_for_gridded_mpole(1) = 0
      Ffield_offset_for_gridded_mpole(1) = 0
      do n = 2,num_mpoles_for_mpole_grid
         Fcoeff_offset_for_gridded_mpole(n) =  &
                                   Fcoeff_offset_for_gridded_mpole(n-1) + &
                                   Mpole_order_for_gridded_mpole(n-1)
         Ffield_offset_for_gridded_mpole(n) =  &
                                   Ffield_offset_for_gridded_mpole(n-1) + &
                                   Ffield_order_for_gridded_mpole(n-1)
      enddo
      tot_siz_Mpole_grid_Fcoeff = &
            Fcoeff_offset_for_gridded_mpole(num_mpoles_for_mpole_grid) + &
            Mpole_order_for_gridded_mpole(num_mpoles_for_mpole_grid)
      tot_siz_Mpole_grid_Ffield = &
            Ffield_offset_for_gridded_mpole(num_mpoles_for_mpole_grid) + &
            Ffield_order_for_gridded_mpole(num_mpoles_for_mpole_grid)
      ! allocate coeff and field arrays
      allocate( Mpole_grid_Fcoeff(tot_siz_Mpole_grid_Fcoeff),stat=ier1)
      allocate( Mpole_grid_Ffield(tot_siz_Mpole_grid_Ffield),stat=ier2)
      if ( ier1 /= 0 .or. ier2 /= 0 )then
         write(out_lun,*)'FH_RECIP_GRIDS_setup_ptrs: allocate fails!'
         stop
      endif
   else
      num_mpoles_for_mpole_grid = 0
   endif !( tot_num_mpole_coeffs > 0 

   ! next sumCH stuff
   if ( tot_num_sumCH_coeffs > 0 )then
      ngp = 0
      do n = 1, num_sites
         lev_mp = sumCH_level_for_site(n)
         if ( lev_mp >= 0 )then
            ngp = ngp + 1
         endif
      enddo
      num_sumCH_for_sumCH_grid = ngp
      if (allocated(Ffield_order_for_gridded_sumCH)) &
              deallocate(Ffield_order_for_gridded_sumCH)
      allocate( Ffield_order_for_gridded_sumCH(num_sumCH_for_sumCH_grid), &
                    stat=ier1)

      if (allocated(Ffield_offset_for_gridded_sumCH)) &
              deallocate(Ffield_offset_for_gridded_sumCH)
      allocate( Ffield_offset_for_gridded_sumCH(num_sumCH_for_sumCH_grid), &
                    stat=ier2)

      if (allocated(site_that_owns_gridded_sumCH)) &
              deallocate(site_that_owns_gridded_sumCH)
      allocate( site_that_owns_gridded_sumCH(num_sumCH_for_sumCH_grid), &
                    stat=ier3)

      if (allocated(sumCH_order_for_gridded_sumCH)) &
              deallocate(sumCH_order_for_gridded_sumCH)
      allocate( sumCH_order_for_gridded_sumCH(num_sumCH_for_sumCH_grid), &
                    stat=ier4)

      if (allocated(Fcoeff_offset_for_gridded_sumCH)) &
              deallocate(Fcoeff_offset_for_gridded_sumCH)
      allocate( Fcoeff_offset_for_gridded_sumCH(num_sumCH_for_sumCH_grid), &
                    stat=ier5)

      if (allocated(Gcoeff_offset_for_gridded_sumCH)) &
              deallocate(Gcoeff_offset_for_gridded_sumCH)
      allocate( Gcoeff_offset_for_gridded_sumCH(num_sumCH_for_sumCH_grid), &
                    stat=ier6)
      if ( ier1 /= 0 .or. ier2 /= 0 .or. ier3 /= 0 )then
          write(out_lun,*)'FH_RECIP_GRIDS_setup_ptrs: allocate fails!'
          stop
      endif
      ! next the site_that_owns_gridded_sumCH
      ngp = 0
      do n = 1, num_sites
         lev_mp = sumCH_level_for_site(n)
         if ( lev_mp >= 0 )then
            ngp = ngp + 1
            site_that_owns_gridded_sumCH(ngp) = n
            sumCH_order_for_gridded_sumCH(ngp) = mp_order_for_level(lev_mp)
            Ffield_order_for_gridded_sumCH(ngp) = mp_order_for_level(lev_mp+1)
            Gcoeff_offset_for_gridded_sumCH(ngp) = sumCH_coeff_off_for_site(n)
         endif
      enddo
      ! get the Fcoeff, Ffield offsets for sumCH
      Fcoeff_offset_for_gridded_sumCH(1) = 0
      Ffield_offset_for_gridded_sumCH(1) = 0
      do n = 2,num_sumCH_for_sumCH_grid
         Fcoeff_offset_for_gridded_sumCH(n) =  &
                                   Fcoeff_offset_for_gridded_sumCH(n-1) + &
                                   sumCH_order_for_gridded_sumCH(n-1)
         Ffield_offset_for_gridded_sumCH(n) =  &
                                   Ffield_offset_for_gridded_sumCH(n-1) + &
                                   Ffield_order_for_gridded_sumCH(n-1)
      enddo
      tot_siz_sumCH_grid_Fcoeff = &
            Fcoeff_offset_for_gridded_sumCH(num_sumCH_for_sumCH_grid) + &
            sumCH_order_for_gridded_sumCH(num_sumCH_for_sumCH_grid)
      tot_siz_sumCH_grid_Ffield = &
            Ffield_offset_for_gridded_sumCH(num_sumCH_for_sumCH_grid) + &
            Ffield_order_for_gridded_sumCH(num_sumCH_for_sumCH_grid)
      ! allocate coeff and field arrays
      if (allocated(sumCH_grid_Fcoeff)) &
              deallocate(sumCH_grid_Fcoeff)
      allocate( sumCH_grid_Fcoeff(tot_siz_sumCH_grid_Fcoeff),stat=ier1)

      if (allocated(sumCH_grid_Ffield)) &
              deallocate(sumCH_grid_Ffield)
      allocate( sumCH_grid_Ffield(tot_siz_sumCH_grid_Ffield),stat=ier2)
      if ( ier1 /= 0 .or. ier2 /= 0 )then
         write(out_lun,*)'FH_RECIP_GRIDS_setup_ptrs: allocate fails!'
         stop
      endif
   else
      num_sumCH_for_sumCH_grid = 0
   endif !( tot_num_sumCH_coeffs > 0 )then

   ! next compact hermite stuff
   if ( tot_num_herm_Cprims > 0 )then
      if (allocated(num_prims_for_HC_grid)) &
              deallocate(num_prims_for_HC_grid)
      allocate( num_prims_for_HC_grid(num_HC_prim_grids),stat=ier1)

      if (allocated(prim_offset_for_HC_grid)) &
              deallocate(prim_offset_for_HC_grid)
      allocate( prim_offset_for_HC_grid(num_HC_prim_grids),stat=ier2)
      if ( ier1 /= 0 .or. ier2 /= 0 )then
         write(out_lun,*)'FH_RECIP_GRIDS_setup_ptrs: allocate fails!'
         stop
      endif
      ! get num_prims_for_HC grid
      num_prims_for_HC_grid = 0
      do np = 1,tot_num_herm_Cprims
         gc = HC_grid_number_for_Cprim(np)
         if ( (gc == 0) .or. (gc > num_HC_prim_grids) )then
            write(out_lun,*)'FH_RECIP_GRIDS_setup_ptrs: hermite Cprim ',np, &
                      'has bad compact grid number!'
            stop
         else
            num_prims_for_HC_grid(gc) = num_prims_for_HC_grid(gc) + 1
         endif
      enddo ! np = 1,tot_num_herm_Cprims
      ! next the HC grid offsets
      prim_offset_for_HC_grid(1) = 0
      do gc = 2,num_HC_prim_grids
         prim_offset_for_HC_grid(gc) = prim_offset_for_HC_grid(gc-1) + &
                                       num_prims_for_HC_grid(gc-1)
      enddo
      tot_num_gridded_HC_prims =  &
                              prim_offset_for_HC_grid(num_HC_prim_grids) + &
                              num_prims_for_HC_grid(num_HC_prim_grids)
      ! allocate pointers for HC grids
      if (allocated(prim_num_for_HC_grid_prim)) &
              deallocate(prim_num_for_HC_grid_prim)
      allocate( prim_num_for_HC_grid_prim(tot_num_gridded_HC_prims), &
                                          stat=ier1)

      if (allocated(herm_order_for_HC_grid_prim)) &
              deallocate(herm_order_for_HC_grid_prim)
      allocate( herm_order_for_HC_grid_prim(tot_num_gridded_HC_prims), &
                                          stat=ier2)

      if (allocated(herm_expon_for_HC_grid_prim)) &
              deallocate(herm_expon_for_HC_grid_prim)
      allocate( herm_expon_for_HC_grid_prim(tot_num_gridded_HC_prims), &
                                          stat=ier3)

      if (allocated(Fcoeff_offset_for_HC_grid_prim)) &
              deallocate(Fcoeff_offset_for_HC_grid_prim)
      allocate( Fcoeff_offset_for_HC_grid_prim(tot_num_gridded_HC_prims), &
                                          stat=ier4)

      if (allocated(Gcoeff_offset_for_HC_grid_prim)) &
              deallocate(Gcoeff_offset_for_HC_grid_prim)
      allocate( Gcoeff_offset_for_HC_grid_prim(tot_num_gridded_HC_prims), &
                                          stat=ier5)

      if (allocated(Ffield_order_for_HC_grid_prim)) &
              deallocate(Ffield_order_for_HC_grid_prim)
      allocate( Ffield_order_for_HC_grid_prim(tot_num_gridded_HC_prims), &
                                          stat=ier6)

      if (allocated(Ffield_offset_for_HC_grid_prim)) &
              deallocate(Fcoeff_offset_for_HC_grid_prim)
      allocate( Ffield_offset_for_HC_grid_prim(tot_num_gridded_HC_prims), &
                                          stat=ier7)

      if (allocated(site_that_owns_HC_grid_prim)) &
              deallocate(site_that_owns_HC_grid_prim)
      allocate( site_that_owns_HC_grid_prim(tot_num_gridded_HC_prims), &
                                          stat=ier8)
      if ( ier1 /= 0 .or. ier2 /= 0 .or. ier3 /= 0 .or. ier4 /= 0 &
                     .or. ier5 /= 0 .or. ier6 /= 0 .or. ier7 /= 0 &
                     .or. ier8 /= 0 )then
         write(out_lun,*)'FH_RECIP_GRIDS_setup_ptrs: allocate fails!'
         stop
      endif
      num_prims_for_HC_grid = 0
      do np = 1,tot_num_herm_Cprims
         gc = HC_grid_number_for_Cprim(np)
         if ( gc > 0 )then
            num_prims_for_HC_grid(gc) = num_prims_for_HC_grid(gc) + 1
            offc = prim_offset_for_HC_grid(gc)
            prim_num_for_HC_grid_prim(offc + num_prims_for_HC_grid(gc)) = np
            lev_prim = hermite_level_for_Cprim(np)
            herm_order_for_HC_grid_prim(offc + num_prims_for_HC_grid(gc)) = &
                     mp_order_for_level(lev_prim)
            Ffield_order_for_HC_grid_prim(offc + num_prims_for_HC_grid(gc)) = &
                     mp_order_for_level(lev_prim+1)
            herm_expon_for_HC_grid_prim(offc + num_prims_for_HC_grid(gc)) = &
                     hermite_rec_expon_for_Cprim(np)
            Gcoeff_offset_for_HC_grid_prim(offc + num_prims_for_HC_grid(gc)) = &
                     herm_coeff_offset_of_Cprim(np)
            n = site_that_owns_Cprim(np)
            site_that_owns_HC_grid_prim(offc + num_prims_for_HC_grid(gc)) = n
         endif ! ( gc > 0 )
      enddo !np = 1,tot_num_herm_Cprims
      ! get the Fcoeff, Ffield offsets for HC grids
      Fcoeff_offset_for_HC_grid_prim(1) = 0
      Ffield_offset_for_HC_grid_prim(1) = 0
      do ngp = 2,tot_num_gridded_HC_prims
         Fcoeff_offset_for_HC_grid_prim(ngp) = &
                            Fcoeff_offset_for_HC_grid_prim(ngp-1) + &
                            herm_order_for_HC_grid_prim(ngp-1)
         Ffield_offset_for_HC_grid_prim(ngp) = &
                            Ffield_offset_for_HC_grid_prim(ngp-1) + &
                            Ffield_order_for_HC_grid_prim(ngp-1)
      enddo
      tot_siz_HC_grid_Fcoeff = &
          Fcoeff_offset_for_HC_grid_prim(tot_num_gridded_HC_prims) + &
          herm_order_for_HC_grid_prim(tot_num_gridded_HC_prims)
      tot_siz_HC_grid_Ffield = &
          Ffield_offset_for_HC_grid_prim(tot_num_gridded_HC_prims) + &
          Ffield_order_for_HC_grid_prim(tot_num_gridded_HC_prims)
      ! allocate coeff and field arrays
      if (allocated(HC_grid_Fcoeff)) &
              deallocate(HC_grid_Fcoeff)
      allocate( HC_grid_Fcoeff(tot_siz_HC_grid_Fcoeff),stat=ier1)
      if (allocated(HC_grid_Ffield)) &
              deallocate(HC_grid_Ffield)
      allocate( HC_grid_Ffield(tot_siz_HC_grid_Ffield),stat=ier2)
      if ( ier1 /= 0 .or. ier2 /= 0 )then
         write(out_lun,*)'FH_RECIP_GRIDS_setup_ptrs: allocate fails!'
         stop
      endif
   else
      tot_num_gridded_HC_prims = 0
   endif !( tot_num_herm_Cprims > 0 )then

   ! next diffuse hermite stuff
   if ( tot_num_herm_Dprims > 0 )then
      if (allocated(num_prims_for_HD_grid)) &
              deallocate(num_prims_for_HD_grid)
      allocate( num_prims_for_HD_grid(num_HD_prim_grids),stat=ier1)
      if (allocated(prim_offset_for_HD_grid)) &
              deallocate(prim_offset_for_HD_grid)
      allocate( prim_offset_for_HD_grid(num_HD_prim_grids),stat=ier2)
      if ( ier1 /= 0 .or. ier2 /= 0 )then
         write(out_lun,*)'FH_RECIP_GRIDS_setup_ptrs: allocate fails!'
         stop
      endif
      ! get num_prims_for_HD grid
      num_prims_for_HD_grid = 0
      do np = 1,tot_num_herm_Dprims
         gd = HD_grid_number_for_Dprim(np)
         if ( (gd == 0) .or. (gd > num_HD_prim_grids) )then
            write(out_lun,*)'FH_RECIP_GRIDS_setup_ptrs: hermite Dprim ',np, &
                      'has bad diffuse grid number!'
            stop
         else
            num_prims_for_HD_grid(gd) = num_prims_for_HD_grid(gd) + 1
         endif
      enddo ! np = 1,tot_num_herm_Dprims
      ! next the HD grid offsets
      prim_offset_for_HD_grid(1) = 0
      do gd = 2,num_HD_prim_grids
         prim_offset_for_HD_grid(gd) = prim_offset_for_HD_grid(gd-1) + &
                                       num_prims_for_HD_grid(gd-1)
      enddo
      tot_num_gridded_HD_prims =  &
                              prim_offset_for_HD_grid(num_HD_prim_grids) + &
                              num_prims_for_HD_grid(num_HD_prim_grids)
      ! allocate pointers for HD grids
      if (allocated(prim_num_for_HD_grid_prim)) &
              deallocate(prim_num_for_HD_grid_prim)
      allocate( prim_num_for_HD_grid_prim(tot_num_gridded_HD_prims), &
                                          stat=ier1)

      if (allocated(herm_order_for_HD_grid_prim)) &
              deallocate(herm_order_for_HD_grid_prim)
      allocate( herm_order_for_HD_grid_prim(tot_num_gridded_HD_prims), &
                                          stat=ier2)

      if (allocated(herm_expon_for_HD_grid_prim)) &
              deallocate(herm_expon_for_HD_grid_prim)
      allocate( herm_expon_for_HD_grid_prim(tot_num_gridded_HD_prims), &
                                          stat=ier3)

      if (allocated(Fcoeff_offset_for_HD_grid_prim)) &
              deallocate(Fcoeff_offset_for_HD_grid_prim)
      allocate( Fcoeff_offset_for_HD_grid_prim(tot_num_gridded_HD_prims), &
                                          stat=ier4)

      if (allocated(Gcoeff_offset_for_HD_grid_prim)) &
              deallocate(Gcoeff_offset_for_HD_grid_prim)
      allocate( Gcoeff_offset_for_HD_grid_prim(tot_num_gridded_HD_prims), &
                                          stat=ier5)

      if (allocated(Ffield_order_for_HD_grid_prim)) &
              deallocate(Ffield_order_for_HD_grid_prim)
      allocate( Ffield_order_for_HD_grid_prim(tot_num_gridded_HD_prims), &
                                          stat=ier6)

      if (allocated(Ffield_offset_for_HD_grid_prim)) &
              deallocate(Ffield_offset_for_HD_grid_prim)
      allocate( Ffield_offset_for_HD_grid_prim(tot_num_gridded_HD_prims), &
                                          stat=ier7)

      if (allocated(site_that_owns_HD_grid_prim)) &
              deallocate(site_that_owns_HD_grid_prim)
      allocate( site_that_owns_HD_grid_prim(tot_num_gridded_HD_prims), &
                                          stat=ier8)
      if ( ier1 /= 0 .or. ier2 /= 0 .or. ier3 /= 0 .or. ier4 /= 0 &
                     .or. ier5 /= 0 .or. ier6 /= 0 .or. ier7 /= 0 &
                     .or. ier8 /= 0 )then
         write(out_lun,*)'FH_RECIP_GRIDS_setup_ptrs: allocate fails!'
         stop
      endif
      num_prims_for_HD_grid = 0
      do np = 1,tot_num_herm_Dprims
         gd = HD_grid_number_for_Dprim(np)
         if ( gd > 0 )then
            num_prims_for_HD_grid(gd) = num_prims_for_HD_grid(gd) + 1
            offd = prim_offset_for_HD_grid(gd)
            prim_num_for_HD_grid_prim(offd + num_prims_for_HD_grid(gd)) = np
            lev_prim = hermite_level_for_Dprim(np)
            herm_order_for_HD_grid_prim(offd + num_prims_for_HD_grid(gd)) = &
                     mp_order_for_level(lev_prim)
            Ffield_order_for_HD_grid_prim(offd + num_prims_for_HD_grid(gd)) = &
                     mp_order_for_level(lev_prim+1)
            herm_expon_for_HD_grid_prim(offd + num_prims_for_HD_grid(gd)) = &
                     hermite_rec_expon_for_Dprim(np)
            Gcoeff_offset_for_HD_grid_prim(offd + num_prims_for_HD_grid(gd)) = &
                     herm_coeff_offset_of_Dprim(np)
            n = site_that_owns_Dprim(np)
            site_that_owns_HD_grid_prim(offd + num_prims_for_HD_grid(gd)) = n
         endif ! ( gc > 0 )
      enddo !np = 1,tot_num_herm_Cprims
      ! get the Fcoeff, Ffield offsets for HD grids
      Fcoeff_offset_for_HD_grid_prim(1) = 0
      Ffield_offset_for_HD_grid_prim(1) = 0
      do ngp = 2,tot_num_gridded_HD_prims
         Fcoeff_offset_for_HD_grid_prim(ngp) = &
                            Fcoeff_offset_for_HD_grid_prim(ngp-1) + &
                            herm_order_for_HD_grid_prim(ngp-1)
         Ffield_offset_for_HD_grid_prim(ngp) = &
                            Ffield_offset_for_HD_grid_prim(ngp-1) + &
                            Ffield_order_for_HD_grid_prim(ngp-1)
      enddo
      tot_siz_HD_grid_Fcoeff = &
          Fcoeff_offset_for_HD_grid_prim(tot_num_gridded_HD_prims) + &
          herm_order_for_HD_grid_prim(tot_num_gridded_HD_prims)
      tot_siz_HD_grid_Ffield = &
          Ffield_offset_for_HD_grid_prim(tot_num_gridded_HD_prims) + &
          Ffield_order_for_HD_grid_prim(tot_num_gridded_HD_prims)
      if (allocated(HD_grid_Fcoeff)) &
              deallocate(HD_grid_Fcoeff)
      allocate( HD_grid_Fcoeff(tot_siz_HD_grid_Fcoeff),stat=ier1)
      if (allocated(HD_grid_Ffield)) &
              deallocate(HD_grid_Ffield)
      allocate( HD_grid_Ffield(tot_siz_HD_grid_Ffield),stat=ier2)
      if ( ier1 /= 0 .or. ier2 /= 0 )then
         write(out_lun,*)'FH_RECIP_GRIDS_setup_ptrs: allocate fails!'
         stop
      endif
   else
      tot_num_gridded_HD_prims = 0
   endif !( tot_num_herm_Dprims > 0 )then
   if ( verbose == 1 )then
      write(out_lun,*)'num_mpoles_for_mpole_grid = ',num_mpoles_for_mpole_grid
      write(out_lun,*)'num_sumCH_for_sumCH_grid = ',num_sumCH_for_sumCH_grid
      write(out_lun,*)'tot_num_gridded_HC_prims = ',tot_num_gridded_HC_prims
      write(out_lun,*)'tot_num_gridded_HD_prims = ',tot_num_gridded_HD_prims
   endif

end subroutine FH_RECIP_GRIDS_setup_ptrs
!-----------------------------------------------------------------
subroutine FH_RECIP_GRIDS_setup_CD_grids(out_lun)

   use user, only : num_HC_prim_grids,num_HD_prim_grids, &
                    grid_type_for_HC_prim_grid, &
                    grid_type_for_HD_prim_grid, &
                    grid_type_for_MPOLES, &
                    grid_type_for_sumCH
   use hermite, only : tot_num_mpole_coeffs, &
                       tot_num_sumCH_coeffs, &
                       tot_num_herm_Cprims, &
                       tot_num_herm_Dprims

   implicit none

   integer, intent(in)  :: out_lun

   integer ier,ier1,ier2
   integer off,g,gt

   ! multipole stuff
   if ( tot_num_mpole_coeffs > 0 )then
      gt = grid_type_for_MPOLES
      tot_siz_MPOLE_grid = 2*nfftdim1_for_gridtype(gt)* &
               nfftdim2_for_gridtype(gt)*nfftdim3_for_gridtype(gt)
      if ( tot_siz_MPOLE_grid > 0 )then
         if (allocated(MC_grid)) &
              deallocate(MC_grid)
         allocate(MC_grid(tot_siz_MPOLE_grid),stat=ier1)
         if (allocated(MD_grid)) &
              deallocate(MD_grid)
         allocate(MD_grid(tot_siz_MPOLE_grid),stat=ier2)
         if ( ier1 /= 0 .or. ier2 /= 0 )then
            write(out_lun,*)'FH_RECIP_GRIDS_setup_CD_grids: allocate fails!'
            stop
         endif
      endif
   else
      tot_siz_MPOLE_grid = 0
   endif

   ! sumCH stuff
   if ( tot_num_sumCH_coeffs > 0 )then
      gt = grid_type_for_sumCH
      tot_siz_sumCH_grid = 2*nfftdim1_for_gridtype(gt)* &
               nfftdim2_for_gridtype(gt)*nfftdim3_for_gridtype(gt)
      if ( tot_siz_sumCH_grid > 0 )then
         if (allocated(sumCH_C_grid)) &
              deallocate(sumCH_C_grid)
         allocate(sumCH_C_grid(tot_siz_sumCH_grid),stat=ier1)
         if (allocated(sumCH_D_grid)) &
              deallocate(sumCH_D_grid)
         allocate(sumCH_D_grid(tot_siz_sumCH_grid),stat=ier2)
         if ( ier1 /= 0 .or. ier2 /= 0 )then
            write(out_lun,*)'FH_RECIP_GRIDS_setup_CD_grids: allocate fails!'
            stop
         endif
      endif
   else
      tot_siz_sumCH_grid = 0
   endif

   ! compact hermite stuff
   if ( tot_num_herm_Cprims > 0 )then
         if (allocated(HC_grid_offset_for_HC_grid)) &
              deallocate(HC_grid_offset_for_HC_grid)
      allocate(HC_grid_offset_for_HC_grid(num_HC_prim_grids),stat=ier)
      if ( ier /= 0 )then
         write(out_lun,*)'FH_RECIP_GRIDS_setup_CD_grids: allocate fails!'
         stop
      endif
      off = 0
      do g = 1,num_HC_prim_grids
         HC_grid_offset_for_HC_grid(g) = off
         gt = grid_type_for_HC_prim_grid(g)
         off = off + 2*nfftdim1_for_gridtype(gt)* &
                       nfftdim2_for_gridtype(gt)* &
                       nfftdim3_for_gridtype(gt)
      enddo
      tot_siz_HC_grid = off
      if ( tot_siz_HC_grid > 0 )then
         if (allocated(HC_grid)) &
              deallocate(HC_grid)
         allocate(HC_grid(tot_siz_HC_grid),stat=ier)
         if ( ier /= 0 )then
            write(out_lun,*)'FH_RECIP_GRIDS_setup_CD_grids: allocate fails!'
            stop
         endif
      endif
   else
      tot_siz_HC_grid = 0
   endif

   ! diffuse hermite stuff
   if ( tot_num_herm_Dprims > 0 )then
         if (allocated(HD_grid_offset_for_HD_grid)) &
              deallocate(HD_grid_offset_for_HD_grid)
      allocate(HD_grid_offset_for_HD_grid(num_HD_prim_grids),stat=ier)
      if ( ier /= 0 )then
         write(out_lun,*)'FH_RECIP_GRIDS_setup_CD_grids: allocate fails!'
         stop
      endif
      off = 0
      do g = 1,num_HD_prim_grids
         HD_grid_offset_for_HD_grid(g) = off
         gt = grid_type_for_HD_prim_grid(g)
         off = off + 2*nfftdim1_for_gridtype(gt)* &
                       nfftdim2_for_gridtype(gt)* &
                       nfftdim3_for_gridtype(gt)
      enddo
      tot_siz_HD_grid = off
      if ( tot_siz_HD_grid > 0 )then
         if (allocated(HD_grid)) &
              deallocate(HD_grid)
         allocate(HD_grid(tot_siz_HD_grid),stat=ier)
         if ( ier /= 0 )then
            write(out_lun,*)'FH_RECIP_GRIDS_setup_CD_grids: allocate fails!'
            stop
         endif
      endif
   else
      tot_siz_HD_grid = 0
   endif

end subroutine FH_RECIP_GRIDS_setup_CD_grids
!-----------------------------------------------------------------
subroutine FH_RECIP_GRIDS_setup_collects(out_lun)

   use user, only : grid_type_for_MPOLES, &
                    grid_type_for_sumCH, &
                    num_HC_prim_grids, &
                    grid_type_for_HC_prim_grid, &
                    num_HD_prim_grids, &
                    grid_type_for_HD_prim_grid, &
                    nfft1_for_grid_type, &
                    nfft2_for_grid_type, &
                    nfft3_for_grid_type
   use hermite, only : tot_num_mpole_coeffs, &
                       tot_num_sumCH_coeffs, &
                       tot_num_herm_Cprims, &
                       tot_num_herm_Dprims

   implicit none

   integer, intent(in)  :: out_lun

   integer g,gt,nfft,nffw,sizfft,sizffw,ier1,ier2,ier3,ier4
   integer, parameter :: TOO_BIG = 1000000


   ! first do the mpole collects and FT potentials
   if ( tot_num_mpole_coeffs > 0 )then
      gt = grid_type_for_MPOLES
      nfft1_for_MPOLE_collect = nfft1_for_grid_type(gt)
      nfft2_for_MPOLE_collect = nfft2_for_grid_type(gt)
      nfft3_for_MPOLE_collect = nfft3_for_grid_type(gt)
      call get_fftdims(nfft1_for_MPOLE_collect, &
                       nfft2_for_MPOLE_collect, &
                       nfft3_for_MPOLE_collect, &
                       nfftdim1_for_MPOLE_collect, &
                       nfftdim2_for_MPOLE_collect, &
                       nfftdim3_for_MPOLE_collect, &
                       nfft,nffw,sizfft,sizffw, &
                       out_lun)
      tot_siz_MPOLE_collect_grid = 2*nfftdim1_for_MPOLE_collect* &
                                  nfftdim2_for_MPOLE_collect* &
                                  nfftdim3_for_MPOLE_collect
      nfft1_for_MPOLE_FT_pot = nfft1_for_MPOLE_collect
      nfft2_for_MPOLE_FT_pot = nfft2_for_MPOLE_collect
      nfft3_for_MPOLE_FT_pot = nfft3_for_MPOLE_collect
      nfftdim1_for_MPOLE_FT_pot = nfftdim1_for_MPOLE_collect
      nfftdim2_for_MPOLE_FT_pot = nfftdim2_for_MPOLE_collect
      nfftdim3_for_MPOLE_FT_pot = nfftdim3_for_MPOLE_collect
      tot_siz_MPOLE_FT_potgrid = tot_siz_MPOLE_collect_grid
      allocate(collect_MC_grid(tot_siz_MPOLE_collect_grid),stat=ier1)
      allocate(collect_MD_grid(tot_siz_MPOLE_collect_grid),stat=ier2)
      allocate(MC_FT_pot_grid(tot_siz_MPOLE_FT_potgrid),stat=ier3)
      allocate(MD_FT_pot_grid(tot_siz_MPOLE_FT_potgrid),stat=ier4)
      if ( ier1 /= 0 .or. ier2 /= 0 .or. ier3 /= 0 .or. ier4 /= 0 )then
         write(out_lun,*)'FH_RECIP_GRIDS_setup_collects: allocate fails!'
         stop
      endif
   else
      tot_siz_MPOLE_collect_grid = 0
      tot_siz_MPOLE_FT_potgrid = tot_siz_MPOLE_collect_grid
      nfft1_for_MPOLE_collect = -1
      nfft2_for_MPOLE_collect = -1
      nfft3_for_MPOLE_collect = -1
      nfft1_for_MPOLE_FT_pot = nfft1_for_MPOLE_collect
      nfft2_for_MPOLE_FT_pot = nfft2_for_MPOLE_collect
      nfft3_for_MPOLE_FT_pot = nfft3_for_MPOLE_collect
   endif

   ! next do the sumCH collects and FT potentials
   if ( tot_num_sumCH_coeffs > 0 )then
      gt = grid_type_for_sumCH
      nfft1_for_sumCH_collect = nfft1_for_grid_type(gt)
      nfft2_for_sumCH_collect = nfft2_for_grid_type(gt)
      nfft3_for_sumCH_collect = nfft3_for_grid_type(gt)
      call get_fftdims(nfft1_for_sumCH_collect, &
                       nfft2_for_sumCH_collect, &
                       nfft3_for_sumCH_collect, &
                       nfftdim1_for_sumCH_collect, &
                       nfftdim2_for_sumCH_collect, &
                       nfftdim3_for_sumCH_collect, &
                       nfft,nffw,sizfft,sizffw, &
                       out_lun)
      tot_siz_sumCH_collect_grid = 2*nfftdim1_for_sumCH_collect* &
                                  nfftdim2_for_sumCH_collect* &
                                  nfftdim3_for_sumCH_collect
      nfft1_for_sumCH_FT_pot = nfft1_for_sumCH_collect
      nfft2_for_sumCH_FT_pot = nfft2_for_sumCH_collect
      nfft3_for_sumCH_FT_pot = nfft3_for_sumCH_collect
      nfftdim1_for_sumCH_FT_pot = nfftdim1_for_sumCH_collect
      nfftdim2_for_sumCH_FT_pot = nfftdim2_for_sumCH_collect
      nfftdim3_for_sumCH_FT_pot = nfftdim3_for_sumCH_collect
      tot_siz_sumCH_FT_potgrid = tot_siz_sumCH_collect_grid
      allocate(collect_sumCH_C_grid(tot_siz_sumCH_collect_grid),stat=ier1)
      allocate(collect_sumCH_D_grid(tot_siz_sumCH_collect_grid),stat=ier2)
      allocate(sumCH_C_FT_pot_grid(tot_siz_sumCH_FT_potgrid),stat=ier3)
      allocate(sumCH_D_FT_pot_grid(tot_siz_sumCH_FT_potgrid),stat=ier4)
      if ( ier1 /= 0 .or. ier2 /= 0 .or. ier3 /= 0 .or. ier4 /= 0 )then
         write(out_lun,*)'FH_RECIP_GRIDS_setup_collects: allocate fails!'
         stop
      endif
   else
      tot_siz_sumCH_collect_grid = 0
      tot_siz_sumCH_FT_potgrid = tot_siz_sumCH_collect_grid
      nfft1_for_sumCH_collect = -1
      nfft2_for_sumCH_collect = -1
      nfft3_for_sumCH_collect = -1
      nfft1_for_sumCH_FT_pot = nfft1_for_sumCH_collect
      nfft2_for_sumCH_FT_pot = nfft2_for_sumCH_collect
      nfft3_for_sumCH_FT_pot = nfft3_for_sumCH_collect
   endif

   ! next do the HC_collects and FT potentials
   ! collective grid dims are minimal over those grids contributing
   if ( tot_num_herm_Cprims > 0 )then
      if ( num_HC_prim_grids == 0 )then
          write(out_lun,*)'compact prims but no compact prim grids!'
          stop
      endif
      if ( num_HC_prim_grids > 0 )then
         nfft1_for_HC_collect = TOO_BIG
         nfft2_for_HC_collect = TOO_BIG
         nfft3_for_HC_collect = TOO_BIG
      endif
      do g = 1,num_HC_prim_grids
         gt = grid_type_for_HC_prim_grid(g)
         if ( nfft1_for_grid_type(gt) < nfft1_for_HC_collect ) &
                nfft1_for_HC_collect = nfft1_for_grid_type(gt)
         if ( nfft2_for_grid_type(gt) < nfft2_for_HC_collect ) &
                nfft2_for_HC_collect = nfft2_for_grid_type(gt)
         if ( nfft3_for_grid_type(gt) < nfft3_for_HC_collect ) &
                nfft3_for_HC_collect = nfft3_for_grid_type(gt)
      enddo
      call get_fftdims(nfft1_for_HC_collect, &
                       nfft2_for_HC_collect, &
                       nfft3_for_HC_collect, &
                       nfftdim1_for_HC_collect, &
                       nfftdim2_for_HC_collect, &
                       nfftdim3_for_HC_collect, &
                       nfft,nffw,sizfft,sizffw, &
                       out_lun)
      tot_siz_collect_HC_grid = 2*nfftdim1_for_HC_collect* &
                                  nfftdim2_for_HC_collect* &
                                  nfftdim3_for_HC_collect
      nfft1_for_HC_FT_pot = nfft1_for_HC_collect
      nfft2_for_HC_FT_pot = nfft2_for_HC_collect
      nfft3_for_HC_FT_pot = nfft3_for_HC_collect
      nfftdim1_for_HC_FT_pot = nfftdim1_for_HC_collect
      nfftdim2_for_HC_FT_pot = nfftdim2_for_HC_collect
      nfftdim3_for_HC_FT_pot = nfftdim3_for_HC_collect
      tot_siz_HC_FT_pot_grid = tot_siz_collect_HC_grid
      allocate(collect_HC_grid(tot_siz_collect_HC_grid),stat=ier1)
      allocate(HC_FT_pot_grid(tot_siz_HC_FT_pot_grid),stat=ier2)
      if ( ier1 /= 0 .or. ier2 /= 0 )then
         write(out_lun,*)'FH_RECIP_GRIDS_setup_collects: allocate fails!'
         stop
      endif
   else
      tot_siz_collect_HC_grid = 0
      tot_siz_HC_FT_pot_grid = tot_siz_collect_HC_grid
      nfft1_for_HC_collect = -1
      nfft2_for_HC_collect = -1
      nfft3_for_HC_collect = -1
      nfft1_for_HC_FT_pot = nfft1_for_HC_collect
      nfft2_for_HC_FT_pot = nfft2_for_HC_collect
      nfft3_for_HC_FT_pot = nfft3_for_HC_collect
   endif !( num_HC_prim_grids > 0 )

   ! next do the HD_collects and FT potentials
   ! collective grid dims are minimal over those grids contributing
   if ( tot_num_herm_Dprims > 0 )then
      if ( num_HD_prim_grids == 0 )then
          write(out_lun,*)'compact prims but no compact prim grids!'
          stop
      endif
      if ( num_HD_prim_grids > 0 )then
         nfft1_for_HD_collect = TOO_BIG
         nfft2_for_HD_collect = TOO_BIG
         nfft3_for_HD_collect = TOO_BIG
      endif
      do g = 1,num_HD_prim_grids
         gt = grid_type_for_HD_prim_grid(g)
         if ( nfft1_for_grid_type(gt) < nfft1_for_HD_collect ) &
                nfft1_for_HD_collect = nfft1_for_grid_type(gt)
         if ( nfft2_for_grid_type(gt) < nfft2_for_HD_collect ) &
                nfft2_for_HD_collect = nfft2_for_grid_type(gt)
         if ( nfft3_for_grid_type(gt) < nfft3_for_HD_collect ) &
                nfft3_for_HD_collect = nfft3_for_grid_type(gt)
      enddo
      call get_fftdims(nfft1_for_HD_collect, &
                       nfft2_for_HD_collect, &
                       nfft3_for_HD_collect, &
                       nfftdim1_for_HD_collect, &
                       nfftdim2_for_HD_collect, &
                       nfftdim3_for_HD_collect, &
                       nfft,nffw,sizfft,sizffw, &
                       out_lun)
      tot_siz_collect_HD_grid = 2*nfftdim1_for_HD_collect* &
                                  nfftdim2_for_HD_collect* &
                                  nfftdim3_for_HD_collect
      nfft1_for_HD_FT_pot = nfft1_for_HD_collect
      nfft2_for_HD_FT_pot = nfft2_for_HD_collect
      nfft3_for_HD_FT_pot = nfft3_for_HD_collect
      nfftdim1_for_HD_FT_pot = nfftdim1_for_HD_collect
      nfftdim2_for_HD_FT_pot = nfftdim2_for_HD_collect
      nfftdim3_for_HD_FT_pot = nfftdim3_for_HD_collect
      tot_siz_HD_FT_pot_grid = tot_siz_collect_HD_grid
      if(allocated(collect_HD_grid)) deallocate(collect_HD_grid)
      allocate(collect_HD_grid(tot_siz_collect_HD_grid),stat=ier1)
      if(allocated(HD_FT_pot_grid)) deallocate(HD_FT_pot_grid)
      allocate(HD_FT_pot_grid(tot_siz_HD_FT_pot_grid),stat=ier2)
      if ( ier1 /= 0 .or. ier2 /= 0 )then
         write(out_lun,*)'FH_RECIP_GRIDS_setup_collects: allocate fails!'
         stop
      endif
   else
      tot_siz_collect_HD_grid = 0
      tot_siz_HD_FT_pot_grid = tot_siz_collect_HD_grid
      nfft1_for_HD_collect = -1
      nfft2_for_HD_collect = -1
      nfft3_for_HD_collect = -1
      nfft1_for_HD_FT_pot = nfft1_for_HD_collect
      nfft2_for_HD_FT_pot = nfft2_for_HD_collect
      nfft3_for_HD_FT_pot = nfft3_for_HD_collect
   endif !( num_HD_prim_grids > 0 )

end subroutine FH_RECIP_GRIDS_setup_collects
!-----------------------------------------------------------------
subroutine FH_RECIP_GRIDS_setup_coulomb(out_lun)

   use user,only : pbc
   use unit_cell, only : recip,sphere

   implicit none

   integer, intent(in)  :: out_lun

   integer :: ncg,nfft,nffw,sizfft,sizffw,ier

   nfft1_for_coulomb = -1
   nfft2_for_coulomb = -1
   nfft3_for_coulomb = -1

   ! coulomb grid dims are maximum over collective grids
   ! mpoles
   if ( nfft1_for_MPOLE_collect > nfft1_for_coulomb ) &
           nfft1_for_coulomb = nfft1_for_MPOLE_collect
   if ( nfft2_for_MPOLE_collect > nfft2_for_coulomb ) &
           nfft2_for_coulomb = nfft2_for_MPOLE_collect
   if ( nfft3_for_MPOLE_collect > nfft3_for_coulomb ) &
           nfft3_for_coulomb = nfft3_for_MPOLE_collect
   ! sumCH
   if ( nfft1_for_sumCH_collect > nfft1_for_coulomb ) &
           nfft1_for_coulomb = nfft1_for_sumCH_collect
   if ( nfft2_for_sumCH_collect > nfft2_for_coulomb ) &
           nfft2_for_coulomb = nfft2_for_sumCH_collect
   if ( nfft3_for_sumCH_collect > nfft3_for_coulomb ) &
           nfft3_for_coulomb = nfft3_for_sumCH_collect
   ! compact hermite
   if ( nfft1_for_HC_collect > nfft1_for_coulomb ) &
           nfft1_for_coulomb = nfft1_for_HC_collect
   if ( nfft2_for_HC_collect > nfft2_for_coulomb ) &
           nfft2_for_coulomb = nfft2_for_HC_collect
   if ( nfft3_for_HC_collect > nfft3_for_coulomb ) &
           nfft3_for_coulomb = nfft3_for_HC_collect
   ! diffuse hermite
   if ( nfft1_for_HD_collect > nfft1_for_coulomb ) &
           nfft1_for_coulomb = nfft1_for_HD_collect
   if ( nfft2_for_HD_collect > nfft2_for_coulomb ) &
           nfft2_for_coulomb = nfft2_for_HD_collect
   if ( nfft3_for_HD_collect > nfft3_for_coulomb ) &
           nfft3_for_coulomb = nfft3_for_HD_collect
   call get_fftdims(nfft1_for_coulomb, &
                    nfft2_for_coulomb, &
                    nfft3_for_coulomb, &
                    nfftdim1_for_coulomb, &
                    nfftdim2_for_coulomb, &
                    nfftdim3_for_coulomb, &
                    nfft,nffw,sizfft,sizffw, &
                    out_lun)
   tot_siz_coulomb = nfft3_for_coulomb* &
                     nfftdim1_for_coulomb* &
                     nfft2_for_coulomb
   
   if(allocated(coulomb_grid)) deallocate(coulomb_grid)
   allocate(coulomb_grid(tot_siz_coulomb),stat=ier)
   if ( ier /= 0 )then
      write(out_lun,*)'FH_RECIP_coulomb_setup: allocate fails!'
      stop
   endif
   call FH_RECIP_GRIDS_get_FT_coulomb( &
                           nfft1_for_coulomb, &
                           nfft2_for_coulomb, &
                           nfft3_for_coulomb, &
                           nfftdim1_for_coulomb, &
                           pbc, sphere,recip,coulomb_grid)

end subroutine FH_RECIP_GRIDS_setup_coulomb
!-----------------------------------------------------------------
subroutine FH_RECIP_GRIDS_setup_diff_adj(out_lun)

   use user,only : diffuse_adjust_exponent
   use unit_cell, only : recip

   implicit none

   integer, intent(in)  :: out_lun

   integer :: ncg,nfft,nffw,sizfft,sizffw,ier
   include 'scale.fh'

   if ( diffuse_adjust_exponent <= small_value )return ! no need for this grid

   nfft1_for_diff_adj = -1
   nfft2_for_diff_adj = -1
   nfft3_for_diff_adj = -1

   ! diff_adj grid dims are maximum over collective grids
   ! mpoles
   if ( nfft1_for_MPOLE_collect > nfft1_for_diff_adj ) &
           nfft1_for_diff_adj = nfft1_for_MPOLE_collect
   if ( nfft2_for_MPOLE_collect > nfft2_for_diff_adj ) &
           nfft2_for_diff_adj = nfft2_for_MPOLE_collect
   if ( nfft3_for_MPOLE_collect > nfft3_for_diff_adj ) &
           nfft3_for_diff_adj = nfft3_for_MPOLE_collect
   ! sumCH
   if ( nfft1_for_sumCH_collect > nfft1_for_diff_adj ) &
           nfft1_for_diff_adj = nfft1_for_sumCH_collect
   if ( nfft2_for_sumCH_collect > nfft2_for_diff_adj ) &
           nfft2_for_diff_adj = nfft2_for_sumCH_collect
   if ( nfft3_for_sumCH_collect > nfft3_for_diff_adj ) &
           nfft3_for_diff_adj = nfft3_for_sumCH_collect
   ! compact hermite
   if ( nfft1_for_HC_collect > nfft1_for_diff_adj ) &
           nfft1_for_diff_adj = nfft1_for_HC_collect
   if ( nfft2_for_HC_collect > nfft2_for_diff_adj ) &
           nfft2_for_diff_adj = nfft2_for_HC_collect
   if ( nfft3_for_HC_collect > nfft3_for_diff_adj ) &
           nfft3_for_diff_adj = nfft3_for_HC_collect
   ! diffuse hermite
   if ( nfft1_for_HD_collect > nfft1_for_diff_adj ) &
           nfft1_for_diff_adj = nfft1_for_HD_collect
   if ( nfft2_for_HD_collect > nfft2_for_diff_adj ) &
           nfft2_for_diff_adj = nfft2_for_HD_collect
   if ( nfft3_for_HD_collect > nfft3_for_diff_adj ) &
           nfft3_for_diff_adj = nfft3_for_HD_collect
   call get_fftdims(nfft1_for_diff_adj, &
                    nfft2_for_diff_adj, &
                    nfft3_for_diff_adj, &
                    nfftdim1_for_diff_adj, &
                    nfftdim2_for_diff_adj, &
                    nfftdim3_for_diff_adj, &
                    nfft,nffw,sizfft,sizffw)
   tot_siz_diff_adj = nfft3_for_diff_adj* &
                      nfftdim1_for_diff_adj* &
                      nfft2_for_diff_adj
   if(allocated(diff_adj_grid)) deallocate(diff_adj_grid)
   allocate(diff_adj_grid(tot_siz_diff_adj),stat=ier)
   if ( ier /= 0 )then
      write(out_lun,*)'FH_RECIP_diff_adj_setup: allocate fails!'
      stop
   endif

   call FH_RECIP_GRIDS_diff_adj_grid( &
                    nfft1_for_diff_adj, &
                    nfft2_for_diff_adj, &
                    nfft3_for_diff_adj, &
                    nfftdim1_for_diff_adj, &
                    diffuse_adjust_exponent, &
                    recip,diff_adj_grid)

end subroutine FH_RECIP_GRIDS_setup_diff_adj
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!-----------------------------------------------------------------
! EVALUATION SUBROUTINES
!-----------------------------------------------------------------
!-----------------------------------------------------------------

!-----------------------------------------------------------------
subroutine FH_RECIP_GRIDS_fft_R_to_K(out_lun)

   use user,only : num_HC_prim_grids, &
                   grid_type_for_HC_prim_grid, &
                   num_HD_prim_grids, &
                   grid_type_for_HD_prim_grid, &
                   grid_type_for_MPOLES, &
                   grid_type_for_sumCH, &
                   struc_fac_method_for_grid_type, &
                   num_pme_grid_types, &
                   num_ffp_grid_types, &
                   nfft1_for_grid_type, &
                   nfft2_for_grid_type, &
                   nfft3_for_grid_type
   use hermite, only : tot_num_mpole_coeffs, &
                       tot_num_sumCH_coeffs, &
                       tot_num_herm_Cprims, &
                       tot_num_herm_Dprims

   implicit none

   integer, intent(in)  :: out_lun

   include 'structure_factor_type.fh'

   integer :: g,gt,off1,off2

   if ( (num_pme_grid_types == 0) .and. &
        (num_ffp_grid_types == 0)  )return ! nothing further to do here

   ! first the MC_grid
   if ( tot_num_mpole_coeffs > 0 )then
      gt = grid_type_for_MPOLES
      if ( (struc_fac_method_for_grid_type(gt) == SF_PME) .OR. &
           (struc_fac_method_for_grid_type(gt) == SF_FFP) )then
         off1 = 0 !just one grid
         off2 = off_fftab_for_gridtype(gt)
         call fft_backrc( &
               MC_grid(off1+1), &
               fftable(off2+1), &
               ffwork, &
               nfft1_for_grid_type(gt), &
               nfft2_for_grid_type(gt), &
               nfft3_for_grid_type(gt), &
               nfftdim1_for_gridtype(gt), &
               nfftdim2_for_gridtype(gt), &
               nfftdim3_for_gridtype(gt), &
               nfftable_for_gridtype(gt), &
               nffwork_for_gridtype(gt), &
               fft_tmpy, &
               fft_alpha, &
               fft_beta, &
               out_lun)
      endif !( (struc_fac_method_for_grid_type(gt) == SF_PME) .OR. &
      ! do the diffuse multipole if SF_FFP (no need if SF_PME)
      if ( struc_fac_method_for_grid_type(gt) == SF_FFP )then
         call fft_backrc( &
               MD_grid(off1+1), &
               fftable(off2+1), &
               ffwork, &
               nfft1_for_grid_type(gt), &
               nfft2_for_grid_type(gt), &
               nfft3_for_grid_type(gt), &
               nfftdim1_for_gridtype(gt), &
               nfftdim2_for_gridtype(gt), &
               nfftdim3_for_gridtype(gt), &
               nfftable_for_gridtype(gt), &
               nffwork_for_gridtype(gt), &
               fft_tmpy, &
               fft_alpha, &
               fft_beta, &
               out_lun)
      endif
   endif !( tot_num_mpole_coeffs > 0 )then
   
   ! next the sumCH_C_grid
   if ( tot_num_sumch_coeffs > 0 )then
      gt = grid_type_for_sumCH
      if ( (struc_fac_method_for_grid_type(gt) == SF_PME) .OR. &
           (struc_fac_method_for_grid_type(gt) == SF_FFP) )then
         off1 = 0 !just one grid
         off2 = off_fftab_for_gridtype(gt)
         call fft_backrc( &
               sumCH_C_grid(off1+1), &
               fftable(off2+1), &
               ffwork, &
               nfft1_for_grid_type(gt), &
               nfft2_for_grid_type(gt), &
               nfft3_for_grid_type(gt), &
               nfftdim1_for_gridtype(gt), &
               nfftdim2_for_gridtype(gt), &
               nfftdim3_for_gridtype(gt), &
               nfftable_for_gridtype(gt), &
               nffwork_for_gridtype(gt), &
               fft_tmpy, &
               fft_alpha, &
               fft_beta, &
               out_lun)
      endif !( (struc_fac_method_for_grid_type(gt) == SF_PME) .OR. &
      ! do the diffuse sumCH if SF_FFP (no need if SF_PME)
      if ( struc_fac_method_for_grid_type(gt) == SF_FFP )then
         call fft_backrc( &
               sumCH_D_grid(off1+1), &
               fftable(off2+1), &
               ffwork, &
               nfft1_for_grid_type(gt), &
               nfft2_for_grid_type(gt), &
               nfft3_for_grid_type(gt), &
               nfftdim1_for_gridtype(gt), &
               nfftdim2_for_gridtype(gt), &
               nfftdim3_for_gridtype(gt), &
               nfftable_for_gridtype(gt), &
               nffwork_for_gridtype(gt), &
               fft_tmpy, &
               fft_alpha, &
               fft_beta, &
               out_lun)
      endif
   endif !( tot_num_sumch_coeffs > 0 )then
   
   ! next the compact hermites
   if ( tot_num_herm_Cprims > 0 )then
      do g = 1,num_HC_prim_grids
         gt = grid_type_for_HC_prim_grid(g)
         if ( (struc_fac_method_for_grid_type(gt) == SF_PME) .OR. &
              (struc_fac_method_for_grid_type(gt) == SF_FFP) )then
            off1 = HC_grid_offset_for_HC_grid(g)
            off2 = off_fftab_for_gridtype(gt)
            call fft_backrc( &
               HC_grid(off1+1), &
               fftable(off2+1), &
               ffwork, &
               nfft1_for_grid_type(gt), &
               nfft2_for_grid_type(gt), &
               nfft3_for_grid_type(gt), &
               nfftdim1_for_gridtype(gt), &
               nfftdim2_for_gridtype(gt), &
               nfftdim3_for_gridtype(gt), &
               nfftable_for_gridtype(gt), &
               nffwork_for_gridtype(gt), &
               fft_tmpy, &
               fft_alpha, &
               fft_beta, &
               out_lun)
         endif !( (struc_fac_method_for_grid_type(gt) == SF_PME) .OR. &
      enddo !g = 1,num_HC_prim_grids
   endif !( tot_num_herm_Cprims > 0 )then

   ! finally the diffuse hermites
   if ( tot_num_herm_Dprims > 0 )then
      do g = 1,num_HD_prim_grids
         gt = grid_type_for_HD_prim_grid(g)
         if ( (struc_fac_method_for_grid_type(gt) == SF_PME) .OR. &
              (struc_fac_method_for_grid_type(gt) == SF_FFP) )then
            off1 = HD_grid_offset_for_HD_grid(g)
            off2 = off_fftab_for_gridtype(gt)
            call fft_backrc( &
               HD_grid(off1+1), &
               fftable(off2+1), &
               ffwork, &
               nfft1_for_grid_type(gt), &
               nfft2_for_grid_type(gt), &
               nfft3_for_grid_type(gt), &
               nfftdim1_for_gridtype(gt), &
               nfftdim2_for_gridtype(gt), &
               nfftdim3_for_gridtype(gt), &
               nfftable_for_gridtype(gt), &
               nffwork_for_gridtype(gt), &
               fft_tmpy,&
               fft_alpha, &
               fft_beta, &
               out_lun)
         endif !( (struc_fac_method_for_grid_type(gt) == SF_PME) .OR. &
      enddo !g = 1,num_HD_prim_grids
   endif !( tot_num_herm_Dprims > 0 )then

end subroutine FH_RECIP_GRIDS_fft_R_to_K
!-----------------------------------------------------------------
subroutine FH_RECIP_GRIDS_fft_K_to_R(out_lun)

   use user,only : num_HC_prim_grids, &
                   grid_type_for_HC_prim_grid, &
                   num_HD_prim_grids, &
                   grid_type_for_HD_prim_grid, &
                   grid_type_for_MPOLES, &
                   grid_type_for_sumCH, &
                   struc_fac_method_for_grid_type, &
                   num_pme_grid_types, &
                   num_ffp_grid_types, &
                   nfft1_for_grid_type, &
                   nfft2_for_grid_type, &
                   nfft3_for_grid_type
   use hermite, only : tot_num_mpole_coeffs, &
                       tot_num_sumCH_coeffs, &
                       tot_num_herm_Cprims, &
                       tot_num_herm_Dprims

   implicit none

   integer, intent(in)  :: out_lun

   include 'structure_factor_type.fh'

   integer :: g,gt,off1,off2

   if ( (num_pme_grid_types == 0) .and. &
        (num_ffp_grid_types == 0)  )return ! nothing further to do here

   ! first the MC_grid
   if ( tot_num_mpole_coeffs > 0 )then
      gt = grid_type_for_MPOLES
      if ( (struc_fac_method_for_grid_type(gt) == SF_PME) .OR. &
           (struc_fac_method_for_grid_type(gt) == SF_FFP) )then
         off1 = 0 !just one grid
         off2 = off_fftab_for_gridtype(gt)
         call fft_forwardrc( &
               MC_grid(off1+1), &
               fftable(off2+1), &
               ffwork, &
               nfft1_for_grid_type(gt), &
               nfft2_for_grid_type(gt), &
               nfft3_for_grid_type(gt), &
               nfftdim1_for_gridtype(gt), &
               nfftdim2_for_gridtype(gt), &
               nfftdim3_for_gridtype(gt), &
               nfftable_for_gridtype(gt), &
               nffwork_for_gridtype(gt), &
               fft_tmpy, &
               fft_alpha, &
               fft_beta, &
               out_lun)
      endif !( (struc_fac_method_for_grid_type(gt) == SF_PME) .OR. &
   endif !( tot_num_mpole_coeffs > 0 )then

   ! next the sumCH_C_grid
   if ( tot_num_sumch_coeffs > 0 )then
      gt = grid_type_for_sumCH
      if ( (struc_fac_method_for_grid_type(gt) == SF_PME) .OR. &
           (struc_fac_method_for_grid_type(gt) == SF_FFP) )then
         off1 = 0 !just one grid
         off2 = off_fftab_for_gridtype(gt)
         call fft_forwardrc( &
               sumCH_C_grid(off1+1), &
               fftable(off2+1), &
               ffwork, &
               nfft1_for_grid_type(gt), &
               nfft2_for_grid_type(gt), &
               nfft3_for_grid_type(gt), &
               nfftdim1_for_gridtype(gt), &
               nfftdim2_for_gridtype(gt), &
               nfftdim3_for_gridtype(gt), &
               nfftable_for_gridtype(gt), &
               nffwork_for_gridtype(gt), &
               fft_tmpy, &
               fft_alpha, &
               fft_beta, &
               out_lun)
      endif !( (struc_fac_method_for_grid_type(gt) == SF_PME) .OR. &
   endif !( tot_num_sumch_coeffs > 0 )then

   ! next the compact hermites
   if ( tot_num_herm_Cprims > 0 )then
      do g = 1,num_HC_prim_grids
         gt = grid_type_for_HC_prim_grid(g)
         if ( (struc_fac_method_for_grid_type(gt) == SF_PME) .OR. &
              (struc_fac_method_for_grid_type(gt) == SF_FFP) )then
            off1 = HC_grid_offset_for_HC_grid(g)
            off2 = off_fftab_for_gridtype(gt)
            call fft_forwardrc( &
               HC_grid(off1+1), &
               fftable(off2+1), &
               ffwork, &
               nfft1_for_grid_type(gt), &
               nfft2_for_grid_type(gt), &
               nfft3_for_grid_type(gt), &
               nfftdim1_for_gridtype(gt), &
               nfftdim2_for_gridtype(gt), &
               nfftdim3_for_gridtype(gt), &
               nfftable_for_gridtype(gt), &
               nffwork_for_gridtype(gt), &
               fft_tmpy, &
               fft_alpha, &
               fft_beta, &
               out_lun)
         endif !( (struc_fac_method_for_grid_type(gt) == SF_PME) .OR. &
      enddo !g = 1,num_HC_prim_grids
   endif !( tot_num_herm_Cprims > 0 )then

   ! finally the diffuse hermites
   if ( tot_num_herm_Dprims > 0 )then
      do g = 1,num_HD_prim_grids
         gt = grid_type_for_HD_prim_grid(g)
         if ( (struc_fac_method_for_grid_type(gt) == SF_PME) .OR. &
              (struc_fac_method_for_grid_type(gt) == SF_FFP) )then
            off1 = HD_grid_offset_for_HD_grid(g)
            off2 = off_fftab_for_gridtype(gt)
            call fft_forwardrc( &
               HD_grid(off1+1), &
               fftable(off2+1), &
               ffwork, &
               nfft1_for_grid_type(gt), &
               nfft2_for_grid_type(gt), &
               nfft3_for_grid_type(gt), &
               nfftdim1_for_gridtype(gt), &
               nfftdim2_for_gridtype(gt), &
               nfftdim3_for_gridtype(gt), &
               nfftable_for_gridtype(gt), &
               nffwork_for_gridtype(gt), &
               fft_tmpy,&
               fft_alpha, &
               fft_beta, &
               out_lun)
         endif !( (struc_fac_method_for_grid_type(gt) == SF_PME) .OR. &
      enddo !g = 1,num_HD_prim_grids
   endif !( tot_num_herm_Dprims > 0 )then

end subroutine FH_RECIP_GRIDS_fft_K_to_R
!-----------------------------------------------------------------
subroutine FH_RECIP_GRIDS_collective(out_lun)

   use user, only : num_HC_prim_grids, &
                    num_HD_prim_grids, &
                    grid_type_for_HC_prim_grid, &
                    grid_type_for_HD_prim_grid, &
                    grid_type_for_MPOLES, &
                    grid_type_for_sumCH, &
                    nfft1_for_grid_type, &
                    nfft2_for_grid_type, &
                    nfft3_for_grid_type
   use hermite, only : tot_num_mpole_coeffs, &
                       tot_num_sumCH_coeffs, &
                       tot_num_herm_Cprims, &
                       tot_num_herm_Dprims

   implicit none

   integer, intent(in)  :: out_lun

   integer g,gt,off

   ! first fill MC collective
   if ( tot_num_mpole_coeffs > 0 )then
      call UTIL_zero_real_array(collect_MC_grid,tot_siz_MPOLE_collect_grid)
      call UTIL_zero_real_array(collect_MD_grid,tot_siz_MPOLE_collect_grid)
      gt = grid_type_for_MPOLES
      off = 0
      call FH_RECIP_GRIDS_accum( &
          nfft1_for_MPOLE_collect, &
          nfft2_for_MPOLE_collect, &
          nfft3_for_MPOLE_collect, &
          nfftdim1_for_MPOLE_collect, &
          nfft1_for_grid_type(gt), &
          nfft2_for_grid_type(gt), &
          nfft3_for_grid_type(gt), &
          nfftdim1_for_gridtype(gt), &
          collect_MC_grid,MC_grid(off+1),out_lun)
      ! next fill MD collective
      call FH_RECIP_GRIDS_accum( &
          nfft1_for_MPOLE_collect, &
          nfft2_for_MPOLE_collect, &
          nfft3_for_MPOLE_collect, &
          nfftdim1_for_MPOLE_collect, &
          nfft1_for_grid_type(gt), &
          nfft2_for_grid_type(gt), &
          nfft3_for_grid_type(gt), &
          nfftdim1_for_gridtype(gt), &
          collect_MD_grid,MD_grid(off+1),out_lun)
   endif !( tot_num_mpole_coeffs > 0 )then

   ! next fill sumCH_C collective
   if ( tot_num_sumCH_coeffs > 0 )then
      call UTIL_zero_real_array(collect_sumCH_C_grid, &
                               tot_siz_sumCH_collect_grid)
      call UTIL_zero_real_array(collect_sumCH_D_grid, &
                               tot_siz_sumCH_collect_grid)
      gt = grid_type_for_sumCH
      off = 0
      call FH_RECIP_GRIDS_accum( &
          nfft1_for_sumCH_collect, &
          nfft2_for_sumCH_collect, &
          nfft3_for_sumCH_collect, &
          nfftdim1_for_sumCH_collect, &
          nfft1_for_grid_type(gt), &
          nfft2_for_grid_type(gt), &
          nfft3_for_grid_type(gt), &
          nfftdim1_for_gridtype(gt), &
          collect_sumCH_C_grid, &
          sumCH_C_grid(off+1),out_lun)
      ! next fill sumCH_D collective
      call FH_RECIP_GRIDS_accum( &
          nfft1_for_sumCH_collect, &
          nfft2_for_sumCH_collect, &
          nfft3_for_sumCH_collect, &
          nfftdim1_for_sumCH_collect, &
          nfft1_for_grid_type(gt), &
          nfft2_for_grid_type(gt), &
          nfft3_for_grid_type(gt), &
          nfftdim1_for_gridtype(gt), &
          collect_sumCH_D_grid, &
          sumCH_D_grid(off+1),out_lun)
   endif !( tot_num_sumCH_Coeffs > 0 )then

   ! next fill HC collective
   if ( tot_num_herm_Cprims > 0 )then
      do g = 1,num_HC_prim_grids
         off = HC_grid_offset_for_HC_grid(g)
         gt = grid_type_for_HC_prim_grid(g)
         call FH_RECIP_GRIDS_accum( &
             nfft1_for_HC_collect, &
             nfft2_for_HC_collect, &
             nfft3_for_HC_collect, &
             nfftdim1_for_HC_collect, &
             nfft1_for_grid_type(gt), &
             nfft2_for_grid_type(gt), &
             nfft3_for_grid_type(gt), &
             nfftdim1_for_gridtype(gt), &
             collect_HC_grid, &
             HC_grid(off+1),out_lun)
      enddo !g = 1,num_HC_prim_grids
   endif !( tot_num_herm_Cprims > 0 )then

   ! next fill HD collective
   if ( tot_num_herm_Dprims > 0 )then
      do g = 1,num_HD_prim_grids
         off = HD_grid_offset_for_HD_grid(g)
         gt = grid_type_for_HD_prim_grid(g)
         call FH_RECIP_GRIDS_accum( &
             nfft1_for_HD_collect, &
             nfft2_for_HD_collect, &
             nfft3_for_HD_collect, &
             nfftdim1_for_HD_collect, &
             nfft1_for_grid_type(gt), &
             nfft2_for_grid_type(gt), &
             nfft3_for_grid_type(gt), &
             nfftdim1_for_gridtype(gt), &
             collect_HD_grid, &
             HD_grid(off+1),out_lun)
      enddo !g = 1,num_HD_prim_grids
   endif !( tot_num_herm_Dprims > 0 )then

end subroutine FH_RECIP_GRIDS_collective
!-----------------------------------------------------------------
subroutine FH_RECIP_GRIDS_mult_coulomb(out_lun)

   use hermite, only : tot_num_mpole_coeffs, &
                       tot_num_sumCH_coeffs, &
                       tot_num_herm_Cprims, &
                       tot_num_herm_Dprims

   implicit none

  integer, intent(in)   :: out_lun

   ! multiply each collective by fourier transform of 1/r (or truncated 1/r)

   ! first the multipoles
   if ( tot_num_mpole_coeffs > 0 )then
      call FH_RECIP_GRIDS_mult_one_coulomb( &
             nfft1_for_MPOLE_collect, &
             nfft2_for_MPOLE_collect, &
             nfft3_for_MPOLE_collect, &
             nfftdim1_for_MPOLE_collect, &
             nfft1_for_coulomb, &
             nfft2_for_coulomb, &
             nfft3_for_coulomb, &
             nfftdim1_for_coulomb, &
             nfft1_for_MPOLE_FT_pot, &
             nfft2_for_MPOLE_FT_pot, &
             nfft3_for_MPOLE_FT_pot, &
             nfftdim1_for_MPOLE_FT_pot, &
             collect_MC_grid,  &
             coulomb_grid, &
             MC_FT_pot_grid,out_lun)
            
      call FH_RECIP_GRIDS_mult_one_coulomb( &
             nfft1_for_MPOLE_collect, &
             nfft2_for_MPOLE_collect, &
             nfft3_for_MPOLE_collect, &
             nfftdim1_for_MPOLE_collect, &
             nfft1_for_coulomb, &
             nfft2_for_coulomb, &
             nfft3_for_coulomb, &
             nfftdim1_for_coulomb, &
             nfft1_for_MPOLE_FT_pot, &
             nfft2_for_MPOLE_FT_pot, &
             nfft3_for_MPOLE_FT_pot, &
             nfftdim1_for_MPOLE_FT_pot, &
             collect_MD_grid,  &
             coulomb_grid, &
             MD_FT_pot_grid,out_lun)
   endif !( tot_num_mpole_coeffs > 0 )then
            
   ! next the sumCH
   if ( tot_num_sumCH_coeffs > 0 )then
      call FH_RECIP_GRIDS_mult_one_coulomb( &
             nfft1_for_sumCH_collect, &
             nfft2_for_sumCH_collect, &
             nfft3_for_sumCH_collect, &
             nfftdim1_for_sumCH_collect, &
             nfft1_for_coulomb, &
             nfft2_for_coulomb, &
             nfft3_for_coulomb, &
             nfftdim1_for_coulomb, &
             nfft1_for_sumCH_FT_pot, &
             nfft2_for_sumCH_FT_pot, &
             nfft3_for_sumCH_FT_pot, &
             nfftdim1_for_sumCH_FT_pot, &
             collect_sumCH_C_grid,  &
             coulomb_grid, &
             sumCH_C_FT_pot_grid,out_lun)
            
      call FH_RECIP_GRIDS_mult_one_coulomb( &
             nfft1_for_sumCH_collect, &
             nfft2_for_sumCH_collect, &
             nfft3_for_sumCH_collect, &
             nfftdim1_for_sumCH_collect, &
             nfft1_for_coulomb, &
             nfft2_for_coulomb, &
             nfft3_for_coulomb, &
             nfftdim1_for_coulomb, &
             nfft1_for_sumCH_FT_pot, &
             nfft2_for_sumCH_FT_pot, &
             nfft3_for_sumCH_FT_pot, &
             nfftdim1_for_sumCH_FT_pot, &
             collect_sumCH_D_grid,  &
             coulomb_grid, &
             sumCH_D_FT_pot_grid,out_lun)
   endif !( tot_num_sumCH_coeffs > 0 )then
            
   ! next the compact hermite
   if ( tot_num_herm_Cprims > 0 )then
      call FH_RECIP_GRIDS_mult_one_coulomb( &
             nfft1_for_HC_collect, &
             nfft2_for_HC_collect, &
             nfft3_for_HC_collect, &
             nfftdim1_for_HC_collect, &
             nfft1_for_coulomb, &
             nfft2_for_coulomb, &
             nfft3_for_coulomb, &
             nfftdim1_for_coulomb, &
             nfft1_for_HC_FT_pot, &
             nfft2_for_HC_FT_pot, &
             nfft3_for_HC_FT_pot, &
             nfftdim1_for_HC_FT_pot, &
             collect_HC_grid, &
             coulomb_grid, &
             HC_FT_pot_grid,out_lun )
   endif !( tot_num_herm_Cprims > 0 )then
            
   ! finally the diffuse hermite
   if ( tot_num_herm_Dprims > 0 )then
      call FH_RECIP_GRIDS_mult_one_coulomb( &
             nfft1_for_HD_collect, &
             nfft2_for_HD_collect, &
             nfft3_for_HD_collect, &
             nfftdim1_for_HD_collect, &
             nfft1_for_coulomb, &
             nfft2_for_coulomb, &
             nfft3_for_coulomb, &
             nfftdim1_for_coulomb, &
             nfft1_for_HD_FT_pot, &
             nfft2_for_HD_FT_pot, &
             nfft3_for_HD_FT_pot, &
             nfftdim1_for_HD_FT_pot, &
             collect_HD_grid,  &
             coulomb_grid, &
             HD_FT_pot_grid,out_lun)
   endif !( tot_num_herm_Dprims > 0 )then

end subroutine FH_RECIP_GRIDS_mult_coulomb
!-----------------------------------------------------------
subroutine FH_RECIP_GRIDS_virial(out_lun)

! need two virial subroutines to be able to access the grids in correct order

   implicit none

   integer, intent(in)  :: out_lun

   ! mpoles first
   call FH_RECIP_GRIDS_virial_one(nfft1_for_MPOLE_collect,&
                                  nfft2_for_MPOLE_collect,&
                              nfft3_for_MPOLE_collect,&
                              nfftdim1_for_MPOLE_collect,collect_MC_grid,&
                              collect_MD_grid,1,out_lun)
   ! sumCH next
   call FH_RECIP_GRIDS_virial_one(nfft1_for_sumCH_collect,&
                                  nfft2_for_sumCH_collect,&
                                  nfft3_for_sumCH_collect,&
                                  nfftdim1_for_sumCH_collect,&
                                  collect_sumCH_C_grid,&
                                  collect_sumCH_D_grid,1,out_lun)
   ! next the hermites
   ! compacts first
   call FH_RECIP_GRIDS_virial_one(nfft1_for_HC_collect,nfft2_for_HC_collect,&
                              nfft3_for_HC_collect,nfftdim1_for_HC_collect,&
                              collect_HC_grid,collect_HC_grid,0,out_lun)
   ! diffuse next
   call FH_RECIP_GRIDS_virial_one(nfft1_for_HD_collect,nfft2_for_HD_collect,&
                              nfft3_for_HD_collect,nfftdim1_for_HD_collect,&
                              collect_HD_grid,collect_HD_grid,0,out_lun)

   return
end subroutine FH_RECIP_GRIDS_virial
!-----------------------------------------------------------
subroutine FH_RECIP_GRIDS_virial_one(nfft1,nfft2,nfft3,nfftdim1,&
                                 grid1,grid2,mp_grid,out_lun)
   ! This subroutine does a scalar sum to calculate the reciprocal virial
   ! similar to am_recip_scalar_sum

   ! GAC: this is similar structure as get_FT_coulomb, but changed to
   !      calculate virial from the SIX collective grids:
   !      collect_MC_grid, collect_MD_grid, collect_sumCH_C_grid,
   !      collect_sumCH_D_grid, collect_HC_grid and collect_HD_grid
   use unit_cell, only : recip

   implicit none
                   
   integer, intent(in) :: nfft1,nfft2,nfft3,nfftdim1,mp_grid
   double precision, intent(in) :: grid1(2,nfft3,nfftdim1,nfft2)
   double precision, intent(in) :: grid2(2,nfft3,nfftdim1,nfft2)
   integer, intent(in) :: out_lun

   integer :: k1,k2,k3,m1,m2,m3,k10,nf1,nf2,nf3

   double precision :: pi,mhat1,mhat2,mhat3,msq,sq_msq
   double precision :: struc2,eterm,vterm
   double precision :: mult_grid1, mult_grid2
   double precision :: tmp1,tmp2,mult,fac
   double precision :: vxx, vxy, vxz, vyy, vyz, vzz

   pi = 3.14159265358979323846d0
   !rec_virial(:,:)=0.d0
   ! for vterm I don't need fac: 
   ! fac = pi**2 / ewald_coeff**2 !DON'T NEED IT BECAUSE GRIDS MULTIPLIED BY
                                  !FACTORS in FH_PME_RECIP_grid_mult
                                  ! -> FH_PME_RECIP_grid_mult_orthog
   ! don't need eterm because all grids get multiplied by pme_multipliers
   ! in FH_RECIP_get_FT_density-> FH_PME_RECIP_FT_density ->
   !                              FH_PME_RECIP_fix_one_FT_density

!have to calculate virial for each collective grid separately
!note we only need 4 calculations because MPOLE and SUM_CH have same nffts

   !print *,nfft1,nfft2,nfft3,nfftdim1,mp_grid
   !print *,grid1(1,1,1,1),grid1(2,1,1,1),maxval(grid1)
   !print *,grid2(1,1,1,1),grid2(2,1,1,1),maxval(grid2)
   !print *,grid1(1,nfft3,nfftdim1,nfft2),grid1(2,nfft3,nfftdim1,nfft2)
   !print *,grid2(1,nfft3,nfftdim1,nfft2),grid2(2,nfft3,nfftdim1,nfft2)


   vxx = 0.d0
   vxy = 0.d0
   vxz = 0.d0
   vyy = 0.d0
   vyz = 0.d0
   vzz = 0.d0
   nf1 = nfft1/2
   if ( 2*nf1 < nfft1)nf1 = nf1+1
   nf2 = nfft2/2
   if ( 2*nf2 < nfft2)nf2 = nf2+1
   nf3 = nfft3/2
   if ( 2*nf3 < nfft3)nf3 = nf3+1

! first the case for mhat = (0,0,0)

   do k2 = 1, nfft2
      m2 = k2 - 1
      if ( k2 > nf2 )m2 = k2 - 1 - nfft2
      do k3 = 1,nfft3
         m3 = k3 - 1
         if ( k3 > nf3 )m3 = k3 - 1 - nfft3
         k10 = 1
         ! already did (1,1,1)
         if(k3+k2 == 2) k10 = 2
         do k1 = k10, nf1+1
            m1 = k1 - 1
            if ( k1 > nf1 )m1 = k1 - 1 - nfft1
            mhat1 = recip(1,1)*m1+recip(1,2)*m2+recip(1,3)*m3
            mhat2 = recip(2,1)*m1+recip(2,2)*m2+recip(2,3)*m3
            mhat3 = recip(3,1)*m1+recip(3,2)*m2+recip(3,3)*m3
            msq = mhat1*mhat1+mhat2*mhat2+mhat3*mhat3
            vterm = 2.d0 * (msq + 1.d0) / msq
            ! if mults or sum_ch then we have two grids for each
            ! (compact and diffuse in grid1 and grid2 respectively)
            if (mp_grid == 1) then
               mult_grid1 = (grid1(1,k3,k1,k2)+&
                            grid2(1,k3,k1,k2))*&
                            (grid1(1,k3,k1,k2)+&
                            grid2(1,k3,k1,k2))
               mult_grid2 = (grid1(2,k3,k1,k2)+&
                            grid2(2,k3,k1,k2))*&
                            (grid1(2,k3,k1,k2)+&
                            grid2(2,k3,k1,k2))
            ! if hermites then we only have one grid
            else if (mp_grid == 0) then
               mult_grid1 = grid1(1,k3,k1,k2)*grid1(1,k3,k1,k2)
               mult_grid2 = grid1(2,k3,k1,k2)*grid1(2,k3,k1,k2)
            else
               write(out_lun,*)'Error in input for virial calculation!!'
            endif
!            if (mp_grid == 1) then
!               mult_grid1 = grid1(1,k3,k1,k2)*&
!                            grid1(1,k3,k1,k2)
!               mult_grid2 = grid2(2,k3,k1,k2)*&
!                            grid2(2,k3,k1,k2)
!            else if (mp_grid == 0) then
!               mult_grid1 = grid1(1,k3,k1,k2)
!               mult_grid2 = grid1(2,k3,k1,k2)
!            else
!               write(out_lun,*)'Error in input for virial calculation!!'
!            endif

            struc2 = mult_grid1+mult_grid2
            tmp1 = struc2 ! removed eterm here because grids multiplied by
                          ! pme_multipliers but may need to fix this
            tmp2 = tmp1 * vterm
            vxx = vxx + tmp2 * mhat1 * mhat1 - tmp1
            vxy = vxy + tmp2 * mhat1 * mhat2
            vxz = vxz + tmp2 * mhat1 * mhat3
            vyy = vyy + tmp2 * mhat2 * mhat2 - tmp1
            vyz = vyz + tmp2 * mhat2 * mhat3
            vzz = vzz + tmp2 * mhat3 * mhat3 - tmp1
         enddo
      enddo
   enddo
   rec_virial(1, 1) = rec_virial(1, 1) + 0.5d0 * vxx
   rec_virial(1, 2) = rec_virial(1, 2) + 0.5d0 * vxy
   rec_virial(2, 1) = rec_virial(2, 1) + 0.5d0 * vxy
   rec_virial(1, 3) = rec_virial(1, 3) + 0.5d0 * vxz
   rec_virial(3, 1) = rec_virial(3, 1) + 0.5d0 * vxz
   rec_virial(2, 2) = rec_virial(2, 2) + 0.5d0 * vyy
   rec_virial(2, 3) = rec_virial(2, 3) + 0.5d0 * vyz
   rec_virial(3, 2) = rec_virial(3, 2) + 0.5d0 * vyz
   rec_virial(3, 3) = rec_virial(3, 3) + 0.5d0 * vzz

   !print *,rec_virial(1, 1),rec_virial(1, 2),rec_virial(1, 3) 
   !print *,rec_virial(2, 1),rec_virial(2, 2),rec_virial(2, 3) 
   !print *,rec_virial(3, 1),rec_virial(3, 2),rec_virial(3, 3) 

   return

end subroutine FH_RECIP_GRIDS_virial_one
!-----------------------------------------------------------------
subroutine FH_RECIP_GRIDS_FT_ene( &
                                 recip_mp_mp_ene, &
                                 recip_mp_ch_ene, &
                                 recip_ch_mp_ene, &
                                 recip_ch_ch_ene, &
                                 recip_mp_dh_ene, &
                                 recip_dh_mp_ene, &
                                 recip_ch_dh_ene, &
                                 recip_dh_ch_ene, &
                                 recip_dh_dh_ene, &
                                 out_lun)

   use hermite, only : tot_num_mpole_coeffs, &
                       tot_num_sumCH_coeffs, &
                       tot_num_herm_Cprims, &
                       tot_num_herm_Dprims
   use unit_cell, only : volume
   use user, only : diffuse_adjust_exponent

   implicit none

   include 'scale.fh'

   double precision, intent(out) :: recip_mp_mp_ene, &
                                    recip_mp_ch_ene, &
                                    recip_ch_mp_ene, &
                                    recip_ch_ch_ene, &
                                    recip_mp_dh_ene, &
                                    recip_dh_mp_ene, &
                                    recip_ch_dh_ene, &
                                    recip_dh_ch_ene, &
                                    recip_dh_dh_ene
   integer, intent(in)           :: out_lun
   double precision energy

   ! mpole mpole energy
   ! just the MC_MD energy this is the ewald compensating ene
   if ( tot_num_mpole_coeffs > 0 )then
      call FH_RECIP_GRIDS_2grid_ene( &
                nfft1_for_MPOLE_collect, &
                nfft2_for_MPOLE_collect, &
                nfft3_for_MPOLE_collect, &
                nfftdim1_for_MPOLE_collect, &
                nfft1_for_MPOLE_FT_pot, &
                nfft2_for_MPOLE_FT_pot, &
                nfft3_for_MPOLE_FT_pot, &
                nfftdim1_for_MPOLE_FT_pot, &
                collect_MC_grid, &
                MD_FT_pot_grid, &
                volume, &
                energy,out_lun) 
       recip_mp_mp_ene = 0.5d0*energy
   endif !( tot_num_mpole_coeffs > 0 )then

   ! sumCH sumCH energy
   ! just the sumCH_C_sumCH_D energy this is the ewald compensating ene
   if ( tot_num_sumCH_coeffs > 0 )then
      call FH_RECIP_GRIDS_2grid_ene( &
                nfft1_for_sumCH_collect, &
                nfft2_for_sumCH_collect, &
                nfft3_for_sumCH_collect, &
                nfftdim1_for_sumCH_collect, &
                nfft1_for_sumCH_FT_pot, &
                nfft2_for_sumCH_FT_pot, &
                nfft3_for_sumCH_FT_pot, &
                nfftdim1_for_sumCH_FT_pot, &
                collect_sumCH_C_grid, &
                sumCH_D_FT_pot_grid, &
                volume, &
                energy,out_lun) 
       recip_ch_ch_ene = 0.5d0*energy
   endif !( tot_num_sumCH_coeffs > 0 )then

   ! sumCH mpole energy
   if ( (tot_num_mpole_coeffs > 0) .and. (tot_num_sumCH_coeffs > 0) )then
      call FH_RECIP_GRIDS_2grid_ene( &
                nfft1_for_mpole_collect, &
                nfft2_for_mpole_collect, &
                nfft3_for_mpole_collect, &
                nfftdim1_for_mpole_collect, &
                nfft1_for_sumCH_FT_pot, &
                nfft2_for_sumCH_FT_pot, &
                nfft3_for_sumCH_FT_pot, &
                nfftdim1_for_sumCH_FT_pot, &
                collect_MC_grid, &
                sumCH_D_FT_pot_grid, &
                volume, &
                energy,out_lun) 
      recip_mp_ch_ene = -0.5d0*energy ! negative for nuclear-electron
      call FH_RECIP_GRIDS_2grid_ene( &
                nfft1_for_sumCH_collect, &
                nfft2_for_sumCH_collect, &
                nfft3_for_sumCH_collect, &
                nfftdim1_for_sumCH_collect, &
                nfft1_for_mpole_FT_pot, &
                nfft2_for_mpole_FT_pot, &
                nfft3_for_mpole_FT_pot, &
                nfftdim1_for_mpole_FT_pot, &
                collect_sumCH_C_grid, &
                MD_FT_pot_grid, &
                volume,  &
                energy,out_lun) 
      recip_ch_mp_ene = -0.5d0*energy ! negative for nuclear-electron
   endif !( (tot_num_mpole_coeffs > 0) .and. (tot_num_sumCH_coeffs > 0) 

   ! mpole diffuse hermite energy
   if ( (tot_num_mpole_coeffs > 0) .and. (tot_num_herm_Dprims > 0) )then
      call FH_RECIP_GRIDS_2grid_ene( &
                nfft1_for_HD_FT_pot, &
                nfft2_for_HD_FT_pot, &
                nfft3_for_HD_FT_pot, &
                nfftdim1_for_HD_FT_pot, &
                nfft1_for_MPOLE_collect, &
                nfft2_for_MPOLE_collect, &
                nfft3_for_MPOLE_collect, &
                nfftdim1_for_MPOLE_collect, &
                collect_HD_grid, &
                MC_FT_pot_grid, &
                volume,  &
                energy,out_lun) 
      recip_dh_mp_ene = -0.5d0*energy ! negative for nuclear-electron
      call FH_RECIP_GRIDS_2grid_ene( &
                nfft1_for_MPOLE_collect, &
                nfft2_for_MPOLE_collect, &
                nfft3_for_MPOLE_collect, &
                nfftdim1_for_MPOLE_collect, &
                nfft1_for_HD_FT_pot, &
                nfft2_for_HD_FT_pot, &
                nfft3_for_HD_FT_pot, &
                nfftdim1_for_HD_FT_pot, &
                collect_MC_grid, &
                HD_FT_pot_grid, &
                volume,  &
                energy,out_lun) 
      recip_mp_dh_ene = -0.5d0*energy ! negative for nuclear-electron
   endif !( (tot_num_mpole_coeffs > 0) .and. (tot_num_herm_Dprims > 0) )then

   ! compact hermite - diffuse hermite energy
   if ( (tot_num_herm_Cprims > 0) .and. (tot_num_herm_Dprims > 0) )then
      ! first the HC_HD energy
      call FH_RECIP_GRIDS_2grid_ene( &
                nfft1_for_HD_collect, &
                nfft2_for_HD_collect, &
                nfft3_for_HD_collect, &
                nfftdim1_for_HD_collect, &
                nfft1_for_HC_FT_pot, &
                nfft2_for_HC_FT_pot, &
                nfft3_for_HC_FT_pot, &
                nfftdim1_for_HC_FT_pot, &
                collect_HD_grid, &
                HC_FT_pot_grid, &
                volume,  &
                energy,out_lun) 
      recip_ch_dh_ene = 0.5d0*energy ! positive for electron-electron
      ! next the HD_HC energy
      call FH_RECIP_GRIDS_2grid_ene( &
                nfft1_for_HC_collect, &
                nfft2_for_HC_collect, &
                nfft3_for_HC_collect, &
                nfftdim1_for_HC_collect, &
                nfft1_for_HD_FT_pot, &
                nfft2_for_HD_FT_pot, &
                nfft3_for_HD_FT_pot, &
                nfftdim1_for_HD_FT_pot, &
                collect_HC_grid, &
                HD_FT_pot_grid, &
                volume,  &
                energy,out_lun) 
      recip_dh_ch_ene = 0.5d0*energy ! positive for electron-electron
   endif !( (tot_num_herm_Cprims > 0) .and. (tot_num_herm_Dprims > 0) )then

   ! finally the diffuse-diffuse hermite energy
   if ( tot_num_herm_Dprims > 0 )then
      if ( diffuse_adjust_exponent < small_value )then
         call FH_RECIP_GRIDS_2grid_ene( &
                nfft1_for_HD_collect, &
                nfft2_for_HD_collect, &
                nfft3_for_HD_collect, &
                nfftdim1_for_HD_collect, &
                nfft1_for_HD_FT_pot, &
                nfft2_for_HD_FT_pot, &
                nfft3_for_HD_FT_pot, &
                nfftdim1_for_HD_FT_pot, &
                collect_HD_grid, &
                HD_FT_pot_grid, &
                volume, &
                energy,out_lun) 
      else
         call FH_RECIP_LOC_diffadj_2grid_ene( &
                nfft1_for_HD_collect, &
                nfft2_for_HD_collect, &
                nfft3_for_HD_collect, &
                nfftdim1_for_HD_collect, &
                nfft1_for_HD_FT_pot, &
                nfft2_for_HD_FT_pot, &
                nfft3_for_HD_FT_pot, &
                nfftdim1_for_HD_FT_pot, &
                nfft1_for_diff_adj, &
                nfft2_for_diff_adj, &
                nfft3_for_diff_adj, &
                nfftdim1_for_diff_adj, &
                collect_HD_grid, &
                HD_FT_pot_grid, &
                diff_adj_grid, &
                volume,  &
                energy,out_lun)
      endif !( diffuse_adjust_exponent < small_value )then
      recip_dh_dh_ene = 0.5d0*energy ! positive for electron-electron
   endif !( tot_num_herm_Dprims > 0 )then

end subroutine FH_RECIP_GRIDS_FT_ene
!-----------------------------------------------------------------
subroutine RECIP_GRIDS_load_FT_pot(out_lun)

   use user, only : nfft1_for_grid_type, &
                    nfft2_for_grid_type, &
                    nfft3_for_grid_type, &
                    grid_type_for_MPOLES, &
                    grid_type_for_sumCH, &
                    diffuse_adjust_exponent, &
                    num_HC_prim_grids, &
                    grid_type_for_HC_prim_grid, &
                    num_HD_prim_grids, &
                    grid_type_for_HD_prim_grid
   use unit_cell, only : volume
   use hermite, only : tot_num_mpole_coeffs, &
                       tot_num_sumCH_coeffs, &
                       tot_num_herm_Cprims, &
                       tot_num_herm_Dprims

   implicit none

   integer, intent(in)  :: out_lun

   include 'scale.fh'

   integer g,gt,off,signum

   ! first clear the coefficient arrays
   if ( tot_num_mpole_coeffs > 0 )then
      call UTIL_zero_real_array(MC_grid,tot_siz_MPOLE_grid)
   endif
   if ( tot_num_sumCH_coeffs > 0 )then
      call UTIL_zero_real_array(sumCH_C_grid,tot_siz_sumCH_grid)
   endif
   if ( tot_num_herm_Cprims > 0 )then
      call UTIL_zero_real_array(HC_grid,tot_siz_HC_grid)
   endif
   if ( tot_num_herm_Dprims > 0 )then
      call UTIL_zero_real_array(HD_grid,tot_siz_HD_grid)
   endif

   ! now accumulate the potentials
   ! first the multipole grid to get forces  on multipoles
   if ( tot_num_mpole_coeffs > 0 )then
      gt = grid_type_for_MPOLES
      ! add the diffuse mpole potential for MD_MC interactions
      signum = 1  ! mpole mpole are positive
      call FH_RECIP_GRIDS_add_to_FT_pot( &
                                 nfft1_for_MPOLE_FT_pot, &
                                 nfft2_for_MPOLE_FT_pot, &
                                 nfft3_for_MPOLE_FT_pot, &
                                 nfftdim1_for_MPOLE_FT_pot, &
                                 nfft1_for_grid_type(gt), &
                                 nfft2_for_grid_type(gt), &
                                 nfft3_for_grid_type(gt), &
                                 nfftdim1_for_gridtype(gt), &
                                 volume, &
                                 MD_FT_pot_grid, &
                                 MC_grid, &
                                 signum,out_lun)
      ! add the diffuse sumCH potential for sumCH_D_MC interactions
      if ( tot_num_sumCH_coeffs > 0 )then
         signum = -1  ! mpole sumCH are negative e.g. nuclear-electron
         call FH_RECIP_GRIDS_add_to_FT_pot( &
                                 nfft1_for_sumCH_FT_pot, &
                                 nfft2_for_sumCH_FT_pot, &
                                 nfft3_for_sumCH_FT_pot, &
                                 nfftdim1_for_sumCH_FT_pot, &
                                 nfft1_for_grid_type(gt), &
                                 nfft2_for_grid_type(gt), &
                                 nfft3_for_grid_type(gt), &
                                 nfftdim1_for_gridtype(gt), &
                                 volume, &
                                 sumCH_D_FT_pot_grid, &
                                 MC_grid, &
                                 signum,out_lun)
      endif !( tot_num_sumCH_coeffs > 0 )then
      ! finally add the diffuse hermite potential for HD_MC
      if ( tot_num_herm_Dprims > 0 )then
         signum = -1  ! mpole sumCH are negative e.g. nuclear-electron
         call FH_RECIP_GRIDS_add_to_FT_pot( &
                                 nfft1_for_HD_FT_pot, &
                                 nfft2_for_HD_FT_pot, &
                                 nfft3_for_HD_FT_pot, &
                                 nfftdim1_for_HD_FT_pot, &
                                 nfft1_for_grid_type(gt), &
                                 nfft2_for_grid_type(gt), &
                                 nfft3_for_grid_type(gt), &
                                 nfftdim1_for_gridtype(gt), &
                                 volume, &
                                 HD_FT_pot_grid, &
                                 MC_grid, &
                                 signum,out_lun)
      endif !( tot_num_herm_Dprims > 0 )then
   endif !( tot_num_mpole_coeffs > 0 )then

   ! next do the forces on sumCH
   if ( tot_num_sumCH_coeffs > 0 )then
      gt = grid_type_for_sumCH
      ! add the diffuse mpole potential for MD_sumCH interactions
      if ( tot_num_mpole_coeffs > 0 )then
         signum = -1  ! mpole sumCH are negative e.g. nuclear-electron
         call FH_RECIP_GRIDS_add_to_FT_pot( &
                                 nfft1_for_MPOLE_FT_pot, &
                                 nfft2_for_MPOLE_FT_pot, &
                                 nfft3_for_MPOLE_FT_pot, &
                                 nfftdim1_for_MPOLE_FT_pot, &
                                 nfft1_for_grid_type(gt), &
                                 nfft2_for_grid_type(gt), &
                                 nfft3_for_grid_type(gt), &
                                 nfftdim1_for_gridtype(gt), &
                                 volume, &
                                 MD_FT_pot_grid, &
                                 sumCH_C_grid, &
                                 signum,out_lun)
      endif !( tot_num_mpole_coeffs > 0 )then
      ! add the diffuse sumCH potential for sumCH_D_sumCH_C interactions
      signum = 1  ! sumCH sumCH are positive e.g. electron-electron
      call FH_RECIP_GRIDS_add_to_FT_pot( &
                                 nfft1_for_sumCH_FT_pot, &
                                 nfft2_for_sumCH_FT_pot, &
                                 nfft3_for_sumCH_FT_pot, &
                                 nfftdim1_for_sumCH_FT_pot, &
                                 nfft1_for_grid_type(gt), &
                                 nfft2_for_grid_type(gt), &
                                 nfft3_for_grid_type(gt), &
                                 nfftdim1_for_gridtype(gt), &
                                 volume, &
                                 sumCH_D_FT_pot_grid, &
                                 sumCH_C_grid, &
                                 signum,out_lun)
   endif !( tot_num_sumCH_coeffs > 0 )then

   ! next the diffuse hermite acting on compact hermites
   ! the effect of mpoles and compact hermites themselves on compact hermites
   ! already accounted for in sumCH_C_grid
   if ( (tot_num_herm_Cprims > 0) .and. (tot_num_herm_Dprims > 0) )then
      signum = 1  ! hermite hermite are positive
      do g = 1,num_HC_prim_grids
         off = HC_grid_offset_for_HC_grid(g)
         gt = grid_type_for_HC_prim_grid(g)
         ! add the diffuse hermite potential for HD_HC interactions
         call FH_RECIP_GRIDS_add_to_FT_pot( &
                                 nfft1_for_HD_FT_pot, &
                                 nfft2_for_HD_FT_pot, &
                                 nfft3_for_HD_FT_pot, &
                                 nfftdim1_for_HD_FT_pot, &
                                 nfft1_for_grid_type(gt), &
                                 nfft2_for_grid_type(gt), &
                                 nfft3_for_grid_type(gt), &
                                 nfftdim1_for_gridtype(gt), &
                                 volume, &
                                 HD_FT_pot_grid, &
                                 HC_grid(off+1), &
                                 signum,out_lun)
      enddo !g = 1,num_HC_prim_grids
   endif !( (tot_num_herm_Cprims > 0) .and. (tot_num_herm_Dprims > 0) )then

   ! finally accumulate onto HD_grid
   if ( tot_num_herm_Dprims > 0 )then
      do g = 1,num_HD_prim_grids
         off = HD_grid_offset_for_HD_grid(g)
         gt = grid_type_for_HC_prim_grid(g)
         ! add the compact MPOLE potential, for MC_HD interactions
         signum = -1  ! mpole hermite are negative
         call FH_RECIP_GRIDS_add_to_FT_pot( &
                                 nfft1_for_MPOLE_FT_pot, &
                                 nfft2_for_MPOLE_FT_pot, &
                                 nfft3_for_MPOLE_FT_pot, &
                                 nfftdim1_for_MPOLE_FT_pot, &
                                 nfft1_for_grid_type(gt), &
                                 nfft2_for_grid_type(gt), &
                                 nfft3_for_grid_type(gt), &
                                 nfftdim1_for_gridtype(gt), &
                                 volume, &
                                 MC_FT_pot_grid, &
                                 HD_grid(off+1), &
                                 signum,out_lun)
         ! add the compact hermite potential, for HC_HD interactions
         signum = 1  ! hermite hermite are positive
         call FH_RECIP_GRIDS_add_to_FT_pot( &
                                 nfft1_for_HC_FT_pot, &
                                 nfft2_for_HC_FT_pot, &
                                 nfft3_for_HC_FT_pot, &
                                 nfftdim1_for_HC_FT_pot, &
                                 nfft1_for_grid_type(gt), &
                                 nfft2_for_grid_type(gt), &
                                 nfft3_for_grid_type(gt), &
                                 nfftdim1_for_gridtype(gt), &
                                 volume, &
                                 HC_FT_pot_grid, &
                                 HD_grid(off+1), &
                                 signum,out_lun)
         ! add the diffuse hermite potential, for HD_HD interactions
         signum = 1  ! hermite hermite are positive
         if ( diffuse_adjust_exponent < small_value )then
            call FH_RECIP_GRIDS_add_to_FT_pot( &
                                 nfft1_for_HD_FT_pot, &
                                 nfft2_for_HD_FT_pot, &
                                 nfft3_for_HD_FT_pot, &
                                 nfftdim1_for_HD_FT_pot, &
                                 nfft1_for_grid_type(gt), &
                                 nfft2_for_grid_type(gt), &
                                 nfft3_for_grid_type(gt), &
                                 nfftdim1_for_gridtype(gt), &
                                 volume, &
                                 HD_FT_pot_grid, &
                                 HD_grid(off+1), &
                                 signum,out_lun)
         else
            ! this one hard-wired to be positive
            call FH_RECIP_GRIDS_da_add_to_FT_pot(  &
                                 nfft1_for_HD_FT_pot, &
                                 nfft2_for_HD_FT_pot, &
                                 nfft3_for_HD_FT_pot, &
                                 nfftdim1_for_HD_FT_pot, &
                                 nfft1_for_grid_type(gt), &
                                 nfft2_for_grid_type(gt), &
                                 nfft3_for_grid_type(gt), &
                                 nfftdim1_for_gridtype(gt), &
                                 nfft1_for_diff_adj, &
                                 nfft2_for_diff_adj, &
                                 nfft3_for_diff_adj, &
                                 nfftdim1_for_diff_adj, &
                                 volume, &
                                 HD_FT_pot_grid, &
                                 HD_grid(off+1), &
                                 diff_adj_grid, out_lun)
         endif !( diffuse_adjust_exponent < small_value )then
      enddo !g = 1,num_HD_prim_grids
   endif !( tot_num_herm_Dprims > 0 )then

end subroutine RECIP_GRIDS_load_FT_pot
!-----------------------------------------------------------------
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!-----------------------------------------------------------------
! DEALLOCATION SUBROUTINES
!-----------------------------------------------------------------
!-----------------------------------------------------------------

!-----------------------------------------------------------------
subroutine FH_RECIP_GRIDS_deallocate()

   implicit none

   call FH_RECIP_GRIDS_deallocate_fft()
   call FH_RECIP_GRIDS_deallocate_ptrs()
   call FH_RECIP_GRIDS_deall_CD_grids()
   call FH_RECIP_GRIDS_deall_collect()
   call FH_RECIP_GRIDS_deall_coulomb()
   call FH_RECIP_GRIDS_deall_diff_adj()

end subroutine FH_RECIP_GRIDS_deallocate
!-----------------------------------------------------------------
subroutine FH_RECIP_GRIDS_deallocate_fft()

   implicit none

   if ( allocated(nfftdim1_for_gridtype) )deallocate(nfftdim1_for_gridtype)
   if ( allocated(nfftdim2_for_gridtype) )deallocate(nfftdim2_for_gridtype)
   if ( allocated(nfftdim3_for_gridtype) )deallocate(nfftdim3_for_gridtype)
   if ( allocated(nfftable_for_gridtype) )deallocate(nfftable_for_gridtype)
   if ( allocated(nffwork_for_gridtype) )deallocate(nffwork_for_gridtype)
   if ( allocated(sizfftab_for_gridtype) )deallocate(sizfftab_for_gridtype)
   if ( allocated(off_fftab_for_gridtype) )deallocate(off_fftab_for_gridtype)
   if ( allocated(sizffwrk_for_gridtype) )deallocate(sizffwrk_for_gridtype)
   if ( allocated(fftable) )deallocate(fftable)
   if ( allocated(ffwork) )deallocate(ffwork)
   if ( allocated(fft_tmpy) )deallocate(fft_tmpy)
   if ( allocated(fft_alpha) )deallocate(fft_alpha)
   if ( allocated(fft_beta) )deallocate(fft_beta)

end subroutine FH_RECIP_GRIDS_deallocate_fft
!-----------------------------------------------------------------
subroutine FH_RECIP_GRIDS_deallocate_ptrs()

   implicit none

   ! mpole stuff
   if ( allocated(Mpole_order_for_gridded_mpole) ) &
          deallocate(Mpole_order_for_gridded_mpole)
   if ( allocated(Fcoeff_offset_for_gridded_mpole) ) &
          deallocate(Fcoeff_offset_for_gridded_mpole)
   if ( allocated(Gcoeff_offset_for_gridded_mpole) ) &
          deallocate(Gcoeff_offset_for_gridded_mpole)
   if ( allocated(Ffield_order_for_gridded_mpole) ) &
          deallocate(Ffield_order_for_gridded_mpole)
   if ( allocated(Ffield_offset_for_gridded_mpole) ) &
          deallocate(Ffield_offset_for_gridded_mpole)
   if ( allocated(site_that_owns_gridded_mpole) ) &
          deallocate(site_that_owns_gridded_mpole)
   if ( allocated(Mpole_grid_Fcoeff) )deallocate(Mpole_grid_Fcoeff)
   if ( allocated(Mpole_grid_Ffield) )deallocate(Mpole_grid_Ffield)

   ! sumCH stuff
   if ( allocated(sumCH_order_for_gridded_sumCH) ) &
          deallocate(sumCH_order_for_gridded_sumCH)
   if ( allocated(Fcoeff_offset_for_gridded_sumCH) ) &
          deallocate(Fcoeff_offset_for_gridded_sumCH)
   if ( allocated(Gcoeff_offset_for_gridded_sumCH) ) &
          deallocate(Gcoeff_offset_for_gridded_sumCH)
   if ( allocated(Ffield_order_for_gridded_sumCH) ) &
          deallocate(Ffield_order_for_gridded_sumCH)
   if ( allocated(Ffield_offset_for_gridded_sumCH) ) &
          deallocate(Ffield_offset_for_gridded_sumCH)
   if ( allocated(site_that_owns_gridded_sumCH) ) &
          deallocate(site_that_owns_gridded_sumCH)
   if ( allocated(sumCH_grid_Fcoeff) )deallocate(sumCH_grid_Fcoeff)
   if ( allocated(sumCH_grid_Ffield) )deallocate(sumCH_grid_Ffield)

   ! compact hermite stuff
   if ( allocated(num_prims_for_HC_grid) )deallocate(num_prims_for_HC_grid)
   if ( allocated(prim_offset_for_HC_grid) ) &
             deallocate(prim_offset_for_HC_grid)
   if ( allocated(prim_num_for_HC_grid_prim) ) &
              deallocate(prim_num_for_HC_grid_prim)
   if ( allocated(herm_order_for_HC_grid_prim) ) &
          deallocate(herm_order_for_HC_grid_prim)
   if ( allocated(Ffield_order_for_HC_grid_prim) ) &
          deallocate(Ffield_order_for_HC_grid_prim)
   if ( allocated(herm_expon_for_HC_grid_prim) ) &
          deallocate(herm_expon_for_HC_grid_prim)
   if ( allocated(Fcoeff_offset_for_HC_grid_prim) ) &
          deallocate(Fcoeff_offset_for_HC_grid_prim)
   if ( allocated(Ffield_offset_for_HC_grid_prim) ) &
          deallocate(Ffield_offset_for_HC_grid_prim)
   if ( allocated(Gcoeff_offset_for_HC_grid_prim) ) &
          deallocate(Gcoeff_offset_for_HC_grid_prim)
   if ( allocated(site_that_owns_HC_grid_prim) ) &
          deallocate(site_that_owns_HC_grid_prim)
   if ( allocated(HC_grid_Fcoeff) )deallocate(HC_grid_Fcoeff)
   if ( allocated(HC_grid_Ffield) )deallocate(HC_grid_Ffield)

   ! diffuse hermite stuff
   if ( allocated(num_prims_for_HD_grid) )deallocate(num_prims_for_HD_grid)
   if ( allocated(prim_offset_for_HD_grid) ) &
             deallocate(prim_offset_for_HD_grid)
   if ( allocated(prim_num_for_HD_grid_prim) ) &
              deallocate(prim_num_for_HD_grid_prim)
   if ( allocated(herm_order_for_HD_grid_prim) ) &
          deallocate(herm_order_for_HD_grid_prim)
   if ( allocated(Ffield_order_for_HD_grid_prim) ) &
          deallocate(Ffield_order_for_HD_grid_prim)
   if ( allocated(herm_expon_for_HD_grid_prim) ) &
          deallocate(herm_expon_for_HD_grid_prim)
   if ( allocated(Fcoeff_offset_for_HD_grid_prim) ) &
          deallocate(Fcoeff_offset_for_HD_grid_prim)
   if ( allocated(Ffield_offset_for_HD_grid_prim) ) &
          deallocate(Ffield_offset_for_HD_grid_prim)
   if ( allocated(Gcoeff_offset_for_HD_grid_prim) ) &
          deallocate(Gcoeff_offset_for_HD_grid_prim)
   if ( allocated(site_that_owns_HD_grid_prim) ) &
          deallocate(site_that_owns_HD_grid_prim)
   if ( allocated(HD_grid_Fcoeff) )deallocate(HD_grid_Fcoeff)
   if ( allocated(HD_grid_Ffield) )deallocate(HD_grid_Ffield)

end subroutine FH_RECIP_GRIDS_deallocate_ptrs
!-----------------------------------------------------------------
subroutine FH_RECIP_GRIDS_deall_CD_grids()

   implicit none

   ! mpoles
   if ( allocated(MC_grid) )deallocate(MC_grid)
   if ( allocated(MD_grid) )deallocate(MD_grid)
   ! sumCH
   if ( allocated(sumCH_C_grid) )deallocate(sumCH_C_grid)
   if ( allocated(sumCH_D_grid) )deallocate(sumCH_D_grid)
   ! compact hermites
   if ( allocated(HC_grid_offset_for_HC_grid) ) &
           deallocate(HC_grid_offset_for_HC_grid)
   if ( allocated(HC_grid) )deallocate(HC_grid)
   ! diffuse hermites
   if ( allocated(HD_grid_offset_for_HD_grid) ) &
           deallocate(HD_grid_offset_for_HD_grid)
   if ( allocated(HD_grid) )deallocate(HD_grid)

end subroutine FH_RECIP_GRIDS_deall_CD_grids
!-----------------------------------------------------------------
subroutine FH_RECIP_GRIDS_deall_collect()

   implicit none

   ! mpoles
   if ( allocated(collect_MC_grid) )deallocate(collect_MC_grid)
   if ( allocated(collect_MD_grid) )deallocate(collect_MD_grid)
   if ( allocated(MC_FT_pot_grid) )deallocate(MC_FT_pot_grid)
   if ( allocated(MD_FT_pot_grid) )deallocate(MD_FT_pot_grid)
   ! sumCH
   if ( allocated(collect_sumCH_C_grid) )deallocate(collect_sumCH_C_grid)
   if ( allocated(collect_sumCH_D_grid) )deallocate(collect_sumCH_D_grid)
   if ( allocated(sumCH_C_FT_pot_grid) )deallocate(sumCH_C_FT_pot_grid)
   if ( allocated(sumCH_D_FT_pot_grid) )deallocate(sumCH_D_FT_pot_grid)
   ! compact hermites
   if ( allocated(collect_HC_grid) )deallocate(collect_HC_grid)
   if ( allocated(HC_FT_pot_grid) )deallocate(HC_FT_pot_grid)
   ! diffuse hermites
   if ( allocated(collect_HD_grid) )deallocate(collect_HD_grid)
   if ( allocated(HD_FT_pot_grid) )deallocate(HD_FT_pot_grid)

end subroutine FH_RECIP_GRIDS_deall_collect
!-----------------------------------------------------------------
subroutine FH_RECIP_GRIDS_deall_coulomb()

   implicit none

   if ( allocated(coulomb_grid) )deallocate(coulomb_grid)
end subroutine FH_RECIP_GRIDS_deall_coulomb
!-----------------------------------------------------------------
subroutine FH_RECIP_GRIDS_deall_diff_adj()

   implicit none

   if ( allocated(diff_adj_grid) )deallocate(diff_adj_grid)
end subroutine FH_RECIP_GRIDS_deall_diff_adj
!-----------------------------------------------------------------
end module recip_grids
!-----------------------------------------------------------------
!-----------------------------------------------------------------
! NON-MODULE (LOCAL) SUBROUTINES
! THESE HANDLE ALL ARGUMENTS THROUGH ARGUMENT LIST
!-----------------------------------------------------------------
!-----------------------------------------------------------------
subroutine FH_RECIP_GRIDS_get_FT_coulomb( &
                           nfft1,nfft2,nfft3,nfftdim1,PBC, &
                           sphere,recip,FT_coulomb)

   implicit none

   integer, intent(in) :: nfft1,nfft2,nfft3,nfftdim1
   integer, intent(in) :: PBC
   double precision, intent(in) :: sphere,recip(3,3)
   double precision, intent(out) :: FT_coulomb(nfft3,nfftdim1,nfft2)
                   
   integer :: k1,k2,k3,m1,m2,m3,k10,nf1,nf2,nf3
   double precision :: pi,mhat1,mhat2,mhat3,msq,sq_msq

   pi = 3.14159265358979323846d0

   nf1 = nfft1/2
   if ( 2*nf1 < nfft1 )nf1 = nf1+1
   nf2 = nfft2/2
   if ( 2*nf2 < nfft2 )nf2 = nf2+1
   nf3 = nfft3/2
   if ( 2*nf3 < nfft3 )nf3 = nf3+1

! first the case for mhat = (0,0,0)
   if ( PBC == 1 )then
      FT_Coulomb(1,1,1) = 0.d0 !actually this is the singularity---
                               !need neutral cell
   else
      FT_Coulomb(1,1,1) = 2.d0*pi*sphere**2
   endif

   do k2 = 1, nfft2
      m2 = k2 - 1
      if ( k2 > nf2 )m2 = k2 - 1 - nfft2
      do k3 = 1,nfft3
         m3 = k3 - 1
         if ( k3 > nf3 )m3 = k3 - 1 - nfft3
         k10 = 1
         ! already did (1,1,1)
         if(k3+k2 == 2) k10 = 2
         do k1 = k10, nf1+1
            m1 = k1 - 1
            if ( k1 > nf1 )m1 = k1 - 1 - nfft1
            mhat1 = recip(1,1)*m1+recip(1,2)*m2+recip(1,3)*m3
            mhat2 = recip(2,1)*m1+recip(2,2)*m2+recip(2,3)*m3
            mhat3 = recip(3,1)*m1+recip(3,2)*m2+recip(3,3)*m3
            msq = mhat1*mhat1+mhat2*mhat2+mhat3*mhat3
            sq_msq = sqrt(msq)
            if ( PBC == 1 )then
               FT_Coulomb(k3,k1,k2) = 1.d0 / (pi*msq)
            else
               FT_coulomb(k3,k1,k2) =  &
                   (1.d0 - cos(2.d0*pi*sq_msq*sphere)) / (pi*msq)
            endif
         enddo
      enddo
   enddo
end subroutine FH_RECIP_GRIDS_get_FT_coulomb
!-----------------------------------------------------------
subroutine FH_RECIP_GRIDS_diff_adj_grid( &
               nfft1,nfft2,nfft3,nfftdim1, &
               expon,recip,diffuse_adjust_grid)

   implicit none

   integer, intent(in) :: nfft1,nfft2,nfft3,nfftdim1
   double precision, intent(in) :: expon, recip(3,3)
   double precision, intent(out) :: diffuse_adjust_grid(nfft3,nfftdim1,nfft2)

   integer :: k1,k2,k3,m1,m2,m3,k10,nf1,nf2,nf3
   double precision :: pi,fac,mhat1,mhat2,mhat3,msq,term

   pi = 3.14159265358979323846d0
   fac = 2.d0*pi*pi / expon ! note the factor of 2 (need to fix 2 interacting
                            ! diffuse densities that have both been adjusted

   nf1 = nfft1/2
   if ( 2*nf1 < nfft1 )nf1 = nf1+1
   nf2 = nfft2/2
   if ( 2*nf2 < nfft2 )nf2 = nf2+1
   nf3 = nfft3/2
   if ( 2*nf3 < nfft3 )nf3 = nf3+1

   do k2 = 1, nfft2
      m2 = k2 - 1
      if ( k2 > nf2 )m2 = k2 - 1 - nfft2
      do k3 = 1,nfft3
         m3 = k3 - 1
         if ( k3 > nf3 )m3 = k3 - 1 - nfft3
         k10 = 1
         ! need (1,1,1) case also
         !if(k3+k2 == 2) k10 = 2
         do k1 = k10, nf1+1
            m1 = k1 - 1
            if ( k1 > nf1 )m1 = k1 - 1 - nfft1
            mhat1 = recip(1,1)*m1+recip(1,2)*m2+recip(1,3)*m3
            mhat2 = recip(2,1)*m1+recip(2,2)*m2+recip(2,3)*m3
            mhat3 = recip(3,1)*m1+recip(3,2)*m2+recip(3,3)*m3
            msq = mhat1*mhat1+mhat2*mhat2+mhat3*mhat3
            diffuse_adjust_grid(k3,k1,k2) = exp(-fac*msq)
         enddo
      enddo
   enddo
end subroutine FH_RECIP_GRIDS_diff_adj_grid
!-----------------------------------------------------------
subroutine FH_RECIP_GRIDS_accum( &
                           nfft1_1,nfft2_1,nfft3_1,nfftdim1_1, &
                           nfft1_2,nfft2_2,nfft3_2,nfftdim1_2, &
                           grid1,grid2,out_lun)

   implicit none

   integer, intent(in) :: nfft1_1,nfft2_1,nfft3_1,nfftdim1_1, &
                          nfft1_2,nfft2_2,nfft3_2,nfftdim1_2
   double precision, intent(in) :: grid2(2,nfft3_2,nfftdim1_2,nfft2_2)
   double precision, intent(inout) :: grid1(2,nfft3_1,nfftdim1_1,nfft2_1)
   integer, intent(in)  :: out_lun

   integer :: k1,k2,k3,k10,nf1

   ! note that grid2 exceeds grid1 in dimensions, so let grid1 dims control
   ! i.e. particular grid (grid2) is bigger than collective grid (grid1)
   ! check dims
   if ( (nfft1_1 > nfft1_2) .or. (nfft2_1 > nfft2_2) .or. &
        (nfft3_1 > nfft3_2) )then
      write(out_lun,*)'FH_RECIP_GRIDS_accum: grid1 bigger than grid2!!'
      stop
   endif
   nf1 = nfft1_1/2
   if ( 2*nf1 < nfft1_1 )nf1 = nf1+1
   do k2 = 1, nfft2_1
      do k3 = 1,nfft3_1
         k10 = 1
         ! need (1,1,1) case also
         !if ( k3+k2 == 2 )k10 = 2
         do k1 = k10, nf1+1
            grid1(1,k3,k1,k2) = grid1(1,k3,k1,k2) + grid2(1,k3,k1,k2)
            grid1(2,k3,k1,k2) = grid1(2,k3,k1,k2) + grid2(2,k3,k1,k2)
         enddo
      enddo
   enddo
end subroutine FH_RECIP_GRIDS_accum
!-----------------------------------------------------------
subroutine FH_RECIP_GRIDS_mult_one_coulomb( &
                           nfft1_1,nfft2_1,nfft3_1,nfftdim1_1, &
                           nfft1_2,nfft2_2,nfft3_2,nfftdim1_2, &
                           nfft1_3,nfft2_3,nfft3_3,nfftdim1_3, &
                           collective_grid,FT_coulomb,potential_grid,out_lun)

   implicit none

   integer, intent(in) :: nfft1_1,nfft2_1,nfft3_1,nfftdim1_1, &
                          nfft1_2,nfft2_2,nfft3_2,nfftdim1_2, &
                          nfft1_3,nfft2_3,nfft3_3,nfftdim1_3
   double precision, intent(in) ::  &
                           collective_grid(2,nfft3_1,nfftdim1_1,nfft2_1), &
                           FT_coulomb(nfft3_2,nfftdim1_2,nfft2_2)
   double precision, intent(out) :: &
                           potential_grid(2,nfft3_3,nfftdim1_3,nfft2_3)
   integer, intent(in)  :: out_lun

   integer :: k1,k2,k3,k10,nf1

   ! note that grid2 exceeds grid1 in dimensions, so let grid1 dims control
   ! i.e. particular grid (grid2) is bigger than collective grid (grid1)
   ! check dims
   if ( (nfft1_1 > nfft1_2) .or. (nfft2_1 > nfft2_2) .or. &
        (nfft3_1 > nfft3_2) )then
      write(out_lun,*)'FH_RECIP_LOC_FT_potential: grid1 bigger than grid2!!'
      stop
   endif
   ! check that collective and potential have same dims
   if ( (nfft1_1 .ne. nfft1_3) .or. (nfft2_1 .ne. nfft2_3) .or. &
        (nfft3_1 .ne. nfft3_3) .or. (nfftdim1_1 .ne. nfftdim1_3) )then
      write(out_lun,*)'FH_RECIP_LOC_FT_potential: grid1 diff size than grid3!!'
      stop
   endif
   nf1 = nfft1_1/2
   if ( 2*nf1 < nfft1_1 )nf1 = nf1+1
   do k2 = 1, nfft2_1
      do k3 = 1,nfft3_1
         k10 = 1
         ! need (1,1,1) case also
         !if ( k3+k2 == 2 )k10 = 2
         do k1 = k10, nf1+1
            potential_grid(1,k3,k1,k2) = collective_grid(1,k3,k1,k2) * &
                                          FT_coulomb(k3,k1,k2)
            potential_grid(2,k3,k1,k2) = collective_grid(2,k3,k1,k2) * &
                                          FT_coulomb(k3,k1,k2)
         enddo
      enddo
   enddo
end subroutine FH_RECIP_GRIDS_mult_one_coulomb
!-----------------------------------------------------------
subroutine FH_RECIP_GRIDS_2grid_ene( &
                           nfft1_1,nfft2_1,nfft3_1,nfftdim1_1, &
                           nfft1_2,nfft2_2,nfft3_2,nfftdim1_2, &
                           grid1,grid2,volume,energy,out_lun)

   implicit none

   integer, intent(in) :: nfft1_1,nfft2_1,nfft3_1,nfftdim1_1, &
                          nfft1_2,nfft2_2,nfft3_2,nfftdim1_2
   double precision,intent(in) :: grid1(2,nfft3_1,nfftdim1_1,nfft2_1), &
                                  grid2(2,nfft3_2,nfftdim1_2,nfft2_2)
   double precision, intent(in) :: volume
   double precision,intent(out) :: energy
   integer, intent(in)          :: out_lun
   
   integer :: k1,k2,k3,k10,nf1
   double precision :: mult,ene

   ene = 0.d0
   ! note that nfft1_2 >= nfft1_1 etc. i.e. tot_FT_pot bigger than FT_potential
   ! so let nfft1_1 etc drive loops
   ! check dims
   if ( (nfft1_1 > nfft1_2) .or. (nfft2_1 > nfft2_2) .or. &
        (nfft3_1 > nfft3_2) )then
      write(out_lun,*)'FH_RECIP_LOC_FT_potential: grid1 bigger than grid2!!'
      stop
   endif
   nf1 = nfft1_1/2
   if ( 2*nf1 < nfft1_1 )nf1 = nf1+1
   do k2 = 1, nfft2_1
      do k3 = 1,nfft3_1
         k10 = 1
         ! need (1,1,1) case also
         do k1 = k10, nf1+1
            if ( k1 > 1 )then
               mult = 2.d0*volume
            else
               mult = 1.d0*volume
            endif
            ene = ene + mult * ( grid1(1,k3,k1,k2)*grid2(1,k3,k1,k2) + &
                              grid1(2,k3,k1,k2)*grid2(2,k3,k1,k2) ) 
         enddo
      enddo
   enddo
   energy = ene
end subroutine FH_RECIP_GRIDS_2grid_ene
!-----------------------------------------------------------
subroutine FH_RECIP_LOC_diffadj_2grid_ene( &
                           nfft1_1,nfft2_1,nfft3_1,nfftdim1_1, &
                           nfft1_2,nfft2_2,nfft3_2,nfftdim1_2, &
                           nfft1_3,nfft2_3,nfft3_3,nfftdim1_3, &
                           grid1,grid2,grid3,volume,energy,out_lun)

   implicit none

   integer, intent(in) :: nfft1_1,nfft2_1,nfft3_1,nfftdim1_1, &
                          nfft1_2,nfft2_2,nfft3_2,nfftdim1_2, &
                          nfft1_3,nfft2_3,nfft3_3,nfftdim1_3    
   double precision,intent(in) :: grid1(2,nfft3_1,nfftdim1_1,nfft2_1), &
                                  grid2(2,nfft3_2,nfftdim1_2,nfft2_2), &
                                  grid3(nfft3_3,nfftdim1_3,nfft2_3)
   double precision, intent(in) :: volume
   double precision,intent(out) :: energy
   integer, intent(in)          :: out_lun
   
   integer :: k1,k2,k3,k10,nf1
   double precision :: mult,ene

   ene = 0.d0
   ! note that nfft1_2 >= nfft1_1 etc. i.e. tot_FT_pot bigger than FT_potential
   ! also nfft1_3 >= nfft1_1
   ! so let nfft1_1 etc drive loops
   ! check dims
   if ( (nfft1_1 > nfft1_2) .or. (nfft2_1 > nfft2_2) .or. &
        (nfft3_1 > nfft3_2) )then
      write(out_lun,*)&
        'FH_RECIP_LOC_diffadj_2grid_ene: grid1 bigger than grid2!!'
      stop
   endif
   if ( (nfft1_1 > nfft1_3) .or. (nfft2_1 > nfft2_3) .or. &
        (nfft3_1 > nfft3_3) )then
      write(out_lun,*)&
        'FH_RECIP_LOC_diffadj_2grid_ene: grid1 bigger than grid3!!'
      stop
   endif
   nf1 = nfft1_1/2
   if ( 2*nf1 < nfft1_1 )nf1 = nf1+1
   do k2 = 1, nfft2_1
      do k3 = 1,nfft3_1
         k10 = 1
         ! need (1,1,1) case also
         do k1 = k10, nf1+1
            if ( k1 > 1 )then
               mult = 2.d0*volume
            else
               mult = 1.d0*volume
            endif
            ene = ene + mult * grid3(k3,k1,k2) * &
                               ( grid1(1,k3,k1,k2)*grid2(1,k3,k1,k2) + &
                                 grid1(2,k3,k1,k2)*grid2(2,k3,k1,k2) )
         enddo
      enddo
   enddo
   energy = ene
end subroutine FH_RECIP_LOC_diffadj_2grid_ene
!-----------------------------------------------------------
subroutine FH_RECIP_GRIDS_add_to_FT_pot(  &
                           nfft1_1,nfft2_1,nfft3_1,nfftdim1_1, &
                           nfft1_2,nfft2_2,nfft3_2,nfftdim1_2, &
                           volume,FT_potential,tot_FT_pot,signum,out_lun)

   implicit none

   integer, intent(in) :: nfft1_1,nfft2_1,nfft3_1,nfftdim1_1, &
                          nfft1_2,nfft2_2,nfft3_2,nfftdim1_2
   double precision, intent(in) :: volume
   double precision, intent(in) :: FT_potential(2,nfft3_1,nfftdim1_1,nfft2_1)
   double precision, intent(inout) :: tot_FT_pot(2,nfft3_2,nfftdim1_2,nfft2_2)
   integer, intent(in)  :: signum
   integer, intent(in)  :: out_lun
                           
   integer :: k1,k2,k3,k10,nf1

   ! note that nfft1_2 >= nfft1_1 etc. i.e. tot_FT_pot bigger than FT_potential
   ! so let nfft1_1 etc drive loops
   ! check dims
   if ( (nfft1_1 > nfft1_2) .or. (nfft2_1 > nfft2_2) .or. &
        (nfft3_1 > nfft3_2) )then
      write(out_lun,*)'FH_RECIP_LOC_FT_potential: grid1 bigger than grid2!!'
      stop
   endif
   nf1 = nfft1_1/2
   if ( 2*nf1 < nfft1_1 )nf1 = nf1+1
   do k2 = 1, nfft2_1
      do k3 = 1,nfft3_1
         k10 = 1
         ! need (1,1,1) case also
         !if ( k3+k2 == 2 )k10 = 2
         do k1 = k10, nf1+1
            !mult by volume for plancherel
            tot_FT_pot(1,k3,k1,k2) = tot_FT_pot(1,k3,k1,k2) + &
                                 signum*volume*FT_potential(1,k3,k1,k2)
            tot_FT_pot(2,k3,k1,k2) = tot_FT_pot(2,k3,k1,k2) + &
                                 signum*volume*FT_potential(2,k3,k1,k2)
         enddo
      enddo
   enddo
end subroutine FH_RECIP_GRIDS_add_to_FT_pot
!-----------------------------------------------------------
subroutine FH_RECIP_GRIDS_da_add_to_FT_pot(  &
                           nfft1_1,nfft2_1,nfft3_1,nfftdim1_1, &
                           nfft1_2,nfft2_2,nfft3_2,nfftdim1_2, &
                           nfft1_3,nfft2_3,nfft3_3,nfftdim1_3, &
                           volume,FT_potential, tot_FT_pot, diff_adj_grid, &
                           out_lun )
   implicit none

   integer, intent(in) :: nfft1_1,nfft2_1,nfft3_1,nfftdim1_1, &
                          nfft1_2,nfft2_2,nfft3_2,nfftdim1_2, &
                          nfft1_3,nfft2_3,nfft3_3,nfftdim1_3
   double precision, intent(in) :: volume
   double precision, intent(in) :: FT_potential(2,nfft3_1,nfftdim1_1,nfft2_1)
   double precision, intent(inout) :: tot_FT_pot(2,nfft3_2,nfftdim1_2,nfft2_2)
   double precision, intent(in) :: diff_adj_grid(nfft3_3,nfftdim1_3,nfft2_3)
   integer, intent(in)  :: out_lun
                           
   integer :: k1,k2,k3,k10,nf1

   ! note that nfft1_2 >= nfft1_1 etc. i.e. tot_FT_pot bigger than FT_potential
   ! so let nfft1_1 etc drive loops
   ! check dims
   if ( (nfft1_1 > nfft1_2) .or. (nfft2_1 > nfft2_2) .or. &
        (nfft3_1 > nfft3_2) )then
      write(out_lun,*)&
        'FH_RECIP_GRIDS_da_add_to_FT_pot: grid1 bigger than grid2!!'
      stop
   endif
   if ( (nfft1_1 > nfft1_3) .or. (nfft2_1 > nfft2_3) .or. &
        (nfft3_1 > nfft3_3) )then
      write(out_lun,*)&
        'FH_RECIP_GRIDS_da_add_to_FT_pot: grid1 bigger than grid3!!'
      stop
   endif
   nf1 = nfft1_1/2
   if ( 2*nf1 < nfft1_1 )nf1 = nf1+1
   do k2 = 1, nfft2_1
      do k3 = 1,nfft3_1
         k10 = 1
         ! need (1,1,1) case also
         !if ( k3+k2 == 2 )k10 = 2
         do k1 = k10, nf1+1
            !mult by volume for plancherel
            tot_FT_pot(1,k3,k1,k2) = tot_FT_pot(1,k3,k1,k2) + &
                                     volume * diff_adj_grid(k3,k1,k2) * &
                                     FT_potential(1,k3,k1,k2)
            tot_FT_pot(2,k3,k1,k2) = tot_FT_pot(2,k3,k1,k2) + &
                                     volume * diff_adj_grid(k3,k1,k2) * &
                                     FT_potential(2,k3,k1,k2)
         enddo
      enddo
   enddo
end subroutine FH_RECIP_GRIDS_da_add_to_FT_pot
!-----------------------------------------------------------
