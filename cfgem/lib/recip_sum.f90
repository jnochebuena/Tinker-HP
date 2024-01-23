module recip_sum
   implicit none
   private

   public FH_RECIP_setup,FH_RECIP_eval,FH_RECIP_deallocate

contains
!-----------------------------------------------------------------
!-----------------------------------------------------------------
! SETUP SUBROUTINES
!-----------------------------------------------------------------
!-----------------------------------------------------------------
subroutine FH_RECIP_setup(out_lun)

   use user, only : do_coulomb, &
                    do_overlap, &
                    coulomb_use_recip, &
                    overlap_use_recip
   use recip_grids, only : FH_RECIP_GRIDS_setup, rec_virial
   use pme_recip, only : FH_PME_RECIP_setup
   use ewald_recip, only : FH_EWALD_RECIP_setup

   implicit none

   integer, intent(in)  :: out_lun
   
   integer :: proceed
   proceed = 0

  ! We stick zero'ing the recip virial here to INSURE it happens, regardless
  ! of whether recip grids get used..

   rec_virial(:,:) = 0.d0

   ! check if doing recip at all
   if ( do_coulomb == 1 .and. coulomb_use_recip == 1 )proceed = 1
   if ( do_overlap == 1 .and. overlap_use_recip == 1 )proceed = 1
   if ( proceed == 0 )return

   call FH_RECIP_GRIDS_setup(out_lun)
   call FH_PME_RECIP_setup(out_lun)
   call FH_EWALD_RECIP_setup(out_lun)

end subroutine FH_recip_setup
!-----------------------------------------------------------------
!-----------------------------------------------------------------
! EVALUATION SUBROUTINES
!-----------------------------------------------------------------
!-----------------------------------------------------------------
subroutine FH_RECIP_eval( &
                       recip_mp_mp_ene, &
                       recip_mp_ch_ene, &
                       recip_ch_mp_ene, &
                       recip_ch_ch_ene, &
                       recip_mp_dh_ene, &
                       recip_dh_mp_ene, &
                       recip_ch_dh_ene, &
                       recip_dh_ch_ene, &
                       recip_dh_dh_ene,factor,virial, out_lun)

   use user, only : verbose, &
                    do_coulomb, &
                    do_overlap, &
                    coulomb_use_recip, &
                    overlap_use_recip
   use recip_grids, only : FH_RECIP_GRIDS_FT_ene,rec_virial

   implicit none

   double precision, intent(in) :: factor
   double precision, intent(inout) :: virial(3,3)
   double precision, intent(out) :: recip_mp_mp_ene, &
                                    recip_mp_ch_ene, &
                                    recip_ch_mp_ene, &
                                    recip_ch_ch_ene, &
                                    recip_mp_dh_ene, &
                                    recip_dh_mp_ene, &
                                    recip_ch_dh_ene, &
                                    recip_dh_ch_ene, &
                                    recip_dh_dh_ene
   integer, intent(in)  :: out_lun

   double precision time1,time2,time3,energy
   integer count,rate,i,j

   integer :: proceed
   proceed = 0

   ! check if doing recip at all
   if ( do_coulomb == 1 .and. coulomb_use_recip == 1 )proceed = 1
   if ( do_overlap == 1 .and. overlap_use_recip == 1 )proceed = 1
   if ( proceed == 0 )return

!  call system_clock( COUNT=count, COUNT_RATE=rate)
!  time1 = dble(count)/dble(rate)
   ! first get the fourier transformed gridded densities
   call FH_RECIP_get_FT_density(out_lun)
   ! get the coulomb potentials in Fourier space
   rec_virial(:,:)=0.d0
   call FH_RECIP_FT_potential(out_lun)

 ! reciprocal virial calculated in FT_potential, add rec_virial to total
 ! virial
  !write(6,*)'---------------GEM Reciprocal VIRIAL-------------------'
  !do i = 1,3
  !   write(6,*)(rec_virial(i,j),j=1,3)
  !enddo
  !write(6,*)'---------------GEM Reciprocal VIRIAL-------------------'

  !write(6,*)'---------------GEM SUM RECIP VIRIAL-------------------'
  !do i = 1,3
  !   do j = 1,3
  !      virial(i,j)=virial(i,j)+rec_virial(i,j)
  !   enddo
  !enddo
  !do i = 1,3
  !   write(6,*)(virial(i,j),j=1,3)
  !enddo
  !write(6,*)'---------------GEM SUM RECIP VIRIAL-------------------'

   ! get the energy components
   call FH_RECIP_GRIDS_FT_ene(recip_mp_mp_ene, &
                              recip_mp_ch_ene, &
                              recip_ch_mp_ene, &
                              recip_ch_ch_ene, &
                              recip_mp_dh_ene, &
                              recip_dh_mp_ene, &
                              recip_ch_dh_ene, &
                              recip_dh_ch_ene, &
                              recip_dh_dh_ene, &
                              out_lun)
!  if ( verbose > 1 )then
!     write(6,*)'recip_mp_mp_ene = ',recip_mp_mp_ene
!     write(6,*)'recip_mp_ch_ene = ',recip_mp_ch_ene
!     write(6,*)'recip_ch_mp_ene = ',recip_ch_mp_ene
!     write(6,*)'recip_ch_ch_ene = ',recip_ch_ch_ene
!     write(6,*)'recip_mp_dh_ene = ',recip_mp_dh_ene
!     write(6,*)'recip_dh_mp_ene = ',recip_dh_mp_ene
!     write(6,*)'recip_ch_dh_ene = ',recip_ch_dh_ene
!     write(6,*)'recip_dh_ch_ene = ',recip_dh_ch_ene
!     write(6,*)'recip_dh_dh_ene = ',recip_dh_dh_ene
!  endif !( verbose > 1 )then

   call FH_RECIP_ene_force(energy, out_lun)

!  call system_clock( COUNT=count, COUNT_RATE=rate)

!  time2 = dble(count)/dble(rate)
!  write(6,*)'recip_ene = ',energy
!  write(6,*)'tot of FT_ene = ',recip_mp_mp_ene + &
!                               recip_mp_ch_ene + &
!                               recip_ch_mp_ene + &
!                               recip_ch_ch_ene + &
!                               recip_mp_dh_ene + &
!                               recip_dh_mp_ene + &
!                               recip_ch_dh_ene + &
!                               recip_dh_ch_ene + &
!                               recip_dh_dh_ene
!  write(6,*)'DONE RECIP!!!!!!!!!', ' time = ',time2-time1

end subroutine FH_RECIP_eval
!-----------------------------------------------------------------
subroutine FH_RECIP_get_FT_density(out_lun)

   use recip_grids, only : FH_RECIP_GRIDS_fft_R_to_K
   use pme_recip, only : FH_PME_RECIP_global_to_frac, &
                         FH_PME_RECIP_Bspline_fill, &
                         FH_PME_RECIP_fill_grids, & 
                         FH_PME_RECIP_grid_mult, &
                         FH_PME_RECIP_FT_density
   use ffp_recip, only : FH_FFP_RECIP_fill_grids, &
                         FH_FFP_RECIP_FT_density
   use ewald_recip, only : FH_EWALD_RECIP_FT_density

   implicit none

   integer, intent(in)  :: out_lun

   !-------------------------------------------------
   ! first fill the grids that need fft transforming
   !-------------------------------------------------
   ! PME grids
   !-------------------------------------------
      ! transform the hermite coeffs
      call FH_PME_RECIP_global_to_frac()
      ! fill the b-splines
      call FH_PME_RECIP_Bspline_fill(out_lun)
      ! grid fill
      call FH_PME_RECIP_fill_grids()

   !---------------------------------------------
   ! FFP grids
   !---------------------------------------------
      ! grid fill
      call FH_FFP_RECIP_fill_grids(out_lun)
   !-------------------------------------------------
   ! next fft the pme and ffp grids to k-space
   !-------------------------------------------------
      call FH_RECIP_GRIDS_fft_R_to_K(out_lun)
   !---------------------------------------------------------------
   !  the ffp grids are already the FT's of their density
   !  except for volume factor
   !  multiply the pme grids by their reciprocal gaussian times prefac
   !----------------------------------------------------------------
      ! get the pme_grid multipliers
      call FH_PME_RECIP_grid_mult(out_lun)
      ! multiply to get structure factors
      call FH_PME_RECIP_FT_density()
      ! also fix the ffp density (scale by volume)
      call FH_FFP_RECIP_FT_density()
   !------------------------------------------------------------
   ! finally do the ewald grids 
   !------------------------------------------------------------
      call FH_EWALD_RECIP_FT_density()

end subroutine FH_RECIP_get_FT_density
!-----------------------------------------------------------------
subroutine FH_RECIP_FT_potential(out_lun)

   use recip_grids, only : FH_RECIP_GRIDS_collective,&
                           FH_RECIP_GRIDS_virial,&
                           FH_RECIP_GRIDS_mult_coulomb

   implicit none

   integer, intent(in)  :: out_lun

   call FH_RECIP_GRIDS_collective(out_lun)

! Need to calculate virial for each grid separately, Multipole and sumCH have
! two collective grids each, hermites only compact and diffuse collective
! last argument is "1" if it is two separate grids (mtp, sumCH) or "0"
! for one grid (herms)

   ! This subroutine calls FH_RECIP_GRIDS_virial_one four times, one for
   ! each type of grid
   !call FH_RECIP_GRIDS_virial(out_lun)

   call FH_RECIP_GRIDS_mult_coulomb(out_lun)
end subroutine FH_RECIP_FT_potential
!-----------------------------------------------------------------
subroutine FH_RECIP_ene_force(energy, out_lun)

   use recip_grids, only : RECIP_GRIDS_load_FT_pot, &
                           FH_RECIP_GRIDS_fft_K_to_R
   use pme_recip, only : FH_PME_RECIP_fix_FT_phi, &
                         FH_PME_RECIP_ene_force_field
   use ffp_recip, only : FH_FFP_RECIP_fix_FT_phi, &
                         FH_FFP_RECIP_ene_force_field
   use ewald_recip, only : FH_EWALD_RECIP_ene_force_field

   implicit none

   double precision, intent(out) :: energy

   integer, intent(in)  :: out_lun

   ! first load the potentials onto the grids, to prepare for xform back to real
   call RECIP_GRIDS_load_FT_pot(out_lun)

   ! next fix the pme recip phi--to put in B-spline correction
   call FH_PME_RECIP_fix_FT_phi()
   ! next fix the ffp recip phi--volume normalization
   call FH_FFP_RECIP_fix_FT_phi()

   ! next fft the pme and ffp grids back to real space
   call FH_RECIP_GRIDS_fft_K_to_R(out_lun)

   ! finally the energy and forces
   energy = 0.d0
   ! first the pme grids
   call FH_PME_RECIP_ene_force_field(energy, out_lun)
   ! next the ffp grids
   call FH_FFP_RECIP_ene_force_field(energy, out_lun)
   ! finally the ewald grids
   call FH_EWALD_RECIP_ene_force_field(energy)

end subroutine FH_RECIP_ene_force
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!-----------------------------------------------------------------
! DEALLOCATE SUBROUTINES
!-----------------------------------------------------------------
!-----------------------------------------------------------------
subroutine FH_RECIP_deallocate()
   use recip_grids, only : FH_RECIP_GRIDS_deallocate
   use pme_recip, only : FH_PME_RECIP_deallocate
   use ewald_recip, only : FH_EWALD_RECIP_deallocate

   implicit none

   call FH_RECIP_GRIDS_deallocate()
   call FH_PME_RECIP_deallocate()
   call FH_EWALD_RECIP_deallocate()

end subroutine FH_RECIP_deallocate
!-----------------------------------------------------------------
!-----------------------------------------------------------------
end module recip_sum
