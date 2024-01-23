module coulomb_site_site

   implicit none

   private
          
   public FH_COULOMB_SITE_SITE_CC, &
          FH_COULOMB_SITE_SITE_self_CC, &
          FH_COULOMB_SITE_SITE_CD_DD, &
          FH_COULOMB_SITE_SITE_self_CD_DD
          
   contains
!--------------------------------------------------------------------

!-----------------------------------------------------------------
!-----------------------------------------------------------------
! EVALUATION SUBROUTINES
!-----------------------------------------------------------------
!-----------------------------------------------------------------

!-------------------------------------------------------------------------
subroutine FH_COULOMB_SITE_SITE_CC( &
                                     correcting_energies, &
                                     use_cutoff, &
                                     full_interactions, &
                                     max_site_site_num_nghbrs, &
                                     site_site_num_nghbrs, &
                                     site_site_list_offset, &
                                     site_site_list, &
                                     site_site_ucell_transind, &
                                     mp_mp_ene, &
                                     mp_ch_ene, &
                                     ch_mp_ene, &
                                     ch_ch_ene,factor,virial)

   use user, only : pbc, &
                    CD_split_expon, &
                    extent_of_compact_hermites
   use sites, only : num_sites, &
                     map_site_crd, &
                     site_frc
   use unit_cell, only : ucell_translate
   use hermite, only : tot_num_mpole_coeffs, &
                       Global_multipole_coeff, &
                       Global_multipole_field, &
                       mpole_level_for_site, &
                       mpole_coeff_off_for_site, &
                       tot_num_sumCH_coeffs, &
                       Global_sumCH_coeff, &
                       Global_sumCH_field, &
                       sumCH_level_for_site, &
                       sumCH_coeff_off_for_site, &
                       tot_num_herm_Cprims, &
                       Global_Chermite_coeff, &
                       Global_Chermite_field, &
                       gauss_extent_for_Cprim, &
                       hermite_level_for_Cprim, &
                       num_herm_Cprim_for_site, &
                       off_herm_Cprim_for_site, &
                       herm_coeff_offset_of_Cprim, &
                       hermite_expon_for_Cprim

   implicit none

   integer, intent(in) :: correcting_energies, &
                          use_cutoff, &
                          full_interactions
   integer, intent(in)    ::  max_site_site_num_nghbrs, &
                              site_site_num_nghbrs(*), &
                              site_site_list(*), &
                              site_site_ucell_transind(*), &
                              site_site_list_offset(*)
   double precision, intent(in) :: factor 
   double precision, intent(inout) :: virial(3,3)
   double precision, intent(inout) :: mp_mp_ene, &
                                      mp_ch_ene, &
                                      ch_mp_ene, &
                                      ch_ch_ene

   include "interact_type.fh"

   integer :: buf_count, neighbor_buf(max_site_site_num_nghbrs)
   integer site1,num_site_list,self_flag, &
           signum,overall_signum,interaction_signum, &
           do_mpole_mpole,do_mpole_cherm,do_cherm_cherm, &
           do_mpole_sumch,do_sumch_sumch

   double precision :: delx(max_site_site_num_nghbrs), &
                       dely(max_site_site_num_nghbrs), &
                       delz(max_site_site_num_nghbrs), &
                       delr2(max_site_site_num_nghbrs)
                       

   integer jj

   self_flag = 0
   if ( tot_num_mpole_coeffs > 0 )then
      do_mpole_mpole = 1
   else
      do_mpole_mpole = 0
   endif

   if ( (tot_num_mpole_coeffs > 0) .and. (tot_num_herm_Cprims > 0) )then
      if ( (full_interactions == 1) .or. (use_cutoff == 1) )then
         do_mpole_cherm = 1
      else
         do_mpole_cherm = 0
      endif
   else
      do_mpole_cherm = 0
   endif !( (tot_num_mpole_coeffs > 0) .and. (tot_num_herm_Cprims > 0) )

   if ( tot_num_herm_Cprims > 0 )then
      if ( (full_interactions == 1) .or. (use_cutoff == 1) )then
         do_cherm_cherm = 1
      else
         do_cherm_cherm = 0
      endif
   else
      do_cherm_cherm = 0
   endif ! ( tot_num_herm_Cprims > 0 )

   if ( (tot_num_mpole_coeffs > 0) .and. (tot_num_sumCH_coeffs > 0) )then
      if ( (full_interactions == 0) .or. (use_cutoff == 1) )then
         do_mpole_sumch = 1
      else
         do_mpole_sumch = 0
      endif
   else
      do_mpole_sumch = 0
   endif !( (tot_num_mpole_coeffs > 0) .and. (tot_num_sumCH_coeffs > 0) )

   if ( tot_num_sumCH_coeffs > 0 )then
      if ( (full_interactions == 0) .or. (use_cutoff == 1) )then
         do_sumch_sumch = 1
      else
         do_sumch_sumch = 0
      endif
   else
      do_sumch_sumch = 0
   endif !( tot_num_sumCH_coeffs > 0 )

   ! overall signum (are we accumulating interactions or correcting them?)
   if ( correcting_energies == 1 )then
      overall_signum = -1
   else
      overall_signum = 1
   endif

   !print *,'num_sites',num_sites
   do site1 = 1,num_sites-1
      call FH_SITE_SITE_LIST_filter(  &
                               site1, &
                               site_site_num_nghbrs, &
                               site_site_list_offset, &
                               site_site_list, &
                               site_site_ucell_transind, &
                               map_site_crd, &
                               ucell_translate, &
                               extent_of_compact_hermites, &
                               pbc, &
                               use_cutoff, &
                               buf_count, &
                               neighbor_buf, &
                               delx,dely,delz,delr2 )
      num_site_list = buf_count
      ! 1st the mpole-mpole
      if ( do_mpole_mpole == 1 )then
         interaction_signum = 1 !e.g. nuclear-nuclear
         signum = overall_signum * interaction_signum
         call FH_COULOMB_SITE_SITE_mp_mp( &
                site1, &
                num_site_list, &
                mpole_level_for_site, &
                mpole_coeff_off_for_site, &
                mpole_level_for_site, &
                mpole_coeff_off_for_site, &
                extent_of_compact_hermites, &
                CD_split_expon, &
                full_interactions, &
                use_cutoff, &
                self_flag, &
                signum, &
                neighbor_buf, &
                delx, &
                dely, &
                delz, &
                delr2, &
                Global_multipole_coeff, &
                Global_multipole_field, &
                Global_multipole_coeff, &
                Global_multipole_field, &
                mp_mp_ene, &
                site_frc,factor,virial )
         !print *,'mpole-mpole do_mpole_mpole',mp_mp_ene
      endif !( do_mpole_mpole == 1 )then
      ! next the mpole - compact hermite
      if ( do_mpole_cherm == 1 )then
         interaction_signum = -1 !e.g. nuclear-electron
         signum = overall_signum * interaction_signum
         call FH_COULOMB_SITE_SITE_mp_herm( &
                site1, &
                num_site_list, &
                mpole_level_for_site, &
                mpole_coeff_off_for_site, &
                hermite_level_for_Cprim, &
                herm_coeff_offset_of_Cprim, &
                num_herm_Cprim_for_site, &
                off_herm_Cprim_for_site, &
                gauss_extent_for_Cprim, &
                hermite_expon_for_Cprim, &
                use_cutoff, &
                signum, &
                neighbor_buf,  &
                delx, &
                dely, &
                delz, &
                delr2, &
                Global_multipole_coeff, &
                Global_multipole_field, &
                Global_Chermite_coeff, &
                Global_Chermite_field, &
                mp_ch_ene, &
                site_frc,factor,virial )
         !print *,'mpole-cherm do_mpole_cherm',mp_ch_ene
         ! next the compact hermite - mpole
         call FH_COULOMB_SITE_SITE_herm_mp( &
                site1, &
                num_site_list, &
                mpole_level_for_site, &
                mpole_coeff_off_for_site, &
                hermite_level_for_Cprim, &
                herm_coeff_offset_of_Cprim, &
                num_herm_Cprim_for_site, &
                off_herm_Cprim_for_site, &
                gauss_extent_for_Cprim, &
                hermite_expon_for_Cprim, &
                use_cutoff, &
                signum, &
                neighbor_buf, &
                delx, &
                dely, &
                delz, &
                delr2, &
                Global_multipole_coeff, &
                Global_multipole_field, &
                Global_Chermite_coeff, &
                Global_Chermite_field, &
                ch_mp_ene, &
                site_frc,factor,virial )
         !print *,'cherm-mpole do_mpole_cherm',ch_mp_ene
      endif !( do_mpole_cherm == 1 )then
      ! next the compact herms with compact herms
      if ( do_cherm_cherm == 1 )then
         interaction_signum = 1 !e.g. electron-electron
         signum = overall_signum * interaction_signum
         call FH_COULOMB_SITE_SITE_herm_herm( &
                site1, &
                num_site_list, &
                hermite_level_for_Cprim, &
                herm_coeff_offset_of_Cprim, &
                num_herm_Cprim_for_site, &
                off_herm_Cprim_for_site, &
                gauss_extent_for_Cprim, &
                hermite_expon_for_Cprim, &
                hermite_level_for_Cprim, &
                herm_coeff_offset_of_Cprim, &
                num_herm_Cprim_for_site, &
                off_herm_Cprim_for_site, &
                gauss_extent_for_Cprim, &
                hermite_expon_for_Cprim, &
                use_cutoff, &
                signum, &
                neighbor_buf, &
                factor, &
                delx, &
                dely, &
                delz, &
                delr2, &
                Global_Chermite_coeff, &
                Global_Chermite_field, &
                Global_Chermite_coeff, &
                Global_Chermite_field, &
                ch_ch_ene, &
                site_frc,virial)
         !print *,'ch-ch do_cherm_cherm',ch_ch_ene
      endif !( do_cherm_cherm == 1 )then
      if ( do_mpole_sumch == 1 )then
         ! mpole-sumch
         interaction_signum = -1 !e.g. nuclear-electron
         signum = overall_signum * interaction_signum
         call FH_COULOMB_SITE_SITE_mp_mp( &
                site1, &
                num_site_list, &
                mpole_level_for_site, &
                mpole_coeff_off_for_site, &
                sumCH_level_for_site, &
                sumCH_coeff_off_for_site, &
                extent_of_compact_hermites, &
                CD_split_expon, &
                full_interactions, &
                use_cutoff, &
                self_flag, &
                signum, &
                neighbor_buf, &
                delx, &
                dely, &
                delz, &
                delr2, &
                Global_multipole_coeff, &
                Global_multipole_field, &
                Global_sumCH_coeff, &
                Global_sumCH_field, &
                mp_ch_ene, &
                site_frc,factor,virial )
         !print *,'mp-ch_nuc_elec do_mpole_sumch',mp_ch_ene
         ! sumch-mpole
         call FH_COULOMB_SITE_SITE_mp_mp( &
                site1, &
                num_site_list, &
                sumCH_level_for_site, &
                sumCH_coeff_off_for_site, &
                mpole_level_for_site, &
                mpole_coeff_off_for_site, &
                extent_of_compact_hermites, &
                CD_split_expon, &
                full_interactions, &
                use_cutoff, &
                self_flag, &
                signum, &
                neighbor_buf, &
                delx, &
                dely, &
                delz, &
                delr2, &
                Global_sumCH_coeff, &
                Global_sumCH_field, &
                Global_multipole_coeff, &
                Global_multipole_field, &
                ch_mp_ene, &
                site_frc,factor,virial )
         !print *,'ch-mp_nuc_elec do_mpole_sumch',ch_mp_ene
      endif !( do_mpole_sumch == 1 )then
      ! finally the sumch-sumch
      if ( do_sumch_sumch == 1 )then
         interaction_signum = 1 !e.g. electron-electron
         signum = overall_signum * interaction_signum
         call FH_COULOMB_SITE_SITE_mp_mp( &
                site1, &
                num_site_list, &
                sumCH_level_for_site, &
                sumCH_coeff_off_for_site, &
                sumCH_level_for_site, &
                sumCH_coeff_off_for_site, &
                extent_of_compact_hermites, &
                CD_split_expon, &
                full_interactions, &
                use_cutoff, &
                self_flag, &
                signum, &
                neighbor_buf, &
                delx, &
                dely, &
                delz, &
                delr2, &
                Global_sumCH_coeff, &
                Global_sumCH_field, &
                Global_sumCH_coeff, &
                Global_sumCH_field, &
                ch_ch_ene, &
                site_frc,factor,virial )
      endif !( do_sumch_sumch == 1 )then
   enddo !site1 = 1,num_sites-1

   return

end subroutine FH_COULOMB_SITE_SITE_CC
!--------------------------------------------------------------
subroutine FH_COULOMB_SITE_SITE_self_CC( &
                                     correcting_energies, &
                                     mp_mp_ene, &
                                     mp_ch_ene, &
                                     ch_mp_ene, &
                                     ch_ch_ene,factor,virial)

   use sites, only : num_sites, &
                     site_frc
   use user, only : CD_split_expon, &
                    extent_of_compact_hermites
   use hermite, only : tot_num_mpole_coeffs, &
                       Global_multipole_coeff, &
                       Global_multipole_field, &
                       mpole_level_for_site, &
                       mpole_coeff_off_for_site, &
                       tot_num_sumCH_coeffs, &
                       Global_sumCH_coeff, &
                       Global_sumCH_field, &
                       sumCH_level_for_site, &
                       sumCH_coeff_off_for_site, &
                       tot_num_herm_Cprims, &
                       Global_Chermite_coeff, &
                       Global_Chermite_field, &
                       gauss_extent_for_Cprim, &
                       hermite_level_for_Cprim, &
                       num_herm_Cprim_for_site, &
                       off_herm_Cprim_for_site, &
                       herm_coeff_offset_of_Cprim, &
                       hermite_expon_for_Cprim

   implicit none

   integer, intent(in) :: correcting_energies
   double precision, intent(in) :: factor 
   double precision, intent(inout) :: virial(3,3)
   double precision, intent(inout) :: mp_mp_ene, &
                                      mp_ch_ene, &
                                      ch_mp_ene, &
                                      ch_ch_ene

   integer :: buf_count, neighbor_buf(1)
   integer :: site1,num_site_list,self_flag, &
              signum,overall_signum,interaction_signum, &
              do_mpole_mpole,do_mpole_sumch,do_sumch_sumch, &
              use_cutoff,full_interactions
   double precision :: delx(1), &
                       dely(1), &
                       delz(1), &
                       delr2(1)
                      

   if ( tot_num_mpole_coeffs > 0 )then
      do_mpole_mpole = 1
   else
      do_mpole_mpole = 0
   endif

   if ( (tot_num_mpole_coeffs > 0) .and. (tot_num_sumCH_coeffs > 0) )then
      do_mpole_sumch = 1
   else
      do_mpole_sumch = 0
   endif !( (tot_num_mpole_coeffs > 0) .and. (tot_num_sumCH_coeffs > 0) )

   if ( tot_num_sumCH_coeffs > 0 )then
      do_sumch_sumch = 1
   else
      do_sumch_sumch = 0
   endif !( tot_num_sumCH_coeffs > 0 )

   ! overall signum (are we accumulating interactions or correcting them?)
   if ( correcting_energies == 1 )then
      overall_signum = -1
   else
      overall_signum = 1
   endif

   use_cutoff = 0
   full_interactions = 0
   do site1 = 1,num_sites
      num_site_list = 1
      delx(1) = 0.d0
      dely(1) = 0.d0
      delz(1) = 0.d0
      delr2(1) = 0.d0
      neighbor_buf(1) = site1
      ! 1st the mpole-mpole
      if ( do_mpole_mpole == 1 )then
         interaction_signum = 1 ! e.g.nuclear-nuclear
         signum = overall_signum * interaction_signum
         self_flag = 1
         call FH_COULOMB_SITE_SITE_mp_mp( &
                site1, &
                num_site_list, &
                mpole_level_for_site, &
                mpole_coeff_off_for_site, &
                mpole_level_for_site, &
                mpole_coeff_off_for_site, &
                extent_of_compact_hermites, &
                CD_split_expon, &
                full_interactions, &
                use_cutoff, &
                self_flag, &
                signum, &
                neighbor_buf, &
                delx, &
                dely, &
                delz, &
                delr2, &
                Global_multipole_coeff, &
                Global_multipole_field, &
                Global_multipole_coeff, &
                Global_multipole_field, &
                mp_mp_ene, &
                site_frc,factor,virial )
      endif !( do_mpole_mpole == 1 )then
      ! next the mpole-sumch
      if ( do_mpole_sumch == 1 )then
         interaction_signum = -1 ! e.g.nuclear-electron
         signum = overall_signum * interaction_signum
         self_flag = 0
         ! next the sumch-mpole
         call FH_COULOMB_SITE_SITE_mp_mp( &
                site1, &
                num_site_list, &
                sumCH_level_for_site, &
                sumCH_coeff_off_for_site, &
                mpole_level_for_site, &
                mpole_coeff_off_for_site, &
                extent_of_compact_hermites, &
                CD_split_expon, &
                full_interactions, &
                use_cutoff, &
                self_flag, &
                signum, &
                neighbor_buf, &
                delx, &
                dely, &
                delz, &
                delr2, &
                Global_sumCH_coeff, &
                Global_sumCH_field, &
                Global_multipole_coeff, &
                Global_multipole_field, &
                ch_mp_ene, &
                site_frc,factor,virial )
      endif !( do_mpole_sumch == 1 )then
      ! finally the sumch-sumch
      if ( do_sumch_sumch == 1 )then
         interaction_signum = 1 ! e.g.electron-electron
         signum = overall_signum * interaction_signum
         self_flag = 1
         call FH_COULOMB_SITE_SITE_mp_mp( &
                site1, &
                num_site_list, &
                sumCH_level_for_site, &
                sumCH_coeff_off_for_site, &
                sumCH_level_for_site, &
                sumCH_coeff_off_for_site, &
                extent_of_compact_hermites, &
                CD_split_expon, &
                full_interactions, &
                use_cutoff, &
                self_flag, &
                signum, &
                neighbor_buf, &
                delx, &
                dely, &
                delz, &
                delr2, &
                Global_sumCH_coeff, &
                Global_sumCH_field, &
                Global_sumCH_coeff, &
                Global_sumCH_field, &
                ch_ch_ene, &
                site_frc,factor,virial )
      endif !( do_sumch_sumch == 1 )then
   enddo !site1 = 1,num_sites-1
   mp_mp_ene = 0.5d0*mp_mp_ene
   ch_mp_ene = 0.5d0*ch_mp_ene
   mp_ch_ene = ch_mp_ene
   ch_ch_ene = 0.5d0*ch_ch_ene

end subroutine FH_COULOMB_SITE_SITE_self_CC
!--------------------------------------------------------------
subroutine FH_COULOMB_SITE_SITE_CD_DD( &
                                     correcting_energies, &
                                     max_site_site_num_nghbrs, &
                                     site_site_num_nghbrs, &
                                     site_site_list_offset, &
                                     site_site_list, &
                                     site_site_ucell_transind, &
                                     mp_dh_ene, &
                                     dh_mp_ene, &
                                     ch_dh_ene, &
                                     dh_ch_ene, &
                                     dh_dh_ene,factor,virial)

   use user, only : pbc, &
                    extent_of_compact_hermites
   use sites, only : num_sites, &
                     map_site_crd, &
                     site_frc
   use unit_cell, only : ucell_translate
   use hermite, only :  &
                       ! mpoles
                       tot_num_mpole_coeffs, &
                       Global_multipole_coeff, &
                       Global_multipole_field, &
                       mpole_level_for_site, &
                       mpole_coeff_off_for_site, &
                       ! compact hermite
                       tot_num_herm_Cprims, &
                       Global_Chermite_coeff, &
                       Global_Chermite_field, &
                       gauss_extent_for_Cprim, &
                       hermite_level_for_Cprim, &
                       num_herm_Cprim_for_site, &
                       off_herm_Cprim_for_site, &
                       herm_coeff_offset_of_Cprim, &
                       hermite_expon_for_Cprim, &
                       ! diffuse hermite
                       tot_num_herm_Dprims, &
                       Global_Dhermite_coeff, &
                       Global_Dhermite_field, &
                       gauss_extent_for_Dprim, &
                       hermite_level_for_Dprim, &
                       num_herm_Dprim_for_site, &
                       off_herm_Dprim_for_site, &
                       herm_coeff_offset_of_Dprim, &
                       hermite_expon_for_Dprim

   implicit none

   integer, intent(in) :: correcting_energies
   integer, intent(in)    ::  max_site_site_num_nghbrs, &
                              site_site_num_nghbrs(*), &
                              site_site_list(*), &
                              site_site_ucell_transind(*), &
                              site_site_list_offset(*)
   double precision, intent(in) :: factor 
   double precision, intent(inout) :: virial(3,3)
   double precision, intent(inout) :: mp_dh_ene, &
                                      dh_mp_ene, &
                                      ch_dh_ene, &
                                      dh_ch_ene, &
                                      dh_dh_ene

   include "interact_type.fh"

   integer :: buf_count, neighbor_buf(max_site_site_num_nghbrs)
   integer site1,num_site_list,self_flag,use_cutoff, &
           signum,overall_signum,interaction_signum, &
           do_mpole_dherm,do_cherm_dherm,do_dherm_dherm

   double precision :: delx(max_site_site_num_nghbrs), &
                       dely(max_site_site_num_nghbrs), &
                       delz(max_site_site_num_nghbrs), &
                       delr2(max_site_site_num_nghbrs)

   self_flag = 0
   use_cutoff = 0  !here  we are either emulating recip sum or correcting it

   if ( (tot_num_mpole_coeffs > 0) .and. (tot_num_herm_Dprims > 0) )then
      do_mpole_dherm = 1
   else
      do_mpole_dherm = 0
   endif !( (tot_num_mpole_coeffs > 0) .and. (tot_num_herm_Dprims > 0) )

   if ( (tot_num_herm_Cprims > 0) .and. (tot_num_herm_Dprims > 0) )then
      do_cherm_dherm = 1
   else
      do_cherm_dherm = 0
   endif !( (tot_num_herm_Cprims > 0) .and. (tot_num_herm_Cprims > 0) )

   if ( tot_num_herm_Dprims > 0 )then
      do_dherm_dherm = 1
   else
      do_dherm_dherm = 0
   endif

   ! overall signum (are we accumulating interactions or correcting them?)
   if ( correcting_energies == 1 )then
      overall_signum = -1
   else
      overall_signum = 1
   endif

   do site1 = 1,num_sites-1
      call FH_SITE_SITE_LIST_filter(  &
                               site1, &
                               site_site_num_nghbrs, &
                               site_site_list_offset, &
                               site_site_list, &
                               site_site_ucell_transind, &
                               map_site_crd, &
                               ucell_translate, &
                               extent_of_compact_hermites, &
                               pbc, &
                               use_cutoff, &
                               buf_count, &
                               neighbor_buf, &
                               delx,dely,delz,delr2 )
      num_site_list = buf_count

      ! first the mpoles with diffuse herms
      if ( do_mpole_dherm == 1 )then
         interaction_signum = -1 !e.g. nuclear-electron
         signum = overall_signum * interaction_signum
         call FH_COULOMB_SITE_SITE_mp_herm( &
                site1, &
                num_site_list, &
                mpole_level_for_site, &
                mpole_coeff_off_for_site, &
                hermite_level_for_Dprim, &
                herm_coeff_offset_of_Dprim, &
                num_herm_Dprim_for_site, &
                off_herm_Dprim_for_site, &
                gauss_extent_for_Dprim, &
                hermite_expon_for_Dprim, &
                use_cutoff, &
                signum, &
                neighbor_buf,  &
                delx, &
                dely, &
                delz, &
                delr2, &
                Global_multipole_coeff, &
                Global_multipole_field, &
                Global_Dhermite_coeff, &
                Global_Dhermite_field, &
                mp_dh_ene, &
                site_frc,factor,virial )
         ! next the diffuse hermite - mpole
         call FH_COULOMB_SITE_SITE_herm_mp( &
                site1, &
                num_site_list, &
                mpole_level_for_site, &
                mpole_coeff_off_for_site, &
                hermite_level_for_Dprim, &
                herm_coeff_offset_of_Dprim, &
                num_herm_dprim_for_site, &
                off_herm_Dprim_for_site, &
                gauss_extent_for_Dprim, &
                hermite_expon_for_Dprim, &
                use_cutoff, &
                signum, &
                neighbor_buf, &
                delx, &
                dely, &
                delz, &
                delr2, &
                Global_multipole_coeff, &
                Global_multipole_field, &
                Global_Dhermite_coeff, &
                Global_Dhermite_field, &
                dh_mp_ene, &
                site_frc,factor,virial )
      endif !( do_mpole_dherm == 1 )then

      ! next the compact herms with diffuse herms
      if ( do_cherm_dherm == 1 )then
         interaction_signum = 1 !e.g. electron-electron
         signum = overall_signum * interaction_signum
         call FH_COULOMB_SITE_SITE_herm_herm( &
                site1, &
                num_site_list, &
                hermite_level_for_Cprim, &
                herm_coeff_offset_of_Cprim, &
                num_herm_Cprim_for_site, &
                off_herm_Cprim_for_site, &
                gauss_extent_for_Cprim, &
                hermite_expon_for_Cprim, &
                hermite_level_for_Dprim, &
                herm_coeff_offset_of_Dprim, &
                num_herm_Dprim_for_site, &
                off_herm_Dprim_for_site, &
                gauss_extent_for_Dprim, &
                hermite_expon_for_Dprim, &
                use_cutoff, &
                signum, &
                neighbor_buf, &
                factor, &
                delx, &
                dely, &
                delz, &
                delr2, &
                Global_Chermite_coeff, &
                Global_Chermite_field, &
                Global_Dhermite_coeff, &
                Global_Dhermite_field, &
                ch_dh_ene, &
                site_frc,virial)
         call FH_COULOMB_SITE_SITE_herm_herm( &
                site1, &
                num_site_list, &
                hermite_level_for_Dprim, &
                herm_coeff_offset_of_Dprim, &
                num_herm_Dprim_for_site, &
                off_herm_Dprim_for_site, &
                gauss_extent_for_Dprim, &
                hermite_expon_for_Dprim, &
                hermite_level_for_Cprim, &
                herm_coeff_offset_of_Cprim, &
                num_herm_Cprim_for_site, &
                off_herm_Cprim_for_site, &
                gauss_extent_for_Cprim, &
                hermite_expon_for_Cprim, &
                use_cutoff, &
                signum, &
                neighbor_buf, &
                factor, &
                delx, &
                dely, &
                delz, &
                delr2, &
                Global_Dhermite_coeff, &
                Global_Dhermite_field, &
                Global_Chermite_coeff, &
                Global_Chermite_field, &
                dh_ch_ene, &
                site_frc,virial)
      endif !( do_cherm_dherm == 1 )then

      ! finally the diffuse herms with diffuse herms
      if ( do_dherm_dherm == 1 )then
         interaction_signum = 1 !e.g. electron-electron
         signum = overall_signum * interaction_signum
         call FH_COULOMB_SITE_SITE_herm_herm( &
                site1, &
                num_site_list, &
                hermite_level_for_Dprim, &
                herm_coeff_offset_of_Dprim, &
                num_herm_Dprim_for_site, &
                off_herm_Dprim_for_site, &
                gauss_extent_for_Dprim, &
                hermite_expon_for_Dprim, &
                hermite_level_for_Dprim, &
                herm_coeff_offset_of_Dprim, &
                num_herm_Dprim_for_site, &
                off_herm_Dprim_for_site, &
                gauss_extent_for_Dprim, &
                hermite_expon_for_Dprim, &
                use_cutoff, &
                signum, &
                neighbor_buf, &
                factor, &
                delx, &
                dely, &
                delz, &
                delr2, &
                Global_Dhermite_coeff, &
                Global_Dhermite_field, &
                Global_Dhermite_coeff, &
                Global_Dhermite_field, &
                dh_dh_ene, &
                site_frc,virial)
      endif !( do_dherm_dherm == 1 )then

   enddo !site1 = 1,num_sites-1

end subroutine FH_COULOMB_SITE_SITE_CD_DD
!--------------------------------------------------------------
subroutine FH_COULOMB_SITE_SITE_self_CD_DD( &
                                     correcting_energies, &
                                     mp_dh_ene, &
                                     dh_mp_ene, &
                                     ch_dh_ene, &
                                     dh_ch_ene, &
                                     dh_dh_ene,factor,virial)

   use sites, only : num_sites, &
                     site_frc
   use hermite, only : tot_num_mpole_coeffs, &
                       Global_multipole_coeff, &
                       Global_multipole_field, &
                       mpole_level_for_site, &
                       mpole_coeff_off_for_site, &
                       ! compact hermite
                       tot_num_herm_Cprims, &
                       Global_Chermite_coeff, &
                       Global_Chermite_field, &
                       gauss_extent_for_Cprim, &
                       hermite_level_for_Cprim, &
                       num_herm_Cprim_for_site, &
                       off_herm_Cprim_for_site, &
                       herm_coeff_offset_of_Cprim, &
                       hermite_expon_for_Cprim, &
                       ! diffuse hermite
                       tot_num_herm_Dprims, &
                       Global_Dhermite_coeff, &
                       Global_Dhermite_field, &
                       gauss_extent_for_Dprim, &
                       hermite_level_for_Dprim, &
                       num_herm_Dprim_for_site, &
                       off_herm_Dprim_for_site, &
                       herm_coeff_offset_of_Dprim, &
                       hermite_expon_for_Dprim

   implicit none

   integer, intent(in) :: correcting_energies
   double precision, intent(in) :: factor 
   double precision, intent(inout) :: virial(3,3)
   double precision, intent(inout) :: mp_dh_ene, &
                                      dh_mp_ene, &
                                      ch_dh_ene, &
                                      dh_ch_ene, &
                                      dh_dh_ene

   integer :: buf_count, neighbor_buf(1)
   integer :: site1,num_site_list, &
              signum,overall_signum,interaction_signum, &
              do_mpole_dherm,do_cherm_dherm,do_dherm_dherm, &
              use_cutoff
   double precision :: delx(1), &
                       dely(1), &
                       delz(1), &
                       delr2(1)

   if ( (tot_num_mpole_coeffs > 0) .and. (tot_num_herm_Dprims > 0) )then
      do_mpole_dherm = 1
   else
      do_mpole_dherm = 0
   endif !( (tot_num_mpole_coeffs > 0) .and. (tot_num_herm_Dprims > 0) )

   if ( (tot_num_herm_Cprims > 0) .and. (tot_num_herm_Dprims > 0) )then
      do_cherm_dherm = 1
   else
      do_cherm_dherm = 0
   endif !( (tot_num_herm_Cprims > 0) .and. (tot_num_herm_Cprims > 0) )

   if ( tot_num_herm_Dprims > 0 )then
      do_dherm_dherm = 1
   else
      do_dherm_dherm = 0
   endif

   ! overall signum (are we accumulating interactions or correcting them?)
   if ( correcting_energies == 1 )then
      overall_signum = -1
   else
      overall_signum = 1
   endif

   use_cutoff = 0
   do site1 = 1,num_sites
      num_site_list = 1
      delx(1) = 0.d0
      dely(1) = 0.d0
      delz(1) = 0.d0
      delr2(1) = 0.d0
      neighbor_buf(1) = site1

      ! first the mpoles with diffuse herms
      if ( do_mpole_dherm == 1 )then
         interaction_signum = -1 !e.g. nuclear-electron
         signum = overall_signum * interaction_signum
         call FH_COULOMB_SITE_SITE_mp_herm( &
                site1, &
                num_site_list, &
                mpole_level_for_site, &
                mpole_coeff_off_for_site, &
                hermite_level_for_Dprim, &
                herm_coeff_offset_of_Dprim, &
                num_herm_Dprim_for_site, &
                off_herm_Dprim_for_site, &
                gauss_extent_for_Dprim, &
                hermite_expon_for_Dprim, &
                use_cutoff, &
                signum, &
                neighbor_buf,  &
                delx, &
                dely, &
                delz, &
                delr2, &
                Global_multipole_coeff, &
                Global_multipole_field, &
                Global_Dhermite_coeff, &
                Global_Dhermite_field, &
                mp_dh_ene, &
                site_frc,factor,virial )
      endif !( do_mpole_dherm == 1 )then

      ! next the compact herms with diffuse herms
      if ( do_cherm_dherm == 1 )then
         interaction_signum = 1 !e.g. electron-electron
         signum = overall_signum * interaction_signum
         call FH_COULOMB_SITE_SITE_herm_herm( &
                site1, &
                num_site_list, &
                hermite_level_for_Cprim, &
                herm_coeff_offset_of_Cprim, &
                num_herm_Cprim_for_site, &
                off_herm_Cprim_for_site, &
                gauss_extent_for_Cprim, &
                hermite_expon_for_Cprim, &
                hermite_level_for_Dprim, &
                herm_coeff_offset_of_Dprim, &
                num_herm_Dprim_for_site, &
                off_herm_Dprim_for_site, &
                gauss_extent_for_Dprim, &
                hermite_expon_for_Dprim, &
                use_cutoff, &
                signum, &
                neighbor_buf, &
                factor, &
                delx, &
                dely, &
                delz, &
                delr2, &
                Global_Chermite_coeff, &
                Global_Chermite_field, &
                Global_Dhermite_coeff, &
                Global_Dhermite_field, &
                ch_dh_ene, &
                site_frc,virial)
      endif !( do_cherm_dherm == 1 )then

      ! finally the diffuse herm self energy
      if ( do_dherm_dherm == 1 )then
         interaction_signum = 1 !e.g. electron-electron
         signum = overall_signum * interaction_signum
         call FH_COULOMB_SITE_SITE_self_herm( &
                                   site1, &
                                   signum, &
                                   hermite_level_for_Dprim, &
                                   herm_coeff_offset_of_Dprim, &
                                   num_herm_Dprim_for_site, &
                                   off_herm_Dprim_for_site, &
                                   hermite_expon_for_Dprim, &
                                   Global_Dhermite_coeff, &
                                   Global_Dhermite_field, &
                                   dh_dh_ene, &
                                   site_frc,factor,virial )
      endif !( do_dherm_dherm == 1 )then
   enddo !site1 = 1,num_sites

   mp_dh_ene = 0.5d0*mp_dh_ene
   dh_mp_ene = mp_dh_ene
   ch_dh_ene = 0.5d0*ch_dh_ene
   dh_ch_ene = ch_dh_ene
   dh_dh_ene = 0.5d0*dh_dh_ene
end subroutine FH_COULOMB_SITE_SITE_self_CD_DD
!--------------------------------------------------------------
!--------------------------------------------------------------
end module coulomb_site_site
!--------------------------------------------------------------

!-----------------------------------------------------------------
!-----------------------------------------------------------------
! NON-MODULE (LOCAL) SUBROUTINES
! THESE HANDLE ALL ARGUMENTS THROUGH ARGUMENT LIST
!-----------------------------------------------------------------

!-----------------------------------------------------------------
subroutine FH_COULOMB_SITE_SITE_mp_mp( &
                                   site1, &
                                   num_site_list, &
                                   mpole1_level, &
                                   mpole1_offset, &
                                   mpole2_level, &
                                   mpole2_offset, &
                                   extent_of_compact_hermites, &
                                   CD_split_expon, &
                                   full_interactions, &
                                   use_cutoff, &
                                   self_flag, &
                                   signum, &
                                   site_list, &
                                   dx,dy,dz,dr2, &
                                   GMpole1_coeff, &
                                   GMpole1_field, &
                                   GMpole2_coeff, &
                                   GMpole2_field, &
                                   energy, &
                                   site_frc,factor,virial )

   implicit none

   integer,intent(in) :: site1,num_site_list
   integer, intent(in) :: mpole1_level(*),mpole1_offset(*)
   integer, intent(in) :: mpole2_level(*),mpole2_offset(*)
   double precision,intent(in) :: extent_of_compact_hermites, &
                                  CD_split_expon
   integer,intent(in) :: full_interactions, &
                         use_cutoff, &
                         self_flag, &
                         signum, &
                         site_list(*)
   double precision, intent(in) :: factor 
   double precision,intent(in) :: dx(*),dy(*),dz(*),dr2(*), &
                                  GMpole1_coeff(*), GMpole2_coeff(*)
   double precision,intent(inout) :: GMpole1_field(*), &
                                     GMpole2_field(*), &
                                     energy,site_frc(3,*)
   double precision, intent(inout) :: virial(3,3)
   include "direct_pointers.fh"
   include "interact_type.fh"
   include 'scale.fh'

   double precision R(NCACHE*(MAXRLEV+1)**4)
   double precision expon2(NCACHE,0:MAXMPLEV)
   double precision delx(NCACHE,0:MAXMPLEV),dely(NCACHE,0:MAXMPLEV), &
                    delz(NCACHE,0:MAXMPLEV),delr2(NCACHE,0:MAXMPLEV), &
                    delr2inv(NCACHE,0:MAXMPLEV)
   double precision field1(NCACHE*MAXMPG), &
                    field2(NCACHE*MAXMP), &
                    hermite_coeff2(NCACHE,MAXMP,0:MAXMPLEV)
   integer numlist(0:MAXMPLEV),site_pointer(NCACHE,0:MAXMPLEV), &
           field_offset2(NCACHE,0:MAXMPLEV)

   integer :: level1,order1,field_order1,offset1, &
              site2,level2,order2,offset2, &
              top_level,nlist,j,n,proceed,mcmur_dav_init_recur_rule
   double precision :: expon1,sum_extent
   double precision :: ewald_expon,ewald_extent


   ewald_expon = 0.5d0*CD_split_expon
   ! the extent goes as the inverse sqrt of exponent
   ewald_extent = sqrt(2.d0)*extent_of_compact_hermites
   if ( use_cutoff == 1 )then
      mcmur_dav_init_recur_rule = MPOLE_COULOMB_MPOLE_MINUS_HERM
   else
      if ( full_interactions == 1 )then
         mcmur_dav_init_recur_rule = MPOLE_COULOMB_MPOLE
      else
         mcmur_dav_init_recur_rule = MPOLE_COULOMB_HERM
      endif
   endif !( use_cutoff == 1 )

   numlist = 0
   level1 = mpole1_level(site1)
   if ( level1 < 0 )return
   order1 = hermite_order_from_level(level1)
   field_order1 = hermite_order_from_level(level1+1)
   offset1 = mpole1_offset(site1)
   expon1 = big_value
   do n = 1,num_site_list
      site2 = site_list(n)
      level2 = mpole2_level(site2)
      sum_extent = ewald_extent
      proceed = 0
      if ( level2 >= 0 )then
         if ( (use_cutoff == 0) .or. &
              ((use_cutoff == 1) .and. &
                (dr2(n) < sum_extent*sum_extent)) )then
              proceed = 1
         endif
      endif !( level2 >= 0 )then
      if ( proceed == 1 )then
         order2 = hermite_order_from_level(level2)
         top_level = level1 + level2 + 1
         numlist(level2) = numlist(level2) + 1
         nlist = numlist(level2)
         expon2(nlist,level2) = ewald_expon
         site_pointer(nlist,level2) = site2
         delx(nlist,level2) = dx(n)
         dely(nlist,level2) = dy(n)
         delz(nlist,level2) = dz(n)
         delr2(nlist,level2) = dr2(n)
         delr2inv(nlist,level2) = 1.d0 / dr2(n)
         offset2 = mpole2_offset(site2)
         field_offset2(nlist,level2) = offset2
         do j = 1,order2
            hermite_coeff2(nlist,j,level2) = GMpole2_coeff(offset2+j)
         enddo
         if ( nlist == NCACHE )then !process full list
            call FH_MCMUR_DAV_recur(nlist,top_level, &
                    mcmur_dav_init_recur_rule, &
                    expon1,expon2(1,level2),R, &
                    delx(1,level2),dely(1,level2), &
                    delz(1,level2),delr2(1,level2),delr2inv(1,level2))
            call FH_MCMUR_DAV_fill_fields(nlist,ncache,top_level, &
                    field_order1,order1,order2, &
                    field1,field2,GMpole1_coeff(offset1+1), &
                    hermite_coeff2(1,1,level2),R)
            call FH_MCMUR_DAV_update_Global_fields( &
                    nlist,ncache,signum,factor,self_flag, &
                    order1,offset1,order2,field_offset2(1,level2), &
                    GMpole1_field,GMpole2_field,field1,field2)
            call FH_MCMUR_DAV_ene_frc(nlist,ncache,signum,factor, &
                    order1,site1,site_pointer(1,level2), &
                    delx(1,level2),dely(1,level2),delz(1,level2),&
                    GMpole1_coeff(offset1+1),field1,energy,site_frc,virial)
            numlist(level2) = 0
         endif ! nlist == NCACHE
      endif !(proceed ==1)
   enddo ! n = 1,num_site_list
   ! process remaining lists
   do level2 = 0,MAXMPLEV
      nlist = numlist(level2)
      top_level = level1 + level2 + 1
      order2 = hermite_order_from_level(level2)
      if ( nlist > 0 )then
         call FH_MCMUR_DAV_recur(nlist,top_level, &
                 mcmur_dav_init_recur_rule, &
                 expon1,expon2(1,level2),R, &
                 delx(1,level2),dely(1,level2), &
                 delz(1,level2),delr2(1,level2),delr2inv(1,level2))
         call FH_MCMUR_DAV_fill_fields(nlist,ncache,top_level, &
                 field_order1,order1,order2, &
                 field1,field2,GMpole1_coeff(offset1+1), &
                 hermite_coeff2(1,1,level2),R)
         call FH_MCMUR_DAV_update_Global_fields( &
                 nlist,ncache,signum,factor,self_flag, &
                 order1,offset1,order2,field_offset2(1,level2), &
                 GMpole1_field,GMpole2_field,field1,field2)
         call FH_MCMUR_DAV_ene_frc(nlist,ncache,signum,factor, &
                 order1,site1,site_pointer(1,level2), &
                 delx(1,level2),dely(1,level2),delz(1,level2),&
                 GMpole1_coeff(offset1+1),field1,energy,site_frc,virial)
      endif ! nlist > 0
   enddo !level2 = 0,MAXMPLEV
end subroutine FH_COULOMB_SITE_SITE_mp_mp
!-----------------------------------------------------------------------
subroutine FH_COULOMB_SITE_SITE_mp_herm( &
                                   site1, &
                                   num_site_list, &
                                   mpole_level, &
                                   mpole_offset, &
                                   hermite_level, &
                                   hermite_offset, &
                                   num_prims_for_site, &
                                   prim_offset_for_site, &
                                   prim_extent, &
                                   prim_expon, &
                                   use_cutoff, &
                                   signum, &
                                   site_list, &
                                   dx, &
                                   dy, &
                                   dz, &
                                   dr2, &
                                   GMpole_coeff, &
                                   GMpole_field, &
                                   GHerm_coeff, &
                                   GHerm_field, &
                                   energy, &
                                   site_frc,factor,virial )

   implicit none

   integer,intent(in) :: site1, &
                         num_site_list, &
                         mpole_level(*), &
                         mpole_offset(*), &
                         hermite_level(*), &
                         hermite_offset(*), &
                         num_prims_for_site(*), &
                         prim_offset_for_site(*)
   double precision,intent(in) :: prim_extent(*), &
                                  prim_expon(*)
   integer,intent(in) :: use_cutoff, &
                         signum, &
                         site_list(*)
   double precision,intent(in) :: factor
   double precision,intent(in) :: dx(*),dy(*),dz(*),dr2(*), &
                                  GMpole_coeff(*), &
                                  GHerm_coeff(*)
   double precision,intent(inout) :: GMpole_field(*), &
                                     GHerm_field(*), &
                                     energy, &
                                     site_frc(3,*)
   double precision, intent(inout) :: virial(3,3)
   include "direct_pointers.fh"
   include "interact_type.fh"
   include 'scale.fh'

   double precision R(NCACHE*(MAXRLEV+1)**4)
   double precision expon2(NCACHE,0:MAXMPLEV)
   double precision delx(NCACHE,0:MAXMPLEV),dely(NCACHE,0:MAXMPLEV), &
                    delz(NCACHE,0:MAXMPLEV),delr2(NCACHE,0:MAXMPLEV), &
                    delr2inv(NCACHE,0:MAXMPLEV)
   double precision field1(NCACHE*MAXMPG), &
                    field2(NCACHE*MAXMP), &
                    hermite_coeff2(NCACHE,MAXMP,0:MAXMPLEV)
   integer numlist(0:MAXMPLEV),site_pointer(NCACHE,0:MAXMPLEV), &
           field_offset2(NCACHE,0:MAXMPLEV)

   integer :: level1,order1,field_order1,offset1, &
              site2,level2,order2,offset2, &
              top_level,nlist,j,n,num2,off2,kp2,np2, &
              proceed,mcmur_dav_init_recur_rule,self_flag
   double precision :: expon1,sum_extent

   if ( use_cutoff == 1 )then
      mcmur_dav_init_recur_rule = MPOLE_COULOMB_HERM_MINUS_MPOLE
   else
      mcmur_dav_init_recur_rule = MPOLE_COULOMB_HERM
   endif !( use_cutoff == 1 )

   level1 = mpole_level(site1)
   if ( level1 < 0 )return
   numlist = 0
   self_flag = 0
   order1 = hermite_order_from_level(level1)
   field_order1 = hermite_order_from_level(level1+1)
   offset1 = mpole_offset(site1)
   expon1 = big_value
   do n = 1,num_site_list
      site2 = site_list(n)
      num2 = num_prims_for_site(site2)
      off2 = prim_offset_for_site(site2)
      do kp2 = 1,num2
         np2 = off2 + kp2
         level2 = hermite_level(np2)
         sum_extent = prim_extent(np2)
         proceed = 0
         if ( level2 >= 0 )then
            if ( (use_cutoff == 0) .or. &
                 ((use_cutoff == 1) .and. &
                   (dr2(n) < sum_extent*sum_extent)) )then
                 proceed = 1
            endif
         endif !( level2 >= 0 )
         if ( (proceed ==1) )then
            order2 = hermite_order_from_level(level2)
            top_level = level1 + level2 + 1
            numlist(level2) = numlist(level2) + 1
            nlist = numlist(level2)
            expon2(nlist,level2) = prim_expon(np2)
            site_pointer(nlist,level2) = site2
            delx(nlist,level2) = dx(n)
            dely(nlist,level2) = dy(n)
            delz(nlist,level2) = dz(n)
            delr2(nlist,level2) = dr2(n)
            delr2inv(nlist,level2) = 1.d0 / dr2(n)
            offset2 = hermite_offset(np2)
            field_offset2(nlist,level2) = offset2
            do j = 1,order2
               hermite_coeff2(nlist,j,level2) = GHerm_coeff(offset2+j)
            enddo
            if ( nlist == NCACHE )then !process full list
               call FH_MCMUR_DAV_recur(nlist,top_level, &
                       mcmur_dav_init_recur_rule, &
                       expon1,expon2(1,level2),R, &
                       delx(1,level2),dely(1,level2), &
                       delz(1,level2),delr2(1,level2), &
                       delr2inv(1,level2))
               call FH_MCMUR_DAV_fill_fields(nlist,ncache,top_level, &
                       field_order1,order1,order2, &
                       field1,field2,GMpole_coeff(offset1+1), &
                       hermite_coeff2(1,1,level2),R)
               call FH_MCMUR_DAV_update_Global_fields( &
                       nlist,ncache,signum,factor,self_flag, &
                       order1,offset1,order2,field_offset2(1,level2), &
                       GMpole_field,GHerm_field,field1,field2)
               call FH_MCMUR_DAV_ene_frc(nlist,ncache,signum,factor, &
                       order1,site1,site_pointer(1,level2), &
                       delx(1,level2),dely(1,level2),delz(1,level2),&
                       GMpole_coeff(offset1+1),field1,energy,site_frc,virial)
               numlist(level2) = 0
            endif ! nlist == NCACHE
         endif !(proceed ==1)
      enddo !np = 1,num
   enddo ! n = 1,num_site_list
   ! process remaining lists
   do level2 = 0,MAXMPLEV
      nlist = numlist(level2)
      top_level = level1 + level2 + 1
      order2 = hermite_order_from_level(level2)
      if ( nlist > 0 )then
         call FH_MCMUR_DAV_recur(nlist,top_level, &
                 mcmur_dav_init_recur_rule, &
                 expon1,expon2(1,level2),R, &
                 delx(1,level2),dely(1,level2), &
                 delz(1,level2),delr2(1,level2),delr2inv(1,level2))
         call FH_MCMUR_DAV_fill_fields(nlist,ncache,top_level, &
                 field_order1,order1,order2, &
                 field1,field2,GMpole_coeff(offset1+1), &
                 hermite_coeff2(1,1,level2),R)
         call FH_MCMUR_DAV_update_Global_fields( &
                 nlist,ncache,signum,factor,self_flag, &
                 order1,offset1,order2,field_offset2(1,level2), &
                 GMpole_field,GHerm_field,field1,field2)
         call FH_MCMUR_DAV_ene_frc(nlist,ncache,signum,factor, &
                 order1,site1,site_pointer(1,level2), &
                 delx(1,level2),dely(1,level2),delz(1,level2),&
                 GMpole_coeff(offset1+1),field1,energy,site_frc,virial)
      endif ! nlist > 0
   enddo !level2 = 0,MAXMPLEV
end subroutine FH_COULOMB_SITE_SITE_mp_herm
!-----------------------------------------------------------------------
subroutine FH_COULOMB_SITE_SITE_herm_mp( &
                                   site1, &
                                   num_site_list, &
                                   mpole_level, &
                                   mpole_offset, &
                                   hermite_level, &
                                   hermite_offset, &
                                   num_prims_for_site, &
                                   prim_offset_for_site, &
                                   prim_extent, &
                                   prim_expon, &
                                   use_cutoff, &
                                   signum, &
                                   site_list, &
                                   dx, &
                                   dy, &
                                   dz, &
                                   dr2, &
                                   GMpole_coeff, &
                                   GMpole_field, &
                                   GHerm_coeff, &
                                   GHerm_field, &
                                   energy, &
                                   site_frc,factor,virial )

   implicit none

   integer,intent(in) :: site1, &
                         num_site_list, &
                         mpole_level(*), &
                         mpole_offset(*), &
                         hermite_level(*), &
                         hermite_offset(*), &
                         num_prims_for_site(*), &
                         prim_offset_for_site(*)
   double precision,intent(in) :: prim_extent(*), &
                                  prim_expon(*)
   integer,intent(in) :: use_cutoff, &
                         signum, &
                         site_list(*)
   double precision,intent(in) :: factor
   double precision,intent(in) :: dx(*),dy(*),dz(*),dr2(*), &
                                  GMpole_coeff(*), &
                                  GHerm_coeff(*)
   double precision,intent(inout) :: GMpole_field(*), &
                                     GHerm_field(*), &
                                     energy, &
                                     site_frc(3,*)
   double precision, intent(inout) :: virial(3,3)
   include "direct_pointers.fh"
   include "interact_type.fh"
   include 'scale.fh'

   double precision R(NCACHE*(MAXRLEV+1)**4)
   double precision expon2(NCACHE,0:MAXMPLEV)
   double precision delx(NCACHE,0:MAXMPLEV),dely(NCACHE,0:MAXMPLEV), &
                    delz(NCACHE,0:MAXMPLEV),delr2(NCACHE,0:MAXMPLEV), &
                    delr2inv(NCACHE,0:MAXMPLEV)
   double precision field1(NCACHE*MAXMPG), &
                    field2(NCACHE*MAXMP), &
                    hermite_coeff2(NCACHE,MAXMP,0:MAXMPLEV)
   integer numlist(0:MAXMPLEV),site_pointer(NCACHE,0:MAXMPLEV), &
           field_offset2(NCACHE,0:MAXMPLEV)

   integer :: level1,order1,field_order1,offset1, &
              site2,level2,order2,offset2, &
              top_level,nlist,j,n,num1,off1,kp1,np1, &
              proceed1,proceed2,mcmur_dav_init_recur_rule,self_flag
   double precision :: expon1,sum_extent

   if ( use_cutoff == 1 )then
      mcmur_dav_init_recur_rule = HERM_MINUS_MPOLE_COULOMB_MPOLE
   else
      mcmur_dav_init_recur_rule = HERM_COULOMB_MPOLE
   endif !( use_cutoff == 1 )

   num1 = num_prims_for_site(site1)
   off1 = prim_offset_for_site(site1)
   self_flag = 0
   do kp1 = 1,num1
      np1 = off1 + kp1
      proceed1 = 1 !all cases pass this hurdle 
                   !(compact and diffuse handled separately now)
      if ( (proceed1 == 1) )then
         level1 = hermite_level(np1)
         order1 = hermite_order_from_level(level1)
         field_order1 = hermite_order_from_level(level1+1)
         offset1 = hermite_offset(np1)
         expon1 = prim_expon(np1)
         numlist = 0
         do n = 1,num_site_list
            site2 = site_list(n)
            level2 = mpole_level(site2)
            sum_extent = prim_extent(np1)
            proceed2 = 0
            if ( level2 >= 0 )then
               if ( (use_cutoff == 0) .or. &
                    ((use_cutoff == 1) .and. &
                      (dr2(n) < sum_extent*sum_extent)) )then
                    proceed2 = 1
               endif
            endif !( level2 >= 0 )then
            if ( (proceed2 ==1) )then
               order2 = hermite_order_from_level(level2)
               top_level = level1 + level2 + 1
               numlist(level2) = numlist(level2) + 1
               nlist = numlist(level2)
               expon2(nlist,level2) = big_value
               site_pointer(nlist,level2) = site2
               delx(nlist,level2) = dx(n)
               dely(nlist,level2) = dy(n)
               delz(nlist,level2) = dz(n)
               delr2(nlist,level2) = dr2(n)
               delr2inv(nlist,level2) = 1.d0 / dr2(n)
               offset2 = mpole_offset(site2)
               field_offset2(nlist,level2) = offset2
               do j = 1,order2
                  hermite_coeff2(nlist,j,level2) =  &
                  GMpole_coeff(offset2+j)
               enddo
               if ( nlist == NCACHE )then !process full list
                  call FH_MCMUR_DAV_recur(nlist, &
                          top_level,mcmur_dav_init_recur_rule, &
                          expon1,expon2(1,level2),R, &
                          delx(1,level2),dely(1,level2), &
                          delz(1,level2),delr2(1,level2), &
                          delr2inv(1,level2))
                  call FH_MCMUR_DAV_fill_fields(nlist,ncache, &
                          top_level, &
                          field_order1,order1,order2, &
                          field1,field2,GHerm_coeff(offset1+1), &
                          hermite_coeff2(1,1,level2),R)
                  call FH_MCMUR_DAV_update_Global_fields( &
                          nlist,ncache,signum,factor,self_flag, &
                          order1,offset1,order2, &
                          field_offset2(1,level2), &
                          GHerm_field,GMpole_field,field1,field2)
                  call FH_MCMUR_DAV_ene_frc(nlist,ncache,signum, &
                          factor,order1,site1, &
                          site_pointer(1,level2), &
                          delx(1,level2),dely(1,level2),delz(1,level2),&
                          GHerm_coeff(offset1+1),field1, &
                          energy,site_frc,virial)
                  numlist(level2) = 0
               endif ! nlist == NCACHE
            endif !(proceed2 ==1)
         enddo !n = 1,num_site_list
         do level2 = 0,MAXMPLEV
            nlist = numlist(level2)
            top_level = level1 + level2 + 1
            order2 = hermite_order_from_level(level2)
            if ( nlist > 0 )then
               call FH_MCMUR_DAV_recur(nlist, &
                       top_level,mcmur_dav_init_recur_rule, &
                       expon1,expon2(1,level2),R, &
                       delx(1,level2),dely(1,level2), &
                       delz(1,level2),delr2(1,level2), &
                       delr2inv(1,level2))
               call FH_MCMUR_DAV_fill_fields(nlist,ncache,top_level, &
                       field_order1,order1,order2, &
                       field1,field2,GHerm_coeff(offset1+1), &
                       hermite_coeff2(1,1,level2),R)
               call FH_MCMUR_DAV_update_Global_fields( &
                       nlist,ncache,signum,factor,self_flag, &
                       order1,offset1,order2,field_offset2(1,level2), &
                       GHerm_field,GMpole_field,field1,field2)
               call FH_MCMUR_DAV_ene_frc(nlist,ncache,signum,factor, &
                       order1,site1,site_pointer(1,level2), &
                       delx(1,level2),dely(1,level2),delz(1,level2),&
                       GHerm_coeff(offset1+1),field1,energy,site_frc,virial)
            endif ! nlist > 0
         enddo !level2 = 0,MAXMPLEV
      endif !(proceed1 == 1)
   enddo !kp1 = 1,num1
   ! process remaining lists
end subroutine FH_COULOMB_SITE_SITE_herm_mp
!-----------------------------------------------------------------------
subroutine FH_COULOMB_SITE_SITE_herm_herm( &
                                   site1, &
                                   num_site_list, &
                                   hermite1_level, &
                                   hermite1_offset, &
                                   num_prim1_for_site, &
                                   prim1_offset_for_site, &
                                   prim1_extent, &
                                   prim1_expon, &
                                   hermite2_level, &
                                   hermite2_offset, &
                                   num_prim2_for_site, &
                                   prim2_offset_for_site, &
                                   prim2_extent, &
                                   prim2_expon, &
                                   use_cutoff, &
                                   signum, &
                                   site_list, &
                                   factor, &
                                   dx, &
                                   dy, &
                                   dz, &
                                   dr2, &
                                   GHerm1_coeff, &
                                   GHerm1_field, &
                                   GHerm2_coeff, &
                                   GHerm2_field, &
                                   energy, &
                                   site_frc,virial )

   implicit none

   integer,intent(in) :: site1,num_site_list, &
                         hermite1_level(*), &
                         hermite1_offset(*), &
                         num_prim1_for_site(*), &
                         prim1_offset_for_site(*)
   double precision,intent(in) :: prim1_extent(*), &
                                  prim1_expon(*)
   integer,intent(in) :: hermite2_level(*), &
                         hermite2_offset(*), &
                         num_prim2_for_site(*), &
                         prim2_offset_for_site(*)
   double precision,intent(in) :: prim2_extent(*), &
                                  prim2_expon(*)
   integer,intent(in) :: use_cutoff, &
                         signum, &
                         site_list(*)
   double precision,intent(in) :: factor,dx(*),dy(*),dz(*),dr2(*), &
                                  GHerm1_coeff(*),GHerm2_coeff(*)
   double precision,intent(inout) :: GHerm1_field(*), &
                                     GHerm2_field(*), &
                                     energy,site_frc(3,*)
   double precision, intent(inout) :: virial(3,3)
   include "direct_pointers.fh"
   include "interact_type.fh"

   double precision R(NCACHE*(MAXRLEV+1)**4)
   double precision expon2(NCACHE,0:MAXMPLEV)
   double precision delx(NCACHE,0:MAXMPLEV),dely(NCACHE,0:MAXMPLEV), &
                    delz(NCACHE,0:MAXMPLEV),delr2(NCACHE,0:MAXMPLEV), &
                    delr2inv(NCACHE,0:MAXMPLEV)
   double precision field1(NCACHE*MAXMPG), &
                    field2(NCACHE*MAXMP), &
                    hermite_coeff2(NCACHE,MAXMP,0:MAXMPLEV)
   integer numlist(0:MAXMPLEV),site_pointer(NCACHE,0:MAXMPLEV), &
           field_offset2(NCACHE,0:MAXMPLEV)

   integer :: level1,order1,field_order1,offset1, &
              site2,level2,order2,offset2, &
              top_level,nlist,j,n,off1,num1,off2,num2,kp1,kp2,np1,np2, &
              proceed1,proceed2,mcmur_dav_init_recur_rule,self_flag
   double precision :: expon1,sum_extent

   if ( use_cutoff == 1 )then
      mcmur_dav_init_recur_rule = H_COULOMB_H_MINUS_MP_COULOMB_MP
   else
      mcmur_dav_init_recur_rule = HERM_COULOMB_HERM
   endif !( use_cutoff == 1 )

   num1 = num_prim1_for_site(site1)
   off1 = prim1_offset_for_site(site1)
   self_flag = 0
   do kp1 = 1,num1
      np1 = off1 + kp1
      proceed1 = 1 !filtering here happened above
      if ( (proceed1 == 1) )then
         level1 = hermite1_level(np1)
         order1 = hermite_order_from_level(level1)
         field_order1 = hermite_order_from_level(level1+1)
         offset1 = hermite1_offset(np1)
         expon1 = prim1_expon(np1)
         numlist = 0
         do n = 1,num_site_list
            site2 = site_list(n)
            num2 = num_prim2_for_site(site2)
            off2 = prim2_offset_for_site(site2)
            do kp2 = 1,num2
               np2 = off2 + kp2
               level2 = hermite2_level(np2)
               sum_extent = prim1_extent(np1) + prim2_extent(np2)
               proceed2 = 0
               if ( level1 >= 0 .and. level2 >= 0 )then
                  if ( (use_cutoff == 0) .or. &
                       ((use_cutoff == 1) .and. &
                         (dr2(n) < sum_extent*sum_extent)) )then
                       proceed2 = 1
                  endif
               endif !( level1 >= 0 .and. level2 >= 0 )then
               if ( (proceed2 == 1) )then
                  order2 = hermite_order_from_level(level2)
                  top_level = level1 + level2 + 1
                  numlist(level2) = numlist(level2) + 1
                  nlist = numlist(level2)
                  expon2(nlist,level2) = prim2_expon(np2)
                  site_pointer(nlist,level2) = site2
                  delx(nlist,level2) = dx(n)
                  dely(nlist,level2) = dy(n)
                  delz(nlist,level2) = dz(n)
                  delr2(nlist,level2) = dr2(n)
                  delr2inv(nlist,level2) = 1.d0 / dr2(n)
                  offset2 = hermite2_offset(np2)
                  field_offset2(nlist,level2) = offset2
                  do j = 1,order2
                     hermite_coeff2(nlist,j,level2) =  &
                     GHerm2_coeff(offset2+j)
                  enddo
                  if ( nlist == NCACHE )then !process full list
                     call FH_MCMUR_DAV_recur( &
                          nlist,top_level, &
                          mcmur_dav_init_recur_rule, &
                          expon1,expon2(1,level2),R, &
                          delx(1,level2),dely(1,level2), &
                          delz(1,level2),delr2(1,level2), &
                          delr2inv(1,level2))
                     call FH_MCMUR_DAV_fill_fields(nlist,ncache,top_level, &
                          field_order1,order1,order2, &
                          field1,field2,GHerm1_coeff(offset1+1), &
                          hermite_coeff2(1,1,level2),R)
                     call FH_MCMUR_DAV_update_Global_fields( &
                          nlist,ncache,signum,factor,self_flag, &
                          order1,offset1,order2,field_offset2(1,level2), &
                          GHerm1_field,GHerm2_field,field1,field2)
                     call FH_MCMUR_DAV_ene_frc(nlist,ncache,signum,factor, &
                          order1,site1,site_pointer(1,level2), &
                          delx(1,level2),dely(1,level2),delz(1,level2),&
                          GHerm1_coeff(offset1+1),field1,energy,site_frc,virial)
                     numlist(level2) = 0
                  endif ! nlist == NCACHE
               endif ! (proceed2 == 1)
            enddo !kp2 = 1,num2
         enddo !n = 1,num_site_list
         ! process remaining lists
         do level2 = 0,MAXMPLEV
            nlist = numlist(level2)
            top_level = level1 + level2 + 1
            order2 = hermite_order_from_level(level2)
            if ( nlist > 0 )then
               call FH_MCMUR_DAV_recur(nlist,top_level, &
                       mcmur_dav_init_recur_rule, &
                       expon1,expon2(1,level2),R, &
                       delx(1,level2),dely(1,level2), &
                       delz(1,level2),delr2(1,level2),delr2inv(1,level2))
               call FH_MCMUR_DAV_fill_fields(nlist,ncache,top_level, &
                       field_order1,order1,order2, &
                       field1,field2,GHerm1_coeff(offset1+1), &
                       hermite_coeff2(1,1,level2),R)
               call FH_MCMUR_DAV_update_Global_fields( &
                       nlist,ncache,signum,factor,self_flag, &
                       order1,offset1,order2,field_offset2(1,level2), &
                       GHerm1_field,GHerm2_field,field1,field2)
               call FH_MCMUR_DAV_ene_frc(nlist,ncache,signum,factor, &
                       order1,site1,site_pointer(1,level2), &
                       delx(1,level2),dely(1,level2),delz(1,level2),&
                       GHerm1_coeff(offset1+1),field1,energy,site_frc,virial)
            endif ! nlist > 0
         enddo !level2 = 0,MAXMPLEV
      endif ! (proceed1 == 1) )
   enddo ! kp1 = 1,num1
end subroutine FH_COULOMB_SITE_SITE_herm_herm
!-----------------------------------------------------------------------
subroutine FH_COULOMB_SITE_SITE_self_herm( &
                                   site1, &
                                   signum, &
                                   hermite_level, &
                                   hermite_offset, &
                                   num_prims_for_site, &
                                   prim_offset_for_site, &
                                   prim_expon, &
                                   GHerm_coeff, &
                                   GHerm_field, &
                                   ene_self, &
                                   site_frc,factor,virial )

   implicit none

   integer,intent(in) :: site1, &
                         signum, &
                         hermite_level(*), &
                         hermite_offset(*), &
                         num_prims_for_site(*), &
                         prim_offset_for_site(*)
   double precision,intent(in) :: factor
   double precision,intent(in) :: prim_expon(*)
   double precision,intent(in) :: GHerm_coeff(*)
   double precision,intent(inout) :: GHerm_field(*), &
                                     ene_self, &
                                     site_frc(3,*)
   double precision, intent(inout) :: virial(3,3)
   include "direct_pointers.fh"
   include "interact_type.fh"

   double precision R(NCACHE*(MAXRLEV+1)**4)
   double precision expon2(NCACHE,0:MAXMPLEV)
   double precision delx(NCACHE,0:MAXMPLEV),dely(NCACHE,0:MAXMPLEV), &
                    delz(NCACHE,0:MAXMPLEV),delr2(NCACHE,0:MAXMPLEV), &
                    delr2inv(NCACHE,0:MAXMPLEV)
   double precision field1(NCACHE*MAXMPG), &
                    field2(NCACHE*MAXMP), &
                    hermite_coeff2(NCACHE,MAXMP,0:MAXMPLEV)
   integer numlist(0:MAXMPLEV),site_pointer(NCACHE,0:MAXMPLEV), &
           field_offset2(NCACHE,0:MAXMPLEV)

   integer :: level1,order1,field_order1,offset1, &
              site2,level2,order2,offset2, &
              top_level,nlist,j,n,off1,num1,off2,num2,kp1,kp2,np1,np2, &
              proceed1,proceed2,self_flag,mcmur_dav_init_recur_rule
   double precision :: expon1,energy1,energy2

   ! setup interaction type
   mcmur_dav_init_recur_rule = HERM_COULOMB_HERM

   ! first the self terms
   num1 = num_prims_for_site(site1)
   off1 = prim_offset_for_site(site1)
   self_flag = 1 ! to prevent double filling of global hermite field
   energy1 = 0.d0
   energy2 = 0.d0

   do kp1 = 1,num1
      np1 = off1 + kp1
      numlist = 0
      level1 = hermite_level(np1)
      if ( level1 >= 0 )then
         order1 = hermite_order_from_level(level1)
         level2 = level1
         order2 = order1
         field_order1 = hermite_order_from_level(level1+1)
         offset1 = hermite_offset(np1)
         offset2 = offset1
         expon1 = prim_expon(np1)
         top_level = level1 + level2 + 1
         numlist(level2) = 1
         nlist = numlist(level2)
         expon2(nlist,level2) = prim_expon(np1)
         site_pointer(nlist,level2) = site1
         delx(nlist,level2) = 0.d0
         dely(nlist,level2) = 0.d0
         delz(nlist,level2) = 0.d0
         delr2(nlist,level2) = 0.d0
         delr2inv(nlist,level2) = 0.d0
         offset2 = hermite_offset(np1)
         field_offset2(nlist,level2) = offset2
         do j = 1,order2
            hermite_coeff2(nlist,j,level2) = GHerm_coeff(offset2+j)
         enddo
         ! process remaining lists
         do level2 = 0,MAXMPLEV
            nlist = numlist(level2)
            top_level = level1 + level2 + 1
            order2 = hermite_order_from_level(level2)
            if ( nlist > 0 )then
               call FH_MCMUR_DAV_recur(nlist,top_level, &
                       mcmur_dav_init_recur_rule, &
                       expon1,expon2(1,level2),R, &
                       delx(1,level2),dely(1,level2), &
                       delz(1,level2),delr2(1,level2),delr2inv(1,level2))
               call FH_MCMUR_DAV_fill_fields(nlist,ncache,top_level, &
                       field_order1,order1,order2, &
                       field1,field2,GHerm_coeff(offset1+1), &
                       hermite_coeff2(1,1,level2),R)
               call FH_MCMUR_DAV_update_Global_fields( &
                       nlist,ncache,signum,factor,self_flag, &
                       order1,offset1,order2,field_offset2(1,level2), &
                       GHerm_field,GHerm_field,field1,field2)
               call FH_MCMUR_DAV_ene_frc(nlist,ncache,signum,factor, &
                       order1,site1,site_pointer(1,level2), &
                       delx(1,level2),dely(1,level2),delz(1,level2),&
                       GHerm_coeff(offset1+1),field1,energy1,site_frc,virial)
            endif ! nlist > 0
         enddo !level2 = 0,MAXMPLEV
      endif !( level1 >= 0 )then
   enddo ! kp1 = 1,num1
   self_flag = 0 ! Global field now different offsets
   ! now the remaining terms
   do kp1 = 1,num1 - 1
      np1 = off1 + kp1
      level1 = hermite_level(np1)
      order1 = hermite_order_from_level(level1)
      field_order1 = hermite_order_from_level(level1+1)
      offset1 = hermite_offset(np1)
      expon1 = prim_expon(np1)
      numlist = 0
      do kp2 = kp1+1,num1
         np2 = off1 + kp2
         level2 = hermite_level(np2)
         proceed2 = 0
         ! note the pair_filter  is ALL_but_COMPACT_with_COMPACT
         if ( (level1 >= 0) .and. (level2 >= 0 ) )then 
              proceed2 = 1
         endif
         if ( (proceed2 == 1) )then
            order2 = hermite_order_from_level(level2)
            top_level = level1 + level2 + 1
            numlist(level2) = numlist(level2) + 1
            nlist = numlist(level2)
            expon2(nlist,level2) = prim_expon(np2)
            site_pointer(nlist,level2) = site1
            delx(nlist,level2) = 0.d0
            dely(nlist,level2) = 0.d0
            delz(nlist,level2) = 0.d0
            delr2(nlist,level2) = 0.d0
            delr2inv(nlist,level2) = 0.d0
            offset2 = hermite_offset(np2)
            field_offset2(nlist,level2) = offset2
            do j = 1,order2
               hermite_coeff2(nlist,j,level2) = GHerm_coeff(offset2+j)
            enddo
            if ( nlist == NCACHE )then !process full list
               call FH_MCMUR_DAV_recur(nlist,top_level, &
                    mcmur_dav_init_recur_rule, &
                    expon1,expon2(1,level2),R, &
                    delx(1,level2),dely(1,level2), &
                    delz(1,level2),delr2(1,level2),delr2inv(1,level2))
               call FH_MCMUR_DAV_fill_fields(nlist,ncache,top_level, &
                    field_order1,order1,order2, &
                    field1,field2,GHerm_coeff(offset1+1), &
                    hermite_coeff2(1,1,level2),R)
               call FH_MCMUR_DAV_update_Global_fields( &
                    nlist,ncache,signum,factor,self_flag, &
                    order1,offset1,order2,field_offset2(1,level2), &
                    GHerm_field,GHerm_field,field1,field2)
               call FH_MCMUR_DAV_ene_frc(nlist,ncache,signum,factor, &
                    order1,site1,site_pointer(1,level2), &
                    delx(1,level2),dely(1,level2),delz(1,level2),&
                    GHerm_coeff(offset1+1),field1,energy2,site_frc,virial)
               numlist(level2) = 0
            endif ! nlist == NCACHE
         endif ! (proceed2 == 1)
      enddo !kp2 = kp1+1,num1
      ! process remaining lists
      do level2 = 0,MAXMPLEV
         nlist = numlist(level2)
         top_level = level1 + level2 + 1
         order2 = hermite_order_from_level(level2)
         if ( nlist > 0 )then
            call FH_MCMUR_DAV_recur(nlist,top_level, &
                    mcmur_dav_init_recur_rule, &
                    expon1,expon2(1,level2),R, &
                    delx(1,level2),dely(1,level2), &
                    delz(1,level2),delr2(1,level2),delr2inv(1,level2))
            call FH_MCMUR_DAV_fill_fields(nlist,ncache,top_level, &
                    field_order1,order1,order2, &
                    field1,field2,GHerm_coeff(offset1+1), &
                    hermite_coeff2(1,1,level2),R)
            call FH_MCMUR_DAV_update_Global_fields( &
                    nlist,ncache,signum,factor,self_flag, &
                    order1,offset1,order2,field_offset2(1,level2), &
                    GHerm_field,GHerm_field,field1,field2)
            call FH_MCMUR_DAV_ene_frc(nlist,ncache,signum,factor, &
                    order1,site1,site_pointer(1,level2), &
                    delx(1,level2),dely(1,level2),delz(1,level2),&
                    GHerm_coeff(offset1+1),field1,energy2,site_frc,virial)
         endif ! nlist > 0
      enddo !level2 = 0,MAXMPLEV
   enddo !kp1 = 1,num1

   ene_self = ene_self + energy1 + 2.d0 * energy2
end subroutine FH_COULOMB_SITE_SITE_self_herm
!-----------------------------------------------------------------------
