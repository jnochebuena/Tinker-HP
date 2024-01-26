module exchange_site_site
   implicit none
   private
          
   public FH_EXCHANGE_SITE_SITE_CC, &
          !FH_EXCHANGE_SITE_SITE_self_CC, &
          FH_EXCHANGE_SITE_SITE_CD_DD, &
          FH_EXCHANGE_SITE_SITE_self_CD_DD
          
   contains
!--------------------------------------------------------------------

!-----------------------------------------------------------------
!-----------------------------------------------------------------
! EVALUATION SUBROUTINES
!-----------------------------------------------------------------
!-----------------------------------------------------------------

!-------------------------------------------------------------------------
subroutine FH_EXCHANGE_SITE_SITE_CC( &
                                     correcting_energies, &
                                     use_cutoff, &
                                     full_interactions, &
                                     max_site_site_num_nghbrs, &
                                     site_site_num_nghbrs, &
                                     site_site_list_offset, &
                                     site_site_list, &
                                     site_site_ucell_transind, &
                                     ch_ch_ene,factor,exch_cutoff,virial)

   use user, only : pbc, &
                    CD_split_expon, &
                    extent_of_compact_hermites
   use sites, only : num_sites, &
                     map_site_crd, &
                     site_frc
   use unit_cell, only : ucell_translate
   use hermite, only : tot_num_herm_Cprims, &
                       Global_Chermite_coeff, &
                       Global_Chermite_field_EXCH, &
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
   double precision, intent(in) :: factor, exch_cutoff
   double precision, intent(inout) :: virial(3,3)
   double precision, intent(inout) :: ch_ch_ene

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
   double precision :: cutoff_dist
                       

   integer jj

   self_flag = 0
   
      ! NO MPOLES, ONLY HERMITES BECAUSE OF OVERLAPS
      do_mpole_mpole = 0
      do_mpole_cherm = 0
      do_mpole_sumch = 0

   if ( tot_num_herm_Cprims > 0 )then
      if ( (full_interactions == 1) .or. (use_cutoff == 1) )then
         do_cherm_cherm = 1
      else
         do_cherm_cherm = 0
      endif
   else
      do_cherm_cherm = 0
   endif ! ( tot_num_herm_Cprims > 0 )

         do_sumch_sumch = 0

   ! overall signum (are we accumulating interactions or correcting them?)
   if ( correcting_energies == 1 )then
      overall_signum = -1
   else
      overall_signum = 1
   endif

   ! TRANSFORM CUTOFF TO BOHR AND SQUARE TO COMPARE TO DISTANCE
   cutoff_dist = (exch_cutoff/0.529177249d0)*(exch_cutoff/0.529177249d0)

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

      ! only compact herms with compact herms
      if ( do_cherm_cherm == 1 )then
         interaction_signum = 1 !e.g. electron-electron
         signum = overall_signum * interaction_signum
         call FH_EXCHANGE_SITE_SITE_herm_herm( &
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
                cutoff_dist, &
                delx, &
                dely, &
                delz, &
                delr2, &
                Global_Chermite_coeff, &
                Global_Chermite_field_EXCH, &
                Global_Chermite_coeff, &
                Global_Chermite_field_EXCH, &
                ch_ch_ene, &
                site_frc,virial)
      endif !( do_cherm_cherm == 1 )then
   enddo !site1 = 1,num_sites-1

end subroutine FH_EXCHANGE_SITE_SITE_CC
!--------------------------------------------------------------
!!!subroutine FH_EXCHANGE_SITE_SITE_self_CC( &
!!!OJO: NO SITE_SITE_self for exchange because no mpoles!
!--------------------------------------------------------------
!--------------------------------------------------------------
subroutine FH_EXCHANGE_SITE_SITE_CD_DD( &
                                     correcting_energies, &
                                     max_site_site_num_nghbrs, &
                                     site_site_num_nghbrs, &
                                     site_site_list_offset, &
                                     site_site_list, &
                                     site_site_ucell_transind, &
                                     ch_dh_ene, &
                                     dh_ch_ene, &
                                     dh_dh_ene,factor,exch_cutoff,virial)

   use user, only : pbc, &
                    extent_of_compact_hermites
   use sites, only : num_sites, &
                     map_site_crd, &
                     site_frc
   use unit_cell, only : ucell_translate
   use hermite, only :  &
                       ! compact hermite
                       tot_num_herm_Cprims, &
                       Global_Chermite_coeff, &
                       Global_Chermite_field_EXCH, &
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
   double precision, intent(in) :: factor,exch_cutoff
   double precision, intent(inout) :: virial(3,3)
   double precision, intent(inout) :: ch_dh_ene, &
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
   double precision :: cutoff_dist

   self_flag = 0
   
   use_cutoff = 0  !here  we are either emulating recip sum or correcting it


      do_mpole_dherm = 0

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

   ! TRANSFORM CUTOFF TO BOHR AND SQUARE TO COMPARE TO DISTANCE
   cutoff_dist = (exch_cutoff/0.529177249d0)*(exch_cutoff/0.529177249d0)

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

      ! 1st compact herms with diffuse herms
      if ( do_cherm_dherm == 1 )then
         interaction_signum = 1 !e.g. electron-electron
         signum = overall_signum * interaction_signum
         call FH_EXCHANGE_SITE_SITE_herm_herm( &
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
                cutoff_dist, &
                delx, &
                dely, &
                delz, &
                delr2, &
                Global_Chermite_coeff, &
                Global_Chermite_field_EXCH, &
                Global_Dhermite_coeff, &
                Global_Dhermite_field, &
                ch_dh_ene, &
                site_frc,virial)
         call FH_EXCHANGE_SITE_SITE_herm_herm( &
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
                cutoff_dist, &
                delx, &
                dely, &
                delz, &
                delr2, &
                Global_Dhermite_coeff, &
                Global_Dhermite_field, &
                Global_Chermite_coeff, &
                Global_Chermite_field_EXCH, &
                dh_ch_ene, &
                site_frc,virial)
      endif !( do_cherm_dherm == 1 )then

      ! finally the diffuse herms with diffuse herms
      if ( do_dherm_dherm == 1 )then
         interaction_signum = 1 !e.g. electron-electron
         signum = overall_signum * interaction_signum
         call FH_EXCHANGE_SITE_SITE_herm_herm( &
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
                cutoff_dist, &
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

end subroutine FH_EXCHANGE_SITE_SITE_CD_DD
!--------------------------------------------------------------
subroutine FH_EXCHANGE_SITE_SITE_self_CD_DD( &
                                     correcting_energies, &
                                     ch_dh_ene, &
                                     dh_ch_ene, &
                                     dh_dh_ene,factor,exch_cutoff,virial)

   use sites, only : num_sites, &
                     site_frc
                       ! compact hermite
   use hermite, only : tot_num_herm_Cprims, &
                       Global_Chermite_coeff, &
                       Global_Chermite_field_EXCH, &
                       gauss_extent_for_Cprim, &
                       hermite_level_for_Cprim, &
                       num_herm_Cprim_for_site, &
                       off_herm_Cprim_for_site, &
                       herm_coeff_offset_of_Cprim, &
                       hermite_expon_for_Cprim, &
                       ! diffuse hermite
                       tot_num_herm_Dprims, &
                       Global_Dhermite_coeff, &
                       Global_Dhermite_field_EXCH, &
                       gauss_extent_for_Dprim, &
                       hermite_level_for_Dprim, &
                       num_herm_Dprim_for_site, &
                       off_herm_Dprim_for_site, &
                       herm_coeff_offset_of_Dprim, &
                       hermite_expon_for_Dprim

   implicit none

   integer, intent(in) :: correcting_energies
   double precision, intent(in) :: factor, exch_cutoff
   double precision, intent(inout) :: virial(3,3)
   double precision, intent(inout) :: ch_dh_ene, &
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
   double precision :: cutoff_dist


      do_mpole_dherm = 0

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
   
   ! TRANSFORM CUTOFF TO BOHR AND SQUARE TO COMPARE TO DISTANCE
   cutoff_dist = (exch_cutoff/0.529177249d0)*(exch_cutoff/0.529177249d0)

   do site1 = 1,num_sites
      num_site_list = 1
      delx(1) = 0.d0
      dely(1) = 0.d0
      delz(1) = 0.d0
      delr2(1) = 0.d0
      neighbor_buf(1) = site1

      ! first the compact herms with diffuse herms
      if ( do_cherm_dherm == 1 )then
         interaction_signum = 1 !e.g. electron-electron
         signum = overall_signum * interaction_signum
         call FH_EXCHANGE_SITE_SITE_herm_herm( &
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
                cutoff_dist, &
                delx, &
                dely, &
                delz, &
                delr2, &
                Global_Chermite_coeff, &
                Global_Chermite_field_EXCH, &
                Global_Dhermite_coeff, &
                Global_Dhermite_field_EXCH, &
                ch_dh_ene, &
                site_frc,virial)
      endif !( do_cherm_dherm == 1 )then

      ! finally the diffuse herm self energy
      if ( do_dherm_dherm == 1 )then
         interaction_signum = 1 !e.g. electron-electron
         signum = overall_signum * interaction_signum
         call FH_EXCHANGE_SITE_SITE_self_herm( &
                                   site1, &
                                   signum, &
                                   hermite_level_for_Dprim, &
                                   herm_coeff_offset_of_Dprim, &
                                   num_herm_Dprim_for_site, &
                                   off_herm_Dprim_for_site, &
                                   hermite_expon_for_Dprim, &
                                   Global_Dhermite_coeff, &
                                   Global_Dhermite_field_EXCH, &
                                   dh_dh_ene, &
                                   site_frc,factor,virial )
      endif !( do_dherm_dherm == 1 )then
   enddo !site1 = 1,num_sites

   ch_dh_ene = 0.5d0*ch_dh_ene
   dh_ch_ene = ch_dh_ene
   dh_dh_ene = 0.5d0*dh_dh_ene
end subroutine FH_EXCHANGE_SITE_SITE_self_CD_DD
!--------------------------------------------------------------
!--------------------------------------------------------------
end module exchange_site_site
!--------------------------------------------------------------

!-----------------------------------------------------------------
!-----------------------------------------------------------------
! NON-MODULE (LOCAL) SUBROUTINES
! THESE HANDLE ALL ARGUMENTS THROUGH ARGUMENT LIST
!-----------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine FH_EXCHANGE_SITE_SITE_herm_herm( &
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
                                   cutoff_dist, &
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
   double precision,intent(in) :: factor,cutoff_dist,&
                                  dx(*),dy(*),dz(*),dr2(*), &
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
              proceed1,proceed2,mcmur_dav_init_recur_rule,self_flag!, &
              !tot_ints, eval_ints
   double precision :: expon1,sum_extent

   
      mcmur_dav_init_recur_rule = HERM_OVERLAP_HERM

   num1 = num_prim1_for_site(site1)
   off1 = prim1_offset_for_site(site1)
   self_flag = 0
   !tot_ints = 0
   !eval_ints = 0
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
!  CHECK HERE IF DISTANCE IS LESS THAN 4.0 ANGSTROM, IF IT IS THEN
!  DO INTEGRALS, IF NOT, SKIP!
          !tot_ints = tot_ints + 1
          if (dr2(n) < cutoff_dist) then
            !eval_ints = eval_ints + 1
            site2 = site_list(n)
            num2 = num_prim2_for_site(site2)
            off2 = prim2_offset_for_site(site2)
            do kp2 = 1,num2
               np2 = off2 + kp2
               level2 = hermite2_level(np2)
               proceed2 = 1 !DBG
!DBG           sum_extent = prim1_extent(np1) + prim2_extent(np2)
!DBG           proceed2 = 0
!DBG           if ( level1 >= 0 .and. level2 >= 0 )then
!DBG              if ( (use_cutoff == 0) .or. &
!DBG                   ((use_cutoff == 1) .and. &
!DBG                     (dr2(n) < sum_extent*sum_extent)) )then
!DBG                   proceed2 = 1
!DBG              endif
!DBG           endif !( level1 >= 0 .and. level2 >= 0 )then
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
           endif ! TO CHECK CENTER DISTANCE TO DO INTEGRALS OR NOT
         enddo !n = 1,num_site_list
         ! process remaining lists
         do level2 = 0,MAXMPLEV
            nlist = numlist(level2)
            top_level = level1 + level2 + 1
            order2 = hermite_order_from_level(level2)
            if ( nlist > 0 )then
!  CHECK HERE IF DISTANCE IS LESS THAN 3.5 ANGSTROM IF IT IS LESS THEN
!  DO INTEGRALS, IF NOT, SKIP!
              !tot_ints = tot_ints + 1
              if (delr2(1,level2) < cutoff_dist) then
               !eval_ints = eval_ints + 1
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
              endif ! TO CHECK CENTER DISTANCE TO DO INTEGRALS OR NOT
            endif ! nlist > 0
         enddo !level2 = 0,MAXMPLEV
      endif ! (proceed1 == 1) )
   enddo ! kp1 = 1,num1

   return

end subroutine FH_EXCHANGE_SITE_SITE_herm_herm
!-----------------------------------------------------------------------
subroutine FH_EXCHANGE_SITE_SITE_self_herm( &
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
   double precision, intent(inout) :: virial(3,3)
   double precision,intent(inout) :: GHerm_field(*), &
                                     ene_self, &
                                     site_frc(3,*)
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
   mcmur_dav_init_recur_rule = HERM_OVERLAP_HERM

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
end subroutine FH_EXCHANGE_SITE_SITE_self_herm
!-----------------------------------------------------------------------
