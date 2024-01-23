!--------------------------------------------------------------------
! these routines build the lists of pairs of sites that interact in
! various circumstances. 4 lists are possible:
! 1) the unmasked_pairs list---used for finite clusters to calculate all
!    intermolecular (unmasked) pairs directly i.e.
!    with no decomposition or splitting--this is triggered in user input 
!    by setting CD_split_expon to 0
! 2) the cutoff_direct list---used in finite clusters or periodic boundary
!    conditions to calculate all unmasked pairs out ot a cutoff using a
!    damped (erfc) coulomb interaction
! 3) the emulate_recip list. this is used in finite clusters to emulate
!    the calculation done by reciprocal space techniques -note that
!    masked as well as unmasked pairs are calculated as undamped 
!    coulomb interactions - however the self interactions calculated in
!    reciprocal space are emulated separately in a self energy routine
!    this latter is necessary due to need to treat self terms specially
! 4) the correct_recip list --this is used to subtract the masked terms 
!    in reciprocal space or emulated as above---used for periodic as well
!    as non-periodic simulations
!--------------------------------------------------------------------
module site_site_list

   implicit none

   private

   ! lists for NO_SPLIT calculations
   integer,save :: max_unmasked_pairs_list, &
                   num_unmasked_pairs_list, &
                   max_unmasked_pairs_num_nghbrs
   integer,save, allocatable :: unmasked_pairs_list(:), &
                                unmasked_pairs_ucell_transind(:), &
                                unmasked_pairs_num_nghbrs(:), &
                                unmasked_pairs_list_offset(:)

   ! lists for direct sum compact-compact interactions for 
   ! SPLIT_DIRECT & SPLIT_RECIP calculations
   integer,save :: max_cutoff_direct_list, &
                   num_cutoff_direct_list, &
                   max_cutoff_direct_num_nghbrs
   integer,save, allocatable :: cutoff_direct_list(:), &
                                cutoff_direct_ucell_transind(:), &
                                cutoff_direct_num_nghbrs(:), &
                                cutoff_direct_list_offset(:)

   ! lists for direct sum compact-diffuse and diffuse-diffuse interactions for 
   ! SPLIT_DIRECT calculations
   integer,save :: max_emulate_recip_list, &
                   num_emulate_recip_list, &
                   max_emulate_recip_num_nghbrs
   integer,save, allocatable :: emulate_recip_list(:), &
                                emulate_recip_ucell_transind(:), &
                                emulate_recip_num_nghbrs(:), &
                                emulate_recip_list_offset(:)

   ! lists for direct sum correction to reciprocal interactions for 
   ! SPLIT_RECIP calculations
   integer,save :: max_correct_recip_list, &
                   num_correct_recip_list, &
                   max_correct_recip_num_nghbrs
   integer,save, allocatable :: correct_recip_list(:), &
                                correct_recip_ucell_transind(:), &
                                correct_recip_num_nghbrs(:), &
                                correct_recip_list_offset(:)

   ! flag variables with an initial value within a subroutine create a sort
   ! re-entrancy issue (simple) - you can't re-do the setup if it these are
   ! buried in the code..

   integer, save        :: first_time

   public FH_SITE_SITE_LIST_init, &
          FH_SITE_SITE_LIST_eval, &
          FH_SITE_SITE_LIST_deallocate, &
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
          
   contains
!--------------------------------------------------------------------

!-----------------------------------------------------------------
!-----------------------------------------------------------------
! SETUP SUBROUTINES
!-----------------------------------------------------------------

subroutine FH_SITE_SITE_LIST_init

  implicit none

  first_time = 1

  return

end subroutine FH_SITE_SITE_LIST_init

!-----------------------------------------------------------------

!-----------------------------------------------------------------
!-----------------------------------------------------------------
!-----------------------------------------------------------------
! EVALUATION SUBROUTINES
!-----------------------------------------------------------------
!-----------------------------------------------------------------

!--------------------------------------------------------------------
subroutine FH_SITE_SITE_LIST_eval(out_lun)

   use sites,only : num_sites, &
                    num_sites_masked_by_site, &
                    masked_sitelist_offset_of_site, &
                    masked_sitelist, &
                    frac_site_crd, &
                    map_site_crd
   use user, only : pbc, &
                    cut_extra, &
                    extent_of_compact_hermites, &
                    calculation_method, &
                    verbose
   use unit_cell, only : ucell, &
                         ucell_translate, &
                         extended_ucell, &
                         num_translates

   implicit none

! Formal arguments:

   integer, intent(in)  :: out_lun

! Local variables:

   include "interact_type.fh"
   include "masking_rules.fh"

   ! NOTE: You MUST call FH_SITE_SITE_LIST_init() at start of each cfgem_coul()
   !       and cfgem_exch() call to initialize first_time..

   integer :: siz_buffer
   integer :: masking_rule, &
              use_cutoff
   integer :: ier1,ier2,ier3,ier4

   ! siz_buffer is the largest num neighbors any site could have
   if ( pbc == 0 )then
      siz_buffer = num_sites
   elseif ( pbc == 1 .and. extended_ucell == 0 )then
      siz_buffer = num_sites
   elseif ( pbc == 1 .and. extended_ucell > 0 )then
      siz_buffer = num_translates*num_sites
   endif
   if ( calculation_method == NO_SPLIT )then !do all unmasked pairs
      masking_rule = SKIP_MASKED_PAIRS
      use_cutoff = 1
      if ( first_time == 1 )then
         call FH_SITE_SITE_LIST_count( &
                 num_sites, &
                 num_sites_masked_by_site, &
                 masked_sitelist_offset_of_site, &
                 masked_sitelist, &
                 masking_rule, &
                 use_cutoff, &
                 pbc, &
                 extended_ucell, &
                 num_translates, &
                 cut_extra, &
                 ucell, &
                 ucell_translate, &
                 extent_of_compact_hermites, &
                 frac_site_crd, &
                 map_site_crd, &
                 num_unmasked_pairs_list, &
                 out_lun)
         if ( verbose == 1 )  &
         write(out_lun,*)'num_unmasked_pairs_list = ',num_unmasked_pairs_list
         max_unmasked_pairs_list = 1.5*num_unmasked_pairs_list
         allocate(unmasked_pairs_list(max_unmasked_pairs_list),stat=ier1)
         allocate(unmasked_pairs_ucell_transind(max_unmasked_pairs_list), &
                  stat=ier2)
         allocate(unmasked_pairs_num_nghbrs(num_sites),stat=ier3)
         allocate(unmasked_pairs_list_offset(num_sites),stat=ier4)
         if ( ier1 /= 0 .or. ier2 /= 0 .or. ier3 /= 0 .or. ier4 /= 0 )then
            write(out_lun,*) &
               'FH_SITE_SITE_lists: trouble allocating unmasked_pairs'
            stop
         endif
      endif !( first_time == 1 )then
      call FH_SITE_SITE_LIST_build(& 
              num_sites, &
              num_sites_masked_by_site, &
              masked_sitelist_offset_of_site, &
              masked_sitelist, &
              masking_rule, &
              use_cutoff, &
              pbc, &
              extended_ucell, &
              num_translates, &
              cut_extra, &
              ucell, &
              ucell_translate, &
              extent_of_compact_hermites, &
              frac_site_crd, &
              map_site_crd, &
              siz_buffer, &
              max_unmasked_pairs_list, &
              num_unmasked_pairs_list, &
              unmasked_pairs_list, &
              unmasked_pairs_num_nghbrs, &
              unmasked_pairs_list_offset, &
              unmasked_pairs_ucell_transind, &
              max_unmasked_pairs_num_nghbrs, &
              out_lun)
   elseif ( (calculation_method == SPLIT_DIRECT) .or. &
            (calculation_method == SPLIT_RECIP) )then
      ! (I)  load the compact-compact direct
      masking_rule = SKIP_MASKED_PAIRS
      use_cutoff = 1
      if ( first_time == 1 )then
         call FH_SITE_SITE_LIST_count( &
                 num_sites, &
                 num_sites_masked_by_site, &
                 masked_sitelist_offset_of_site, &
                 masked_sitelist, &
                 masking_rule, &
                 use_cutoff, &
                 pbc, &
                 extended_ucell, &
                 num_translates, &
                 cut_extra, &
                 ucell, &
                 ucell_translate, &
                 extent_of_compact_hermites, &
                 frac_site_crd, &
                 map_site_crd, &
                 num_cutoff_direct_list, &
                 out_lun)
         if ( verbose == 1 )  &
         write(out_lun,*)'num_cutoff_direct_list = ',num_cutoff_direct_list
         max_cutoff_direct_list = 1.5*num_cutoff_direct_list
         allocate(cutoff_direct_list(max_cutoff_direct_list), &
                  stat=ier1)
         allocate(cutoff_direct_ucell_transind(max_cutoff_direct_list), &
                  stat=ier2)
         allocate(cutoff_direct_num_nghbrs(num_sites),stat=ier3)
         allocate(cutoff_direct_list_offset(num_sites),stat=ier4)
         if ( ier1 /= 0 .or. ier2 /= 0 .or. ier3 /= 0 .or. ier4 /= 0 )then
            write(out_lun,*) &
               'FH_SITE_SITE_lists: trouble allocating cutoff_direct'
            stop
         endif
      endif !( first_time == 1 )then
      call FH_SITE_SITE_LIST_build(& 
              num_sites, &
              num_sites_masked_by_site, &
              masked_sitelist_offset_of_site, &
              masked_sitelist, &
              masking_rule, &
              use_cutoff, &
              pbc, &
              extended_ucell, &
              num_translates, &
              cut_extra, &
              ucell, &
              ucell_translate, &
              extent_of_compact_hermites, &
              frac_site_crd, &
              map_site_crd, &
              siz_buffer, &
              max_cutoff_direct_list, &
              num_cutoff_direct_list, &
              cutoff_direct_list, &
              cutoff_direct_num_nghbrs, &
              cutoff_direct_list_offset, &
              cutoff_direct_ucell_transind, &
              max_cutoff_direct_num_nghbrs, &
              out_lun)
      ! (2) load the correct recip list: only masked pairs 
      masking_rule = ONLY_MASKED_PAIRS
      use_cutoff = 0
      if ( first_time == 1 )then
         call FH_SITE_SITE_LIST_count( &
                 num_sites, &
                 num_sites_masked_by_site, &
                 masked_sitelist_offset_of_site, &
                 masked_sitelist, &
                 masking_rule, &
                 use_cutoff, &
                 pbc, &
                 extended_ucell, &
                 num_translates, &
                 cut_extra, &
                 ucell, &
                 ucell_translate, &
                 extent_of_compact_hermites, &
                 frac_site_crd, &
                 map_site_crd, &
                 num_correct_recip_list, &
                 out_lun)
         if ( verbose == 1 )  &
         write(out_lun,*)'num_correct_recip_list = ',num_correct_recip_list
         max_correct_recip_list = 1.5*num_correct_recip_list
         allocate(correct_recip_list(max_correct_recip_list), &
                  stat=ier1)
         allocate(correct_recip_ucell_transind(max_correct_recip_list), &
                  stat=ier2)
         allocate(correct_recip_num_nghbrs(num_sites),stat=ier3)
         allocate(correct_recip_list_offset(num_sites),stat=ier4)
         if ( ier1 /= 0 .or. ier2 /= 0 .or. ier3 /= 0 .or. ier4 /= 0 )then
            write(out_lun,*) &
               'FH_SITE_SITE_lists: trouble allocating correct_recip'
            stop
         endif
      endif !( first_time == 1 )then
      call FH_SITE_SITE_LIST_build(& 
              num_sites, &
              num_sites_masked_by_site, &
              masked_sitelist_offset_of_site, &
              masked_sitelist, &
              masking_rule, &
              use_cutoff, &
              pbc, &
              extended_ucell, &
              num_translates, &
              cut_extra, &
              ucell, &
              ucell_translate, &
              extent_of_compact_hermites, &
              frac_site_crd, &
              map_site_crd, &
              siz_buffer, &
              max_correct_recip_list, &
              num_correct_recip_list, &
              correct_recip_list, &
              correct_recip_num_nghbrs, &
              correct_recip_list_offset, &
              correct_recip_ucell_transind, &
              max_correct_recip_num_nghbrs, &
              out_lun)
   endif !( calculation_method == NO_SPLIT )

   if ( calculation_method == SPLIT_DIRECT  )then
      ! load the list for direct space version of recip calculation
      ! this should get every pair except self interactions, as
      ! the recip calculation does
      masking_rule = INCLUDE_MASKED_PAIRS
      use_cutoff = 0
      if ( first_time == 1 )then
         call FH_SITE_SITE_LIST_count( &
                 num_sites, &
                 num_sites_masked_by_site, &
                 masked_sitelist_offset_of_site, &
                 masked_sitelist, &
                 masking_rule, &
                 use_cutoff, &
                 pbc, &
                 extended_ucell, &
                 num_translates, &
                 cut_extra, &
                 ucell, &
                 ucell_translate, &
                 extent_of_compact_hermites, &
                 frac_site_crd, &
                 map_site_crd, &
                 num_emulate_recip_list, &
                 out_lun)
         if ( verbose == 1 )  &
         write(out_lun,*)'num_emulate_recip_list = ',num_emulate_recip_list
         max_emulate_recip_list = 1.5*num_emulate_recip_list
         allocate(emulate_recip_list(max_emulate_recip_list), &
                  stat=ier1)
         allocate(emulate_recip_ucell_transind(max_emulate_recip_list), &
                  stat=ier2)
         allocate(emulate_recip_num_nghbrs(num_sites),stat=ier3)
         allocate(emulate_recip_list_offset(num_sites),stat=ier4)
         if ( ier1 /= 0 .or. ier2 /= 0 .or. ier3 /= 0 .or. ier4 /= 0 )then
            write(out_lun,*) &
               'FH_SITE_SITE_lists: trouble allocating emulate_recip'
            stop
         endif
      endif !( first_time == 1 )then
      call FH_SITE_SITE_LIST_build(& 
              num_sites, &
              num_sites_masked_by_site, &
              masked_sitelist_offset_of_site, &
              masked_sitelist, &
              masking_rule, &
              use_cutoff, &
              pbc, &
              extended_ucell, &
              num_translates, &
              cut_extra, &
              ucell, &
              ucell_translate, &
              extent_of_compact_hermites, &
              frac_site_crd, &
              map_site_crd, &
              siz_buffer, &
              max_emulate_recip_list, &
              num_emulate_recip_list, &
              emulate_recip_list, &
              emulate_recip_num_nghbrs, &
              emulate_recip_list_offset, &
              emulate_recip_ucell_transind, &
              max_emulate_recip_num_nghbrs, &
              out_lun)
   endif! ( calculation_method == SPLIT_DIRECT  )
   first_time = 0
end subroutine FH_SITE_SITE_LIST_eval
!--------------------------------------------------------------

!-----------------------------------------------------------------
!-----------------------------------------------------------------
! DEALLOCATE SUBROUTINES
!-----------------------------------------------------------------
!-----------------------------------------------------------------

!-----------------------------------------------------------------
subroutine FH_SITE_SITE_LIST_deallocate()

   implicit none

   ! lists for NO_SPLIT calculations
   if ( allocated(unmasked_pairs_list) )deallocate(unmasked_pairs_list)
   if ( allocated(unmasked_pairs_ucell_transind) ) &
             deallocate(unmasked_pairs_ucell_transind)
   if ( allocated(unmasked_pairs_num_nghbrs) ) &
             deallocate(unmasked_pairs_num_nghbrs)
   if ( allocated(unmasked_pairs_list_offset) ) &
             deallocate(unmasked_pairs_list_offset)

   ! lists for direct sum compact-compact interactions for 
   ! SPLIT_DIRECT & SPLIT_RECIP calculations
   if ( allocated(cutoff_direct_list) )deallocate(cutoff_direct_list)
   if ( allocated(cutoff_direct_ucell_transind) ) &
             deallocate(cutoff_direct_ucell_transind)
   if ( allocated(cutoff_direct_num_nghbrs) ) &
             deallocate(cutoff_direct_num_nghbrs)
   if ( allocated(cutoff_direct_list_offset) ) &
             deallocate(cutoff_direct_list_offset)

   ! lists for direct sum compact-diffuse and diffuse-diffuse interactions for 
   ! SPLIT_DIRECT calculations
   if ( allocated(emulate_recip_list) )deallocate(emulate_recip_list)
   if ( allocated(emulate_recip_ucell_transind) ) &
             deallocate(emulate_recip_ucell_transind)
   if ( allocated(emulate_recip_num_nghbrs) ) &
             deallocate(emulate_recip_num_nghbrs)
   if ( allocated(emulate_recip_list_offset) ) &
             deallocate(emulate_recip_list_offset)

   ! lists for direct sum compact-diffuse and diffuse-diffuse correction
   ! (to reciprocal) interactions for  SPLIT_RECIP calculations
   if ( allocated(correct_recip_list) )deallocate(correct_recip_list)
   if ( allocated(correct_recip_ucell_transind) ) &
             deallocate(correct_recip_ucell_transind)
   if ( allocated(correct_recip_num_nghbrs) ) &
             deallocate(correct_recip_num_nghbrs)
   if ( allocated(correct_recip_list_offset) ) &
             deallocate(correct_recip_list_offset)
end subroutine FH_SITE_SITE_LIST_deallocate
!-----------------------------------------------------------------
!-----------------------------------------------------------------
end module site_site_list
!-----------------------------------------------------------------
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!-----------------------------------------------------------------
! NON-MODULE (LOCAL) SUBROUTINES
! THESE HANDLE ALL ARGUMENTS THROUGH ARGUMENT LIST
!-----------------------------------------------------------------
!-----------------------------------------------------------------

!-----------------------------------------------------------------

subroutine FH_SITE_SITE_LIST_count( &
                    num_sites, &
                    num_sites_masked_by_site, &
                    masked_sitelist_offset_of_site, &
                    masked_sitelist, &
                    masking_rule, &
                    use_cutoff, &
                    pbc, &
                    extended_ucell, &
                    num_translates, &
                    cut_extra, &
                    ucell, &
                    ucell_translate, &
                    extent_of_compact_hermites, &
                    frac_site_crd, &
                    map_site_crd, &
                    num_list, &
                    out_lun)

   use user, only : do_overlap, exchange_cutoff

   implicit none

! Formal arguments:

   integer, intent(in)           :: num_sites
   integer, intent(in)           :: num_sites_masked_by_site(num_sites), &
                                    masked_sitelist_offset_of_site(num_sites), &
                                    masked_sitelist(*)
   integer, intent(in)           :: masking_rule, use_cutoff
   integer, intent(in)           :: pbc, extended_ucell, num_translates
   double precision, intent(in)  :: cut_extra
   double precision, intent(in)  :: ucell(3,3), ucell_translate(3,*), &
                                    extent_of_compact_hermites, &
                                    frac_site_crd(3,*), map_site_crd(3,*)
   integer,intent(out)           :: num_list

   integer, intent(in)           :: out_lun

! Local variables:

   include "masking_rules.fh"

   double precision, parameter  :: bohr = 0.52917721092d0
   double precision, parameter  :: bohr_per_angstrom = 1.d0 / bohr

   integer :: site1,site2,j,k,off,num,scratch_array(num_sites)
   integer :: num_per_site,start,proceed
   double precision :: cut,cut2,df(3)
   double precision :: xf,yf,zf,xm,ym,zm,dx0,dy0,dz0,dx,dy,dz,dr2,dr
   
   ! check sanity of parameter combinations
   if ( (use_cutoff == 1) .and. &
        (  (masking_rule==INCLUDE_MASKED_PAIRS) .or. &
           (masking_rule==ONLY_MASKED_PAIRS)  )  )then
       write(out_lun,*)'FH_SITE_SITE_LIST_count: ', &
                 'use_cutoff=1 incompatible with INCLUDE_MASKED_PAIRS', &
                 ' or ONLY_MASKED_PAIRS'
       stop
   endif
   ! clear the scratch array
   call UTIL_zero_int_array(scratch_array,num_sites)

   if ( use_cutoff == 1 )then
      if (do_overlap .ne. 0) then
        cut = exchange_cutoff * bohr_per_angstrom + cut_extra
        cut2 = cut*cut
      else
        cut = 2.d0*extent_of_compact_hermites + cut_extra
        cut2 = cut*cut
      end if
   endif
   num_list = 0
   do site1 = 1,num_sites
      off = masked_sitelist_offset_of_site(site1)
      num = num_sites_masked_by_site(site1)
      num_per_site = 0
      start = site1+1
      do j = 1,num
         k = masked_sitelist(off+j)
         scratch_array(k) = site1
      enddo
      xm = map_site_crd(1,site1)
      ym = map_site_crd(2,site1)
      zm = map_site_crd(3,site1)
      xf = frac_site_crd(1,site1)
      yf = frac_site_crd(2,site1)
      zf = frac_site_crd(3,site1)
      if ( use_cutoff == 1 )then
         if ( pbc == 1 .and.  extended_ucell == 0 )then!minimum image only
            do site2 = start,num_sites
               ! masking rule must be SKIP_MASKED_PAIRS
               if ( scratch_array(site2) /= site1 )then
                  ! in this case use fractional coords
                  df(1) = frac_site_crd(1,site2) - xf
                  df(2) = frac_site_crd(2,site2) - yf
                  df(3) = frac_site_crd(3,site2) - zf
                  df(1) = df(1) - dnint(df(1))
                  df(2) = df(2) - dnint(df(2))
                  df(3) = df(3) - dnint(df(3))
                  dx = ucell(1,1)*df(1) + ucell(1,2)*df(2) + &
                       ucell(1,3)*df(3)
                  dy = ucell(2,1)*df(1) + ucell(2,2)*df(2) + &
                       ucell(2,3)*df(3)
                  dz = ucell(3,1)*df(1) + ucell(3,2)*df(2) + &
                       ucell(3,3)*df(3)
                  dr2 = dx*dx + dy*dy + dz*dz
                  if ( dr2 < cut2 )then
                     num_per_site = num_per_site + 1
                  endif
               endif !( scratch_array(site2) /= site1 )then
            enddo !site2 = start,num_sites
         elseif ( pbc == 1 .and. extended_ucell > 0 )then !use translates
            do site2 = start,num_sites
               ! masking rule must be SKIP_MASKED_PAIRS
               if ( scratch_array(site2) /= site1 )then
                  ! in this case crd are actually mapped into central ucell
                  dx0 = map_site_crd(1,site2) - xm
                  dy0 = map_site_crd(2,site2) - ym
                  dz0 = map_site_crd(3,site2) - zm
                  do k = 1,num_translates
                     dx = dx0 + ucell_translate(1,k)
                     dy = dy0 + ucell_translate(2,k)
                     dz = dz0 + ucell_translate(3,k)
                     dr2 = dx*dx + dy*dy + dz*dz
                     if ( dr2 < cut2 )then
                        num_per_site = num_per_site + 1
                     endif
                  enddo !k = 1,num_translates
               endif !( scratch_array(site2) == site1 )
            enddo !site2 = site1+1,num_sites
         elseif ( pbc == 0 )then ! use reg crds and unimaged dx,dy,dz
            do site2 = start,num_sites
               ! masking rule must be SKIP_MASKED_PAIRS
               if ( scratch_array(site2) /= site1 )then
                  dx = map_site_crd(1,site2) - xm
                  dy = map_site_crd(2,site2) - ym
                  dz = map_site_crd(3,site2) - zm
                  dr2 = dx*dx + dy*dy + dz*dz
                  if ( dr2 < cut2 )then
                     num_per_site = num_per_site + 1
                  endif
               endif!( scratch_array(site2) /= site1 )
            enddo !site2 = start,num_sites
         endif !( pbc == 1 .and. extended_ucell == 0 )
      elseif ( use_cutoff == 0 )then
         ! extended_ucell must be 0, so minimum image--one pair per site2
         do site2 = start,num_sites
            proceed = 0
            if ( (masking_rule==SKIP_MASKED_PAIRS) .and. &
                 (scratch_array(site2) /= site1) )then
               proceed = 1
            elseif ( (masking_rule==ONLY_MASKED_PAIRS) .and. &
                 (scratch_array(site2) == site1) )then
               proceed = 1
            elseif ( (masking_rule==INCLUDE_MASKED_PAIRS) )then
               proceed = 1
            endif
            if ( proceed == 1 )then
               num_per_site = num_per_site + 1
            endif
         enddo !site2 = start,num_sites
      endif !( use_cutoff == 1 )then
      num_list = num_list + num_per_site
   enddo !site1 = 1,num_sites-1
end subroutine FH_SITE_SITE_LIST_count
!-----------------------------------------------------------------------
subroutine FH_SITE_SITE_LIST_build(& 
                    num_sites, &
                    num_sites_masked_by_site, &
                    masked_sitelist_offset_of_site, &
                    masked_sitelist, &
                    masking_rule, &
                    use_cutoff, &
                    pbc, &
                    extended_ucell, &
                    num_translates, &
                    cut_extra, &
                    ucell, &
                    ucell_translate, &
                    extent_of_compact_hermites, &
                    frac_site_crd, &
                    map_site_crd, &
                    siz_buffer, &
                    max_list, &
                    num_list, &
                    list, &
                    num_neighbors, &
                    list_offset, &
                    ucell_translate_index_list, &
                    max_num_neighbors, &
                    out_lun)

   use user, only : do_overlap, exchange_cutoff

   implicit none

   integer,intent(in) :: num_sites
   integer,intent(in) :: num_sites_masked_by_site(num_sites), &
                         masked_sitelist_offset_of_site(num_sites), &
                         masked_sitelist(*)
   integer,intent(in) :: max_list
   integer,intent(in) :: masking_rule,use_cutoff
   integer,intent(in) :: pbc,extended_ucell,num_translates
   double precision,intent(in) :: cut_extra
   double precision,intent(in) :: ucell(3,3),ucell_translate(3,*), &
                                  extent_of_compact_hermites, &
                                  frac_site_crd(3,*),map_site_crd(3,*)
   integer,intent(in) :: siz_buffer
   integer,intent(out) :: num_list,list(*),num_neighbors(num_sites), &
                          list_offset(num_sites), &
                          ucell_translate_index_list(*), &
                          max_num_neighbors
   integer,intent(in)   :: out_lun

   include 'masking_rules.fh'

   double precision, parameter  :: bohr = 0.52917721092d0
   double precision, parameter  :: bohr_per_angstrom = 1.d0 / bohr

   integer :: list_buffer(siz_buffer),trans_ind_buffer(siz_buffer)
   integer :: site1,site2,j,k,m1,m2,m3,scratch_array(num_sites),buf_count, &
              off,num,start,proceed
   double precision :: cut,cut2,df(3)
   double precision :: xf,yf,zf,xm,ym,zm,dx0,dy0,dz0,dx,dy,dz,dr2,dr
   
   ! check sanity of parameter combinations
   if ( (use_cutoff == 1) .and. &
        (  (masking_rule==INCLUDE_MASKED_PAIRS) .or. &
           (masking_rule==ONLY_MASKED_PAIRS)  )  )then
       write(out_lun,*)'FH_SITE_SITE_LIST_build: ', &
                 'use_cutoff=1 incompatible with INCLUDE_MASKED_PAIRS', &
                 ' or ONLY_MASKED_PAIRS'
       stop
   endif

   ! clear the scratch array
   call UTIL_zero_int_array(scratch_array,num_sites)

   num_list = 0
   if ( use_cutoff == 1 )then
      if (do_overlap .ne. 0) then
        cut = exchange_cutoff* bohr_per_angstrom  + cut_extra
        cut2 = cut*cut
      else
        cut = 2.d0*extent_of_compact_hermites + cut_extra
        cut2 = cut*cut
      end if
   endif
   do site1 = 1,num_sites
      off = masked_sitelist_offset_of_site(site1)
      num = num_sites_masked_by_site(site1)
      start = site1+1
      do j = 1,num
         k = masked_sitelist(off+j)
         scratch_array(k) = site1
      enddo
      xm = map_site_crd(1,site1)
      ym = map_site_crd(2,site1)
      zm = map_site_crd(3,site1)
      xf = frac_site_crd(1,site1)
      yf = frac_site_crd(2,site1)
      zf = frac_site_crd(3,site1)
      buf_count = 0
      if ( use_cutoff == 1 )then
         if ( pbc == 1 .and. extended_ucell == 0 )then!!minimum image only
            do site2 = start,num_sites
               ! masking rule must be SKIP_MASKED_PAIRS
               if ( scratch_array(site2) /= site1 )then
                  ! in this case use fractional coords
                  df(1) = frac_site_crd(1,site2) - xf
                  df(2) = frac_site_crd(2,site2) - yf
                  df(3) = frac_site_crd(3,site2) - zf
                  m1 = -dnint(df(1))
                  m2 = -dnint(df(2))
                  m3 = -dnint(df(3))
                  df(1) = df(1) + m1
                  df(2) = df(2) + m2
                  df(3) = df(3) + m3
                  dx = ucell(1,1)*df(1) + ucell(1,2)*df(2) + &
                       ucell(1,3)*df(3)
                  dy = ucell(2,1)*df(1) + ucell(2,2)*df(2) + &
                       ucell(2,3)*df(3)
                  dz = ucell(3,1)*df(1) + ucell(3,2)*df(2) + &
                       ucell(3,3)*df(3)
                  dr2 = dx*dx + dy*dy + dz*dz
                  if ( dr2 < cut2 )then
                     buf_count = buf_count + 1
                     list_buffer(buf_count) = site2
                     k = 9*(m3+1) + 3*(m2+1) + m1+1
                     trans_ind_buffer(buf_count) = k+1
                  endif
               endif !( scratch_array(site2) /= site1 )then
            enddo !site2 = start,num_sites
         elseif ( pbc == 1 .and. extended_ucell > 0 )then! translates
            do site2 = start,num_sites
               ! masking rule must be SKIP_MASKED_PAIRS
               if ( scratch_array(site2) /= site1 )then
                  ! in this case crd are actually mapped into central ucell
                  dx0 = map_site_crd(1,site2) - xm
                  dy0 = map_site_crd(2,site2) - ym
                  dz0 = map_site_crd(3,site2) - zm
                  do k = 1,num_translates
                     dx = dx0 + ucell_translate(1,k)
                     dy = dy0 + ucell_translate(2,k)
                     dz = dz0 + ucell_translate(3,k)
                     dr2 = dx*dx + dy*dy + dz*dz
                     if ( dr2 < cut2 )then
                        buf_count = buf_count + 1
                        list_buffer(buf_count) = site2
                        trans_ind_buffer(buf_count) = k
                     endif
                  enddo !k = 1,num_translates
               endif !( scratch_array(site2) == site1 )
            enddo !site2 = site1+1,num_sites
         elseif ( pbc == 0 )then ! use reg crds and unimaged dx,dy,dz
            do site2 = start,num_sites
               ! masking rule must be SKIP_MASKED_PAIRS
               if ( scratch_array(site2) /= site1 )then
                  dx = map_site_crd(1,site2) - xm
                  dy = map_site_crd(2,site2) - ym
                  dz = map_site_crd(3,site2) - zm
                  dr2 = dx*dx + dy*dy + dz*dz
                  if ( dr2 < cut2 )then
                     buf_count = buf_count + 1
                     list_buffer(buf_count) = site2
                  endif
               endif!( scratch_array(site2) /= site1 )
            enddo !site2 = start,num_sites
         endif !( pbc == 1 .and. extended_ucell == 0 )
      elseif ( use_cutoff == 0 )then
         ! if pbc minimum image--one pair per site2
         do site2 = start,num_sites
            proceed = 0
            if ( (masking_rule==SKIP_MASKED_PAIRS) .and. &
                 (scratch_array(site2) /= site1) )then
               proceed = 1
            elseif ( (masking_rule==ONLY_MASKED_PAIRS) .and. &
                 (scratch_array(site2) == site1) )then
               proceed = 1
            elseif ( (masking_rule==INCLUDE_MASKED_PAIRS) )then
               proceed = 1
            endif
            if ( proceed == 1 )then
               buf_count = buf_count + 1
               if ( pbc == 1 )then
                  ! in this case use fractional coords
                  df(1) = frac_site_crd(1,site2) - xf
                  df(2) = frac_site_crd(2,site2) - yf
                  df(3) = frac_site_crd(3,site2) - zf
                  m1 = -dnint(df(1))
                  m2 = -dnint(df(2))
                  m3 = -dnint(df(3))
                  list_buffer(buf_count) = site2
                  k = 9*(m3+1) + 3*(m2+1) + m1+1
                  trans_ind_buffer(buf_count) = k+1
               else
                  list_buffer(buf_count) = site2
               endif !( pbc == 1 )then
            endif !( proceed == 1 )then
         enddo !site2 = start,num_sites
      endif !( use_cutoff == 1 )then
      ! now fill regular lists
      num_neighbors(site1) = buf_count
      list_offset(site1) = num_list
      if ( buf_count > max_num_neighbors )max_num_neighbors = buf_count
      if ( num_list + buf_count > max_list )then
         write(out_lun,*)&
            'list exceeded: site1,num_list,buf_count,max_list = ', &
            site1,num_list,buf_count,max_list
         stop
      endif
      if ( pbc == 1 )then 
         do j = 1,buf_count
            list(num_list+j) = list_buffer(j)
            ucell_translate_index_list(num_list+j) = trans_ind_buffer(j)
         enddo
      else
         do j = 1,buf_count
            list(num_list+j) = list_buffer(j)
         enddo
      endif
      num_list = num_list + buf_count
   enddo !site1 = 1,num_sites-1
end subroutine FH_SITE_SITE_LIST_build

!-----------------------------------------------------------------------

subroutine FH_SITE_SITE_LIST_filter(  &
                               site1, &
                               num_site_nghbrs, &
                               list_offset, &
                               list, &
                               ucell_translate_index_list, &
                               map_site_crd, &
                               ucell_translate, &
                               extent_of_compact_hermites, &
                               pbc, &
                               use_cutoff, &
                               buf_count, &
                               neighbor_buf, &
                               delx, &
                               dely, &
                               delz, &
                               delr2 )

   use user, only : do_overlap, exchange_cutoff

   implicit none

   integer,intent(in) :: site1, &
                         num_site_nghbrs(*), &
                         list_offset(*), &
                         list(*), &
                         ucell_translate_index_list(*)
   double precision,intent(in) :: map_site_crd(3,*), &
                                  ucell_translate(3,*), &
                                  extent_of_compact_hermites
   integer,intent(in) :: pbc,use_cutoff
   integer,intent(out) :: buf_count,neighbor_buf(*)
   double precision,intent(out) :: delx(*),dely(*),delz(*),delr2(*)

   double precision, parameter  :: bohr = 0.52917721092d0
   double precision, parameter  :: bohr_per_angstrom = 1.d0 / bohr

   integer :: num,off,site2,ind,j
   double precision :: x,y,z,dx,dy,dz,dr2,cut,cut2

   buf_count = 0
   if ( use_cutoff == 1 )then
      if (do_overlap .ne. 0) then
        cut = exchange_cutoff * bohr_per_angstrom
        cut2 = cut*cut
      else
        cut = 2.d0*extent_of_compact_hermites
        cut2 = cut*cut
      end if
   endif
   num = num_site_nghbrs(site1)
   off = list_offset(site1)
   x = map_site_crd(1,site1)
   y = map_site_crd(2,site1)
   z = map_site_crd(3,site1)
   if ( use_cutoff == 1 .and. pbc == 1 )then
      do j = 1,num
         site2 = list(off+j)
         ind = ucell_translate_index_list(off+j)
         dx = map_site_crd(1,site2) + ucell_translate(1,ind) - x
         dy = map_site_crd(2,site2) + ucell_translate(2,ind) - y
         dz = map_site_crd(3,site2) + ucell_translate(3,ind) - z
         dr2 = dx*dx + dy*dy + dz*dz
         if ( dr2 < cut2 )then
            buf_count = buf_count + 1
            neighbor_buf(buf_count) = site2
            delx(buf_count) = dx
            dely(buf_count) = dy
            delz(buf_count) = dz
            delr2(buf_count) = dr2
         endif !( dr2 < cut2 )
      enddo !j = 1,num
   elseif ( use_cutoff == 1 .and. pbc == 0 )then !cutoff, no image
      do j = 1,num
         site2 = list(off+j)
         dx = map_site_crd(1,site2) - x
         dy = map_site_crd(2,site2) - y
         dz = map_site_crd(3,site2) - z
         dr2 = dx*dx + dy*dy + dz*dz
         if ( dr2 < cut2 )then
            buf_count = buf_count + 1
            neighbor_buf(buf_count) = site2
            delx(buf_count) = dx
            dely(buf_count) = dy
            delz(buf_count) = dz
            delr2(buf_count) = dr2
         endif !( dr2 < cut2 )
      enddo!j = 1,num
   elseif ( use_cutoff == 0 .and. pbc == 0 )then !no cutoff, no image
      do j = 1,num
         site2 = list(off+j)
         dx = map_site_crd(1,site2) - x
         dy = map_site_crd(2,site2) - y
         dz = map_site_crd(3,site2) - z
         dr2 = dx*dx + dy*dy + dz*dz
         buf_count = buf_count + 1
         neighbor_buf(buf_count) = site2
         delx(buf_count) = dx
         dely(buf_count) = dy
         delz(buf_count) = dz
         delr2(buf_count) = dr2
      enddo!j = 1,num
   elseif ( use_cutoff == 0 .and. pbc == 1 )then  !no cutoff, but image
      do j = 1,num
         site2 = list(off+j)
         ind = ucell_translate_index_list(off+j)
         dx = map_site_crd(1,site2) + ucell_translate(1,ind) - x
         dy = map_site_crd(2,site2) + ucell_translate(2,ind) - y
         dz = map_site_crd(3,site2) + ucell_translate(3,ind) - z
         dr2 = dx*dx + dy*dy + dz*dz
         buf_count = buf_count + 1
         neighbor_buf(buf_count) = site2
         delx(buf_count) = dx
         dely(buf_count) = dy
         delz(buf_count) = dz
         delr2(buf_count) = dr2
      enddo!j = 1,num
   endif !( compact_compact == 1 .and. pbc == 1 )
end subroutine FH_SITE_SITE_LIST_filter
!-----------------------------------------------------------------------
