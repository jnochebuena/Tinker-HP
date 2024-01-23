! NOTE! - A number of variables in the namelist are basically meaningless. I
!         have left them in so that using them will not cause an error, BUT
!         any implied effect on functionality does not exist if the variable is
!         listed here.  I also won't issue any warning messages, but may later
!         log a message in some sort of logfile if it becomes an issue. -- RED
!
!         dump_mol_forces - ignored, functionality not available

module user

   implicit none

   private
   
   integer, save :: user_num_site_types

   double precision, save :: gaussian_extent_tol, gaussian_recip_tol

   double precision, save :: erfcinv_gauss_extent_tol

   double precision, save :: exchange_factor

   double precision, save :: exchange_cutoff

   integer, save :: num_HC_prim_grids, num_HD_prim_grids

   integer, save :: user_num_aux_Cprims, user_num_aux_Dprims

   integer, save, dimension(100) :: HC_grid_number_for_Caux, &
                                    HD_grid_number_for_Daux

   integer, save, dimension(100) :: Bspline_order_for_grid_type

   integer, save, dimension(100) :: nfft1_for_grid_type, &
                                    nfft2_for_grid_type, &
                                    nfft3_for_grid_type

   double precision, save, dimension(3) :: unit_cell_length, unit_cell_angle

   integer,save :: pbc

   integer,save :: do_coulomb,do_overlap, &
                   calculation_method, &
                   Coulomb_use_recip,Overlap_use_recip

   integer, save :: orthogonal_ucell

   integer, save, dimension(100) :: grid_type_for_HC_prim_grid, &
                                    grid_type_for_HD_prim_grid
   integer, save                 :: grid_type_for_MPOLES, &
                                    grid_type_for_sumCH

   integer, save :: num_prim_grid_types

   integer, save, dimension(100) :: struc_fac_method_for_grid_type

   integer, save :: num_pme_grid_types, num_ffp_grid_types, &
                    num_ewald_grid_types

   double precision, save :: diffuse_adjust_exponent, cut_extra

   double precision, save :: CD_split_expon, extent_of_compact_hermites

   integer, save :: ffp_auto_setup, pme_auto_setup, reg_ewald_auto_setup

   integer,save :: verbose, dump_atom_forces, dump_mol_forces

   public FH_USER_read, &
          FH_USER_first_check, &
          FH_USER_second_check, &
          user_num_site_types, &
          gaussian_extent_tol, &
          erfcinv_gauss_extent_tol, &
          gaussian_recip_tol, &
          exchange_factor, &
          exchange_cutoff, &
          user_num_aux_Cprims, &
          user_num_aux_Dprims, &
          num_HC_prim_grids, &
          num_HD_prim_grids, &
          HC_grid_number_for_Caux, &
          HD_grid_number_for_Daux, &
          num_prim_grid_types, &
          grid_type_for_HC_prim_grid, &
          grid_type_for_HD_prim_grid, &
          grid_type_for_MPOLES, &
          grid_type_for_sumCH, &
          struc_fac_method_for_grid_type, &
          nfft1_for_grid_type, &
          nfft2_for_grid_type, &
          nfft3_for_grid_type, &
          Bspline_order_for_grid_type,&
          unit_cell_length, &
          unit_cell_angle, &
          pbc, &
          verbose, &
          dump_atom_forces, &
          dump_mol_forces, &
          do_coulomb, &
          do_overlap, &
          calculation_method, &
          Coulomb_use_recip, &
          Overlap_use_recip, &
          CD_split_expon, &
          extent_of_compact_hermites, &
          num_pme_grid_types, &
          num_ffp_grid_types, &
          num_ewald_grid_types, &
          ffp_auto_setup, &
          pme_auto_setup, &
          reg_ewald_auto_setup, &
          diffuse_adjust_exponent, &
          orthogonal_ucell, &
          cut_extra

contains

subroutine FH_USER_read(filename, free_lun, out_lun, pbc_box_size)

   implicit none

! Formal arguments:

   character(len=*), intent(in)                 :: filename
   integer, intent(in)                          :: free_lun
   integer, intent(in)                          :: out_lun
   double precision, optional, intent(in)       :: pbc_box_size(3)

! Namelist:

   namelist /user/ gaussian_extent_tol, &
                   CD_split_expon, &
                   gaussian_recip_tol, &
                   exchange_factor, &
                   exchange_cutoff, &
                   nfft1_for_grid_type, &
                   nfft2_for_grid_type, &
                   nfft3_for_grid_type, &
                   Bspline_order_for_grid_type, &
                   unit_cell_length, &
                   do_coulomb, &
                   do_overlap

! Local variables:

   integer      :: ios

   gaussian_extent_tol = 1.d-8
   CD_split_expon = 0.d0
   gaussian_recip_tol = 1.d-8
   exchange_factor = 6.6899d0
   exchange_cutoff = 6.d0
   user_num_site_types = 0
   user_num_aux_Cprims = 0
   user_num_aux_Dprims = 0
   num_HC_prim_grids = 0
   num_HD_prim_grids = 0
   HC_grid_number_for_Caux = 0
   HD_grid_number_for_Daux = 0
   num_prim_grid_types = 0
   grid_type_for_HC_prim_grid = 0
   grid_type_for_HD_prim_grid = 0
   grid_type_for_MPOLES = 0
   grid_type_for_sumCH = 0
   struc_fac_method_for_grid_type = 0
   nfft1_for_grid_type = 0
   nfft2_for_grid_type = 0
   nfft3_for_grid_type = 0
   Bspline_order_for_grid_type = 0
   unit_cell_length(:) = 0.d0
   unit_cell_angle(1) = 90.d0
   unit_cell_angle(2) = 90.d0
   unit_cell_angle(3) = 90.d0
   pbc = 1
   verbose = 0
   dump_atom_forces = 0
   dump_mol_forces = 0
   do_coulomb = 0
   do_overlap = 0
   Coulomb_use_recip = 1
   Overlap_use_recip = 0
   ffp_auto_setup = 0
   pme_auto_setup = 1
   reg_ewald_auto_setup = 0
   cut_extra = 1.d0             ! Bohr, sigh..

   open(unit=free_lun, status='old', file=filename)

   read(free_lun, nml=user)

   close(free_lun)

   if (present(pbc_box_size)) unit_cell_length(:) = pbc_box_size(:)

   if (verbose == 1) then
      write(out_lun, nml=user)
   endif

   return

end subroutine FH_USER_read

!-------------------------------------------------------------
subroutine FH_USER_first_check(do_overlap, out_lun)

   implicit none

! Formal arguments:

   integer, intent(in) :: do_overlap
   integer, intent(in) :: out_lun

! Local variables:

   include 'scale.fh'
   include 'interact_type.fh'

   ! if not an auto run then need some sanity checks

   if ((ffp_auto_setup==0) .and. (pme_auto_setup==0) .and. &
        (reg_ewald_auto_setup==0)) then

      if (user_num_site_types == 0 )then
         write(out_lun,*)'FH_USER_read: user_num_site_types = 0 !!'
         stop
      endif

   endif

   ! get erfcinv_gauss_extent_tol

   if (gaussian_extent_tol < very_small_value) then
      erfcinv_gauss_extent_tol = big_value ! no cutoffs
   else
      call FH_erfcinv(gaussian_extent_tol,erfcinv_gauss_extent_tol)
   endif

   ! get interact type

   if ( CD_split_expon < small_value )then

      if (do_coulomb .eq. 1) then
         write(out_lun,*)'Doing coulomb calc, but CD_split_expon is 0!'
         stop
      endif

   else

      extent_of_compact_hermites =  &
                 erfcinv_gauss_extent_tol / sqrt(CD_split_expon)

      if (Coulomb_use_recip == 1) then
         calculation_method = SPLIT_RECIP
      else
         calculation_method = SPLIT_DIRECT
      endif

   endif

   if (do_overlap .eq. 1) then
      calculation_method=NO_SPLIT
      extent_of_compact_hermites = big_value
   endif

   if (verbose == 1) then
     if (calculation_method == NO_SPLIT) then 
        write(out_lun,*)'calculation_method is NO_SPLIT'
     elseif ( calculation_method == SPLIT_DIRECT )then 
        write(out_lun,*)'calculation_method is SPLIT_DIRECT'
     elseif ( calculation_method == SPLIT_RECIP )then
        write(out_lun,*)'calculation_method is SPLIT_RECIP'
     endif
   endif

   ! catch all in case we're doing overlap

   if (do_overlap == 1 .and. CD_split_expon .ne. 0.d0) then
      write(out_lun,*)'do_overlap is 1 but CD_split_expon is != 0.d0!'
      write(out_lun,*)'setting CD_split_expon to 0.d0'
      CD_split_expon = 0.d0
      calculation_method = NO_SPLIT
      extent_of_compact_hermites = big_value
   endif

   if ( verbose == 1 )then
      write(out_lun,*)'extent of compact hermites = ',extent_of_compact_hermites
   endif

   return

end subroutine FH_USER_first_check

!-------------------------------------------------------------
subroutine FH_USER_second_check(out_lun)

   implicit none

! Formal arguments:

   integer, intent(in) :: out_lun

! Local variables:

   include 'structure_factor_type.fh'
   include 'mpole_sizes.fh'
   include 'scale.fh'
   include 'interact_type.fh'

   double precision bohr
   integer g,gt,n,j,needit


   if ( calculation_method == SPLIT_RECIP  )then 
      !check on ffp with diffuse_adjust_exponent
      needit = 0
      ! CD_split_expon must be positive (FH_USER_first_check)
      ! check the compact hermites
      do g = 1,num_HC_prim_grids
         gt = grid_type_for_HC_prim_grid(g)
         if ( struc_fac_method_for_grid_type(gt) == SF_FFP)needit = 1
      enddo
      ! check the compact multipoles
      gt = grid_type_for_MPOLES
      if ( struc_fac_method_for_grid_type(gt) == SF_FFP)needit = 1
      ! sumCH must be active
      gt = grid_type_for_sumCH
      if ( struc_fac_method_for_grid_type(gt) == SF_FFP)needit = 1
      if ( needit == 1 ) then
         diffuse_adjust_exponent = 2.d0*CD_split_expon
      else
         diffuse_adjust_exponent = 0.d0
      endif
   endif !( calculation_method == SPLIT_RECIP  )
   if ( (pbc == 1) .and. &
        ((do_coulomb == 1) .and. (coulomb_use_recip == 0)) )then
      write(out_lun,*)'need reciprocal sum for PBC!!'
      stop
   endif 
   ! check sanity of grid types
   if ( num_HC_prim_grids > 0 )then
      do n = 1,user_num_aux_Cprims
         if ( HC_grid_number_for_Caux(n) < 0 .or. &
                HC_grid_number_for_Caux(n) > num_HC_prim_grids )then
            write(out_lun,*)'HC_grid_number_for_aux(n) wrong! ', &
                   'n,HC_grid_number_for_aux(n),num_HC_prim_grids = ',&
                    n,HC_grid_number_for_Caux(n),num_HC_prim_grids
            stop
         endif
      enddo
   endif
   if ( num_HD_prim_grids > 0 )then
      do n = 1,user_num_aux_Dprims
         if ( HD_grid_number_for_Daux(n) < 0 .or. &
                HD_grid_number_for_Daux(n) > num_HD_prim_grids )then
            write(out_lun,*)'HD_grid_number_for_aux(n) wrong! ', &
                   'n,HD_grid_number_for_aux(n),num_HD_prim_grids = ',&
                    n,HD_grid_number_for_Daux(n),num_HD_prim_grids
            stop
         endif
      enddo
   endif
   if ((Coulomb_use_recip == 1) .or. (Overlap_use_recip == 1))then
      if ( num_prim_grid_types <= 0 .or. &
            ! recall there are multipole grids and sumCH_grid
            num_prim_grid_types > num_HC_prim_grids + &
                                  num_HD_prim_grids + 2 )then
         write(out_lun,*)'recip sum  but num_prim_grid_types wrong'
         write(out_lun,*)'num_prim_grid_types,num_HC_prim_grids,', &
                   'num_HD_prim_grids = ',num_prim_grid_types, &
                    num_HC_prim_grids,num_HD_prim_grids
         stop
      endif
      do n = 1,num_HC_prim_grids
         if ( grid_type_for_HC_prim_grid(n) <= 0 .or. &
              grid_type_for_HC_prim_grid(n) > num_prim_grid_types )then
            write(out_lun,*)'grid_type_for_HC_prim_grid(n) wrong! ', &
                      'n,grid_type_for_HC_prim_grid(n),num_prim_grid_types = ',&
                       n,grid_type_for_HC_prim_grid(n),num_prim_grid_types
            stop
         endif
      enddo
      do n = 1,num_HD_prim_grids
         if ( grid_type_for_HD_prim_grid(n) <= 0 .or. &
              grid_type_for_HD_prim_grid(n) > num_prim_grid_types )then
            write(out_lun,*)'grid_type_for_HD_prim_grid(n) wrong! ', &
                      'n,grid_type_for_HD_prim_grid(n),num_prim_grid_types = ',&
                       n,grid_type_for_HD_prim_grid(n),num_prim_grid_types
            stop
         endif
      enddo
      
      if ( grid_type_for_MPOLES <= 0 .or. &
           grid_type_for_MPOLES > num_prim_grid_types )then
            write(out_lun,*)'grid_type_for_MPOLES wrong! ', &
                      'grid_type_for_MPOLES,num_prim_grid_types = ',&
                       grid_type_for_MPOLES,num_prim_grid_types
            stop
      endif

      if ( num_HC_prim_grids > 0 )then ! check sumCH_grid
         if ( grid_type_for_sumCH <= 0 .or. &
              grid_type_for_sumCH > num_prim_grid_types )then
            write(out_lun,*)'grid_type_for_sumCH wrong! ', &
                      'grid_type_for_sumCH,num_prim_grid_types = ',&
                       grid_type_for_sumCH,num_prim_grid_types
            stop
         endif
      endif
   endif
   num_pme_grid_types = 0
   num_ffp_grid_types = 0
   num_ewald_grid_types = 0
   if ( num_prim_grid_types > 0 )then
      do gt = 1,num_prim_grid_types
         if ( struc_fac_method_for_grid_type(gt) == SF_PME )then 
            num_pme_grid_types = num_pme_grid_types + 1
         elseif ( struc_fac_method_for_grid_type(gt) == SF_FFP )then 
            num_ffp_grid_types = num_ffp_grid_types + 1
         elseif ( struc_fac_method_for_grid_type(gt) == SF_EWALD )then 
            num_ewald_grid_types = num_ewald_grid_types + 1
         else
            write(out_lun,*)'bad structure factor type for grid type ',gt
            stop
         endif
      enddo
      if ( verbose == 1 )then
         write(out_lun,*)'user: num pme,ffp,ewald grid types = ', &
            num_pme_grid_types,num_ffp_grid_types,num_ewald_grid_types
      endif
      do gt = 1,num_prim_grid_types
         if ( (struc_fac_method_for_grid_type(gt) == SF_PME) .and. &
              ( Bspline_order_for_grid_type(gt) == 0 .or.  &
                Bspline_order_for_grid_type(gt) > MAXSplineOrder) )then
            write(out_lun,*)'bad spline order for grid type ',gt
            stop
         endif
      enddo
   endif

!  bohr=0.529177249d0
   bohr=0.52917721092d0
   if ( pbc == 1 )then
      do n = 1,3
        if ( unit_cell_length(n) < 1.d-3 )then
           write(out_lun,*)'bad unit cell length'
           stop
        endif
        unit_cell_length(n) = unit_cell_length(n) / bohr
        if ( unit_cell_angle(n) < 1.d-3 .or. &
              unit_cell_angle(n) > 180.d0-1.d-3 )then
           write(out_lun,*)'bad unit cell angle'
           stop
        endif
      enddo
   endif
   if ( (pbc==0) .or. ((pbc==1).and.(unit_cell_angle(1)==90.0).and. &
        (unit_cell_angle(2)==90.0).and.(unit_cell_angle(3)==90.0))  )then
      orthogonal_ucell = 1
   else
      orthogonal_ucell = 0
   endif

   return
   
end subroutine FH_USER_second_check

!-------------------------------------------------------------

end module user
