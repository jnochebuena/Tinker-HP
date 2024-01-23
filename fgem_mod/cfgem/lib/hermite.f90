module hermite

   implicit none

   private

   integer, save, allocatable :: max_xform_order_for_site(:)
   ! multipole stuff
   integer, save              :: tot_num_mpole_coeffs
   integer, save, allocatable :: mpole_order_for_site(:), &
                                 mpole_level_for_site(:), &
                                 mpole_coeff_off_for_site(:)
   double precision, save, allocatable :: Local_multipole_coeff(:), &
                                          Global_multipole_coeff(:), &
                                          Global_multipole_field(:)
   double precision, save     :: compact_mpole_expon, &
                                 diffuse_mpole_expon, &
                                 compact_rec_mpole_expon, &
                                 diffuse_rec_mpole_expon

   ! sum of compact hermites stuff
   integer, save              :: tot_num_sumCH_coeffs
   integer, save, allocatable :: sumCH_order_for_site(:), &
                                 sumCH_level_for_site(:), &
                                 sumCH_coeff_off_for_site(:)
   double precision, save, allocatable :: Local_sumCH_coeff(:), &
                                          Global_sumCH_coeff(:), &
                                          Global_sumCH_field(:)
   double precision, save     :: compact_sumCH_expon, &
                                 diffuse_sumCH_expon, &
                                 compact_rec_sumCH_expon, &
                                 diffuse_rec_sumCH_expon

   ! compact hermite stuff
   integer, save              :: tot_num_herm_Cprims, &
                                 tot_num_Cherm_coeffs
   integer, save, allocatable :: num_herm_Cprim_for_site(:), &
                                 off_herm_Cprim_for_site(:)
   integer, save, allocatable :: hermite_order_for_Cprim(:), &
                                 hermite_level_for_Cprim(:), &
                                 site_that_owns_Cprim(:), &
                                 herm_coeff_offset_of_Cprim(:), &
                                 HC_grid_number_for_Cprim(:)
   double precision, save, allocatable :: hermite_expon_for_Cprim(:), &
                                          hermite_rec_expon_for_Cprim(:), &
                                          gauss_extent_for_Cprim(:), &
                                          Local_Chermite_coeff(:), &
                                          Global_Chermite_coeff(:), &
                                          Global_Chermite_field(:), &
                                          Global_Chermite_field_EXCH(:)
                                          !GAC: for EXCHANGE, the fields
                                          !     need to include the factor

   ! diffuse hermite stuff
   integer, save              :: tot_num_herm_Dprims, &
                                 tot_num_Dherm_coeffs
   integer, save, allocatable :: num_herm_Dprim_for_site(:), &
                                 off_herm_Dprim_for_site(:)
   integer, save, allocatable :: hermite_order_for_Dprim(:), &
                                 hermite_level_for_Dprim(:), &
                                 site_that_owns_Dprim(:), &
                                 herm_coeff_offset_of_Dprim(:), &
                                 HD_grid_number_for_Dprim(:)
   double precision, save, allocatable :: hermite_expon_for_Dprim(:), &
                                          hermite_rec_expon_for_Dprim(:), &
                                          gauss_extent_for_Dprim(:), &
                                          Local_Dhermite_coeff(:), &
                                          Global_Dhermite_coeff(:), &
                                          Global_Dhermite_field(:), &
                                          Global_Dhermite_field_EXCH(:)

   double precision, save,dimension(3) :: lower_density_extent, &
                                          upper_density_extent

   public  FH_hermite_load, &
           FH_hermite_deallocate, &
           FH_hermite_density_extent, &
           FH_hermite_local_to_global, &
           FH_hermite_field_de_drot, &
           FH_hermite_field_de_drot_EXCH, &
           FH_hermite_field_energy, &
           FH_hermite_field_energy_EXCH, &
           lower_density_extent, &
           upper_density_extent, &
           ! mpole stuff
           tot_num_mpole_coeffs, &
           mpole_level_for_site, &
           mpole_order_for_site, &
           mpole_coeff_off_for_site, &
           compact_rec_mpole_expon, &
           diffuse_rec_mpole_expon, &
           Global_multipole_coeff, &
           Global_multipole_field, &
           ! sumCH stuff
           tot_num_sumCH_coeffs, &
           sumCH_level_for_site, &
           sumCH_order_for_site, &
           sumCH_coeff_off_for_site, &
           compact_rec_sumCH_expon, &
           diffuse_rec_sumCH_expon, &
           Global_sumCH_coeff, &
           Global_sumCH_field, &
           ! compact hermite stuff
           tot_num_herm_Cprims, &
           num_herm_Cprim_for_site, &
           off_herm_Cprim_for_site, &
           hermite_level_for_Cprim, &
           hermite_order_for_Cprim, &
           hermite_rec_expon_for_Cprim, &
           herm_coeff_offset_of_Cprim, &
           site_that_owns_Cprim, &
           hermite_expon_for_Cprim, &
           gauss_extent_for_Cprim, &
           HC_grid_number_for_Cprim, &
           Local_Chermite_coeff, &
           Global_Chermite_coeff, &
           Global_Chermite_field, &
           Global_Chermite_field_EXCH, &
           ! diffuse hermite stuff
           tot_num_herm_Dprims, &
           num_herm_Dprim_for_site, &
           off_herm_Dprim_for_site, &
           hermite_level_for_Dprim, &
           hermite_order_for_Dprim, &
           hermite_rec_expon_for_Dprim, &
           herm_coeff_offset_of_Dprim, &
           site_that_owns_Dprim, &
           hermite_expon_for_Dprim, &
           gauss_extent_for_Dprim, &
           HD_grid_number_for_Dprim, &
           Local_Dhermite_coeff, &
           Global_Dhermite_coeff, &
           Global_Dhermite_field, &
           Global_Dhermite_field_EXCH
           
contains
!-----------------------------------------------------------------
!-----------------------------------------------------------------
! SETUP SUBROUTINES
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!----------------------------------------------------------------

subroutine FH_hermite_load(out_lun)

   use user,only : gaussian_extent_tol,HC_grid_number_for_Caux, &
                   HD_grid_number_for_Daux,CD_split_expon,verbose, &
                   diffuse_adjust_exponent

   use sites,only : num_sites,aux_type_of_site,site_crd

   use auxiliary, only : num_aux_site_types, &
                         tot_num_aux_Cprim, &
                         num_aux_Cprim_for_site_type, &
                         off_aux_Cprim_for_site_type, &
                         herm_exp_for_aux_Cprim, &
                         herm_rec_exp_for_aux_Cprim, &
                         herm_ord_for_aux_Cprim, &
                         herm_lev_for_aux_Cprim, &
                         herm_coeff_off_for_aux_Cprim, &
                         aux_Chermite_coeff, &
                         tot_num_aux_Dprim, &
                         num_aux_Dprim_for_site_type, &
                         off_aux_Dprim_for_site_type, &
                         herm_exp_for_aux_Dprim, &
                         herm_rec_exp_for_aux_Dprim, &
                         herm_ord_for_aux_Dprim, &
                         herm_lev_for_aux_Dprim, &
                         herm_coeff_off_for_aux_Dprim, &
                         aux_Dhermite_coeff, &
                         sumCH_order_for_site_type, &
                         sumCH_level_for_site_type, &
                         sumCH_coeff_off_for_site_type, &
                         aux_sumCH_coeff, &
                         mpole_order_for_site_type, &
                         mpole_level_for_site_type, &
                         mpole_coeff_off_for_site_type, &
                         aux_mpole_coeff

   implicit none

! Formal arguments:

   integer, intent(in)  :: out_lun

! Local variables:

   include 'scale.fh'
   double precision :: erfcinv_tol,expon,expon1
   integer :: j,k,n,np,nauxp,offprim,offaux,itype,numprim,order, &
              offcoeff,offauxcoeff,max_order
   integer :: ier1,ier2,ier3,ier4,ier5,ier6,ier7,ier8,ier9,ier10,ier11

   if ( gaussian_extent_tol < 1.d-20 )then
      erfcinv_tol = 1.d10 ! no cutoffs
   else
      call FH_erfcinv(gaussian_extent_tol,erfcinv_tol)
   endif

   allocate(max_xform_order_for_site(num_sites),stat=ier1)
   if ( ier1 /= 0 )then
      write(out_lun,*)'FH_hermite_load: allocate fails!'
      stop
   endif
   ! setup mpoles
   allocate(mpole_order_for_site(num_sites),stat=ier1)
   allocate(mpole_level_for_site(num_sites),stat=ier2)
   allocate(mpole_coeff_off_for_site(num_sites),stat=ier3)
   if ( ier1 /= 0 .or. ier2 /= 0 .or. ier3 /= 0 )then
      write(out_lun,*)'FH_hermite_load: allocate fails!'
      stop
   endif
   do n = 1,num_sites
      itype = aux_type_of_site(n)
      mpole_order_for_site(n) = mpole_order_for_site_type(itype)
      mpole_level_for_site(n) = mpole_level_for_site_type(itype)
      max_xform_order_for_site(n) = mpole_order_for_site(n)
   enddo
   mpole_coeff_off_for_site(1) = 0
   do n = 2,num_sites
      mpole_coeff_off_for_site(n) = mpole_coeff_off_for_site(n-1) + &
                                    mpole_order_for_site(n-1)
   enddo
   tot_num_mpole_coeffs = mpole_coeff_off_for_site(num_sites) + &
                          mpole_order_for_site(num_sites)
   allocate(Local_multipole_coeff(tot_num_mpole_coeffs),stat=ier1)
   allocate(Global_multipole_coeff(tot_num_mpole_coeffs),stat=ier2)
   allocate(Global_multipole_field(tot_num_mpole_coeffs),stat=ier3)
   if ( ier1 /= 0 .or. ier2 /= 0 .or. ier3 /= 0 )then
      write(out_lun,*)'FH_hermite_load: allocate fails!'
      stop
   endif
   do n = 1,num_sites
      itype = aux_type_of_site(n)
      offcoeff = mpole_coeff_off_for_site(n)
      offauxcoeff = mpole_coeff_off_for_site_type(itype)
      order = mpole_order_for_site(n)
      if ( order > 0 )then
         do k = 1,order
            Local_multipole_coeff(offcoeff+k) = aux_mpole_coeff(offauxcoeff+k)
         enddo
      endif
   enddo
   ! next mpole exponents
   compact_mpole_expon = big_value
   diffuse_mpole_expon = 0.5d0*cd_split_expon
   compact_rec_mpole_expon = compact_mpole_expon
   diffuse_rec_mpole_expon = diffuse_mpole_expon
   if ( diffuse_adjust_exponent > 1.d-20 )then
      ! make compact more diffuse
      compact_rec_mpole_expon = diffuse_adjust_exponent
      ! make diffuse less diffuse
      expon = diffuse_rec_mpole_expon
      expon = (expon*diffuse_adjust_exponent) / &
              (diffuse_adjust_exponent - expon)
      diffuse_rec_mpole_expon = expon
   endif
   if ( verbose == 1 )then
      write(out_lun,*)'tot_num_mpole_coeffs = ',tot_num_mpole_coeffs
   endif

   ! set up sumCH if appropriate
   if ( CD_split_expon > small_value )then
      allocate(sumCH_order_for_site(num_sites),stat=ier1)
      allocate(sumCH_level_for_site(num_sites),stat=ier2)
      allocate(sumCH_coeff_off_for_site(num_sites),stat=ier3)
      if ( ier1 /= 0 .or. ier2 /= 0 .or. ier3 /= 0 )then
         write(out_lun,*)'FH_hermite_load: allocate fails!'
         stop
      endif
      do n = 1,num_sites
         itype = aux_type_of_site(n)
         sumCH_order_for_site(n) = sumCH_order_for_site_type(itype)
         sumCH_level_for_site(n) = sumCH_level_for_site_type(itype)
         if ( sumCH_order_for_site(n) > max_xform_order_for_site(n) ) &
             max_xform_order_for_site(n) = sumCH_order_for_site(n)
      enddo
      sumCH_coeff_off_for_site(1) = 0
      do n = 2,num_sites
         sumCH_coeff_off_for_site(n) = sumCH_coeff_off_for_site(n-1) + &
                                       sumCH_order_for_site(n-1)
      enddo
      tot_num_sumCH_coeffs = sumCH_coeff_off_for_site(num_sites) + &
                             sumCH_order_for_site(num_sites)
      allocate(Local_sumCH_coeff(tot_num_sumCH_coeffs),stat=ier1)
      allocate(Global_sumCH_coeff(tot_num_sumCH_coeffs),stat=ier2)
      allocate(Global_sumCH_field(tot_num_sumCH_coeffs),stat=ier3)
      if ( ier1 /= 0 .or. ier2 /= 0 .or. ier3 /= 0 )then
         write(out_lun,*)'FH_hermite_load: allocate fails!'
         stop
      endif
      do n = 1,num_sites
         itype = aux_type_of_site(n)
         offcoeff = sumCH_coeff_off_for_site(n)
         offauxcoeff = sumCH_coeff_off_for_site_type(itype)
         order = sumCH_order_for_site(n)
         if ( order > 0 )then
            do k = 1,order
               Local_sumCH_coeff(offcoeff+k) =  &
                    aux_sumCH_coeff(offauxcoeff+k)
            enddo
         endif
      enddo
      ! next sumCH exponents
      compact_sumCH_expon = big_value
      diffuse_sumCH_expon = 0.5d0*CD_split_expon
      compact_rec_sumCH_expon = compact_sumCH_expon
      diffuse_rec_sumCH_expon = diffuse_sumCH_expon
      if ( diffuse_adjust_exponent > small_value )then
         ! make compact more diffuse
         compact_rec_sumCH_expon = diffuse_adjust_exponent
         ! make diffuse less diffuse
         expon = diffuse_rec_sumCH_expon
         expon = (expon*diffuse_adjust_exponent) / &
              (diffuse_adjust_exponent - expon)
         diffuse_rec_sumCH_expon = expon
      endif
   else
      tot_num_sumCH_coeffs = 0
   endif
   if ( verbose == 1 )then
      write(out_lun,*)'tot_num_sumCH_coeffs = ',tot_num_sumCH_coeffs
   endif

   ! now the compact hermites
   allocate(num_herm_Cprim_for_site(num_sites),stat=ier1)
   allocate(off_herm_Cprim_for_site(num_sites),stat=ier2)
   if ( ier1 /= 0 .or. ier2 /= 0 )then
      write(out_lun,*)'FH_hermite_load: allocate fails!'
      stop
   endif
   do n = 1,num_sites
      itype = aux_type_of_site(n)
      num_herm_Cprim_for_site(n) = num_aux_Cprim_for_site_type(itype)
   enddo
   off_herm_Cprim_for_site(1) = 0
   do n = 2,num_sites
      off_herm_Cprim_for_site(n) = off_herm_Cprim_for_site(n-1) + &
                                  num_herm_Cprim_for_site(n-1)
   enddo
   tot_num_herm_Cprims = off_herm_Cprim_for_site(num_sites) + &
                         num_herm_Cprim_for_site(num_sites)
   if ( tot_num_herm_Cprims > 0 )then
      allocate(hermite_expon_for_Cprim(tot_num_herm_Cprims),stat=ier1)
      allocate(hermite_rec_expon_for_Cprim(tot_num_herm_Cprims),stat=ier2)
      allocate(gauss_extent_for_Cprim(tot_num_herm_Cprims),stat=ier3)
      allocate(hermite_order_for_Cprim(tot_num_herm_Cprims),stat=ier4)
      allocate(hermite_level_for_Cprim(tot_num_herm_Cprims),stat=ier5)
      allocate(site_that_owns_Cprim(tot_num_herm_Cprims),stat=ier6)
      allocate(herm_coeff_offset_of_Cprim(tot_num_herm_Cprims),stat=ier7)
      allocate(HC_grid_number_for_Cprim(tot_num_herm_Cprims),stat=ier8)
      if ( ier1 /= 0 .or. ier2 /= 0 .or. ier3 /= 0 .or. ier4 /= 0 .or. &
           ier5 /= 0 .or. ier6 /= 0 .or. ier7 /= 0 .or. ier8 /= 0  )then
         write(out_lun,*)'FH_hermite_load: allocate fails!'
         stop
      endif
      
      do n = 1,num_sites
         itype = aux_type_of_site(n)
         offprim = off_herm_Cprim_for_site(n)
         offaux = off_aux_Cprim_for_site_type(itype)
         numprim = num_herm_Cprim_for_site(n)
         do k = 1,numprim
            np = offprim + k
            nauxp = offaux + k
            expon = herm_exp_for_aux_Cprim(nauxp)
            expon1 = herm_rec_exp_for_aux_Cprim(nauxp)
            hermite_expon_for_Cprim(np) = expon
            hermite_rec_expon_for_Cprim(np) = expon1
            gauss_extent_for_Cprim(np) = erfcinv_tol / sqrt(expon)
            site_that_owns_Cprim(np) = n
            hermite_order_for_Cprim(np) = herm_ord_for_aux_Cprim(nauxp)
            hermite_level_for_Cprim(np) = herm_lev_for_aux_Cprim(nauxp)
            HC_grid_number_for_Cprim(np) = HC_grid_number_for_Caux(nauxp)
         enddo
      enddo
      herm_coeff_offset_of_Cprim(1) = 0
      do np = 2,tot_num_herm_Cprims
         herm_coeff_offset_of_Cprim(np) = herm_coeff_offset_of_Cprim(np-1) + &
                                          hermite_order_for_Cprim(np-1)
      enddo
      tot_num_Cherm_coeffs =  &
                herm_coeff_offset_of_Cprim(tot_num_herm_Cprims) + &
                hermite_order_for_Cprim(tot_num_herm_Cprims)
      allocate(Local_Chermite_coeff(tot_num_Cherm_coeffs),stat=ier1)
      allocate(Global_Chermite_coeff(tot_num_Cherm_coeffs),stat=ier2)
      allocate(Global_Chermite_field(tot_num_Cherm_coeffs),stat=ier3)
      allocate(Global_Chermite_field_EXCH(tot_num_Cherm_coeffs),stat=ier3)
      if ( ier1 /= 0 .or. ier2 /= 0 .or. ier3 /= 0 )then
         write(out_lun,*)'FH_hermite_load: allocate fails!'
         stop
      endif
      do n = 1,num_sites
         itype = aux_type_of_site(n)
         offprim = off_herm_Cprim_for_site(n)
         offaux = off_aux_Cprim_for_site_type(itype)
         numprim = num_herm_Cprim_for_site(n)
         max_order = -1
         do k = 1,numprim
            np = offprim + k
            nauxp = offaux + k
            order = hermite_order_for_Cprim(np)
            if ( order > max_order )max_order = order
            offcoeff = herm_coeff_offset_of_Cprim(np)
            offauxcoeff = herm_coeff_off_for_aux_Cprim(nauxp)
            do j = 1,order
               Local_Chermite_coeff(offcoeff+j) =  &
                         aux_Chermite_coeff(offauxcoeff+j)
            enddo
         enddo
         if ( max_order > max_xform_order_for_site(n) ) &
               max_xform_order_for_site(n) = max_order
      enddo !n = 1,num_sites
   else
      tot_num_Cherm_coeffs = 0
   endif !( tot_num_herm_Cprims > 0 )
   if ( verbose == 1 )then
      write(out_lun,*)'tot_num_herm_Cprims = ',tot_num_herm_Cprims
      write(out_lun,*)'tot_num_Cherm_coeffs = ',tot_num_Cherm_coeffs
   endif

   ! now the diffuse hermites
   allocate(num_herm_Dprim_for_site(num_sites),stat=ier1)
   allocate(off_herm_Dprim_for_site(num_sites),stat=ier2)
   if ( ier1 /= 0 .or. ier2 /= 0 )then
      write(out_lun,*)'FH_hermite_load: allocate fails!'
      stop
   endif
   do n = 1,num_sites
      itype = aux_type_of_site(n)
      num_herm_Dprim_for_site(n) = num_aux_Dprim_for_site_type(itype)
   enddo
   off_herm_Dprim_for_site(1) = 0
   do n = 2,num_sites
      off_herm_Dprim_for_site(n) = off_herm_Dprim_for_site(n-1) + &
                                  num_herm_Dprim_for_site(n-1)
   enddo
   tot_num_herm_Dprims = off_herm_Dprim_for_site(num_sites) + &
                         num_herm_Dprim_for_site(num_sites)


   if ( tot_num_herm_Dprims > 0 )then
      allocate(hermite_expon_for_Dprim(tot_num_herm_Dprims),stat=ier1)
      allocate(hermite_rec_expon_for_Dprim(tot_num_herm_Dprims),stat=ier2)
      allocate(gauss_extent_for_Dprim(tot_num_herm_Dprims),stat=ier3)
      allocate(hermite_order_for_Dprim(tot_num_herm_Dprims),stat=ier4)
      allocate(hermite_level_for_Dprim(tot_num_herm_Dprims),stat=ier5)
      allocate(site_that_owns_Dprim(tot_num_herm_Dprims),stat=ier6)
      allocate(herm_coeff_offset_of_Dprim(tot_num_herm_Dprims),stat=ier7)
      allocate(HD_grid_number_for_Dprim(tot_num_herm_Dprims),stat=ier8)
      if ( ier1 /= 0 .or. ier2 /= 0 .or. ier3 /= 0 .or. ier4 /= 0 .or. &
           ier5 /= 0 .or. ier6 /= 0 .or. ier7 /= 0 .or. ier8 /= 0  )then
         write(out_lun,*)'FH_hermite_load: allocate fails!'
         stop
      endif
      do n = 1,num_sites
         itype = aux_type_of_site(n)
         offprim = off_herm_Dprim_for_site(n)
         offaux = off_aux_Dprim_for_site_type(itype)
         numprim = num_herm_Dprim_for_site(n)
         do k = 1,numprim
            np = offprim + k
            nauxp = offaux + k
            expon = herm_exp_for_aux_Dprim(nauxp)
            expon1 = herm_rec_exp_for_aux_Dprim(nauxp)
            hermite_expon_for_Dprim(np) = expon
            hermite_rec_expon_for_Dprim(np) = expon1
            gauss_extent_for_Dprim(np) = erfcinv_tol / sqrt(expon)
            site_that_owns_Dprim(np) = n
            hermite_order_for_Dprim(np) = herm_ord_for_aux_Dprim(nauxp)
            hermite_level_for_Dprim(np) = herm_lev_for_aux_Dprim(nauxp)
            HD_grid_number_for_Dprim(np) = HD_grid_number_for_Daux(nauxp)
         enddo
      enddo
      herm_coeff_offset_of_Dprim(1) = 0
      do np = 2,tot_num_herm_Dprims
         herm_coeff_offset_of_Dprim(np) =  &
                   herm_coeff_offset_of_Dprim(np-1) + &
                   hermite_order_for_Dprim(np-1)
      enddo
      tot_num_Dherm_coeffs =  &
                         herm_coeff_offset_of_Dprim(tot_num_herm_Dprims) + &
                         hermite_order_for_Dprim(tot_num_herm_Dprims)
      allocate(Local_Dhermite_coeff(tot_num_Dherm_coeffs),stat=ier1)
      allocate(Global_Dhermite_coeff(tot_num_Dherm_coeffs),stat=ier2)
      allocate(Global_Dhermite_field(tot_num_Dherm_coeffs),stat=ier3)
      allocate(Global_Dhermite_field_EXCH(tot_num_Dherm_coeffs),stat=ier3)
      if ( ier1 /= 0 .or. ier2 /= 0 .or. ier3 /= 0 )then
         write(out_lun,*)'FH_hermite_load: allocate fails!'
         stop
      endif
      do n = 1,num_sites
         itype = aux_type_of_site(n)
         offprim = off_herm_Dprim_for_site(n)
         offaux = off_aux_Dprim_for_site_type(itype)
         numprim = num_herm_Dprim_for_site(n)
         max_order = -1
         do k = 1,numprim
            np = offprim + k
            nauxp = offaux + k
            order = hermite_order_for_Dprim(np)
            if ( order > max_order )max_order = order
            offcoeff = herm_coeff_offset_of_Dprim(np)
            offauxcoeff = herm_coeff_off_for_aux_Dprim(nauxp)
            do j = 1,order
               Local_Dhermite_coeff(offcoeff+j) =  &
                  aux_Dhermite_coeff(offauxcoeff+j)
            enddo
         enddo
         if ( max_order > max_xform_order_for_site(n) ) &
               max_xform_order_for_site(n) = max_order
      enddo 
   else
      tot_num_Dherm_coeffs = 0
   endif !( tot_num_herm_Dprims > 0 )then
   if ( verbose == 1 )then
      write(out_lun,*)'tot_num_herm_Dprims = ',tot_num_herm_Dprims
      write(out_lun,*)'tot_num_Dherm_coeffs = ',tot_num_Dherm_coeffs
   endif

end subroutine FH_hermite_load

!----------------------------------------------------------------
subroutine FH_hermite_density_extent(out_lun)

   use user, only : verbose,calculation_method
   use sites, only : num_sites,site_crd,aux_type_of_site
   use auxiliary, only : max_extent_for_site_type

   implicit none

! Formal arguments:

   integer, intent(in)  :: out_lun

! Local variables:

   integer k,n,itype
   double precision r
   include 'scale.fh'
   include 'interact_type.fh'

   if ( calculation_method /= SPLIT_RECIP )return  !not using fourier methods
   ! get the density extent
   do k = 1,3
      lower_density_extent(k) = big_value
      upper_density_extent(k) = -big_value
   enddo
   ! handle multipoles and sumCH
   do n = 1,num_sites
      itype = aux_type_of_site(n)
      r = max_extent_for_site_type(itype)
      do k = 1,3
         if ( site_crd(k,n) - r < lower_density_extent(k) ) &
              lower_density_extent(k) = site_crd(k,n) - r
         if ( site_crd(k,n) + r > upper_density_extent(k) ) &
              upper_density_extent(k) = site_crd(k,n) + r
      enddo
   enddo
   if ( verbose == 1 )then
      write(out_lun,*)'lower density extent = ',(lower_density_extent(k),k=1,3)
      write(out_lun,*)'upper density extent = ',(upper_density_extent(k),k=1,3)
   endif
end subroutine FH_hermite_density_extent
!----------------------------------------------------------------
!-----------------------------------------------------------------
! EVALUATION SUBROUTINES
!-----------------------------------------------------------------
!-----------------------------------------------------------------
subroutine FH_hermite_local_to_global()

   use sites,only : num_sites, &
                    frame_index_of_site, &
                    frame

   implicit none

   include "mpole_sizes.fh"

   integer j,k,n,np,num,off,max_order,indframe,order,coeff_off,dimxy
   double precision Mpole_xy(MAXMP,MAXMP)

   do n = 1,num_sites
      max_order = max_xform_order_for_site(n)
      if ( max_order > 1 )then
         indframe = frame_index_of_site(n)
         call XFORM_MPOLE_matrix(frame(:,:,indframe),Mpole_xy,max_order)
      endif
      ! first transform the multipoles
      if ( tot_num_mpole_coeffs > 0 )then
         order = mpole_order_for_site(n)
         coeff_off = mpole_coeff_off_for_site(n)
         dimxy = max_order
         ! if order is zero nothing happens
         if ( order == 1 )then ! just copy
            Global_multipole_coeff(coeff_off+1) =  &
            Local_multipole_coeff(coeff_off+1)
         elseif ( order > 1 )then
            call XFORM_MPOLE(Mpole_xy,dimxy, &
                          Local_multipole_coeff(coeff_off+1), &
                          Global_multipole_coeff(coeff_off+1),order)
         endif
      endif !( tot_num_mpole_coeffs > 0 )
      !print *,'MULTIPOLES'
      !print *,'site,global,local'
      !if (order == 1)then
      !   print *,n,Global_multipole_coeff(coeff_off+j),&
      !             Local_multipole_coeff(coeff_off+j)
      !else if (order == 4)then
      !   do j=1,3
      !      print *,n,Global_multipole_coeff(coeff_off+j),&
      !                Local_multipole_coeff(coeff_off+j)
      !   enddo
      !else if (order == 10)then
      !   do j=1,10
      !      print *,n,Global_multipole_coeff(coeff_off+j),&
      !                Local_multipole_coeff(coeff_off+j)
      !   enddo
      !endif

      ! next transform the sumCH if they exist
      if ( tot_num_sumCH_coeffs > 0 )then
         order = sumCH_order_for_site(n)
         coeff_off = sumCH_coeff_off_for_site(n)
         dimxy = max_order
         ! if order is zero nothing happens
         if ( order == 1 )then ! just copy
            Global_sumCH_coeff(coeff_off+1) =  &
            Local_sumCH_coeff(coeff_off+1)
         elseif ( order > 1 )then
            call XFORM_MPOLE(Mpole_xy,dimxy, &
                          Local_sumCH_coeff(coeff_off+1), &
                          Global_sumCH_coeff(coeff_off+1),order)
         endif
      endif !( tot_num_sumCH_coeffs > 0 )

      ! next transform the compact hermites
      if ( tot_num_herm_Cprims > 0 )then
         num = num_herm_Cprim_for_site(n)
         off = off_herm_Cprim_for_site(n)
         do k = 1,num
            np = off+k
            order = hermite_order_for_Cprim(np)
            coeff_off = herm_coeff_offset_of_Cprim(np)
            dimxy = max_order
            if ( order == 1 )then ! just copy
               Global_Chermite_coeff(coeff_off+1) =  &
               Local_Chermite_coeff(coeff_off+1)
            elseif( order > 1 )then
               call XFORM_MPOLE(Mpole_xy,dimxy, &
                             Local_Chermite_coeff(coeff_off+1), &
                             Global_Chermite_coeff(coeff_off+1),order)
            endif
         enddo! k = 1,num
      endif !( tot_num_herm_Cprims > 0 )then
      !print *,'Compact HERMITES'
      !print *,'site,global,local'
      !if (order == 1)then
      !   print *,n,Global_multipole_coeff(coeff_off+j),&
      !             Local_multipole_coeff(coeff_off+j)
      !else if (order == 4)then
      !   do j=1,3
      !      print *,n,Global_multipole_coeff(coeff_off+j),&
      !                Local_multipole_coeff(coeff_off+j)
      !   enddo
      !else if (order == 10)then
      !   do j=1,10
      !      print *,n,Global_multipole_coeff(coeff_off+j),&
      !                Local_multipole_coeff(coeff_off+j)
      !   enddo
      !endif

      ! next transform the diffuse hermites
      if ( tot_num_herm_Dprims > 0 )then
         num = num_herm_Dprim_for_site(n)
         off = off_herm_Dprim_for_site(n)
         do k = 1,num
            np = off+k
            order = hermite_order_for_Dprim(np)
            coeff_off = herm_coeff_offset_of_Dprim(np)
            dimxy = max_order
            if ( order == 1 )then ! just copy
               Global_Dhermite_coeff(coeff_off+1) =  &
               Local_Dhermite_coeff(coeff_off+1)
            elseif( order > 1 )then
               call XFORM_MPOLE(Mpole_xy,dimxy, &
                             Local_Dhermite_coeff(coeff_off+1), &
                             Global_Dhermite_coeff(coeff_off+1),order)
            endif
         enddo! k = 1,num
      endif !( tot_num_herm_Dprims > 0 )then
      !print *,'Compact HERMITES'
      !print *,'site,global,local'
      !if (order == 1)then
      !   print *,n,Global_multipole_coeff(coeff_off+j),&
      !             Local_multipole_coeff(coeff_off+j)
      !else if (order == 4)then
      !   do j=1,3
      !      print *,n,Global_multipole_coeff(coeff_off+j),&
      !                Local_multipole_coeff(coeff_off+j)
      !   enddo
      !else if (order == 10)then
      !   do j=1,10
      !      print *,n,Global_multipole_coeff(coeff_off+j),&
      !                Local_multipole_coeff(coeff_off+j)
      !   enddo
      !endif

   enddo !n = 1,num_sites
   
end subroutine FH_hermite_local_to_global
!----------------------------------------------------------------
subroutine FH_hermite_field_de_drot()

   use sites,only : num_sites,de_drotsite

   implicit none

   include "mpole_sizes.fh"

   double precision DMP_x(MAXMP*MAXMP),DMP_y(MAXMP*MAXMP),  &
                    DMP_z(MAXMP*MAXMP),A_xy(3,3),DA_xy(3,3),  &
                    Tmp_x(MAXMP),Tmp_y(MAXMP),Tmp_z(MAXMP)

   integer :: i,j,k,n,np,num,off,order,coeff_off,dimxy

! to get de_drot we calculate the deriv of mpole wrt infinitesmal
! rotations about x,y and z axis
   do i = 1,3
      do j = 1,3
         A_xy(i,j) = 0.d0
      enddo
      A_xy(i,i) = 1.d0
   enddo
     
! x-axis rotation of dtheta
   do i = 1,3
      do j = 1,3
         DA_xy(i,j) = 0.d0
      enddo
   enddo
   DA_xy(3,2) = 1.d0
   DA_xy(2,3) = -1.d0
! do the maximal order
   call XFORM_MPOLE_deriv_matrix(A_xy,DA_xy,DMP_x,MAXMP)

! y-axis
   do i = 1,3
      do j = 1,3
         DA_xy(i,j) = 0.d0
      enddo
   enddo
   DA_xy(3,1) = -1.d0
   DA_xy(1,3) = 1.d0
   call XFORM_MPOLE_deriv_matrix(A_xy,DA_xy,DMP_y,MAXMP)

! z-axis
   do i = 1,3
      do j = 1,3
         DA_xy(i,j) = 0.d0
      enddo
   enddo
   DA_xy(2,1) = 1.d0
   DA_xy(1,2) = -1.d0
   call XFORM_MPOLE_deriv_matrix(A_xy,DA_xy,DMP_z,MAXMP)
 
   dimxy = MAXMP

   do n = 1,num_sites
      de_drotsite(1,n) = 0.d0
      de_drotsite(2,n) = 0.d0
      de_drotsite(3,n) = 0.d0
   enddo

   ! first do multipoles
   if ( tot_num_mpole_coeffs > 0 )then
      do n = 1,num_sites
         coeff_off = mpole_coeff_off_for_site(n)
         order = mpole_order_for_site(n)
         if ( order > 1 )then ! skip those with no torque
            call XFORM_MPOLE(DMP_x,dimxy, &
                          Global_multipole_coeff(coeff_off+1),Tmp_x,order)
            call XFORM_MPOLE(DMP_y,dimxy, &
                          Global_multipole_coeff(coeff_off+1),Tmp_y,order)
            call XFORM_MPOLE(DMP_z,dimxy, &
                          Global_multipole_coeff(coeff_off+1),Tmp_z,order)
            do j = 1,order
               de_drotsite(1,n) = de_drotsite(1,n) +  &
                               Tmp_x(j)*Global_multipole_field(coeff_off+j)
               de_drotsite(2,n) = de_drotsite(2,n) +  &
                               Tmp_y(j)*Global_multipole_field(coeff_off+j)
               de_drotsite(3,n) = de_drotsite(3,n) +  &
                               Tmp_z(j)*Global_multipole_field(coeff_off+j)
            enddo
         endif!( order > 1 )
      enddo !n = 1,num_sites
   endif !( tot_num_mpole_coeffs > 0 )then

   ! next do sumCH
   if ( tot_num_sumCH_coeffs > 0 )then
      do n = 1,num_sites
         coeff_off = sumCH_coeff_off_for_site(n)
         order = sumCH_order_for_site(n)
         if ( order > 1 )then ! skip those with no torque
            call XFORM_MPOLE(DMP_x,dimxy, &
                          Global_sumCH_coeff(coeff_off+1),Tmp_x,order)
            call XFORM_MPOLE(DMP_y,dimxy, &
                          Global_sumCH_coeff(coeff_off+1),Tmp_y,order)
            call XFORM_MPOLE(DMP_z,dimxy, &
                          Global_sumCH_coeff(coeff_off+1),Tmp_z,order)
            do j = 1,order
               de_drotsite(1,n) = de_drotsite(1,n) +  &
                               Tmp_x(j)*Global_sumCH_field(coeff_off+j)
               de_drotsite(2,n) = de_drotsite(2,n) +  &
                               Tmp_y(j)*Global_sumCH_field(coeff_off+j)
               de_drotsite(3,n) = de_drotsite(3,n) +  &
                               Tmp_z(j)*Global_sumCH_field(coeff_off+j)
            enddo
         endif!( order > 1 )
      enddo !n = 1,num_sites
   endif !( tot_num_mpole_coeffs > 0 )then

   ! next the compact hermites
   if ( tot_num_herm_Cprims > 0 )then
      do n = 1,num_sites
         num = num_herm_Cprim_for_site(n)
         off = off_herm_Cprim_for_site(n)
         do k = 1,num
            np = off+k
            order = hermite_order_for_Cprim(np)
            coeff_off = herm_coeff_offset_of_Cprim(np)
            if ( order > 1 )then ! skip those with no torque
               call XFORM_MPOLE(DMP_x,dimxy, &
                             Global_Chermite_coeff(coeff_off+1),Tmp_x,order)
               call XFORM_MPOLE(DMP_y,dimxy, &
                             Global_Chermite_coeff(coeff_off+1),Tmp_y,order)
               call XFORM_MPOLE(DMP_z,dimxy, &
                             Global_Chermite_coeff(coeff_off+1),Tmp_z,order)
               do j = 1,order
                  de_drotsite(1,n) = de_drotsite(1,n) +  &
                                  Tmp_x(j)*Global_Chermite_field(coeff_off+j)
                  de_drotsite(2,n) = de_drotsite(2,n) +  &
                                  Tmp_y(j)*Global_Chermite_field(coeff_off+j)
                  de_drotsite(3,n) = de_drotsite(3,n) +  &
                                  Tmp_z(j)*Global_Chermite_field(coeff_off+j)
               enddo
            endif!( order > 1 )
         enddo ! k = 1,num
      enddo !n = 1,num_sites
   endif !( tot_num_herm_Cprims > 0 )then

   ! next the diffuse hermites
   if ( tot_num_herm_Dprims > 0 )then
      do n = 1,num_sites
         num = num_herm_Dprim_for_site(n)
         off = off_herm_Dprim_for_site(n)
         do k = 1,num
            np = off+k
            order = hermite_order_for_Dprim(np)
            coeff_off = herm_coeff_offset_of_Dprim(np)
            if ( order > 1 )then ! skip those with no torque
               call XFORM_MPOLE(DMP_x,dimxy, &
                             Global_Dhermite_coeff(coeff_off+1),Tmp_x,order)
               call XFORM_MPOLE(DMP_y,dimxy, &
                             Global_Dhermite_coeff(coeff_off+1),Tmp_y,order)
               call XFORM_MPOLE(DMP_z,dimxy, &
                             Global_Dhermite_coeff(coeff_off+1),Tmp_z,order)
               do j = 1,order
                  de_drotsite(1,n) = de_drotsite(1,n) +  &
                                  Tmp_x(j)*Global_Dhermite_field(coeff_off+j)
                  de_drotsite(2,n) = de_drotsite(2,n) +  &
                                  Tmp_y(j)*Global_Dhermite_field(coeff_off+j)
                  de_drotsite(3,n) = de_drotsite(3,n) +  &
                                  Tmp_z(j)*Global_Dhermite_field(coeff_off+j)
               enddo
            endif!( order > 1 )
         enddo ! k = 1,num
      enddo !n = 1,num_sites
   endif !( tot_num_herm_Cprims > 0 )then

end subroutine FH_hermite_field_de_drot
!----------------------------------------------------------------
subroutine FH_hermite_field_energy(energy)

   use sites,only : num_sites

   implicit none

   double precision, intent(out) :: energy

   integer :: n,np,j,order,off

   energy = 0.d0
   ! first mpoles
   if ( tot_num_mpole_coeffs > 0 )then
      do n = 1,num_sites
         off = mpole_coeff_off_for_site(n)
         order = mpole_order_for_site(n)
         do j = 1,order
            energy = energy +  &
              Global_multipole_coeff(off+j)*Global_multipole_field(off+j)
         enddo
      enddo
   endif

   ! next sumCH
   if ( tot_num_sumCH_coeffs > 0 )then
      do n = 1,num_sites
         off = sumCH_coeff_off_for_site(n)
         order = sumCH_order_for_site(n)
         do j = 1,order
            energy = energy +  &
              Global_sumCH_coeff(off+j)*Global_sumCH_field(off+j)
         enddo
      enddo
   endif

   ! next compact hermites
   if ( tot_num_herm_Cprims > 0 )then
      do np = 1,tot_num_herm_Cprims
         order = hermite_order_for_Cprim(np)
         off = herm_coeff_offset_of_Cprim(np)
         do j = 1,order
            energy = energy +  &
              Global_Chermite_coeff(off+j)*Global_Chermite_field(off+j)
         enddo
      enddo
   endif !( tot_num_herm_Cprims > 0 )then

   ! next diffuse hermites
   if ( tot_num_herm_Dprims > 0 )then
      do np = 1,tot_num_herm_Dprims
         order = hermite_order_for_Dprim(np)
         off = herm_coeff_offset_of_Dprim(np)
         do j = 1,order
            energy = energy +  &
              Global_Dhermite_coeff(off+j)*Global_Dhermite_field(off+j)
         enddo
      enddo
   endif !( tot_num_herm_Dprims > 0 )then

   energy = 0.5d0*energy
end subroutine FH_hermite_field_energy
!----------------------------------------------- -----------------
!----------------------------------------------------------------
subroutine FH_hermite_field_de_drot_EXCH()

   use sites,only : num_sites,de_drotsite

   implicit none

   include "mpole_sizes.fh"

   double precision DMP_x(MAXMP*MAXMP),DMP_y(MAXMP*MAXMP),  &
                    DMP_z(MAXMP*MAXMP),A_xy(3,3),DA_xy(3,3),  &
                    Tmp_x(MAXMP),Tmp_y(MAXMP),Tmp_z(MAXMP)

   integer :: i,j,k,n,np,num,off,order,coeff_off,dimxy

! to get de_drot we calculate the deriv of mpole wrt infinitesmal
! rotations about x,y and z axis
   do i = 1,3
      do j = 1,3
         A_xy(i,j) = 0.d0
      enddo
      A_xy(i,i) = 1.d0
   enddo
     
! x-axis rotation of dtheta
   do i = 1,3
      do j = 1,3
         DA_xy(i,j) = 0.d0
      enddo
   enddo
   DA_xy(3,2) = 1.d0
   DA_xy(2,3) = -1.d0
! do the maximal order
   call XFORM_MPOLE_deriv_matrix(A_xy,DA_xy,DMP_x,MAXMP)

! y-axis
   do i = 1,3
      do j = 1,3
         DA_xy(i,j) = 0.d0
      enddo
   enddo
   DA_xy(3,1) = -1.d0
   DA_xy(1,3) = 1.d0
   call XFORM_MPOLE_deriv_matrix(A_xy,DA_xy,DMP_y,MAXMP)

! z-axis
   do i = 1,3
      do j = 1,3
         DA_xy(i,j) = 0.d0
      enddo
   enddo
   DA_xy(2,1) = 1.d0
   DA_xy(1,2) = -1.d0
   call XFORM_MPOLE_deriv_matrix(A_xy,DA_xy,DMP_z,MAXMP)
 
   dimxy = MAXMP

   do n = 1,num_sites
      de_drotsite(1,n) = 0.d0
      de_drotsite(2,n) = 0.d0
      de_drotsite(3,n) = 0.d0
   enddo

   ! NO multipoles or sumCH for EXCH because only Hermites here!
 
   ! next the compact hermites
   if ( tot_num_herm_Cprims > 0 )then
      do n = 1,num_sites
         num = num_herm_Cprim_for_site(n)
         off = off_herm_Cprim_for_site(n)
         do k = 1,num
            np = off+k
            order = hermite_order_for_Cprim(np)
            coeff_off = herm_coeff_offset_of_Cprim(np)
            if ( order > 1 )then ! skip those with no torque
               call XFORM_MPOLE(DMP_x,dimxy, &
                             Global_Chermite_coeff(coeff_off+1),Tmp_x,order)
               call XFORM_MPOLE(DMP_y,dimxy, &
                             Global_Chermite_coeff(coeff_off+1),Tmp_y,order)
               call XFORM_MPOLE(DMP_z,dimxy, &
                             Global_Chermite_coeff(coeff_off+1),Tmp_z,order)
               do j = 1,order
                  de_drotsite(1,n) = de_drotsite(1,n) +  &
                                  Tmp_x(j)*Global_Chermite_field_EXCH(coeff_off+j)
                  de_drotsite(2,n) = de_drotsite(2,n) +  &
                                  Tmp_y(j)*Global_Chermite_field_EXCH(coeff_off+j)
                  de_drotsite(3,n) = de_drotsite(3,n) +  &
                                  Tmp_z(j)*Global_Chermite_field_EXCH(coeff_off+j)
               enddo
            endif!( order > 1 )
         enddo ! k = 1,num
      enddo !n = 1,num_sites
   endif !( tot_num_herm_Cprims > 0 )then

   ! next the diffuse hermites
   if ( tot_num_herm_Dprims > 0 )then
      do n = 1,num_sites
         num = num_herm_Dprim_for_site(n)
         off = off_herm_Dprim_for_site(n)
         do k = 1,num
            np = off+k
            order = hermite_order_for_Dprim(np)
            coeff_off = herm_coeff_offset_of_Dprim(np)
            if ( order > 1 )then ! skip those with no torque
               call XFORM_MPOLE(DMP_x,dimxy, &
                             Global_Dhermite_coeff(coeff_off+1),Tmp_x,order)
               call XFORM_MPOLE(DMP_y,dimxy, &
                             Global_Dhermite_coeff(coeff_off+1),Tmp_y,order)
               call XFORM_MPOLE(DMP_z,dimxy, &
                             Global_Dhermite_coeff(coeff_off+1),Tmp_z,order)
               do j = 1,order
                  de_drotsite(1,n) = de_drotsite(1,n) +  &
                                  Tmp_x(j)*Global_Dhermite_field_EXCH(coeff_off+j)
                  de_drotsite(2,n) = de_drotsite(2,n) +  &
                                  Tmp_y(j)*Global_Dhermite_field_EXCH(coeff_off+j)
                  de_drotsite(3,n) = de_drotsite(3,n) +  &
                                  Tmp_z(j)*Global_Dhermite_field_EXCH(coeff_off+j)
               enddo
            endif!( order > 1 )
         enddo ! k = 1,num
      enddo !n = 1,num_sites
   endif !( tot_num_herm_Cprims > 0 )then

end subroutine FH_hermite_field_de_drot_EXCH
!----------------------------------------------------------------
subroutine FH_hermite_field_energy_EXCH(energy)

   use sites,only : num_sites

   implicit none

   double precision, intent(out) :: energy

   integer :: n,np,j,order,off

   energy = 0.d0
   ! NO CONTRIBUTION from mpoles or sumCH because only overlap for EXCH

   ! next compact hermites
   if ( tot_num_herm_Cprims > 0 )then
      do np = 1,tot_num_herm_Cprims
         order = hermite_order_for_Cprim(np)
         off = herm_coeff_offset_of_Cprim(np)
         do j = 1,order
            energy = energy +  &
              Global_Chermite_coeff(off+j)*Global_Chermite_field_EXCH(off+j)
         enddo
      enddo
   endif !( tot_num_herm_Cprims > 0 )then

   ! next diffuse hermites
   if ( tot_num_herm_Dprims > 0 )then
      do np = 1,tot_num_herm_Dprims
         order = hermite_order_for_Dprim(np)
         off = herm_coeff_offset_of_Dprim(np)
         do j = 1,order
            energy = energy +  &
              Global_Dhermite_coeff(off+j)*Global_Dhermite_field_EXCH(off+j)
         enddo
      enddo
   endif !( tot_num_herm_Dprims > 0 )then

   energy = 0.5d0*energy
end subroutine FH_hermite_field_energy_EXCH
!----------------------------------------------------------------

!----------------------------------------------------------------
!-----------------------------------------------------------------
! DEALLOCATE SUBROUTINES
!-----------------------------------------------------------------
!-----------------------------------------------------------------

!----------------------------------------------------------------
subroutine FH_hermite_deallocate()

   implicit none

   integer ii
   ii = 1
   if ( allocated(max_xform_order_for_site) ) &
               deallocate(max_xform_order_for_site)
   ! multipole stuff
   if ( allocated(mpole_order_for_site) ) &
               deallocate(mpole_order_for_site)
   if ( allocated(mpole_level_for_site) ) &
               deallocate(mpole_level_for_site)
   if ( allocated(mpole_coeff_off_for_site) ) &
               deallocate(mpole_coeff_off_for_site)
   if ( allocated(mpole_coeff_off_for_site) ) &
               deallocate(mpole_coeff_off_for_site)
   if ( allocated(Local_multipole_coeff) ) &
               deallocate(Local_multipole_coeff)
   if ( allocated(Global_multipole_coeff) ) &
               deallocate(Global_multipole_coeff)
   if ( allocated(Global_multipole_field) ) &
               deallocate(Global_multipole_field)
   ! sum of compact hermites stuff
   if ( allocated(sumCH_order_for_site) ) &
               deallocate(sumCH_order_for_site)
   if ( allocated(sumCH_level_for_site) ) &
               deallocate(sumCH_level_for_site)
   if ( allocated(sumCH_coeff_off_for_site) ) &
               deallocate(sumCH_coeff_off_for_site)
   if ( allocated(Local_sumCH_coeff) ) &
               deallocate(Local_sumCH_coeff)
   if ( allocated(Global_sumCH_coeff) ) &
               deallocate(Global_sumCH_coeff)
   !if ( ii == 1 )stop
   if ( allocated(Global_sumCH_field) ) &
               deallocate(Global_sumCH_field)
   ! compact hermite stuff
   if ( allocated(num_herm_Cprim_for_site) ) &
               deallocate(num_herm_Cprim_for_site)
   if ( allocated(off_herm_Cprim_for_site) ) &
               deallocate(off_herm_Cprim_for_site)
   if ( allocated(hermite_order_for_Cprim) ) &
               deallocate(hermite_order_for_Cprim)
   if ( allocated(hermite_level_for_Cprim) ) &
               deallocate(hermite_level_for_Cprim)
   if ( allocated(site_that_owns_Cprim) ) &
               deallocate(site_that_owns_Cprim)
   if ( allocated(herm_coeff_offset_of_Cprim) ) &
               deallocate(herm_coeff_offset_of_Cprim)
   if ( allocated(HC_grid_number_for_Cprim) ) &
               deallocate(HC_grid_number_for_Cprim)
   if ( allocated(hermite_expon_for_Cprim) ) &
               deallocate(hermite_expon_for_Cprim)
   if ( allocated(hermite_rec_expon_for_Cprim) ) &
               deallocate(hermite_rec_expon_for_Cprim)
   if ( allocated(gauss_extent_for_Cprim) ) &
               deallocate(gauss_extent_for_Cprim)
   if ( allocated(Local_Chermite_coeff) ) &
               deallocate(Local_Chermite_coeff)
   if ( allocated(Global_Chermite_coeff) ) &
               deallocate(Global_Chermite_coeff)
   if ( allocated(Global_Chermite_field) ) &
               deallocate(Global_Chermite_field)
   if ( allocated(Global_Chermite_field_EXCH) ) &
               deallocate(Global_Chermite_field_EXCH)
   ! diffuse hermite stuff
   if ( allocated(num_herm_Dprim_for_site) ) &
               deallocate(num_herm_Dprim_for_site)
   if ( allocated(off_herm_Dprim_for_site) ) &
               deallocate(off_herm_Dprim_for_site)
   if ( allocated(hermite_order_for_Dprim) ) &
               deallocate(hermite_order_for_Dprim)
   if ( allocated(hermite_level_for_Dprim) ) &
               deallocate(hermite_level_for_Dprim)
   if ( allocated(site_that_owns_Dprim) ) &
               deallocate(site_that_owns_Dprim)
   if ( allocated(herm_coeff_offset_of_Dprim) ) &
               deallocate(herm_coeff_offset_of_Dprim)
   if ( allocated(HD_grid_number_for_Dprim) ) &
               deallocate(HD_grid_number_for_Dprim)
   if ( allocated(hermite_expon_for_Dprim) ) &
               deallocate(hermite_expon_for_Dprim)
   if ( allocated(hermite_rec_expon_for_Dprim) ) &
               deallocate(hermite_rec_expon_for_Dprim)
   if ( allocated(gauss_extent_for_Dprim) ) &
               deallocate(gauss_extent_for_Dprim)
   if ( allocated(Local_Dhermite_coeff) ) &
               deallocate(Local_Dhermite_coeff)
   if ( allocated(Global_Dhermite_coeff) ) &
               deallocate(Global_Dhermite_coeff)
   if ( allocated(Global_Dhermite_field) ) &
               deallocate(Global_Dhermite_field)
   if ( allocated(Global_Dhermite_field_EXCH)) &
               deallocate(Global_Dhermite_field_EXCH)
end subroutine FH_hermite_deallocate
!----------------------------------------------------------------
end module hermite
