! build arrays of various attributes of auxiliary primitives with pointers 
! from array of site-types
! build array of hermite coefficients with pointers from auxiliary primitives

module auxiliary

   implicit none

   private

   integer,save :: num_aux_site_types, &
                   tot_num_aux_Cprim, &
                   tot_num_aux_Dprim, &
                   tot_num_aux_Cherm_coeff, &
                   tot_num_aux_Dherm_coeff, &
                   tot_num_aux_mpole_coeff, &
                   tot_num_aux_sumCH_coeff


   integer, save, allocatable ::  &
                                 num_aux_Cprim_for_site_type(:), &
                                 num_aux_Dprim_for_site_type(:), &
                                 off_aux_Cprim_for_site_type(:), &
                                 off_aux_Dprim_for_site_type(:), &
                                 nuclear_charge_for_site_type(:), &
                                 sumCH_order_for_site_type(:), &
                                 sumCH_level_for_site_type(:), &
                                 sumCH_coeff_off_for_site_type(:), &
                                 mpole_order_for_site_type(:), &
                                 mpole_level_for_site_type(:), &
                                 mpole_coeff_off_for_site_type(:)
   double precision, save, allocatable :: aux_mpole_coeff(:)
   double precision, save, allocatable :: aux_sumCH_coeff(:)

   integer, save, allocatable :: herm_ord_for_aux_Cprim(:), &
                                 herm_lev_for_aux_Cprim(:), &
                                 herm_coeff_off_for_aux_Cprim(:)
   double precision, save, allocatable :: herm_exp_for_aux_Cprim(:)
   double precision, save, allocatable :: herm_rec_exp_for_aux_Cprim(:)
   double precision, save, allocatable :: aux_Chermite_coeff(:)

   integer, save, allocatable :: herm_ord_for_aux_Dprim(:), &
                                 herm_lev_for_aux_Dprim(:), &
                                 herm_coeff_off_for_aux_Dprim(:)
   double precision, save, allocatable :: herm_exp_for_aux_Dprim(:)
   double precision, save, allocatable :: herm_rec_exp_for_aux_Dprim(:)
   double precision, save, allocatable :: aux_Dhermite_coeff(:)

   double precision, save, allocatable :: max_extent_for_site_type(:)

   public FH_AUX_read, &
          FH_AUX_deallocate, &
          FH_AUX_get_rec_expon, &
          FH_AUX_recip_auto_setup, &
          num_aux_site_types, &
          max_extent_for_site_type, &
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
          
contains

subroutine FH_AUX_read(filename, free_lun, out_lun)

   use user, only : user_num_site_types, &
                    user_num_aux_Cprims, &
                    user_num_aux_Dprims, &
                    verbose,CD_split_expon, &
                    erfcinv_gauss_extent_tol, &
                    extent_of_compact_hermites, &
                    pme_auto_setup, &
                    ffp_auto_setup, &
                    reg_ewald_auto_setup

   implicit none

! Formal arguments:

   character(len=*), intent(in) :: filename
   integer,          intent(in) :: free_lun
   integer,          intent(in) :: out_lun

! Local variables:

   include "interact_type.fh"
   include "mpole_sizes.fh"
   include "scale.fh"

   character(len=100) :: line
   integer :: ios,ier1,ier2,ier3,ier4,ier5,ier6
   integer :: j,k,m,n,itype,nprim,ncprim,ndprim,nuclear_cg,primtype,iprim
   integer :: np,order,level,max_Corder,max_Clevel,offc,offd,herm_off, &
              mp_order
   double precision :: expon,min_expon,coeff,sum_coeff(MAXMP),mp_coeff

   open(unit=free_lun,file=filename,status='old')
   read(free_lun,'(a)',iostat=ios)line
   if ( ios /= 0 .or. line(1:9) /= "num_types" )then
      write(out_lun,*)'FH_AUX_read: bad line(num_types): ',line(1:25)
      stop
   endif ! ( ios /= 0 )
   read(line(10:),*)num_aux_site_types
   ! if not-auto ffp,pme or ewald 1st check against user info
   if ( (ffp_auto_setup==0) .and. (pme_auto_setup==0) .and. &
        (reg_ewald_auto_setup==0) )then
      if ( user_num_site_types /= num_aux_site_types )then
         write(out_lun,*) &
           'FH_AUX_read: mismatch in number of types with user input!'
         stop
      endif
   endif
   allocate(num_aux_Cprim_for_site_type(num_aux_site_types),stat=ier1)
   allocate(num_aux_Dprim_for_site_type(num_aux_site_types),stat=ier2)
   allocate(off_aux_Cprim_for_site_type(num_aux_site_types),stat=ier3)
   allocate(off_aux_Dprim_for_site_type(num_aux_site_types),stat=ier4)
   allocate(nuclear_charge_for_site_type(num_aux_site_types),stat=ier5)
   allocate(max_extent_for_site_type(num_aux_site_types),stat=ier6)
   if ( ier1 /= 0 .or. ier2 /= 0 .or. ier3 /= 0 .or. ier4 /= 0 .or. &
          ier5 /= 0 .or. ier6 /= 0 )then
      write(out_lun,*)'FH_AUX_read: allocate fails!'
      stop
   endif
   allocate(mpole_order_for_site_type(num_aux_site_types),stat=ier1)
   allocate(mpole_level_for_site_type(num_aux_site_types),stat=ier2)
   allocate(mpole_coeff_off_for_site_type(num_aux_site_types),stat=ier3)
   if ( ier1 /= 0 .or. ier2 /= 0 .or. ier3 /= 0 )then
      write(out_lun,*)'FH_AUX_read: allocate fails!'
      stop
   endif
   if ( CD_split_expon > small_value )then
      allocate(sumCH_order_for_site_type(num_aux_site_types),stat=ier1)
      allocate(sumCH_level_for_site_type(num_aux_site_types),stat=ier2)
      allocate(sumCH_coeff_off_for_site_type(num_aux_site_types),stat=ier3)
   endif
   ! initialize things
   do n = 1,num_aux_site_types
      num_aux_Cprim_for_site_type(n) = 0
      num_aux_Dprim_for_site_type(n) = 0
      nuclear_charge_for_site_type(n) = 0
      max_extent_for_site_type(n) = big_value
      mpole_order_for_site_type(n) = 0
      mpole_level_for_site_type(n) = -1
      if ( CD_split_expon > small_value )then
         sumCH_order_for_site_type(n) = 0
         sumCH_level_for_site_type(n) = -1
      endif
   enddo !n = 1,num_aux_site_types

   do n = 1,num_aux_site_types
      read(free_lun,'(a)',iostat=ios)line
      if ( ios /= 0 .or. line(1:4) /= "type" )then
         write(out_lun,*)'FH_AUX_read: bad line(type): ',line(1:25)
         stop
      endif ! ( ios /= 0 )
      read(free_lun,*,iostat=ios)itype,nprim,nuclear_cg
      if ( ios /= 0 )then
         write(out_lun,*)'FH_AUX_read: bad read in type #',n
         stop
      endif ! ( ios /= 0 )
      if ( n /= itype )then
         write(out_lun,*)'FH_AUX_read: type out of order!'
         stop
      endif
      nuclear_charge_for_site_type(n) = nuclear_cg
      ! read multipole info
      read(free_lun,'(a)',iostat=ios)line
      if ( ios /= 0 .or. line(1:15) /= "multipole_order" )then
         write(out_lun,*)'FH_AUX_read: bad line(multipole_order): ',line(1:25)
         stop
      endif ! ( ios /= 0 )
      read(free_lun,*,iostat=ios)mp_order,mp_coeff
      if ( ios /= 0 )then
         write(out_lun,*)'FH_AUX_read: bad read in type ',n,' mpole header'
         stop
      endif ! ( ios /= 0 )
      mpole_order_for_site_type(n) = mp_order
      if ( mp_order == 0 )then
           mpole_level_for_site_type(n) = -1
      elseif ( mp_order == 1 )then
           mpole_level_for_site_type(n) = 0
      elseif ( mp_order == 4 )then
           mpole_level_for_site_type(n) = 1
      elseif ( mp_order == 10 )then
           mpole_level_for_site_type(n) = 2
      elseif ( mp_order == 20 )then
           mpole_level_for_site_type(n) = 3
      elseif ( mp_order == 35 )then
           mpole_level_for_site_type(n) = 4
      else
           write(out_lun,*) &
             'bad multipole order',mp_order,mpole_level_for_site_type(n)
           stop
      endif
      do j = 2,mp_order
         read(free_lun,*,iostat=ios)mp_coeff
         if ( ios /= 0 )then
            write(out_lun,*)'FH_AUX_read: bad read in type,mpole coeff #',n,j
            stop
         endif ! ( ios /= 0 )
      enddo
      ncprim = 0
      ndprim = 0
      min_expon = big_value
      if ( nprim > 0 )then
         read(free_lun,'(a)',iostat=ios)line
         if ( ios /= 0 .or. line(1:8) /= "prim_num" )then
            write(out_lun,*)'FH_AUX_read: bad line(prim_num): ',line(1:25)
            stop
         endif ! ( ios /= 0 )
         do k = 1,nprim
            read(free_lun,*,iostat=ios)m,expon,order,coeff
            if ( ios /= 0 )then
               write(out_lun,*)'FH_AUX_read: bad read in type,prim #',n,k
               stop
            endif ! ( ios /= 0 )
            if ( expon < min_expon )min_expon = expon
            ! note if CD_split_expon = 0 all prims are compact
            if ( expon >= CD_split_expon )then
               ncprim = ncprim + 1
            else
               ndprim = ndprim + 1
            endif
            do j = 2,order
               read(free_lun,*,iostat=ios)coeff
               if ( ios /= 0 )then
                  write(out_lun,*) &
                    'FH_AUX_read: bad read in type,prim,coeff #',n,k,j
                  stop
               endif ! ( ios /= 0 )
            enddo
         enddo !k = 1,nprim
      endif !nprim > 0
      num_aux_Cprim_for_site_type(n) = ncprim
      num_aux_Dprim_for_site_type(n) = ndprim
      max_extent_for_site_type(n) =  &
              erfcinv_gauss_extent_tol / sqrt(min_expon)
      if ( max_extent_for_site_type(n) < extent_of_compact_hermites ) &
           max_extent_for_site_type(n) = extent_of_compact_hermites
   enddo !n = 1,num_aux_site_types

   ! get offsets
   mpole_coeff_off_for_site_type(1) = 0
   off_aux_Cprim_for_site_type(1) = 0
   off_aux_Dprim_for_site_type(1) = 0
   do n = 2,num_aux_site_types
      off_aux_Cprim_for_site_type(n) = off_aux_Cprim_for_site_type(n-1) + &
                                       num_aux_Cprim_for_site_type(n-1)
      off_aux_Dprim_for_site_type(n) = off_aux_Dprim_for_site_type(n-1) + &
                                       num_aux_Dprim_for_site_type(n-1)
      mpole_coeff_off_for_site_type(n) = mpole_coeff_off_for_site_type(n-1) + &
                                    mpole_order_for_site_type(n-1)
   enddo

   tot_num_aux_mpole_coeff = &
            mpole_coeff_off_for_site_type(num_aux_site_types)+ &
            mpole_order_for_site_type(num_aux_site_types)
   if ( verbose == 1 )then
      write(out_lun,*)'tot_num_aux_mpole_coeff = ',tot_num_aux_mpole_coeff
   endif
   allocate(aux_mpole_coeff(tot_num_aux_mpole_coeff),stat=ier1)
   if ( ier1 /= 0 )then
      write(out_lun,*)'FH_AUX_read: allocate fails!'
      stop
   endif
   tot_num_aux_Cprim = off_aux_Cprim_for_site_type(num_aux_site_types)+ &
                       num_aux_Cprim_for_site_type(num_aux_site_types)
   tot_num_aux_Dprim = off_aux_Dprim_for_site_type(num_aux_site_types)+ &
                       num_aux_Dprim_for_site_type(num_aux_site_types)
   ! again if not auto check against user input
   if ( (ffp_auto_setup==0) .and. (pme_auto_setup==0) .and. &
        (reg_ewald_auto_setup==0) )then
      if ( tot_num_aux_Cprim /= user_num_aux_Cprims .or. &
           tot_num_aux_Dprim /= user_num_aux_Dprims  )then
         write(out_lun,*) &
           'FH_AUX_read: mismatch with user_input in number of auxprims!'
         write(out_lun,*) &
           ' tot_num_auxCprim,user_num_aux_Cprims = ', &
             tot_num_aux_Cprim,user_num_aux_Cprims
         write(out_lun,*) &
           ' tot_num_auxDprim,user_num_aux_Dprims = ', &
             tot_num_aux_Dprim,user_num_aux_Dprims
         stop
      endif
   endif
   !allocate
   if ( tot_num_aux_Cprim > 0 )then
      allocate(herm_exp_for_aux_Cprim(tot_num_aux_Cprim),stat=ier1)
      allocate(herm_rec_exp_for_aux_Cprim(tot_num_aux_Cprim),stat=ier2)
      allocate(herm_ord_for_aux_Cprim(tot_num_aux_Cprim),stat=ier3)
      allocate(herm_lev_for_aux_Cprim(tot_num_aux_Cprim),stat=ier4)
      allocate(herm_coeff_off_for_aux_Cprim(tot_num_aux_Cprim),stat=ier5)
      if ( ier1 /= 0 .or. ier2 /= 0 .or. ier3 /= 0 .or. ier4 /= 0 .or. &
           ier5 /= 0 )then
         write(out_lun,*)'FH_AUX_read: allocate fails!'
         stop
      endif
   endif !( tot_num_aux_Cprim > 0 )then

   if ( tot_num_aux_Dprim > 0 )then
      allocate(herm_exp_for_aux_Dprim(tot_num_aux_Dprim),stat=ier1)
      allocate(herm_rec_exp_for_aux_Dprim(tot_num_aux_Dprim),stat=ier2)
      allocate(herm_ord_for_aux_Dprim(tot_num_aux_Dprim),stat=ier3)
      allocate(herm_lev_for_aux_Dprim(tot_num_aux_Dprim),stat=ier4)
      allocate(herm_coeff_off_for_aux_Dprim(tot_num_aux_Dprim),stat=ier5)
      if ( ier1 /= 0 .or. ier2 /= 0 .or. ier3 /= 0 .or. ier4 /= 0 .or. &
           ier5 /= 0 )then
         write(out_lun,*)'FH_AUX_read: allocate fails!'
         stop
      endif
   endif !( tot_num_aux_Dprim > 0 )then

   ! pass two load prim arrays and get num coeffs
   rewind(free_lun)
   read(free_lun,'(a)',iostat=ios)line
   if ( ios /= 0 .or. line(1:9) /= "num_types" )then
      write(out_lun,*)'FH_AUX_read: bad line(num_types): ',line(1:25)
      stop
   endif ! ( ios /= 0 )
   do n = 1,num_aux_site_types
      read(free_lun,'(a)',iostat=ios)line
      if ( ios /= 0 .or. line(1:4) /= "type" )then
         write(out_lun,*)'FH_AUX_read: bad line(type): ',line(1:25)
         stop
      endif ! ( ios /= 0 )
      read(free_lun,*,iostat=ios)itype,nprim,nuclear_cg
      if ( ios /= 0 )then
         write(out_lun,*)'FH_AUX_read: bad read in type #',n
         stop
      endif ! ( ios /= 0 )
      ! read multipole info
      read(free_lun,'(a)',iostat=ios)line
      if ( ios /= 0 .or. line(1:15) /= "multipole_order" )then
         write(out_lun,*)'FH_AUX_read: bad line(multipole_order): ',line(1:25)
         stop
      endif ! ( ios /= 0 )
      read(free_lun,*,iostat=ios)mp_order,mp_coeff
      do j = 2,mp_order
         read(free_lun,*,iostat=ios)mp_coeff
         if ( ios /= 0 )then
            write(out_lun,*)'FH_AUX_read: bad read in type,mpole coeff #',n,j
            stop
         endif ! ( ios /= 0 )
      enddo
      if ( ios /= 0 )then
         write(out_lun,*)'FH_AUX_read: bad read in type ',n,' mpole header'
         stop
      endif ! ( ios /= 0 )
      ncprim = 0
      ndprim = 0
      if ( nprim > 0 )then
         read(free_lun,'(a)',iostat=ios)line
         if ( ios /= 0 .or. line(1:8) /= "prim_num" )then
            write(out_lun,*)'FH_AUX_read: bad line(prim_num): ',line(1:25)
            stop
         endif ! ( ios /= 0 )
         max_Corder = -1
         max_Clevel = -1
         offc = off_aux_Cprim_for_site_type(n)
         offd = off_aux_Dprim_for_site_type(n)
         do k = 1,nprim
            read(free_lun,*,iostat=ios)m,expon,order,coeff
            if ( ios /= 0 )then
               write(out_lun,*)'FH_AUX_read: bad read in type #',n
               stop
            endif ! ( ios /= 0 )
            if ( order == 1 )then
               level = 0
            elseif ( order == 4 )then
               level = 1
            elseif ( order == 10 )then
               level = 2
            elseif ( order == 20 )then
               level = 3
            elseif ( order == 35 )then
               level = 4
            else
               write(out_lun,*)'bad hermite order'
               stop
            endif
            if ( expon >= CD_split_expon )then
               ncprim = ncprim + 1
               np = offc + ncprim
               herm_exp_for_aux_Cprim(np) = expon
               herm_ord_for_aux_Cprim(np) = order
               herm_lev_for_aux_Cprim(np) = level
               if ( CD_split_expon > small_value ) then
                  if ( order > max_Corder )max_Corder = order
                  if ( level > max_Clevel )max_Clevel = level
               endif
            else
               ndprim = ndprim + 1
               np = offd + ndprim
               herm_exp_for_aux_Dprim(np) = expon
               herm_ord_for_aux_Dprim(np) = order
               herm_lev_for_aux_Dprim(np) = level
            endif
            do j = 2,order
               read(free_lun,*,iostat=ios)coeff
               if ( ios /= 0 )then
                  write(out_lun,*)'FH_AUX_read: bad read in type #',n
                  stop
               endif ! ( ios /= 0 )
            enddo
         enddo !k = 1,nprim
         if ( CD_split_expon > small_value ) then
            sumCH_order_for_site_type(n) = max_Corder
            sumCH_level_for_site_type(n) = max_Clevel
         endif
      endif !nprim > 0
   enddo !n = 1,num_aux_site_types
   ! get sumCH offsets
   if ( CD_split_expon > small_value ) then
      sumCH_coeff_off_for_site_type(1) = 0
      do n = 2,num_aux_site_types
         sumCH_coeff_off_for_site_type(n) =  &
                                    sumCH_coeff_off_for_site_type(n-1) + &
                                    sumCH_order_for_site_type(n-1)
      enddo
      tot_num_aux_sumCH_coeff = &
            sumCH_coeff_off_for_site_type(num_aux_site_types)+ &
            sumCH_order_for_site_type(num_aux_site_types)
      if ( verbose == 1 )then
         write(out_lun,*)'tot_num_aux_sumCH_coeff = ',tot_num_aux_sumCH_coeff
      endif
      allocate(aux_sumCH_coeff(tot_num_aux_sumCH_coeff),stat=ier1)
      if ( ier1 /= 0 )then
         write(out_lun,*)'FH_AUX_read: allocate fails!'
         stop
      endif
   else
      tot_num_aux_sumCH_coeff = 0
   endif

   ! get hermite offsets
   if ( tot_num_aux_Cprim > 0 )then
     herm_coeff_off_for_aux_Cprim(1) = 0
      do np = 2,tot_num_aux_Cprim
         herm_coeff_off_for_aux_Cprim(np) =  &
                                       herm_coeff_off_for_aux_Cprim(np-1) + &
                                       herm_ord_for_aux_Cprim(np-1)
      enddo
      tot_num_aux_Cherm_coeff =  &
                           herm_coeff_off_for_aux_Cprim(tot_num_aux_Cprim) + &
                           herm_ord_for_aux_Cprim(tot_num_aux_Cprim)
   else
      tot_num_aux_Cherm_coeff = 0
   endif !( tot_num_aux_Cprim > 0 )then
   if ( verbose == 1 )then
      write(out_lun,*)'tot_num_aux_Cherm_coeff = ',tot_num_aux_Cherm_coeff
   endif
   ! allocate coeff array
   if ( tot_num_aux_Cprim > 0 )then
      allocate(aux_Chermite_coeff(tot_num_aux_Cherm_coeff),stat=ier1)
      if ( ier1 /= 0 )then
         write(out_lun,*)'FH_AUX_read: bad allocate of herm coeff array'
         stop
      endif
   endif !( tot_num_aux_Cprim > 0 )then
   if ( tot_num_aux_Dprim > 0 )then
      herm_coeff_off_for_aux_Dprim(1) = 0
      do np = 2,tot_num_aux_Dprim
         herm_coeff_off_for_aux_Dprim(np) =  &
                                       herm_coeff_off_for_aux_Dprim(np-1) + &
                                       herm_ord_for_aux_Dprim(np-1)
      enddo
      tot_num_aux_Dherm_coeff =  &
                           herm_coeff_off_for_aux_Dprim(tot_num_aux_Dprim) + &
                           herm_ord_for_aux_Dprim(tot_num_aux_Dprim)
   else
      tot_num_aux_Dherm_coeff = 0
   endif !( tot_num_aux_Dprim > 0 )then
   if ( verbose == 1 )then
      write(out_lun,*)'tot_num_aux_Dherm_coeff = ',tot_num_aux_Dherm_coeff
   endif
   ! allocate coeff array
   if ( tot_num_aux_Dprim > 0 )then
      allocate(aux_Dhermite_coeff(tot_num_aux_Dherm_coeff),stat=ier1)
      if ( ier1 /= 0 )then
         write(out_lun,*)'FH_AUX_read: bad allocate of herm coeff array'
         stop
      endif
   endif !( tot_num_aux_Dprim > 0 )then

   ! pass 3 fill coeff array
   rewind(free_lun)
   read(free_lun,'(a)',iostat=ios)line
   if ( ios /= 0 .or. line(1:9) /= "num_types" )then
      write(out_lun,*)'FH_AUX_read: bad line(num_types): ',line(1:25)
      stop
   endif ! ( ios /= 0 )
   do n = 1,num_aux_site_types
      read(free_lun,'(a)',iostat=ios)line
      if ( ios /= 0 .or. line(1:4) /= "type" )then
         write(out_lun,*)'FH_AUX_read: bad line(type): ',line(1:25)
         stop
      endif ! ( ios /= 0 )
      read(free_lun,*,iostat=ios)itype,nprim,nuclear_cg
      if ( ios /= 0 )then
         write(out_lun,*)'FH_AUX_read: bad read in type #',n
         stop
      endif ! ( ios /= 0 )
      ! read multipole info
      read(free_lun,'(a)',iostat=ios)line
      if ( ios /= 0 .or. line(1:15) /= "multipole_order" )then
         write(out_lun,*)'FH_AUX_read: bad line(multipole_order): ',line(1:25)
         stop
      endif ! ( ios /= 0 )
      read(free_lun,*,iostat=ios)mp_order,mp_coeff
      offc = mpole_coeff_off_for_site_type(n)
      if ( mp_order > 0 )then
         aux_mpole_coeff(offc+1) = mp_coeff
      endif
      do j = 2,mp_order
         read(free_lun,*,iostat=ios)mp_coeff
         if ( ios /= 0 )then
            write(out_lun,*)'FH_AUX_read: bad read in type,mpole coeff #',n,j
            stop
         endif ! ( ios /= 0 )
         aux_mpole_coeff(offc+j) = mp_coeff
      enddo
      ncprim = 0
      ndprim = 0
      if ( CD_split_expon > small_value ) then
         !init the compact sum coeffs
         do j = 1,sumCH_order_for_site_type(n)
            sum_coeff(j) = 0.d0
         enddo
      endif
      if ( nprim > 0 )then
      read(free_lun,'(a)',iostat=ios)line
         if ( ios /= 0 .or. line(1:8) /= "prim_num" )then
            write(out_lun,*)'FH_AUX_read: bad line(prim_num): ',line(1:25)
            stop
         endif ! ( ios /= 0 )
         offc = off_aux_Cprim_for_site_type(n)
         offd = off_aux_Dprim_for_site_type(n)
         do k = 1,nprim
            read(free_lun,*,iostat=ios)m,expon,order,coeff
            if ( ios /= 0 .or. m /= k )then
               write(out_lun,*)'FH_AUX_read: bad read in type #',n
               stop
            endif ! ( ios /= 0 )
            if ( expon >= CD_split_expon )then
               ncprim = ncprim + 1
               np = offc + ncprim
               herm_off = herm_coeff_off_for_aux_Cprim(np)
               aux_Chermite_coeff(herm_off+1) = coeff
               do j = 2,order
                  read(free_lun,*,iostat=ios)coeff
                  if ( ios /= 0 )then
                     write(out_lun,*)'FH_AUX_read: bad read in type #',n
                     stop
                  endif ! ( ios /= 0 )
                  aux_Chermite_coeff(herm_off+j) = coeff
               enddo
               if ( CD_split_expon > small_value ) then
                  do j = 1,order
                     sum_coeff(j) = sum_coeff(j) +  &
                                    aux_Chermite_coeff(herm_off+j)
                  enddo
               endif
            else
               ndprim = ndprim + 1
               np = offd + ndprim
               herm_off = herm_coeff_off_for_aux_Dprim(np)
               aux_Dhermite_coeff(herm_off+1) = coeff
               do j = 2,order
                  read(free_lun,*,iostat=ios)coeff
                  if ( ios /= 0 )then
                     write(out_lun,*)'FH_AUX_read: bad read in type #',n
                     stop
                  endif ! ( ios /= 0 )
                  aux_Dhermite_coeff(herm_off+j) = coeff
               enddo
            endif
         enddo !k = 1,nprim
         ! fill in the sumCH coeffs
         if ( CD_split_expon > small_value ) then
            herm_off = sumCH_coeff_off_for_site_type(n)
            do j = 1,sumCH_order_for_site_type(n)
               aux_sumCH_coeff(herm_off+j) = sum_coeff(j)
            enddo
         endif
      endif !nprim > 0
   enddo !n = 1,num_aux_site_types

   close(free_lun)

end subroutine FH_AUX_read
!--------------------------------------------------------------------
subroutine FH_AUX_recip_auto_setup(out_lun)

   use user, only : Coulomb_use_recip, &
                    pme_auto_setup, &
                    ffp_auto_setup, &
                    reg_ewald_auto_setup, &
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
                    num_pme_grid_types, &
                    num_ffp_grid_types, &
                    num_ewald_grid_types, &
                    CD_split_expon, &
                    verbose

   implicit none

! Formal arguments:

   integer, intent(in)  :: out_lun

! Local variables:

   include 'structure_factor_type.fh'
   include 'interact_type.fh'
   integer off,num,kp,np,order,num_cs,num_cbuf,num_dbuf, &
           new,g,num_multiple,itype
   double precision expon,compact_expon_buf(tot_num_aux_Cprim), &
                          diffuse_expon_buf(tot_num_aux_Dprim)

   if ( Coulomb_use_recip == 0 )return !nothing to do here

   if ( pme_auto_setup + ffp_auto_setup + reg_ewald_auto_setup > 1 )then
      write(out_lun,*)'more than one recip auto setup is chosen!'
      stop
   endif

   if ( pme_auto_setup == 1 )then
      do np = 1,tot_num_aux_Cprim
            HC_grid_number_for_Caux(np) = 0
      enddo
      do np = 1,tot_num_aux_Dprim
            HD_grid_number_for_Daux(np) = 0
      enddo
      ! first get the compact non-ewald counterion s-gaussians on each site
      num_multiple = 0
      do itype = 1,num_aux_site_types
         num = num_aux_Cprim_for_site_type(itype)
         off = off_aux_Cprim_for_site_type(itype)
         ! first pass through itypes aux_herms to see if any compact s-gauss
         num_cs = 0
         do kp = 1,num
            np = off + kp
            if ( herm_ord_for_aux_Cprim(np)==1 )then
               num_cs = num_cs + 1
            endif
         enddo
         if ( num_cs > 0 )then
            num_multiple = num_multiple + 1
            do kp = 1,num
               np = off + kp
               if ( herm_ord_for_aux_Cprim(np)==1 )then
                  HC_grid_number_for_Caux(np) = num_multiple
               endif
            enddo
         endif
      enddo
      ! next get the remaining grids
      num_cbuf = 0
      do np = 1,tot_num_aux_Cprim
         if ( HC_grid_number_for_Caux(np) == 0 )then !not already assigned
            expon = herm_exp_for_aux_Cprim(np)
            new = 1
            do kp = 1,num_cbuf
               if ( abs(expon-compact_expon_buf(kp)) < 1.d-5 )then
                  new = 0
                  HC_grid_number_for_Caux(np) = num_multiple+kp
               endif 
            enddo
            if ( new == 1 )then
               num_cbuf = num_cbuf + 1
               compact_expon_buf(num_cbuf) = expon
               HC_grid_number_for_Caux(np) = num_multiple+num_cbuf
            endif
         endif !( HC_grid_number_for_aux(np) == 0 )
      enddo 
      num_dbuf = 0
      do np = 1,tot_num_aux_Dprim
         expon = herm_exp_for_aux_Dprim(np)
         new = 1
         do kp = 1,num_dbuf
            if ( abs(expon-diffuse_expon_buf(kp)) < 1.d-5 )then
               new = 0
               HD_grid_number_for_Daux(np) = kp
            endif 
         enddo
         if ( new == 1 )then
            num_dbuf = num_dbuf + 1
            diffuse_expon_buf(num_dbuf) = expon
            HD_grid_number_for_Daux(np) = num_dbuf
         endif
      enddo 
      num_HC_prim_grids = num_multiple+num_cbuf
      num_HD_prim_grids = num_dbuf
      num_prim_grid_types = 1
      do g = 1,num_HC_prim_grids
         grid_type_for_HC_prim_grid(g) = 1
      enddo
      do g = 1,num_HD_prim_grids
         grid_type_for_HD_prim_grid(g) = 1
      enddo
      grid_type_for_MPOLES = 1
      if ( num_HC_prim_grids > 0 )then
         grid_type_for_sumCH = 1
      else
         grid_type_for_sumCH = 0
      endif
      struc_fac_method_for_grid_type(1) = SF_PME
   endif

   if ( ffp_auto_setup == 1 )then
      num_prim_grid_types = 1
      grid_type_for_MPOLES = 1
      struc_fac_method_for_grid_type(1) = SF_FFP
      if ( tot_num_aux_Cprim > 0 )then
         num_HC_prim_grids = 1
         grid_type_for_HC_prim_grid(1) = 1
         grid_type_for_sumCH = 1
         do np = 1,tot_num_aux_Cprim
            HC_grid_number_for_Caux(np) = 1
         enddo
      else
         num_HC_prim_grids = 0
         grid_type_for_sumCH = 0
      endif
      if ( tot_num_aux_Dprim > 0 )then
         num_HD_prim_grids = 1
         grid_type_for_HD_prim_grid(1) = 1
         do np = 1,tot_num_aux_Dprim
            HD_grid_number_for_Daux(np) = 1
         enddo
      else
         num_HD_prim_grids = 0
      endif
   endif

   if ( reg_ewald_auto_setup == 1 )then
      num_prim_grid_types = 1
      grid_type_for_MPOLES = 1
      struc_fac_method_for_grid_type(1) = SF_EWALD
      if ( tot_num_aux_Cprim > 0 )then
         num_HC_prim_grids = 1
         grid_type_for_HC_prim_grid(1) = 1
         grid_type_for_sumCH = 1
         do np = 1,tot_num_aux_Cprim
            HC_grid_number_for_Caux(np) = 1
         enddo
      else
         num_HC_prim_grids = 0
         grid_type_for_sumCH = 0
      endif
      if ( tot_num_aux_Dprim > 0 )then
         num_HD_prim_grids = 1
         grid_type_for_HD_prim_grid(1) = 1
         do np = 1,tot_num_aux_Dprim
            HD_grid_number_for_Daux(np) = 1
         enddo
      else
         num_HD_prim_grids = 0
      endif
   endif
   if ( verbose == 1 )then
      write(out_lun,*)'tot_num_aux_Cprim = ',tot_num_aux_Cprim
      write(out_lun,*)'tot_num_aux_Dprim = ',tot_num_aux_Dprim
      write(out_lun,*)'num_HC_prim_grids,num_HD_prim_grids = ', &
                 num_HC_prim_grids,num_HD_prim_grids
      write(out_lun,'(a,30i2)')'HC_grid_number_for_Caux = ', &
                (HC_grid_number_for_Caux(np),np=1,tot_num_aux_Cprim)
      write(out_lun,'(a,30i2)')'HD_grid_number_for_Daux = ', &
                (HD_grid_number_for_Daux(np),np=1,tot_num_aux_Dprim)
   endif

end subroutine FH_AUX_recip_auto_setup
!--------------------------------------------------------------------
subroutine FH_AUX_get_rec_expon(out_lun)

   use user, only : diffuse_adjust_exponent,do_coulomb,do_overlap, &
                    Coulomb_use_recip,Overlap_use_recip
   implicit none

! Formal arguments:

   integer, intent(in)  :: out_lun

! Local variables:

   integer np
   double precision expon
 
   do np = 1,tot_num_aux_Cprim
      herm_rec_exp_for_aux_Cprim(np) = herm_exp_for_aux_Cprim(np)
   enddo
   do np = 1,tot_num_aux_Dprim
      herm_rec_exp_for_aux_Dprim(np) = herm_exp_for_aux_Dprim(np)
   enddo

   if ( diffuse_adjust_exponent > 0.d0 )then
      ! make the compact more diffuse
      do np = 1,tot_num_aux_Cprim
         expon = herm_rec_exp_for_aux_Cprim(np)
         expon = (expon*diffuse_adjust_exponent) / &
                                 (diffuse_adjust_exponent + expon)
         herm_rec_exp_for_aux_Cprim(np) = expon
      enddo
      ! make the diffuse corespondingly less diffuse
      do np = 1,tot_num_aux_Dprim
         expon = herm_rec_exp_for_aux_Dprim(np)
         if ( expon >= diffuse_adjust_exponent )then
             write(out_lun,*) &
                'non-compact expon bigger than diffuse_adjust_exponent!'
             write(out_lun,*) &
                'for diffuse prim ',np,' with expon = ',expon
             write(out_lun,*)'do_coulomb,Coulomb_use_recip = ', &
                 do_coulomb,Coulomb_use_recip
             write(out_lun,*)'do_overlap,Overlap_use_recip = ', &
                  do_overlap,Overlap_use_recip
             stop
         endif
         expon = (expon*diffuse_adjust_exponent) / &
                           (diffuse_adjust_exponent - expon)
         herm_rec_exp_for_aux_Dprim(np) = expon
      enddo
   endif

end subroutine FH_AUX_get_rec_expon
!--------------------------------------------------------------------
subroutine FH_AUX_deallocate()

   implicit none

   if ( allocated(max_extent_for_site_type) ) &
               deallocate(max_extent_for_site_type)

   if ( allocated(num_aux_Cprim_for_site_type) ) &
               deallocate(num_aux_Cprim_for_site_type)
   if ( allocated(num_aux_Dprim_for_site_type) ) &
               deallocate(num_aux_Dprim_for_site_type)
   if ( allocated(off_aux_Cprim_for_site_type) ) &
               deallocate(off_aux_Cprim_for_site_type)
   if ( allocated(off_aux_Dprim_for_site_type) ) &
               deallocate(off_aux_Dprim_for_site_type)
   if ( allocated(nuclear_charge_for_site_type) ) &
               deallocate(nuclear_charge_for_site_type)
   if ( allocated(sumCH_order_for_site_type) ) &
               deallocate(sumCH_order_for_site_type)
   if ( allocated(sumCH_level_for_site_type) ) &
               deallocate(sumCH_level_for_site_type)
   if ( allocated(sumCH_coeff_off_for_site_type) ) &
               deallocate(sumCH_coeff_off_for_site_type)
   if ( allocated(mpole_order_for_site_type) ) &
               deallocate(mpole_order_for_site_type)
   if ( allocated(mpole_level_for_site_type) ) &
               deallocate(mpole_level_for_site_type)
   if ( allocated(mpole_coeff_off_for_site_type) ) &
               deallocate(mpole_coeff_off_for_site_type)
   if ( allocated(aux_mpole_coeff) ) deallocate(aux_mpole_coeff)
   if ( allocated(aux_sumCH_coeff) ) deallocate(aux_sumCH_coeff)

   if ( allocated(herm_ord_for_aux_Cprim) ) &
               deallocate(herm_ord_for_aux_Cprim)
   if ( allocated(herm_lev_for_aux_Cprim) ) &
               deallocate(herm_lev_for_aux_Cprim)
   if ( allocated(herm_coeff_off_for_aux_Cprim) ) &
               deallocate(herm_coeff_off_for_aux_Cprim)
   if ( allocated(herm_exp_for_aux_Cprim) ) &
               deallocate(herm_exp_for_aux_Cprim)
   if ( allocated(herm_rec_exp_for_aux_Cprim) ) &
               deallocate(herm_rec_exp_for_aux_Cprim)
   if ( allocated(aux_Chermite_coeff) ) deallocate(aux_Chermite_coeff)

   if ( allocated(herm_ord_for_aux_Dprim) ) &
               deallocate(herm_ord_for_aux_Dprim)
   if ( allocated(herm_lev_for_aux_Dprim) ) &
               deallocate(herm_lev_for_aux_Dprim)
   if ( allocated(herm_coeff_off_for_aux_Dprim) ) &
               deallocate(herm_coeff_off_for_aux_Dprim)
   if ( allocated(herm_exp_for_aux_Dprim) ) &
               deallocate(herm_exp_for_aux_Dprim)
   if ( allocated(herm_rec_exp_for_aux_Dprim) ) &
               deallocate(herm_rec_exp_for_aux_Dprim)
   if ( allocated(aux_Dhermite_coeff) ) deallocate(aux_Dhermite_coeff)
end subroutine FH_AUX_deallocate
!--------------------------------------------------------------------
end module auxiliary
