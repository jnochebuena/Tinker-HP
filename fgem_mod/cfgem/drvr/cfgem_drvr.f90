program cfgem_drvr

   use cfgem_libmod

   implicit none

   integer              :: iargc
   integer              :: argc, arg
   character(len=256)   :: argv
   character(len=80)    :: aux_file, crd_file, list_file
   character(len=80)    :: user_coul_file, user_exch_file

   double precision, allocatable, save  :: atom_crds(:,:) ! in angstrom

   double precision     :: x, y, z      ! tmp coord input

   double precision     :: tot_mpole_mpole_ene, &
                           tot_mpole_herm_ene, &
                           tot_herm_herm_ene

   integer              :: i,j,k
   integer              :: atom_cnt
   character*5          :: atom_name
   integer              :: mol_num      ! trash input..
   integer              :: io_stat
   integer              :: alloc_failed


   ! If necessary, move the lun to something unused.  Should not present
   ! problems here..

   integer, parameter   :: crd_file_lun = 10

   integer, parameter   :: free_lun = 11
   integer, parameter   :: out_lun = 6

   ! Energy and force outputs of cfgem()

   double precision     :: cfgem_tot_coul_ene ! ALL coulomb

   double precision     :: cfgem_coul_ene_fld ! (for torque)

   type(cfgem_ene_rec)  :: cfgem_coul_recip_ene
   type(cfgem_ene_rec)  :: cfgem_coul_direct_ene
   type(cfgem_ene_rec)  :: cfgem_coul_adj_recip_ene
   type(cfgem_ene_rec)  :: cfgem_coul_self_ene
   type(cfgem_ene_rec)  :: cfgem_coul_wigner_ene

   double precision     :: cfgem_exch_ch_ch_ene

   double precision     :: cfgem_exch_ene_fld ! (for torque)

   double precision, allocatable, save :: cfgem_coul_frc(:,:)
   double precision, allocatable, save :: cfgem_exch_frc(:,:)

   aux_file = ""
   list_file = ""
   user_coul_file = ""
   user_exch_file = ""

   arg = 1
   argc = iargc()

   do while (arg <= argc)

      call getarg(arg,argv)

      if (argv == '-aux') then
         arg = arg + 1
         call getarg(arg,aux_file)
      else if (argv == '-crd') then
         arg = arg + 1
         call getarg(arg,crd_file)
      else if (argv == '-list') then
         arg = arg + 1
         call getarg(arg,list_file)
      else if (argv == '-user_coul') then
         arg = arg + 1
         call getarg(arg,user_coul_file)
      else if (argv == '-user_exch') then
         arg = arg + 1
         call getarg(arg,user_exch_file)
      endif

      arg = arg + 1

   end do

   if (aux_file == '') then
      write(out_lun,*)'Auxiliary file not specified'
      stop
   endif  

   if (crd_file == '') then
      write(out_lun,*)'Coordinate file not specified'
      stop
   endif  

   if (list_file == '') then
      write(out_lun,*)'List file not specified'
      stop
   endif  

   if (user_coul_file == '') then
      write(out_lun,*)'User coulomb file not specified'
      stop
   endif  

   if (user_exch_file == '') then
      write(out_lun,*)'User exchange file not specified'
      stop
   endif  

! The coordinate reading stuff, inlined for simplicity, since only used here:

   open(unit=crd_file_lun, file=crd_file, status='old')

   read(crd_file_lun, *, iostat=io_stat) atom_cnt

   if (io_stat .ne. 0) then
      write(out_lun, *)'Bad coordinate file name', crd_file
      stop
   end if

   allocate(atom_crds(3,atom_cnt), &
            cfgem_coul_frc(3,atom_cnt), &
            cfgem_exch_frc(3,atom_cnt), &
            stat=alloc_failed)

   if (alloc_failed .ne. 0) then
     write(out_lun, *)'Allocation of arrays failed!'
     stop
   end if

   do i = 1, atom_cnt

      read(crd_file_lun, *, iostat=io_stat) atom_name, mol_num, x, y, z

      if (io_stat .ne. 0) then
         write(out_lun,*)'Error reading crds for atom number ', i
         stop
      end if

      atom_crds(1,i) = x
      atom_crds(2,i) = y
      atom_crds(3,i) = z

   end do

   close(crd_file_lun)

! End coordinate reading stuff

   ! Atom coordinate input is in angstrom...

   call cfgem_coul(aux_file, list_file, user_coul_file, &
                   atom_cnt, atom_crds, free_lun, out_lun, &
                   cfgem_tot_coul_ene, &
                   cfgem_coul_ene_fld, &
                   cfgem_coul_recip_ene, &
                   cfgem_coul_direct_ene, &
                   cfgem_coul_adj_recip_ene, &
                   cfgem_coul_self_ene, &
                   cfgem_coul_wigner_ene, &
                   cfgem_coul_frc)

   tot_mpole_mpole_ene = cfgem_coul_direct_ene%mp_mp + &
                         cfgem_coul_recip_ene%mp_mp + &
                         cfgem_coul_wigner_ene%mp_mp + &
                         cfgem_coul_adj_recip_ene%mp_mp + &
                         cfgem_coul_self_ene%mp_mp

   tot_mpole_herm_ene  = cfgem_coul_direct_ene%mp_ch + &
                         cfgem_coul_direct_ene%ch_mp + &
                         cfgem_coul_recip_ene%mp_ch + &
                         cfgem_coul_recip_ene%mp_dh + &
                         cfgem_coul_recip_ene%ch_mp + &
                         cfgem_coul_recip_ene%dh_mp + &
                         cfgem_coul_wigner_ene%mp_ch + &
                         cfgem_coul_wigner_ene%mp_dh + &
                         cfgem_coul_wigner_ene%ch_mp + &
                         cfgem_coul_wigner_ene%dh_mp + &
                         cfgem_coul_adj_recip_ene%mp_ch + &
                         cfgem_coul_adj_recip_ene%mp_dh + &
                         cfgem_coul_adj_recip_ene%ch_mp + &
                         cfgem_coul_adj_recip_ene%dh_mp + &
                         cfgem_coul_self_ene%mp_ch + &
                         cfgem_coul_self_ene%mp_dh + &
                         cfgem_coul_self_ene%ch_mp + &
                         cfgem_coul_self_ene%dh_mp

   tot_herm_herm_ene =   cfgem_coul_direct_ene%ch_ch + &
                         cfgem_coul_recip_ene%ch_ch + &
                         cfgem_coul_recip_ene%dh_ch + &
                         cfgem_coul_recip_ene%ch_dh + &
                         cfgem_coul_recip_ene%dh_dh + &
                         cfgem_coul_wigner_ene%ch_ch + &
                         cfgem_coul_wigner_ene%ch_dh + &
                         cfgem_coul_wigner_ene%dh_ch + &
                         cfgem_coul_wigner_ene%dh_dh + &
                         cfgem_coul_adj_recip_ene%ch_ch + &
                         cfgem_coul_adj_recip_ene%ch_dh + &
                         cfgem_coul_adj_recip_ene%dh_ch + &
                         cfgem_coul_adj_recip_ene%dh_dh + &
                         cfgem_coul_self_ene%ch_ch + &
                         cfgem_coul_self_ene%ch_dh + &
                         cfgem_coul_self_ene%dh_ch + &
                         cfgem_coul_self_ene%dh_dh

   write(out_lun,*)'----------------------------------------'
   write(out_lun,*)'tot_mpole_mpole_ene = ',tot_mpole_mpole_ene
   write(out_lun,*)'tot_mpole_herm_ene = ',tot_mpole_herm_ene
   write(out_lun,*)'tot_herm_herm_ene = ',tot_herm_herm_ene
   write(out_lun,*)'tot_coulomb_ene = ',cfgem_tot_coul_ene,' in kcals = ', &
                627.509474277d0*cfgem_tot_coul_ene
   write(out_lun,*)

   write(out_lun,*)'field ene is ', cfgem_coul_ene_fld,  &
         'in kcals ',627.509474277d0 * cfgem_coul_ene_fld

!---
   write(out_lun,*)'----------------------------------------'
   write(out_lun,*)'mpole-mpole components: '

   write(out_lun,*)'cutoff_mpole_mpole_ene = ', &
                   cfgem_coul_direct_ene%mp_mp

   write(out_lun,*)'recip_mpole_mpole_ene = ', &
                   cfgem_coul_recip_ene%mp_mp

   write(out_lun,*)'wigner_mpole_mpole_ene = ', &
                   cfgem_coul_wigner_ene%mp_mp

   write(out_lun,*)'correct_recip_mpole_mpole_ene = ', &
                   cfgem_coul_adj_recip_ene%mp_mp

   write(out_lun,*)'self_mpole_mpole_ene = ', &
                   cfgem_coul_self_ene%mp_mp

   write(out_lun,*)'----------------------------------------'

   write(out_lun,*)'mpole-herm components: '

   write(out_lun,*)'cutoff_mpole_herm_ene = ', &
                   cfgem_coul_direct_ene%mp_ch + &
                   cfgem_coul_direct_ene%ch_mp + &
                   cfgem_coul_direct_ene%mp_dh + &
                   cfgem_coul_direct_ene%dh_mp

   write(out_lun,*)'recip_mpole_herm_ene = ', &
                   cfgem_coul_recip_ene%mp_ch + &
                   cfgem_coul_recip_ene%ch_mp + &
                   cfgem_coul_recip_ene%mp_dh + &
                   cfgem_coul_recip_ene%dh_mp

   write(out_lun,*)'wigner_mpole_herm_ene = ', &
                   cfgem_coul_wigner_ene%mp_ch + &
                   cfgem_coul_wigner_ene%ch_mp + &
                   cfgem_coul_wigner_ene%mp_dh + &
                   cfgem_coul_wigner_ene%dh_mp

   write(out_lun,*)'correct_recip_mpole_herm_ene = ', &
                   cfgem_coul_adj_recip_ene%mp_ch + &
                   cfgem_coul_adj_recip_ene%ch_mp + &
                   cfgem_coul_adj_recip_ene%mp_dh + &
                   cfgem_coul_adj_recip_ene%dh_mp

   write(out_lun,*)'self_mpole_herm_ene = ', &
                   cfgem_coul_self_ene%mp_ch + &
                   cfgem_coul_self_ene%ch_mp + &
                   cfgem_coul_self_ene%mp_dh + &
                   cfgem_coul_self_ene%dh_mp

   write(out_lun,*)'----------------------------------------'

   write(out_lun,*)'herm-herm components: '

   write(out_lun,*)'cutoff_herm_herm_ene = ', &
                   cfgem_coul_direct_ene%ch_ch + &
                   cfgem_coul_direct_ene%dh_ch + &
                   cfgem_coul_direct_ene%ch_dh + &
                   cfgem_coul_direct_ene%dh_dh

   write(out_lun,*)'recip_herm_herm_ene = ', &
                   cfgem_coul_recip_ene%ch_ch + &
                   cfgem_coul_recip_ene%dh_ch + &
                   cfgem_coul_recip_ene%ch_dh + &
                   cfgem_coul_recip_ene%dh_dh

   write(out_lun,*)'wigner_herm_herm_ene = ', &
                   cfgem_coul_wigner_ene%ch_ch + &
                   cfgem_coul_wigner_ene%dh_ch + &
                   cfgem_coul_wigner_ene%ch_dh + &
                   cfgem_coul_wigner_ene%dh_dh

   write(out_lun,*)'correct_recip_herm_herm_ene = ', &
                   cfgem_coul_adj_recip_ene%ch_ch + &
                   cfgem_coul_adj_recip_ene%dh_ch + &
                   cfgem_coul_adj_recip_ene%ch_dh + &
                   cfgem_coul_adj_recip_ene%dh_dh

   write(out_lun,*)'self_herm_herm_ene = ', &
                   cfgem_coul_self_ene%ch_ch + &
                   cfgem_coul_self_ene%dh_ch + &
                   cfgem_coul_self_ene%ch_dh + &
                   cfgem_coul_self_ene%dh_dh

   call cfgem_exch(aux_file, list_file, user_exch_file, &
                   atom_cnt, atom_crds, free_lun, out_lun, &
                   cfgem_exch_ch_ch_ene, &
                   cfgem_exch_ene_fld, &
                   cfgem_exch_frc)

   write(out_lun,*)'----------------------------------------'
   write(out_lun,*)'tot_exchange_ene = ',cfgem_exch_ch_ch_ene,' in kcals = ', &
             627.509474277d0*cfgem_exch_ch_ch_ene

   write(out_lun,*)'field ene is ', cfgem_exch_ene_fld,  &
      'in kcals ',627.509474277d0*(cfgem_exch_ene_fld)

   ! Deallocate local allocations here:

   if (allocated(cfgem_coul_frc)) deallocate(cfgem_coul_frc)
   if (allocated(cfgem_exch_frc)) deallocate(cfgem_exch_frc)
   if (allocated(atom_crds)) deallocate(atom_crds)

   return

end program cfgem_drvr
