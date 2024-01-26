!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!     University of North Texas
!
!
!     subroutine geminit: initialize GEM, read in data
!
      subroutine geminit

      use cfgem_libmod

      use atmtyp
      use atoms
      use bath
      use bound
      use boxes
      use couple
      use cutoff
      use deriv
      use domdec
      use files
      use gemstuff
      use keys
      use freeze
      use inform
      use iounit
      use langevin
      use mdstuf
      use molcul
      use moldyn
      use mpole
      use neigh
      use polpot
      use units

      implicit none
   
      integer              :: freeunit ! GAC to call tinker function
   
      !atom_crds need to be in MODgemstuff
      !double precision, allocatable, save  :: atom_crds(:,:) ! in angstrom
      double precision, allocatable        :: atm_crd_tmp(:,:) ! in angstrom
   
      double precision     :: x_gem, y_gem, z_gem      ! tmp coord input
   
      double precision     :: tot_mpole_mpole_ene, &
                           tot_mpole_herm_ene, &
                           tot_herm_herm_ene
   
      integer              :: i,j,k
      integer              :: atom_cnt_tmp, tem1, tem2
      integer              :: mol_num      ! trash input..
      integer              :: io_stat
      integer              :: alloc_failed
      integer              :: next ! GAC to read keyfile
      character*5          :: atom_name
      character*20         :: keyword
      character*120        :: record
      character*120        :: string

      integer              :: out_lun

      ! Energy and force outputs of cfgem()

      double precision     :: cfgem_tot_coul_ene ! ALL coulomb
      double precision     :: cfgem_coul_ene_fld ! (for torque)

      !type(cfgem_ene_rec)  :: cfgem_coul_recip_ene
      !type(cfgem_ene_rec)  :: cfgem_coul_direct_ene
      !type(cfgem_ene_rec)  :: cfgem_coul_adj_recip_ene
      !type(cfgem_ene_rec)  :: cfgem_coul_self_ene
      !type(cfgem_ene_rec)  :: cfgem_coul_wigner_ene

      double precision     :: cfgem_exch_ch_ch_ene

      double precision     :: cfgem_exch_ene_fld ! (for torque)


      aux_file = ""
      crd_file = ""
      list_file = ""
      user_coul_file = ""
      user_exch_file = ""

      first_time_gem = 1

      ! GAC get unit number using tinker function
      crd_file_lun = freeunit ()
      free_lun = freeunit ()

      ! make out_lun the stdout unit from tinker in iounit module
      out_lun = iout

      ! put number of atoms in atom_cnt
      atom_cnt = n

      ! put box info in gem box array
      pbc_box_size(1) = xbox
      pbc_box_size(2) = ybox
      pbc_box_size(3) = zbox
      !print *,'box in gem ',pbc_box_size(1),pbc_box_size(2)

      ! Have to read GEM file names from keyword file
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:120)
         if (keyword(1:7) .eq. 'GEMAUX ') then
            call getword (record,aux_file,next)
            !call upcase (aux_file)
         else if (keyword(1:11) .eq. 'GEMAUXEXCH ') then
            call getword (record,aux_exch_file,next)
            !print *,'aux_exch_file = ',trim(aux_exch_file)
            !call upcase (aux_exch_file)
         else if (keyword(1:7) .eq. 'GEMCRD ') then
            call getword (record,crd_file,next)
            !print *,'crd_file = ',trim(crd_file)
            !call upcase (crd_file)
         else if (keyword(1:8) .eq. 'GEMLIST ') then
            call getword (record,list_file,next)
            !call upcase (list_file)
         else if (keyword(1:8) .eq. 'GEMCOUL ') then
            call getword (record,user_coul_file,next)
            !call upcase (user_coul_file)
         else if (keyword(1:8) .eq. 'GEMEXCH ') then
            call getword (record,user_exch_file,next)
            !call upcase (user_exch_file)
         endif
      end do 

      if (aux_exch_file == '') then
         aux_exch_file = aux_file ! If no Exch Aux use Coulomb
      endif  

      if (aux_file == '') then
         write(out_lun,*)'Auxiliary file not specified'
         call fatal
      endif  

      if (crd_file == '') then
         write(out_lun,*)'Coordinate file not specified'
         call fatal
      endif  

      if (list_file == '') then
         write(out_lun,*)'List file not specified'
         call fatal
      endif  

      if (user_exch_file == '') then
         write(out_lun,*)'User exchange file not specified'
         call fatal
      endif  

      if (user_coul_file == '') then
         write(out_lun,*)'User coulomb file not specified'
         call fatal
      endif  

! The coordinate reading stuff, inlined for simplicity, since only used here:

      open(unit=crd_file_lun, file=trim(crd_file), status='old')

      read(crd_file_lun, *, iostat=io_stat) atom_cnt_tmp, tem1, tem2
      if (io_stat .ne. 0) then
         write(out_lun, *)'Bad coordinate file name', crd_file
         call fatal
      end if

      ! consistency check on number of atoms
      if (atom_cnt .ne. atom_cnt_tmp) then
         write(out_lun,*)'Number of atoms in GEM different than TINKER'
         call fatal
      endif
      !print *,'after check crd_file',atom_cnt,atom_cnt_tmp

      if (allocated(atom_crds)) deallocate(atom_crds)
      if (allocated(cfgem_coul_frc)) deallocate(cfgem_coul_frc)
      if (allocated(cfgem_exch_frc)) deallocate(cfgem_exch_frc)
      if (allocated(atm_crd_tmp)) deallocate(atm_crd_tmp)
      allocate(atom_crds(3,atom_cnt), &
            cfgem_coul_frc(3,atom_cnt), &
            cfgem_exch_frc(3,atom_cnt), &
            atm_crd_tmp(3,atom_cnt), &
            stat=alloc_failed)

      if (alloc_failed .ne. 0) then
        write(out_lun, *)'Allocation of arrays failed!'
        call fatal
      end if

      !!allocate temporary array to check atom crd inputs
      !allocate(atm_crd_tmp(3,atom_cnt), stat=alloc_failed)
!
!      if (alloc_failed .ne. 0) then
!        write(out_lun, *)'Allocation of atm_crd_tmp failed!',atom_cnt
!        call fatal
!      end if

      do i = 1, atom_cnt

         read(crd_file_lun, *, iostat=io_stat) atom_name, mol_num, &
             x_gem, y_gem, z_gem

         if (io_stat .ne. 0) then
            write(out_lun,*)'Error reading crds for atom number ', i
            write(out_lun,*)atom_name,mol_num,x_gem,y_gem,z_gem
            call fatal
         end if

         atm_crd_tmp(1,i) = x_gem
         atm_crd_tmp(2,i) = y_gem
         atm_crd_tmp(3,i) = z_gem

         ! here is where we put the TINKER crds into GEM arrays
         atom_crds(1,:) = x!*0.52917721092d0
         atom_crds(2,:) = y!*0.52917721092d0
         atom_crds(3,:) = z!*0.52917721092d0

         ! sanity check (to be removed later) for coordinates
         !if ( (atm_crd_tmp(1,i) .ne. atom_crds(1,i)) .or. &
         !     (atm_crd_tmp(2,i) .ne. atom_crds(2,i)) .or. &
         !     (atm_crd_tmp(3,i) .ne. atom_crds(3,i)) ) then
         !  write(out_lun, *)'Coordinates are different !'
         !  call fatal
         !end if

      end do

      close(crd_file_lun)

      return

      end subroutine geminit

