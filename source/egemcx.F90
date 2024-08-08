!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!     University of North Texas
!
!
!     subroutine egemcx: calculate GEM Coulomb and Exchange energy/forces, 
!
      subroutine egemcx
      use cfgem_libmod
      use gemstuff
      use atmtyp
      use atoms
      use bath
      use bound
      use couple
      use cutoff
      use deriv
      use domdec
      use energi
      use files
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
   

      !integer              :: iargc
      !integer              :: argc, arg
      integer              :: freeunit ! GAC to call tinker function
      !character(len=256)   :: argv
      !character(len=80)    :: aux_file, aux_exch_file, crd_file
      !character(len=80)    :: user_coul_file, user_exch_file
      !character(len=80)    :: list_file
   
      !atom_crds need to be in MOD (si a huevo)
      !double precision, allocatable, save  :: atom_crds(:,:) ! in angstrom
      double precision, allocatable        :: atm_crd_tmp(:,:) ! in angstrom
   
      double precision     :: x_gem, y_gem, z_gem      ! tmp coord input
      double precision     :: au_to_kcal, au_to_kcal_ang
   
      double precision     :: tot_mpole_mpole_ene, tot_mpole_herm_ene
      double precision     :: tot_herm_herm_ene
   
      integer              :: i,j,k,h
      !integer              :: atom_cnt
      integer              :: mol_num      ! trash input..
      integer              :: io_stat
      integer              :: alloc_failed
      integer              :: next ! GAC to read keyfile
      character*5          :: atom_name
      character*20         :: keyword


      ! If necessary, move the lun to something unused.  Should not present
      ! problems here..

      !integer, parameter   :: crd_file_lun = 10
      !integer, parameter   :: free_lun = 11
      !integer, parameter   :: out_lun = 6

      !integer              :: crd_file_lun 
      !integer              :: free_lun 
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


      ! GAC the TINKER crds in GEM arrays were put in BEEMAN
         atom_cnt = n

      ! zero out coulomb energy and derivs
      em = 0.d0
      dem = 0d0
      demrec = 0d0
      cfgem_tot_coul_ene = 0.d0
      cfgem_coul_ene_fld = 0.d0
      !cfgem_coul_recip_ene = 0d0
      !cfgem_coul_direct_ene = 0d0
      !cfgem_coul_adj_recip_ene = 0d0
      !cfgem_coul_self_ene = 0d0
      !cfgem_coul_wigner_ene = 0d0
      cfgem_coul_frc = 0d0
      cfgem_exch_ch_ch_ene = 0d0
      cfgem_exch_ene_fld = 0d0
      cfgem_exch_frc = 0d0
      out_lun = iout

      ! constants to convert E and force
      au_to_kcal = 627.509474277d0
      au_to_kcal_ang = au_to_kcal/0.52917721092d0

         !allocate(cfgem_coul_frc(3,n),cfgem_exch_frc(3,n),status=alloc_fail)

      ! Atom coordinate input is in angstrom...

      !do h = 1, atom_cnt
      !   print *,atom_crds(1,h),atom_crds(2,h),atom_crds(3,h)
      !enddo
      ! print *,aux_file, list_file, user_coul_file
      ! print *,atom_cnt, free_lun, out_lun

      call cfgem_coul(aux_file, list_file, user_coul_file, &
                   atom_cnt, atom_crds, free_lun, out_lun, &
                   cfgem_tot_coul_ene, &
                   cfgem_coul_ene_fld, &
                   cfgem_coul_frc,&
                   pbc_box_size)
      

      ! don't need all the _ene structs for detailed info
      !call cfgem_coul(aux_file, list_file, user_coul_file, &
      !             atom_cnt, atom_crds, free_lun, out_lun, &
      !             cfgem_tot_coul_ene, &
      !             cfgem_coul_ene_fld, &
      !             !cfgem_coul_recip_ene, &
      !             !cfgem_coul_direct_ene, &
      !             !cfgem_coul_adj_recip_ene, &
      !             !cfgem_coul_self_ene, &
      !             !cfgem_coul_wigner_ene, &
      !             cfgem_coul_frc) ! ME FALTA EL PINCHE PBC_BOX_SIZE!!
!
      ! update TINKER-HP Energy and Force Arrays with GEM Coulomb
      !S:JORGE
      if (.not. onlyHal) then
         em = em + (cfgem_tot_coul_ene*au_to_kcal)
      end if
      !E:JORGE

      !dem(:,:) = dem(:,:) + (cfgem_coul_frc(:,:)*au_to_kcal_ang)
      !dem(:,:) = dem(:,:) + (cfgem_coul_frc(:,:)*au_to_kcal_ang*-1.d0)
      do k=1,3
       do h = 1,atom_cnt
         dem(k,h)=dem(k,h)+(cfgem_coul_frc(k,h)*au_to_kcal_ang*(-1.d0))
       enddo 
      enddo 

      !tot_mpole_mpole_ene = cfgem_coul_direct_ene%mp_mp + &
      !                   cfgem_coul_recip_ene%mp_mp + &
      !                   cfgem_coul_wigner_ene%mp_mp + &
      !                   cfgem_coul_adj_recip_ene%mp_mp + &
      !                   cfgem_coul_self_ene%mp_mp

      !tot_mpole_herm_ene  = cfgem_coul_direct_ene%mp_ch + &
      !                   cfgem_coul_direct_ene%ch_mp + &
      !                   cfgem_coul_recip_ene%mp_ch + &
      !                   cfgem_coul_recip_ene%mp_dh + &
      !                   cfgem_coul_recip_ene%ch_mp + &
      !                   cfgem_coul_recip_ene%dh_mp + &
      !                   cfgem_coul_wigner_ene%mp_ch + &
      !                   cfgem_coul_wigner_ene%mp_dh + &
      !                   cfgem_coul_wigner_ene%ch_mp + &
      !                   cfgem_coul_wigner_ene%dh_mp + &
      !                   cfgem_coul_adj_recip_ene%mp_ch + &
      !                   cfgem_coul_adj_recip_ene%mp_dh + &
      !                   cfgem_coul_adj_recip_ene%ch_mp + &
      !                   cfgem_coul_adj_recip_ene%dh_mp + &
      !                   cfgem_coul_self_ene%mp_ch + &
      !                   cfgem_coul_self_ene%mp_dh + &
      !                   cfgem_coul_self_ene%ch_mp + &
      !                   cfgem_coul_self_ene%dh_mp

      !tot_herm_herm_ene =   cfgem_coul_direct_ene%ch_ch + &
      !                   cfgem_coul_recip_ene%ch_ch + &
      !                   cfgem_coul_recip_ene%dh_ch + &
      !                   cfgem_coul_recip_ene%ch_dh + &
      !                   cfgem_coul_recip_ene%dh_dh + &
      !                   cfgem_coul_wigner_ene%ch_ch + &
      !                   cfgem_coul_wigner_ene%ch_dh + &
      !                   cfgem_coul_wigner_ene%dh_ch + &
      !                   cfgem_coul_wigner_ene%dh_dh + &
      !                   cfgem_coul_adj_recip_ene%ch_ch + &
      !                   cfgem_coul_adj_recip_ene%ch_dh + &
      !                   cfgem_coul_adj_recip_ene%dh_ch + &
      !                   cfgem_coul_adj_recip_ene%dh_dh + &
      !                   cfgem_coul_self_ene%ch_ch + &
      !                   cfgem_coul_self_ene%ch_dh + &
      !                   cfgem_coul_self_ene%dh_ch + &
      !                   cfgem_coul_self_ene%dh_dh

      !write(iout,*)'--------------------------------------------------'
      !write(out_lun,*)'tot_mpole_mpole_ene = ',tot_mpole_mpole_ene
      !write(out_lun,*)'tot_mpole_herm_ene = ',tot_mpole_herm_ene
      !write(out_lun,*)'tot_herm_herm_ene = ',tot_herm_herm_ene
      !write(iout,*)'tot_coulomb_ene in kcals = ', &
      !          627.509474277d0*cfgem_tot_coul_ene
      !write(iout,*)
      !write(iout,*)'field ene in kcals ',627.509474277d0*(cfgem_coul_ene_fld)

      call cfgem_exch(aux_exch_file, list_file, user_exch_file, &
                      atom_cnt, atom_crds, free_lun, out_lun, &
                      cfgem_exch_ch_ch_ene, &
                      cfgem_exch_ene_fld, &
                      cfgem_exch_frc, &
                      pbc_box_size)

      ! update TINKER-HP Energy and Force Arrays with GEM Exchange
      ! NOTE WE'RE ADDING EXCHANGE-REPULSION IN COULOMB FOR NOW!!
      !S:JORGE
      if (.not. onlyHal) then
         em = em + (cfgem_exch_ch_ch_ene*au_to_kcal)
      end if
      !E:JORGE

      !dem(:,:) = dem(:,:) + (cfgem_exch_frc(:,:)*au_to_kcal_ang)
      !dem(:,:) = dem(:,:) + (cfgem_exch_frc(:,:)*au_to_kcal_ang*-1.d0)
      do k=1,3
       do h = 1,atom_cnt
         dem(k,h)=dem(k,h)+(cfgem_exch_frc(k,h)*au_to_kcal_ang*(-1.d0))
       enddo 
      enddo 


      !write(iout,*)'tot_exchange_ene  in kcals = ', &
      !                 627.509474277d0*cfgem_exch_ch_ch_ene
      !write(iout,*)'field ene is ', cfgem_exch_ene_fld, &
      !   'in kcals ',627.509474277d0*(cfgem_exch_ene_fld)
      !do h = 1,atom_cnt
      !   write(iout,*)'der ',(dem(h,k),k=1,3)
      !enddo 
      !write(iout,*)'--------------------------------------------------'

      !S:JORGE
      ene_gem_coul = 627.509474277d0*cfgem_tot_coul_ene
      ene_gem_exch = 627.509474277d0*cfgem_exch_ch_ch_ene
      !E:JORGE

      if (first_time_gem == 1) first_time_gem = 0 ! to avoid reading every time
      if (allocated(cfgem_coul_frc)) deallocate(cfgem_coul_frc)
      if (allocated(cfgem_exch_frc)) deallocate(cfgem_exch_frc)
      if (allocated(atom_crds)) deallocate(atom_crds)

      return

      end subroutine egemcx

!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!     University of North Texas
!
!
!     subroutine egemcx: calculate GEM Coulomb and Exchange energy/forces, 
!
      subroutine egemcx0
      use cfgem_libmod
      use gemstuff
      use atmtyp
      use atoms
      use bath
      use bound
      use couple
      use cutoff
      use domdec
      use energi
      use files
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

      double precision, allocatable        :: atm_crd_tmp(:,:) ! in angstrom
   
      double precision     :: x_gem, y_gem, z_gem      ! tmp coord input
      double precision     :: au_to_kcal, au_to_kcal_ang
   
      double precision     :: tot_mpole_mpole_ene, tot_mpole_herm_ene
      double precision     :: tot_herm_herm_ene
   
      integer              :: i,j,k,h
      !integer              :: atom_cnt
      integer              :: mol_num      ! trash input..
      integer              :: io_stat
      integer              :: alloc_failed
      integer              :: next ! GAC to read keyfile
      character*5          :: atom_name
      character*20         :: keyword


      integer              :: out_lun

      ! Energy and force outputs of cfgem()

      double precision     :: cfgem_tot_coul_ene ! ALL coulomb

      double precision     :: cfgem_coul_ene_fld ! (for torque)
      double precision,allocatable  :: tmp_coul_frc(:,:) ! (for torque)

      double precision     :: cfgem_exch_ch_ch_ene

      double precision     :: cfgem_exch_ene_fld ! (for torque)
      double precision,allocatable  :: tmp_exch_frc(:,:) ! (for torque)


      ! GAC the TINKER crds in GEM arrays were put in BEEMAN
         atom_cnt = n


      allocate(tmp_coul_frc(3,atom_cnt),tmp_exch_frc(3,atom_cnt)) 
      ! zero out coulomb energy and derivs
      em = 0.d0
      cfgem_tot_coul_ene = 0.d0
      cfgem_coul_ene_fld = 0.d0
      cfgem_exch_ch_ch_ene = 0d0
      cfgem_exch_ene_fld = 0d0
      tmp_coul_frc = 0d0
      tmp_exch_frc = 0d0
      out_lun = iout

      ! constants to convert E and force
      au_to_kcal = 627.509474277d0
      au_to_kcal_ang = au_to_kcal/0.52917721092d0

      call cfgem_coul(aux_file, list_file, user_coul_file, &
                   atom_cnt, atom_crds, free_lun, out_lun, &
                   cfgem_tot_coul_ene, &
                   cfgem_coul_ene_fld, &
                   tmp_coul_frc,&
                   pbc_box_size)
      em = em + (cfgem_tot_coul_ene*au_to_kcal)

      call cfgem_exch(aux_exch_file, list_file, user_exch_file, &
                      atom_cnt, atom_crds, free_lun, out_lun, &
                      cfgem_exch_ch_ch_ene, &
                      cfgem_exch_ene_fld, &
                      tmp_exch_frc,&
                      pbc_box_size)

      if(allocated(tmp_coul_frc)) deallocate(tmp_coul_frc)
      if(allocated(tmp_exch_frc)) deallocate(tmp_exch_frc)
      if (allocated(atom_crds)) deallocate(atom_crds)

      em = em + (cfgem_exch_ch_ch_ene*au_to_kcal)

      return

      end subroutine egemcx0

