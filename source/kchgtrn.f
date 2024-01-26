c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine kchgtrn  --  charge transfer term assignment  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "kchgtrn" assigns charge magnitude and damping parameters for
c     charge transfer interactions and processes any new or changed
c     values for these parameters
c
c
      subroutine kchgtrn
      use atoms
      use atmtyp
      use chgpen
      use chgtrn
      use inform
      use iounit
      use kctrn
      use keys
      use mplpot
      use mpole
      use polar
      use polpot
      use potent
      use sizes
      implicit none
      integer i,k
      integer ia,ic,next
      real*8 chtrn,actrn
      logical header
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     process keywords containing charge transfer parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:7) .eq. 'CHGTRN ') then
            k = 0
            chtrn = 0.0d0
            actrn = 0.0d0
            call getnumb (record,k,next)
            string = record(next:240)
            read (string,*,err=10,end=10)  chtrn,actrn
   10       continue
            if (k .gt. 0) then
               if (header .and. .not.silent) then
                  header = .false.
                  write (iout,20)
   20             format (/,' Additional Charge Transfer',
     &                       ' Parameters :',
     &                    //,5x,'Atom Class',13x,'Charge',11x,'Damp',/)
               end if
               if (k .le. maxclass) then
                  ctchg(k) = chtrn
                  ctdmp(k) = actrn
                  if (.not. silent) then
                     write (iout,30)  k,chtrn,actrn
   30                format (6x,i6,7x,f15.4,f15.4)
                  end if
               else
                  write (iout,40)
   40             format (/,' KCHGTRN  --  Too many Charge',
     &                       ' Transfer Parameters')
                  abort = .true.
               end if
            end if
         end if
      end do
c
c     perform dynamic allocation of some global arrays
c
      call alloc_shared_chgct
c
c      if (allocated(chgct))  deallocate (chgct)
c      if (allocated(dmpct))  deallocate (dmpct)
c      allocate (chgct(n))
c      allocate (dmpct(n))
c
c     assign the charge transfer charge and alpha parameters 
c     
      do i = 1, n
         ic = class(i)
         chgct(i) = ctchg(ic)
         dmpct(i) = ctdmp(ic)
      end do
c
c     process keywords containing atom specific charge transfer
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:7) .eq. 'CHGTRN ') then
            ia = 0
            chtrn = 0.0d0
            actrn = 0.0d0
            string = record(next:240)
            read (string,*,err=70,end=70)  ia,chtrn,actrn
            if (ia.lt.0 .and. ia.ge.-n) then
               ia = -ia
               if (header .and. .not.silent) then
                  header = .false.
                  write (iout,50)
   50             format (/,' Additional Charge Transfer Values',
     &                       ' for Specific Atoms :',
     &                    //,8x,'Atom',16x,'Charge',11x,'Damp',/)
               end if
               if (.not. silent) then
                  write (iout,60)  ia,chtrn,actrn
   60             format (6x,i6,7x,f15.4,f15.4)
               end if
               chgct(ia) = chtrn
               dmpct(ia) = actrn
            end if
   70       continue
         end if
      end do
c
c     remove zero or undefined electrostatic sites from the list
c
      npole = 0
      ncp = 0
      npolar = 0
      nct = 0
      do i = 1, n
c         if (polarity(i) .eq. 0.0d0)  douind(i) = .false.
         if (polsiz(i).ne.0 .or. polarity(i).ne.0.0d0 .or.
     &          chgct(i).ne. 0.0d0 .or. dmpct(i).ne.0.0d0) then
            npole = npole + 1
            ipole(npole) = i
            pollist(i) = npole
            zaxis(npole) = zaxis(i)
            xaxis(npole) = xaxis(i)
            yaxis(npole) = yaxis(i)
            polaxe(npole) = polaxe(i)
            do k = 1, maxpole
               pole(k,npole) = pole(k,i)
            end do
            if (palpha(i) .ne. 0.0d0)  ncp = ncp + 1
            pcore(npole) = pcore(i)
            pval(npole) = pval(i)
            palpha(npole) = palpha(i)
            if (polarity(i) .ne. 0.0d0) then
               npolar = npolar + 1
c               ipolar(npolar) = npole
c               douind(i) = .true.
            end if
            polarity(npole) = polarity(i)
            thole(npole) = thole(i)
c            dirdamp(npole) = dirdamp(i)
            if (chgct(i).ne.0.0d0 .or. dmpct(i).ne.0.0d0) then
               nct = nct + 1
            end if
            chgct(npole) = chgct(i)
            dmpct(npole) = dmpct(i)
         end if
      end do
c
c     test multipoles at chiral sites and invert if necessary
c
      call chkpole(.true.)
c
c     turn off individual electrostatic potentials if not used
c
      if (npole .eq. 0)  use_mpole = .false.
c      if (ncp .ne. 0)  use_chgpen = .true.
      if (npolar .eq. 0)  use_polar = .false.
      if (use_polar) then
         do i = 1, npole
            if (thole(i) .ne. 0.0d0) then
               use_thole = .true.
               goto 80
            end if
         end do
   80    continue
         do i = 1, npole
            if (dirdamp(i) .ne. 0.0d0) then
               use_dirdamp = .true.
               goto 90
            end if
         end do
   90    continue
      end if
      if (nct .eq. 0)  use_chgtrn = .false.
      return
      end
c
c     subroutine alloc_shared_chgct : allocate shared memory pointers for charge transfer
c     parameter arrays
c
      subroutine alloc_shared_chgct
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
      use sizes
      use atoms
      use domdec
      use chgtrn
      use mpi
      implicit none
      integer :: win,win2
      INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
      INTEGER :: disp_unit,ierr,total
      TYPE(C_PTR) :: baseptr
      integer :: arrayshape(1),arrayshape2(2)
c
      if (associated(chgct)) then
        CALL MPI_Win_shared_query(winchgct, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winchgct,ierr)
      end if
      if (associated(dmpct)) then
        CALL MPI_Win_shared_query(windmpct, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(windmpct,ierr)
      end if
c
c     chgct
c
      arrayshape=(/n/)
      if (hostrank == 0) then
        windowsize = int(n,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winchgct, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winchgct, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,chgct,arrayshape)
c
c     dmpct
c
      arrayshape=(/n/)
      if (hostrank == 0) then
        windowsize = int(n,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, windmpct, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(windmpct, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,dmpct,arrayshape)
c
      return
      end
