c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine angles  --  locate and store bond angles  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "angles" finds the total number of bond angles and stores
c     the atom numbers of the atoms defining each angle; for
c     each angle to a trivalent central atom, the third bonded
c     atom is stored for use in out-of-plane bending
c
c
      subroutine angles(init)
      use angle
      use atmlst
      use atoms
      use couple
      use domdec
      use iounit
      implicit none
      integer i,j,k,m,iglob,jglob,kglob
      logical init,docompute
      real*8 pos(3,4)
c
      if (init) then
c
c     loop over all atoms, storing the atoms in each bond angle
c
        nangle = 0
        do i = 1, n
           do j = 1, n12(i)-1
              do k = j+1, n12(i)
                nangle = nangle + 1
                if (nangle .gt. 4*n) then
                   if (rank.eq.0) write (iout,10)
   10              format (/,' ANGLES  --  Too many Bond Angles;',
     &                        ' Increase MAXANG')
                   call fatal
                end if
              end do
           end do
        end do
c
c       allocate arrays
c
        call alloc_shared_angles
c
        nangleloc = 0
        do i = 1, n
           m = 0
           do j = 1, n12(i)-1
              do k = j+1, n12(i)
                nangleloc = nangleloc + 1
                m = m + 1
                anglist(m,i) = nangleloc
                iang(1,nangleloc) = i12(j,i)
                iang(2,nangleloc) = i
                iang(3,nangleloc) = i12(k,i)
                iang(4,nangleloc) = 0
              end do
           end do
c
c       set the out-of-plane atom for angles at trivalent centers
c
           if (n12(i) .eq. 3) then
              iang(4,nangleloc) = i12(1,i)
              iang(4,nangleloc-1) = i12(2,i)
              iang(4,nangleloc-2) = i12(3,i)
           end if
        end do
        if (allocated(angleloc)) deallocate(angleloc)
        allocate (angleloc(nangle))
      end if
      if (allocated(angleglob)) deallocate(angleglob)
      allocate (angleglob(4*nbloc))
      nangleloc = 0
c      do i = 1, nbloc
c        m = 0
c        iglob = glob(i)
c        pos(1,1) = x(iglob)
c        pos(2,1) = y(iglob)
c        pos(3,1) = z(iglob)
c        do j = 1, n12(iglob)-1
c          jglob = i12(j,iglob)
c          pos(1,2) = x(jglob)
c          pos(2,2) = y(jglob)
c          pos(3,2) = z(jglob)
c          do k = j+1, n12(iglob)
c            m = m + 1
c            kglob = i12(k,iglob)
c            pos(1,3) = x(kglob)
c            pos(2,3) = y(kglob)
c            pos(3,3) = z(kglob)
c            call midpointgroup(pos,3,docompute)
c            if (.not.(docompute)) cycle
c            nangleloc = nangleloc + 1
c            angleglob(nangleloc) = anglist(m,iglob)
c          end do
c        end do
c      end do
      do i = 1, nloc
        m = 0
        iglob = glob(i)
        do j = 1, n12(iglob)-1
          do k = j+1, n12(iglob)
            nangleloc = nangleloc + 1
            m = m + 1
            angleglob(nangleloc) = anglist(m,iglob)
            angleloc(anglist(m,iglob)) = nangleloc
          end do
        end do
      end do
      return
      end
c
c
c     subroutine alloc_shared_angles : allocate shared memory pointers for angles
c     parameter arrays
c
      subroutine alloc_shared_angles
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
      use angang
      use angle
      use angpot
      use atmlst
      use atoms
      use bitor
      use domdec
      use improp
      use imptor
      use opbend
      use opdist
      use strbnd
      use urey
      use mpi
      implicit none
      integer :: win
      INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
      INTEGER :: disp_unit,ierr
      TYPE(C_PTR) :: baseptr
      integer :: arrayshape(1),arrayshape2(2)
c
c      if (associated(anglist)) deallocate(anglist)
c      if (associated(iang)) deallocate (iang)
c      if (associated(ak)) deallocate (ak)
c      if (associated(anat)) deallocate (anat)
c      if (associated(afld)) deallocate (afld)
c      if (associated(angtyp)) deallocate (angtyp)
c      if (associated(isb)) deallocate (isb)
c      if (associated(sbk)) deallocate (sbk)
c      if (associated(uk)) deallocate (uk)
c      if (associated(ul)) deallocate (ul)
c      if (associated(iury)) deallocate (iury)
c      if (associated(nbbitors)) deallocate (nbbitors)
c      if (associated(nbstrbnd)) deallocate (nbstrbnd)
c      if (associated(nburey)) deallocate (nburey)
c      if (associated(nbangang)) deallocate (nbangang)
c      if (associated(nbopbend)) deallocate (nbopbend)
c      if (associated(nbopdist)) deallocate (nbopdist)
c      if (associated(nbimprop)) deallocate (nbimprop)
c      if (associated(nbimptor)) deallocate (nbimptor)

      if (associated(anglist)) then
        CALL MPI_Win_shared_query(winanglist, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winanglist,ierr)
      end if
      if (associated(iang)) then
        CALL MPI_Win_shared_query(winiang, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winiang,ierr)
      end if
      if (associated(ak)) then
        CALL MPI_Win_shared_query(winak, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winak,ierr)
      end if
      if (associated(anat)) then
        CALL MPI_Win_shared_query(winanat, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winanat,ierr)
      end if
      if (associated(afld)) then
        CALL MPI_Win_shared_query(winafld, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winafld,ierr)
      end if
      if (associated(angtyp)) then
        CALL MPI_Win_shared_query(winangtyp, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winangtyp,ierr)
      end if
      if (associated(isb)) then
        CALL MPI_Win_shared_query(winisb, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winisb,ierr)
      end if
      if (associated(sbk)) then
        CALL MPI_Win_shared_query(winsbk, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winsbk,ierr)
      end if
      if (associated(uk)) then
        CALL MPI_Win_shared_query(winuk, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winuk,ierr)
      end if
      if (associated(ul)) then
        CALL MPI_Win_shared_query(winul, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winul,ierr)
      end if
      if (associated(iury)) then
        CALL MPI_Win_shared_query(winiury, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winiury,ierr)
      end if
      if (associated(nbbitors)) then
        CALL MPI_Win_shared_query(winnbbitors, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winnbbitors,ierr)
      end if
      if (associated(nbstrbnd)) then
        CALL MPI_Win_shared_query(winnbstrbnd, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winnbstrbnd,ierr)
      end if
      if (associated(nburey)) then
        CALL MPI_Win_shared_query(winnburey, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winnburey,ierr)
      end if
      if (associated(nbangang)) then
        CALL MPI_Win_shared_query(winnbangang, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winnbangang,ierr)
      end if
      if (associated(nbopbend)) then
        CALL MPI_Win_shared_query(winnbopbend, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winnbopbend,ierr)
      end if
      if (associated(nbopdist)) then
        CALL MPI_Win_shared_query(winnbopdist, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winnbopdist,ierr)
      end if
      if (associated(nbimprop)) then
        CALL MPI_Win_shared_query(winnbimprop, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winnbimprop,ierr)
      end if
      if (associated(nbimptor)) then
        CALL MPI_Win_shared_query(winnbimptor, 0, windowsize, disp_unit,
     $  baseptr, ierr)
        CALL MPI_Win_free(winnbimptor,ierr)
      end if
c
c     anglist
c
      arrayshape2=(/16,n/)
      if (hostrank == 0) then
        windowsize = int(16*n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winanglist, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winanglist, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,anglist,arrayshape2)
c
c     iang
c
      arrayshape2=(/4,nangle/)
      if (hostrank == 0) then
        windowsize = int(4*nangle,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winiang, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winiang, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,iang,arrayshape2)
c
c     ak
c
      arrayshape=(/nangle/)
      if (hostrank == 0) then
        windowsize = int(nangle,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winak, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winak, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,ak,arrayshape)
c
c     anat
c
      arrayshape=(/nangle/)
      if (hostrank == 0) then
        windowsize = int(nangle,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winanat, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winanat, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,anat,arrayshape)
c
c     afld
c
      arrayshape=(/nangle/)
      if (hostrank == 0) then
        windowsize = int(nangle,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winafld, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winafld, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,afld,arrayshape)
c
c     angtyp
c
      arrayshape=(/nangle/)
      if (hostrank == 0) then
        windowsize = int(nangle,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winangtyp, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winangtyp, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,angtyp,arrayshape)
c
c     isb
c
      arrayshape2=(/3,nangle/)
      if (hostrank == 0) then
        windowsize = int(3*nangle,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winisb, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winisb, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,isb,arrayshape2)
c
c     sbk
c
      arrayshape2=(/2,nangle/)
      if (hostrank == 0) then
        windowsize = int(2*nangle,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winsbk, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winsbk, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,sbk,arrayshape2)
c
c    uk
c
      arrayshape=(/nangle/)
      if (hostrank == 0) then
        windowsize = int(nangle,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winuk, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winuk, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,uk,arrayshape)
c
c    ul
c
      arrayshape=(/nangle/)
      if (hostrank == 0) then
        windowsize = int(nangle,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winul, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winul, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,ul,arrayshape)
c
c    iury
c
      arrayshape2=(/3,nangle/)
      if (hostrank == 0) then
        windowsize = int(3*nangle,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winiury, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winiury, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,iury,arrayshape2)
c
c     nbbitors
c
      arrayshape=(/nangle/)
      if (hostrank == 0) then
        windowsize = int(nangle,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winnbbitors, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winnbbitors, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,nbbitors,arrayshape)
c
c     nbstrbnd
c
      arrayshape=(/nangle/)
      if (hostrank == 0) then
        windowsize = int(nangle,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winnbstrbnd, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winnbstrbnd, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,nbstrbnd,arrayshape)
c
c     nburey
c
      arrayshape=(/nangle/)
      if (hostrank == 0) then
        windowsize = int(nangle,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winnburey, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winnburey, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,nburey,arrayshape)
c
c     nbangang
c
      arrayshape=(/n/)
      if (hostrank == 0) then
        windowsize = int(n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winnbangang, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winnbangang, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,nbangang,arrayshape)
c
c     nbopbend
c
      arrayshape=(/nangle/)
      if (hostrank == 0) then
        windowsize = int(nangle,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winnbopbend, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winnbopbend, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,nbopbend,arrayshape)
c
c     nbopdist
c
      arrayshape=(/n/)
      if (hostrank == 0) then
        windowsize = int(n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winnbopdist, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winnbopdist, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,nbopdist,arrayshape)
c
c     nbimprop
c
      arrayshape=(/n/)
      if (hostrank == 0) then
        windowsize = int(n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winnbimprop, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winnbimprop, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,nbimprop,arrayshape)
c
c     nbimptor
c
      arrayshape=(/n/)
      if (hostrank == 0) then
        windowsize = int(n,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      else
        windowsize = 0_MPI_ADDRESS_KIND
      end if
      disp_unit = 1
c
c    allocation
c
      CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL,
     $  hostcomm, baseptr, winnbimptor, ierr)
      if (hostrank /= 0) then
        CALL MPI_Win_shared_query(winnbimptor, 0, windowsize, disp_unit,
     $  baseptr, ierr)
      end if
c
c    association with fortran pointer
c
      CALL C_F_POINTER(baseptr,nbimptor,arrayshape)
      return
      end
