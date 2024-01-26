c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine kimptor  --  improper torsion parameters  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "kimptor" assigns torsional parameters to each improper
c     torsion in the structure and processes any changed values
c
c
      subroutine kimptor(init)
      use atmlst
      use atmtyp
      use atoms
      use couple
      use domdec
      use keys
      use imptor
      use inform
      use iounit
      use kitors
      use math
      use potent
      use tors
      use utils
      implicit none
      integer i,j,k,nti
      integer ia,ib,ic,id
      integer ita,itb,itc,itd
      integer ita1,itb1,itc1,itd1
      integer iglob,imptorcount,nitorsloc1
      integer size,next,zero
      integer ft(3)
      integer::isys=0
      integer*8 ipt
      real*8 angle,symm
      real*8 vt(3),st(3)
      logical header,done
      character*4 pa,pb,pc,pd
      character*4 zeros
      integer*8 pt(6),pti,pt0,pt1,pt2,pt3
c      character*16 blank,pti
c      character*16 pt0,pt1
c      character*16 pt2,pt3
c      character*16 pt(6)
      character*20 keyword
      character*120 record
      character*120 string
      logical init
c
c      blank = '                '
c      zeros = '0000'
      if (init) then
c
c     process keywords containing improper torsion parameters
c
        header = .true.
        do i = 1, nkey
           next = 1
           record = keyline(i)
           call gettext (record,keyword,next)
           call upcase (keyword)
           if (keyword(1:8) .eq. 'IMPTORS ') then
              ia = 0
              ib = 0
              ic = 0
              id = 0
              do j = 1, 3
                 vt(j) = 0.0d0
                 st(j) = 0.0d0
                 ft(j) = 0
              end do
              string = record(next:120)
              read (string,*,err=10,end=10)  ia,ib,ic,id,
     &                                       (vt(j),st(j),ft(j),j=1,3)
   10         continue
c              size = 4
c              call numeral (ia,pa,size)
c              call numeral (ib,pb,size)
c              call numeral (ic,pc,size)
c              call numeral (id,pd,size)
c              pti = pa//pb//pc//pd
              call front_convert_base(ia,ib,ic,id,0,ipt)
              call torphase (ft,vt,st)
              if (.not. silent) then
                 if (header) then
                    header = .false.
                    if (rank.eq.0) write (iout,20)
   20             format (/,' Additional Improper Torsion Parameters :',
     &                    //,5x,'Atom Classes',15x,'1-Fold',12x,
     &                       '2-Fold',12x,'3-Fold',/)
                 end if
                 if (rank.eq.0) write (iout,30)  ia,ib,ic,id,
     &               (vt(j),st(j),j=1,3)
   30            format (4x,4i4,2x,3(f11.3,f7.1))
              end if
              do j = 1, maxnti
                 if (kti(j).eq.-1 .or. kti(j).eq.ipt) then
                    kti(j) = pti
                    ti1(1,j) = vt(1)
                    ti1(2,j) = st(1)
                    ti2(1,j) = vt(2)
                    ti2(2,j) = st(2)
                    ti3(1,j) = vt(3)
                    ti3(2,j) = st(3)
                    goto 50
                 end if
              end do
              if (rank.eq.0) write (iout,40)
   40         format (/,' KIMPTOR  --  Too many Improper Torsion',
     &                   ' Parameters')
              abort = .true.
   50         continue
           end if
        end do
c
c       determine the total number of forcefield parameters
c
        nti = maxnti
        do i = maxnti, 1, -1
           if (kti(i) .eq. -1)  nti = i - 1
        end do
c
c       assign improper torsional parameters for each improper torsion;
c       multiple symmetrical parameters are given partial weights
c
        nitors = 0
        if (nti .ne. 0) then
           do i = 1, n
              if (n12(i) .eq. 3) then
                 ia = i12(1,i)
                 ib = i12(2,i)
                 ic = i
                 id = i12(3,i)
                 nbimptor(i) = nitors
                 ita = class(ia)
                 itb = class(ib)
                 itc = class(ic)
                 itd = class(id)
c                 size = 4
c                 call numeral (ita,pa,size)
c                 call numeral (itb,pb,size)
c                 call numeral (itc,pc,size)
c                 call numeral (itd,pd,size)
c                 pt(1) = pa//pb//pc//pd
                 call front_convert_base(ita,itb,itc,itd,0,pt(1))
c                 pt(2) = pb//pa//pc//pd
                 call front_convert_base(itb,ita,itc,itd,0,pt(2))
c                 pt(3) = pa//pd//pc//pb
                 call front_convert_base(ita,itd,itc,itb,0,pt(3))
c                 pt(4) = pd//pa//pc//pb
                 call front_convert_base(itd,ita,itc,itb,0,pt(4))
c                 pt(5) = pb//pd//pc//pa
                 call front_convert_base(itb,itd,itc,ita,0,pt(5))
c                 pt(6) = pd//pb//pc//pa
                 call front_convert_base(itd,itb,itc,ita,0,pt(6))
c                 pt3 = zeros//zeros//pc//pd
                 call front_convert_base(0,0,itc,itd,0,pt3)
c                 pt2 = zeros//zeros//pc//pb
                 call front_convert_base(0,0,itc,itb,0,pt2)
c                 pt1 = zeros//zeros//pc//pa
                 call front_convert_base(0,0,itc,ita,0,pt1)
c                 pt0 = zeros//zeros//pc//zeros
                 call front_convert_base(0,0,itc,0,0,pt0)
                 symm = 1.0d0
c                 if (pa.eq.pb .or. pa.eq.pd .or. pb.eq.pd)  symm = 2.0d0
c                 if (pa.eq.pb .and. pa.eq.pd .and. pb.eq.pd)  symm = 6.0d0
                 if (ita.eq.itb .or. ita.eq.itd .or. itb.eq.itd)  
     $               symm = 2.0d0
                 if (ita.eq.itb .and. ita.eq.itd .and. itb.eq.itd)  
     $               symm = 6.0d0
                 done = .false.
                 do j = 1, nti
                   call back_convert_base(ita1,itb1,itc1,itd1,zero,
     $                 kti(j))
c                    if (kti(j)(9:12) .eq. pc) then
                    if (itc1 .eq. itc) then
                       do k = 1, 6
                          if (kti(j) .eq. pt(k)) then
                             nitors = nitors + 1
                             iitors(3,nitors) = ic
                             if (k .eq. 1) then
                                iitors(1,nitors) = ia
                                iitors(2,nitors) = ib
                                iitors(4,nitors) = id
                             else if (k .eq. 2) then
                                iitors(1,nitors) = ib
                                iitors(2,nitors) = ia
                                iitors(4,nitors) = id
                             else if (k .eq. 3) then
                                iitors(1,nitors) = ia
                                iitors(2,nitors) = id
                                iitors(4,nitors) = ib
                             else if (k .eq. 4) then
                                iitors(1,nitors) = id
                                iitors(2,nitors) = ia
                                iitors(4,nitors) = ib
                             else if (k .eq. 5) then
                                iitors(1,nitors) = ib
                                iitors(2,nitors) = id
                                iitors(4,nitors) = ia
                             else if (k .eq. 6) then
                                iitors(1,nitors) = id
                                iitors(2,nitors) = ib
                                iitors(4,nitors) = ia
                             end if
                             itors1(1,nitors) = ti1(1,j) / symm
                             itors1(2,nitors) = ti1(2,j)
                             itors2(1,nitors) = ti2(1,j) / symm
                             itors2(2,nitors) = ti2(2,j)
                             itors3(1,nitors) = ti3(1,j) / symm
                             itors3(2,nitors) = ti3(2,j)
                             done = .true.
                             if (.not.is_find8(kti_sys(1),isys,kti(j)))
     $                          then
                                isys = isys + 1
                                kti_sys(isys) = pt(k)
                             end if
                          end if
                       end do
                    end if
                 end do
                 if (.not. done) then
                    do j = 1, nti
                       if (kti(j) .eq. pt1) then
                          symm = 3.0d0
                          do k = 1, 3
                             nitors = nitors + 1
                             iitors(3,nitors) = ic
                             if (k .eq. 1) then
                                iitors(1,nitors) = ia
                                iitors(2,nitors) = ib
                                iitors(4,nitors) = id
                             else if (k .eq. 2) then
                                iitors(1,nitors) = ib
                                iitors(2,nitors) = id
                                iitors(4,nitors) = ia
                             else if (k .eq. 3) then
                                iitors(1,nitors) = id
                                iitors(2,nitors) = ia
                                iitors(4,nitors) = ib
                             end if
                             itors1(1,nitors) = ti1(1,j) / symm
                             itors1(2,nitors) = ti1(2,j)
                             itors2(1,nitors) = ti2(1,j) / symm
                             itors2(2,nitors) = ti2(2,j)
                             itors3(1,nitors) = ti3(1,j) / symm
                             itors3(2,nitors) = ti3(2,j)
                          end do
                          done = .true.
                          if (.not.is_find8(kti_sys(1),isys,kti(j)))
     $                       then
                             isys = isys + 1
                             kti_sys(isys) = pt1
                          end if
                       else if (kti(j) .eq. pt2) then
                          symm = 3.0d0
                          do k = 1, 3
                             nitors = nitors + 1
                             iitors(3,nitors) = ic
                             if (k .eq. 1) then
                                iitors(1,nitors) = ia
                                iitors(2,nitors) = ib
                                iitors(4,nitors) = id
                             else if (k .eq. 2) then
                                iitors(1,nitors) = ib
                                iitors(2,nitors) = id
                                iitors(4,nitors) = ia
                             else if (k .eq. 3) then
                                iitors(1,nitors) = id
                                iitors(2,nitors) = ia
                                iitors(4,nitors) = ib
                             end if
                             itors1(1,nitors) = ti1(1,j) / symm
                             itors1(2,nitors) = ti1(2,j)
                             itors2(1,nitors) = ti2(1,j) / symm
                             itors2(2,nitors) = ti2(2,j)
                             itors3(1,nitors) = ti3(1,j) / symm
                             itors3(2,nitors) = ti3(2,j)
                          end do
                          done = .true.
                          if (.not.is_find8(kti_sys(1),isys,kti(j)))
     $                       then
                             isys = isys + 1
                             kti_sys(isys) = pt2
                          end if
                       else if (kti(j) .eq. pt3) then
                          symm = 3.0d0
                          do k = 1, 3
                             nitors = nitors + 1
                             iitors(3,nitors) = ic
                             if (k .eq. 1) then
                                iitors(1,nitors) = ia
                                iitors(2,nitors) = ib
                                iitors(4,nitors) = id
                             else if (k .eq. 2) then
                                iitors(1,nitors) = ib
                                iitors(2,nitors) = id
                                iitors(4,nitors) = ia
                             else if (k .eq. 3) then
                                iitors(1,nitors) = id
                                iitors(2,nitors) = ia
                                iitors(4,nitors) = ib
                             end if
                             itors1(1,nitors) = ti1(1,j) / symm
                             itors1(2,nitors) = ti1(2,j)
                             itors2(1,nitors) = ti2(1,j) / symm
                             itors2(2,nitors) = ti2(2,j)
                             itors3(1,nitors) = ti3(1,j) / symm
                             itors3(2,nitors) = ti3(2,j)
                          end do
                          done = .true.
                          if (.not.is_find8(kti_sys(1),isys,kti(j)))
     $                       then
                             isys = isys + 1
                             kti_sys(isys) = pt3
                          end if
                       end if
                    end do
                 end if
                 if (.not. done) then
                    do j = 1, nti
                       if (kti(j) .eq. pt0) then
                          symm = 3.0d0
                          do k = 1, 3
                             nitors = nitors + 1
                             iitors(3,nitors) = ic
                             if (k .eq. 1) then
                                iitors(1,nitors) = ia
                                iitors(2,nitors) = ib
                                iitors(4,nitors) = id
                             else if (k .eq. 2) then
                                iitors(1,nitors) = ib
                                iitors(2,nitors) = id
                                iitors(4,nitors) = ia
                             else if (k .eq. 3) then
                                iitors(1,nitors) = id
                                iitors(2,nitors) = ia
                                iitors(4,nitors) = ib
                             end if
                             itors1(1,nitors) = ti1(1,j) / symm
                             itors1(2,nitors) = ti1(2,j)
                             itors2(1,nitors) = ti2(1,j) / symm
                             itors2(2,nitors) = ti2(2,j)
                             itors3(1,nitors) = ti3(1,j) / symm
                             itors3(2,nitors) = ti3(2,j)
                          end do
                          if (.not.is_find8(kti_sys(1),isys,kti(j)))
     $                       then
                             isys = isys + 1
                             kti_sys(isys) = pt0
                          end if
                       end if
                    end do
                 end if
              end if
           end do
        end if
c
c       find the cosine and sine of the phase angle for each torsion
c
        do i = 1, nitors
           angle = itors1(2,i) / radian
           itors1(3,i) = cos(angle)
           itors1(4,i) = sin(angle)
           angle = itors2(2,i) / radian
           itors2(3,i) = cos(angle)
           itors2(4,i) = sin(angle)
           angle = itors3(2,i) / radian
           itors3(3,i) = cos(angle)
           itors3(4,i) = sin(angle)
        end do
c
c       turn off the improper torsional potential if it is not used
c
        if (nitors .eq. 0)  use_imptor = .false.
      end if
      nti = maxnti
      do i = maxnti, 1, -1
         if (kti(i) .eq. -1)  nti = i - 1
      end do
      if (allocated(imptorglob)) deallocate(imptorglob)
      allocate (imptorglob(6*nbloc))
      nitorsloc = 0
      if (nti .ne. 0) then
         do i = 1, nloc
            iglob = glob(i)
            imptorcount = nbimptor(iglob)
            nitorsloc1 = 0
            if (n12(iglob) .eq. 3) then
               ia = i12(1,iglob)
               ib = i12(2,iglob)
               ic = iglob
               id = i12(3,iglob)
               ita = class(ia)
               itb = class(ib)
               itc = class(ic)
               itd = class(id)
c               size = 4
c               call numeral (ita,pa,size)
c               call numeral (itb,pb,size)
c               call numeral (itc,pc,size)
c               call numeral (itd,pd,size)
c               pt(1) = pa//pb//pc//pd
               call front_convert_base(ita,itb,itc,itd,0,pt(1))
c               pt(2) = pb//pa//pc//pd
               call front_convert_base(itb,ita,itc,itd,0,pt(2))
c               pt(3) = pa//pd//pc//pb
               call front_convert_base(ita,itd,itc,itb,0,pt(3))
c               pt(4) = pd//pa//pc//pb
               call front_convert_base(itd,ita,itc,itb,0,pt(4))
c               pt(5) = pb//pd//pc//pa
               call front_convert_base(itb,itd,itc,ita,0,pt(5))
c               pt(6) = pd//pb//pc//pa
               call front_convert_base(itd,itb,itc,ita,0,pt(6))
c               pt3 = zeros//zeros//pc//pd
               call front_convert_base(0,0,itc,ita,0,pt3)
c               pt2 = zeros//zeros//pc//pb
               call front_convert_base(0,0,itc,itb,0,pt2)
c               pt1 = zeros//zeros//pc//pa
               call front_convert_base(0,0,itc,ita,0,pt1)
c               pt0 = zeros//zeros//pc//zeros
               call front_convert_base(0,0,itc,0,0,pt0)
               symm = 1.0d0
               done = .false.
               do j = 1, nti
                  call back_convert_base(ita1,itb1,itc1,itd1,zero,
     $               kti(j))
c                  if (kti(j)(9:12) .eq. pc) then
                  if (itc1 .eq. itc) then
                     do k = 1, 6
                        if (kti(j) .eq. pt(k)) then
                           nitorsloc = nitorsloc + 1
                           nitorsloc1 = nitorsloc1 + 1
                           imptorglob(nitorsloc)=imptorcount +nitorsloc1
                           done = .true.
                        end if
                     end do
                  end if
               end do
               if (.not. done) then
                  do j = 1, nti
                     if (kti(j) .eq. pt1) then
                        do k = 1, 3
                           nitorsloc = nitorsloc + 1
                           nitorsloc1 = nitorsloc1 + 1
                           imptorglob(nitorsloc)=imptorcount +nitorsloc1
                           done = .true.
                        end do
                        done = .true.
                     else if (kti(j) .eq. pt2) then
                        do k = 1, 3
                           nitorsloc = nitorsloc + 1
                           nitorsloc1 = nitorsloc1 + 1
                           imptorglob(nitorsloc)=imptorcount +nitorsloc1
                        end do
                        done = .true.
                     else if (kti(j) .eq. pt3) then
                        do k = 1, 3
                           nitorsloc = nitorsloc + 1
                           nitorsloc1 = nitorsloc1 + 1
                           imptorglob(nitorsloc)=imptorcount +nitorsloc1
                        end do
                        done = .true.
                     end if
                  end do
               end if
               if (.not. done) then
                  do j = 1, nti
                     if (kti(j) .eq. pt0) then
                        do k = 1, 3
                           nitorsloc = nitorsloc + 1
                           nitorsloc1 = nitorsloc1 + 1
                           imptorglob(nitorsloc)=imptorcount +nitorsloc1
                        end do
                     end if
                  end do
               end if
            end if
         end do
      end if
      return
      end
