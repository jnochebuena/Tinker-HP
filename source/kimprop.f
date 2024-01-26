c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine kimprop  --  improper dihedral parameters  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "kimprop" assigns potential parameters to each improper
c     dihedral in the structure and processes any changed values
c
c
      subroutine kimprop(init)
      use atmlst
      use atmtyp
      use atoms
      use couple
      use domdec
      use improp
      use inform
      use iounit
      use keys
      use kiprop
      use potent
      use tors
      use utils
      implicit none
      integer i,j,k,ndi
      integer ia,ib,ic,id
      integer ita,itb,itc,itd
      integer ita1,itb1,itc1,itd1
      integer iglob,impropcount,niproploc1
      integer size,next,zero
      integer::isys=0
      integer*8 ipt
      real*8 tk,tv,symm
      logical header,done
c      character*4 pa,pb,pc,pd
c      character*12 zeros
c      character*16 blank
c      character*16 pt0,pt1
c      character*16 pt2,pt3
c      character*16 pt(6)
      integer*8 pt(6),pt0,pt1,pt2,pt3
      character*20 keyword
      character*120 record
      character*120 string
      logical init
c
c      blank = '                '
c      zeros = '000000000000'
      if (init) then
c
c     process keywords containing improper dihedral parameters
c
        header = .true.
        do i = 1, nkey
           next = 1
           record = keyline(i)
           call gettext (record,keyword,next)
           call upcase (keyword)
           if (keyword(1:9) .eq. 'IMPROPER ') then
              ia = 0
              ib = 0
              ic = 0
              id = 0
              tk = 0.0d0
              tv = 0.0d0
              string = record(next:120)
              read (string,*,err=10,end=10)  ia,ib,ic,id,tk,tv
   10         continue
c              size = 4
c              call numeral (ia,pa,size)
c              call numeral (ib,pb,size)
c              call numeral (ic,pc,size)
c              call numeral (id,pd,size)
c              pti = pa//pb//pc//pd
              call front_convert_base(ia,ib,ic,id,0,ipt)
              if (.not. silent) then
                 if (header) then
                    header = .false.
                    if (rank.eq.0) write (iout,20)
   20               format (/,' Additional Improper Dihedral',
     &                         ' Parameters :',
     &                      //,5x,'Atom Classes',20x,'K(ID)',
     &                         7x,'Angle',/)
                 end if
                 if (rank.eq.0) write (iout,30)  ia,ib,ic,id,tk,tv
   30            format (4x,4i4,10x,2f12.3)
              end if
              do j = 1, maxndi
c                 if (kdi(j).eq.blank .or. kdi(j).eq.pti) then
                 if (kdi(j).eq.-1 .or. kdi(j).eq.ipt) then
                    kdi(j) = ipt
                    dcon(j) = tk
                    tdi(j) = tv
                    goto 50
                 end if
              end do
              if (rank.eq.0) write (iout,40)
   40         format (/,' KIMPROP  --  Too many Improper Dihedral',
     &                   ' Parameters')
              abort = .true.
   50         continue
           end if
        end do
c
c       determine the total number of forcefield parameters
c
        ndi = maxndi
        do i = maxndi, 1, -1
           if (kdi(i) .eq. -1)  ndi = i - 1
        end do
c
c       assign improper dihedral parameters for each improper angle;
c       multiple symmetrical parameters are given partial weights
c
        niprop = 0
        if (ndi .ne. 0) then
           do i = 1, n
              if (n12(i) .eq. 3) then
                 ia = i
                 ib = i12(1,i)
                 ic = i12(2,i)
                 id = i12(3,i)
                 nbimprop(i) = niprop
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
c                 pt(2) = pa//pb//pd//pc
                 call front_convert_base(ita,itb,itd,itc,0,pt(2))
c                 pt(3) = pa//pc//pb//pd
                 call front_convert_base(ita,itc,itb,itd,0,pt(3))
c                 pt(4) = pa//pc//pd//pb
                 call front_convert_base(ita,itc,itd,itb,0,pt(4))
c                 pt(5) = pa//pd//pb//pc
                 call front_convert_base(ita,itd,itb,itc,0,pt(5))
c                 pt(6) = pa//pd//pc//pb
                 call front_convert_base(ita,itd,itc,itb,0,pt(6))
c                 pt3 = pa//pb//zeros//zeros
                 call front_convert_base(ita,itb,0,0,0,pt3)
c                 pt2 = pa//pc//zeros//zeros
                 call front_convert_base(ita,itc,0,0,0,pt2)
c                 pt1 = pa//pd//zeros//zeros
                 call front_convert_base(ita,itd,0,0,0,pt1)
                 call front_convert_base(ita,0,0,0,0,pt0)
c                 pt0 = pa//zeros
                 symm = 1.0d0
c                 if (pb.eq.pc .or. pb.eq.pd .or. pc.eq.pd)  symm = 2.0d0
c                 if (pb.eq.pc .and. pb.eq.pd .and. pc.eq.pd)  symm = 6.0d0
                 if (itb.eq.itc .or. itb.eq.itd .or. itc.eq.itd)  
     $             symm = 2.0d0
                 if (itb.eq.itc .and. itb.eq.itd .and. itc.eq.itd)  
     $             symm = 6.0d0
                 done = .false.
                 do j = 1, ndi
                   call back_convert_base(ita1,itb1,itc1,itd1,zero,
     $                 kdi(j))
                    if (ita1 .eq. ita) then
                       do k = 1, 6
                          if (kdi(j) .eq. pt(k)) then
                             niprop = niprop + 1
                             iiprop(1,niprop) = ia
                             if (k .eq. 1) then
                                iiprop(2,niprop) = ib
                                iiprop(3,niprop) = ic
                                iiprop(4,niprop) = id
                             else if (k .eq. 2) then
                                iiprop(2,niprop) = ib
                                iiprop(3,niprop) = id
                                iiprop(4,niprop) = ic
                             else if (k .eq. 3) then
                                iiprop(2,niprop) = ic
                                iiprop(3,niprop) = ib
                                iiprop(4,niprop) = id
                             else if (k .eq. 4) then
                                iiprop(2,niprop) = ic
                                iiprop(3,niprop) = id
                                iiprop(4,niprop) = ib
                             else if (k .eq. 5) then
                                iiprop(2,niprop) = id
                                iiprop(3,niprop) = ib
                                iiprop(4,niprop) = ic
                             else if (k .eq. 6) then
                                iiprop(2,niprop) = id
                                iiprop(3,niprop) = ic
                                iiprop(4,niprop) = ib
                             end if
                             kprop(niprop) = dcon(j) / symm
                             vprop(niprop) = tdi(j)
                             done = .true.
                             if (.not.is_find8(kdi_sys(1),isys,kdi(j)))
     $                          then
                                isys = isys + 1
                                kdi_sys(isys) = pt(k)
                             end if
                          end if
                       end do
                    end if
                 end do
                 if (.not. done) then
                    do j = 1, ndi
                       if (kdi(j) .eq. pt1) then
                          symm = 3.0d0
                          do k = 1, 3
                             niprop = niprop + 1
                             iiprop(1,niprop) = ia
                             if (k .eq. 1) then
                                iiprop(2,niprop) = ib
                                iiprop(3,niprop) = ic
                                iiprop(4,niprop) = id
                             else if (k .eq. 2) then
                                iiprop(2,niprop) = ic
                                iiprop(3,niprop) = id
                                iiprop(4,niprop) = ib
                             else if (k .eq. 3) then
                                iiprop(2,niprop) = id
                                iiprop(3,niprop) = ib
                                iiprop(4,niprop) = ic
                             end if
                             kprop(niprop) = dcon(j) / symm
                             vprop(niprop) = tdi(j)
                          end do
                          done = .true.
                          if (.not.is_find8(kdi_sys(1),isys,kdi(j)))
     $                       then
                             isys = isys + 1
                             kdi_sys(isys) = pt1
                          end if
                       else if (kdi(j) .eq. pt2) then
                          symm = 3.0d0
                          do k = 1, 3
                             niprop = niprop + 1
                             iiprop(1,niprop) = ia
                             if (k .eq. 1) then
                                iiprop(2,niprop) = ib
                                iiprop(3,niprop) = ic
                                iiprop(4,niprop) = id
                             else if (k .eq. 2) then
                                iiprop(2,niprop) = ic
                                iiprop(3,niprop) = id
                                iiprop(4,niprop) = ib
                             else if (k .eq. 3) then
                                iiprop(2,niprop) = id
                                iiprop(3,niprop) = ib
                                iiprop(4,niprop) = ic
                             end if
                             kprop(niprop) = dcon(j) / symm
                             vprop(niprop) = tdi(j)
                          end do
                          done = .true.
                          if (.not.is_find8(kdi_sys(1),isys,kdi(j)))
     $                       then
                             isys = isys + 1
                             kdi_sys(isys) = pt2
                          end if
                       else if (kdi(j) .eq. pt3) then
                          symm = 3.0d0
                          do k = 1, 3
                             niprop = niprop + 1
                             iiprop(1,niprop) = ia
                             if (k .eq. 1) then
                                iiprop(2,niprop) = ib
                                iiprop(3,niprop) = ic
                                iiprop(4,niprop) = id
                             else if (k .eq. 2) then
                                iiprop(2,niprop) = ic
                                iiprop(3,niprop) = id
                                iiprop(4,niprop) = ib
                             else if (k .eq. 3) then
                                iiprop(2,niprop) = id
                                iiprop(3,niprop) = ib
                                iiprop(4,niprop) = ic
                             end if
                             kprop(niprop) = dcon(j) / symm
                             vprop(niprop) = tdi(j)
                          end do
                          done = .true.
                          if (.not.is_find8(kdi_sys(1),isys,kdi(j)))
     $                       then
                             isys = isys + 1
                             kdi_sys(isys) = pt3
                          end if
                       end if
                    end do
                 end if
                 if (.not. done) then
                    do j = 1, ndi
                       if (kdi(j) .eq. pt0) then
                          symm = 3.0d0
                          do k = 1, 3
                             niprop = niprop + 1
                             iiprop(1,niprop) = ia
                             if (k .eq. 1) then
                                iiprop(2,niprop) = ib
                                iiprop(3,niprop) = ic
                                iiprop(4,niprop) = id
                             else if (k .eq. 2) then
                                iiprop(2,niprop) = ic
                                iiprop(3,niprop) = id
                                iiprop(4,niprop) = ib
                             else if (k .eq. 3) then
                                iiprop(2,niprop) = id
                                iiprop(3,niprop) = ib
                                iiprop(4,niprop) = ic
                             end if
                             kprop(niprop) = dcon(j) / symm
                             vprop(niprop) = tdi(j)
                          end do
                          if (.not.is_find8(kdi_sys(1),isys,kdi(j)))
     $                       then
                             isys = isys + 1
                             kdi_sys(isys) = pt0
                          end if
                       end if
                    end do
                 end if
              end if
              kdi_sys(0) = isys
           end do
        end if
c
c       turn off the improper dihedral potential if it is not used
c
        if (niprop .eq. 0)  use_improp = .false.
      end if
      ndi = maxndi
      do i = maxndi, 1, -1
         if (kdi(i) .eq. -1)  ndi = i - 1
      end do
      if (allocated(impropglob)) deallocate(impropglob)
      allocate (impropglob(6*nbloc))
      niproploc = 0
      if (ndi .ne. 0) then
         do i = 1, nloc
            iglob = glob(i)
            impropcount = nbimprop(iglob)
            if (n12(iglob) .eq. 3) then
               ia = iglob
               ib = i12(1,iglob)
               ic = i12(2,iglob)
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
c               pt(2) = pa//pb//pd//pc
               call front_convert_base(ita,itb,itd,itc,0,pt(2))
c               pt(3) = pa//pc//pb//pd
               call front_convert_base(ita,itc,itb,itd,0,pt(3))
c               pt(4) = pa//pc//pd//pb
               call front_convert_base(ita,itc,itd,itb,0,pt(4))
c               pt(5) = pa//pd//pb//pc
               call front_convert_base(ita,itd,itb,itc,0,pt(5))
c               pt(6) = pa//pd//pc//pb
               call front_convert_base(ita,itd,itc,itb,0,pt(6))
c               pt3 = pa//pb//zeros//zeros
               call front_convert_base(ita,itb,0,0,0,pt3)
c               pt2 = pa//pc//zeros//zeros
               call front_convert_base(ita,itc,0,0,0,pt2)
c               pt1 = pa//pd//zeros//zeros
               call front_convert_base(ita,itd,0,0,0,pt1)
c               pt0 = pa//zeros
               call front_convert_base(ita,0,0,0,0,pt0)
               done = .false.
               niproploc1 = 0
               do j = 1, ndi
                  call back_convert_base(ita1,itb1,itc1,itd1,zero,
     $               kdi(j))
                  if (ita1 .eq. ita) then
                     do k = 1, 6
                        if (kdi(j) .eq. pt(k)) then
                           niproploc = niproploc + 1
                           niproploc1 = niproploc1 + 1
                           impropglob(niproploc)=impropcount +niproploc1
                           done = .true.
                        end if
                     end do
                  end if
               end do
               if (.not. done) then
                  do j = 1, ndi
                     if (kdi(j) .eq. pt1) then
                        do k = 1, 3
                           niproploc = niproploc + 1
                           niproploc1 = niproploc1 + 1
                           impropglob(niproploc)=impropcount +niproploc1
                        end do
                        done = .true.
                     else if (kdi(j) .eq. pt2) then
                        do k = 1, 3
                           niproploc = niproploc + 1
                           niproploc1 = niproploc1 + 1
                           impropglob(niproploc)=impropcount +niproploc1
                        end do
                        done = .true.
                     else if (kdi(j) .eq. pt3) then
                        do k = 1, 3
                           niproploc = niproploc + 1
                           niproploc1 = niproploc1 + 1
                           impropglob(niproploc)=impropcount +niproploc1
                        end do
                        done = .true.
                     end if
                  end do
               end if
               if (.not. done) then
                  do j = 1, ndi
                     if (kdi(j) .eq. pt0) then
                        do k = 1, 3
                           niproploc = niproploc + 1
                           niproploc1 = niproploc1 + 1
                           impropglob(niproploc)=impropcount +niproploc1
                        end do
                     end if
                  end do
               end if
            end if
         end do
      end if
c
      return
      end
