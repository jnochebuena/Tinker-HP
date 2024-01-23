c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine echgtrn3  --  charge transfer energy & analysis  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "echgtrn3" calculates the charge transfer energy; also partitions
c     the energy among the atoms
c
c
      subroutine echgtrn3
      implicit none
c
c
c     choose method for summing over charge transfer interactions
c
      call echgtrn3c
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine echgtrn3c  --  neighbor list chgtrn analysis  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "echgtrn3c" calculates the charge transfer interaction energy
c     and also partitions the energy among the atoms using a pairwise
c     neighbor list
c
c
      subroutine echgtrn3c
      use action
      use analyz
      use atmtyp
      use atmlst
      use atoms
      use bound
      use chgpot
      use chgtrn
      use cell
      use couple
      use ctrpot
      use domdec
      use energi
      use group
      use inform
      use inter
      use iounit
      use molcul
      use mplpot
      use mpole
      use neigh
      use shunt
      use usage
      implicit none
      integer i,j,k,iipole,iglob,kglob
      integer ii,kk,kkk,kbis,kkpole
      real*8 e,f,fgrp
      real*8 r,r2,r3
      real*8 r4,r5
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 chgi,chgk
      real*8 chgik
      real*8 alphai,alphak
      real*8 expi,expk
      real*8 expik
      real*8 taper
      real*8, allocatable :: mscale(:)
      logical proceed,usei
      logical header,huge
      character*10 mode
c
c     zero out the charge transfer energy and partitioning terms
c
      nect = 0
      ect = 0.0d0
      aect = 0.0d0
      if (nct .eq. 0)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
c
c     initialize connected atom exclusion coefficients
c
      mscale = 1.0d0
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'CHGTRN'
      call switch (mode)
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. nct.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Charge Transfer Interactions :',
     &           //,' Type',14x,'Atom Names',15x,'Distance',
     &              8x,'Energy',/)
      end if
c
c     compute the charge transfer energy
c
      do ii = 1, npolelocnl
         iipole = poleglobnl(ii)
         iglob = ipole(iipole)
         i = loc(iglob)
         xi = x(iglob)
         yi = y(iglob)
         zi = z(iglob)
         chgi = chgct(iipole)
         alphai = dmpct(iipole)
         if (alphai .eq. 0.0d0)  alphai = 100.0d0
         usei = use(iglob)
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(iglob)
            mscale(i12(j,iglob)) = m2scale
         end do
         do j = 1, n13(iglob)
            mscale(i13(j,iglob)) = m3scale
         end do
         do j = 1, n14(iglob)
            mscale(i14(j,iglob)) = m4scale
         end do
         do j = 1, n15(iglob)
            mscale(i15(j,iglob)) = m5scale
         end do
c
c     evaluate all sites within the cutoff distance
c
         do kkk = 1, nelst(ii)
            kkpole = elst(kkk,ii)
            kglob = ipole(kkpole)
            kbis = loc(kglob)
            proceed = (usei .or. use(kglob))
            if (proceed) then
               xr = x(kglob) - xi
               yr = y(kglob) - yi
               zr = z(kglob) - zi
               if (use_bounds)  call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  chgk = chgct(kkpole)
                  alphak = dmpct(kkpole)
                  if (alphak .eq. 0.0d0)  alphak = 100.0d0
                  if (ctrntyp .eq. 'SEPARATE') then
                     expi = exp(-alphai*r)
                     expk = exp(-alphak*r)
                     e = -chgi*expk - chgk*expi
                  else
                     chgik = sqrt(abs(chgi*chgk))
                     expik = exp(-0.5d0*(alphai+alphak)*r)
                     e = -chgik * expik
                  end if
                  e = f * e * mscale(kglob)
c
c     use energy switching if near the cutoff distance
c
                  if (r2 .gt. cut2) then
                     r3 = r2 * r
                     r4 = r2 * r2
                     r5 = r2 * r3
                     taper = c5*r5 + c4*r4 + c3*r3
     &                          + c2*r2 + c1*r + c0
                     e = e * taper
                  end if
c     
c     increment the overall charge transfer energy components
c
                  if (e .ne. 0.0d0) then
                     nect = nect + 1
                     ect = ect + e
                     aect(i) = aect(i) + 0.5d0*e
                     aect(kbis) = aect(kbis) + 0.5d0*e
                  end if
c
c     increment the total intermolecular energy
c
                  if (molcule(iglob) .ne. molcule(kglob)) then
                     einter = einter + e
                  end if
c
c     print a message if the energy of this interaction is large
c
                  huge = (e .gt. 10.0d0)
                  if ((debug.and.e.ne.0.0d0)
     &                  .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,20)
   20                   format (/,' Individual Charge Transfer',
     &                             ' Interactions :',
     &                          //,' Type',14x,'Atom Names',
     &                             15x,'Distance',8x,'Energy',/)
                     end if
                     write (iout,30)  iglob,name(iglob),kglob,
     &                 name(kglob),r,e
   30                format (' ChgTrn',4x,2(i7,'-',a3),
     &                          9x,f10.4,2x,f12.4)
                  end if
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(iglob)
            mscale(i12(j,iglob)) = 1.0d0
         end do
         do j = 1, n13(iglob)
            mscale(i13(j,iglob)) = 1.0d0
         end do
         do j = 1, n14(iglob)
            mscale(i14(j,iglob)) = 1.0d0
         end do
         do j = 1, n15(iglob)
            mscale(i15(j,iglob)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      return
      end
