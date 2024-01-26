c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     "empole1" : driver for calculation of the multipole and dipole polarization
c     energy and derivatives with respect to Cartesian coordinates
c
c
      subroutine empoleGEM
      use domdec
      use energi
      use potent
      use mpi
      implicit none
      real*8 time0, time1
c
c     choose the method for summing over multipole interactions
c
      time0 = mpi_wtime()
      call empole1cGEM
      time1 = mpi_wtime()
c      write(*,*) 'time empole1 = ',time1-time0
c
c     zero out energy and derivative terms which are not in use
c
      em = 0.0d0

      return
      end
c
c
c     "empole1c" calculates the multipole energy and derivatives
c     with respect to Cartesian coordinates using particle mesh Ewald
c     summation and a neighbor list
c
c
      subroutine empole1cGEM
      use sizes
      use atmlst
      use atoms
      use boxes
      use chgpot
      use deriv
      use domdec
      use energi
      use ewald
      use math
      use mpole
      use potent
      use timestat
      use virial
      use mpi
      implicit none
      integer i,j,ii
      integer iipole,iglob,ierr
      real*8 e,f
      real*8 term,fterm
      real*8 cii,dii,qii
      real*8 xd,yd,zd
      real*8 xq,yq,zq
      real*8 xv,yv,zv,vterm
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 xdfield,ydfield
      real*8 zdfield
      real*8 trq(3),frcx(3)
      real*8 frcy(3),frcz(3)
      real*8 time0,time1
c
c
c     zero out the atomic multipole energy and derivatives
c
      em = 0.0d0
      dem = 0.0d0
      if (npole .eq. 0)  return
c
c     set the energy unit conversion factor
c
      f = electric / dielec
c
c     check the sign of multipole components at chiral sites
c
      call chkpole(.false.)
c
c     rotate the multipole components into the global frame
c
      call rotpole

      return
      end
