subroutine get_fftdims(nfft1,nfft2,nfft3, &
      nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork, &
      sizfftab,sizffwrk,out_lun)

   implicit none

   integer nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3, &
         nfftable,nffwork,sizfftab,sizffwrk
   integer n,nfftmax
   integer out_lun

   nfftdim1 = nfft1
   n = nfft1/2
   if ( nfft1 == 2*n )then
      nfftdim1 = nfft1/2+2
   else
      write(out_lun,*)'GET_FFTDIMS (ew_fft.f)', &
               'For RealComplex FFT nfft1 must be even'
      stop
   end if
   nfftdim2 = nfft2
   n = nfft2/2
   if ( nfft2 == 2*n )nfftdim2 = nfft2+1
   nfftdim3 = nfft3
   n = nfft3/2
   if ( nfft3 == 2*n )nfftdim3 = nfft3+1
   nfftmax = max(nfft1,nfft2,nfft3)

   !-----------------------------------------------------------

   nfftable = 4*nfftmax + 15
   sizfftab = 3*nfftable

   !       space work for the individual 1D ffts

   nffwork = 2*nfftdim2

   nffwork = nffwork + nfftdim1*nfft3*2*(nfft2+ 1)
   
   sizffwrk  = nffwork

end subroutine get_fftdims
!-----------------------------------------------------------
subroutine fft_setup_gem(array,fftable,ffwork, &
      nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3, &
      nfftable,nffwork,out_lun)

   implicit none
   
   double precision  array(*),fftable(*),ffwork(*)
   integer nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3
   integer nfftable,nffwork
   integer out_lun

   double precision  alpha,beta,tmpy
   integer isign
   double precision  scale
   
   isign = 0
   scale = 1.d0
   call fft3d0rc(isign,nfft1,nfft2,nfft3,scale,array, &
         nfftdim1,nfftdim2,fftable, &
         ffwork,tmpy,alpha,beta,out_lun )

   return

end subroutine fft_setup_gem
!-----------------------------------------------------------
!-------------------------------------------------------------------
!     --- FFT_BACK RC ---


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine fft_backrc here]
subroutine fft_backrc(array,fftable,ffwork, &
      nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3, &
      nfftable,nffwork, tmpy, alpha,beta, out_lun)

   implicit none

   double precision  array(*),fftable(*),ffwork(*)
   integer nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3
   integer nfftable,nffwork
   double precision  tmpy(*), &
         alpha(*),beta(*)
   integer out_lun

   integer isign
   double precision  scale
   
   isign = -1
   scale = 1.d0
   call fft3d0rc(isign,nfft1,nfft2,nfft3,scale,array, &
         nfftdim1,nfftdim2,fftable, &
         ffwork,tmpy,alpha,beta,out_lun )

   return
end subroutine fft_backrc 

!-------------------------------------------------------------------

!     --- FFT_FORWARD RC---


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine fft_forwardrc here]
subroutine fft_forwardrc(array,fftable,ffwork, &
      nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3, &
      nfftable,nffwork,tmpy,alpha,beta,out_lun)

   implicit none
   
   double precision  array(*),fftable(*),ffwork(*)
   integer nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3
   integer nfftable,nffwork
   double precision  tmpy(*), &
         alpha(*),beta(*)
   integer out_lun

   integer isign
   double precision  scale
   
   isign = 1
   scale = 1.d0
   call fft3d_zxyrc(isign,nfft1,nfft2,nfft3,scale,array, &
         nfftdim1,nfftdim2,array,nfftdim1,nfftdim2,fftable, &
         ffwork,tmpy,alpha,beta,out_lun )
   
end subroutine fft_forwardrc
!-----------------------------------------------------------
subroutine fft3d0rc(isign,n1,n2,n3, scale, &
      x, ldx,ldx2, table, work, tmpy, alpha, beta, out_lun )

   implicit none

   double precision,parameter :: one = 1.d0

   integer isign
   integer ldx,ldx2, n1, n2, n3
   double precision  x(*)
   double precision  table(*), work(*), scale
   double precision  tmpy(*),alpha(*),beta(*)
   integer           out_lun

   integer i, k, k0, ks
   integer j, ja, ja00, jz, jz0, jz00, jj, jidx, jtask

   integer mxyslabs,mxzslabs,ntxyslab,ntxzslab
   integer indz

   integer my_first_slab_xy
   integer my_first_slab_xz
   integer n1x

   
   !--------------------------------------------------------------
   !         startup: initialize tables for all three dimensions
   

   mxyslabs=n3
   mxzslabs=n2
   ntxyslab = (ldx) * ldx2 *2
   ntxzslab = (ldx) * n3 *2
   my_first_slab_xy = 0
   my_first_slab_xz = 0
   indz = 2*n2
   n1x = n1/2
   scale = one

   if(isign == 0)then

      call fft2drc(0,n1,n2,one, x,ldx, table,work, &
            tmpy,alpha,beta,out_lun)
      call cfgem_cffti(n3,table(4*n1+15 + 4*n2+15 + 1))

      return
   end if

   !-----------------------------------------------------------------------
   ! each PE should do their 2D ffts now
   !-----------------------------------------------------------------------
   do j = 1, mxyslabs
      jj = (j-1)*ntxyslab+1
      call fft2drc(isign, n1, n2, one, &
            x(jj),ldx, table, work, tmpy, alpha, beta,out_lun)
   end do

   !mfc ????  WHAT is this for???
   call UTIL_zero_real_array(work(indz+1),n3*ldx*mxzslabs)
   !mfc ????  WHAT was that for???

   do ks = 0,mxzslabs-1
      ja00 = (my_first_slab_xz+ks)*ldx*2
      jz00 = ks*ntxzslab + 2*my_first_slab_xy
      do j = 0,mxyslabs-1
         jz = jz00 + 2*j +1
         ja = ja00 + j*ntxyslab+1
         do i = 0,n1x
            work(indz+jz+i*n3*2) = x(ja+i*2)
            work(indz+jz+1+i*n3*2) = x(ja+i*2+1)
         end do
      end do
   end do

   !-----------------------------------------------------------------------
   !         END of TRANSPOSE
   !-----------------------------------------------------------------------
   !    Now do Z-FFTs
   !-----------------------------------------------------------------------
   do k = 0,mxzslabs-1
      k0 = k*ntxzslab
      do j = 0, n1x
         jidx=k0 + j*n3*2 +1
         call cfgem_cfftf(n3,work(indz+jidx),table(4*n1+15 + 4*n2+15 +1))
      end do
   end do
   !*************************************************
   !     Leave it   DISTRIBUTED AS Z-X SLABS
   !*************************************************
   do k = 1,ntxzslab*mxzslabs
      x(k) = work(indz+k)
   end do

   return
end subroutine fft3d0rc 
!----------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine fft3d_zxyrc here]

subroutine fft3d_zxyrc(isign,n1,n2,n3, scale, &
      x, ldx,ldx2, y, ldy, ldy2, &
      table, work,tmpy,alpha,beta,out_lun )

   implicit none

   !****************************************************************
   
   !  isign = 1 for forward fft
   !        =-1 for reverse fft
   !        = 0 initializes table for all three dimensions
   !  n1, n2, n3   dimensions of the fft to be performed
   !  scale  data will be scaled by this number on return
   !         see ccfgem_cfft3d from Cray for details of how to use this one
   !  x, ldx, ldx2  complex 3-d array
   !                the input array with declared dimensions ldx and ldx2
   !                in the first two dimensions
   !  y, ldy, ldy2  complex 3-d array output 3D array
   !  table  real size 2*(n1+n2+n3) (for Cray CCFFT)
   !  work   workspace size 4*( max(n1,n2,n3) ) (for Cray CCFFT)
   !                     + ldx * ldx2 * 2 * (n3/numtasks + 1)
   !                     + ldx * n3   * 2 * (n2/numtasks + 1)
   !  isys   use 0 (read cray docs on ccfft); currently unused.
   !*************************************************************

   integer ldx,ldx2, ldy, ldy2, n1, n2, n3
   double precision  x(*), y(*)
   double precision  table(*), work(*), scale
   double precision  alpha(*),beta(*),tmpy(*)
   integer out_lun
   integer k, k0, ks
   integer ja, ja00, jz, jz0, jz00, jj, jidx, jtask
   integer  i, j, isign, ndim

   integer mxyslabs,mxzslabs,ntxyslab,ntxzslab
   integer indz,mytaskid
   double precision, parameter :: one=1.d0

   integer my_first_slab_xy
   integer my_first_slab_xz
   integer n1x, loop_limit

   !--------------------------------------------------------------
   !         startup: no initialization possible in this routine
   
   if(isign == 0)then
      write(out_lun,*)"fft3d_zxy ERROR, cannot do an isign of 0, I QUIT"
      stop
   end if

   if (isign /= 1 ) then
      write(out_lun,*)'isign for 2nd fft should be 1'
      stop
   end if

   mxyslabs=n3
   mxzslabs=n2
   ntxyslab = ldx * ldx2 *2
   ntxzslab = ldx * n3 *2
   my_first_slab_xy = 0
   my_first_slab_xz = 0
   indz = 2*n2
   n1x=n1/2

   !***********************************************************************
   !       DISTRIBUTED AS Z-X SLABS, put the data into z area of work
   !***********************************************************************
   loop_limit = ntxzslab*mxzslabs
   do k = 1,loop_limit
      work(indz+k) = x(k)
   end do
   !***********************************************************************
   !**** DO Z FFTs NOW **************************************************
   !***********************************************************************
   do k = 0,mxzslabs-1
      k0 = k*ntxzslab
      !           do j = 0, n1-1
      do j = 0, n1x
         jidx=k0 + j*n3*2 +1

         call cfgem_cfftb(n3,work(indz+jidx),table(4*n1+15 + 4*n2+15 +1))
      end do
   end do

   do ks = 0,mxzslabs-1
      ja00 = (my_first_slab_xz+ks)*ldx*2
      jz00 = ks*ntxzslab + 2*my_first_slab_xy
      do j = 0,mxyslabs-1
         jz = jz00 + 2*j +1
         ja = ja00 + j*ntxyslab+1
         do i = 0,n1x
            x(ja+i*2) = work(indz+jz+i*n3*2)
            x(ja+i*2+1) = work(indz+jz+1+i*n3*2)
         end do
      end do
   end do


   !-----------------------------------------------------------------------
   ! each PE should do their 2D ffts now
   !-----------------------------------------------------------------------
   
   do j = 1, mxyslabs
      jj = (j-1)*ntxyslab+1
      call fft2drc(isign, n1, n2, one, &
            x(jj),ldx, table, work, tmpy, alpha, beta,out_lun)
   end do

end subroutine fft3d_zxyrc 
!----------------------------------------------------------------
!************ 2D  FFT real-complex-real ***********************

!**************************************************************
!     pubfft implementation for REAL-to-COMPLEX fft
!**************************************************************
!      subroutine fft2drc
!   Complex version May 22 1995
!**************************************************************


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine fft2drc here]
subroutine fft2drc(isign,n1, n2, scale, x, ldx, &
      table, work, tmpy, alpha, beta, out_lun)
   
   !**************************************************************
   !     1 PE non-distributed 2 dim. FFT
   !         calls a 1 dim FFT (CCFFT for the 1)
   
   !         Author: Michael F. Crowley
   !                 Pittsburgh Supercomputing Center
   !                 Oct 20, 1994
   !**************************************************************
   
   implicit none

   integer isign,n1,n2,ldx
   double precision  scale
   double precision  x(2, 0:ldx-1, 0:n2-1)
   double precision  work(*), table(*)
   double precision  tmpy(2,0:ldx-1)
   double precision  alpha(0:n1),beta(0:n1)
   double precision, parameter :: pi=3.14159265358979323846d0
   double precision, parameter :: TWOPI=2.d0*pi
   double precision, parameter :: half = 0.5d0, zero=0.d0
   integer, intent(in)         :: out_lun

   integer i, idx, idy, j, jy
   integer n1rc,kr,kkr,ki,kki,idt,n1x
   integer j1,j2,j3,j4,j1rl,j1im,j2rl,j2im,j3rl,j3im,j4rl,j4im,jdx
   integer istar
   double precision  a,b,c,d,pi2n,theta

   !---------------------------------------------------------------

   pi2n=TWOPI/n1
   n1x=n1/2

   !=================================================
   !             initialize fft tables
   
   if(isign == 0)then
      if(mod(n1,2) /= 0)then
         write(out_lun,*)" NEED factor 2 value for nfft1 for RC fft"
         stop
      end if

      call cfgem_cffti(n1x,table)
      call cfgem_cffti(n2,table(4*n1+15 + 1))

      return
   end if

   !-----------------------------------------------------
   do i=0,n1x-1
      theta=pi2n*i
      alpha(i) = cos(theta)
      beta(i)  = sin(theta)
   end do
   !---------------------------------------------------------------
   !       Backward fft real to complex
   !---------------------------------------------------------------
   if(isign == -1)then
      !---------------------------------------------------------------
      !  First the x direction, the data is already contiguous
      
      !-----------------------------------------------------------------------
      do j = 0,n2-1
         do i = 0, n1x-1
            tmpy(1,i)=x(1,i,j)
            tmpy(2,i)=x(2,i,j)
         end do
         call cfgem_cfftf(n1x, tmpy(1,0), table)

         do i = 1, n1x-1
            a =  half*(tmpy(1,i)+tmpy(1,n1x-i)) ! Real F even
            b =  half*(tmpy(2,i)-tmpy(2,n1x-i)) ! Imag F even
            c =  half*(tmpy(2,i)+tmpy(2,n1x-i)) ! Real F odd
            d = -half*(tmpy(1,i)-tmpy(1,n1x-i)) ! Imag F odd
            x(1,i,j) = a + alpha(i)*c + beta(i)*d
            x(2,i,j) = b + alpha(i)*d - beta(i)*c
         end do
         !--------------------------------------------
         !     DC and nyquist
         x(1,0,j)  =tmpy(1,0)+tmpy(2,0)
         x(2,0,j)  =zero
         x(1,n1x,j)=tmpy(1,0)-tmpy(2,0)
         x(2,n1x,j)=zero
      end do
      
      !-----------------------------------------------------------------------
      !     Now in the y direction, the data is in y now and
      !     we will put it into a contiguous 1D array first, transform,
      !     then put it back.
      !     ibegw should be adjusted to be thesize of the work
      !     area necessary for the machine specific fft.
      !     for pubfft, there is no work area used so ibegw=0
      
      do j1 = 0, n1x
         i=2*j1+1
         istar=2*(n1-j1)
         do j2 = 0, n2-1
            j=2*j2+1
            work(j)   = x(1,j1,j2)
            work(j+1) = x(2,j1,j2)
         end do
         call cfgem_cfftf(n2, work(1), table(4*n1+16))
         do j2 = 0, n2-1
            j=2*j2+1
            x(1,j1,j2) = work(j)
            x(2,j1,j2) = work(j+1)
         end do
      end do
      !---------------------------------------------------------------
   else
      !---------------------------------------------------------------
      !     Forward fft complex to real
      !---------------------------------------------------------------
      !-----------------------------------------------------------------------
      !  Now in the y direction, the data is in y now and
      !     we will put it into a contiguous 1D array first, transform,
      !           then put it back.
      !     ibegw should be adjusted to be thesize of the work
      !           area necessary for the machine specific fft.
      !           for pubfft, there is no work area used so ibegw=0
      
      
      
      do i = 0, n1x
         do j = 0, n2-1

            work(2*j+1) = x(1,i,j)
            work(2*j+2) = x(2,i,j)
         end do
         call cfgem_cfftb(n2, work(1), table(4*n1+16))
         do j = 0, n2-1
            x(1,i,j) = work(2*j+1)
            x(2,i,j) = work(2*j+2)
         end do
      end do

      do j = 0,n2-1
         do i = 1, n1x-1
            a =  (x(1,i,j)+x(1,n1x-i,j)) ! Real F even
            b =  (x(2,i,j)-x(2,n1x-i,j)) ! Imag F even
            c =  (x(2,i,j)+x(2,n1x-i,j)) ! F odd contrib
            d =  (x(1,i,j)-x(1,n1x-i,j)) ! F odd contrib
            tmpy(1,i) = a - alpha(i)*c - beta(i)*d
            tmpy(2,i) = b + alpha(i)*d - beta(i)*c
         end do
         tmpy(1,0) = (x(1,0,j)+x(1,n1x,j))
         tmpy(2,0) = (x(1,0,j)-x(1,n1x,j))
         call cfgem_cfftb(n1x, tmpy(1,0), table)
         do i = 0, n1x-1
            x(1,i,j)=tmpy(1,i)
            x(2,i,j)=tmpy(2,i)
         end do
      end do
   end if  ! (isign == -1)

end subroutine fft2drc 
!--------------------------------------------------------------------
