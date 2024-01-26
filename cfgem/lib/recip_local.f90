subroutine FH_RECIP_LOC_pme_grid_mult( &
                           nfft1,nfft2,nfft3,nfftdim1, &
                           expon,recip,volume,prefac1,prefac2,prefac3, &
                           pme_grid_multiplier)

   implicit none

   integer, intent(in) :: nfft1,nfft2,nfft3,nfftdim1
   double precision, intent(in) :: expon, recip(3,3),volume, &
                                   prefac1(2,*), prefac2(2,*),prefac3(2,*)
   double precision, intent(out) :: pme_grid_multiplier(2,nfft3,nfftdim1,nfft2)

   integer :: k1,k2,k3,m1,m2,m3,k10,nf1,nf2,nf3
   double precision :: pi,fac,mhat1,mhat2,mhat3,msq, &
                       term,tmp1r,tmp1i,tmp2r,tmp2i

   pi = 3.14159265358979323846d0
   fac = pi*pi / expon

   nf1 = nfft1/2
   if ( 2*nf1 < nfft1 )nf1 = nf1+1
   nf2 = nfft2/2
   if ( 2*nf2 < nfft2 )nf2 = nf2+1
   nf3 = nfft3/2
   if ( 2*nf3 < nfft3 )nf3 = nf3+1

   do k2 = 1, nfft2
      m2 = k2 - 1
      if ( k2 > nf2 )m2 = k2 - 1 - nfft2
      do k3 = 1,nfft3
         m3 = k3 - 1
         if ( k3 > nf3 )m3 = k3 - 1 - nfft3
         k10 = 1
         ! need (1,1,1) case also
         !if(k3+k2 == 2) k10 = 2
         do k1 = k10, nf1+1
            m1 = k1 - 1
            if ( k1 > nf1 )m1 = k1 - 1 - nfft1
            mhat1 = recip(1,1)*m1+recip(1,2)*m2+recip(1,3)*m3
            mhat2 = recip(2,1)*m1+recip(2,2)*m2+recip(2,3)*m3
            mhat3 = recip(3,1)*m1+recip(3,2)*m2+recip(3,3)*m3
            msq = mhat1*mhat1+mhat2*mhat2+mhat3*mhat3
            term = exp(-fac*msq) / volume
            tmp1r = prefac1(1,k1)*prefac2(1,k2) - prefac1(2,k1)*prefac2(2,k2)
            tmp1i = prefac1(1,k1)*prefac2(2,k2) + prefac1(2,k1)*prefac2(1,k2)
            tmp2r = tmp1r*prefac3(1,k3) - tmp1i*prefac3(2,k3)
            tmp2i = tmp1r*prefac3(2,k3) + tmp1i*prefac3(1,k3)
            pme_grid_multiplier(1,k3,k1,k2) = tmp2r*term
            pme_grid_multiplier(2,k3,k1,k2) = tmp2i*term
         enddo
      enddo
   enddo
end subroutine FH_RECIP_LOC_pme_grid_mult
!-----------------------------------------------------------
subroutine dump_FT_density(nf,nfft1,nfft2,nfft3,nfftdim1, &
                           grid)

   implicit none

   integer, intent(in) :: nf,nfft1,nfft2,nfft3,nfftdim1
   double precision, intent(in) :: grid(2,nfft3,nfftdim1,nfft2)

   integer :: k1,k2,k3,k10,nf1
   nf1 = nfft1/2
   if ( 2*nf1 < nfft1 )nf1 = nf1+1
   do k2 = 1, nfft2
      do k3 = 1,nfft3
         k10 = 1
         ! need (1,1,1) case also
         !if ( k3+k2 == 2 )k10 = 2
         do k1 = k10, nf1+1
            write(nf,'(3i4,1x,2e20.12)') &
              k3,k1,k2,grid(1,k3,k1,k2),grid(2,k3,k1,k2)
         enddo
      enddo
   enddo
end subroutine dump_FT_density
!-----------------------------------------------------------
