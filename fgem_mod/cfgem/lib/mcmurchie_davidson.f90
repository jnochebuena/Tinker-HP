
!-----------------------------------------------------------------
!-----------------------------------------------------------------
! SETUP SUBROUTINES
!-----------------------------------------------------------------
!-----------------------------------------------------------------

!----------------------------------------------------------
subroutine FH_MCMUR_DAV_init()

   implicit none

   include "interact_type.fh"
   ! initialize pointers
   call FH_MCMUR_DAV_tuv_parts()
   call FH_MCMUR_DAV_order_from_level()
   call FH_MCMUR_DAV_fill_grad_ptrs()
end subroutine FH_MCMUR_DAV_init
!----------------------------------------------------------
subroutine FH_MCMUR_DAV_tuv_parts()

  implicit none

  include "mpole_index.fh"
  include "direct_pointers.fh"

! charge; no derivs
  t_part(IND_000) = 0; u_part(IND_000) = 0; v_part(IND_000) = 0
! dipoles
  t_part(IND_100) = 1; u_part(IND_100) = 0; v_part(IND_100) = 0
  t_part(IND_010) = 0; u_part(IND_010) = 1; v_part(IND_010) = 0
  t_part(IND_001) = 0; u_part(IND_001) = 0; v_part(IND_001) = 1
! quadrupoles
  t_part(IND_200) = 2; u_part(IND_200) = 0; v_part(IND_200) = 0
  t_part(IND_020) = 0; u_part(IND_020) = 2; v_part(IND_020) = 0
  t_part(IND_002) = 0; u_part(IND_002) = 0; v_part(IND_002) = 2
  t_part(IND_110) = 1; u_part(IND_110) = 1; v_part(IND_110) = 0
  t_part(IND_101) = 1; u_part(IND_101) = 0; v_part(IND_101) = 1
  t_part(IND_011) = 0; u_part(IND_011) = 1; v_part(IND_011) = 1
! octupoles
  t_part(IND_300) = 3; u_part(IND_300) = 0; v_part(IND_300) = 0
  t_part(IND_030) = 0; u_part(IND_030) = 3; v_part(IND_030) = 0
  t_part(IND_003) = 0; u_part(IND_003) = 0; v_part(IND_003) = 3
  t_part(IND_210) = 2; u_part(IND_210) = 1; v_part(IND_210) = 0
  t_part(IND_201) = 2; u_part(IND_201) = 0; v_part(IND_201) = 1
  t_part(IND_120) = 1; u_part(IND_120) = 2; v_part(IND_120) = 0
  t_part(IND_021) = 0; u_part(IND_021) = 2; v_part(IND_021) = 1
  t_part(IND_102) = 1; u_part(IND_102) = 0; v_part(IND_102) = 2
  t_part(IND_012) = 0; u_part(IND_012) = 1; v_part(IND_012) = 2
  t_part(IND_111) = 1; u_part(IND_111) = 1; v_part(IND_111) = 1
! hexadecapoles
  t_part(IND_400) = 4; u_part(IND_400) = 0; v_part(IND_400) = 0
  t_part(IND_040) = 0; u_part(IND_040) = 4; v_part(IND_040) = 0
  t_part(IND_004) = 0; u_part(IND_004) = 0; v_part(IND_004) = 4
  t_part(IND_310) = 3; u_part(IND_310) = 1; v_part(IND_310) = 0
  t_part(IND_301) = 3; u_part(IND_301) = 0; v_part(IND_301) = 1
  t_part(IND_130) = 1; u_part(IND_130) = 3; v_part(IND_130) = 0
  t_part(IND_031) = 0; u_part(IND_031) = 3; v_part(IND_031) = 1
  t_part(IND_103) = 1; u_part(IND_103) = 0; v_part(IND_103) = 3
  t_part(IND_013) = 0; u_part(IND_013) = 1; v_part(IND_013) = 3
  t_part(IND_220) = 2; u_part(IND_220) = 2; v_part(IND_220) = 0
  t_part(IND_202) = 2; u_part(IND_202) = 0; v_part(IND_202) = 2
  t_part(IND_022) = 0; u_part(IND_022) = 2; v_part(IND_022) = 2
  t_part(IND_211) = 2; u_part(IND_211) = 1; v_part(IND_211) = 1
  t_part(IND_121) = 1; u_part(IND_121) = 2; v_part(IND_121) = 1
  t_part(IND_112) = 1; u_part(IND_112) = 1; v_part(IND_112) = 2
!  5th order
  t_part(IND_500) = 5; u_part(IND_500) = 0; v_part(IND_500) = 0
  t_part(IND_050) = 0; u_part(IND_050) = 5; v_part(IND_050) = 0
  t_part(IND_005) = 0; u_part(IND_005) = 0; v_part(IND_005) = 4
  t_part(IND_410) = 4; u_part(IND_410) = 1; v_part(IND_410) = 0
  t_part(IND_401) = 4; u_part(IND_401) = 0; v_part(IND_401) = 1
  t_part(IND_140) = 1; u_part(IND_140) = 4; v_part(IND_140) = 0
  t_part(IND_041) = 0; u_part(IND_041) = 4; v_part(IND_041) = 1
  t_part(IND_104) = 1; u_part(IND_104) = 0; v_part(IND_104) = 4
  t_part(IND_014) = 0; u_part(IND_014) = 1; v_part(IND_014) = 4
  t_part(IND_320) = 3; u_part(IND_320) = 2; v_part(IND_320) = 0
  t_part(IND_302) = 3; u_part(IND_302) = 0; v_part(IND_302) = 2
  t_part(IND_230) = 2; u_part(IND_230) = 3; v_part(IND_230) = 0
  t_part(IND_032) = 0; u_part(IND_032) = 3; v_part(IND_032) = 2
  t_part(IND_203) = 2; u_part(IND_203) = 0; v_part(IND_203) = 3
  t_part(IND_023) = 0; u_part(IND_023) = 2; v_part(IND_023) = 3
  t_part(IND_311) = 3; u_part(IND_311) = 1; v_part(IND_311) = 1
  t_part(IND_131) = 1; u_part(IND_131) = 3; v_part(IND_131) = 1
  t_part(IND_113) = 1; u_part(IND_113) = 1; v_part(IND_113) = 3
  t_part(IND_221) = 2; u_part(IND_221) = 2; v_part(IND_221) = 1
  t_part(IND_212) = 2; u_part(IND_212) = 1; v_part(IND_212) = 2
  t_part(IND_122) = 1; u_part(IND_122) = 2; v_part(IND_122) = 2
   
end subroutine FH_MCMUR_DAV_tuv_parts
!------------------------------------------------------------------
subroutine FH_MCMUR_DAV_order_from_level()

   implicit none

   include "direct_pointers.fh"

   hermite_order_from_level(0) = 1
   hermite_order_from_level(1) = 4
   hermite_order_from_level(2) = 10
   hermite_order_from_level(3) = 20
   hermite_order_from_level(4) = 35
   hermite_order_from_level(5) = 56
end subroutine FH_MCMUR_DAV_order_from_level
!------------------------------------------------------------------
subroutine FH_MCMUR_DAV_fill_grad_ptrs()

   implicit none

   include "direct_pointers.fh"
   include "mpole_index.fh"

  p_xder(IND_000)=IND_100; p_yder(IND_000)=IND_010; p_zder(IND_000)=IND_001
  p_xder(IND_100)=IND_200; p_yder(IND_100)=IND_110; p_zder(IND_100)=IND_101
  p_xder(IND_010)=IND_110; p_yder(IND_010)=IND_020; p_zder(IND_010)=IND_011
  p_xder(IND_001)=IND_101; p_yder(IND_001)=IND_011; p_zder(IND_001)=IND_002
  p_xder(IND_200)=IND_300; p_yder(IND_200)=IND_210; p_zder(IND_200)=IND_201
  p_xder(IND_020)=IND_120; p_yder(IND_020)=IND_030; p_zder(IND_020)=IND_021
  p_xder(IND_002)=IND_102; p_yder(IND_002)=IND_012; p_zder(IND_002)=IND_003
  p_xder(IND_110)=IND_210; p_yder(IND_110)=IND_120; p_zder(IND_110)=IND_111
  p_xder(IND_101)=IND_201; p_yder(IND_101)=IND_111; p_zder(IND_101)=IND_102
  p_xder(IND_011)=IND_111; p_yder(IND_011)=IND_021; p_zder(IND_011)=IND_012
  p_xder(IND_300)=IND_400; p_yder(IND_300)=IND_310; p_zder(IND_300)=IND_301
  p_xder(IND_030)=IND_130; p_yder(IND_030)=IND_040; p_zder(IND_030)=IND_031
  p_xder(IND_003)=IND_103; p_yder(IND_003)=IND_013; p_zder(IND_003)=IND_004
  p_xder(IND_210)=IND_310; p_yder(IND_210)=IND_220; p_zder(IND_210)=IND_211
  p_xder(IND_201)=IND_301; p_yder(IND_201)=IND_211; p_zder(IND_201)=IND_202
  p_xder(IND_120)=IND_220; p_yder(IND_120)=IND_130; p_zder(IND_120)=IND_121
  p_xder(IND_021)=IND_121; p_yder(IND_021)=IND_031; p_zder(IND_021)=IND_022
  p_xder(IND_102)=IND_202; p_yder(IND_102)=IND_112; p_zder(IND_102)=IND_103
  p_xder(IND_012)=IND_112; p_yder(IND_012)=IND_022; p_zder(IND_012)=IND_013
  p_xder(IND_111)=IND_211; p_yder(IND_111)=IND_121; p_zder(IND_111)=IND_112
  p_xder(IND_400)=IND_500; p_yder(IND_400)=IND_410; p_zder(IND_400)=IND_401
  p_xder(IND_040)=IND_140; p_yder(IND_040)=IND_050; p_zder(IND_040)=IND_041
  p_xder(IND_004)=IND_104; p_yder(IND_004)=IND_014; p_zder(IND_004)=IND_005
  p_xder(IND_310)=IND_410; p_yder(IND_310)=IND_320; p_zder(IND_310)=IND_311
  p_xder(IND_301)=IND_401; p_yder(IND_301)=IND_311; p_zder(IND_301)=IND_302
  p_xder(IND_130)=IND_230; p_yder(IND_130)=IND_140; p_zder(IND_130)=IND_131
  p_xder(IND_031)=IND_131; p_yder(IND_031)=IND_041; p_zder(IND_031)=IND_032
  p_xder(IND_103)=IND_203; p_yder(IND_103)=IND_113; p_zder(IND_103)=IND_104
  p_xder(IND_013)=IND_113; p_yder(IND_013)=IND_023; p_zder(IND_013)=IND_014
  p_xder(IND_220)=IND_320; p_yder(IND_220)=IND_230; p_zder(IND_220)=IND_221
  p_xder(IND_202)=IND_302; p_yder(IND_202)=IND_212; p_zder(IND_202)=IND_203
  p_xder(IND_022)=IND_122; p_yder(IND_022)=IND_032; p_zder(IND_022)=IND_023
  p_xder(IND_211)=IND_311; p_yder(IND_211)=IND_221; p_zder(IND_211)=IND_212
  p_xder(IND_121)=IND_221; p_yder(IND_121)=IND_131; p_zder(IND_121)=IND_122
  p_xder(IND_112)=IND_212; p_yder(IND_112)=IND_122; p_zder(IND_112)=IND_113
end subroutine FH_MCMUR_DAV_fill_grad_ptrs
!----------------------------------------------------------

!-----------------------------------------------------------------
!-----------------------------------------------------------------
!-----------------------------------------------------------------
! EVALUATION SUBROUTINES
!-----------------------------------------------------------------
!-----------------------------------------------------------------

!-----------------------------------------------------------------
subroutine FH_MCMUR_DAV_recur(numlist,toplevel, &
                           mcmur_dav_init_recur_rule, &
                           expon1,expon2,R, &
                           delx,dely,delz,delr2,delr2inv)

   implicit none

   integer,intent(in) :: numlist,toplevel,mcmur_dav_init_recur_rule
   double precision,intent(in) :: expon1,expon2(numlist), &
                                  delx(numlist),dely(numlist),delz(numlist), &
                                  delr2(numlist),delr2inv(numlist)
   double precision, intent(inout) ::  &
                    R(numlist,0:toplevel,0:toplevel,0:toplevel,0:toplevel)
   include "interact_type.fh"
   include "mpole_sizes.fh"
   include 'scale.fh'

   double precision alpha,pi,x,y,two_over_root_pi
   integer j,t,u,v,n
   double precision term(NCACHE),term2(NCACHE), &
                    B(NCACHE,0:MAXRLEV),C(NCACHE), &
                    B2(NCACHE,0:MAXRLEV),C2(NCACHE)

   pi = 3.14159265358979323846d0
   two_over_root_pi = 2.d0 / sqrt(pi)

   ! first fill R(j,0,0,0,n) using the rule
   ! R(j,0,0,0,n+1) = (1/r*d/dr) R(j,0,0,0,n)

   if ( mcmur_dav_init_recur_rule == HERM_OVERLAP_HERM )then
      do j = 1,numlist
         alpha = expon1*expon2(j) / (expon1+expon2(j))
         B(j,0) = (alpha/pi)*sqrt(alpha/pi)*exp(-alpha*delr2(j))
         term(j) = -2.d0*alpha
      enddo
      do n = 1,toplevel
         ! rule here is R(j,0,0,0,n+1) = 1/r*d/dr[ R(j,0,0,0,n) ]
         do j = 1,numlist
            B(j,n) = term(j)*B(j,n-1)
         enddo
      enddo
      do n = 0,toplevel
         do j = 1,numlist
            R(j,0,0,0,n) = B(j,n)
         enddo
      enddo
   elseif ( mcmur_dav_init_recur_rule == HERM_COULOMB_HERM )then
      do j = 1,numlist
         alpha = expon1*expon2(j) / (expon1+expon2(j))
         x = alpha*delr2(j)
         if ( delr2(j) > small_value )then
            call FH_erfcfun(sqrt(x),y)!y = erfc(sqrt(x));1-y = erf(sqrt(x))
            B(j,0) = (1.d0-y)*sqrt(delr2(j))*delr2inv(j)
            C(j) = two_over_root_pi*sqrt(alpha)*exp(-x)
            term(j) = -2.d0*alpha
         else ! assume delr2(j) = 0
            term(j) = -2.d0*alpha
            C(j) = two_over_root_pi*sqrt(alpha)
            B(j,0) = C(j)
         endif
      enddo
      do n = 1,toplevel
         ! rule here is B(j,n) = 1/r*d/dr B(j,n-1)
         do j = 1,numlist
            if ( delr2(j) > small_value )then
               B(j,n) = (-(2.d0*dble(n)-1.d0)*B(j,n-1)+C(j))*delr2inv(j)
               C(j) = term(j)*C(j)
            else ! assume delr2(j) = 0
               C(j) = term(j)*C(j)
               B(j,n) = C(j) / (2.d0*dble(n) + 1.d0)
            endif
         enddo
      enddo
      do n = 0,toplevel
         do j = 1,numlist
            R(j,0,0,0,n) = B(j,n)
         enddo
      enddo
   elseif ( mcmur_dav_init_recur_rule == H_COULOMB_H_MINUS_MP_COULOMB_MP )then
      
      do j = 1,numlist
         alpha = expon1*expon2(j) / (expon1+expon2(j))
         x = alpha*delr2(j)
         call FH_erfcfun(sqrt(x),y) !y = erfc(sqrt(x));1-y=erf(sqrt(x))
         B(j,0) = (1.d0-y)*sqrt(delr2(j))*delr2inv(j)
         C(j) = two_over_root_pi*sqrt(alpha)*exp(-x)
         term(j) = -2.d0*alpha
     
         B2(j,0) = sqrt(delr2(j))*delr2inv(j)
      enddo
      do n = 1,toplevel
         ! rule here is B(j,n+1) = 1/r*d/dr B(j,n)
         do j = 1,numlist
            B(j,n) = (-(2.d0*n-1.d0)*B(j,n-1)+C(j))*delr2inv(j)
            C(j) = term(j)*C(j)
            B2(j,n) = (-(2.d0*n-1.d0)*B2(j,n-1))*delr2inv(j)
         enddo
      enddo
      do n = 0,toplevel
         do j = 1,numlist
            R(j,0,0,0,n) = B(j,n) - B2(j,n)
         enddo
      enddo
   elseif ( mcmur_dav_init_recur_rule == MPOLE_COULOMB_MPOLE )then
      do j = 1,numlist
         B(j,0) = sqrt(delr2(j))*delr2inv(j)
      enddo
      do n = 1,toplevel
         ! rule here is B(j,n+1) = 1/r*d/dr B(j,n)
         do j = 1,numlist
            B(j,n) = (-(2.d0*n-1.d0)*B(j,n-1))*delr2inv(j)
         enddo
      enddo
      do n = 0,toplevel
         do j = 1,numlist
            R(j,0,0,0,n) = B(j,n)
         enddo
      enddo
   elseif ( mcmur_dav_init_recur_rule == MPOLE_COULOMB_MPOLE_MINUS_HERM )then
      do j = 1,numlist
         B(j,0) = sqrt(delr2(j))*delr2inv(j)
         alpha = expon2(j)  ! multipole is hermite with infinite exponent
         x = alpha*delr2(j)
         call FH_erfcfun(sqrt(x),y) !y = erfc(sqrt(x));1-y=erf(sqrt(x))
         B2(j,0) = (1.d0-y)*sqrt(delr2(j))*delr2inv(j)
         C2(j) = two_over_root_pi*sqrt(alpha)*exp(-x)
         term2(j) = -2.d0*alpha
      enddo
      do n = 1,toplevel
         do j = 1,numlist
            B(j,n) = (-(2.d0*n-1.d0)*B(j,n-1))*delr2inv(j)
            B2(j,n) = (-(2.d0*n-1.d0)*B2(j,n-1)+C2(j))*delr2inv(j)
            C2(j) = term2(j)*C2(j)
         enddo
      enddo
      do n = 0,toplevel
         do j = 1,numlist
            R(j,0,0,0,n) = B(j,n) - B2(j,n)
         enddo
      enddo
   elseif ( mcmur_dav_init_recur_rule == MPOLE_COULOMB_HERM )then
      do j = 1,numlist
         alpha = expon2(j)  ! multipole is hermite with infinite exponent
         x = alpha*delr2(j)
         if ( delr2(j) > small_value )then
            call FH_erfcfun(sqrt(x),y)!y = erfc(sqrt(x));1-y = erf(sqrt(x))
            B(j,0) = (1.d0-y)*sqrt(delr2(j))*delr2inv(j)
            C(j) = two_over_root_pi*sqrt(alpha)*exp(-x)
            term(j) = -2.d0*alpha
         else ! assume delr2(j) = 0
            term(j) = -2.d0*alpha
            C(j) = two_over_root_pi*sqrt(alpha)
            B(j,0) = C(j)
         endif
      enddo
      do n = 1,toplevel
         do j = 1,numlist
            if ( delr2(j) > small_value )then
               B(j,n) = (-(2.d0*dble(n)-1.d0)*B(j,n-1)+C(j))*delr2inv(j)
               C(j) = term(j)*C(j)
            else ! assume delr2(j) = 0
               C(j) = term(j)*C(j)
               B(j,n) = C(j) / (2.d0*dble(n) + 1.d0)
            endif
         enddo
      enddo
      do n = 0,toplevel
         do j = 1,numlist
            R(j,0,0,0,n) = B(j,n)
         enddo
      enddo
   elseif ( mcmur_dav_init_recur_rule == HERM_COULOMB_MPOLE )then
      do j = 1,numlist
         alpha = expon1  ! multipole is hermite with infinite exponent
         x = alpha*delr2(j)
         if ( delr2(j) > small_value )then
            call FH_erfcfun(sqrt(x),y) !y = erfc(sqrt(x));1-y=erf(sqrt(x))
            B(j,0) = (1.d0-y)*sqrt(delr2(j))*delr2inv(j)
            C(j) = two_over_root_pi*sqrt(alpha)*exp(-x)
            term(j) = -2.d0*alpha
         else ! assume delr2(j) = 0
            term(j) = -2.d0*alpha
            C(j) = two_over_root_pi*sqrt(alpha)
            B(j,0) = C(j)
         endif
      enddo
      do n = 1,toplevel
         do j = 1,numlist
            if ( delr2(j) > small_value )then
               B(j,n) = (-(2.d0*dble(n)-1.d0)*B(j,n-1)+C(j))*delr2inv(j)
               C(j) = term(j)*C(j)
            else ! assume delr2(j) = 0
               C(j) = term(j)*C(j)
               B(j,n) = C(j) / (2.d0*dble(n) + 1.d0)
            endif
         enddo
      enddo
      do n = 0,toplevel
         do j = 1,numlist
            R(j,0,0,0,n) = B(j,n)
         enddo
      enddo
   elseif ( mcmur_dav_init_recur_rule == MPOLE_COULOMB_HERM_MINUS_MPOLE )then
      do j = 1,numlist
         alpha = expon2(j)  ! multipole is hermite with infinite exponent
         x = alpha*delr2(j)
         call FH_erfcfun(sqrt(x),y) !y = erfc(sqrt(x));1-y=erf(sqrt(x))
         B(j,0) = (1.d0-y)*sqrt(delr2(j))*delr2inv(j)
         C(j) = two_over_root_pi*sqrt(alpha)*exp(-x)
         term(j) = -2.d0*alpha
         B2(j,0) = sqrt(delr2(j))*delr2inv(j)
      enddo
      do n = 1,toplevel
         do j = 1,numlist
            B(j,n) = (-(2.d0*n-1.d0)*B(j,n-1)+C(j))*delr2inv(j)
            C(j) = term(j)*C(j)
            B2(j,n) = (-(2.d0*n-1.d0)*B2(j,n-1))*delr2inv(j)
         enddo
      enddo
      do n = 0,toplevel
         do j = 1,numlist
            R(j,0,0,0,n) = B(j,n) - B2(j,n)
         enddo
      enddo
   elseif ( mcmur_dav_init_recur_rule == HERM_MINUS_MPOLE_COULOMB_MPOLE )then
      do j = 1,numlist
         alpha = expon1  ! multipole is hermite with infinite exponent
         x = alpha*delr2(j)
         call FH_erfcfun(sqrt(x),y) !y = erfc(sqrt(x));1-y=erf(sqrt(x))
         B(j,0) = (1.d0-y)*sqrt(delr2(j))*delr2inv(j)
         C(j) = two_over_root_pi*sqrt(alpha)*exp(-x)
         term(j) = -2.d0*alpha
         B2(j,0) = sqrt(delr2(j))*delr2inv(j)
      enddo
      do n = 1,toplevel
         do j = 1,numlist
            B(j,n) = (-(2.d0*n-1.d0)*B(j,n-1)+C(j))*delr2inv(j)
            C(j) = term(j)*C(j)
            B2(j,n) = (-(2.d0*n-1.d0)*B2(j,n-1))*delr2inv(j)
         enddo
      enddo
      do n = 0,toplevel
         do j = 1,numlist
            R(j,0,0,0,n) = B(j,n) - B2(j,n)
         enddo
      enddo
   endif
   ! next use standard recursion for all types of mcmur_dav_init_recur_rules
   v = 1
   do n = 0,toplevel - v
      do j = 1,numlist
         R(j,0,0,v,n) = delz(j)*R(j,0,0,0,n+1)
      enddo
   enddo
   do v = 2,toplevel
      do n = 0,toplevel - v
         do j = 1,numlist
            R(j,0,0,v,n) = delz(j)*R(j,0,0,v-1,n+1)+(v-1)*R(j,0,0,v-2,n+1)
         enddo
      enddo
   enddo
   u = 1
   do v = 0,toplevel - u
      do n = 0,toplevel - (u+v)
         do j = 1,numlist
            R(j,0,u,v,n) = dely(j)*R(j,0,0,v,n+1)
         enddo
      enddo
   enddo
   do u = 2,toplevel
      do v = 0,toplevel - u
         do n = 0,toplevel - (u+v)
            do j = 1,numlist
               R(j,0,u,v,n) = dely(j)*R(j,0,u-1,v,n+1)+(u-1)*R(j,0,u-2,v,n+1)
            enddo
         enddo
      enddo
   enddo
   t = 1
   do u = 0,toplevel - t
      do v = 0,toplevel - (t+u)
         do n = 0,toplevel - (t+u+v)
            do j = 1,numlist
               R(j,t,u,v,n) = delx(j)*R(j,0,u,v,n+1)
            enddo
         enddo
      enddo
   enddo
   do t = 2,toplevel
      do u = 0,toplevel - t
         do v = 0,toplevel - (t+u)
            do n = 0,toplevel - (t+u+v)
               do j = 1,numlist
                  R(j,t,u,v,n) = delx(j)*R(j,t-1,u,v,n+1)+(t-1)*R(j,t-2,u,v,n+1)
               enddo
            enddo
         enddo
      enddo
   enddo
end subroutine FH_MCMUR_DAV_recur
!-----------------------------------------------------
subroutine FH_MCMUR_DAV_fill_fields(numlist,dim1,toplevel, &
                            field_order1,hermite_order1,hermite_order2, &
                            field1,field2,hermite_coeff1,hermite_coeff2,R)

   implicit none

   integer,intent(in) ::  numlist,toplevel,dim1,field_order1, &
                          hermite_order1,hermite_order2
   double precision,intent(inout) :: field1(dim1,field_order1), &
                                     field2(dim1,hermite_order2)
   double precision, intent(in) ::  hermite_coeff1(hermite_order1), &
                                    hermite_coeff2(dim1,hermite_order2)
   double precision, intent(in) :: &
                R(numlist,0:toplevel,0:toplevel,0:toplevel,0:toplevel)
   include "direct_pointers.fh"

   integer j,fld,mp,tau,mu,nu,t,u,v,t_tau,u_mu,v_nu,mult,signum(0:9)
   double precision term

   signum(0) = 1
   do j = 1,9
      signum(j) = -signum(j-1)
   enddo
   ! clear fields
   do fld = 1,field_order1
      do j = 1,numlist
         field1(j,fld) = 0.d0
      enddo
   enddo
   do fld = 1,hermite_order2
      do j = 1,numlist
         field2(j,fld) = 0.d0
      enddo
   enddo

   do fld = 1,field_order1
      tau = t_part(fld)
      mu = u_part(fld)
      nu = v_part(fld)
      mult = signum(tau+mu+nu)
      do mp = 1,hermite_order2
         t = t_part(mp)
         u = u_part(mp)
         v = v_part(mp)
         t_tau = t + tau
         u_mu = u + mu
         v_nu = v + nu
         do j = 1,numlist
            field1(j,fld) = field1(j,fld) + &
                  mult*hermite_coeff2(j,mp)*R(j,t_tau,u_mu,v_nu,0) 
         enddo
      enddo
   enddo
   do fld = 1,hermite_order2
      t = t_part(fld)
      u = u_part(fld)
      v = v_part(fld)
      do mp = 1,hermite_order1
         tau = t_part(mp)
         mu = u_part(mp)
         nu = v_part(mp)
         t_tau = t + tau
         u_mu = u + mu
         v_nu = v + nu
         mult = signum(tau+mu+nu)
         term = hermite_coeff1(mp)
         do j = 1,numlist
            field2(j,fld) = field2(j,fld) + &
                  mult*term*R(j,t_tau,u_mu,v_nu,0) 
         enddo
      enddo
   enddo
      
end subroutine FH_MCMUR_DAV_fill_fields
!----------------------------------------------------------
subroutine FH_MCMUR_DAV_update_Global_fields( &
                           numlist,dim1,signum,factor,self_flag, &
                           order1,offset1,order2,field_offset2, &
                           Global_field1,Global_field2,field1,field2)

   implicit none 

   integer,intent(in) :: numlist,dim1,signum,self_flag, &
                         order1,offset1,order2
   double precision,intent(in) :: factor
   integer,intent(in) :: field_offset2(*)
   double precision,intent(inout) :: Global_field1(*),Global_field2(*)
   double precision,intent(in) :: field1(dim1,*),field2(dim1,*)

   integer fld,j,offset2

   do fld = 1,order1
      do j = 1,numlist
         Global_field1(offset1+fld) = Global_field1(offset1+fld) +  &
                                    signum*factor*field1(j,fld) 
      enddo
   enddo
   if ( self_flag == 1 )return !field2 is double counting
   ! now for field2. 
   do j = 1,numlist
      offset2 = field_offset2(j)
      do fld = 1,order2
         Global_field2(offset2+fld) = Global_field2(offset2+fld) +  &
                                    signum*factor*field2(j,fld) 
      enddo
   enddo
end subroutine FH_MCMUR_DAV_update_Global_fields
!------------------------------------------------------------------
subroutine FH_MCMUR_DAV_ene_frc(numlist,dim1,signum,factor,order1,site1, &
                             site_pointer,delx,dely,delz, &
                             hermite_coeff1,field1,energy,site_frc,virial)

   implicit none

   integer,intent(in) :: numlist,dim1,signum,order1,site1
   double precision,intent(in) :: factor
   integer,intent(in) :: site_pointer(*)
   double precision,intent(in) :: hermite_coeff1(*)
   double precision,intent(in) :: field1(dim1,*)
   double precision, intent(in)    :: delx(numlist)
   double precision, intent(in)    :: dely(numlist)
   double precision, intent(in)    :: delz(numlist)
   double precision,intent(inout) :: energy,site_frc(3,*)
   double precision, intent(inout) :: virial(3,3)
   
   include "direct_pointers.fh"

   integer j,mp,kx,ky,kz,site2
   double precision term,dfx,dfy,dfz,sumx,sumy,sumz,ene
   double precision      :: vxx, vxy, vxz, vyx, vyy, vyz, vzx, vzy, vzz
   sumx = 0.d0
   sumy = 0.d0
   sumz = 0.d0
   ene = 0.d0
   vxx = 0.d0
   vxy = 0.d0
   vxz = 0.d0
   vyx = 0.d0
   vyy = 0.d0
   vyz = 0.d0
   vzx = 0.d0
   vzy = 0.d0
   vzz = 0.d0

   do mp = 1,order1
      kx = p_xder(mp)
      ky = p_yder(mp)
      kz = p_zder(mp)
      term = signum*factor*hermite_coeff1(mp)
      do j = 1,numlist
         ene = ene + term*field1(j,mp)
         dfx = term*field1(j,kx)
         dfy = term*field1(j,ky)
         dfz = term*field1(j,kz)
         sumx = sumx + dfx
         sumy = sumy + dfy
         sumz = sumz + dfz
         site2 = site_pointer(j)
         site_frc(1,site2) = site_frc(1,site2) + dfx
         site_frc(2,site2) = site_frc(2,site2) + dfy
         site_frc(3,site2) = site_frc(3,site2) + dfz
         vxx = vxx - delx(site2) * dfx
         vxy = vxy - delx(site2) * dfy
         vxz = vxz - delx(site2) * dfz
         vyx = vyx - dely(site2) * dfx
         vyy = vyy - dely(site2) * dfy
         vyz = vyz - dely(site2) * dfz
         vzx = vzx - delz(site2) * dfx
         vzy = vzy - delz(site2) * dfy
         vzz = vzz - delz(site2) * dfz
      enddo
   enddo
   energy = energy + ene
   site_frc(1,site1) = site_frc(1,site1) - sumx
   site_frc(2,site1) = site_frc(2,site1) - sumy
   site_frc(3,site1) = site_frc(3,site1) - sumz
   virial(1, 1) = virial(1, 1) + vxx
   virial(1, 2) = virial(1, 2) + 0.5d0 * (vxy + vyx)
   virial(1, 3) = virial(1, 3) + 0.5d0 * (vxz + vzx)
   virial(2, 1) = virial(2, 1) + 0.5d0 * (vxy + vyx)
   virial(2, 2) = virial(2, 2) + vyy
   virial(2, 3) = virial(2, 3) + 0.5d0 * (vyz + vzy)
   virial(3, 1) = virial(3, 1) + 0.5d0 * (vxz + vzx)
   virial(3, 2) = virial(3, 2) + 0.5d0 * (vyz + vzy)
   virial(3, 3) = virial(3, 3) + vzz
   
end subroutine FH_MCMUR_DAV_ene_frc
!------------------------------------------------------------------
