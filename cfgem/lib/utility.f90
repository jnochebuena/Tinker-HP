!------------------------------------------------------------------
! UTILITY routines
!----------------------------------------------------------------
subroutine UTIL_zero_real_array(array,num)

  implicit none

  double precision,intent(out) :: array(*)
  integer,intent(in) :: num

  integer i
  do i = 1,num
    array(i) = 0.d0
  end do
end subroutine UTIL_zero_real_array
!----------------------------------------------------------------
subroutine UTIL_zero_int_array(array,num)

  implicit none

  integer,intent(out) :: array(*)
  integer,intent(in) :: num

  integer i
  do i = 1,num
    array(i) = 0
  end do
end subroutine UTIL_zero_int_array
!----------------------------------------------------------
subroutine UTIL_copy_real_array(a,b,num)

   implicit none

   double precision, intent(in) :: a(*)
   double precision, intent(out) :: b(*)
  integer,intent(in) :: num

  integer i
  do i = 1,num
    b(i) = a(i)
  end do
end subroutine UTIL_copy_real_array
!----------------------------------------------------------
subroutine dot(v1,v2,result)

  implicit none

  double precision,intent(in) :: v1(3),v2(3)
  double precision,intent(out) :: result
  result = v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3)
  return
end subroutine dot
!--------------------------------------------------------
subroutine cross(v1,v2,v12)

  implicit none
!
!    v12 is cross product of v1 and v2
!
  double precision,intent(in) :: v1(3),v2(3)
  double precision,intent(out) :: v12(3)
  v12(1) = v1(2)*v2(3)-v1(3)*v2(2)
  v12(2) = v1(3)*v2(1)-v1(1)*v2(3)
  v12(3) = v1(1)*v2(2)-v1(2)*v2(1)
  return
end subroutine cross
!----------------------------------------------------------
subroutine FH_erfcinv(tol,erfcinv_tol)

   implicit none

   double precision,intent(in) :: tol
   double precision,intent(out) :: erfcinv_tol

   double precision xlo,xhi,x,delta,diff,y
   delta = 1.d-12
   xlo = 1.d0
   xhi = 30.d0
   diff = xhi-xlo
   do while( diff > delta )
      x = 0.5d0*(xhi+xlo)
      call FH_erfcfun(x,y)
      if ( y > tol )then
         xlo = x
      else
         xhi = x
      endif
      diff = xhi - xlo 
   enddo
   erfcinv_tol = x
end subroutine FH_erfcinv
!----------------------------------------------------------------
subroutine FH_erfcfun(x,erfc)

   implicit none

   double precision,intent(in) ::  x
   double precision,intent(out) :: erfc

   double precision absx, c, p, q, nonexperfc, erf

   absx=abs(x)
   if (x > 26.d0) then
      erfc = 0.d0

   else if (x < -5.5d0) then
      erfc = 2.0d0

   else if (absx <= 0.5d0) then
      c = x * x
      p=((-0.356098437018154d-1*c+0.699638348861914d1)*c+   &
            0.219792616182942d2)*c+0.242667955230532d3
      q=((c+0.150827976304078d2)*c+0.911649054045149d2)*c+  &
            0.215058875869861d3
      erf = x*p/q
      erfc = 1.d0-erf

   else if (absx < 4.d0) then
      c=absx
      p=((((((-0.136864857382717d-6*c+0.564195517478974d0)*c+    &
         0.721175825088309d1)*c+0.431622272220567d2)*c+    &
         0.152989285046940d3)*c+0.339320816734344d3)*c+    &
         0.451918953711873d3)*c+0.300459261020162d3
      q=((((((c+0.127827273196294d2)*c+0.770001529352295d2)*c+    &
         0.277585444743988d3)*c+0.638980264465631d3)*c+    &
         0.931354094850610d3)*c+0.790950925327898d3)*c+    &
         0.300459260956983d3
      if ( x > 0.d0 ) then
         nonexperfc = p/q
      else
         nonexperfc = 2.d0*exp(x*x) - p/q
      end if
      erfc = exp(-absx*absx)*nonexperfc
      if (x < 0.d0) erfc = 2.d0- erfc

   else
      c=1.d0/(x*x)
      p=(((0.223192459734185d-1*c+0.278661308609648d0)*c+      &
         0.226956593539687d0)*c+0.494730910623251d-1)*c+      &
         0.299610707703542d-2
      q=(((c+0.198733201817135d1)*c+0.105167510706793d1)*c+      &
         0.191308926107830d0)*c+0.106209230528468d-1
      c=(-c*p/q + 0.564189583547756d0)/absx
      if( x > 0.d0 ) then
         nonexperfc = c
      else
         nonexperfc = 2.d0*exp(x*x) - c
      end if
      erfc = exp(-absx*absx)*nonexperfc
      if (x < 0.d0) erfc = 2.d0- erfc
   end if
   return
end subroutine FH_erfcfun
!----------------------------------------------------------------
