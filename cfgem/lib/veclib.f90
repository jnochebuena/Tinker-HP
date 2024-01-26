subroutine vdexp( n, x, y )

   implicit none

   integer  n
   double precision :: x(n), y(n)
   
      y(1:n) = exp(x(1:n))
end subroutine vdexp
