c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c
c     Convert forward a series of number to base n
c
      subroutine front_convert_base3 (n2,n1,n0,number)
      implicit none
      integer  ,intent(in) :: n2,n1,n0
      integer*8,intent(out):: number
      integer*8,parameter  :: base=1000

      if (n2.gt.base.or.n1.gt.base.or.n0.gt.base) then
         print*, 'ERROR in front_convert_base3 :: ',
     &           'call arguments are \n',n2,n1,n0 
      end if
      number = int(n0,8) + int(n1,8)*base + int(n2,8)*base**2
      end subroutine
c
      subroutine front_convert_base5 (n4,n3,n2,n1,n0,number)
      implicit none
      integer  ,intent(in) :: n4,n3,n2,n1,n0
      integer*8,intent(out):: number
      integer*8,parameter  :: base=1000

      if (n4.gt.base.or.n3.gt.base.or.n2.gt.base.or.
     &    n1.gt.base.or.n0.gt.base) then
         print*, 'ERROR in front_convert_base5 :: ',
     &           'call arguments are \n',n4,n3,n2,n1,n0 
      end if
      number = int(n0,8) + int(n1,8)*base + int(n2,8)*base**2 + 
     &         int(n3,8)*base**3 + int(n4,8)*base**4
      end subroutine
c
c     Decomposition of base-n number
c
      subroutine back_convert_base5 (n4,n3,n2,n1,n0,number)
      implicit none
      integer*8,intent(in) :: number
      integer  ,intent(out):: n4,n3,n2,n1,n0
      integer*8,parameter  :: base=1000
      integer*8 :: cnum

      cnum=number
      if (cnum.lt.0) then
         print*,"ERROR ! Can't convert back a negative number",
     &          number
      end if

      n0   = mod(cnum,base)
      cnum = cnum/base
      n1   = mod(cnum,base)
      cnum = cnum/base
      n2   = mod(cnum,base)
      cnum = cnum/base
      n3   = mod(cnum,base)
      cnum = cnum/base
      n4   = mod(cnum,base)
      cnum = cnum/base

      if (cnum.ne.0)
     &   print*,"ERROR ! convert unfinished",cnum
      end subroutine
c
