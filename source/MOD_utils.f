c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###################################################################
c     ##                                                               ##
c     ##  module utils  --  General Stuff                              ##
c     ##                                                               ##
c     ###################################################################
c
c
      module utils
      implicit none
      interface front_convert_base
        subroutine front_convert_base3(n2,n1,n0,number)
          integer  ,intent(in) :: n2,n1,n0
          integer*8,intent(out):: number
        end subroutine
        subroutine front_convert_base5(n4,n3,n2,n1,n0,number)
          integer  ,intent(in) :: n4,n3,n2,n1,n0
          integer*8,intent(out):: number
        end subroutine
      end interface
      interface back_convert_base
        subroutine back_convert_base5(n4,n3,n2,n1,n0,number)
          integer  ,intent(out):: n4,n3,n2,n1,n0
          integer*8,intent(in) :: number
        end subroutine
      end interface
      interface is_find
        module procedure is_find8
        module procedure is_find4
      end interface

      contains
      function is_find8(array,n,number) 
      implicit none
      integer*8,intent(in ):: array(*)
      integer*8,intent(in ):: number
      integer  ,intent(in ):: n
      logical  :: is_find8
      integer i 

      is_find8 = .false.
      do i = 1, n
         if (array(i).eq.number) then
            is_find8 = .true.
            return
         end if
      end do
      end
      function is_find4(array,n,number)
      implicit none
      integer  ,intent(in ):: array(*)
      integer  ,intent(in ):: number
      integer  ,intent(in ):: n
      logical  :: is_find4
      integer i 

      is_find4 = .false.
      do i = 1, n
         if (array(i).eq.number) then
            is_find4 = .true.
            return
         end if
      end do
      end

      end module
