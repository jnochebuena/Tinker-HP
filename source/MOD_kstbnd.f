c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  module kstbnd  --  forcefield parameters for stretch-bend  ##
c     ##                                                             ##
c     #################################################################
c
c
c     maxnsb   maximum number of stretch-bend parameter entries
c
c     stbn     force constant parameters for stretch-bend terms
c     ksb      string of atom classes for stretch-bend terms
c
c
      module kstbnd
      implicit none
      integer maxnsb
      parameter (maxnsb=2000)
      integer*8 ksb(maxnsb)
      integer*8 ksb_sys(0:maxnsb)
      real*8 stbn(2,maxnsb)
c      character*12 ksb(maxnsb)
      save
      end
