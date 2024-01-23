c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine image  --  compute the minimum image distance  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "image" takes the components of pairwise distance between
c     two points in a periodic box and converts to the components
c     of the minimum image distance
c
c
c      subroutine image (xr,yr,zr)
c      use sizes
c      use boxes
c      use cell
c      implicit none
c      real*8 xr,yr,zr
cc
cc
cc     for orthogonal lattice, find the desired image directly
cc
c      if (orthogonal) then
c         do while (abs(xr) .gt. xcell2)
c            xr = xr - sign(xcell,xr)
c         end do
c         do while (abs(yr) .gt. ycell2)
c            yr = yr - sign(ycell,yr)
c         end do
c         do while (abs(zr) .gt. zcell2)
c            zr = zr - sign(zcell,zr)
c         end do
cc
cc     for monoclinic lattice, convert "xr" and "zr" to
cc     fractional coordinates, find desired image and then
cc     translate fractional coordinates back to Cartesian
cc
c      else if (monoclinic) then
c         zr = zr / beta_sin
c         xr = xr - zr*beta_cos
c         do while (abs(xr) .gt. xcell2)
c            xr = xr - sign(xcell,xr)
c         end do
c         do while (abs(yr) .gt. ycell2)
c            yr = yr - sign(ycell,yr)
c         end do
c         do while (abs(zr) .gt. zcell2)
c            zr = zr - sign(zcell,zr)
c         end do
c         xr = xr + zr*beta_cos
c         zr = zr * beta_sin
cc
cc     for triclinic lattice, convert pairwise components to
cc     fractional coordinates, find desired image and then
cc     translate fractional coordinates back to Cartesian
cc
c      else if (triclinic) then
c         zr = zr / gamma_term
c         yr = (yr - zr*beta_term) / gamma_sin
c         xr = xr - yr*gamma_cos - zr*beta_cos
c         do while (abs(xr) .gt. xcell2)
c            xr = xr - sign(xcell,xr)
c         end do
c         do while (abs(yr) .gt. ycell2)
c            yr = yr - sign(ycell,yr)
c         end do
c         do while (abs(zr) .gt. zcell2)
c            zr = zr - sign(zcell,zr)
c         end do
c         xr = xr + yr*gamma_cos + zr*beta_cos
c         yr = yr*gamma_sin + zr*beta_term
c         zr = zr * gamma_term
cc
cc     for truncated octahedron, use orthogonal box equations,
cc     then perform extra tests to remove corner pieces
cc
c      else if (octahedron) then
c         do while (abs(xr) .gt. xbox2)
c            xr = xr - sign(xbox,xr)
c         end do
c         do while (abs(yr) .gt. ybox2)
c            yr = yr - sign(ybox,yr)
c         end do
c         do while (abs(zr) .gt. zbox2)
c            zr = zr - sign(zbox,zr)
c         end do
c         if (abs(xr)+abs(yr)+abs(zr) .gt. box34) then
c            xr = xr - sign(xbox2,xr)
c            yr = yr - sign(ybox2,yr)
c            zr = zr - sign(zbox2,zr)
c         end if
c      end if
c      return
c      end

c     #################################################################
c     ##                                                             ##
c     ##  subroutine image2  --  compute the minimum image distance  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "image" takes the components of pairwise distance between
c     two points in a periodic box and converts to the components
c     of the minimum image distance
c
c
      subroutine image2 (xr,yr,zr)
      use sizes
      use boxes
      use cell
      implicit none
      real*8 xr,yr,zr,cel
c
c
c     for orthogonal lattice, find the desired image directly
c
      if (orthogonal) then
         if (abs(xr) .gt. xcell2) then
            cel      = sign(xcell,xr)
            xr = xr - cel*floor(abs(xr)/xcell)
            if ((abs(xr)) .gt. xcell2)
     &         xr = xr-sign(xcell,xr)
         end if
         if (abs(yr) .gt. xcell2) then
            cel      = sign(ycell,yr)
            yr = yr - cel*floor(abs(yr)/ycell)
            if ((abs(yr)) .gt. ycell2)
     &         yr = yr-sign(ycell,yr)
         end if
         if (abs(zr) .gt. zcell2) then
            cel      = sign(zcell,zr)
            zr = zr - cel*floor(abs(zr)/zcell)
            if ((abs(zr)) .gt. zcell2)
     &         zr = zr-sign(zcell,zr)
         end if
c         do while (abs(xr) .gt. xcell2)
c            xr = xr - sign(xcell,xr)
c         end do
c         do while (abs(yr) .gt. ycell2)
c            yr = yr - sign(ycell,yr)
c         end do
c         do while (abs(zr) .gt. zcell2)
c            zr = zr - sign(zcell,zr)
c         end do
c
c     for monoclinic lattice, convert "xr" and "zr" to
c     fractional coordinates, find desired image and then
c     translate fractional coordinates back to Cartesian
c
      else if (monoclinic) then
         zr = zr / beta_sin
         xr = xr - zr*beta_cos
         do while (abs(xr) .gt. xcell2)
            xr = xr - sign(xcell,xr)
         end do
         do while (abs(yr) .gt. ycell2)
            yr = yr - sign(ycell,yr)
         end do
         do while (abs(zr) .gt. zcell2)
            zr = zr - sign(zcell,zr)
         end do
         xr = xr + zr*beta_cos
         zr = zr * beta_sin
c
c     for triclinic lattice, convert pairwise components to
c     fractional coordinates, find desired image and then
c     translate fractional coordinates back to Cartesian
c
      else if (triclinic) then
         zr = zr / gamma_term
         yr = (yr - zr*beta_term) / gamma_sin
         xr = xr - yr*gamma_cos - zr*beta_cos
         do while (abs(xr) .gt. xcell2)
            xr = xr - sign(xcell,xr)
         end do
         do while (abs(yr) .gt. ycell2)
            yr = yr - sign(ycell,yr)
         end do
         do while (abs(zr) .gt. zcell2)
            zr = zr - sign(zcell,zr)
         end do
         xr = xr + yr*gamma_cos + zr*beta_cos
         yr = yr*gamma_sin + zr*beta_term
         zr = zr * gamma_term
c
c     for truncated octahedron, use orthogonal box equations,
c     then perform extra tests to remove corner pieces
c
      else if (octahedron) then
         do while (abs(xr) .gt. xbox2)
            xr = xr - sign(xbox,xr)
         end do
         do while (abs(yr) .gt. ybox2)
            yr = yr - sign(ybox,yr)
         end do
         do while (abs(zr) .gt. zbox2)
            zr = zr - sign(zbox,zr)
         end do
         if (abs(xr)+abs(yr)+abs(zr) .gt. box34) then
            xr = xr - sign(xbox2,xr)
            yr = yr - sign(ybox2,yr)
            zr = zr - sign(zbox2,zr)
         end if
      end if
      return
      end
c
      subroutine image(xr,yr,zr) 
      use cell
      implicit none
      real*8,intent(inout):: xr,yr,zr

      if (abs(xr) .gt. xcell2)
     &   xr  = xr - sign(xcell,xr)*floor((abs(xr)+xcell2)/xcell)
      if (abs(yr) .gt. ycell2)
     &   yr  = yr - sign(ycell,yr)*floor((abs(yr)+ycell2)/ycell)
      if (abs(zr) .gt. zcell2)
     &   zr  = zr - sign(zcell,zr)*floor((abs(zr)+zcell2)/zcell)
c
c!     Much slower on GPU
c     xr =  xr - int( (abs(xr) - xcell2) / xcell + 1.0_ti_p )
c    &         * sign (xcell,xr)
c     yr =  yr - int( (abs(yr) - ycell2) / ycell + 1.0_ti_p )
c    &         * sign (ycell,yr)
c     zr =  zr - int( (abs(zr) - zcell2) / zcell + 1.0_ti_p )
c    &         * sign (zcell,xr)
      end
