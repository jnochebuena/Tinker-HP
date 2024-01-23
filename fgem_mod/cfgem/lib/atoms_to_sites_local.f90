!----------------------------------------------------------------------
subroutine ATMSITE_atom_to_sitecrds(numatoms,crd,sitecrd,site_index_of_atom)

  implicit none

  integer,intent(in) :: numatoms,site_index_of_atom(numatoms)
  double precision,intent(in) :: crd(3,*)
  double precision,intent(out) :: sitecrd(3,*)
  integer n,m
  if ( numatoms == 0 )return
  do n = 1,numatoms
    m = site_index_of_atom(n)  ! site this atom is mapped to
    sitecrd(1,m) = crd(1,n)
    sitecrd(2,m) = crd(2,n)
    sitecrd(3,m) = crd(3,n)
  enddo
  return
end subroutine ATMSITE_atom_to_sitecrds
!----------------------------------------------------------------------
subroutine ATMSITE_build_centers(numatoms,num_cenlist,cen_list,cen_wt,sitecrd)

  implicit none

  integer,intent(in) :: numatoms,num_cenlist
  integer,intent(in) :: cen_list(2,num_cenlist)
  double precision,intent(in) :: cen_wt(num_cenlist)
  double precision,intent(out) :: sitecrd(3,*)

  integer k,m,n
  double precision w

  if ( num_cenlist == 0 )return
  do n = 1,num_cenlist
    k = cen_list(1,n)  !contributing atom site
    m = cen_list(2,n)  !resulting extra center
    w = cen_wt(n)  ! weight of atom k contribution
    sitecrd(1,m) = sitecrd(1,m) + w*sitecrd(1,k)
    sitecrd(2,m) = sitecrd(2,m) + w*sitecrd(2,k)
    sitecrd(3,m) = sitecrd(3,m) + w*sitecrd(3,k)
  enddo
  return
end subroutine ATMSITE_build_centers
!----------------------------------------------------------------------
subroutine ATMSITE_build_frame_def_pts(   &
                 numfrdeflist,framedeflist,sitecrd, &
                 nframes,frdefpt,frame_def_pt_markflag)

   implicit none

   integer,intent(in) :: numfrdeflist,nframes
   integer,intent(in) :: framedeflist(5,numfrdeflist)
   double precision,intent(in) :: sitecrd(3,*)
   double precision,intent(out) :: frdefpt(3,2,nframes)
   integer,intent(out) :: frame_def_pt_markflag(nframes)

   integer n,i,j,k,l,m
   double precision dx,dy,dz,wt

  if ( numfrdeflist == 0 )return

   do n = 1,numfrdeflist
      i = framedeflist(1,n)   ! center i gives tail of vector
      j = framedeflist(2,n)   ! center j gives head of vector
      k = framedeflist(3,n)   ! the frame number this vector helps define
      l = framedeflist(4,n)   ! the frame defpoint (1 or 2) within frame
      m = framedeflist(5,n)   ! m is number of such (unit) vectors 
                                ! averaged to get the frame defpoint k
      dx = sitecrd(1,j) - sitecrd(1,i)
      dy = sitecrd(2,j) - sitecrd(2,i)
      dz = sitecrd(3,j) - sitecrd(3,i)
      wt = m*sqrt(dx*dx+dy*dy+dz*dz) ! divide by length of ij times num pairs
      frdefpt(1,l,k) = frdefpt(1,l,k) + dx/wt
      frdefpt(2,l,k) = frdefpt(2,l,k) + dy/wt
      frdefpt(3,l,k) = frdefpt(3,l,k) + dz/wt
      frame_def_pt_markflag(k) = 1
   enddo
end subroutine ATMSITE_build_frame_def_pts
!----------------------------------------------------------------------
subroutine ATMSITE_defpoints_to_frames(nframes,frame_axis,frdefpt,frame, &
                           frame_def_pt_markflag,frame_doneflag)

   implicit none

   integer,intent(in) :: nframes,frame_axis(3)
   double precision,intent(in) :: frdefpt(3,2,nframes)
   double precision,intent(out) :: frame(3,3,nframes)
   integer,intent(in) :: frame_def_pt_markflag(nframes)
   integer,intent(inout) :: frame_doneflag(nframes)
  

   integer i,n,k1,k2,k3
   double precision u(3),v(3),w(3),siz,dot

   k1 = frame_axis(1)
   k2 = frame_axis(2)
   k3 = frame_axis(3)
   do n = 1,nframes
      ! check if defpts are nonzero and frame has not previously been done
      if ( frame_def_pt_markflag(n) == 1 .and. frame_doneflag(n) == 0 )then
         ! u is unit vector in direction of primary def pt
         do i = 1,3
            u(i) = frdefpt(i,1,n)
         enddo
         siz = sqrt(u(1)*u(1)+u(2)*u(2)+u(3)*u(3))
         do i = 1,3
            u(i) = u(i) / siz
         enddo
         ! v is unit vector given by component of secondary pt orthog to u
         do i = 1,3
            v(i) = frdefpt(i,2,n)
         enddo
         dot = u(1)*v(1)+u(2)*v(2)+u(3)*v(3)
         do i = 1,3
            v(i) = v(i) - dot*u(i)
         enddo
         siz = sqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
         do i = 1,3
            v(i) = v(i) / siz
         enddo
         ! w is u cross v
         w(1) = u(2)*v(3) - u(3)*v(2)
         w(2) = u(3)*v(1) - u(1)*v(3)
         w(3) = u(1)*v(2) - u(2)*v(1)
         !  now define frame
         do i = 1,3
            frame(i,k1,n) = u(i)
            frame(i,k2,n) = v(i)
            frame(i,k3,n) = w(i)
         enddo
         frame_doneflag(n) = 1
      endif
  enddo
end subroutine ATMSITE_defpoints_to_frames
!----------------------------------------------------------------------
subroutine ATMSITE_build_extrasites(num_extra_site_list, &
                 extra_site_list,loc_vector,frame,sitecrd)

  implicit none

  integer,intent(in) ::  num_extra_site_list
  integer,intent(in) ::  extra_site_list(3,num_extra_site_list)
  double precision,intent(in) :: loc_vector(3,num_extra_site_list), &
                                 frame(3,3,*)
  double precision,intent(inout) :: sitecrd(3,*)

  integer k,l,m,n
  double precision u,v,w
  if ( num_extra_site_list == 0 )return

  do n = 1,num_extra_site_list
    k = extra_site_list(1,n)   !site number for this site
    l = extra_site_list(2,n)   !origin number for this site
    m = extra_site_list(3,n)   !frame number for local coords
    u = loc_vector(1,n)
    v = loc_vector(2,n)
    w = loc_vector(3,n)
    sitecrd(1,k) = sitecrd(1,l) + u*frame(1,1,m)+v*frame(1,2,m)+w*frame(1,3,m)
    sitecrd(2,k) = sitecrd(2,l) + u*frame(2,1,m)+v*frame(2,2,m)+w*frame(2,3,m)
    sitecrd(3,k) = sitecrd(3,l) + u*frame(3,1,m)+v*frame(3,2,m)+w*frame(3,3,m)
  enddo
end subroutine ATMSITE_build_extrasites
!-----------------------------------------------------------------------------
subroutine ATMSITE_accum_de_dframe_rot(  &
               nsites,nframes,indframe,frame_axis,    &
               de_drotsite,frame,          &
               frdefpt,                    &
               de_drotframe)
!----------------------------------------------------------------------

   implicit none

   integer,intent(in) :: nsites,nframes,indframe(nsites),frame_axis(3)
   double precision,intent(in) :: de_drotsite(3,nsites),frame(3,3,nframes)
   double precision,intent(in) :: frdefpt(3,2,nframes)
   double precision,intent(out) :: de_drotframe(3,nframes)

   integer n,j,k,k1,k2,k3
   double precision p2unit(3),siz
   do n = 1,nframes   !clear the derivs wrt frame rotation
     do j = 1,3
        de_drotframe(j,n) = 0.d0
     enddo
   enddo
   k1 = frame_axis(1)
   k2 = frame_axis(2)
   k3 = frame_axis(3)
   ! deriv of energy with respect to rotation about unit vectors along
   ! p1 p2 and their cross product
   ! note that unit vector along p1 corresponds to k1st frame axis
   ! and unit vector in p1 x p2 direction corresponds to k3rd frame axis
   ! the energy derivative with respect to rotation about any unit vector
   ! is for each mpole given by the dot product of de_drotpole 
   ! (which is negative of torque due to that mpole) with the unit vector
   ! a frame may define several mpoles; hence the need for  pointers
   do n = 1,nsites
      k = indframe(n)
      if ( k > 0 )then   !skip those with no frame pointers
         siz = sqrt(frdefpt(1,2,k)**2+frdefpt(2,2,k)**2+frdefpt(3,2,k)**2)
         do j = 1,3
            p2unit(j) = frdefpt(j,2,k) / siz
         enddo
         de_drotframe(k1,k) = de_drotframe(k1,k) +   &
                              de_drotsite(1,n)*frame(1,k1,k) +  &
                              de_drotsite(2,n)*frame(2,k1,k) +  &
                              de_drotsite(3,n)*frame(3,k1,k)
         de_drotframe(k2,k) = de_drotframe(k2,k) +   &
                              de_drotsite(1,n)*p2unit(1) +   &
                              de_drotsite(2,n)*p2unit(2) +   &
                              de_drotsite(3,n)*p2unit(3)
         de_drotframe(k3,k) = de_drotframe(k3,k) +    &
                              de_drotsite(1,n)*frame(1,k3,k) +  &
                              de_drotsite(2,n)*frame(2,k3,k) +  &
                              de_drotsite(3,n)*frame(3,k3,k)
      endif
   enddo
end subroutine ATMSITE_accum_de_dframe_rot
!----------------------------------------------------------------------
subroutine ATMSITE_accum_de_ddefpts(   &
           nframes,frame_axis,   &
           de_drotframe,frame,defpts,  &
           de_ddefpt)

   implicit none

   integer,intent(in) :: nframes,frame_axis(3)
   double precision,intent(in) :: de_drotframe(3,nframes),frame(3,3,nframes)
   double precision,intent(in) :: defpts(3,2,nframes)
   double precision,intent(out) :: de_ddefpt(3,2,nframes)
   ! get the derivs of energy with respect to movement of defpoints
! expressed in the local frame coord system
   integer n,k,j,k1,k2,k3
   double precision :: p1(3),p2(3),p2unit(3),p2perp1(3),p1perp2(3),  &
                       u(3),v(3),w(3),dotu,dotv,dotw,  &
                       sizp1perp2,sizp2perp1,dot12,dot21,  &
                       sizp1,sizp2,dedrotp1,dedrotp2,dedu,dedv,dedw, &
                       de_drotu,de_drotv,de_drotw

   k1 = frame_axis(1)
   k2 = frame_axis(2)
   k3 = frame_axis(3)
   do n = 1,nframes
      do j = 1,3
         p1(j) = defpts(j,1,n)
         p2(j) = defpts(j,2,n)
         u(j) = frame(j,k1,n)
         v(j) = frame(j,k2,n)
         w(j) = frame(j,k3,n)
      enddo
      de_drotu = de_drotframe(k1,n)
      de_drotv = de_drotframe(k2,n)
      de_drotw = de_drotframe(k3,n)

      sizp1 = sqrt( p1(1)**2 + p1(2)**2 + p1(3)**2 )
      sizp2 = sqrt( p2(1)**2 + p2(2)**2 + p2(3)**2 )
      do j = 1,3
         p2unit(j) = p2(j) / sizp2
         !p1unit(j) = u(j) so no need to recalculate
      enddo
      dot21 = u(1)*p2(1) + u(2)*p2(2) + u(3)*p2(3)
      dot12 = p1(1)*p2unit(1) + p1(2)*p2unit(2) + p1(3)*p2unit(3)
      do j = 1,3
         p2perp1(j) = p2(j) - dot21*u(j)
         p1perp2(j) = p1(j) - dot12*p2unit(j)
      enddo
      sizp2perp1 = sqrt(p2perp1(1)**2 + p2perp1(2)**2 + p2perp1(3)**2)
      sizp1perp2 = sqrt(p1perp2(1)**2 + p1perp2(2)**2 + p1perp2(3)**2)
      ! def point one is along axis one. movement du parallel to that axis does
      ! not rotate the frame..so deriv is zero
      !    dedu = 0.d0
      ! movement dv in v-axis direction corresponds to rotation about local 
      ! w-axis of dtheta = dv/sizp1; thus a change in energy of 
      !    dE = dedrotw*dtheta = dedrotw*dv/sizp1
      !    dE/dv = dedrotw /sizp1
      dedv = de_drotw/sizp1
      ! movement dw in w-axis direction corresponds to rotation about p2unit
      ! of dtheta = -dw/sizp1perp2 (clockwise rotation) 
      dedw = -de_drotv/sizp1perp2
      ! Now convert to derivs wrt x,y,z. u = p.u = x*u(1)+y*u(2)+z*u(3)
      ! thus dudx = u(1). similaryl dvdx = v(1)
      de_ddefpt(1,1,n) = dedv*v(1) + dedw*w(1)
      de_ddefpt(2,1,n) = dedv*v(2) + dedw*w(2)
      de_ddefpt(3,1,n) = dedv*v(3) + dedw*w(3)
 
      ! for point 2..any movement in the local uv plane does not affect the 
      ! frame
      !       dedu = 0.d0
      !       dedv = 0.d0
      ! movement dw in w direction corresponds to rotation about local u-axis 
      ! of dtheta = dw/sizp2perpu
      dedw = de_drotu/sizp2perp1
      de_ddefpt(1,2,n) = dedw*w(1)
      de_ddefpt(2,2,n) = dedw*w(2)
      de_ddefpt(3,2,n) = dedw*w(3)
   enddo      
end subroutine ATMSITE_accum_de_ddefpts
!----------------------------------------------------------------------
subroutine ATMSITE_de_ddefpts_to_sites(  &
                numfrdeflist,framedeflist,sitecrd,sitefrc,  &
                nframes,de_ddefpt,virial)

   implicit none

   integer,intent(in) :: numfrdeflist,nframes
   integer,intent(in) :: framedeflist(5,numfrdeflist)
   double precision,intent(in) :: sitecrd(3,*),de_ddefpt(3,2,nframes)
   double precision,intent(inout) :: sitefrc(3,*)
   double precision,intent(inout) :: virial(3,3)

   integer n,i,j,k,l,m,p
   double precision siz,siz2,dx,dy,dz,dux_dx,dux_dy,dux_dz, &
        duy_dy,duy_dz,duz_dz,dedux,deduy,deduz,dedx,dedy,dedz,  &
        siz3inv,sizinv
   double precision vxx, vxy, vxz, vyx, vyy, vyz, vzx, vzy, vzz

    vxx = 0.d0
    vxy = 0.d0
    vxz = 0.d0
    vyx = 0.d0
    vyy = 0.d0
    vyz = 0.d0
    vzx = 0.d0
    vzy = 0.d0
    vzz = 0.d0
   do n = 1,numfrdeflist
      i = framedeflist(1,n)   ! site i gives tail of vector
      j = framedeflist(2,n)   ! site j gives head of vector
      k = framedeflist(3,n)   ! the frame number this vector helps define
      l = framedeflist(4,n)   ! the frame defpoint (1 or 2) within frame
      m = framedeflist(5,n)   ! m is number of such (unit) vectors 
                            ! averaged to get the frame defpoint k
      dx = sitecrd(1,j) - sitecrd(1,i)
      dy = sitecrd(2,j) - sitecrd(2,i)
      dz = sitecrd(3,j) - sitecrd(3,i)
      siz2 = dx*dx+dy*dy+dz*dz
      siz = sqrt(siz2)
      siz3inv = 1.d0/(siz2*siz)
      sizinv = 1.d0/siz
      ! ux, uy, uz are given by dx/siz, dy/siz, and dz/siz 
      dux_dx = sizinv - dx*dx*siz3inv
      dux_dy = -dx*dy*siz3inv   ! note duy_dx = dux_dy use this below
      dux_dz = -dx*dz*siz3inv
      duy_dy = sizinv - dy*dy*siz3inv
      duy_dz = -dy*dz*siz3inv
      duz_dz = sizinv - dz*dz*siz3inv
      ! the derivs of E wrt coordinates of unit vector in ij direction are 
      ! given  by (1/m) times the derivs of E wrt coords of def point (l,k)
      ! since the def point is the simple average of m of these unit vectors
      dedux = de_ddefpt(1,l,k) / m
      deduy = de_ddefpt(2,l,k) / m
      deduz = de_ddefpt(3,l,k) / m
      ! now apply chain rule, using symmetry e.g. dux_dy = duy_dx
      dedx = dedux*dux_dx + deduy*dux_dy + deduz*dux_dz
      dedy = dedux*dux_dy + deduy*duy_dy + deduz*duy_dz
      dedz = dedux*dux_dz + deduy*duy_dz + deduz*duz_dz
      ! finally apply forces. note force is negative of deriv of energy 
      ! wrt position
      ! also note e.g. deriv of dx wrt x position of atoms i and j is -1,+1
      !print *,'Orig Forces'
      !print *,sitefrc(1,i),sitefrc(2,i),sitefrc(3,i)
      !print *,'Torques'
      !print *,dedx,dedy,dedz
      sitefrc(1,i) = sitefrc(1,i) + dedx
      sitefrc(2,i) = sitefrc(2,i) + dedy
      sitefrc(3,i) = sitefrc(3,i) + dedz
      sitefrc(1,j) = sitefrc(1,j) - dedx
      sitefrc(2,j) = sitefrc(2,j) - dedy
      sitefrc(3,j) = sitefrc(3,j) - dedz
      vxx = vxx + dedx * dx
      vxy = vxy + dedx * dy
      vxz = vxz + dedx * dz
      vyx = vyx + dedy * dx
      vyy = vyy + dedy * dy
      vyz = vyz + dedy * dz
      vzx = vzx + dedz * dx
      vzy = vzy + dedz * dy
      vzz = vzz + dedz * dz
   enddo

    virial(1,1) = virial(1,1) + vxx
    virial(1,2) = virial(1,2) + 0.5d0 * (vxy + vyx)
    virial(1,3) = virial(1,3) + 0.5d0 * (vxz + vzx)
    virial(2,1) = virial(2,1) + 0.5d0 * (vxy + vyx)
    virial(2,2) = virial(2,2) + vyy
    virial(2,3) = virial(2,3) + 0.5d0 * (vyz + vzy)
    virial(3,1) = virial(3,1) + 0.5d0 * (vxz + vzx)
    virial(3,2) = virial(3,2) + 0.5d0 * (vyz + vzy)
    virial(3,3) = virial(3,3) + vzz

end subroutine ATMSITE_de_ddefpts_to_sites
!-----------------------------------------------------------------------------
subroutine ATMSITE_add_extra_site_contrib(  &
                num_extra_site_list,extra_site_list,sitecrd, &
                sitefrc,de_drotsite)

   implicit none

   integer,intent(in) :: num_extra_site_list, &
                         extra_site_list(3,num_extra_site_list)      
   double precision,intent(in) :: sitecrd(3,*)
   double precision,intent(inout) :: sitefrc(3,*)
   double precision,intent(inout) :: de_drotsite(3,*)
!----------------------------------------------------------------------
! add the derivs of energy due to translation of extra sites in terms
! of frame rotations. The torque is given by r x F, so the term we are looking
! for is given by minus the torque
!----------------------------------------------------------------------
   integer :: k,l,n
   double precision :: r(3),torque(3)
   do n = 1,num_extra_site_list
      k = extra_site_list(1,n)   !site number for this site
      l = extra_site_list(2,n)   !origin number for this site
      ! add the translational contribution
      sitefrc(1,l) = sitefrc(1,l) + sitefrc(1,k)
      sitefrc(2,l) = sitefrc(2,l) + sitefrc(2,k)
      sitefrc(3,l) = sitefrc(3,l) + sitefrc(3,k)
      ! next the torque
      r(1) = sitecrd(1,k) - sitecrd(1,l)
      r(2) = sitecrd(2,k) - sitecrd(2,l)
      r(3) = sitecrd(3,k) - sitecrd(3,l)
      torque(1) = r(2)*sitefrc(3,k) - r(3)*sitefrc(2,k)
      torque(2) = r(3)*sitefrc(1,k) - r(1)*sitefrc(3,k)
      torque(3) = r(1)*sitefrc(2,k) - r(2)*sitefrc(1,k)
      de_drotsite(1,l) = de_drotsite(1,l) - torque(1)
      de_drotsite(2,l) = de_drotsite(2,l) - torque(2)
      de_drotsite(3,l) = de_drotsite(3,l) - torque(3)
   enddo

end subroutine ATMSITE_add_extra_site_contrib
!----------------------------------------------------------------------
subroutine ATMSITE_trans_de_dcens_to_atoms(numatoms,ncenlist,  &
                                       cen_list,cen_wt,sitefrc)

   implicit none

   integer,intent(in) :: numatoms,ncenlist,cen_list(2,*)
   double precision,intent(in) :: cen_wt(*)
   double precision,intent(inout) :: sitefrc(3,*)

   integer :: n,k,m
   double precision :: w
   do n = 1,ncenlist
      k = cen_list(1,n)  !contributing atom site
      m = cen_list(2,n)  !resulting extra center
      w = cen_wt(n)  ! weight of atom k contribution
      sitefrc(1,k) = sitefrc(1,k) + w * sitefrc(1,m)
      sitefrc(2,k) = sitefrc(2,k) + w * sitefrc(2,m)
      sitefrc(3,k) = sitefrc(3,k) + w * sitefrc(3,m)
   enddo

end subroutine ATMSITE_trans_de_dcens_to_atoms
!-----------------------------------------------------------
subroutine ATMSITE_sitefrc_to_frc(numatoms,sitefrc,frc,site_index_of_atom)

   implicit none

   integer,intent(in) :: numatoms,site_index_of_atom(numatoms)
   double precision,intent(in) :: sitefrc(3,*)
   double precision,intent(inout) :: frc(3,*)

   integer n,m
   do n = 1,numatoms
     m = site_index_of_atom(n)  ! site this atom is mapped to
     frc(1,n) = frc(1,n) + sitefrc(1,m)
     frc(2,n) = frc(2,n) + sitefrc(2,m)
     frc(3,n) = frc(3,n) + sitefrc(3,m)
   enddo
end subroutine ATMSITE_sitefrc_to_frc
!----------------------------------------------------------------------
