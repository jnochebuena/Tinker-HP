module sites

   implicit none

   private

   integer num_sites,num_frames,num_cenwt_list,num_reg_frame_list, & 
           num_extrapt_frame_list,num_extrapt_list,size_masked_sitelist
   integer,save, allocatable :: site_index_of_atom(:)
   double precision, save, allocatable :: site_crd(:,:),site_frc(:,:), &
                                          frac_site_crd(:,:), &
                                          map_site_crd(:,:), &
                                          frame_def_pt(:,:,:), &
                                          de_ddef_pt(:,:,:), &
                                          frame(:,:,:), &
                                          de_drotsite(:,:), &
                                          de_drotframe(:,:)
   integer,save, allocatable :: frame_def_pt_markflag(:),frame_doneflag(:)
   integer,save, allocatable :: aux_type_of_site(:),molnum_of_site(:), &
                                frame_index_of_site(:)

   integer,save, allocatable :: center_site_deflist(:,:)
   double precision, save, allocatable :: center_site_weight(:)
   integer, save, dimension(3) :: frame_axis_order
   integer, save, allocatable :: reg_frame_list(:,:)
   integer, save, allocatable :: extrapt_frame_list(:,:)
   integer, save, allocatable :: extra_pt_list(:,:)
   double precision, save, allocatable :: extra_pt_loc_crd(:,:)
   integer, save, allocatable :: num_sites_masked_by_site(:)
   integer, save, allocatable :: masked_sitelist_offset_of_site(:)
   integer, save, allocatable :: masked_sitelist(:)
   public FH_SITES_readfile,FH_SITES_deallocate,FH_SITES_build, &
          num_sites,aux_type_of_site,frame_index_of_site,frame, &
          de_drotsite,FH_SITES_de_drotsite_to_sitefrc, &
          FH_SITES_sitefrc_to_frc,FH_SITES_add_extrapt_torque, &
          num_sites_masked_by_site,masked_sitelist_offset_of_site, &
          masked_sitelist,site_crd,site_frc,frac_site_crd,map_site_crd
contains

subroutine FH_SITES_readfile(filename, free_lun, out_lun)

   use atoms,only : num_atoms
   use user, only : verbose

   implicit none

! Formal arguments:

   character(len=*), intent(in) :: filename
   integer,          intent(in) :: free_lun
   integer,          intent(in) :: out_lun

! Local variables:

   character(len=100)line
   integer j,n,natom,max_site,atom,site,itype,molnum,frameind,num,off
   integer ios,ier,ier1,ier2,ier3,ier4,ier5,ier6
   double precision weight,bohr,xloc,yloc,zloc

   open(unit=free_lun,file=filename,status='old')
   ! atom site index`
   read(free_lun,'(a)',iostat=ios)line
   if ( ios /= 0 .or. line(1:13) /= "atom_sitelist" )then
      write(out_lun,*)'FH_SITES_readfile: bad line: (atom_sitelist)',line(1:20)
      stop
   endif ! ( ios /= 0 )
   read(free_lun,'(a)',iostat=ios)line
   if ( ios /= 0 .or. line(1:9) /= "num_atoms" )then
      write(out_lun,*)'FH_SITES_readfile: bad line: (num_atoms)',line(1:20)
      stop
   endif ! ( ios /= 0 )
   read(line(10:),*)natom
   if ( natom /= num_atoms )then
      write(out_lun,*)'FH_SITES_readfile: mismatch in number of atoms'
      stop
   endif
   allocate(site_index_of_atom(num_atoms),stat=ier)
   if ( ier /= 0 )then
      write(out_lun,*)'FH_SITES_readfile: allocate fails!!'
      stop
   endif
   read(free_lun,'(a)',iostat=ios)line
   if ( ios /= 0 .or. line(1:11) /= "atom_number" )then
      write(out_lun,*)'FH_SITES_readfile: bad line: (atom_number)',line(1:20)
      stop
   endif ! ( ios /= 0 )
   max_site = -1
   do n = 1,num_atoms
      read(free_lun,*)atom,site
      if ( site > max_site )max_site = site
      if ( atom /= n )then
         write(out_lun,*)&
           'FH_SITES_readfile:atom mismatch in atom_site_index read'
         stop
      endif
      site_index_of_atom(n) = site
   enddo
   ! sitelist (type and molnum)
   read(free_lun,'(a)',iostat=ios)line
   if ( ios /= 0 .or. line(1:8) /= "sitelist" )then
      write(out_lun,*)'FH_SITES_readfile: bad line: (sitelist)',line(1:20)
      stop
   endif ! ( ios /= 0 )
   read(free_lun,'(a)',iostat=ios)line
   if ( ios /= 0 .or. line(1:8) /= "numsites" )then
      write(out_lun,*)'FH_SITES_readfile: bad line: (numsites)',line(1:20)
      stop
   endif ! ( ios /= 0 )
   read(line(9:),*)num_sites
   if ( max_site > num_sites )then
      write(out_lun,*)'num_sites too small for maximal atom_site_index'
      stop
   endif
   read(free_lun,'(a)',iostat=ios)line
   if ( ios /= 0 .or. line(1:9) /= "numframes" )then
      write(out_lun,*)'FH_SITES_readfile: bad line: (numframes)',line(1:20)
      stop
   endif ! ( ios /= 0 )
   read(line(10:),*)num_frames
   if ( verbose == 1 )then
      write(out_lun,*)'num_sites,frames = ',num_sites,num_frames
   endif

   ! following 4 needed later
   allocate(site_crd(3,num_sites),stat=ier1)
   allocate(frac_site_crd(3,num_sites),stat=ier2)
   allocate(map_site_crd(3,num_sites),stat=ier3)
   if ( ier1 /= 0 .or. ier2 /= 0 .or. ier3 /= 0 )then
      write(out_lun,*)'FH_SITES_readfile: allocate fails!'
      stop
   endif
   allocate(site_frc(3,num_sites),stat=ier1)
   allocate(frame_def_pt(3,2,num_frames),stat=ier2)
   allocate(frame(3,3,num_frames),stat=ier3)
   allocate(frame_def_pt_markflag(num_frames),stat=ier4)
   allocate(frame_doneflag(num_frames),stat=ier5)
   if ( ier1 /= 0 .or. ier2 /= 0 .or. ier3 /= 0 .or. ier4 /= 0 .or. &
        ier5 /= 0 )then
      write(out_lun,*)'FH_SITES_readfile: allocate fails!'
      stop
   endif
   ! these read in now
   allocate(aux_type_of_site(num_sites),stat=ier1)
   allocate(molnum_of_site(num_sites),stat=ier2)
   allocate(frame_index_of_site(num_sites),stat=ier3)
   allocate(de_drotsite(3,num_sites),stat=ier4)
   allocate(de_drotframe(3,num_frames),stat=ier5)
   allocate(de_ddef_pt(3,2,num_frames),stat=ier6)
   if ( ier1 /= 0 .or. ier2 /= 0 .or. ier3 /= 0 .or. ier4 /= 0 .or. &
        ier5 /= 0 .or. ier6 /= 0 )then
      write(out_lun,*)'FH_SITES_readfile: allocate fails!'
      stop
   endif
   read(free_lun,'(a)',iostat=ios)line
   if ( ios /= 0 .or. line(1:11) /= "site_number" )then
      write(out_lun,*)'FH_SITES_readfile: bad line: (site_number)',line(1:20)
      stop
   endif ! ( ios /= 0 )
   do n = 1,num_sites
      read(free_lun,*)site,itype,molnum,frameind
      if ( site /= n )then
         write(out_lun,*)'FH_SITES_readfile: reading sitelist; mismatch at #',n
         stop
      endif
      aux_type_of_site(n) = itype
      molnum_of_site(n) = molnum
      frame_index_of_site(n) = frameind
      if ( frameind > num_frames )then
         write(out_lun,*)'frame > num_frames: ',frame,num_frames
         stop
      endif
   enddo

   ! center list and weights
   read(free_lun,'(a)',iostat=ios)line
   if ( ios /= 0 .or. line(1:18) /= "center_weight_list" )then
      write(out_lun,*)&
        'FH_SITES_readfile: bad line: (center_weight_list)',line(1:20)
      stop
   endif ! ( ios /= 0 )
   read(free_lun,'(a)',iostat=ios)line
   if ( ios /= 0 .or. line(1:7) /= "numlist" )then
      write(out_lun,*)'FH_SITES_readfile: bad line: (numlist)',line(1:20)
      stop
   endif ! ( ios /= 0 )
   read(line(8:),*)num_cenwt_list
   if ( verbose == 1 )then
      write(out_lun,*)'num_cenwt_list = ',num_cenwt_list
   endif
   if ( num_cenwt_list > 0 )then
      allocate(center_site_deflist(2,num_cenwt_list),stat=ier1)
      allocate(center_site_weight(num_cenwt_list),stat=ier2)
      if ( ier1 /= 0 .or. ier2 /= 0 )then
         write(out_lun,*)'FH_SITES_readfile: allocate fails!'
         stop
      endif
      read(free_lun,'(a)',iostat=ios)line
      if ( ios /= 0 .or. line(1:8) /= "atom_num" )then
         write(out_lun,*)'FH_SITES_readfile: bad line: (atom_num)',line(1:20)
         stop
      endif ! ( ios /= 0 )
      do n = 1,num_cenwt_list
         read(free_lun,*)atom,site,weight
         center_site_deflist(1,n) = atom
         center_site_deflist(2,n) = site
         center_site_weight(n) = weight
      enddo
   endif
   ! frame axes
   read(free_lun,'(a)',iostat=ios)line
   if ( ios /= 0 .or. line(1:10) /= "frame_axes" )then
      write(out_lun,*)'FH_SITES_readfile: bad line: (frame_axes)',line(1:20)
      stop
   endif ! ( ios /= 0 )
   read(line(11:),*)frame_axis_order(1),frame_axis_order(2),frame_axis_order(3)
   ! regular framelist
   read(free_lun,'(a)',iostat=ios)line
   if ( ios /= 0 .or. line(1:17) /= "regular_framelist" )then
      write(out_lun,*)&
        'FH_SITES_readfile: bad line: (regular_framelist)',line(1:20)
      stop
   endif ! ( ios /= 0 )
   read(free_lun,'(a)',iostat=ios)line
   if ( ios /= 0 .or. line(1:13) /= "num_framelist" )then
      write(out_lun,*)'FH_SITES_readfile: bad line: (num_framelist)',line(1:20)
      stop
   endif ! ( ios /= 0 )
   read(line(14:),*)num_reg_frame_list
   if ( verbose == 1 )then
      write(out_lun,*)'num_reg_frame_list = ',num_reg_frame_list
   endif
   if ( num_reg_frame_list > 0 )then
      allocate(reg_frame_list(5,num_reg_frame_list),stat=ier)
      if ( ier /= 0 )then
         write(out_lun,*)'FH_SITES_readfile: allocate fails!'
         stop
      endif
      read(free_lun,'(a)',iostat=ios)line
      if ( ios /= 0 .or. line(1:4) /= "tail" )then
         write(out_lun,*)'FH_SITES_readfile: bad line: (tail)',line(1:20)
         stop
      endif ! ( ios /= 0 )
      do n = 1,num_reg_frame_list
         read(free_lun,*)(reg_frame_list(j,n),j=1,5)
      enddo
   endif
   ! extra_point framelist
   ! this is separate since the extrapts are defined in terms of reg frames
   ! hence chicken egg problem---first build reg frames, then extrapts
   ! finall the extrapt frames
   read(free_lun,'(a)',iostat=ios)line
   if ( ios /= 0 .or. line(1:17) /= "extrapt_framelist" )then
      write(out_lun,*)&
        'FH_SITES_readfile: bad line: (extrapt_framelist)',line(1:20)
      stop
   endif ! ( ios /= 0 )
   read(free_lun,'(a)',iostat=ios)line
   if ( ios /= 0 .or. line(1:13) /= "num_framelist" )then
      write(out_lun,*)'FH_SITES_readfile: bad line: (num_framelist)',line(1:20)
      stop
   endif ! ( ios /= 0 )
   read(line(14:),*)num_extrapt_frame_list
   if ( verbose == 1 )then
      write(out_lun,*)'num_extrapt_frame_list = ',num_extrapt_frame_list
   endif
   if ( num_extrapt_frame_list > 0 )then
      allocate(extrapt_frame_list(5,num_extrapt_frame_list),stat=ier)
      if ( ier /= 0 )then
         write(out_lun,*)'FH_SITES_readfile: allocate fails!'
         stop
      endif
      read(free_lun,'(a)',iostat=ios)line
      if ( ios /= 0 .or. line(1:4) /= "tail" )then
         write(out_lun,*)'FH_SITES_readfile: bad line: (tail)',line(1:20)
         stop
      endif ! ( ios /= 0 )
      do n = 1,num_extrapt_frame_list
         read(free_lun,*)(extrapt_frame_list(j,n),j=1,5)
      enddo
   endif
   ! extra points list and coords
   read(free_lun,'(a)',iostat=ios)line
   if ( ios /= 0 .or. line(1:14) /= "extra_pts_list" )then
      write(out_lun,*)'FH_SITES_readfile: bad line: (extra_pts_list)',line(1:20)
      stop
   endif ! ( ios /= 0 )
   read(free_lun,'(a)',iostat=ios)line
   if ( ios /= 0 .or. line(1:8) /= "num_list" )then
      write(out_lun,*)'FH_SITES_readfile: bad line: (num_list)',line(1:20)
      stop
   endif ! ( ios /= 0 )
   read(line(9:),*)num_extrapt_list
   if ( verbose == 1 )then
      write(out_lun,*)'num_extrapt_list = ',num_extrapt_list
   endif
   if ( num_extrapt_list > 0 )then
      read(free_lun,'(a)',iostat=ios)line
      if ( ios /= 0 .or. line(1:8) /= "site_num" )then
         write(out_lun,*)'FH_SITES_readfile: bad line: (site_num)',line(1:20)
         stop
      endif ! ( ios /= 0 )
      allocate(extra_pt_list(3,num_extrapt_list),stat=ier1)
      allocate(extra_pt_loc_crd(3,num_extrapt_list),stat=ier2)
      if ( ier1 /= 0 .or. ier2 /= 0 )then
         write(out_lun,*)'FH_SITES_readfile: allocate fails!'
         stop
      endif
      ! conversion since loc crds in file are in angstroms
!     bohr=0.529177249d0
      bohr=0.52917721092d0
      do n = 1,num_extrapt_list
         read(free_lun,*)(extra_pt_list(j,n),j=1,3),xloc,yloc,zloc
         extra_pt_loc_crd(1,n) = xloc / bohr
         extra_pt_loc_crd(2,n) = yloc / bohr
         extra_pt_loc_crd(3,n) = zloc / bohr
      enddo
   endif
   ! masked list for direct sum interaction list
   read(free_lun,'(a)',iostat=ios)line
   if ( ios /= 0 .or. line(1:14) /= "site_mask_list" )then
      write(out_lun,*)'FH_SITES_readfile: bad line: (site_mask_list)',line(1:20)
      stop
   endif ! ( ios /= 0 )
   allocate(num_sites_masked_by_site(num_sites),stat=ier1)
   allocate(masked_sitelist_offset_of_site(num_sites),stat=ier2)
   if ( ier1 /= 0 .or. ier2 /= 0 )then
      write(out_lun,*)'FH_SITES_readfile: allocate fails!'
      stop
   endif
   do n = 1,num_sites
      read(free_lun,*)j,num_sites_masked_by_site(j)
      if ( n /= j )then
         write(out_lun,*)'FH_SITES_readfile: masklist site out of order'
         stop
      endif
   enddo 
   masked_sitelist_offset_of_site(1) = 0
   do n = 2,num_sites
      masked_sitelist_offset_of_site(n) =  &
          masked_sitelist_offset_of_site(n-1) + num_sites_masked_by_site(n-1)
   enddo
   size_masked_sitelist = masked_sitelist_offset_of_site(num_sites) + &
                          num_sites_masked_by_site(num_sites)
   if ( verbose == 1 )then
      write(out_lun,*)'size_masked_sitelist = ',size_masked_sitelist
   endif
   allocate(masked_sitelist(size_masked_sitelist),stat=ier)
   if ( ier /= 0 )then
      write(out_lun,*)'FH_SITES_readfile: allocate fails!'
      stop
   endif
   do n = 1,num_sites
      num = num_sites_masked_by_site(n)
      off = masked_sitelist_offset_of_site(n)
      read(free_lun,*)(masked_sitelist(off+j),j=1,num)
   enddo

   close(free_lun)

   return

end subroutine FH_SITES_readfile
!----------------------------------------------------------------
subroutine FH_SITES_build()

   use user, only : pbc
   use atoms, only : num_atoms,atomic_crds

   implicit none

   integer j,k,n,nmol
   double precision x(3)
   !  clear the sitecrd array
   call UTIL_zero_real_array(site_crd,3*num_sites)
   !  first copy the atomic coords into site coord array
   call ATMSITE_atom_to_sitecrds(num_atoms,atomic_crds,site_crd,  &
                                site_index_of_atom)
   ! next define the centers..weighted averages of atomic positions
   call ATMSITE_build_centers(num_atoms,num_cenwt_list, &
                              center_site_deflist,center_site_weight,site_crd)
   ! clear frame arrays
   call UTIL_zero_real_array(frame_def_pt,3*2*num_sites)
   call UTIL_zero_real_array(frame,3*3*num_sites)
   call UTIL_zero_int_array(frame_def_pt_markflag,num_sites)
   call UTIL_zero_int_array(frame_doneflag,num_sites)
   ! fill the regular framelist
   if ( num_reg_frame_list > 0 )then
      ! first the def points
      call ATMSITE_build_frame_def_pts( &
         num_reg_frame_list,reg_frame_list,site_crd, &
         num_sites,frame_def_pt,frame_def_pt_markflag)
      ! next complete the frames
      call ATMSITE_defpoints_to_frames(num_sites,frame_axis_order, &
               frame_def_pt,frame,frame_def_pt_markflag,frame_doneflag)
  endif
  ! fill any extra sites
   if ( num_extrapt_list > 0 )then
      call ATMSITE_build_extrasites(num_extrapt_list,extra_pt_list, &
                                    extra_pt_loc_crd,frame,site_crd)
   endif
   ! fill the extrapts framelist
   if ( num_extrapt_frame_list > 0 )then
      ! first the def points
      call ATMSITE_build_frame_def_pts( &
         num_extrapt_frame_list,extrapt_frame_list,site_crd, &
         num_sites,frame_def_pt,frame_def_pt_markflag)
      ! next complete the frames
      call ATMSITE_defpoints_to_frames(num_sites,frame_axis_order, &
               frame_def_pt,frame,frame_def_pt_markflag,frame_doneflag)
   endif
    ! print out ideal water monomer geoms overlaying actual
    nmol = num_sites / 5 !MIDPOINTS HERE!!!
    do k = 1,nmol
       n = 5*(k-1) + 1
       !write(18,'(a2,1x,i4,3f12.6)')'O ',k,0.529177249d0*site_crd(1,n), &
                                    !0.529177249d0*site_crd(2,n), &
                                    !0.529177249d0*site_crd(3,n)
       do j = 1,3
          x(j) = 0.529177249d0*site_crd(j,n) + 0.584463d0*frame(j,3,n) + &
                                           0.764099d0*frame(j,1,n)
       enddo 
       !write(18,'(a2,1x,i4,3f12.6)')'H1',k,x(1),x(2),x(3)
       do j = 1,3
          x(j) = 0.529177249d0*site_crd(j,n) + 0.584463d0*frame(j,3,n) - &
                                           0.764099d0*frame(j,1,n)
       enddo 
       !write(18,'(a2,1x,i4,3f12.6)')'H2',k,x(1),x(2),x(3)
    enddo
    !do n = 1,num_sites
      !write(17,'(i4,3(1x,f14.9))')n,0.529177249d0*site_crd(1,n), &
                                    !0.529177249d0*site_crd(2,n), &
                                    !0.529177249d0*site_crd(3,n)
   !enddo
end subroutine FH_SITES_build
!----------------------------------------------------------------
subroutine FH_SITES_de_drotsite_to_sitefrc(pass,virial)

   implicit none

   integer,intent(in) :: pass
   double precision, intent(inout)    :: virial(3,3)

   ! convert de dsite_rot to effect of rotating about frame defining axes
   call ATMSITE_accum_de_dframe_rot(  &
               num_sites,num_frames, &
               frame_index_of_site,frame_axis_order,  &
               de_drotsite,frame, &
               frame_def_pt,   &
               de_drotframe)
   ! from these find gradient wrt change in def point positions
   call ATMSITE_accum_de_ddefpts(   &
               num_frames,frame_axis_order,   &
               de_drotframe,frame,frame_def_pt,  &
               de_ddef_pt)

   if ( pass == 1 )then
      call ATMSITE_de_ddefpts_to_sites(  &
                   num_reg_frame_list,reg_frame_list, &
                   site_crd,site_frc,  &
                   num_frames,de_ddef_pt, virial)
      call ATMSITE_de_ddefpts_to_sites(  &
                   num_extrapt_frame_list,extrapt_frame_list, &
                   site_crd,site_frc,  &
                   num_frames,de_ddef_pt, virial)
   else
      call ATMSITE_de_ddefpts_to_sites(  &
                   num_reg_frame_list,reg_frame_list, &
                   site_crd,site_frc,  &
                   num_frames,de_ddef_pt, virial)
   endif
end subroutine FH_SITES_de_drotsite_to_sitefrc
!----------------------------------------------------------------
subroutine FH_SITES_add_extrapt_torque()

   implicit none

   ! first clear de_drotsite
   de_drotsite = 0.d0
   call ATMSITE_add_extra_site_contrib(  &
                num_extrapt_list,extra_pt_list,site_crd, &
                site_frc,de_drotsite)

end subroutine FH_SITES_add_extrapt_torque
!----------------------------------------------------------------
subroutine FH_SITES_sitefrc_to_frc()

   use atoms, only : num_atoms,atomic_frcs

   implicit none

   ! weighted center contribution
   call ATMSITE_trans_de_dcens_to_atoms( &
                               num_atoms,num_cenwt_list,  &
                               center_site_deflist,center_site_weight, &
                               site_frc)
   ! pullback to atoms
   atomic_frcs = 0.d0
   call ATMSITE_sitefrc_to_frc(num_atoms,site_frc,atomic_frcs, &
                               site_index_of_atom)
end subroutine FH_SITES_sitefrc_to_frc
!----------------------------------------------------------------
subroutine FH_SITES_deallocate()

   implicit none

   if ( allocated(site_index_of_atom) )deallocate(site_index_of_atom)
   if ( allocated(site_crd) )deallocate(site_crd)
   if ( allocated(frac_site_crd) )deallocate(frac_site_crd)
   if ( allocated(map_site_crd) )deallocate(map_site_crd)
   if ( allocated(site_frc) )deallocate(site_frc)
   if ( allocated(frame_def_pt) )deallocate(frame_def_pt)
   if ( allocated(frame) )deallocate(frame)
   if ( allocated(frame_def_pt_markflag) )deallocate(frame_def_pt_markflag)
   if ( allocated(frame_doneflag) )deallocate(frame_doneflag)
   if ( allocated(aux_type_of_site) )deallocate(aux_type_of_site)
   if ( allocated(molnum_of_site) )deallocate(molnum_of_site)
   if ( allocated(frame_index_of_site) )deallocate(frame_index_of_site)
   if ( allocated(de_drotsite) )deallocate(de_drotsite)
   if ( allocated(de_drotframe) )deallocate(de_drotframe)
   if ( allocated(de_ddef_pt) )deallocate(de_ddef_pt)
   if ( allocated(center_site_deflist) )deallocate(center_site_deflist)
   if ( allocated(center_site_weight) )deallocate(center_site_weight)
   if ( allocated(reg_frame_list) )deallocate(reg_frame_list)
   if ( allocated(extrapt_frame_list) )deallocate(extrapt_frame_list)
   if ( allocated(extra_pt_list) )deallocate(extra_pt_list)
   if ( allocated(extra_pt_loc_crd) )deallocate(extra_pt_loc_crd)
   if ( allocated(num_sites_masked_by_site) ) &
              deallocate(num_sites_masked_by_site)
   if ( allocated(masked_sitelist_offset_of_site) ) &
              deallocate(masked_sitelist_offset_of_site)
   if ( allocated(masked_sitelist) )deallocate(masked_sitelist)
end subroutine FH_SITES_deallocate
!----------------------------------------------------------------
end module sites
