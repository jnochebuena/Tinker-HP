c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine unitcell  --  get periodic boundary conditions  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "unitcell" gets the periodic boundary box size and related
c     values from an external keyword file
c
c
      subroutine unitcell
      use bound
      use boxes
      use cutoff
      ! GAC
      use gemstuff
      use mdstuf, only: useGEM
      ! GAC
      use keys
      use iounit
      implicit none
      integer i,next
      real*8 boxmax,rmax
      logical nosymm
      character*20 keyword
      character*120 record
      character*120 string
c
c
c     set the default values for periodic boundary conditions
c
      use_bounds = .false.
c      use_replica = .false.
c
c     set the default values for the unitcell variables
c
      xbox = 0.0d0
      ybox = 0.0d0
      zbox = 0.0d0
      alpha = 0.0d0
      beta = 0.0d0
      gamma = 0.0d0
      orthogonal = .false.
      monoclinic = .false.
      triclinic = .false.
      octahedron = .false.
      spacegrp = '          '
      nosymm = .false.
c     S: JORGE
      useGEM = .false. 
      usarGEM = .false.
      onlyHal = .false.
      onlyXC = .false.
c     E:JORGE

c
c     get keywords containing crystal lattice dimensions
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:120)
         if (keyword(1:7) .eq. 'X-AXIS ') then
            read (string,*,err=10,end=10)  xbox
         else if (keyword(1:7) .eq. 'Y-AXIS ') then
            read (string,*,err=10,end=10)  ybox
         else if (keyword(1:7) .eq. 'Z-AXIS ') then
            read (string,*,err=10,end=10)  zbox
         else if (keyword(1:7) .eq. 'A-AXIS ') then
            read (string,*,err=10,end=10)  xbox
         else if (keyword(1:7) .eq. 'B-AXIS ') then
            read (string,*,err=10,end=10)  ybox
         else if (keyword(1:7) .eq. 'C-AXIS ') then
            read (string,*,err=10,end=10)  zbox
         else if (keyword(1:6) .eq. 'ALPHA ') then
            read (string,*,err=10,end=10)  alpha
         else if (keyword(1:5) .eq. 'BETA ') then
            read (string,*,err=10,end=10)  beta
         else if (keyword(1:6) .eq. 'GAMMA ') then
            read (string,*,err=10,end=10)  gamma
c     S:JORGE
         else if (keyword(1:8) .eq. 'USE-GEM ') then
            useGEM = .true.
            usarGEM = .true.
         !else if (keyword(1:9) .eq. 'USE-QMMM ') then
           ! QMMM = .true.
         else if (keyword(1:13) .eq. 'ONLY-HALGREN ') then
            onlyHal = .true.
         else if (keyword(1:8) .eq. 'ONLY-XC ') then
            onlyXC = .true.
c     E:JORGE
         else if (keyword(1:11) .eq. 'OCTAHEDRON ') then
            octahedron = .true.
         else if (keyword(1:11) .eq. 'SPACEGROUP ') then
            call getword (record,spacegrp,next)
         else if (keyword(1:12) .eq. 'NO-SYMMETRY ') then
            nosymm = .true.
         end if
   10    continue
      end do
c
c     use periodic boundary conditions if a cell was defined
c
      boxmax = max(xbox,ybox,zbox)
      if (boxmax .ne. 0.0d0)  use_bounds = .true.
c
c     set unspecified periodic boundary box lengths and angles
c
      if (use_bounds) then
         if (xbox .eq. 0.0d0)  xbox = boxmax
         if (ybox .eq. 0.0d0)  ybox = xbox
         if (zbox .eq. 0.0d0)  zbox = xbox
         if (alpha .eq. 0.0d0)  alpha = 90.0d0
         if (beta .eq. 0.0d0)  beta = 90.0d0
         if (gamma .eq. 0.0d0)  gamma = 90.0d0
      
         ! GAC
         if (useGEM) then
            pbc_box_size(1) = xbox
            pbc_box_size(2) = ybox
            pbc_box_size(3) = zbox
         endif
         ! GAC
c
c     determine the general periodic boundary lattice type
c
         if (nosymm) then
            triclinic = .true.
         else if (alpha.eq.90.0d0 .and. beta.eq.90.0d0
     &               .and. gamma.eq.90.0d0) then
            orthogonal = .true.
         else if (alpha.eq.90.0d0 .and. gamma.eq.90.0d0) then
            monoclinic = .true.
         else
            triclinic = .true.
         end if
      end if
c
c     check for proper use of truncated octahedron boundary
c
      if (octahedron) then
         if (xbox.eq.ybox .and. xbox.eq.zbox .and. orthogonal) then
            orthogonal = .false.
            monoclinic = .false.
            triclinic = .false.
         else
            write (iout,20)
   20       format (/,' UNITCELL  --  Truncated Octahedron',
     &                 ' Incompatible with Defined Cell')
            call fatal
         end if
      end if
   
      ! GAC fail safe to make sure box is orthogonal
      if (useGEM) then
         if (.not. orthogonal) then
            write(iout,*) 'GEM can only use rectangular boxes, exiting'
            call fatal
         endif
      endif
c
c     set the system extent for nonperiodic Ewald summation
c
      if (.not. use_bounds) then
         call extent (rmax)
         xbox = 2.0d0 * (rmax+ewaldcut)
         ybox = xbox
         zbox = xbox
         alpha = 90.0d0
         beta = 90.0d0
         gamma = 90.0d0
         orthogonal = .true.
c         call lattice
c         boundary = 'NONE'
c         dens = 0.75d0
      end if
c
      return
      end
