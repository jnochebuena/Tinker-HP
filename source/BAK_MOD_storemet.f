c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #####################################################################
c     ##                                                                 ##
c     ##  module storemet   --  stored variables for HMC                 ##
c     ##                                                                 ##
c     #####################################################################
c
c
c     vold : array of velocities
c     aold : array of accelerations
c     
c      globold
c      locold
c      repartold
c      repartrecold
c      domlenold
c      domlenrecold number of reciprocal atoms in the reciprocal domains
c      domlenpoleold
c      domlenpolerecold
c      globrecold
c      locrecold
c      globrec1old
c      locrec1old
c      bufbegrecold
c      bufbegpoleold
c      bufbegold
c     buflen1old,buflen2old,buf1old,buf2old,bufbeg1old,bufbeg2old explicit direct-reciprocal atomic correspondance, 
c      nlocold
c      nblocold
c      nlocrecold
c      nblocrecold local + neighbors reciprocal number of atoms
c      nlocnlold local nl number of atoms
c      nblocrecdirold local + neighbors direct+reciprocal number of atoms
c
c      molecule 
c      nmolelocold,molculeglobold
c
c      VDW :
c        nvdwblocold,vdwglobold,nvdwlocnlold,vdwglobnlold
c        nvlstold,vlstold
c
c      BONDS:
c        nbondlocold,bndglobold
c
c      STRETCH-BEND:
c        nstrbndlocold,strbndglobold
c
c      ANGLE-ANGLE:
c        nanganglocold,angangglobold
c
c      OP-BENDING:
c        nopbendlocold,opbendglobold
c
c      OP-DIST:
c        nopdistlocold,opdistglobold
c
c      IMPROP:
c        niproplocold,impropglobold
c
c      IMPTOR:
c        nitorslocold,imptorglobold 
c
c      TORSION:
c        ntorslocold,torsglobold 
c
c      PITORSION:
c        npitorslocold,pitorsglobold 
c
c      STRETCH-TORSION:
c        nstrtorlocold,strtorglobold 
c
c      TORSION-TORSION:
c        ntortorlocold,tortorglobold 
c
c      ANGLE:
c        nanglelocold,angleglobold 
c
c      CHARGE:
c        nionreclocold,chgrecglobold 
c        nionlocold,chgglobold
c        nionlocnlold,chgglobnlold
c        nelstold,elstold
c
c      MULTIPOLE:
c        npolereclocold,polerecglobold 
c        npolelocold,poleglobold,polelocold
c        npolelocnlold,poleglobnlold
c        
c      POLARIZATION:
c        udaltold,upaltold,nualtold
c        uindold,uinpold
c      STAT:
c        istepmet,nstepmet
c
c        
      module storemet
      implicit none
      real*8 :: epotold,etotold
c
c     stored quantities
c
      integer :: nlocold,nblocold,nlocrecold
      integer :: nblocrecold,nlocnlold,nblocrecdirold
      integer, allocatable :: globold(:),locold(:),repartold(:)
      integer, allocatable :: ineignlold(:),locnlold(:) 
      integer, allocatable :: repartrecold(:), domlenold(:)
      integer, allocatable :: domlenrecold(:),domlenpoleold(:)
      integer, allocatable :: domlenpolerecold(:),globrecold(:)
      integer, allocatable :: locrecold(:),globrec1old(:)
      integer, allocatable :: locrec1old(:)
      integer, allocatable :: bufbegrecold(:)
      integer, allocatable :: bufbegpoleold(:),bufbegold(:)
      integer, allocatable :: buflen1old(:),buflen2old(:)
      integer, allocatable :: buf1old(:),buf2old(:)
      integer, allocatable :: bufbeg1old(:),bufbeg2old(:)
      integer, allocatable :: nmolelocold(:),molculeglobold(:)
      real*8, allocatable :: vold(:,:)
      real*8, allocatable :: aold(:,:)
c
c     BOND-STRETCHING
c
      integer :: nbondlocold
      integer, allocatable :: bndglobold(:)
c
c     STRETCH-BENDING
c
      integer :: nstrbndlocold
      integer, allocatable :: strbndglobold(:)
c
c     UREY-BRADLEY
c
      integer :: nureylocold
      integer, allocatable :: ureyglobold(:)
c
c     ANGLE-ANGLE
c
      integer :: nanganglocold
      integer, allocatable :: angangglobold(:)
c
c     OP-BENDING
c
      integer :: nopbendlocold
      integer, allocatable :: opbendglobold(:)
c
c     OP-DIST
c
      integer :: nopdistlocold
      integer, allocatable :: opdistglobold(:)
c
c     IMPROP
c
      integer :: niproplocold
      integer, allocatable :: impropglobold(:)
c
c     IMPTOR
c
      integer :: nitorslocold
      integer, allocatable :: imptorglobold(:)
c
c     TORSION
c
      integer :: ntorslocold
      integer, allocatable :: torsglobold(:)
c
c     PITORSION
c
      integer :: npitorslocold
      integer, allocatable :: pitorsglobold(:)
c
c     STRETCH-TORSION
c
      integer :: nstrtorlocold
      integer, allocatable :: strtorglobold(:)
c
c     TORSION-TORSION
c
      integer :: ntortorlocold
      integer, allocatable :: tortorglobold(:)
c
c     ANGLE
c
      integer :: nanglelocold
      integer, allocatable :: angleglobold(:)
c
c     CHARGE
c
      integer :: nionreclocold,nionlocold,nionlocnlold
      integer, allocatable :: chgrecglobold(:)
      integer, allocatable :: chgglobold(:)
      integer, allocatable :: chgglobnlold(:)
      integer, allocatable :: nelstold(:),elstold(:,:)
c
c     MULTIPOLE
c
      integer :: npolereclocold,npolelocold,npoleblocold,npolelocnlold
      integer, allocatable :: polerecglobold(:)
      integer, allocatable :: poleglobold(:)
      integer, allocatable :: polelocold(:)
      integer, allocatable :: poleglobnlold(:)
c 
c      POLARIZATION
c
      integer :: nualtold
      real*8, allocatable :: udaltold(:,:,:),upaltold(:,:,:)
      real*8, allocatable :: uindold(:,:),uinpold(:,:)
c
c     VDW
c
      integer :: nvdwblocold,nvdwlocnlold
      integer, allocatable :: vdwglobold(:)
      integer, allocatable :: vdwglobnlold(:)
      integer, allocatable :: nvlstold(:),vlstold(:,:)
c
c     STATS
c
      integer :: istepmet,nstepmet
      integer :: naccept,ntrial
      save

      contains
      subroutine allocmetro
      use angle
      use atoms
      use bitor
      use domdec
      use molcul
      use neigh
      use sizes
      use pitors
      use potent
      use tors
      use uprior
      implicit none
      real*8 etemp,energy
c
      if (allocated(globold)) deallocate (globold)
      allocate (globold(n))
      if (allocated(locold)) deallocate (locold)
      allocate (locold(n))
      if (allocated(repartold)) deallocate (repartold)
      allocate (repartold(n))
      if (allocated(repartrecold)) deallocate (repartrecold)
      allocate (repartrecold(n))
      if (allocated(domlenold)) deallocate (domlenold)
      allocate (domlenold(nproc))
      if (allocated(domlenrecold)) deallocate (domlenrecold)
      allocate (domlenrecold(nproc))
      if (allocated(domlenpoleold)) deallocate (domlenpoleold)
      allocate (domlenpoleold(nproc))
      if (allocated(domlenpolerecold)) deallocate (domlenpolerecold)
      allocate (domlenpolerecold(nproc))
      if (allocated(globrecold)) deallocate (globrecold)
      allocate (globrecold(n))
      if (allocated(locrecold)) deallocate (locrecold)
      allocate (locrecold(n))
      if (allocated(globrec1old)) deallocate (globrec1old)
      allocate (globrec1old(n))
      if (allocated(locrec1old)) deallocate (locrec1old)
      allocate (locrec1old(n))
      if (allocated(bufbegrecold)) deallocate (bufbegrecold)
      allocate (bufbegrecold(nproc))
      if (allocated(bufbegpoleold)) deallocate (bufbegpoleold)
      allocate (bufbegpoleold(nproc))
      if (allocated(bufbegold)) deallocate (bufbegold)
      allocate (bufbegold(nproc))
      if (allocated(buflen1old)) deallocate (buflen1old)
      allocate (buflen1old(nproc))
      if (allocated(buflen2old)) deallocate (buflen2old)
      allocate (buflen2old(nproc))
      if (allocated(bufbeg1old)) deallocate (bufbeg1old)
      allocate (bufbeg1old(nproc))
      if (allocated (buf1old)) deallocate(buf1old)
      allocate (buf1old(n))
      if (allocated (buf2old)) deallocate(buf2old)
      allocate (buf2old(n))
      if (allocated (bufbeg1old)) deallocate(bufbeg1old)
      allocate (bufbeg1old(nproc))
      if (allocated(bufbeg2old)) deallocate (bufbeg2old)
      allocate (bufbeg2old(nproc))
      if (allocated(vold)) deallocate (vold)
      allocate (vold(3,n))
      if (allocated(aold)) deallocate (aold)
      allocate (aold(3,n))
      if (allocated(molculeglobold)) deallocate (molculeglobold)
      allocate (molculeglobold(nmol))
      if (allocated(ineignlold)) deallocate (ineignlold)
      allocate (ineignlold(n))
      if (allocated(locnlold)) deallocate (locnlold)
      allocate (locnlold(n))



      if (use_vdw) then
        if (allocated(vdwglobold)) deallocate (vdwglobold)
        allocate (vdwglobold(n))
        if (allocated(vdwglobnlold)) deallocate (vdwglobnlold)
        allocate (vdwglobnlold(n))
      end if

      if (use_bond) then
       if (allocated(bndglobold)) deallocate(bndglobold)
       allocate (bndglobold(8*n))
      end if

      if (use_strbnd) then
        if (allocated(strbndglobold)) deallocate(strbndglobold)
        allocate (strbndglobold(nangle))
      end if
      
      if (use_urey) then
        if (allocated(ureyglobold)) deallocate(ureyglobold)
        allocate (ureyglobold(nangle))
      end if

      if (use_angang) then
        if (allocated(angangglobold)) deallocate(angangglobold)
        allocate (angangglobold(ntors))
      end if

      if (use_opbend) then
        if (allocated(opbendglobold)) deallocate(opbendglobold)
        allocate (opbendglobold(nangle))
      end if

      if (use_opdist) then
        if (allocated(opdistglobold)) deallocate(opdistglobold)
        allocate (opdistglobold(n))
      end if

      if (use_improp) then
        if (allocated(impropglobold)) deallocate(impropglobold)
        allocate (impropglobold(6*n))
      end if

      if (use_imptor) then
        if (allocated(imptorglobold)) deallocate(imptorglobold)
        allocate (imptorglobold(6*n))
      end if

      if (use_tors) then
        if (allocated(torsglobold)) deallocate(torsglobold)
        allocate (torsglobold(6*n))
      end if

      if (use_pitors) then 
        if (allocated(pitorsglobold)) deallocate(pitorsglobold)
        allocate (pitorsglobold(ntors))
      end if

      if (use_strtor) then
        if (allocated(strtorglobold)) deallocate(strtorglobold)
        allocate (strtorglobold(ntors))
      end if

      if (use_tortor) then
        if (allocated(tortorglobold)) deallocate(tortorglobold)
        allocate (tortorglobold(nbitor))
      end if

      if (use_angle) then
        if (allocated(angleglobold)) deallocate(angleglobold)
        allocate (angleglobold(4*n))
      end if

      if (use_charge) then
        if (allocated(chgrecglobold)) deallocate(chgrecglobold)
        allocate (chgrecglobold(n))
        if (allocated(chgglobold)) deallocate(chgglobold)
        allocate (chgglobold(n))
        if (allocated(chgglobnlold)) deallocate(chgglobnlold)
        allocate (chgglobnlold(n))
      end if

      if (use_polar) then
        if (allocated(udaltold)) deallocate(udaltold)
        allocate (udaltold(maxualt,3,n))
        if (allocated(upaltold)) deallocate(upaltold)
        allocate (upaltold(maxualt,3,n))
        if (allocated(uindold)) deallocate(uindold)
        allocate (uindold(3,n))
        if (allocated(uinpold)) deallocate(uinpold)
        allocate (uinpold(3,n))
      end if

      if (use_mpole) then
        if (allocated(polerecglobold)) deallocate(polerecglobold)
        allocate (polerecglobold(n))
        if (allocated(poleglobold)) deallocate(poleglobold)
        allocate (poleglobold(n))
        if (allocated(polelocold)) deallocate(polelocold)
        allocate (polelocold(n))
        if (allocated(poleglobnlold)) deallocate(poleglobnlold)
        allocate (poleglobnlold(n))
      end if
c

      if (use_vdw) then
        if (allocated(nvlstold)) deallocate (nvlstold)
        allocate (nvlstold(nlocnl))
        if (allocated(vlstold)) deallocate (vlstold)
        allocate (vlstold(maxvlst,nlocnl))
      end if
      if (use_mpole.or.use_charge) then
        if (allocated(nelstold)) deallocate (nelstold)
        allocate (nelstold(nlocnl))
        if (allocated(elstold)) deallocate (elstold)
        allocate (elstold(maxelst,nlocnl))
      end if

      ntrial = 0
      naccept = 0
      istepmet = 1
      etemp = energy()
      call reduceen(etemp)
      etemp = -1000
      call storemetro(etemp,etemp)


      end subroutine

      subroutine storemetro(epot,etot)
      use angang
      use angle
      use atoms
      use atmlst
      use atmtyp
      use bond
      use charge
      use domdec
      use energi
      use improp
      use imptor
      use moldyn
      use mpole
      use neigh
      use opbend
      use opdist
      use pitors
      use polar
      use potent
      use sizes
      use strbnd
      use strtor
      use timestat
      use tors
      use tortor
      use uprior
      use urey
      use vdw
      implicit none
      real*8 epot,etot

c
c     parallelism
c
      nlocold = nloc
      nblocold = nbloc
      nlocrecold = nlocrec
      nblocrecold = nblocrec
      nlocnlold = nlocnl
      nblocrecdirold = nblocrecdir
      globold = glob
      locold = loc
      ineignlold = ineignl
c      locnlold = locnl
      repartold = repart
      repartrecold = repartrec
      domlenold = domlen
      domlenrecold = domlenrec
      domlenpoleold = domlenpole
      domlenpolerecold = domlenpolerec
      globrecold = globrec
      locrecold = locrec
      globrec1old = globrec1
      locrec1old = locrec1
      bufbegrecold = bufbegrec
      bufbegpoleold = bufbegpole
      bufbegold = bufbeg
      buflen1old = buflen1
      buf1old = buf1
      buflen2old = buflen2
      buf2old = buf2
      bufbeg1old =  bufbeg1
      bufbeg2old = bufbeg2

      epotold=epot
      etotold=etot
c
c     positions, speed, mass
c
      vold = v
      aold = a
      xold = x
      yold = y
      zold = z
c
c     VDW 
c
      if (use_vdw) then
        nvdwblocold = nvdwbloc
        vdwglobold = vdwglob
        vdwglobnlold = vdwglobnl
        nvdwlocnlold = nvdwlocnl
          nvlstold = nvlst
          vlstold = vlst
      end if
c
c     BONDS
c
      if (use_bond) then
        bndglobold = bndglob
        nbondlocold = nbondloc
      end if
c
c     STRETCH-BEND
c
      if (use_strbnd) then
        strbndglobold = strbndglob
        nstrbndlocold= nstrbndloc
      end if
c
c     UREY-BRADLEY
c
      if (use_urey) then
        ureyglobold = ureyglob
        nureylocold = nureyloc
      end if
c
c     ANGlE-ANGLE
c
      if (use_angang) then
        angangglobold = angangglob
        nanganglocold = nangangloc
      end if
c
c     OP-BENDING
c
      if (use_opbend) then
        opbendglobold = opbendglob
        nopbendlocold = nopbendloc
      end if
c
c     OP-DIST
c
      if (use_opdist) then
        opdistglobold = opdistglob
        nopdistlocold = nopdistloc
      end if
c
c     IMPROP
c
      if (use_improp) then
        impropglobold = impropglob
        niproplocold = niproploc
      end if
c
c     IMPTOR
c
      if (use_imptor) then
        imptorglobold = imptorglob
        nitorslocold = nitorsloc
      end if
c
c     TORSION
c
      if (use_tors) then
        torsglobold = torsglob
        ntorslocold = ntorsloc
      end if
c
c     PITORSION
c
      if (use_pitors) then
        pitorsglobold = pitorsglob
        npitorslocold = npitorsloc
      end if
c
c     STRETCH-TORSION
c
      if (use_strtor) then
        strtorglobold = strtorglob
        nstrtorlocold = nstrtorloc
      end if
c
c     TORSION-TORSION
c
      if (use_tortor) then
        tortorglobold = tortorglob
        ntortorlocold = ntortorloc
      end if
c
c     ANGLE
c
      if (use_angle) then
        angleglobold = angleglob
        nanglelocold = nangleloc
      end if
c
c     CHARGE
c
      if (use_charge) then
        chgrecglobold = chgrecglob
        nionreclocold = nionrecloc
        chgglobold = chgglob
        nionlocold = nionloc
        chgglobnlold = chgglobnl
        nionlocnlold = nionlocnl
      end if
      if (use_charge.or.use_mpole) then
        nelstold = nelst
        elstold = elst
      end if
c
c     MULTIPOLE
c
      if (use_mpole) then
        polerecglobold = polerecglob
        npolereclocold= npolerecloc
        poleglobold = poleglob
        polelocold = poleloc
        npolelocold = npoleloc
        npoleblocold = npolebloc
        poleglobnlold = poleglobnl
        npolelocnlold = npolelocnl
      end if
c
c     POLARIZATION
c
      if (use_polar) then
        nualtold= nualt
        udaltold = udalt
        upaltold = upalt
        uindold = uind
        uinpold = uinp
      end if

c
      end subroutine
c
c
      subroutine recovermetro
      use angang
      use angle
      use atoms
      use atmlst
      use atmtyp
      use bond
      use charge
      use domdec
      use energi
      use improp
      use imptor
      use moldyn
      use mpole
      use neigh
      use opbend
      use opdist
      use pitors
      use polar
      use potent
      use sizes
      use strbnd
      use strtor
      use timestat
      use tors
      use tortor
      use uprior
      use urey
      use vdw
      implicit none

c
c     parallelism
c
      nloc = nlocold
      nbloc = nblocold
      nlocrec = nlocrecold
      nblocrec = nblocrecold
      nlocnl = nlocnlold
      nblocrecdir = nblocrecdirold
      glob = globold
      loc = locold
      ineignl = ineignlold
      repart = repartold
      repartrec = repartrecold
      domlen = domlenold
      domlenrec = domlenrecold
      domlenpole = domlenpoleold
      domlenpolerec = domlenpolerecold
      globrec = globrecold
      locrec = locrecold
      globrec1 = globrec1old
      locrec1 = locrec1old
      bufbegrec = bufbegrecold
      bufbegpole = bufbegpoleold
      bufbeg = bufbegold
      buflen1 = buflen1old
      buf1 = buf1old
      buflen2 = buflen2old
      buf2 = buf2old
      bufbeg1 = bufbeg1old
      bufbeg2 = bufbeg2old

c      epot=epotold
c
c     positions, speed, mass
c
      x = xold
      y = yold
      z = zold 
c
c     switch momentum variables
c
      v = -vold
      a = aold
c
c     VDW
c
      if (use_vdw) then
        nvdwbloc = nvdwblocold
        vdwglob = vdwglobold
        vdwglobnl = vdwglobnlold
        nvdwlocnl = nvdwlocnlold
          nvlst = nvlstold
          vlst = vlstold
      end if
c
c     BOND
c
      if (use_bond) then
        nbondloc = nbondlocold
        bndglob = bndglobold
      end if
c
c     STRETCH-BEND
c
      if (use_strbnd) then
        nstrbndloc = nstrbndlocold
        strbndglob = strbndglobold
      end if
c
c     UREY-BRADLEY
c
      if (use_urey) then
        nureyloc = nureylocold
        ureyglob = ureyglobold
      end if
c
c     ANGLE-ANGLE
c
      if (use_angang) then
        nangangloc = nanganglocold
        angangglob = angangglobold
      end if
c
c     OP-BENDING
c
      if (use_opbend) then
        nopbendloc = nopbendlocold
        opbendglob = opbendglobold
      end if
c
c     OP-DIST
c
      if (use_opdist) then
        nopdistloc = nopdistlocold
        opdistglob = opdistglobold
      end if
c
c     IMPROP
c
      if (use_improp) then
        niproploc = niproplocold
        impropglob = impropglobold
      end if
c
c     IMPTOR
c
      if (use_imptor) then
        nitorsloc = nitorslocold
        imptorglob = imptorglobold
      end if
c
c     TORSION
c
      if (use_tors) then
        ntorsloc = ntorslocold
        torsglob = torsglobold
      end if
c
c     PITORSION
c
      if (use_pitors) then
        npitorsloc = npitorslocold
        pitorsglob = pitorsglobold
      end if
c
c     STRETCH-TORSION
c
      if (use_strtor) then
        nstrtorloc = nstrtorlocold
        strtorglob = strtorglobold
      end if
c
c     TORSION-TORSION
c
      if (use_tortor) then
        ntortorloc = ntortorlocold
        tortorglob = tortorglobold
      end if
c
c     ANGLE
c
      if (use_angle) then
        nangleloc = nanglelocold
        angleglob = angleglobold
      end if
c
c     CHARGE
c
      if (use_charge) then
        nionrecloc = nionreclocold
        chgrecglob = chgrecglobold
        nionloc = nionlocold
        chgglob = chgglobold
        chgglobnl = chgglobnlold
        nionlocnl = nionlocnlold
      end if
      if (use_charge.or.use_mpole) then
        nelst = nelstold
        elst = elstold
      end if
c
c     MULTIPOLE
c
      if (use_mpole) then
        npolerecloc = npolereclocold
        polerecglob = polerecglobold
        npoleloc = npolelocold
        npolebloc = npoleblocold
        poleglob = poleglobold
        poleloc = polelocold
        poleglobnl = poleglobnlold
        npolelocnl = npolelocnlold
      end if
c
c     POLARIZATION
c
      if (use_polar) then
        nualt = nualtold
        udalt = udaltold
        upalt = upaltold
        uind = uindold
        uinp = uinpold
      end if
c
      end subroutine
      end
c
