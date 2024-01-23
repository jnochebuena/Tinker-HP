c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ########################################################################################
c     ##                                                                                    ##
c     ##  subroutine metro  --  HYBRID Metropolis-BAOAB Langevin molecular dynamics step    ##
c     ##                                                                                    ##
c     ########################################################################################
c
c
c
c     literature reference:
c
c     Efficient molecular dynamics using geodesic integration 
c     and solvent-solute splitting B. Leimkuhler and C. Matthews,
c     Proceedings of the Royal Society A, 472: 20160138, 2016
c
c
      subroutine metro (istep,dt)
      use atmtyp
      use atoms
      use bath
      use cutoff
      use domdec
      use energi
      use freeze
      use langevin
      use mdstuf
      use moldyn
      use storemet
      use timestat
      use units
      use usage
      use mpi
      implicit none
      integer i,j,istep,iglob,ierr
      real*8 dt,dt_2,factor
      real*8 etot,eksum,epot,eold
      real*8 temp,pres
      real*8 part1,part2
      real*8 a1,a2,normal
      real*8 ekin(3,3)
      real*8 stress(3,3)
      real*8 time0,time1
      real*8, allocatable :: derivs(:,:)
      real*8 random,ratio,proba,test
      logical accept
     
      nstepmet = 1
c
c     set some time values for the dynamics integration
c
      dt_2 = 0.5d0 * dt
c
c     set time values and coefficients for BAOAB integration
c
      a1 = exp(-gamma*dt)
      a2 = sqrt((1-a1**2)*boltzmann*kelvin)
c
c     find quarter step velocities and half step positions via BAOAB recursion
c
      do i = 1, nloc
         iglob = glob(i)
         if (use(iglob)) then
            do j = 1, 3
               v(j,iglob) = v(j,iglob) + dt_2*a(j,iglob)
            end do
         end if
      end do
c
      do i = 1, nloc
        iglob = glob(i)
        if (use(iglob)) then
          x(iglob) = x(iglob) + v(1,iglob)*dt_2
          y(iglob) = y(iglob) + v(2,iglob)*dt_2
          z(iglob) = z(iglob) + v(3,iglob)*dt_2
        end if
      end do
c
c
c
c     compute random part
c
      deallocate (Rn)
      allocate (Rn(3,nloc))
      do i = 1, nloc
        do j = 1, 3
          Rn(j,i) = normal()
        end do
      end do
      do i = 1, nloc
         iglob = glob(i)
         if (use(iglob)) then
            do j = 1, 3
               v(j,iglob) = a1*v(j,iglob) + 
     $            a2*Rn(j,i)/sqrt(mass(iglob))
            end do
         end if
      end do
c
c
c     find full step positions via BAOAB recursion
c
      do i = 1, nloc
         iglob = glob(i)
         if (use(iglob)) then
            x(iglob) = x(iglob) + v(1,iglob)*dt_2
            y(iglob) = y(iglob) + v(2,iglob)*dt_2
            z(iglob) = z(iglob) + v(3,iglob)*dt_2
         end if
      end do
c
c
c
c     make half-step temperature and pressure corrections
c
c      call temper (dt,eksum,ekin,temp)
      call pressure2 (epot,temp)
c
c     Reassign the particules that have changed of domain
c
c     -> real space
c
      time0 = mpi_wtime()
c
      call reassign
c
c     -> reciprocal space
c
      call reassignpme(.false.)
      time1 = mpi_wtime()
      timereneig = timereneig + time1 - time0
c
c     communicate positions
c
      time0 = mpi_wtime()
      call commpos
      call commposrec
      time1 = mpi_wtime()
      timecommstep = timecommstep + time1 - time0
c
      allocate (derivs(3,nbloc))
      derivs = 0d0
c
      call reinitnl(istep)
c
      time0 = mpi_wtime()
      call mechanicstep(istep)
      time1 = mpi_wtime()
c
      timeparam = timeparam + time1 - time0
c
      time0 = mpi_wtime()
      call allocstep
      time1 = mpi_wtime()
      timeclear = timeclear  + time1 - time0
c
c     rebuild the neighbor lists
c
      if (use_list) call nblist(istep)
c
c     get the potential energy and atomic forces
c
      call gradient (epot,derivs)
c
c     MPI : get total energy
c
      call reduceen(epot)
c
c     communicate forces
c
      call commforces(derivs)
c
c     use Newton's second law to get the next accelerations;
c     find the full-step velocities using the BAOAB recursion
c
      do i = 1, nloc
         iglob = glob(i)
         if (use(iglob)) then
            do j = 1, 3
               a(j,iglob) = -convert * derivs(j,i)/mass(iglob)
               v(j,iglob) = v(j,iglob) + dt_2*a(j,iglob)
            end do
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (derivs)
cc
cc     metropolis-hastings criterion
cc
c      if (rank.eq.0) then
c        write(*,*) 'epot = ',epot,'epotold = ',epotold
cc        ratio = exp(-boltzmann*kelvin*(epot-epotold))
c        ratio = exp(-boltzmann*kelvin*(etot-etotold))
c        proba = min(1d0,ratio)
c        test = random()
c        write(*,*) 'proba = ',proba,test
c        write(*,*) 'istepmet = ',istepmet
c        accept = (test.le.proba)
c      end if
c      if (rank.eq.0) write(*,*) 'tx accept =',dble(naccept)/dble(ntrial)
c      call MPI_BCAST(accept,1,MPI_LOGICAL,0,COMM_TINKER,ierr)
cc      if (accept) then
cc        if (rank.eq.0) write(*,*) 'ACCEPT'
ccc
ccc     store datas of this configuration
ccc
cc        call storemetro(epot)
cc        istepmet = 1
cc      end if
c      if (istepmet.eq.nstepmet) then
c        ntrial = ntrial + 1
c        if (accept) then
c          if (rank.eq.0) write(*,*) 'ACCEPT',epot
c          naccept = naccept + 1
c          call storemetro(epot,etot)
c          istepmet = 1
c        else
c            if (rank.eq.0) write(*,*) 'REJECT'
c            call recovermetro
c            istepmet = 1
c        end if
c      else
c        istepmet = istepmet + 1
c      end if
c
c     find the constraint-corrected full-step velocities
c
c
      call temper (dt,eksum,ekin,temp)
      call pressure (dt,epot,ekin,temp,pres,stress,istep)
c
c     total energy is sum of kinetic and potential energies
c
      etot = eksum + esum
c
c     metropolis-hastings criterion
c
      if (rank.eq.0) then
        write(*,*) 'epot = ',epot,'epotold = ',epotold
c        ratio = exp(-boltzmann*kelvin*(epot-epotold))
        ratio = exp(-boltzmann*kelvin*(etot-etotold))
        proba = min(1d0,ratio)
        test = random()
        write(*,*) 'proba = ',proba,test
        write(*,*) 'istepmet = ',istepmet
        accept = (test.le.proba)
      end if
      if (rank.eq.0) write(*,*) 'tx accept =',dble(naccept)/dble(ntrial)
      call MPI_BCAST(accept,1,MPI_LOGICAL,0,COMM_TINKER,ierr)
c      if (accept) then
c        if (rank.eq.0) write(*,*) 'ACCEPT'
cc
cc     store datas of this configuration
cc
c        call storemetro(epot)
c        istepmet = 1
c      end if
      if (istepmet.eq.nstepmet) then
        ntrial = ntrial + 1
        if (accept) then
          if (rank.eq.0) write(*,*) 'ACCEPT',epot
          naccept = naccept + 1
          call storemetro(epot,etot)
          istepmet = 1
        else
            if (rank.eq.0) write(*,*) 'REJECT'
            call recovermetro
            istepmet = 1
        end if
      else
        istepmet = istepmet + 1
      end if
c
c     compute statistics and save trajectory for this step
c
c      call mdstat (istep,dt,etot,epot,eksum,temp,pres)
c      call mdsave (istep,dt,epot)
      call mdrest (istep)
      return
      end
