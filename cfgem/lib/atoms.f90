module atoms_gem

   implicit none

   private

   integer, save                        :: num_atoms
   double precision, allocatable, save  :: atomic_crds(:,:)
   double precision, allocatable, save  :: atomic_frcs(:,:)
   !double precision, allocatable        :: atomic_crds(:,:)
   !double precision, allocatable        :: atomic_frcs(:,:)

   public       FH_ATOMS_init, FH_ATOMS_deallocate

   public       num_atoms, atomic_crds, atomic_frcs

contains

subroutine FH_ATOMS_init(atom_cnt, atom_crds, out_lun)

   implicit none

! Formal arguments:

  integer,          intent(in)  :: atom_cnt
  double precision, intent(in)  :: atom_crds(3, atom_cnt)
  integer,          intent(in)  :: out_lun                ! printed output..

! Local variables:

   integer      :: alloc_failed

   ! NOTE units in use for this code - these all come from 2010 CODATA, and
   ! were adjusted to that standard to exactly match pmemd.gem..

   double precision, parameter  :: bohr = 0.52917721092d0
   double precision, parameter  :: bohr_per_angstrom = 1.d0 / bohr

   ! Convert coordinates from Angstrom to atomic units (Bohr)..

   num_atoms = atom_cnt ! Save atom count to module variable..

   if(allocated(atomic_crds)) deallocate(atomic_crds)
   if(allocated(atomic_frcs)) deallocate(atomic_frcs)
   allocate(atomic_crds(3, num_atoms), &
            atomic_frcs(3, num_atoms), &
            stat = alloc_failed)

   if (alloc_failed .ne. 0) then
      write(out_lun,*)'FH_ATOMS_init: allocate fails!!'
      stop
   endif
   
   atomic_crds(:,:) = atom_crds(:,:) * bohr_per_angstrom

   return

end subroutine FH_ATOMS_init

!-----------------------------------------------------------------------
subroutine FH_ATOMS_deallocate()

   implicit none

   if (allocated(atomic_crds)) deallocate(atomic_crds)
   if (allocated(atomic_frcs)) deallocate(atomic_frcs)

   return

end subroutine FH_ATOMS_deallocate

end module atoms_gem
