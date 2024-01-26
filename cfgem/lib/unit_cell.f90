module unit_cell

   implicit none

   private

   integer, save                       :: extended_ucell, num_translates

   double precision, save              :: ucell(3,3), recip(3,3), volume, sphere

   double precision, allocatable, save :: ucell_translate(:,:)

   public       FH_UNIT_CELL_setup, FH_UNIT_CELL_deallocate

   public       ucell, recip, volume, sphere, ucell_translate, &
                extended_ucell, num_translates

contains

subroutine FH_UNIT_CELL_setup(out_lun)

   use user,only : unit_cell_length, &
                   unit_cell_angle, &
                   pbc, &
                   verbose, &
                   do_coulomb, &
                   do_overlap, &
                   extent_of_compact_hermites, &
                   exchange_cutoff, &
                   cut_extra

   use hermite, only : lower_density_extent,upper_density_extent

   implicit none

! Formal arguments:

   integer, intent(in)  :: out_lun

! Local variables:

   double precision     :: a, b, c, alpha, beta, gamma, R

   double precision, parameter  :: bohr = 0.52917721092d0

   integer              :: k, k1, k2, k3, m1, m2, m3, ier

   include 'interact_type.fh'

   a = unit_cell_length(1) 
   b = unit_cell_length(2)
   c = unit_cell_length(3)
   alpha = unit_cell_angle(1)
   beta = unit_cell_angle(2)
   gamma = unit_cell_angle(3)

   call FH_UNIT_CELL_fill(a,b,c,alpha,beta,gamma,ucell,recip, &
                                volume,sphere)
   if ( verbose .eq. 1 )then
      write(out_lun,*)'a,b,c,alpha,beta,gamma = ',a,b,c,alpha,beta,gamma
      write(out_lun,*)'sphere radius = ',sphere
   endif

   if (do_coulomb .eq. 1) then
     if (2.d0 * extent_of_compact_hermites + cut_extra >= sphere) then
      write(out_lun,*)'coulomb direct cut too large for unit cell!'
      stop
     end if
   else if (do_overlap .eq. 1) then
     if (exchange_cutoff/bohr + cut_extra >= sphere) then
      write(out_lun,*)'exchange direct cut too large for unit cell!'
      stop
     end if
   endif

   extended_ucell = 0

   ! get the translates

   allocate(ucell_translate(3,27),stat=ier)
   if ( ier /= 0 )then
     write(out_lun,*)'FH_UNIT_CELL_setup: cannot allocate ucell_translate'
     stop
   endif

   do k = 1,27
      k3 = k - 1
      m3 = k3 / 9
      k2 = k3 - 9*m3
      m2 = k2 / 3
      k1 = k2 - 3*m2
      m3 = m3 - 1
      m2 = m2 - 1
      m1 = k1 - 1
      ucell_translate(1,k) = ucell(1,1)*m1 + ucell(1,2)*m2 + &
                             ucell(1,3)*m3
      ucell_translate(2,k) = ucell(2,1)*m1 + ucell(2,2)*m2 + &
                             ucell(2,3)*m3
      ucell_translate(3,k) = ucell(3,1)*m1 + ucell(3,2)*m2 + &
                             ucell(3,3)*m3
   enddo

   return

end subroutine FH_UNIT_CELL_setup
!-----------------------------------------------------------
subroutine FH_UNIT_CELL_deallocate()

   implicit none

   if ( allocated(ucell_translate) )deallocate(ucell_translate)
end subroutine FH_UNIT_CELL_deallocate
!----------------------------------------------------------------

end module unit_cell
!-----------------------------------------------------------
subroutine FH_UNIT_CELL_fill(a,b,c,alpha,beta,gamma,ucell,recip, &
                               volume,sphere)
  implicit none 

  double precision,intent(in) :: a,b,c,alpha,beta,gamma
  double precision,intent(out) :: ucell(3,3),recip(3,3),volume,sphere

  double precision factor,distance,result,u23(3),u31(3),u12(3),reclng(3)
  integer i,j

  factor = 3.14159265358979323846d0/180.d0
  ucell(1,1) = a
  ucell(2,1) = 0.d0
  ucell(3,1) = 0.d0
  ucell(1,2) = b*cos(factor*gamma)
  ucell(2,2) = b*sin(factor*gamma)
  ucell(3,2) = 0.d0
  ucell(1,3) = c*cos(factor*beta)
  ucell(2,3) = (b*c*cos(factor*alpha)-ucell(1,3)*ucell(1,2))/ucell(2,2)
  ucell(3,3) = sqrt( c*c - ucell(1,3)*ucell(1,3) - ucell(2,3)*ucell(2,3) )

!  now get reciprocal vectors

  call cross(ucell(1,2),ucell(1,3),u23)
  call cross(ucell(1,3),ucell(1,1),u31)
  call cross(ucell(1,1),ucell(1,2),u12)
  call dot(ucell(1,1),u23,volume)
  do j = 1,3
    recip(j,1) = u23(j)/volume
    recip(j,2) = u31(j)/volume
    recip(j,3) = u12(j)/volume
  end do
  reclng(1) = 1.d0/sqrt( recip(1,1)*recip(1,1) +    &
                         recip(2,1)*recip(2,1) +   &
                         recip(3,1)*recip(3,1) )
  reclng(2) = 1.d0/sqrt( recip(1,2)*recip(1,2) +  &
                         recip(2,2)*recip(2,2) +  &
                         recip(3,2)*recip(3,2) )
  reclng(3) = 1.d0/sqrt( recip(1,3)*recip(1,3) +   &
                         recip(2,3)*recip(2,3) +   &
                         recip(3,3)*recip(3,3) )
! interfacial distances given by dot of direct,recip
! sphere is radius of largest sphere inscribed in unit cell
! the minimum image cutoff must be less than or equal to this

  sphere = a+b+c
  do i = 1,3
    call dot(recip(1,i),ucell(1,i),result)
    distance = result*reclng(i)
    if ( distance .lt. sphere )sphere = distance
  end do
  sphere = 0.5d0*sphere

end subroutine FH_UNIT_CELL_fill
!-----------------------------------------------------------
