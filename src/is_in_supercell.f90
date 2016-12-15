module is_in_supercell_module
contains

! ----------------------------------------
! checks if pos is in the supercell
! ----------------------------------------
logical function is_in_supercell(pos, super_lattice)
  use constants,      only : dp
  use linear_algebra, only : inv_33
  implicit none
  
  real(dp), intent(in) :: pos(3)
  real(dp), intent(in) :: super_lattice(3,3)
  real(dp)             :: trans_super_lattice(3,3)
  real(dp)             :: inv_super_lattice(3,3)
  real(dp)             :: frac_pos(3)
  real(dp)             :: t
  real(dp)             :: a
  real(dp)             :: b
  real(dp)             :: c
  real(dp)             :: f1
  real(dp)             :: f2
  real(dp)             :: f3
  real(dp), parameter  :: tol=1.d-2
  
  trans_super_lattice = transpose(super_lattice)
  call inv_33(trans_super_lattice, inv_super_lattice)
  frac_pos(1:3) = pos(1)*inv_super_lattice(1:3,1)&
               &+ pos(2)*inv_super_lattice(1:3,2)&
               &+ pos(3)*inv_super_lattice(1:3,3)
  
  if (frac_pos(1)>-tol .and. frac_pos(1)<(1.d0-tol)) then
    if (frac_pos(2)>-tol .and. frac_pos(2)<(1.d0-tol)) then
      is_in_supercell = (frac_pos(3)>-tol .and. frac_pos(3)<(1.d0-tol))
    else
      is_in_supercell = .false.
    endif
  else
    is_in_supercell = .false.
  endif
end function
end module
