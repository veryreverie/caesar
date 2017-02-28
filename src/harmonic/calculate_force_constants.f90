! This program determines the atomic displacements required to construct the
! matrix of force constants

! The symmetry operations to consider are nicely outlined here:
! http://www.homepages.ucl.ac.uk/~ucfbdxa/phon/node4.html

module calculate_force_constants_module
  implicit none
contains

function calculate_force_constants(structure,structure_sc) &
   & result(force_constants)
  use constants,      only : dp
  use linear_algebra, only : invert
  use structure_module
  use string_module
  implicit none
  
  ! Inputs
  type(StructureData), intent(in) :: structure
  type(StructureData), intent(in) :: structure_sc
  integer, allocatable            :: force_constants(:,:)
  
  ! Input variables
  real(dp),allocatable :: offset(:,:)
  
  ! Parameters
  real(dp), parameter :: x(3) = (/ 1.0_dp, 0.0_dp, 0.0_dp /)
  real(dp), parameter :: tol=1.0e-5_dp
  real(dp), parameter :: tol_mod=1.0e-8_dp
  
  ! Working variables
  real(dp)             :: rot_pos(3),reduced_pos(3),rot_pos_frac(3)
  integer, allocatable :: related(:)
  integer              :: y_related
  integer              :: z_related
  integer              :: no_constants
  
  real(dp), allocatable :: frac_atoms(:,:)
  
  ! Temporary variables
  integer :: i,j,k
  
  ! Transform offsets to primitive cell coordinates
  allocate(offset(3,structure_sc%no_symmetries))
  offset = matmul(structure_sc%supercell%recip_supercell,structure_sc%offsets) &
         & / structure_sc%supercell%sc_size
  
  allocate(frac_atoms(3,structure_sc%no_atoms))
  frac_atoms = matmul(structure%recip_lattice,structure_sc%atoms)
  
  ! Check which Cartesian directions are needed
  ! y_related = 1 if any rotation maps (0,1,0) onto (1,0,0)
  ! z_related = 1 if any rotation maps (0,0,1) onto (1,0,0)
  y_related = 0
  z_related = 0
  do i=1,structure_sc%no_symmetries
    if (all(abs(structure_sc%rotation_matrices(:,2,i)-x)<tol)) then
      y_related = 1
    endif
    if (all(abs(structure_sc%rotation_matrices(:,3,i)-x)<tol)) then
      z_related = 1
    endif
  enddo

  ! Apply point group to atomic positions
  allocate(related(structure_sc%no_atoms))
  related = 0
  do i=1,structure_sc%no_atoms
    do_j : do j=1,structure_sc%no_symmetries
      ! Work out the rotated position of atom i
      rot_pos=matmul( structure_sc%rotation_matrices(:,:,j), &
                    & structure_sc%atoms(:,i))
      rot_pos_frac = matmul(structure%recip_lattice,rot_pos) + offset(:,j)
      do k=1,i-1
        ! Calculate the distance between atom k and the rotated atom i
        reduced_pos = rot_pos_frac(:)-frac_atoms(:,k)-tol_mod
        ! Check if they are at the same position (possibly in different cells)
        if (all(dabs(reduced_pos-nint(reduced_pos))<tol)) then
          related(i) = 1
          exit do_j
        endif 
      enddo
    enddo do_j
  enddo
  
  no_constants = (3-y_related-z_related)*(size(related)-sum(related))
  allocate(force_constants(2,no_constants))
  j = 1
  do i=1,structure_sc%no_atoms
    if (related(i) == 0) then
      force_constants(:,j) = (/ i, 1 /)
      j = j + 1
      
      if (y_related == 0) then
        force_constants(:,j) = (/ i, 2 /)
        j = j + 1
      endif
      
      if (z_related == 0) then
        force_constants(:,j) = (/ i, 3 /)
        j = j + 1
      endif
    endif
  enddo
end function
end module
