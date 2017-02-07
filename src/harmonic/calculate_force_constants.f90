! This program determines the atomic displacements required to construct the
! matrix of force constants

! The symmetry operations to consider are nicely outlined here:
! http://www.homepages.ucl.ac.uk/~ucfbdxa/phon/node4.html

module calculate_force_constants_module
  implicit none
contains

function calculate_force_constants(structure,supercell_int,structure_sc) &
   & result(force_constants)
  use constants,        only : dp
  use linear_algebra,   only : inv_33
  use file_module,      only : open_read_file, open_write_file
  use structure_module, only : StructureData, read_structure_file, drop
  use string_module
  implicit none
  
  ! Inputs
  type(StructureData), intent(in) :: structure
  integer,             intent(in) :: supercell_int(3,3)
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
  logical, allocatable :: related(:)
  integer              :: y_related
  integer              :: z_related
  real(dp)             :: supercell(3,3)
  
  ! Temporary variables
  integer :: i,j,k
  
  supercell = supercell_int
  supercell = inv_33(transpose(supercell))
  
  ! Transform offsets to primitive cell coordinates
  allocate(offset(3,structure_sc%no_symmetries))
  offset = matmul(supercell,structure_sc%offsets)  

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
  related = .false.
  do i=1,structure_sc%no_atoms
    do_j : do j=1,structure_sc%no_symmetries
      rot_pos=matmul( structure_sc%rotation_matrices(:,:,j), &
                    & structure_sc%atoms(:,i))
      rot_pos_frac = matmul(structure%recip_lattice,rot_pos) + offset(:,j)
      do k=1,i-1
        reduced_pos = rot_pos_frac(:)-structure_sc%frac_atoms(:,k)-tol_mod
        if (all(dabs(reduced_pos-nint(reduced_pos))<tol)) then
          related(i)=.true.
          exit do_j
        endif 
      enddo
    enddo do_j
  enddo
  
  allocate(force_constants(2,3-y_related-z_related))
  j = 0
  do i=1,structure_sc%no_atoms
    if (.not. related(i)) then
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
