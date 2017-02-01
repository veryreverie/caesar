! This program determines the atomic displacements required to construct the
! matrix of force constants

! The symmetry operations to consider are nicely outlined here:
! http://www.homepages.ucl.ac.uk/~ucfbdxa/phon/node4.html

module construct_matrix_force_cnsts_module
  implicit none
contains

subroutine construct_matrix_force_cnsts(structure,supercell_int,structure_sc, &
   & force_constants_filename)
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
  type(String),        intent(in) :: force_constants_filename
  
  ! Input variables
  real(dp),allocatable :: offset(:,:)
  
  ! Parameters
  real(dp), parameter :: x(3) = (/ 1.d0, 0.d0, 0.d0 /)
  real(dp), parameter :: tol=1.d-5
  real(dp), parameter :: tol_mod=1.d-8
  
  ! Working variables
  integer :: i,j,k
  real(dp) :: rot_pos(3),reduced_pos(3),rot_pos_frac(3)
  logical :: related,yrelated,zrelated
  real(dp) :: supercell(3,3)
  
  ! file units
  integer :: force_constants_file
  
  supercell = supercell_int
  supercell = inv_33(transpose(supercell))
  
  ! Transform offsets to primitive cell coordinates
  allocate(offset(3,structure_sc%no_symmetries))
  offset = matmul(supercell,structure_sc%offsets)  

  ! Check which Cartesian directions are needed
  yrelated=.false.
  zrelated=.false.
  do i=1,structure_sc%no_symmetries
    if (all(abs(structure_sc%rotation_matrices(:,2,i)-x)<tol)) then
      yrelated = .true.
    endif
    if (all(abs(structure_sc%rotation_matrices(:,3,i)-x)<tol)) then
      zrelated = .true.
    endif
  enddo ! i
  
  ! Prepare output file
  force_constants_file = open_write_file(force_constants_filename)

  ! Apply point group to atomic positions
  do i=1,structure_sc%no_atoms
    related=.false.
    do_j : do j=1,structure_sc%no_symmetries
      rot_pos=matmul( structure_sc%rotation_matrices(:,:,j), &
                    & structure_sc%atoms(:,i))
      rot_pos_frac = matmul(structure%recip_lattice,rot_pos) + offset(:,j)
      do k=1,i-1
        reduced_pos = rot_pos_frac(:)-structure_sc%frac_atoms(:,k)-tol_mod
        if (all(dabs(reduced_pos-nint(reduced_pos))<tol)) then
          related=.true.
          exit do_j
        endif 
      enddo
    enddo do_j
    
    if (.not. related) then
      write(force_constants_file,*) i, 1 
      if (.not. yrelated) write(force_constants_file,*) i, 2
      if (.not. zrelated) write(force_constants_file,*) i, 3
    endif
  enddo
  
  close(force_constants_file)
end subroutine
end module
