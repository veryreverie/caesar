module atom_mapping_module
contains

! ----------------------------------------------------------------------
! Finds a mapping between the atoms in two versions of the same structure.
! ----------------------------------------------------------------------
function atom_mapping(structure_a,structure_b) result(output)
  use constants, only : dp
  use string_module
  use structure_module
  use group_module
  use err_module
  implicit none
  
  type(StructureData), intent(in) :: structure_a
  type(StructureData), intent(in) :: structure_b
  type(Group)                     :: output
  
  ! Fractional atomic positions.
  real(dp), allocatable :: frac_pos_a(:,:)
  real(dp), allocatable :: frac_pos_b(:,:)
  real(dp)              :: difference(3)
  
  ! Output, but in array form.
  integer, allocatable :: mapping(:)
  
  ! Tolerance
  real(dp), parameter :: tol=1.0e-10_dp
  
  ! Temporary variables.
  integer :: i,j
  
  ! Check no_atoms is the same.
  if (structure_b%no_atoms/=structure_a%no_atoms) then
    call print_line('Atom counts do not match')
    call err()
  endif
  
  ! Check lattices are the same.
  if (.not. all(structure_b%lattice-structure_a%lattice < tol)) then
    call print_line('Lattices do not match.')
  endif
  
  ! Calculate fractional atomic positions.
  allocate(frac_pos_b(3,structure_b%no_atoms))
  allocate(frac_pos_a(3,structure_b%no_atoms))
  frac_pos_b = matmul(structure_b%recip_lattice,structure_b%atoms)
  frac_pos_a = matmul(structure_a%recip_lattice,structure_a%atoms)
  
  ! Find mapping between b and a atoms.
  allocate(mapping(structure_b%no_atoms))
  mapping = 0
  do i=1,structure_a%no_atoms
    do j=1,structure_b%no_atoms
      ! Check if a atom i and b atom j are the same, down to translation
      !    by lattice vectors.
      difference = frac_pos_a(:,i) - frac_pos_b(:,j)
      if (all(abs(difference-nint(difference)) < tol)) then
        if (mapping(i)/=0) then
          call print_line('Duplicate atom: atom '//i//' in a matches &
             &atom '//mapping(i)//' and atom '//j//' in b.')
          call err()
        endif
        mapping(i) = j
      endif
    enddo
  enddo
  
  do i=1,structure_b%no_atoms
    if (mapping(i)==0) then
      call print_line('Atom '//i//' in a has no equivalent in b.')
      call err()
    endif
    
    do j=1,i-1
      if (mapping(i)==mapping(j)) then
        call print_line('Duplicate atom: atom '//mapping(i)//' in b &
           &matches atom '//i//' and atom '//j//' in a.')
      endif
    enddo
  enddo
  
  output = mapping
end function
end module