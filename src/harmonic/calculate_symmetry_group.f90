module calculate_symmetry_group_module
  use constants_module, only : dp
  use string_module
  use io_module
contains

function calculate_symmetry_group(structure) result(output)
  use structure_module
  use group_module
  use linear_algebra_module
  implicit none
  
  type(StructureData), intent(in) :: structure
  type(Group), allocatable        :: output(:)
  
  integer,  allocatable :: operations(:,:)
  
  ! Objects in fractional supercell co-ordinates.
  !   n.b. in these co-ordintates, supercell lattice vectors are unit vectors.
  !   This is a different convention to the scaled co-ordinates used elsewhere.
  type(RealVector), allocatable :: atom_pos_frac(:)
  type(RealVector)              :: transformed_pos_frac
  
  ! Distances between atoms and transformed atoms.
  type(RealVector)      :: delta
  real(dp), allocatable :: distances(:)
  
  ! Temporary variables.
  integer :: i,j,k
  
  allocate(atom_pos_frac(structure%no_atoms))
  allocate(operations(structure%no_atoms,structure%no_symmetries))
  allocate(distances(structure%no_atoms))
  
  ! Transform atom positions into fractional supercell co-ordinates.
  atom_pos_frac = structure%recip_lattice * structure%atoms
  
  ! Work out which atoms map to which atoms under each symmetry operation.
  do i=1,structure%no_symmetries
    do j=1,structure%no_atoms
      ! Calculate the position of the transformed atom.
      transformed_pos_frac = structure%rotations(i) * atom_pos_frac(j) &
                         & + structure%translations(i)
      
      ! Identify which atom is closest to the transformed position,
      !    modulo supercell lattice vectors.
      do k=1,structure%no_atoms
        delta = transformed_pos_frac - atom_pos_frac(k)
        distances(k) = l2_norm(delta-vec(nint(dble(delta))))
      enddo
      operations(j,i) = minloc(distances,1)
      
      ! Check that the transformed atom is acceptably close to its image.
      if (distances(operations(j,i))>1.0e-10_dp) then
        call err()
      endif
    enddo
  enddo
  
  ! Check that each symmetry is one-to-one, and that mapped atoms are of the
  !    same species.
  do i=1,structure%no_symmetries
    do j=1,structure%no_atoms
      if (count(operations(:,i)==j)/=1) then
        call print_line('Error: symmetry operation not one-to-one.')
        call err()
      endif
      
      if (structure%species(operations(j,i))/=structure%species(j)) then
        call print_line('Error: symmetry operation between different species.')
        call err()
      endif
    enddo
  enddo
  
  allocate(output(size(operations,2)))
  do i=1,size(output)
    output(i) = operations(:,i)
  enddo
end function
end module
