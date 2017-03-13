module calculate_symmetry_group_module
contains

function calculate_symmetry_group(structure) result(output)
  use constants, only : dp
  use structure_module
  use group_module
  use string_module
  implicit none
  
  type(StructureData), intent(in) :: structure
  type(Group), allocatable        :: output(:)
  
  integer,  allocatable :: operations(:,:)
  
  ! Objects in fractional supercell co-ordinates.
  !   n.b. in these co-ordintates, supercell lattice vectors are unit vectors.
  !   This is a different convention to the scaled co-ordinates used elsewhere.
  real(dp), allocatable :: atom_pos_frac(:,:)
  real(dp)              :: rotation_frac(3,3)
  real(dp)              :: transformed_pos_frac(3)
  
  ! Distances between atoms and transformed atoms.
  real(dp)              :: delta(3)
  real(dp), allocatable :: distances(:)
  
  ! Temporary variables.
  integer :: i,j,k
  
  allocate(atom_pos_frac(3,structure%no_atoms))
  allocate(operations(structure%no_atoms,structure%no_symmetries))
  allocate(distances(structure%no_atoms))
  
  ! Transform atom positions into fractional supercell co-ordinates.
  atom_pos_frac = matmul(structure%recip_lattice,structure%atoms)
  
  ! Work out which atoms map to which atoms under each symmetry operation.
  do i=1,structure%no_symmetries
    ! Transform the rotations into fractional supercell co-ordinates.
    rotation_frac = matmul(matmul(                                      &
                                 & structure%recip_lattice,             &
                                 & structure%rotation_matrices(:,:,i)), &
                                 & transpose(structure%lattice))
    do j=1,structure%no_atoms
      ! Calculate the position of the transformed atom.
      transformed_pos_frac = matmul(rotation_frac,atom_pos_frac(:,j)) &
                         & + structure%offsets(:,i)
      
      ! Identify which atom is closest to the transformed position,
      !    modulo supercell lattice vectors.
      do k=1,structure%no_atoms
        delta = transformed_pos_frac - atom_pos_frac(:,k)
        distances(k) = norm2(delta-nint(delta))
      enddo
      operations(j,i) = minloc(distances,1)
    enddo
  enddo
  
  ! Check that each symmetry is one-to-one, and that mapped atoms are of the
  !    same species.
  do i=1,structure%no_symmetries
    do j=1,structure%no_atoms
      if (count(operations(:,i)==j)/=1) then
        call print_line('Error: symmetry operation not one-to-one.')
        stop
      endif
      
      if (structure%species(operations(j,i))/=structure%species(j)) then
        call print_line('Error: symmetry operation between different species.')
        stop
      endif
    enddo
  enddo
  
  allocate(output(size(operations,2)))
  do i=1,size(output)
    output(i) = operations(:,i)
  enddo
end function
end module
