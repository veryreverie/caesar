submodule (caesar_unique_directions_module) caesar_unique_directions_submodule
  use caesar_dynamical_matrices_module
contains

module procedure new_UniqueDirection
  output%atom_id             = atom_id
  output%direction           = direction
  output%atomic_displacement = atomic_displacement
end procedure

module procedure calculate_unique_directions
  ! A description of the symmetries which map an atom onto itself.
  integer, parameter :: &
     & no_symmetry = 0, & ! Symmetries don't map between directions.
     & y_to_z      = 1, & ! Symmetries only map y to z.
     & x_to_one    = 2, & ! Symmetries only map x to one of y or z.
     & x_to_all    = 3    ! Symmetries map x to both y and z.
  
  ! Unique atom variables.
  integer              :: no_unique_atoms
  integer, allocatable :: unique_atoms(:)
  type(AtomData)       :: atom
  
  ! Unit vectors.
  type(IntVector) :: ex,ey,ez
  type(IntVector) :: unit_vectors(3)
  
  ! Displacement variables.
  type(String)     :: directions(6)
  type(RealVector) :: displacements(6)
  logical          :: direction_required(6)
  
  ! Symmetry information.
  integer                :: point_group
  type(SymmetryOperator) :: symmetry
  type(IntVector)        :: transformed_vectors(3)
  integer                :: dot_prod(3,3)
  type(IntVector)        :: transformed_x
  
  ! Unique direction variables.
  integer :: no_unique_directions
  
  ! Temporary variables
  integer :: i,j,k,l,ialloc
  
  ex = vec([1,0,0])
  ey = vec([0,1,0])
  ez = vec([0,0,1])
  unit_vectors = [ex, ey, ez]
  
  ! --------------------------------------------------
  ! Construct displacements of magnitude harmonic_displacement along the
  !    positive and negative directions along each lattice vector.
  ! --------------------------------------------------
  directions = [ str('+dx'), &
               & str('-dx'), &
               & str('+dy'), &
               & str('-dy'), &
               & str('+dz'), &
               & str('-dz')  ]
  do i=1,3
    ! Construct displacements in the positive directions.
    displacements(2*i-1) = transpose(structure%lattice) * unit_vectors(i)
    displacements(2*i-1) = displacements(2*i-1)  &
                       & * harmonic_displacement &
                       & / l2_norm(displacements(2*i-1))
    ! Construct displacements in the negative directions.
    displacements(2*i) = -displacements(2*i-1)
  enddo
  
  ! --------------------------------------------------
  ! Identify a minimal set of atoms from which the others can be constructed
  !    by symmetry.
  ! --------------------------------------------------
  allocate(unique_atoms(structure%no_atoms), stat=ialloc); call err(ialloc)
  no_unique_atoms = 0
  do_i : do i=1,structure%no_atoms
    ! Search for an already chosen atom, k, which is a copy of atom i.
    do j=1,size(structure%symmetries)
      do k=1,no_unique_atoms
        if (structure%symmetries(j)%atom_group*i == unique_atoms(k)) then
          ! If such an atom is found, atom i is a copy of an already chosen
          !    atom, and should not itself be chosen.
          cycle do_i
        endif
      enddo
    enddo
    ! If no such atom is found, atom i not a copy of an already chosen
    !    atom, and should be chosen.
    no_unique_atoms = no_unique_atoms + 1
    unique_atoms(no_unique_atoms) = i
  enddo do_i
  unique_atoms = unique_atoms(:no_unique_atoms)
  
  ! --------------------------------------------------
  ! For each atom in unique_atoms, identify a minimal set of displacements
  !   from which a complete basis can be constructed by symmetry.
  ! --------------------------------------------------
  allocate(output(no_unique_atoms*6), stat=ialloc); call err(ialloc)
  no_unique_directions = 0
  do i=1,no_unique_atoms
    atom = structure%atoms(unique_atoms(i))
    direction_required = .true.
    
    ! Identify the properties of the point-group.
    point_group = no_symmetry
    do j=1,size(structure%symmetries)
      symmetry = structure%symmetries(j)
      
      ! Ignore symmetries which do not map this atom to itself.
      if (symmetry%atom_group * atom%id() /= atom%id()) then
        cycle
      endif
      
      ! Calculate the dot products of the transformed vectors with
      !    other vectors.
      do k=1,3
        transformed_vectors(k) = symmetry%tensor * unit_vectors(k)
        do l=1,3
          dot_prod(l,k) = unit_vectors(l) * transformed_vectors(k)
        enddo
      enddo
      
      ! Check if any symmetries map a direction onto minus itself.
      ! If so, the displacement is only needed in one direction.
      do k=1,3
        if (dot_prod(k,k)==-1) then
          ! Set the displacement in the negative direction to not required.
          direction_required(2*k) = .false.
        endif
      enddo
      
      if (point_group == x_to_all) then
        ! The point group is already known to map x onto all directions.
        ! No further checks are required.
        cycle
      else
        if (dot_prod(2,1)/=0 .or. dot_prod(3,1)/=0) then
          ! The symmetry maps x onto another direction.
          if (point_group == x_to_one) then
            ! There is already a symmetry mapping x onto another direction.
            ! Check if x and its two transformed copies are
            !    linearly independent.
            if ( triple_product( unit_vectors(1),         &
               &                 transformed_x,           &
               &                 transformed_vectors(1) ) &
               & /= 0) then
              point_group = x_to_all
            endif
          else
            ! This is the first symmetry mapping x onto another direction.
            point_group = x_to_one
            transformed_x = transformed_vectors(1)
          endif
        elseif (dot_prod(3,2)/=0) then
          ! The symmetry maps y onto z.
          ! This is only relevant if no other mappings have been found.
          if (point_group == no_symmetry) then
            point_group = y_to_z
          endif
        endif
      endif
    enddo
    
    ! Update which directions are required based on the point group.
    if (point_group==no_symmetry) then
      continue
    elseif (point_group==y_to_z) then
      ! There is a symmetry mapping y onto z,
      !    so displacements in the z direction are not needed.
      direction_required(5:6) = .false.
    elseif (point_group==x_to_one) then
      if (transformed_x*vec([0,0,1])/=0) then
        ! There is a symmetry mapping x onto z,
        !    so displacements in the z direction are not needed.
        direction_required(5:6) = .false.
      else
        ! There is a symmetry mapping x onto y,
        !    so displacements in the y direction are not needed.
        direction_required(3:4) = .false.
      endif
    elseif (point_group==x_to_all) then
      ! There are symmetries mapping x onto y and z,
      !    so displacements are not needed in the y or z directions.
      direction_required(3:6) = .false.
    else
      call print_line(CODE_ERROR//': Illegal enum for point_group.')
      call err()
    endif
    
    ! Construct output.
    do j=1,6
      if (direction_required(j)) then
        no_unique_directions = no_unique_directions + 1
        output(no_unique_directions) = UniqueDirection( atom%id(),     &
                                                      & directions(j), &
                                                      & displacements(j) )
      endif
    enddo
  enddo
  output = output(:no_unique_directions)
  
  ! --------------------------------------------------
  ! Check output.
  ! --------------------------------------------------
  call check_unique_directions( output,           &
                              & structure,        &
                              & displacements(1), &
                              & displacements(3), &
                              & displacements(5))
end procedure

module procedure check_unique_directions
  type(AtomData) :: atom_1
  type(AtomData) :: atom_1p
  
  type(RealVector)              :: x
  type(RealMatrix), allocatable :: xx(:)
  
  real(dp) :: min_determinant
  
  integer :: i,j,ialloc
  
  min_determinant = 1.0e-7_dp * abs(triple_product(a,b,c))**2
  
  ! Construct xx = sum[ (S.x)^(S.x) ].
  allocate(xx(structure%no_atoms), stat=ialloc); call err(ialloc)
  xx = dblemat(zeroes(3,3))
  do i=1,size(structure%symmetries)
    do j=1,size(unique_directions)
      atom_1 = structure%atoms(unique_directions(j)%atom_id)
      atom_1p = structure%atoms( structure%symmetries(i)%atom_group &
                             & * atom_1%id())
      x = structure%symmetries(i)%cartesian_tensor &
      & * unique_directions(j)%atomic_displacement
      xx(atom_1p%id()) = xx(atom_1p%id()) + outer_product(x,x)
    enddo
  enddo
  
  ! Check that xx can be inverted.
  do i=1,structure%no_atoms
    if (determinant(xx(i))<min_determinant) then
      call print_line(CODE_ERROR//": error generating unique directions; &
         &x'^x' can't be inverted.")
      call print_lines(xx(i))
      call err()
    endif
  enddo
end procedure

module procedure new_CartesianDisplacement_UniqueDirection
  type(RealVector), allocatable :: displacements(:)
  
  integer :: ialloc
  
  allocate(displacements(structure%no_atoms), stat=ialloc); call err(ialloc)
  displacements = dblevec(zeroes(3))
  displacements(input%atom_id) = input%atomic_displacement
  
  this = CartesianDisplacement(displacements)
end procedure

module procedure read_UniqueDirection
  type(String), allocatable :: line(:)
  
  integer          :: atom_id
  type(String)     :: direction
  type(RealVector) :: atomic_displacement
  
  select type(this); type is(UniqueDirection)
    if (size(input)/=3) then
      call print_line(ERROR//': Unable to parse UniqueDirection from strings:')
      call print_lines(input)
      call err()
    endif
    
    ! Read in atom ID.
    line = split_line(input(1))
    atom_id = int(line(3))
    
    ! Read in displacement direction.
    line = split_line(input(2))
    direction = line(3)
    
    ! Read in displacement.
    line = split_line(input(3))
    atomic_displacement = vec(dble(line(3:5)))
    
    this = UniqueDirection(atom_id,direction,atomic_displacement)
  class default
    call err()
  end select
end procedure

module procedure write_UniqueDirection
  select type(this); type is(UniqueDirection)
    output = [ 'Atom         : '//this%atom_id,     &
             & 'Direction    : '//this%direction,   &
             & 'Displacement : '//this%atomic_displacement ]
  class default
    call err()
  end select
end procedure

module procedure new_UniqueDirection_Strings
  call this%read(input)
end procedure

module procedure new_UniqueDirection_StringArray
  this = UniqueDirection(str(input))
end procedure
end submodule
