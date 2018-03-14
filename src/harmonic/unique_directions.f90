! ======================================================================
! Determines which directions are related by symmetry,
!    and thus which atoms should be displaced in which directions 
!    in order to construct the matrix of force constants.
! ======================================================================
! The symmetry operations to consider are nicely outlined here:
! http://www.homepages.ucl.ac.uk/~ucfbdxa/phon/node4.html
module unique_directions_module
  use common_module
  implicit none
  
  ! A direction in which an atom needs to be perturbed in order to map out
  !    the harmonic Born-Oppenheimer surface.
  type UniqueDirection
    ! Atom id.
    integer :: atom_id
    
    ! The direction in which the atom is displaced: '+dx', '-dx', ..., '-dz'.
    type(String) :: direction
    
    ! The vector through which the atom is displaced.
    type(RealVector) :: displacement
  end type
  
  interface UniqueDirection
    module procedure new_UniqueDirection_contents
  end interface
contains

! ----------------------------------------------------------------------
! Constructor
! ----------------------------------------------------------------------
function new_UniqueDirection_contents(atom_id,direction,displacement) &
   & result(output)
  implicit none
  
  integer,          intent(in) :: atom_id
  type(String),     intent(in) :: direction
  type(Realvector), intent(in) :: displacement
  type(UniqueDirection)        :: output
  
  output%atom_id      = atom_id
  output%direction    = direction
  output%displacement = displacement
end function

! ----------------------------------------------------------------------
! Reads an array of UniqueDirection from a file.
! ----------------------------------------------------------------------
function read_unique_directions_file(filename) result(output)
  implicit none
  
  type(String), intent(in)           :: filename
  type(UniqueDirection), allocatable :: output(:)
  
  type(IFile) :: unique_directions_file
  
  integer :: no_directions
  
  type(String), allocatable :: line(:)
  integer                   :: i,ialloc
  
  unique_directions_file = filename
  
  no_directions = size(unique_directions_file)-1
  
  allocate(output(no_directions), stat=ialloc); call err(ialloc)
  
  do i=1,no_directions
    line = unique_directions_file%split_line(i+1)
    output(i)%atom_id = int(line(2))
    output(i)%direction = line(4)
    output(i)%displacement = dble(line(6:8))
  enddo
end function

! ----------------------------------------------------------------------
! Writes a UniqueDirections to a file.
! ----------------------------------------------------------------------
subroutine write_unique_directions_file(unique_directions,filename)
  implicit none
  
  type(UniqueDirection), intent(in) :: unique_directions(:)
  type(String),          intent(in) :: filename
  
  type(OFile) :: unique_directions_file
  
  integer :: i
  
  unique_directions_file = filename
  call unique_directions_file%print_line('Atoms to be perturbed in order to &
     &map out the harmonic Born-Oppenheimer surface:')
  do i=1,size(unique_directions)
    call unique_directions_file%print_line(                &
       & 'atom: '//unique_directions(i)%atom_id         // &
       & ' direction: '//unique_directions(i)%direction // &
       & ' displacement: '//unique_directions(i)%displacement)
  enddo
end subroutine

! ----------------------------------------------------------------------
! Works out which atoms need to be perturbed in which directions
!    in order to map out the harmonic Born-Oppenheimer surface.
! ----------------------------------------------------------------------
function calculate_unique_directions(structure,harmonic_displacement) &
   & result(output)
  implicit none
  
  type(StructureData), intent(in)    :: structure
  real(dp),            intent(in)    :: harmonic_displacement
  type(UniqueDirection), allocatable :: output(:)
  
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
  type(IntVector)        :: rotated_vectors(3)
  integer                :: dot_prod(3,3)
  type(IntVector)        :: rotated_x
  
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
      
      ! Calculate the dot products of the rotated vectors with other vectors.
      do k=1,3
        rotated_vectors(k) = symmetry%rotation * unit_vectors(k)
        do l=1,3
          dot_prod(l,k) = unit_vectors(l) * rotated_vectors(k)
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
            ! Check if x and its two rotated copies are linearly independent.
            if ( triple_product( unit_vectors(1),    &
               &                 rotated_x,          &
               &                 rotated_vectors(1)) &
               & /= 0) then
              point_group = x_to_all
            endif
          else
            ! This is the first symmetry mapping x onto another direction.
            point_group = x_to_one
            rotated_x = rotated_vectors(1)
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
      if (rotated_x*vec([0,0,1])/=0) then
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
end function

! ----------------------------------------------------------------------
! Check that symmetry can be used to reconstruct all directions from
!    the unique directions.
! ----------------------------------------------------------------------
subroutine check_unique_directions(unique_directions,structure, &
   & a,b,c)
  implicit none
  
  type(UniqueDirection), intent(in) :: unique_directions(:)
  type(StructureData),   intent(in) :: structure
  type(RealVector),      intent(in) :: a
  type(RealVector),      intent(in) :: b
  type(RealVector),      intent(in) :: c
  
  type(AtomData) :: atom_1
  type(AtomData) :: atom_1p
  
  type(RealVector)              :: x
  type(RealMatrix), allocatable :: xx(:)
  
  real(dp) :: min_determinant
  
  integer :: i,j,ialloc
  
  min_determinant = 1.0e-7_dp * abs(triple_product(a,b,c))
  
  ! Construct xx = sum[ (S.x)^(S.x) ].
  allocate(xx(structure%no_atoms), stat=ialloc); call err(ialloc)
  xx = dblemat(zeroes(3,3))
  do i=1,size(structure%symmetries)
    do j=1,size(unique_directions)
      atom_1 = structure%atoms(unique_directions(j)%atom_id)
      atom_1p = structure%atoms( structure%symmetries(i)%atom_group &
                             & * atom_1%id())
      x = structure%symmetries(i)%cartesian_rotation &
      & * unique_directions(j)%displacement
      xx(atom_1p%id()) = xx(atom_1p%id()) + outer_product(x,x)
    enddo
  enddo
  
  ! Check that xx can be inverted.
  do i=1,structure%no_atoms
    if (determinant(xx(i))<min_determinant) then
      call print_line(CODE_ERROR//": error generating unique directions; &
         &x'^x' can't be inverted.")
      call err()
    endif
  enddo
end subroutine
end module
