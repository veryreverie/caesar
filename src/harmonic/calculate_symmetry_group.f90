! ======================================================================
! Calculates various properties of a structure's symmetries.
! ======================================================================
module calculate_symmetry_group_module
  use constants_module, only : dp
  use string_module
  use io_module
contains

! ----------------------------------------------------------------------
! Calculates how symmetries map atoms onto other atoms.
! ----------------------------------------------------------------------
! If symmetry i maps atom j to atom k then output(i)*j=k.
function calculate_atom_symmetry_group(structure) result(output)
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
  do i=1,structure%no_atoms
    atom_pos_frac(i) = structure%recip_lattice * structure%atoms(i)
  enddo
  
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

! ----------------------------------------------------------------------
! Calculates how symmetries map symmetries onto other symmetries.
! ----------------------------------------------------------------------
! If symmetry i * symmetry j = symmetry k then output(i)*j=k.
function calculate_operator_symmetry_group(structure,atom_symmetry_group) &
   & result(output)
  use structure_module
  use group_module
  use linear_algebra_module
  implicit none
  
  type(StructureData), intent(in) :: structure
  type(Group),         intent(in) :: atom_symmetry_group(:)
  type(Group), allocatable        :: output(:)
  
  type(IntMatrix) :: rotation_k
  type(Group)     :: operation_k
  
  ! Temporary variables.
  integer              :: i,j,k,ialloc
  integer, allocatable :: temp_operator_group(:)
  
  allocate( temp_operator_group(structure%no_symmetries), &
          & output(structure%no_symmetries),              &
          & stat=ialloc); call err(ialloc)
  do i=1,structure%no_symmetries
    do_j : do j=1,structure%no_symmetries
      rotation_k = structure%rotations(i)*structure%rotations(j)
      operation_k = atom_symmetry_group(i)*atom_symmetry_group(j)
      do k=1,structure%no_symmetries
        if ( rotation_k==structure%rotations(k) .and. &
           & operation_k==atom_symmetry_group(k)) then
          temp_operator_group(j) = k
          cycle do_j
        endif
      enddo
      
      call print_line('Error: symmetry '//i//' times symmetry '//j//' is not &
         &itself a symmetry.')
      call err()
    enddo do_j
    
    output(i) = temp_operator_group
  enddo
end function

! ----------------------------------------------------------------------
! Calculates the order of each symmetry operation.
!    (if x^n=I, then the order of x is n.)
! ----------------------------------------------------------------------
function calculate_symmetry_orders(structure,operator_symmetry_group) &
   & result(output)
  use structure_module
  use group_module
  implicit none
  
  type(StructureData), intent(in) :: structure
  type(Group), intent(in)         :: operator_symmetry_group(:)
  integer, allocatable            :: output(:)
  
  integer :: temp
  
  integer :: i,j,ialloc
  
  allocate( output(size(operator_symmetry_group)), &
          & stat=ialloc); call err(ialloc)
  
  ! Loop over groups, identifying the order of each.
  ! This uses: if x^(n+1)=x then x^n=I  =>  the order of x is n.
  do_i : do i=1,size(operator_symmetry_group)
    temp = i
    do j=1,max(6,structure%sc_size)
      temp = operator_symmetry_group(i) * temp
      if (temp==i) then
        output(i) = j
        cycle do_i
      endif
    enddo
    
    call print_line('Error: Symmetry order too high.')
    call err()
  enddo do_i
end function

! ----------------------------------------------------------------------
! Finds an irreducible basis of symmetries from which the rest can be
!    constructed.
! ----------------------------------------------------------------------
! Returns an array of the IDs of the relevant symmetries.
function calculate_irreducible_symmetries(structure,operator_symmetry_group, &
   &orders) result(output)
  use constants_module, only : identity
  use structure_module
  use group_module
  use linear_algebra_module
  implicit none
  
  type(StructureData), intent(in) :: structure
  type(Group),         intent(in) :: operator_symmetry_group(:)
  integer,             intent(in) :: orders(:)
  integer, allocatable            :: output(:)
  
  integer              :: no_symmetries
  
  ! Base variables.
  integer, allocatable :: base(:)
  integer, allocatable :: base_operators(:)
  logical, allocatable :: is_translation(:)
  
  ! Sort variables.
  integer, allocatable :: sort_values(:)
  integer, allocatable :: sort_base(:)
  integer, allocatable :: sort(:)
  
  logical, allocatable :: is_irreducible(:)
  
  ! Temporary variables.
  integer              :: i,j,k,ialloc
  
  integer              :: temp
  integer              :: total
  integer              :: op
  
  no_symmetries = size(operator_symmetry_group)
  if (no_symmetries /= size(orders)) then
    call print_line('Error: arrays of different sizes passed to symmetry &
       &finder.')
    call err()
  endif
  
  ! Identify which operators are simply powers of another operator.
  allocate( is_irreducible(no_symmetries), &
          & stat=ialloc); call err(ialloc)
  is_irreducible = .true.
  allocate(base(no_symmetries), stat=ialloc); call err(ialloc)
  do i=1,size(base)
    base(i) = i
  enddo
  do i=1,no_symmetries
    ! Ignore operators which are already powers of another operator.
    if (base(i)/=i) then
      cycle
    endif
    
    ! Find the operators which are powers of operator i.
    temp = i
    do j=2,orders(i)-1
      temp = operator_symmetry_group(i) * temp
      if (orders(base(temp)) >= orders(i)) then
        base(temp) = i
        is_irreducible(temp) = .false.
      endif
    enddo
  enddo
  
  total = 0
  do i=1,no_symmetries
    if (base(i)==i) then
      total = total+1
    endif
  enddo
  allocate(base_operators(total), stat=ialloc); call err(ialloc)
  total = 0
  do i=1,no_symmetries
    if (base(i)==i) then
      total = total+1
      base_operators(total) = i
    endif
  enddo
  
  ! Identify the pure translations.
  allocate(is_translation(size(base_operators)), stat=ialloc); call err(ialloc)
  do i=1,size(base_operators)
    op = base_operators(i)
    if (structure%rotations(op)==mat(identity) .and. orders(op)>1) then
      is_translation(i) = .true.
    else
      is_translation(i) = .false.
    endif
  enddo
  
  ! Order base operations in terms of:
  !    Pure translations first.
  !    Then in order of descending "order" (6-fold,4-fold,3-fold,2-fold,I).
  allocate(sort_values(size(base_operators)), stat=ialloc); call err(ialloc)
  do i=1,size(base_operators)
    op = base_operators(i)
    if (is_translation(i)) then
      sort_values(i) = 0
    else
      sort_values(i) = orders(op)
    endif
  enddo
  
  allocate(sort_base(size(base_operators)), stat=ialloc); call err(ialloc)
  sort_base = 0
  do i=1,size(base_operators)
    j = maxloc(sort_values,1,sort_base==0)
    sort_base(j) = i
  enddo
  
  allocate(sort(no_symmetries), stat=ialloc); call err(ialloc)
  sort = 0
  do i=1,size(base_operators)
    sort(base_operators(i)) = sort_base(i)
  enddo
  
  ! Identify a minimal set of operators.
  do i=1,no_symmetries
    do j=1,no_symmetries
      k = operator_symmetry_group(j) * i
      
      if (sort(base(i))>sort(base(k)) .and. sort(base(j))>sort(base(k))) then
        is_irreducible(k) = .false.
      endif
    enddo
  enddo
  
  allocate(output(count(is_irreducible)), stat=ialloc); call err(ialloc)
  j = 0
  do i=1,size(is_irreducible)
    if (is_irreducible(i)) then
      j = j+1
      output(j) = i
    endif
  enddo
  
  ! Check that the number of operations is correct.
  total = 1
  do i=1,size(output)
    total = total*orders(output(i))
  enddo
  if (total/=no_symmetries) then
    call print_line(base_operators)
    call print_line(sort)
    call print_line(total)
    call print_line(no_symmetries)
    call print_line('Error:')
    !call err()
  endif
end function
end module
