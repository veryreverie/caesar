module lift_degeneracies_module
  use constants_module, only : dp
  use string_module
  use io_module
  implicit none
contains

function lift_degeneracies(input,structure) result(output)
  use normal_mode_module
  use structure_module
  use linear_algebra_module
  use group_module
  use symmetry_module
  implicit none
  
  type(ComplexMode),   intent(in) :: input(:)
  type(StructureData), intent(in) :: structure
  type(ComplexMode), allocatable  :: output(:)
  
  type(SymmetryOperator), allocatable :: relevant_symmetries(:)
  type(ComplexMatrix), allocatable :: symmetries_nm_basis(:)
  
  output = input
  
  !relevant_symmetries = identify_relevant_symmetries(supercell%symmetries)
  
  !symmetries_nm_basis = calculate_normal_mode_symmetries( &
  !                                 & relevant_symmetries, &
  !                                 & normal_modes)
end function

! ----------------------------------------------------------------------
! Calculates which symmetries commute with translation by all G-vectors,
!    but are not themselves a translation by a lattice vector.
! ----------------------------------------------------------------------
function identify_relevant_symmetries(symmetries) result(output)
  use linear_algebra_module
  use structure_module
  use group_module
  use symmetry_module
  implicit none
  
  type(SymmetryOperator), intent(in)  :: symmetries(:)
  type(SymmetryOperator), allocatable :: output(:)
  
  logical, allocatable :: translations(:)
  logical, allocatable :: relevant(:)
  
  type(IntMatrix) :: identity_matrix
  type(Group)     :: identity_group
  
  integer :: i,j,ialloc
  
  identity_matrix = make_identity_matrix(3)
  identity_group = make_identity_group(size(symmetries(1)%atom_group))
  
  ! Identify purely translational symmetries.
  allocate(translations(size(symmetries)), stat=ialloc); call err(ialloc)
  translations = .false.
  do i=1,size(symmetries)
    if (symmetries(i)%rotation==identity_matrix) then
      translations(i) = .true.
      exit
    endif
  enddo
  
  ! Identify all relevant symmetries.
  allocate(relevant(size(symmetries)), stat=ialloc); call err(ialloc)
  relevant = .false.
  do i=1,size(symmetries)
    ! Ignore pure translations.
    if (translations(i)) then
      cycle
    endif
    
    ! Check if the symmetry commutes with each translation.
    ! n.b. the rotations will always commute, since R_j=I -> R_i*R_j=R_j*R_i.
    ! Thus only atom mappings need to be checked.
    do j=1,size(symmetries)
      if (translations(j)) then
        if (  symmetries(i)%atom_group*symmetries(j)%atom_group &
         & == symmetries(j)%atom_group*symmetries(i)%atom_group ) then
          relevant(i) = .true.
        endif
      endif
    enddo
  enddo
  
  ! Produce output.
  allocate( output(count(relevant)), stat=ialloc); call err(ialloc)
  j = 0
  do i=1,size(symmetries)
    if (relevant(i)) then
      j = j+1
      output(j) = symmetries(i)
    endif
  enddo
  
  if (j/=size(output)) then
    call err()
  endif
end function

! ----------------------------------------------------------------------
! Calculates the eigenvectors of a symmetry operation.
! ----------------------------------------------------------------------
function calculate_eigenvectors(rotation,symmetry_group) result(output)
  use linear_algebra_module
  use group_module
  implicit none
  
  type(IntMatrix), intent(in) :: rotation
  type(Group),     intent(in) :: symmetry_group
  integer                     :: output
  
  integer              :: order
  integer, allocatable :: phases(:)
  
  integer :: i,ialloc
  
  ! order = n
  order = calculate_symmetry_order(rotation,symmetry_group)
  
  ! phase=j corresponds to e^(2*pi*i*j/n)
  allocate(phases(order), stat=ialloc); call err(ialloc)
  do i=1,order
    phases(i) = i-1
  enddo
end function

! ----------------------------------------------------------------------
! Calculates the order of each symmetry operation.
! The order is the smallest integer i>0 s.t. x^i=I, where I is the identity.
! ----------------------------------------------------------------------
function calculate_symmetry_order(rotation,symmetry_group) &
   & result(output)
  use linear_algebra_module
  use group_module
  implicit none
  
  type(IntMatrix), intent(in) :: rotation
  type(Group),     intent(in) :: symmetry_group
  integer                     :: output
  
  type(IntMatrix) :: identity_matrix
  type(Group)     :: identity_group
  
  type(IntMatrix) :: power_matrix
  type(Group)     :: power_group
  
  integer :: i,j,ialloc
  
  identity_matrix = make_identity_matrix(3)
  identity_group = make_identity_group(size(symmetry_group))
  
  power_matrix = identity_matrix
  power_group = identity_group
  
  do i=1,huge(0)
    power_matrix = rotation * power_matrix
    power_group = symmetry_group * power_group
    
    if (power_matrix==identity_matrix .and. power_group==identity_group) then
      output = i
      exit
    endif
  enddo
end function

! ----------------------------------------------------------------------
! Finds an irreducible basis of symmetries from which the rest can be
!    constructed.
! ----------------------------------------------------------------------
! Returns an array of the IDs of the relevant symmetries.
function calculate_irreducible_symmetries(structure,operator_symmetry_group, &
   &orders) result(output)
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
    if ( structure%symmetries(op)%rotation==make_identity_matrix(3) .and. &
       & orders(op)>1) then
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
