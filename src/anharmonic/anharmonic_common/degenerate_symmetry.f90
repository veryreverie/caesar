! ======================================================================
! Symmetries in normal mode co-ordinates.
! ======================================================================
! Since symmetries only map between modes in the same degenerate subspace,
!    the symmetry is stored as an array of matrices, one for each subspace.
! The ids of the modes corresponding to the elements of the symmetry are
!    stored separately.
! e.g. if there are two subspaces, with IDs 4 and 6, where subspace 4 contains
!    modes with IDs 1,3 and 7, and subspace 6 contains modes with IDs 5 and 8,
!    then subspace_ids_ = [4,6], and mode_ids_ = [[1,3,7],[5,8]],
!    and single_mode_symmetries contains two matrices, the first being a 3x3
!    matrix with elements corresponding to modes 1, 3 and 7, and the second
!    being a 2x2 matrix with elements corresponding to modes 5 and 8.
module degenerate_symmetry_module
  use common_module
  
  use mode_monomial_module
  implicit none
  
  private
  
  public :: DegenerateSymmetry
  
  type :: DegenerateSymmetry
    integer,             allocatable, private :: subspace_ids_(:)
    type(IntArray1D),    allocatable, private :: mode_ids_(:)
    type(ComplexMatrix), allocatable, private :: single_mode_symmetries_(:)
  contains
    generic,   public  :: calculate_symmetry =>             &
                        & calculate_symmetry_ModeMonomials, &
                        & calculate_symmetry_ComplexMonomials
    procedure, private :: calculate_symmetry_ModeMonomials
    procedure, private :: calculate_symmetry_ComplexMonomials
  end type
  
  interface DegenerateSymmetry
    module procedure new_DegenerateSymmetry
  end interface
contains

function new_DegenerateSymmetry(symmetry,degenerate_subspaces,modes,qpoints, &
   & logfile) result(output)
  implicit none
  
  type(SymmetryOperator),   intent(in)    :: symmetry
  type(DegenerateSubspace), intent(in)    :: degenerate_subspaces(:)
  type(ComplexMode),        intent(in)    :: modes(:)
  type(QpointData),         intent(in)    :: qpoints(:)
  type(OFile),              intent(inout) :: logfile
  type(DegenerateSymmetry)                :: output
  
  type(ComplexMode), allocatable :: degenerate_modes(:)
  type(QpointData),  allocatable :: degenerate_qpoints(:)
  
  integer :: max_mode_id
  
  integer :: i,j,ialloc
  
  do i=1,size(degenerate_subspaces)
    max_mode_id = max(max_mode_id, maxval(degenerate_subspaces(i)%mode_ids))
  enddo
  
  allocate( output%subspace_ids_(size(degenerate_subspaces)),           &
          & output%mode_ids_(size(degenerate_subspaces)),               &
          & output%single_mode_symmetries_(size(degenerate_subspaces)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(degenerate_subspaces)
    degenerate_modes = degenerate_subspaces(i)%modes(modes)
    degenerate_qpoints = degenerate_subspaces(i)%qpoints(modes,qpoints)
    output%single_mode_symmetries_(i) =                                &
       & calculate_symmetry_in_normal_coordinates( degenerate_modes,   &
       &                                           degenerate_qpoints, &
       &                                           symmetry,           &
       &                                           logfile)
    output%subspace_ids_(i) = degenerate_subspaces(i)%id
    output%mode_ids_(i) = degenerate_subspaces(i)%mode_ids
  enddo
end function

! Calculate the symmetry in a basis of mode monomials.
function calculate_symmetry_ModeMonomials(this,input,modes) result(output)
  implicit none
  
  class(DegenerateSymmetry), intent(in) :: this
  type(ModeMonomial),        intent(in) :: input(:)
  type(ComplexMode),         intent(in) :: modes(:)
  type(ComplexMatrix)                   :: output
  
  type(ComplexMode)             :: mode
  integer                       :: subspace_position
  type(IntArray1D), allocatable :: subspace_positions(:)
  type(IntArray1D), allocatable :: mode_positions(:)
  
  complex(dp), allocatable :: matrix(:,:)
  complex(dp), allocatable :: single_mode_symmetry(:,:)
  
  integer :: i,j,k,ialloc
  
  if (size(input)==0) then
    output = cmplxmat(zeroes(0,0))
    return
  endif
  
  ! Locate the position of the element corresponding to each mode in each
  !    monomial.
  allocate( subspace_positions(size(input)), &
          & mode_positions(size(input)),     &
          & stat=ialloc); call err(ialloc)
  do i=1,size(input)
    allocate( subspace_positions(i)%i(size(input(i))), &
            & mode_positions(i)%i(size(input(i))),     &
            & stat=ialloc); call err(ialloc)
    do j=1,size(input(i))
      mode = modes(first(modes%id==input(i)%mode_ids(j)))
      subspace_position = first(this%subspace_ids_==mode%subspace_id)
      subspace_positions(i)%i(j) = subspace_position
      mode_positions(i)%i(j) = &
         & first(this%mode_ids_(subspace_position)%i==mode%id)
    enddo
  enddo
  
  ! Construct the symmetry as the product of one-mode symmetries.
  allocate(matrix(size(input), size(input)), stat=ialloc); call err(ialloc)
  do i=1,size(input)
    do j=1,size(input)
      if (size(input(i))/=size(input(j))) then
        matrix(j,i) = 0
      elseif (subspace_positions(i)/=subspace_positions(j)) then
        matrix(j,i) = 0
      else
        matrix(j,i) = 1
        do k=1,size(input(i))
          single_mode_symmetry = &
             & cmplx(this%single_mode_symmetries_(subspace_positions(i)%i(j)))
          matrix(j,i) = matrix(j,i)                                   &
                    & * single_mode_symmetry( mode_positions(j)%i(k), &
                    &                         mode_positions(i)%i(k)  )
        enddo
      endif
    enddo
  enddo
  
  output = matrix
end function

! Calculate the symmetry in a basis of complex monomials.
! If include_symmetry is .false. then the symmetry will not include the
!    coefficients of the basis monomials.
! If include_symmetry is .true., then the coefficients will be included.
function calculate_symmetry_ComplexMonomials(this,input,modes, &
   & include_coefficients) result(output)
  implicit none
  
  class(DegenerateSymmetry), intent(in) :: this
  type(ComplexMonomial),     intent(in) :: input(:)
  type(ComplexMode),         intent(in) :: modes(:)
  logical,                   intent(in) :: include_coefficients
  type(ComplexMatrix)                   :: output
  
  integer              :: max_mode_id
  integer, allocatable :: subspace_pos(:)
  integer, allocatable :: mode_pos(:)
  
  type(IntArray2D), allocatable :: permutations(:)
  
  integer, allocatable :: mode_ids_i(:)
  integer, allocatable :: mode_ids_j(:)
  integer, allocatable :: subspace_pos_i(:)
  integer, allocatable :: subspace_pos_j(:)
  integer, allocatable :: mode_pos_i(:)
  integer, allocatable :: mode_pos_j(:)
  
  complex(dp), allocatable :: single_mode_symmetry(:,:)
  complex(dp), allocatable :: symmetry(:,:)
  complex(dp), allocatable :: element
  
  integer :: i,j,k,l,ialloc
  
  if (size(input)==0) then
    output = cmplxmat(zeroes(0,0))
    return
  endif
  
  ! For each mode, locate the corresponding symmetry matrix and the element
  !    within that matrix.
  max_mode_id = maxval(modes%id)
  allocate( subspace_pos(max_mode_id), &
          & mode_pos(max_mode_id),     &
          & stat=ialloc); call err(ialloc)
  do i=1,max_mode_id
    j = first(modes%id==i, default=0)
    if (j==0) then
      subspace_pos(i) = 0
      mode_pos(i) = 0
    else
      subspace_pos(i) = first( this%subspace_ids_==modes(j)%subspace_id, &
                             & default=0)
      if (subspace_pos(i)==0) then
        mode_pos(i) = 0
      else
        mode_pos(i) = first( this%mode_ids_(subspace_pos(i))%i==modes(j)%id, &
                           & default=0)
      endif
    endif
  enddo
  
  ! TODO: This can be accelerated by considering q-points.
  ! Symmetry S takes q-points q_i to q-points Sq_i.
  ! If the q-points corresponding to the modes of the first permutation of a
  !    monomial are [q_1, q_2, ... ], then the only permutations which need
  !    considering have q-points [Sq_1, Sq_2, ...].
  ! END TODO
  
  ! Identify all permutations of the modes in each monomial.
  permutations = generate_permutations(input)
  
  ! Construct the symmetry as the product of one-mode symmetries.
  allocate(symmetry(size(input), size(input)), stat=ialloc); call err(ialloc)
  symmetry = 0
  do i=1,size(input)
    do j=1,size(input)
      ! Average over all permutations of monomial i.
      mode_ids_j = permutations(j)%i(1)%i
      do k=1,size(permutations(i))
        mode_ids_i = permutations(i)%i(k)%i
        subspace_pos_i = subspace_pos(mode_ids_i)
        subspace_pos_j = subspace_pos(mode_ids_j)
        mode_pos_i = mode_pos(mode_ids_i)
        mode_pos_j = mode_pos(mode_ids_j)
        if (any(mode_pos_i==0) .or. any(mode_pos_j==0)) then
          call print_line(ERROR//': Mode not present in symmetry.')
          call err()
        elseif (size(subspace_pos_i)/=size(subspace_pos_j)) then
          cycle
        elseif (any(subspace_pos_i/=subspace_pos_j)) then
          cycle
        else
          element = 1
          do l=1,size(mode_ids_i)
            single_mode_symmetry = &
               & cmplx(this%single_mode_symmetries_(subspace_pos_i(l)))
            element = element                        &
                  & * single_mode_symmetry( mode_pos_j(l), &
                  &                         mode_pos_i(l)  )
          enddo
          symmetry(j,i) = symmetry(j,i) + element
        endif
      enddo
      ! TODO: Is this right?
      symmetry(j,i) = symmetry(j,i)         &
                  & * size(permutations(i)) &
                  & / size(permutations(j))
      if (include_coefficients) then
        symmetry(j,i) = symmetry(j,i)        &
                    & * input(j)%coefficient &
                    & / input(i)%coefficient
      endif
    enddo
  enddo
  
  output = symmetry
end function

impure elemental function generate_permutations(input) result(output)
  implicit none
  
  type(ComplexMonomial), intent(in) :: input
  type(IntArray2D)                  :: output
  
  type(IntArray1D) :: root
  
  root = array([integer::])
  
  output%i = generate_permutations_helper(input,root)
end function

recursive function generate_permutations_helper(input,root) result(output)
  implicit none
  
  type(ComplexMonomial), intent(in) :: input
  type(IntArray1D),      intent(in) :: root
  type(IntArray1D), allocatable     :: output(:)
  
  type(ComplexMonomial) :: new_input
  type(IntArray1D)      :: new_root
  
  integer :: i
  
  if (input%total_power()==0) then
    output = [root]
    return
  endif
  
  output = [IntArray1D::]
  
  do i=1,size(input)
    if (input%modes(i)%id==input%modes(i)%paired_id) then
      if (input%modes(i)%power>0) then
        new_input = input
        new_input%modes(i)%power = new_input%modes(i)%power - 1
        new_input%modes(i)%paired_power = new_input%modes(i)%paired_power - 1
        new_root = [root%i, new_input%modes(i)%id]
        output = [output, generate_permutations_helper(new_input,new_root)]
      endif
    else
      if (input%modes(i)%power>0) then
        new_input = input
        new_input%modes(i)%power = new_input%modes(i)%power - 1
        new_root = [root%i, new_input%modes(i)%id]
        output = [output, generate_permutations_helper(new_input,new_root)]
      endif
      
      if (input%modes(i)%paired_power>0) then
        new_input = input
        new_input%modes(i)%paired_power = new_input%modes(i)%paired_power - 1
        new_root = [root%i, new_input%modes(i)%paired_id]
        output = [output, generate_permutations_helper(new_input,new_root)]
      endif
    endif
  enddo
end function
end module
