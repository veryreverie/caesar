module degenerate_symmetry_module
  use common_module
  
  use mode_monomial_module
  implicit none
  
  private
  
  public :: DegenerateSymmetry
  
  type :: DegenerateSymmetry
    type(ComplexMatrix), allocatable, private :: single_mode_symmetries_(:)
    integer,             allocatable, private :: subspace_ids_(:)
    type(IntArray1D),    allocatable, private :: mode_ids_(:)
    integer,             allocatable, private :: subspace_pos_(:)
    integer,             allocatable, private :: mode_pos_(:)
  contains
    procedure, public :: calculate_symmetry
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
  
  allocate( output%single_mode_symmetries_(size(degenerate_subspaces)), &
          & output%subspace_ids_(size(degenerate_subspaces)),           &
          & output%mode_ids_(size(degenerate_subspaces)),               &
          & output%subspace_pos_(max_mode_id),                          &
          & output%mode_pos_(max_mode_id),                              &
          & stat=ialloc); call err(ialloc)
  output%subspace_pos_ = 0
  output%mode_pos_ = 0
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
    do j=1,size(output%mode_ids_(i))
      output%subspace_pos_(output%mode_ids_(i)%i(j)) = i
      output%mode_pos_(output%mode_ids_(i)%i(j)) = j
    enddo
  enddo
end function

function calculate_symmetry(this,mode_monomials) result(output)
  implicit none
  
  class(DegenerateSymmetry), intent(in) :: this
  type(ModeMonomial),        intent(in) :: mode_monomials(:)
  type(ComplexMatrix)                   :: output
  
  complex(dp), allocatable :: matrix(:,:)
  complex(dp), allocatable :: single_mode_symmetry(:,:)
  
  integer :: no_mode_monomials
  
  type(IntArray1D), allocatable :: mode_ids(:)
  type(IntArray1D), allocatable :: subspace_positions(:)
  type(IntArray1D), allocatable :: mode_positions(:)
  integer,          allocatable :: sort_key(:)
  
  integer :: i,j,k,ialloc
  
  if (size(mode_monomials)==0) then
    output = cmplxmat(zeroes(0,0))
    return
  endif
  
  no_mode_monomials = size(mode_monomials(1))
  do i=1,size(mode_monomials)
    if (size(mode_monomials(i))/=no_mode_monomials) then
      call err()
    endif
  enddo
  
  allocate( mode_ids(size(mode_monomials)),                   &
          & subspace_positions(size(mode_monomials)),         &
          & mode_positions(size(mode_monomials)),             &
          & stat=ialloc); call err(ialloc)
  do i=1,size(mode_monomials)
    mode_ids(i) = mode_monomials(i)%ids
    subspace_positions(i) = this%subspace_pos_(mode_ids(i)%i)
    mode_positions(i) = this%mode_pos_(mode_ids(i)%i)
    
    sort_key = sort(subspace_positions(i)%i)
    
    mode_ids(i) = mode_ids(i)%i(sort_key)
    subspace_positions(i) = subspace_positions(i)%i(sort_key)
    mode_positions(i) = mode_positions(i)%i(sort_key)
  enddo
  
  do i=2,size(mode_monomials)
    if (subspace_positions(i)/=subspace_positions(1)) then
      call err()
    endif
  enddo
  
  allocate( matrix(size(mode_monomials),size(mode_monomials)), &
          & stat=ialloc); call err(ialloc)
  matrix = 1
  do i=1,no_mode_monomials
    single_mode_symmetry = &
       & cmplx(this%single_mode_symmetries_(subspace_positions(1)%i(i)))
    do j=1,size(mode_monomials)
      do k=1,size(mode_monomials)
        matrix(k,j) = matrix(k,j)                                   &
                  & * single_mode_symmetry( mode_positions(k)%i(i), &
                                          & mode_positions(j)%i(i))
      enddo
    enddo
  enddo
  
  output = matrix
end function
end module
