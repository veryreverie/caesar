! ======================================================================
! Symmetries in normal mode co-ordinates.
! ======================================================================
module degenerate_symmetry_module
  use common_module
  implicit none
  
  private
  
  public :: DegenerateSymmetry
  
  ! The mode u_i, and the action of the symmetry on u_i, S.u_i, given in terms
  !    of other modes {u_j}.
  ! S.u_i = sum_j a_j u_j
  ! symmetric_mode_ids gives the set of j for which a_j is non-zero,
  ! and symmetric_mode_coefficients are the corresponding a_j.
  type, extends(Stringable) :: SingleModeSymmetry
    integer                  :: mode_id
    integer,     allocatable :: symmetric_mode_ids(:)
    complex(dp), allocatable :: symmetric_mode_coefficients(:)
  contains
    procedure, public :: read  => read_SingleModeSymmetry
    procedure, public :: write => write_SingleModeSymmetry
  end type
  
  interface SingleModeSymmetry
    module procedure new_SingleModeSymmetry
    module procedure new_SingleModeSymmetry_String
  end interface
  
  ! A mode-by-mode list of SingleModeSymmetry.
  type, extends(Stringsable) :: DegenerateSymmetry
    integer                                        :: symmetry_id
    type(SingleModeSymmetry), allocatable, private :: symmetries_(:)
  contains
    procedure, public :: calculate_symmetry
    
    procedure, public :: transform_monomial
    
    ! I/O.
    procedure, public :: read  => read_DegenerateSymmetry
    procedure, public :: write => write_DegenerateSymmetry
  end type
  
  interface DegenerateSymmetry
    module procedure new_DegenerateSymmetry
    module procedure new_DegenerateSymmetry_Strings
    module procedure new_DegenerateSymmetry_StringArray
  end interface
contains

function new_SingleModeSymmetry(mode_id,symmetric_mode_ids, &
   & symmetric_mode_coefficients) result(this)
  implicit none
  
  integer,     intent(in)  :: mode_id
  integer,     intent(in)  :: symmetric_mode_ids(:)
  complex(dp), intent(in)  :: symmetric_mode_coefficients(:)
  type(SingleModeSymmetry) :: this
  
  if (size(symmetric_mode_ids)/=size(symmetric_mode_coefficients)) then
    call err()
  endif
  
  this%mode_id = mode_id
  this%symmetric_mode_ids = symmetric_mode_ids
  this%symmetric_mode_coefficients = symmetric_mode_coefficients
end function

function new_DegenerateSymmetry(symmetry,subspaces,modes,qpoints) result(this)
  implicit none
  
  type(SymmetryOperator),   intent(in)    :: symmetry
  type(DegenerateSubspace), intent(in)    :: subspaces(:)
  type(ComplexMode),        intent(in)    :: modes(:)
  type(QpointData),         intent(in)    :: qpoints(:)
  type(DegenerateSymmetry)                :: this
  
  type(ComplexMatrix), allocatable :: subspace_symmetries(:)
  
  type(ComplexMode), allocatable :: degenerate_modes(:)
  type(QpointData),  allocatable :: degenerate_qpoints(:)
  
  integer :: max_mode_id
  
  type(SingleModeSymmetry), allocatable :: symmetries(:)
  type(ComplexMode)                     :: mode
  type(QpointData)                      :: qpoint
  type(QpointData)                      :: transformed_qpoint
  integer                               :: subspace_position
  type(DegenerateSubspace)              :: subspace
  type(ComplexMode),        allocatable :: symmetric_modes(:)
  integer                               :: mode_i_position
  integer,                  allocatable :: symmetric_mode_positions(:)
  complex(dp),              allocatable :: symmetric_mode_coefficients(:)
  complex(dp),              allocatable :: subspace_symmetry(:,:)
  
  integer :: i,j,ialloc
  
  this%symmetry_id = symmetry%id
  
  ! Calculate symmetries subspace by subspace.
  allocate(subspace_symmetries(size(subspaces)), stat=ialloc); call err(ialloc)
  do i=1,size(subspaces)
    degenerate_modes = subspaces(i)%modes(modes)
    degenerate_qpoints = subspaces(i)%qpoints(modes,qpoints)
    subspace_symmetries(i) =                                           &
       & calculate_symmetry_in_normal_coordinates( degenerate_modes,   &
       &                                           degenerate_qpoints, &
       &                                           symmetry            )
  enddo
  
  allocate(this%symmetries_(size(modes)), stat=ialloc); call err(ialloc)
  do i=1,size(this%symmetries_)
    mode = modes(i)
    
    ! Identify the q-point of mode i, and the q-point of the set of modes
    !    which the symmetry transfoms mode i to.
    qpoint = qpoints(first(qpoints%id==mode%qpoint_id))
    transformed_qpoint = symmetry*qpoint
    transformed_qpoint = qpoints(first(qpoints==transformed_qpoint))
    
    ! Identify the subspace to which mode i belongs, and its position in
    !    the array subspaces.
    subspace_position = first(subspaces%id==modes(i)%subspace_id)
    subspace = subspaces(subspace_position)
    
    ! Identify the set of modes which the symmetry transforms mode i to.
    ! These modes must be degenerate with mode i, and at the q-point which
    !    the symmetry transforms mode i's q-point to.
    symmetric_modes = subspace%modes(modes)
    symmetric_modes = symmetric_modes(                            &
       & filter(symmetric_modes%qpoint_id==transformed_qpoint%id) )
    
    ! Identify the position of mode i within the subspace's mode_ids,
    !    and the positions of the symmetric modes within the same.
    mode_i_position = first(subspace%mode_ids==mode%id)
    symmetric_mode_positions =                               &
       & [( first(subspace%mode_ids==symmetric_modes(j)%id), &
       & j=1,                                                &
       & size(symmetric_modes)                               )]
    
    ! Identify the coefficients of the symmetry.
    ! If mode i is ui, and symmetric_modes are {uj}, then these are the
    !    coefficients of S.ui in terms of uj.
    subspace_symmetry = cmplx(subspace_symmetries(subspace_position))
    symmetric_mode_coefficients = subspace_symmetry( &
                        & symmetric_mode_positions,  &
                        & mode_i_position            )
    
    ! Construct the symmetry.
    this%symmetries_(i) = SingleModeSymmetry( mode%id,                    &
                                            & symmetric_modes%id,         &
                                            & symmetric_mode_coefficients )
  enddo
end function

! Calculate the symmetry in a basis of complex monomials.
! If include_symmetry is .false. then the symmetry will not include the
!    coefficients of the basis monomials.
! If include_symmetry is .true., then the coefficients will be included.
function calculate_symmetry(this,input,modes,include_coefficients) &
   & result(output)
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
  
  type(ComplexPolynomial) :: transformed_input
  
  if (size(input)==0) then
    output = cmplxmat(zeroes(0,0))
    return
  endif
  
  allocate(symmetry(size(input), size(input)), stat=ialloc); call err(ialloc)
  do i=1,size(input)
    transformed_input = this%transform_monomial(input(i), modes)
    do j=1,size(input)
      k = first_equivalent( transformed_input%terms,   &
                          & input(j),                  &
                          & compare_complex_monomials, &
                          & default = 0                )
      if (k==0) then
        symmetry(j,i) = 0
      else
        symmetry(j,i) = transformed_input%terms(k)%coefficient
      endif
      
      if (include_coefficients) then
        symmetry(j,i) = symmetry(j,i)        &
                    & * input(i)%coefficient &
                    & / input(j)%coefficient
      endif
    enddo
  enddo
  
  output = symmetry
end function

function transform_monomial(this,input,modes) result(output)
  implicit none
  
  class(DegenerateSymmetry), intent(in) :: this
  type(ComplexMonomial),     intent(in) :: input
  type(ComplexMode),         intent(in) :: modes(:)
  type(ComplexPolynomial)               :: output
  
  integer,                  allocatable :: mode_ids(:)
  type(SingleModeSymmetry), allocatable :: symmetry
  type(ComplexMode)                     :: symmetric_mode
  type(ComplexMonomial),    allocatable :: monomials(:)
  
  integer :: i,j,k,ialloc
  
  ! Convert the monomial into a list of modes.
  allocate(mode_ids(input%total_power()), stat=ialloc); call err(ialloc)
  k = 0
  do i=1,size(input)
    do j=1,input%power(i)
      k = k+1
      mode_ids(k) = input%id(i)
    enddo
    
    if (input%id(i)/=input%paired_id(i)) then
      do j=1,input%paired_power(i)
        k = k+1
        mode_ids(k) = input%paired_id(i)
      enddo
    endif
  enddo
  
  do i=1,size(mode_ids)
    symmetry = this%symmetries_(first(this%symmetries_%mode_id==mode_ids(i)))
    allocate( monomials(size(symmetry%symmetric_mode_ids)), &
            & stat=ialloc); call err(ialloc)
    do j=1,size(symmetry%symmetric_mode_ids)
      symmetric_mode = modes(first(modes%id==symmetry%symmetric_mode_ids(j)))
      monomials(j) = ComplexMonomial(                                 &
         & coefficient = symmetry%symmetric_mode_coefficients(j),     &
         & modes       = [ComplexUnivariate(symmetric_mode, power=1)] )
    enddo
    
    if (i==1) then
      output = ComplexPolynomial(monomials)
    else
      output%terms = [(                                          &
         & (output%terms(i)*monomials(j), j=1, size(monomials)), &
         & i=1,                                                  &
         & size(output)                                          )]
    endif
    deallocate(monomials, stat=ialloc); call err(ialloc)
  enddo
  
  call output%simplify()
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_SingleModeSymmetry(this,input)
  implicit none
  
  class(SingleModeSymmetry), intent(out) :: this
  type(String),              intent(in)  :: input
  
  type(String), allocatable :: line(:)
  type(String), allocatable :: symmetric_mode(:)
  
  integer                  :: mode_id
  integer,     allocatable :: symmetric_mode_ids(:)
  complex(dp), allocatable :: symmetric_mode_coefficients(:)
  
  integer :: i,ialloc
  
  select type(this); type is(SingleModeSymmetry)
    ! If the symmetry takes u1 to a*u2+b*u3+... then
    ! Input = 'u1 -> au2 + bu3 + ...'
    ! line = ['u1', '->', 'au2', '+', 'bu3', ...]
    line = split_line(input)
    
    ! Split off the ID of the input mode. (trim the 'u' off the front.)
    mode_id = int(slice(line(1),2,len(line(1))))
    
    ! Split each 'au2', 'bu3' etc. term by the 'u', to leave the coefficient
    !    and mode id of each.
    allocate( symmetric_mode_ids((size(line)-1)/2),          &
            & symmetric_mode_coefficients((size(line)-1)/2), &
            & stat=ialloc); call err(ialloc)
    do i=1,size(symmetric_mode_ids)
      symmetric_mode = split_line(line(2*i+1), delimiter='u')
      symmetric_mode_ids(i) = int(symmetric_mode(2))
      symmetric_mode_coefficients(i) = cmplx(symmetric_mode(1))
    enddo
    
    this = SingleModeSymmetry( mode_id,                    &
                             & symmetric_mode_ids,         &
                             & symmetric_mode_coefficients )
  class default
    call err()
  end select
end subroutine

function write_SingleModeSymmetry(this) result(output)
  implicit none
  
  class(SingleModeSymmetry), intent(in) :: this
  type(String)                          :: output
  
  type(String), allocatable :: symmetric_modes(:)
  
  integer :: i,ialloc
  
  select type(this); type is(SingleModeSymmetry)
    allocate( symmetric_modes(size(this%symmetric_mode_ids)), &
            & stat=ialloc); call err(ialloc)
    do i=1,size(symmetric_modes)
      symmetric_modes(i) = this%symmetric_mode_coefficients(i)// &
                         & 'u'//this%symmetric_mode_ids(i)
    enddo
    output = 'u'//this%mode_id//' -> '//join(symmetric_modes,delimiter=' + ')
  class default
    call err()
  end select
end function

impure elemental function new_SingleModeSymmetry_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(SingleModeSymmetry) :: this
  
  call this%read(input)
end function

subroutine read_DegenerateSymmetry(this,input)
  implicit none
  
  class(DegenerateSymmetry), intent(out) :: this
  type(String),              intent(in)  :: input(:)
  
  type(String), allocatable :: line(:)
  
  select type(this); type is(DegenerateSymmetry)
    line = split_line(input(1))
    this%symmetry_id = int(line(3))
    this%symmetries_ = SingleModeSymmetry(input(2:))
  class default
    call err()
  end select
end subroutine

function write_DegenerateSymmetry(this) result(output)
  implicit none
  
  class(DegenerateSymmetry), intent(in) :: this
  type(String), allocatable             :: output(:)
  
  select type(this); type is(DegenerateSymmetry)
    output = [ 'Symmetry : '//this%symmetry_id, &
             & str(this%symmetries_)            ]
  class default
    call err()
  end select
end function

function new_DegenerateSymmetry_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(DegenerateSymmetry) :: this
  
  call this%read(input)
end function

impure elemental function new_DegenerateSymmetry_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(DegenerateSymmetry)      :: this
  
  this = DegenerateSymmetry(str(input))
end function
end module
