! ======================================================================
! A vector along a single real mode.
! ======================================================================
module real_single_mode_vector_submodule
  use utils_module
  
  use structure_module
  
  use cartesian_vector_submodule
  use real_mode_submodule
  implicit none
  
  private
  
  public :: RealSingleModeVector
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: operator(-)
  
  type, extends(Stringable) :: RealSingleModeVector
    ! The id of the mode.
    integer :: id
    
    ! The magnitude of vector along the mode.
    real(dp) :: magnitude
  contains
    ! Convert to cartesian co-ordinates.
    procedure, public :: cartesian_vector => &
       & cartesian_vector_RealSingleModeVector
    
    ! I/O.
    procedure, public :: read  => read_RealSingleModeVector
    procedure, public :: write => write_RealSingleModeVector
  end type
  
  interface RealSingleModeVector
    module procedure new_RealSingleModeVector
    module procedure new_RealSingleModeVector_CartesianVector
    module procedure new_RealSingleModeVector_String
  end interface
  
  interface operator(*)
    module procedure multiply_real_RealSingleModeVector
    module procedure multiply_RealSingleModeVector_real
  end interface
  
  interface operator(/)
    module procedure divide_RealSingleModeVector_real
  end interface
  
  interface operator(+)
    module procedure add_RealSingleModeVector_RealSingleModeVector
  end interface
  
  interface operator(-)
    module procedure negative_RealSingleModeVector
    module procedure subtract_RealSingleModeVector_RealSingleModeVector
  end interface
contains

! Constructor.
function new_RealSingleModeVector(id,magnitude) result(this)
  implicit none
  
  integer,  intent(in)       :: id
  real(dp), intent(in)       :: magnitude
  type(RealSingleModeVector) :: this
  
  this%id        = id
  this%magnitude = magnitude
end function

! ----------------------------------------------------------------------
! Arithmetic.
! ----------------------------------------------------------------------
impure elemental function multiply_real_RealSingleModeVector(this,that) &
   & result(output)
  implicit none
  
  real(dp),                   intent(in) :: this
  type(RealSingleModeVector), intent(in) :: that
  type(RealSingleModeVector)             :: output
  
  output = RealSingleModeVector( id        = that%id, &
                               & magnitude = this*that%magnitude)
end function

impure elemental function multiply_RealSingleModeVector_real(this,that) &
   & result(output)
  implicit none
  
  type(RealSingleModeVector), intent(in) :: this
  real(dp),                   intent(in) :: that
  type(RealSingleModeVector)             :: output
  
  output = RealSingleModeVector( id        = this%id, &
                               & magnitude = this%magnitude*that)
end function

impure elemental function divide_RealSingleModeVector_real(this,that) &
   & result(output)
  implicit none
  
  type(RealSingleModeVector), intent(in) :: this
  real(dp),                   intent(in) :: that
  type(RealSingleModeVector)             :: output
  
  output = RealSingleModeVector( id        = this%id, &
                               & magnitude = this%magnitude/that)
end function

impure elemental function add_RealSingleModeVector_RealSingleModeVector(this, &
   & that)  result(output)
  implicit none
  
  type(RealSingleModeVector), intent(in) :: this
  type(RealSingleModeVector), intent(in) :: that
  type(RealSingleModeVector)             :: output
  
  if (this%id/=that%id) then
    call print_line(ERROR//': Trying to add vectors along different modes.')
    call err()
  endif
  
  output = RealSingleModeVector( id        = this%id, &
                               & magnitude = this%magnitude+that%magnitude)
end function

impure elemental function negative_RealSingleModeVector(this) result(output)
  implicit none
  
  type(RealSingleModeVector), intent(in) :: this
  type(RealSingleModeVector)             :: output
  
  output = RealSingleModeVector(id=this%id, magnitude=-this%magnitude)
end function

impure elemental function subtract_RealSingleModeVector_RealSingleModeVector( &
   & this,that) result(output)
  implicit none
  
  type(RealSingleModeVector), intent(in) :: this
  type(RealSingleModeVector), intent(in) :: that
  type(RealSingleModeVector)             :: output
  
  if (this%id/=that%id) then
    call print_line(ERROR//': Trying to subtract vectors along different &
       &modes.')
    call err()
  endif
  
  output = RealSingleModeVector( id        = this%id, &
                               & magnitude = this%magnitude-that%magnitude)
end function

! ----------------------------------------------------------------------
! Conversions to and from cartesian co-ordinates.
! ----------------------------------------------------------------------
! Constructs the CartesianVector corresponding to this vector.
function cartesian_vector_RealSingleModeVector(this,real_mode,structure, &
   & qpoint) result(output)
  implicit none
  
  class(RealSingleModeVector), intent(in) :: this
  type(RealMode),              intent(in) :: real_mode
  type(StructureData),         intent(in) :: structure
  type(QpointData),            intent(in) :: qpoint
  type(CartesianVector)                   :: output
  
  if (real_mode%id/=this%id) then
    call print_line(CODE_ERROR//': Mode and vector incompatible.')
    call err()
  endif
  
  output = this%magnitude * real_mode%cartesian_vector(structure,qpoint)
end function

! Constructs the vector corresponding to the component of a cartesian
!    vector along this mode.
function new_RealSingleModeVector_CartesianVector(mode, &
   & vector,structure,qpoint) result(this)
  implicit none
  
  type(RealMode),        intent(in) :: mode
  type(CartesianVector), intent(in) :: vector
  type(StructureData),   intent(in) :: structure
  type(QpointData),      intent(in) :: qpoint
  type(RealSingleModeVector)        :: this
  
  type(CartesianVector) :: mode_vector
  
  real(dp) :: numerator
  real(dp) :: denominator
  
  integer :: i
  
  ! Modes are orthogonal in mass-reduced co-ordinates,
  !    but normalised in cartesian co-ordinates.
  ! If M is the mass-weighting matrix, M_ab = 1/sqrt(m_a*m_b),
  !    where m_a and m_b are the masses of atoms a and b respectively, then
  !
  ! u_i and u_j are the cartesian representations of modes i and j.
  ! r is the cartesian vector.
  !
  ! u_i.u_i=1   for all i (Modes are normalised in cartesian co-ordinates.)
  ! u_i.M.u_j=0 for i/=j  (Modes are orthogonal in mass-reduced co-ordinates.
  !
  ! r = sum_j a_j*u_j
  ! => u_i.M.r = sum_j a_j*u_i.M.u_j = a_i*u_i.M.u_i
  ! => a_j = u_i.M.r / u_i.M.u_i
  
  mode_vector = mode%cartesian_vector(structure,qpoint)
  
  numerator = 0
  denominator = 0
  do i=1,structure%no_atoms
    numerator = numerator               &
            & + vector%vectors(i) &
            & * mode_vector%vectors(i)  &
            & / structure%atoms(i)%mass()
    denominator = denominator            &
              & + mode_vector%vectors(i) &
              & * mode_vector%vectors(i) &
              & / structure%atoms(i)%mass()
  enddo
  
  this = RealSingleModeVector( id        = mode%id, &
                             & magnitude = numerator/denominator)
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_RealSingleModeVector(this,input)
  implicit none
  
  class(RealSingleModeVector), intent(out) :: this
  type(String),                intent(in)  :: input
  
  type(String), allocatable :: split_string(:)
  integer                   :: id
  real(dp)                  :: magnitude
  
  select type(this); type is(RealSingleModeVector)
    split_string = split_line(input)
    if (size(split_string)/=3) then
      call print_line(ERROR//': unable to parse real single mode vector &
         &from string: '//input)
      call err()
    endif
    
    ! If e.g. id=3 and power=2.1 then split_string = ["u3","=","2.1"]
    ! The 'u' needs stripping off the first element to give the id.
    id = int(slice(split_string(1),2,len(split_string(1))))
    magnitude = dble(split_string(3))
    
    this = RealSingleModeVector(id,magnitude)
  end select
end subroutine

function write_RealSingleModeVector(this) result(output)
  implicit none
  
  class(RealSingleModeVector), intent(in) :: this
  type(String)                            :: output
  
  select type(this); type is(RealSingleModeVector)
    output = 'u'//this%id//' = '//this%magnitude
  end select
end function

impure elemental function new_RealSingleModeVector_String(input) &
   & result(this)
  implicit none
  
  type(String), intent(in)   :: input
  type(RealSingleModeVector) :: this
  
  this = input
end function
end module
