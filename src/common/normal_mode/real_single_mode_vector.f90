! ======================================================================
! A vector along a single real mode, in mass-weighted co-ordinates.
! ======================================================================
module real_single_mode_vector_submodule
  use utils_module
  
  use structure_module
  
  use cartesian_vector_submodule
  use mass_weighted_vector_submodule
  use real_mode_submodule
  implicit none
  
  private
  
  public :: RealSingleModeVector
  public :: MassWeightedVector
  public :: CartesianVector
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: operator(-)
  
  type, extends(Stringable) :: RealSingleModeVector
    ! The id of the mode.
    integer :: id
    
    ! The magnitude of vector along the mode, in mass-weighted co-ordinates.
    real(dp) :: magnitude
  contains
    ! I/O.
    procedure, public :: read  => read_RealSingleModeVector
    procedure, public :: write => write_RealSingleModeVector
  end type
  
  interface RealSingleModeVector
    module procedure new_RealSingleModeVector
    module procedure new_RealSingleModeVector_MassWeightedVector
    module procedure new_RealSingleModeVector_CartesianVector
    module procedure new_RealSingleModeVector_String
  end interface
  
  interface MassWeightedVector
    module procedure new_MassWeightedVector_RealSingleModeVector
  end interface
  
  interface CartesianVector
    module procedure new_CartesianVector_RealSingleModeVector
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
! Conversions to and from mass-weighted cartesian co-ordinates.
! ----------------------------------------------------------------------
! Constructs the MassWeightedVector corresponding to this vector.
function new_MassWeightedVector_RealSingleModeVector(this,real_mode, &
   & structure,qpoint) result(output)
  implicit none
  
  class(RealSingleModeVector), intent(in) :: this
  type(RealMode),              intent(in) :: real_mode
  type(StructureData),         intent(in) :: structure
  type(QpointData),            intent(in) :: qpoint
  type(MassWeightedVector)                :: output
  
  if (real_mode%id/=this%id) then
    call print_line(CODE_ERROR//': Mode and vector incompatible.')
    call err()
  endif
  
  output = this%magnitude * MassWeightedVector(real_mode,structure,qpoint)
end function

! Constructs the CartesianVector corresponding to this vector.
function new_CartesianVector_RealSingleModeVector(this,real_mode,structure, &
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
  
  output = CartesianVector(                                 &
     & MassWeightedVector(this,real_mode,structure,qpoint), &
     & structure)
end function

! Constructs the vector corresponding to the component of a mass-weighted
!    vector along this mode.
function new_RealSingleModeVector_MassWeightedVector(mode,vector,structure, &
   & qpoint) result(this)
  implicit none
  
  type(RealMode),           intent(in) :: mode
  type(MassWeightedVector), intent(in) :: vector
  type(StructureData),      intent(in) :: structure
  type(QpointData),         intent(in) :: qpoint
  type(RealSingleModeVector)           :: this
  
  type(MassWeightedVector) :: mode_vector
  
  real(dp) :: magnitude
  
  ! Modes are orthonormal in mass-reduced co-ordinates.
  ! If M is the mass-weighting matrix, M_ab = 1/sqrt(m_a*m_b),
  !    where m_a and m_b are the masses of atoms a and b respectively, then
  !
  ! u_i and u_j are the cartesian representations of modes i and j.
  ! r is the cartesian vector.
  !
  ! u_i.M.u_j = 0 if i/=j
  !           = n if i=j, where n is the number of primitive cells.
  !
  ! r = sum_j[ a_j*u_j ]
  ! => u_i.M.r = sum_j[ a_j*u_i.M.u_j ] = a_i*u_i.M.u_i = a_i*n
  ! => a_j = u_i.M.r / n
  
  mode_vector = MassWeightedVector(mode,structure,qpoint)
  
  magnitude = sum( vector%vectors       &
          &      * mode_vector%vectors) &
          & / structure%sc_size
  
  this = RealSingleModeVector(id=mode%id, magnitude=magnitude)
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
  
  this = RealSingleModeVector( mode,                                 &
                             & MassWeightedVector(vector,structure), &
                             & structure,                            &
                             & qpoint)
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
