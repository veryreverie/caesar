! ======================================================================
! A force in complex mode co-ordinates.
! ======================================================================
module complex_mode_force_submodule
  use utils_module
  
  use structure_module
  
  use complex_mode_submodule
  use complex_single_mode_force_submodule
  implicit none
  
  private
  
  public :: ComplexModeForce
  public :: size
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: sum
  public :: operator(-)
  
  type, extends(Stringsable) :: ComplexModeForce
    type(ComplexSingleForce), allocatable :: vectors(:)
  contains
    ! The component of the force along a given mode.
    generic,   public  :: force =>  &
                        & force_id, &
                        & force_mode
    procedure, private :: force_id
    procedure, private :: force_mode
    
    ! I/O.
    procedure, public :: read  => read_ComplexModeForce
    procedure, public :: write => write_ComplexModeForce
  end type
  
  interface ComplexModeForce
    module procedure new_ComplexModeForce
    module procedure new_ComplexModeForce_ComplexModes
    module procedure new_ComplexModeForce_Strings
    module procedure new_ComplexModeForce_StringArray
  end interface
  
  interface size
    module procedure size_ComplexModeForce
  end interface
  
  interface operator(*)
    module procedure multiply_real_ComplexModeForce
    module procedure multiply_ComplexModeForce_real
    module procedure multiply_complex_ComplexModeForce
    module procedure multiply_ComplexModeForce_complex
  end interface
  
  interface operator(/)
    module procedure divide_ComplexModeForce_complex
  end interface
  
  interface operator(+)
    module procedure add_ComplexModeForce_ComplexModeForce
  end interface
  
  interface sum
    module procedure sum_ComplexModeForces
  end interface
  
  interface operator(-)
    module procedure negative_ComplexModeForce
    module procedure subtract_ComplexModeForce_ComplexModeForce
  end interface
contains

! Constructors and size() function.
function new_ComplexModeForce(forces) result(this)
  implicit none
  
  type(ComplexSingleForce), intent(in) :: forces(:)
  type(ComplexModeForce)               :: this
  
  this%vectors = forces
end function

function new_ComplexModeForce_ComplexModes(modes,forces) result(this)
  implicit none
  
  type(ComplexMode), intent(in) :: modes(:)
  complex(dp),       intent(in) :: forces(:)
  type(ComplexModeForce)        :: this
  
  this = ComplexModeForce(ComplexSingleForce(modes,forces))
end function

function size_ComplexModeForce(this) result(output)
  implicit none
  
  type(ComplexModeForce), intent(in) :: this
  integer                            :: output
  
  output = size(this%vectors)
end function

! ----------------------------------------------------------------------
! Arithmetic.
! ----------------------------------------------------------------------
impure elemental function multiply_real_ComplexModeForce(this,that) &
   & result(output)
  implicit none
  
  real(dp),               intent(in) :: this
  type(ComplexModeForce), intent(in) :: that
  type(ComplexModeForce)             :: output
  
  output = ComplexModeForce(this*that%vectors)
end function

impure elemental function multiply_ComplexModeForce_real(this,that) &
   & result(output)
  implicit none
  
  type(ComplexModeForce), intent(in) :: this
  real(dp),               intent(in) :: that
  type(ComplexModeForce)             :: output
  
  output = ComplexModeForce(this%vectors*that)
end function

impure elemental function multiply_complex_ComplexModeForce(this,that) &
   & result(output)
  implicit none
  
  complex(dp),            intent(in) :: this
  type(ComplexModeForce), intent(in) :: that
  type(ComplexModeForce)             :: output
  
  output = ComplexModeForce(this*that%vectors)
end function

impure elemental function multiply_ComplexModeForce_complex(this,that) &
   & result(output)
  implicit none
  
  type(ComplexModeForce), intent(in) :: this
  complex(dp),            intent(in) :: that
  type(ComplexModeForce)             :: output
  
  output = ComplexModeForce(this%vectors*that)
end function

impure elemental function divide_ComplexModeForce_complex(this,that) &
   & result(output)
  implicit none
  
  type(ComplexModeForce), intent(in) :: this
  complex(dp),            intent(in) :: that
  type(ComplexModeForce)             :: output
  
  output = ComplexModeForce(this%vectors/that)
end function

impure elemental function add_ComplexModeForce_ComplexModeForce(this, &
   & that) result(output)
  implicit none
  
  type(ComplexModeForce), intent(in) :: this
  type(ComplexModeForce), intent(in) :: that
  type(ComplexModeForce)             :: output
  
  integer :: i,j
  
  output = this
  do i=1,size(that)
    j = first(this%vectors%id==that%vectors(i)%id, default=0)
    if (j==0) then
      output%vectors = [output%vectors, that%vectors(i)]
    else
      output%vectors(j) = output%vectors(j) + that%vectors(i)
    endif
  enddo
end function

function sum_ComplexModeForces(this) result(output)
  implicit none
  
  type(ComplexModeForce), intent(in) :: this(:)
  type(ComplexModeForce)             :: output
  
  integer :: i
  
  if (size(this)==0) then
    call print_line(ERROR//': Trying to sum an empty list.')
    call err()
  endif
  
  output = this(1)
  do i=2,size(this)
    output = output + this(i)
  enddo
end function

impure elemental function negative_ComplexModeForce(this) result(output)
  implicit none
  
  type(ComplexModeForce), intent(in) :: this
  type(ComplexModeForce)             :: output
  
  output = ComplexModeForce(-this%vectors)
end function

impure elemental function subtract_ComplexModeForce_ComplexModeForce( &
   & this,that) result(output)
  implicit none
  
  type(ComplexModeForce), intent(in) :: this
  type(ComplexModeForce), intent(in) :: that
  type(ComplexModeForce)             :: output
  
  output = this + (-that)
end function

! ----------------------------------------------------------------------
! The displacement along a given mode.
! ----------------------------------------------------------------------
impure elemental function force_id(this,id) result(output)
  implicit none
  
  class(ComplexModeForce), intent(in) :: this
  integer,                 intent(in) :: id
  complex(dp)                         :: output
  
  type(ComplexSingleForce) :: force
  
  integer :: i
  
  i = first(this%vectors%id==id, default=0)
  if (i==0) then
    output = 0.0_dp
  else
    output = this%vectors(i)%magnitude
  endif
end function

impure elemental function force_mode(this,mode) result(output)
  implicit none
  
  class(ComplexModeForce), intent(in) :: this
  type(ComplexMode),       intent(in) :: mode
  complex(dp)                         :: output
  
  output = this%force(mode%id)
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_ComplexModeForce(this,input)
  implicit none
  
  class(ComplexModeForce), intent(out) :: this
  type(String),            intent(in)  :: input(:)
  
  select type(this); type is(ComplexModeForce)
    this = ComplexModeForce(ComplexSingleForce(input))
  class default
    call err()
  end select
end subroutine

function write_ComplexModeForce(this) result(output)
  implicit none
  
  class(ComplexModeForce), intent(in) :: this
  type(String), allocatable           :: output(:)
  
  select type(this); type is(ComplexModeForce)
    output = str(this%vectors)
  class default
    call err()
  end select
end function

function new_ComplexModeForce_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(ComplexModeForce)   :: this
  
  call this%read(input)
end function

impure elemental function new_ComplexModeForce_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(ComplexModeForce)        :: this
  
  this = ComplexModeForce(str(input))
end function
end module
