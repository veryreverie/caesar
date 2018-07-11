! ======================================================================
! A displacement in complex mode co-ordinates.
! ======================================================================
module complex_mode_displacement_submodule
  use utils_module
  
  use structure_module
  
  use complex_mode_submodule
  use complex_single_mode_displacement_submodule
  implicit none
  
  private
  
  public :: ComplexModeDisplacement
  public :: size
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: sum
  public :: operator(-)
  
  type, extends(Stringsable) :: ComplexModeDisplacement
    type(ComplexSingleDisplacement), allocatable :: vectors(:)
  contains
    ! The component of the displacement along a given mode.
    generic,   public  :: displacement =>  &
                        & displacement_id, &
                        & displacement_mode
    procedure, private :: displacement_id
    procedure, private :: displacement_mode
    
    ! I/O.
    procedure, public :: read  => read_ComplexModeDisplacement
    procedure, public :: write => write_ComplexModeDisplacement
  end type
  
  interface ComplexModeDisplacement
    module procedure new_ComplexModeDisplacement
    module procedure new_ComplexModeDisplacement_ComplexModes
    module procedure new_ComplexModeDisplacement_Strings
    module procedure new_ComplexModeDisplacement_StringArray
  end interface
  
  interface size
    module procedure size_ComplexModeDisplacement
  end interface
  
  interface operator(*)
    module procedure multiply_real_ComplexModeDisplacement
    module procedure multiply_ComplexModeDisplacement_real
    module procedure multiply_complex_ComplexModeDisplacement
    module procedure multiply_ComplexModeDisplacement_complex
  end interface
  
  interface operator(/)
    module procedure divide_ComplexModeDisplacement_complex
  end interface
  
  interface operator(+)
    module procedure add_ComplexModeDisplacement_ComplexModeDisplacement
  end interface
  
  interface sum
    module procedure sum_ComplexModeDisplacements
  end interface
  
  interface operator(-)
    module procedure negative_ComplexModeDisplacement
    module procedure subtract_ComplexModeDisplacement_ComplexModeDisplacement
  end interface
contains

! Constructors and size() function.
function new_ComplexModeDisplacement(displacements) result(this)
  implicit none
  
  type(ComplexSingleDisplacement), intent(in) :: displacements(:)
  type(ComplexModeDisplacement)               :: this
  
  this%vectors = displacements
end function

function new_ComplexModeDisplacement_ComplexModes(modes,displacements) &
   & result(this)
  implicit none
  
  type(ComplexMode), intent(in) :: modes(:)
  complex(dp),       intent(in) :: displacements(:)
  type(ComplexModeDisplacement) :: this
  
  this = ComplexModeDisplacement(ComplexSingleDisplacement( modes,        &
                                                          & displacements ))
end function

function size_ComplexModeDisplacement(this) result(output)
  implicit none
  
  type(ComplexModeDisplacement), intent(in) :: this
  integer                                   :: output
  
  output = size(this%vectors)
end function

! ----------------------------------------------------------------------
! Arithmetic.
! ----------------------------------------------------------------------
impure elemental function multiply_real_ComplexModeDisplacement(this,that) &
   & result(output)
  implicit none
  
  real(dp),                      intent(in) :: this
  type(ComplexModeDisplacement), intent(in) :: that
  type(ComplexModeDisplacement)             :: output
  
  output = ComplexModeDisplacement(this*that%vectors)
end function

impure elemental function multiply_ComplexModeDisplacement_real(this,that) &
   & result(output)
  implicit none
  
  type(ComplexModeDisplacement), intent(in) :: this
  real(dp),                      intent(in) :: that
  type(ComplexModeDisplacement)             :: output
  
  output = ComplexModeDisplacement(this%vectors*that)
end function

impure elemental function multiply_complex_ComplexModeDisplacement(this,that) &
   & result(output)
  implicit none
  
  complex(dp),                   intent(in) :: this
  type(ComplexModeDisplacement), intent(in) :: that
  type(ComplexModeDisplacement)             :: output
  
  output = ComplexModeDisplacement(this*that%vectors)
end function

impure elemental function multiply_ComplexModeDisplacement_complex(this,that) &
   & result(output)
  implicit none
  
  type(ComplexModeDisplacement), intent(in) :: this
  complex(dp),                   intent(in) :: that
  type(ComplexModeDisplacement)             :: output
  
  output = ComplexModeDisplacement(this%vectors*that)
end function

impure elemental function divide_ComplexModeDisplacement_complex(this,that) &
   & result(output)
  implicit none
  
  type(ComplexModeDisplacement), intent(in) :: this
  complex(dp),                   intent(in) :: that
  type(ComplexModeDisplacement)             :: output
  
  output = ComplexModeDisplacement(this%vectors/that)
end function

impure elemental function add_ComplexModeDisplacement_ComplexModeDisplacement(&
   & this,that) result(output)
  implicit none
  
  type(ComplexModeDisplacement), intent(in) :: this
  type(ComplexModeDisplacement), intent(in) :: that
  type(ComplexModeDisplacement)             :: output
  
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

function sum_ComplexModeDisplacements(this) result(output)
  implicit none
  
  type(ComplexModeDisplacement), intent(in) :: this(:)
  type(ComplexModeDisplacement)             :: output
  
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

impure elemental function negative_ComplexModeDisplacement(this) result(output)
  implicit none
  
  type(ComplexModeDisplacement), intent(in) :: this
  type(ComplexModeDisplacement)             :: output
  
  output = ComplexModeDisplacement(-this%vectors)
end function

impure elemental function                                                &
   & subtract_ComplexModeDisplacement_ComplexModeDisplacement(this,that) &
   & result(output)
  implicit none
  
  type(ComplexModeDisplacement), intent(in) :: this
  type(ComplexModeDisplacement), intent(in) :: that
  type(ComplexModeDisplacement)             :: output
  
  output = this + (-that)
end function

! ----------------------------------------------------------------------
! The displacement along a given mode.
! ----------------------------------------------------------------------
impure elemental function displacement_id(this,id) result(output)
  implicit none
  
  class(ComplexModeDisplacement), intent(in) :: this
  integer,                        intent(in) :: id
  complex(dp)                                :: output
  
  type(ComplexSingleDisplacement) :: displacement
  
  integer :: i
  
  i = first(this%vectors%id==id, default=0)
  if (i==0) then
    output = 0.0_dp
  else
    output = this%vectors(i)%magnitude
  endif
end function

impure elemental function displacement_mode(this,mode) result(output)
  implicit none
  
  class(ComplexModeDisplacement), intent(in) :: this
  type(ComplexMode),              intent(in) :: mode
  complex(dp)                                :: output
  
  output = this%displacement(mode%id)
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_ComplexModeDisplacement(this,input)
  implicit none
  
  class(ComplexModeDisplacement), intent(out) :: this
  type(String),                   intent(in)  :: input(:)
  
  select type(this); type is(ComplexModeDisplacement)
    this = ComplexModeDisplacement(ComplexSingleDisplacement(input))
  class default
    call err()
  end select
end subroutine

function write_ComplexModeDisplacement(this) result(output)
  implicit none
  
  class(ComplexModeDisplacement), intent(in) :: this
  type(String), allocatable                  :: output(:)
  
  select type(this); type is(ComplexModeDisplacement)
    output = str(this%vectors)
  class default
    call err()
  end select
end function

function new_ComplexModeDisplacement_Strings(input) result(this)
  implicit none
  
  type(String), intent(in)      :: input(:)
  type(ComplexModeDisplacement) :: this
  
  call this%read(input)
end function

impure elemental function new_ComplexModeDisplacement_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(ComplexModeDisplacement) :: this
  
  this = ComplexModeDisplacement(str(input))
end function
end module
