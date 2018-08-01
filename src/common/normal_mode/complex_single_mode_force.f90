! ======================================================================
! A force along a single complex mode.
! ======================================================================
module complex_single_mode_force_submodule
  use utils_module
  
  use complex_mode_submodule
  implicit none
  
  private
  
  public :: ComplexSingleForce
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: operator(-)
  public :: select_mode
  public :: select_modes
  
  type, extends(Stringable) :: ComplexSingleForce
    ! The id of the mode.
    integer :: id
    
    ! The magnitude of the force along the mode.
    complex(dp) :: magnitude
  contains
    procedure, public :: read  => read_ComplexSingleForce
    procedure, public :: write => write_ComplexSingleForce
  end type
  
  interface ComplexSingleForce
    module procedure new_ComplexSingleForce
    module procedure new_ComplexSingleForce_ComplexMode
    module procedure new_ComplexSingleForce_String
  end interface
  
  interface operator(*)
    module procedure multiply_real_ComplexSingleForce
    module procedure multiply_ComplexSingleForce_real
    module procedure multiply_complex_ComplexSingleForce
    module procedure multiply_ComplexSingleForce_complex
  end interface
  
  interface operator(/)
    module procedure divide_ComplexSingleForce_complex
  end interface
  
  interface operator(+)
    module procedure add_ComplexSingleForce_ComplexSingleForce
  end interface
  
  interface operator(-)
    module procedure negative_ComplexSingleForce
    module procedure subtract_ComplexSingleForce_ComplexSingleForce
  end interface
  
  interface select_mode
    module procedure select_mode_ComplexSingleForce
  end interface
  
  interface select_modes
    module procedure select_modes_ComplexSingleForces
  end interface
contains

! Constructors.
impure elemental function new_ComplexSingleForce(id,magnitude) result(this)
  implicit none
  
  integer,     intent(in)  :: id
  complex(dp), intent(in)  :: magnitude
  type(ComplexSingleForce) :: this
  
  this%id        = id
  this%magnitude = magnitude
end function

impure elemental function new_ComplexSingleForce_ComplexMode(mode,magnitude) &
   & result(this)
  implicit none
  
  type(ComplexMode), intent(in) :: mode
  complex(dp),       intent(in) :: magnitude
  type(ComplexSingleForce)      :: this
  
  this = ComplexSingleForce(id=mode%id, magnitude=magnitude)
end function

! Arithmetic.
impure elemental function multiply_real_ComplexSingleForce(this,that) &
   & result(output)
  implicit none
  
  real(dp),                 intent(in) :: this
  type(ComplexSingleForce), intent(in) :: that
  type(ComplexSingleForce)             :: output
  
  output = ComplexSingleForce( id        = that%id,            &
                             & magnitude = this*that%magnitude )
end function

impure elemental function multiply_ComplexSingleForce_real(this,that) &
   & result(output)
  implicit none
  
  type(ComplexSingleForce), intent(in) :: this
  real(dp),                 intent(in) :: that
  type(ComplexSingleForce)             :: output
  
  output = ComplexSingleForce( id        = this%id,            &
                             & magnitude = this%magnitude*that )
end function

impure elemental function multiply_complex_ComplexSingleForce(this,that) &
   & result(output)
  implicit none
  
  complex(dp),              intent(in) :: this
  type(ComplexSingleForce), intent(in) :: that
  type(ComplexSingleForce)             :: output
  
  output = ComplexSingleForce( id        = that%id,            &
                             & magnitude = this*that%magnitude )
end function

impure elemental function multiply_ComplexSingleForce_complex(this,that) &
   & result(output)
  implicit none
  
  type(ComplexSingleForce), intent(in) :: this
  complex(dp),              intent(in) :: that
  type(ComplexSingleForce)             :: output
  
  output = ComplexSingleForce( id        = this%id,            &
                             & magnitude = this%magnitude*that )
end function

impure elemental function divide_ComplexSingleForce_complex(this,that) &
   & result(output)
  implicit none
  
  type(ComplexSingleForce), intent(in) :: this
  complex(dp),              intent(in) :: that
  type(ComplexSingleForce)             :: output
  
  output = ComplexSingleForce( id        = this%id,            &
                             & magnitude = this%magnitude/that )
end function

impure elemental function add_ComplexSingleForce_ComplexSingleForce( &
   & this,that) result(output)
  implicit none
  
  type(ComplexSingleForce), intent(in) :: this
  type(ComplexSingleForce), intent(in) :: that
  type(ComplexSingleForce)             :: output
  
  if (this%id/=that%id) then
    call print_line(ERROR//': Trying to add vectors along different modes.')
    call err()
  endif
  
  output = ComplexSingleForce( id        = this%id,                      &
                             & magnitude = this%magnitude+that%magnitude )
end function

impure elemental function negative_ComplexSingleForce(this) result(output)
  implicit none
  
  type(ComplexSingleForce), intent(in) :: this
  type(ComplexSingleForce)             :: output
  
  output = ComplexSingleForce(id=this%id, magnitude=-this%magnitude)
end function

impure elemental function subtract_ComplexSingleForce_ComplexSingleForce( &
   & this,that) result(output)
  implicit none
  
  type(ComplexSingleForce), intent(in) :: this
  type(ComplexSingleForce), intent(in) :: that
  type(ComplexSingleForce)             :: output
  
  if (this%id/=that%id) then
    call print_line(ERROR//': Trying to subtract vectors along different &
       &modes.')
    call err()
  endif
  
  output = ComplexSingleForce( id        = this%id,                      &
                             & magnitude = this%magnitude-that%magnitude )
end function

! Select modes corresponding to a given force or forces.
function select_mode_ComplexSingleForce(force,modes) result(output)
  implicit none
  
  type(ComplexSingleForce), intent(in) :: force
  type(ComplexMode),        intent(in) :: modes(:)
  type(ComplexMode)                    :: output
  
  output = modes(first(modes%id==force%id))
end function

function select_modes_ComplexSingleForces(forces,modes) result(output)
  implicit none
  
  type(ComplexSingleForce), intent(in) :: forces(:)
  type(ComplexMode),        intent(in) :: modes(:)
  type(ComplexMode), allocatable       :: output(:)
  
  integer :: i,ialloc
  
  allocate(output(size(forces)), stat=ialloc); call err(ialloc)
  do i=1,size(forces)
    output(i) = select_mode(forces(i), modes)
  enddo
end function

! I/O.
subroutine read_ComplexSingleForce(this,input)
  implicit none
  
  class(ComplexSingleForce), intent(out) :: this
  type(String),              intent(in)  :: input
  
  type(String), allocatable :: split_string(:)
  integer                   :: id
  complex(dp)               :: magnitude
  
  select type(this); type is(ComplexSingleForce)
    split_string = split_line(input)
    if (size(split_string)/=3) then
      call print_line(ERROR//': unable to parse complex single mode &
         &vector from string: '//input)
      call err()
    endif
    
    ! If e.g. id=3 and power=2.1+1.2i then split_string = ["u3","=","2.1+1.2i"]
    ! The 'u' needs stripping off the first element to give the id.
    id = int(slice(split_string(1),2,len(split_string(1))))
    magnitude = cmplx(split_string(3))
    
    this = ComplexSingleForce(id,magnitude)
  end select
end subroutine

function write_ComplexSingleForce(this) result(output)
  implicit none
  
  class(ComplexSingleForce), intent(in) :: this
  type(String)                          :: output
  
  select type(this); type is(ComplexSingleForce)
    output = 'u'//this%id//' = '//this%magnitude
  end select
end function

impure elemental function new_ComplexSingleForce_String(input) &
   & result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(ComplexSingleForce) :: this
  
  call this%read(input)
end function
end module
