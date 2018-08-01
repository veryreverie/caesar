! ======================================================================
! A displacement along a single complex mode.
! ======================================================================
module complex_single_mode_displacement_submodule
  use utils_module
  
  use complex_mode_submodule
  implicit none
  
  private
  
  public :: ComplexSingleDisplacement
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: operator(-)
  public :: select_mode
  public :: select_modes
  
  type, extends(Stringable) :: ComplexSingleDisplacement
    ! The id of the mode.
    integer :: id
    
    ! The magnitude of the displacement along the mode.
    complex(dp) :: magnitude
  contains
    procedure, public :: read  => read_ComplexSingleDisplacement
    procedure, public :: write => write_ComplexSingleDisplacement
  end type
  
  interface ComplexSingleDisplacement
    module procedure new_ComplexSingleDisplacement
    module procedure new_ComplexSingleDisplacement_ComplexMode
    module procedure new_ComplexSingleDisplacement_String
  end interface
  
  interface operator(*)
    module procedure multiply_real_ComplexSingleDisplacement
    module procedure multiply_ComplexSingleDisplacement_real
    module procedure multiply_complex_ComplexSingleDisplacement
    module procedure multiply_ComplexSingleDisplacement_complex
  end interface
  
  interface operator(/)
    module procedure divide_ComplexSingleDisplacement_complex
  end interface
  
  interface operator(+)
    module procedure add_ComplexSingleDisplacement_ComplexSingleDisplacement
  end interface
  
  interface operator(-)
    module procedure negative_ComplexSingleDisplacement
    module procedure &
       & subtract_ComplexSingleDisplacement_ComplexSingleDisplacement
  end interface
  
  interface select_mode
    module procedure select_mode_ComplexSingleDisplacement
  end interface
  
  interface select_modes
    module procedure select_modes_ComplexSingleDisplacements
  end interface
contains

! Constructors.
impure elemental function new_ComplexSingleDisplacement(id,magnitude) &
   & result(this)
  implicit none
  
  integer,     intent(in)         :: id
  complex(dp), intent(in)         :: magnitude
  type(ComplexSingleDisplacement) :: this
  
  this%id        = id
  this%magnitude = magnitude
end function

impure elemental function new_ComplexSingleDisplacement_ComplexMode(mode, &
   & magnitude) result(this)
  implicit none
  
  type(ComplexMode), intent(in)   :: mode
  complex(dp),       intent(in)   :: magnitude
  type(ComplexSingleDisplacement) :: this
  
  this = ComplexSingleDisplacement(id=mode%id, magnitude=magnitude)
end function

! Arithmetic.
impure elemental function multiply_real_ComplexSingleDisplacement(this,that) &
   & result(output)
  implicit none
  
  real(dp),                        intent(in) :: this
  type(ComplexSingleDisplacement), intent(in) :: that
  type(ComplexSingleDisplacement)             :: output
  
  output = ComplexSingleDisplacement( id        = that%id,            &
                                    & magnitude = this*that%magnitude )
end function

impure elemental function multiply_ComplexSingleDisplacement_real(this,that) &
   & result(output)
  implicit none
  
  type(ComplexSingleDisplacement), intent(in) :: this
  real(dp),                        intent(in) :: that
  type(ComplexSingleDisplacement)             :: output
  
  output = ComplexSingleDisplacement( id        = this%id,            &
                                    & magnitude = this%magnitude*that )
end function

impure elemental function multiply_complex_ComplexSingleDisplacement(this, &
   & that) result(output)
  implicit none
  
  complex(dp),                     intent(in) :: this
  type(ComplexSingleDisplacement), intent(in) :: that
  type(ComplexSingleDisplacement)             :: output
  
  output = ComplexSingleDisplacement( id        = that%id,            &
                                    & magnitude = this*that%magnitude )
end function

impure elemental function multiply_ComplexSingleDisplacement_complex(this, &
   & that) result(output)
  implicit none
  
  type(ComplexSingleDisplacement), intent(in) :: this
  complex(dp),                     intent(in) :: that
  type(ComplexSingleDisplacement)             :: output
  
  output = ComplexSingleDisplacement( id        = this%id,            &
                                    & magnitude = this%magnitude*that )
end function

impure elemental function divide_ComplexSingleDisplacement_complex(this,that) &
   & result(output)
  implicit none
  
  type(ComplexSingleDisplacement), intent(in) :: this
  complex(dp),                     intent(in) :: that
  type(ComplexSingleDisplacement)             :: output
  
  output = ComplexSingleDisplacement( id        = this%id,            &
                                    & magnitude = this%magnitude/that )
end function

impure elemental function                                               &
   & add_ComplexSingleDisplacement_ComplexSingleDisplacement(this,that) &
   & result(output)
  implicit none
  
  type(ComplexSingleDisplacement), intent(in) :: this
  type(ComplexSingleDisplacement), intent(in) :: that
  type(ComplexSingleDisplacement)             :: output
  
  if (this%id/=that%id) then
    call print_line(ERROR//': Trying to add vectors along different modes.')
    call err()
  endif
  
  output = ComplexSingleDisplacement(            &
     & id        = this%id,                      &
     & magnitude = this%magnitude+that%magnitude )
end function

impure elemental function negative_ComplexSingleDisplacement(this) &
   & result(output)
  implicit none
  
  type(ComplexSingleDisplacement), intent(in) :: this
  type(ComplexSingleDisplacement)             :: output
  
  output = ComplexSingleDisplacement(id=this%id, magnitude=-this%magnitude)
end function

impure elemental function                                                    &
   & subtract_ComplexSingleDisplacement_ComplexSingleDisplacement(this,that) &
   & result(output)
  implicit none
  
  type(ComplexSingleDisplacement), intent(in) :: this
  type(ComplexSingleDisplacement), intent(in) :: that
  type(ComplexSingleDisplacement)             :: output
  
  if (this%id/=that%id) then
    call print_line(ERROR//': Trying to subtract vectors along different &
       &modes.')
    call err()
  endif
  
  output = ComplexSingleDisplacement(            &
     & id        = this%id,                      &
     & magnitude = this%magnitude-that%magnitude )
end function

! Select modes corresponding to a given force or forces.
function select_mode_ComplexSingleDisplacement(displacement,modes) &
   & result(output)
  implicit none
  
  type(ComplexSingleDisplacement), intent(in) :: displacement
  type(ComplexMode),               intent(in) :: modes(:)
  type(ComplexMode)                           :: output
  
  output = modes(first(modes%id==displacement%id))
end function

function select_modes_ComplexSingleDisplacements(displacements,modes) &
   & result(output)
  implicit none
  
  type(ComplexSingleDisplacement), intent(in) :: displacements(:)
  type(ComplexMode),               intent(in) :: modes(:)
  type(ComplexMode), allocatable              :: output(:)
  
  integer :: i,ialloc
  
  allocate(output(size(displacements)), stat=ialloc); call err(ialloc)
  do i=1,size(displacements)
    output(i) = select_mode(displacements(i), modes)
  enddo
end function

! I/O.
subroutine read_ComplexSingleDisplacement(this,input)
  implicit none
  
  class(ComplexSingleDisplacement), intent(out) :: this
  type(String),                     intent(in)  :: input
  
  type(String), allocatable :: split_string(:)
  integer                   :: id
  complex(dp)               :: magnitude
  
  select type(this); type is(ComplexSingleDisplacement)
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
    
    this = ComplexSingleDisplacement(id,magnitude)
  end select
end subroutine

function write_ComplexSingleDisplacement(this) result(output)
  implicit none
  
  class(ComplexSingleDisplacement), intent(in) :: this
  type(String)                                 :: output
  
  select type(this); type is(ComplexSingleDisplacement)
    output = 'u'//this%id//' = '//this%magnitude
  end select
end function

impure elemental function new_ComplexSingleDisplacement_String(input) &
   & result(this)
  implicit none
  
  type(String), intent(in)        :: input
  type(ComplexSingleDisplacement) :: this
  
  call this%read(input)
end function
end module
