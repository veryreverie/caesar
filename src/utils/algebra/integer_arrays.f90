! ======================================================================
! A heterogeneous multi-dimensional integer array.
! ======================================================================
module integer_arrays_submodule
  use precision_module
  use io_module
  implicit none
  
  private
  
  public :: IntArray1D
  public :: IntArray2D
  public :: array
  public :: size
  public :: operator(==)
  public :: operator(/=)
  
  ! An array of integers.
  type, extends(Stringable) :: IntArray1D
    integer, public, allocatable :: i(:)
  contains
    generic,   public  :: assignment (= ) => assign_IntArray1D_integers
    procedure, private ::                    assign_IntArray1D_integers
    
    generic, public :: operator (//) => concatenate_IntArray1D_integer,  &
                                      & concatenate_integer_IntArray1D,  &
                                      & concatenate_IntArray1D_integers, &
                                      & concatenate_integers_IntArray1D, &
                                      & concatenate_IntArray1D_IntArray1D
    procedure, private             ::   concatenate_IntArray1D_integer
    procedure, private, pass(that) ::   concatenate_integer_IntArray1D
    procedure, private             ::   concatenate_IntArray1D_integers
    procedure, private, pass(that) ::   concatenate_integers_IntArray1D
    procedure, private             ::   concatenate_IntArray1D_IntArray1D
    
    procedure, public :: read  => read_IntArray1D
    procedure, public :: write => write_IntArray1D
  end type
  
  ! A array of IntArrays.
  type, extends(Stringsable) :: IntArray2D
    type(IntArray1D), allocatable :: i(:)
  contains
    generic,   public  :: assignment (= ) => assign_IntArray2D_IntArray1Ds
    procedure, private ::                    assign_IntArray2D_IntArray1Ds
    
    generic, public :: operator (//) => concatenate_IntArray2D_IntArray1D,  &
                                      & concatenate_IntArray1D_IntArray2D,  &
                                      & concatenate_IntArray2D_IntArray1Ds, &
                                      & concatenate_IntArray1Ds_IntArray2D, &
                                      & concatenate_IntArray2D_IntArray2D
    procedure, private             ::   concatenate_IntArray2D_IntArray1D
    procedure, private, pass(that) ::   concatenate_IntArray1D_IntArray2D
    procedure, private             ::   concatenate_IntArray2D_IntArray1Ds
    procedure, private, pass(that) ::   concatenate_IntArray1Ds_IntArray2D
    procedure, private             ::   concatenate_IntArray2D_IntArray2D
    
    procedure, public :: read  => read_IntArray2D
    procedure, public :: write => write_IntArray2D
  end type
  
  ! Type conversion to array.
  interface array
    module procedure array_IntArray1D_integers
    module procedure array_IntArray2D_IntArray1Ds
  end interface
  
  ! Size.
  interface size
    module procedure size_IntArray1D
    module procedure size_IntArray2D
  end interface
  
  ! Equality and non-equality.
  interface operator(==)
    module procedure equality_IntArray1D_IntArray1D
    module procedure equality_IntArray2D_IntArray2D
  end interface
  
  interface operator(/=)
    module procedure non_equality_IntArray1D_IntArray1D
    module procedure non_equality_IntArray2D_IntArray2D
  end interface
  
contains

! ----------------------------------------------------------------------
! Conversions from integer(:) to IntArray1D and IntArray1D(:) to IntArray2D
! ----------------------------------------------------------------------
subroutine assign_IntArray1D_integers(this,that)
  implicit none
  
  class(IntArray1D), intent(out) :: this
  integer,           intent(in)  :: that(:)
  
  this%i = that
end subroutine

subroutine assign_IntArray2D_IntArray1Ds(this,that)
  implicit none
  
  class(IntArray2D), intent(out) :: this
  type(IntArray1D),  intent(in)  :: that(:)
  
  this%i = that
end subroutine

function array_IntArray1D_integers(input) result(output)
  implicit none
  
  integer, intent(in) :: input(:)
  type(IntArray1D)    :: output
  
  output = input
end function

function array_IntArray2D_IntArray1Ds(input) result(output)
  implicit none
  
  type(IntArray1D), intent(in) :: input(:)
  type(IntArray2D)             :: output
  
  output = input
end function

! ----------------------------------------------------------------------
! Concatenation operations
! ----------------------------------------------------------------------

! IntArray1D = IntArray1D // integer
function concatenate_IntArray1D_integer(this,that) result(output)
  implicit none
  
  class(IntArray1D), intent(in) :: this
  integer,           intent(in) :: that
  type(IntArray1D)              :: output
  
  output = [this%i, that]
end function

! IntArray1D = integer // IntArray1D
function concatenate_integer_IntArray1D(this,that) result(output)
  implicit none
  
  integer,           intent(in) :: this
  class(IntArray1D), intent(in) :: that
  type(IntArray1D)              :: output
  
  output = [this, that%i]
end function

! IntArray1D = IntArray1D // integer(:)
function concatenate_IntArray1D_integers(this,that) result(output)
  implicit none
  
  class(IntArray1D), intent(in) :: this
  integer,           intent(in) :: that(:)
  type(IntArray1D)              :: output
  
  output = [this%i, that]
end function

! IntArray1D = integer(:) // IntArray1D
function concatenate_integers_IntArray1D(this,that) result(output)
  implicit none
  
  integer,           intent(in) :: this(:)
  class(IntArray1D), intent(in) :: that
  type(IntArray1D)              :: output
  
  output = [this, that%i]
end function

! IntArray1D = IntArray1D // IntArray1D
function concatenate_IntArray1D_IntArray1D(this,that) result(output)
  implicit none
  
  class(IntArray1D), intent(in) :: this
  class(IntArray1D), intent(in) :: that
  type(IntArray1D)              :: output
  
  output = [this%i, that%i]
end function

! IntArray2D = IntArray2D // IntArray1D
function concatenate_IntArray2D_IntArray1D(this,that) result(output)
  implicit none
  
  class(IntArray2D), intent(in) :: this
  type(IntArray1D),  intent(in) :: that
  type(IntArray2D)              :: output
  
  output = [this%i, that]
end function

! IntArray2D = IntArray1D // IntArray2D
function concatenate_IntArray1D_IntArray2D(this,that) result(output)
  implicit none
  
  type(IntArray1D),  intent(in) :: this
  class(IntArray2D), intent(in) :: that
  type(IntArray2D)              :: output
  
  output = [this, that%i]
end function

! IntArray2D = IntArray2D // IntArray1D(:)
function concatenate_IntArray2D_IntArray1Ds(this,that) result(output)
  implicit none
  
  class(IntArray2D), intent(in) :: this
  type(IntArray1D),  intent(in) :: that(:)
  type(IntArray2D)              :: output
  
  output = [this%i, that]
end function

! IntArray2D = IntArray1D(:) // IntArray2D
function concatenate_IntArray1Ds_IntArray2D(this,that) result(output)
  implicit none
  
  type(IntArray1D),  intent(in) :: this(:)
  class(IntArray2D), intent(in) :: that
  type(IntArray2D)              :: output
  
  output = [this, that%i]
end function

! IntArray2D = IntArray2D // IntArray2D
function concatenate_IntArray2D_IntArray2D(this,that) result(output)
  implicit none
  
  class(IntArray2D), intent(in) :: this
  class(IntArray2D), intent(in) :: that
  type(IntArray2D)              :: output
  
  output = [this%i, that%i]
end function

! ----------------------------------------------------------------------
! Size function.
! ----------------------------------------------------------------------
function size_IntArray1D(input) result(output)
  implicit none
  
  type(IntArray1D), intent(in) :: input
  integer                      :: output
  
  output = size(input%i)
end function

function size_IntArray2D(input) result(output)
  implicit none
  
  type(IntArray2D), intent(in) :: input
  integer                      :: output
  
  output = size(input%i)
end function

! ----------------------------------------------------------------------
! Equality and non-equality.
! ----------------------------------------------------------------------
impure elemental function equality_IntArray1D_IntArray1D(this,that) &
   & result(output)
  implicit none
  
  type(IntArray1D), intent(in) :: this
  type(IntArray1D), intent(in) :: that
  logical                      :: output
  
  output = all(this%i==that%i)
end function

impure elemental function equality_IntArray2D_IntArray2D(this,that) &
   & result(output)
  implicit none
  
  type(IntArray2D), intent(in) :: this
  type(IntArray2D), intent(in) :: that
  logical                      :: output
  
  output = all(this%i==that%i)
end function

impure elemental function non_equality_IntArray1D_IntArray1D(this,that) &
   & result(output)
  implicit none
  
  type(IntArray1D), intent(in) :: this
  type(IntArray1D), intent(in) :: that
  logical                      :: output
  
  output = .not. this==that
end function

impure elemental function non_equality_IntArray2D_IntArray2D(this,that) &
   & result(output)
  implicit none
  
  type(IntArray2D), intent(in) :: this
  type(IntArray2D), intent(in) :: that
  logical                      :: output
  
  output = .not. this==that
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_IntArray1D(this,input)
  implicit none
  
  class(IntArray1D), intent(out) :: this
  type(String),      intent(in)  :: input
  
  select type(this); type is(IntArray1D)
    this = IntArray1D(int(split(input)))
  end select
end subroutine

function write_IntArray1D(this) result(output)
  implicit none
  
  class(IntArray1D), intent(in) :: this
  type(String)                  :: output
  
  select type(this); type is(IntArray1D)
    output = join(this%i)
  end select
end function

subroutine read_IntArray2D(this,input)
  implicit none
  
  class(IntArray2D), intent(out) :: this
  type(String),      intent(in)  :: input(:)
  
  type(IntArray1D), allocatable :: contents(:)
  
  integer :: i,ialloc
  
  select type(this); type is(IntArray2D)
    allocate(contents(size(input)), stat=ialloc); call err(ialloc)
    do i=1,size(contents)
      contents(i) = input(i)
    enddo
    this = IntArray2D(contents)
  end select
end subroutine

function write_IntArray2D(this) result(output)
  implicit none
  
  class(IntArray2D), intent(in) :: this
  type(String), allocatable     :: output(:)
  
  integer :: i,ialloc
  
  select type(this); type is(IntArray2D)
    allocate(output(size(this)), stat=ialloc); call err(ialloc)
    do i=1,size(this)
      output(i) = str(this%i(i))
    enddo
  end select
end function
end module
