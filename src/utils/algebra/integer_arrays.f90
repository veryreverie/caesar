! ======================================================================
! A heterogeneous multi-dimensional integer array.
! ======================================================================
module caesar_integer_arrays_module
  use caesar_foundations_module
  use caesar_io_module
  implicit none
  
  private
  
  public :: IntArray1D
  public :: IntArray2D
  public :: array
  public :: size
  public :: operator(==)
  public :: operator(/=)
  public :: operator(//)
  
  ! An array of integers.
  type, extends(Stringable) :: IntArray1D
    integer, public, allocatable :: i(:)
  contains
    procedure, public :: read  => read_IntArray1D
    procedure, public :: write => write_IntArray1D
  end type
  
  ! A array of IntArrays.
  type, extends(Stringsable) :: IntArray2D
    type(IntArray1D), allocatable :: i(:)
  contains
    procedure, public :: read  => read_IntArray2D
    procedure, public :: write => write_IntArray2D
  end type
  
  interface IntArray1D
    ! ----------------------------------------------------------------------
    ! Constructors.
    ! ----------------------------------------------------------------------
    module function new_IntArray1D(i) result(this) 
      integer, intent(in) :: i(:)
      type(IntArray1D)    :: this
    end function
  end interface
  
  interface IntArray2D
    module function new_IntArray2D(i) result(this) 
      type(IntArray1D), intent(in) :: i(:)
      type(IntArray2D)             :: this
    end function
  end interface
  
  interface array
    ! ----------------------------------------------------------------------
    ! Conversions from integer(:) to IntArray1D and IntArray1D(:) to IntArray2D
    ! ----------------------------------------------------------------------
    module function array_IntArray1D_integers(input) result(output) 
      integer, intent(in) :: input(:)
      type(IntArray1D)    :: output
    end function
  
    module function array_IntArray2D_IntArray1Ds(input) result(output) 
      type(IntArray1D), intent(in) :: input(:)
      type(IntArray2D)             :: output
    end function
  end interface
  
  interface size
    ! ----------------------------------------------------------------------
    ! Size module function.
    ! ----------------------------------------------------------------------
    module function size_IntArray1D(input) result(output) 
      type(IntArray1D), intent(in) :: input
      integer                      :: output
    end function
  
    module function size_IntArray2D(input) result(output) 
      type(IntArray2D), intent(in) :: input
      integer                      :: output
    end function
  end interface
  
  interface operator(//)
    ! ----------------------------------------------------------------------
    ! Concatenation operations
    ! ----------------------------------------------------------------------
    
    ! IntArray1D = IntArray1D // integer
    module function concatenate_IntArray1D_integer(this,that) result(output) 
      type(IntArray1D), intent(in) :: this
      integer,          intent(in) :: that
      type(IntArray1D)             :: output
    end function
  
    ! IntArray1D = integer // IntArray1D
    module function concatenate_integer_IntArray1D(this,that) result(output) 
      integer,          intent(in) :: this
      type(IntArray1D), intent(in) :: that
      type(IntArray1D)             :: output
    end function
  
    ! IntArray1D = IntArray1D // integer(:)
    module function concatenate_IntArray1D_integers(this,that) result(output) 
      type(IntArray1D), intent(in) :: this
      integer,          intent(in) :: that(:)
      type(IntArray1D)             :: output
    end function
  
    ! IntArray1D = integer(:) // IntArray1D
    module function concatenate_integers_IntArray1D(this,that) result(output) 
      integer,          intent(in) :: this(:)
      type(IntArray1D), intent(in) :: that
      type(IntArray1D)             :: output
    end function
  
    ! IntArray1D = IntArray1D // IntArray1D
    module function concatenate_IntArray1D_IntArray1D(this,that) &
       & result(output) 
      type(IntArray1D), intent(in) :: this
      type(IntArray1D), intent(in) :: that
      type(IntArray1D)             :: output
    end function
  
    ! IntArray2D = IntArray2D // IntArray1D
    module function concatenate_IntArray2D_IntArray1D(this,that) &
       & result(output) 
      type(IntArray2D), intent(in) :: this
      type(IntArray1D), intent(in) :: that
      type(IntArray2D)             :: output
    end function
  
    ! IntArray2D = IntArray1D // IntArray2D
    module function concatenate_IntArray1D_IntArray2D(this,that) &
       & result(output) 
      type(IntArray1D), intent(in) :: this
      type(IntArray2D), intent(in) :: that
      type(IntArray2D)             :: output
    end function
  
    ! IntArray2D = IntArray2D // IntArray1D(:)
    module function concatenate_IntArray2D_IntArray1Ds(this,that) &
       & result(output) 
      type(IntArray2D), intent(in) :: this
      type(IntArray1D), intent(in) :: that(:)
      type(IntArray2D)             :: output
    end function
  
    ! IntArray2D = IntArray1D(:) // IntArray2D
    module function concatenate_IntArray1Ds_IntArray2D(this,that) &
       & result(output) 
      type(IntArray1D), intent(in) :: this(:)
      type(IntArray2D), intent(in) :: that
      type(IntArray2D)             :: output
    end function
  
    ! IntArray2D = IntArray2D // IntArray2D
    module function concatenate_IntArray2D_IntArray2D(this,that) &
       & result(output) 
      type(IntArray2D), intent(in) :: this
      type(IntArray2D), intent(in) :: that
      type(IntArray2D)             :: output
    end function
  end interface
  
  interface operator(==)
    ! ----------------------------------------------------------------------
    ! Equality and non-equality.
    ! ----------------------------------------------------------------------
    impure elemental module function equality_IntArray1D_IntArray1D(this, &
       & that) result(output) 
      type(IntArray1D), intent(in) :: this
      type(IntArray1D), intent(in) :: that
      logical                      :: output
    end function
  
    impure elemental module function equality_IntArray2D_IntArray2D(this, &
       & that) result(output) 
      type(IntArray2D), intent(in) :: this
      type(IntArray2D), intent(in) :: that
      logical                      :: output
    end function
  end interface
  
  interface operator(/=)
    impure elemental module function non_equality_IntArray1D_IntArray1D(this,that) result(output) 
      type(IntArray1D), intent(in) :: this
      type(IntArray1D), intent(in) :: that
      logical                      :: output
    end function
  
    impure elemental module function non_equality_IntArray2D_IntArray2D(this,that) result(output) 
      type(IntArray2D), intent(in) :: this
      type(IntArray2D), intent(in) :: that
      logical                      :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_IntArray1D(this,input) 
      class(IntArray1D), intent(out) :: this
      type(String),      intent(in)  :: input
    end subroutine
  end interface
  
  interface
    module function write_IntArray1D(this) result(output) 
      class(IntArray1D), intent(in) :: this
      type(String)                  :: output
    end function
  end interface
  
  interface IntArray1D
    impure elemental module function new_IntArray1D_String(input) result(this) 
      type(String), intent(in) :: input
      type(IntArray1D)         :: this
    end function
  end interface
  
  interface
    module subroutine read_IntArray2D(this,input) 
      class(IntArray2D), intent(out) :: this
      type(String),      intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_IntArray2D(this) result(output) 
      class(IntArray2D), intent(in) :: this
      type(String), allocatable     :: output(:)
    end function
  end interface
  
  interface IntArray2D
    module function new_IntArray2D_Strings(input) result(this) 
      type(String), intent(in) :: input(:)
      type(IntArray2D)         :: this
    end function
  
    impure elemental module function new_IntArray2D_StringArray(input) &
       & result(this) 
      type(StringArray), intent(in) :: input
      type(IntArray2D)              :: this
    end function
  end interface
end module
