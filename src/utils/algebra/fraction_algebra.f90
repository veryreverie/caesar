!> Provides the [[FractionVector(type)]] and [[FractionMatrix(type)]] classes
!>    and related methods.
module caesar_fraction_algebra_module
  use caesar_foundations_module
  use caesar_io_module
  
  use caesar_fraction_module
  use caesar_linear_algebra_module
  use caesar_algebra_utils_module
  implicit none
  
  private
  
  public :: FractionVector
  public :: FractionMatrix
  public :: vec
  public :: mat
  public :: frac
  public :: fracvec
  public :: fracmat
  public :: intvec
  public :: intmat
  public :: dblevec
  public :: dblemat
  public :: size
  public :: is_int
  public :: operator(==)
  public :: operator(/=)
  public :: transpose
  public :: invert
  public :: operator(+)
  public :: operator(-)
  public :: operator(*)
  public :: operator(/)
  public :: sum
  
  !> A vector whose elements are of type [[IntFraction(type)]].
  type, extends(Stringable) :: FractionVector
    type(IntFraction), allocatable, private :: contents_(:)
  contains
    procedure, public :: read  => read_FractionVector
    procedure, public :: write => write_FractionVector
  end type
  
  !> A matrix whose elements are of type [[IntFraction(type)]].
  type, extends(Stringsable) :: FractionMatrix
    type(IntFraction), allocatable, private :: contents_(:,:)
  contains
    procedure, public :: read  => read_FractionMatrix
    procedure, public :: write => write_FractionMatrix
  end type
  
  interface vec
    !> Conversion from [[IntFraction(type)]] array to [[FractionVector(type)]].
    module function vec_IntFractions(input) result(output) 
      type(IntFraction), intent(in) :: input(:)
      type(FractionVector)          :: output
    end function
  end interface
  
  interface mat
    !> Conversion from a 2D [[IntFraction(type)]] array to
    !>    [[FractionMatrix(type)]].
    module function mat_IntFractions(input) result(output) 
      type(IntFraction), intent(in) :: input(:,:)
      type(FractionMatrix)          :: output
    end function
  
    !> Conversion from an array of [[IntFraction(type)]] matrix elements,
    !>    and the matrix `shape` to a [[FractionMatrix(type)]].
    !> N.B. `input` should be in row-major order.
    module function mat_IntFractions_shape(input,shape) result(output) 
      type(IntFraction), intent(in) :: input(:)
      integer,           intent(in) :: shape(2)
      type(FractionMatrix)          :: output
    end function
  end interface
  
  interface frac
    !> Conversion from [[FractionVector(type)]] to an array of
    !>    [[IntFraction(type)]].
    module function frac_FractionVector(input) result(output) 
      type(FractionVector), intent(in) :: input
      type(IntFraction), allocatable   :: output(:)
    end function
  
    !> Conversion from [[FractionMatrix(type)]] to an array of
    !>    [[IntFraction(type)]].
    module function frac_FractionMatrix(input) result(output) 
      type(FractionMatrix), intent(in) :: input
      type(IntFraction), allocatable   :: output(:,:)
    end function
  end interface
  
  interface fracvec
    !> Conversion from [[IntVector(type)]] to [[FractionVector(type)]].
    impure elemental module function fracvec_IntVector(input) result(output) 
      type(IntVector), intent(in) :: input
      type(FractionVector)        :: output
    end function
  end interface
  
  interface fracmat
    !> Conversion from [[IntMatrix(type)]] to [[FractionMatrix(type)]].
    impure elemental module function fracmat_IntMatrix(input) result(output) 
      type(IntMatrix), intent(in)    :: input
      type(FractionMatrix)           :: output
    end function
  end interface
  
  interface intvec
    !> Conversion from [[FractionVector(type)]] to [[IntVector(type)]].
    impure elemental module function intvec_FractionVector(input) &
       & result(output) 
      type(FractionVector), intent(in) :: input
      type(IntVector)                  :: output
    end function
  end interface
  
  interface intmat
    !> Conversion from [[FractionMatrix(type)]] to [[IntMatrix(type)]].
    impure elemental module function intmat_FractionMatrix(input) &
       & result(output) 
      type(FractionMatrix), intent(in) :: input
      type(IntMatrix)                  :: output
    end function
  end interface
  
  interface dblevec
    !> Conversion from [[FractionVector(type)]] to [[RealVector(type)]].
    impure elemental module function dblevec_FractionVector(input) &
       & result(output) 
      type(FractionVector), intent(in) :: input
      type(RealVector)                 :: output
    end function
  end interface
  
  interface dblemat
    !> Conversion from [[FractionMatrix(type)]] to [[RealMatrix(type)]].
    impure elemental module function dblemat_FractionMatrix(input) &
       & result(output) 
      type(FractionMatrix), intent(in) :: input
      type(RealMatrix)                 :: output
    end function
  end interface
  
  interface size
    !> Returns the number of elements of a [[FractionVector(type)]].
    module function size_FractionVector(this) result(output) 
      type(FractionVector), intent(in) :: this
      integer                          :: output
    end function
  
    !> Returns the number of elements of a [[FractionMatrix(type)]]
    !>    along the given `dim`.
    module function size_FractionMatrix(this,dim) result(output) 
      type(FractionMatrix), intent(in) :: this
      integer,              intent(in) :: dim
      integer                          :: output
    end function
  end interface
  
  interface is_int
    !> Returns whether or not every element of `this` is an integer.
    impure elemental module function is_int_FractionVector(this) &
       & result(output) 
      type(FractionVector), intent(in) :: this
      logical                          :: output
    end function
  
    !> Returns whether or not every element of `this` is an integer.
    impure elemental module function is_int_FractionMatrix(this) &
       & result(output) 
      type(FractionMatrix), intent(in) :: this
      logical                          :: output
    end function
  end interface
  
  interface operator(==)
    !> Equality comparison between [[FractionVector(type)]] and
    !>    [[FractionVector(type)]].
    impure elemental module function equality_FractionVector_FractionVector( &
       & this,that) result(output) 
      type(FractionVector), intent(in) :: this
      type(FractionVector), intent(in) :: that
      logical                          :: output
    end function
  
    !> Equality comparison between [[FractionVector(type)]] and
    !>    [[IntVector(type)]].
    impure elemental module function equality_FractionVector_IntVector(this, &
       & that) result(output) 
      type(FractionVector), intent(in) :: this
      type(IntVector),      intent(in) :: that
      logical                          :: output
    end function
  
    !> Equality comparison between [[IntVector(type)]] and
    !>    [[FractionVector(type)]].
    impure elemental module function equality_IntVector_FractionVector(this, &
       & that) result(output) 
      type(IntVector),      intent(in) :: this
      type(FractionVector), intent(in) :: that
      logical                          :: output
    end function
  
    !> Equality comparison between [[FractionMatrix(type)]] and
    !>    [[FractionMatrix(type)]].
    impure elemental module function equality_FractionMatrix_FractionMatrix( &
       & this,that) result(output) 
      type(FractionMatrix), intent(in) :: this
      type(FractionMatrix), intent(in) :: that
      logical                          :: output
    end function
  
    !> Equality comparison between [[FractionMatrix(type)]] and
    !>    [[IntMatrix(type)]].
    impure elemental module function equality_FractionMatrix_IntMatrix(this, &
       & that) result(output) 
      type(FractionMatrix), intent(in) :: this
      type(IntMatrix),      intent(in) :: that
      logical                          :: output
    end function
  
    !> Equality comparison between [[IntMatrix(type)]] and
    !>    [[FractionMatrix(type)]].
    impure elemental module function equality_IntMatrix_FractionMatrix(this, &
       & that) result(output) 
      type(IntMatrix),      intent(in) :: this
      type(FractionMatrix), intent(in) :: that
      logical                          :: output
    end function
  end interface
  
  interface operator(/=)
    !> Non-equality comparison between [[FractionVector(type)]] and
    !>    [[FractionVector(type)]].
    impure elemental module function &
       & non_equality_FractionVector_FractionVector(this,that) result(output) 
      type(FractionVector), intent(in) :: this
      type(FractionVector), intent(in) :: that
      logical                          :: output
    end function
  
    !> Non-equality comparison between [[FractionVector(type)]] and
    !>    [[IntVector(type)]].
    impure elemental module function non_equality_FractionVector_IntVector( &
       & this,that) result(output) 
      type(FractionVector), intent(in) :: this
      type(IntVector),      intent(in) :: that
      logical                          :: output
    end function
  
    !> Non-equality comparison between [[IntVector(type)]] and
    !>    [[FractionVector(type)]].
    impure elemental module function non_equality_IntVector_FractionVector( &
       & this,that) result(output) 
      type(IntVector),      intent(in) :: this
      type(FractionVector), intent(in) :: that
      logical                          :: output
    end function
  
    !> Non-equality comparison between [[FractionMatrix(type)]] and
    !>    [[FractionMatrix(type)]].
    impure elemental module function &
       & non_equality_FractionMatrix_FractionMatrix(this,that) result(output) 
      type(FractionMatrix), intent(in) :: this
      type(FractionMatrix), intent(in) :: that
      logical                          :: output
    end function
  
    !> Non-equality comparison between [[FractionMatrix(type)]] and
    !>    [[IntMatrix(type)]].
    impure elemental module function &
       & non_equality_FractionMatrix_IntMatrix(this,that) result(output) 
      type(FractionMatrix), intent(in) :: this
      type(IntMatrix),      intent(in) :: that
      logical                          :: output
    end function
  
    !> Non-equality comparison between [[IntMatrix(type)]] and
    !>    [[FractionMatrix(type)]].
    impure elemental module function &
       & non_equality_IntMatrix_FractionMatrix(this,that) result(output) 
      type(IntMatrix),      intent(in) :: this
      type(FractionMatrix), intent(in) :: that
      logical                          :: output
    end function
  end interface
  
  interface transpose
    !> Returns the matrix transpose of a [[FractionMatrix(type)]].
    module function transpose_FractionMatrix(this) result(output) 
      type(FractionMatrix), intent(in) :: this
      type(FractionMatrix)             :: output
    end function
  end interface
  
  interface invert
    !> Returns the inverse of a 3x3 `integer` matrix,
    !>    as a 3x3 [[IntFraction(type)]] matrix.
    module function invert_integers(input) result(output) 
      integer, intent(in) :: input(:,:)
      type(IntFraction)   :: output(3,3)
    end function
  
    !> Returns the inverse of a 3x3 [[IntMatrix(type)]],
    !>    as a 3x3 [[FractionMatrix(type)]].
    module function invert_IntMatrix(input) result(output) 
      type(IntMatrix), intent(in) :: input
      type(FractionMatrix)        :: output
    end function
  end interface
  
  interface operator(+)
    !> Vector addition.
    impure elemental module function add_FractionVector_FractionVector(this, &
       & that) result(output) 
      type(FractionVector), intent(in) :: this
      type(FractionVector), intent(in) :: that
      type(FractionVector)             :: output
    end function
  
    !> Vector addition.
    impure elemental module function add_FractionVector_IntVector(this,that) &
       & result(output) 
      type(IntVector),      intent(in) :: this
      type(FractionVector), intent(in) :: that
      type(FractionVector)             :: output
    end function
  
    !> Vector addition.
    impure elemental module function add_IntVector_FractionVector(this,that) &
       & result(output) 
      type(FractionVector), intent(in) :: this
      type(IntVector),      intent(in) :: that
      type(FractionVector)             :: output
    end function
  
    !> Matrix addition.
    impure elemental module function add_FractionMatrix_FractionMatrix(this, &
       & that) result(output) 
      type(FractionMatrix), intent(in) :: this
      type(FractionMatrix), intent(in) :: that
      type(FractionMatrix)             :: output
    end function
  
    !> Matrix addition.
    impure elemental module function add_FractionMatrix_IntMatrix(this,that) &
       & result(output) 
      type(IntMatrix),      intent(in) :: this
      type(FractionMatrix), intent(in) :: that
      type(FractionMatrix)             :: output
    end function
  
    !> Matrix addition.
    impure elemental module function add_IntMatrix_FractionMatrix(this,that) &
       & result(output) 
      type(FractionMatrix), intent(in) :: this
      type(IntMatrix),      intent(in) :: that
      type(FractionMatrix)             :: output
    end function
  end interface
  
  interface operator(-)
    !> Negative of a vector.
    impure elemental module function negative_FractionVector(this) &
       & result(output) 
      type(FractionVector), intent(in) :: this
      type(FractionVector)             :: output
    end function
  
    !> Negative of a matrix.
    impure elemental module function negative_FractionMatrix(this) &
       & result(output) 
      type(FractionMatrix), intent(in) :: this
      type(FractionMatrix)             :: output
    end function
  
    !> Vector subtraction.
    impure elemental module function subtract_FractionVector_FractionVector( &
       & this,that) result(output) 
      type(FractionVector), intent(in) :: this
      type(FractionVector), intent(in) :: that
      type(FractionVector)             :: output
    end function
  
    !> Vector subtraction.
    impure elemental module function subtract_FractionVector_IntVector(this, &
       & that) result(output) 
      type(FractionVector), intent(in) :: this
      type(IntVector),      intent(in) :: that
      type(FractionVector)             :: output
    end function
  
    !> Vector subtraction.
    impure elemental module function subtract_IntVector_FractionVector(this, &
       & that) result(output) 
      type(IntVector),      intent(in) :: this
      type(FractionVector), intent(in) :: that
      type(FractionVector)             :: output
    end function
  
    !> Matrix subtraction.
    impure elemental module function subtract_FractionMatrix_FractionMatrix( &
       & this,that) result(output) 
      type(FractionMatrix), intent(in) :: this
      type(FractionMatrix), intent(in) :: that
      type(FractionMatrix)             :: output
    end function
  
    !> Matrix subtraction.
    impure elemental module function subtract_FractionMatrix_IntMatrix(this, &
       & that) result(output) 
      type(FractionMatrix), intent(in) :: this
      type(IntMatrix),      intent(in) :: that
      type(FractionMatrix)             :: output
    end function
  
    !> Matrix subtraction.
    impure elemental module function subtract_IntMatrix_FractionMatrix(this, &
       & that) result(output) 
      type(IntMatrix),      intent(in) :: this
      type(FractionMatrix), intent(in) :: that
      type(FractionMatrix)             :: output
    end function
  end interface
  
  interface operator(*)
    !> Multiplication of a vector by a scalar.
    impure elemental module function multiply_FractionVector_integer(this, &
       & that) result(output) 
      type(FractionVector), intent(in) :: this
      integer,              intent(in) :: that
      type(FractionVector)             :: output
    end function
  
    !> Multiplication of a vector by a scalar.
    impure elemental module function multiply_integer_FractionVector(this, &
       & that) result(output) 
      integer,              intent(in) :: this
      type(FractionVector), intent(in) :: that
      type(FractionVector)             :: output
    end function
  
    !> Multiplication of a vector by a scalar.
    impure elemental module function multiply_FractionVector_IntFraction( &
       & this,that) result(output) 
      type(FractionVector), intent(in) :: this
      type(IntFraction),    intent(in) :: that
      type(FractionVector)             :: output
    end function
  
    !> Multiplication of a vector by a scalar.
    impure elemental module function multiply_IntFraction_FractionVector( &
       & this,that) result(output) 
      type(IntFraction),    intent(in) :: this
      type(FractionVector), intent(in) :: that
      type(FractionVector)             :: output
    end function
  
    !> Multiplication of a matrix by a scalar.
    impure elemental module function multiply_FractionMatrix_integer(this, &
       & that) result(output) 
      type(FractionMatrix), intent(in) :: this
      integer,              intent(in) :: that
      type(FractionMatrix)             :: output
    end function
  
    !> Multiplication of a matrix by a scalar.
    impure elemental module function multiply_integer_FractionMatrix(this, &
       & that) result(output) 
      integer,              intent(in) :: this
      type(FractionMatrix), intent(in) :: that
      type(FractionMatrix)             :: output
    end function
  
    !> Multiplication of a matrix by a scalar.
    impure elemental module function multiply_FractionMatrix_IntFraction( &
       & this,that) result(output) 
      type(FractionMatrix), intent(in) :: this
      type(IntFraction),    intent(in) :: that
      type(FractionMatrix)             :: output
    end function
  
    !> Multiplication of a matrix by a scalar.
    impure elemental module function multiply_IntFraction_FractionMatrix( &
       & this,that) result(output) 
      type(IntFraction),    intent(in) :: this
      type(FractionMatrix), intent(in) :: that
      type(FractionMatrix)             :: output
    end function
  
    !> Vector-vector dot product.
    impure elemental module function dot_FractionVector_FractionVector(this, &
       & that) result(output) 
      type(FractionVector), intent(in) :: this
      type(FractionVector), intent(in) :: that
      type(IntFraction)                :: output
    end function
  
    !> Vector-vector dot product.
    impure elemental module function dot_FractionVector_IntVector(this,that) &
       & result(output) 
      type(FractionVector), intent(in) :: this
      type(IntVector),      intent(in) :: that
      type(IntFraction)                :: output
    end function
  
    !> Vector-vector dot product.
    impure elemental module function dot_IntVector_FractionVector(this,that) &
       & result(output) 
      type(IntVector),      intent(in) :: this
      type(FractionVector), intent(in) :: that
      type(IntFraction)                :: output
    end function
  
    !> Matrix-vector dot product.
    impure elemental module function dot_FractionMatrix_FractionVector(this, &
       & that) result(output) 
      type(FractionMatrix), intent(in) :: this
      type(FractionVector), intent(in) :: that
      type(FractionVector)             :: output
    end function
  
    !> Matrix-vector dot product.
    impure elemental module function dot_FractionMatrix_IntVector(this,that) &
       & result(output) 
      type(FractionMatrix), intent(in) :: this
      type(IntVector),      intent(in) :: that
      type(FractionVector)             :: output
    end function
  
    !> Matrix-vector dot product.
    impure elemental module function dot_IntMatrix_FractionVector(this,that) &
       & result(output) 
      type(IntMatrix),      intent(in) :: this
      type(FractionVector), intent(in) :: that
      type(FractionVector)             :: output
    end function
  
    !> Vector-matrix dot product.
    impure elemental module function dot_FractionVector_FractionMatrix(this, &
       & that) result(output) 
      type(FractionVector), intent(in) :: this
      type(FractionMatrix), intent(in) :: that
      type(FractionVector)             :: output
    end function
  
    !> Vector-matrix dot product.
    impure elemental module function dot_FractionVector_IntMatrix(this,that) &
       & result(output) 
      type(FractionVector), intent(in) :: this
      type(IntMatrix),      intent(in) :: that
      type(FractionVector)             :: output
    end function
  
    !> Vector-matrix dot product.
    impure elemental module function dot_IntVector_FractionMatrix(this,that) &
       & result(output) 
      type(IntVector),      intent(in) :: this
      type(FractionMatrix), intent(in) :: that
      type(FractionVector)             :: output
    end function
  
    !> Matrix-matrix dot product.
    impure elemental module function dot_FractionMatrix_FractionMatrix(this, &
       & that) result(output) 
      type(FractionMatrix), intent(in) :: this
      type(FractionMatrix), intent(in) :: that
      type(FractionMatrix)             :: output
    end function
  
    !> Matrix-matrix dot product.
    impure elemental module function dot_FractionMatrix_IntMatrix(this,that) &
       & result(output) 
      type(FractionMatrix), intent(in) :: this
      type(IntMatrix),      intent(in) :: that
      type(FractionMatrix)             :: output
    end function
  
    !> Matrix-matrix dot product.
    impure elemental module function dot_IntMatrix_FractionMatrix(this,that) &
       & result(output) 
      type(IntMatrix),      intent(in) :: this
      type(FractionMatrix), intent(in) :: that
      type(FractionMatrix)             :: output
    end function
  end interface
  
  interface operator(/)
    !> Division of a vector by a scalar.
    impure elemental module function divide_FractionVector_integer(this,that) &
       & result(output) 
      type(FractionVector), intent(in) :: this
      integer,              intent(in) :: that
      type(FractionVector)             :: output
    end function
  
    !> Division of a vector by a scalar.
    impure elemental module function divide_FractionVector_IntFraction(this, &
       & that) result(output) 
      type(FractionVector), intent(in) :: this
      type(IntFraction),    intent(in) :: that
      type(FractionVector)             :: output
    end function
  
    !> Division of a matrix by a scalar.
    impure elemental module function divide_FractionMatrix_integer(this,that) &
       & result(output) 
      type(FractionMatrix), intent(in) :: this
      integer,              intent(in) :: that
      type(FractionMatrix)             :: output
    end function
  
    !> Division of a matrix by a scalar.
    impure elemental module function divide_FractionMatrix_IntFraction(this, &
       & that) result(output) 
      type(FractionMatrix), intent(in) :: this
      type(IntFraction),    intent(in) :: that
      type(FractionMatrix)             :: output
    end function
  end interface
  
  interface sum
    !> Returns the sum of a [[FractionVector(type)]] array.
    module function sum_FractionVector(input) result(output) 
      type(FractionVector), intent(in) :: input(:)
      type(FractionVector)             :: output
    end function
  
    !> Returns the sum of a [[FractionMatrix(type)]] array.
    module function sum_FractionMatrix(input) result(output) 
      type(FractionMatrix), intent(in) :: input(:)
      type(FractionMatrix)             :: output
    end function
  end interface
  
  interface
    !> Conversion from [[String(type)]] to [[FractionVector(type)]].
    module subroutine read_FractionVector(this,input) 
      class(FractionVector), intent(out) :: this
      type(String),          intent(in)  :: input
    end subroutine
  end interface
  
  interface
    !> Conversion from [[FractionVector(type)]] to [[String(type)]].
    module function write_FractionVector(this) result(output) 
      class(FractionVector), intent(in) :: this
      type(String)                      :: output
    end function
  end interface
  
  interface FractionVector
    !> Conversion from [[String(type)]] to [[FractionVector(type)]].
    impure elemental module function new_FractionVector_String(input) &
       & result(this) 
      type(String), intent(in) :: input
      type(FractionVector)     :: this
    end function
  end interface
  
  interface
    !> Conversion from [[String(type)]] array to [[FractionVector(type)]].
    module subroutine read_FractionMatrix(this,input) 
      class(FractionMatrix), intent(out) :: this
      type(String),          intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    !> Conversion from [[FractionVector(type)]] to [[String(type)]] array.
    module function write_FractionMatrix(this) result(output) 
      class(FractionMatrix), intent(in) :: this
      type(String), allocatable         :: output(:)
    end function
  end interface
  
  interface FractionMatrix
    !> Conversion from [[String(type)]] array to [[FractionVector(type)]].
    module function new_FractionMatrix_Strings(input) result(this) 
      type(String), intent(in) :: input(:)
      type(FractionMatrix)     :: this
    end function
  
    !> Conversion from [[StringArray(type)]] to [[FractionVector(type)]].
    impure elemental module function new_FractionMatrix_StringArray(input) &
       & result(this) 
      type(StringArray), intent(in) :: input
      type(FractionMatrix)          :: this
    end function
  end interface
end module
