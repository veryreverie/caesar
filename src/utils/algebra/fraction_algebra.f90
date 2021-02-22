! ======================================================================
! Vectors and Matrices of type IntFraction.
! ======================================================================
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
  public :: exp_2pii
  public :: cos_2pi
  public :: sin_2pi
  
  type, extends(Stringable) :: FractionVector
    type(IntFraction), allocatable, private :: contents_(:)
  contains
    procedure, public :: read  => read_FractionVector
    procedure, public :: write => write_FractionVector
  end type
  
  type, extends(Stringsable) :: FractionMatrix
    type(IntFraction), allocatable, private :: contents_(:,:)
  contains
    procedure, public :: read  => read_FractionMatrix
    procedure, public :: write => write_FractionMatrix
  end type
  
  interface vec
    ! ----------------------------------------------------------------------
    ! Procedures involving contents_
    ! ----------------------------------------------------------------------
    
    ! Conversion to and from vector and matrix types.
    module function vec_IntFractions(input) result(output) 
      type(IntFraction), intent(in) :: input(:)
      type(FractionVector)          :: output
    end function
  end interface
  
  interface mat
    module function mat_IntFractions(input) result(output) 
      type(IntFraction), intent(in) :: input(:,:)
      type(FractionMatrix)          :: output
    end function
  
    module function mat_IntFractions_shape(input,m,n) result(output) 
      type(IntFraction), intent(in) :: input(:)
      integer,           intent(in) :: m
      integer,           intent(in) :: n
      type(FractionMatrix)          :: output
    end function
  end interface
  
  interface frac
    ! Conversion to fraction(:). Effectively getters for contents_.
    module function frac_FractionVector(input) result(output) 
      type(FractionVector), intent(in) :: input
      type(IntFraction), allocatable   :: output(:)
    end function
  
    module function frac_FractionMatrix(input) result(output) 
      type(FractionMatrix), intent(in) :: input
      type(IntFraction), allocatable   :: output(:,:)
    end function
  end interface
  
  interface fracvec
    ! ----------------------------------------------------------------------
    ! Procedures not involving contents_
    ! ----------------------------------------------------------------------
    ! N.B. the number of procedures accessing contents_ directly is intentionally
    !    limited for stability reasons.
    ! The above procedures behave well if contents_ has not been allocated,
    !    and this good behaviour is automatically passed to the procedures below.
    
    ! Conversion to and from vector and matrix types.
    
    impure elemental module function fracvec_IntVector(input) result(output) 
      type(IntVector), intent(in) :: input
      type(FractionVector)        :: output
    end function
  end interface
  
  interface fracmat
    impure elemental module function fracmat_IntMatrix(input) result(output) 
      type(IntMatrix), intent(in)    :: input
      type(FractionMatrix)           :: output
    end function
  end interface
  
  interface intvec
    impure elemental module function intvec_FractionVector(input) &
       & result(output) 
      type(FractionVector), intent(in) :: input
      type(IntVector)                  :: output
    end function
  end interface
  
  interface intmat
    impure elemental module function intmat_FractionMatrix(input) &
       & result(output) 
      type(FractionMatrix), intent(in) :: input
      type(IntMatrix)                  :: output
    end function
  end interface
  
  interface dblevec
    impure elemental module function dblevec_FractionVector(input) &
       & result(output) 
      type(FractionVector), intent(in) :: input
      type(RealVector)                 :: output
    end function
  end interface
  
  interface dblemat
    impure elemental module function dblemat_FractionMatrix(input) &
       & result(output) 
      type(FractionMatrix), intent(in) :: input
      type(RealMatrix)                 :: output
    end function
  end interface
  
  interface size
    ! Properties of the vectors and matrices.
    module function size_FractionVector(this) result(output) 
      type(FractionVector), intent(in) :: this
      integer                          :: output
    end function
  
    module function size_FractionMatrix(this,dim) result(output) 
      type(FractionMatrix), intent(in) :: this
      integer,              intent(in) :: dim
      integer                          :: output
    end function
  end interface
  
  interface is_int
    impure elemental module function is_int_FractionVector(this) &
       & result(output) 
      type(FractionVector), intent(in) :: this
      logical                          :: output
    end function
  
    impure elemental module function is_int_FractionMatrix(this) &
       & result(output) 
      type(FractionMatrix), intent(in) :: this
      logical                          :: output
    end function
  end interface
  
  interface operator(==)
    ! Comparisons.
    impure elemental module function equality_FractionVector_FractionVector(this,that) result(output) 
      type(FractionVector), intent(in) :: this
      type(FractionVector), intent(in) :: that
      logical                          :: output
    end function
  
    impure elemental module function equality_FractionVector_IntVector(this, &
       & that) result(output) 
      type(FractionVector), intent(in) :: this
      type(IntVector),      intent(in) :: that
      logical                          :: output
    end function
  
    impure elemental module function equality_IntVector_FractionVector(this, &
       & that) result(output) 
      type(IntVector),      intent(in) :: this
      type(FractionVector), intent(in) :: that
      logical                          :: output
    end function
  
    impure elemental module function equality_FractionMatrix_FractionMatrix(this,that) result(output) 
      type(FractionMatrix), intent(in) :: this
      type(FractionMatrix), intent(in) :: that
      logical                          :: output
    end function
  
    impure elemental module function equality_FractionMatrix_IntMatrix(this, &
       & that) result(output) 
      type(FractionMatrix), intent(in) :: this
      type(IntMatrix),      intent(in) :: that
      logical                          :: output
    end function
  
    impure elemental module function equality_IntMatrix_FractionMatrix(this, &
       & that) result(output) 
      type(IntMatrix),      intent(in) :: this
      type(FractionMatrix), intent(in) :: that
      logical                          :: output
    end function
  end interface
  
  interface operator(/=)
    impure elemental module function non_equality_FractionVector_FractionVector(this,that) result(output) 
      type(FractionVector), intent(in) :: this
      type(FractionVector), intent(in) :: that
      logical                          :: output
    end function
  
    impure elemental module function non_equality_FractionVector_IntVector(this,that) result(output) 
      type(FractionVector), intent(in) :: this
      type(IntVector),      intent(in) :: that
      logical                          :: output
    end function
  
    impure elemental module function non_equality_IntVector_FractionVector(this,that) result(output) 
      type(IntVector),      intent(in) :: this
      type(FractionVector), intent(in) :: that
      logical                          :: output
    end function
  
    impure elemental module function non_equality_FractionMatrix_FractionMatrix(this,that) result(output) 
      type(FractionMatrix), intent(in) :: this
      type(FractionMatrix), intent(in) :: that
      logical                          :: output
    end function
  
    impure elemental module function non_equality_FractionMatrix_IntMatrix(this,that) result(output) 
      type(FractionMatrix), intent(in) :: this
      type(IntMatrix),      intent(in) :: that
      logical                          :: output
    end function
  
    impure elemental module function non_equality_IntMatrix_FractionMatrix(this,that) result(output) 
      type(IntMatrix),      intent(in) :: this
      type(FractionMatrix), intent(in) :: that
      logical                          :: output
    end function
  end interface
  
  interface transpose
    ! Matrix transpose.
    module function transpose_FractionMatrix(this) result(output) 
      type(FractionMatrix), intent(in) :: this
      type(FractionMatrix)             :: output
    end function
  end interface
  
  interface invert
    ! Inversion of a 3x3 integer matrix.
    module function invert_integers(input) result(output) 
      integer, intent(in) :: input(:,:)
      type(IntFraction)   :: output(3,3)
    end function
  
    module function invert_IntMatrix(input) result(output) 
      type(IntMatrix), intent(in) :: input
      type(FractionMatrix)        :: output
    end function
  end interface
  
  interface operator(+)
    ! Linear algebra.
    impure elemental module function add_FractionVector_FractionVector(this, &
       & that) result(output) 
      type(FractionVector), intent(in) :: this
      type(FractionVector), intent(in) :: that
      type(FractionVector)             :: output
    end function
  
    impure elemental module function add_FractionVector_IntVector(this,that) &
       & result(output) 
      type(IntVector),      intent(in) :: this
      type(FractionVector), intent(in) :: that
      type(FractionVector)             :: output
    end function
  
    impure elemental module function add_IntVector_FractionVector(this,that) &
       & result(output) 
      type(FractionVector), intent(in) :: this
      type(IntVector),      intent(in) :: that
      type(FractionVector)             :: output
    end function
  
    impure elemental module function add_FractionMatrix_FractionMatrix(this, &
       & that) result(output) 
      type(FractionMatrix), intent(in) :: this
      type(FractionMatrix), intent(in) :: that
      type(FractionMatrix)             :: output
    end function
  
    impure elemental module function add_FractionMatrix_IntMatrix(this,that) &
       & result(output) 
      type(IntMatrix),      intent(in) :: this
      type(FractionMatrix), intent(in) :: that
      type(FractionMatrix)             :: output
    end function
  
    impure elemental module function add_IntMatrix_FractionMatrix(this,that) &
       & result(output) 
      type(FractionMatrix), intent(in) :: this
      type(IntMatrix),      intent(in) :: that
      type(FractionMatrix)             :: output
    end function
  end interface
  
  interface operator(-)
    impure elemental module function negative_FractionVector(this) &
       & result(output) 
      type(FractionVector), intent(in) :: this
      type(FractionVector)             :: output
    end function
  
    impure elemental module function negative_FractionMatrix(this) &
       & result(output) 
      type(FractionMatrix), intent(in) :: this
      type(FractionMatrix)             :: output
    end function
  
    impure elemental module function subtract_FractionVector_FractionVector(this,that) result(output) 
      type(FractionVector), intent(in) :: this
      type(FractionVector), intent(in) :: that
      type(FractionVector)             :: output
    end function
  
    impure elemental module function subtract_FractionVector_IntVector(this, &
       & that) result(output) 
      type(FractionVector), intent(in) :: this
      type(IntVector),      intent(in) :: that
      type(FractionVector)             :: output
    end function
  
    impure elemental module function subtract_IntVector_FractionVector(this, &
       & that) result(output) 
      type(IntVector),      intent(in) :: this
      type(FractionVector), intent(in) :: that
      type(FractionVector)             :: output
    end function
  
    impure elemental module function subtract_FractionMatrix_FractionMatrix(this,that) result(output) 
      type(FractionMatrix), intent(in) :: this
      type(FractionMatrix), intent(in) :: that
      type(FractionMatrix)             :: output
    end function
  
    impure elemental module function subtract_FractionMatrix_IntMatrix(this, &
       & that) result(output) 
      type(FractionMatrix), intent(in) :: this
      type(IntMatrix),      intent(in) :: that
      type(FractionMatrix)             :: output
    end function
  
    impure elemental module function subtract_IntMatrix_FractionMatrix(this, &
       & that) result(output) 
      type(IntMatrix),      intent(in) :: this
      type(FractionMatrix), intent(in) :: that
      type(FractionMatrix)             :: output
    end function
  end interface
  
  interface operator(*)
    impure elemental module function multiply_FractionVector_integer(this, &
       & that) result(output) 
      type(FractionVector), intent(in) :: this
      integer,              intent(in) :: that
      type(FractionVector)             :: output
    end function
  
    impure elemental module function multiply_integer_FractionVector(this, &
       & that) result(output) 
      integer,              intent(in) :: this
      type(FractionVector), intent(in) :: that
      type(FractionVector)             :: output
    end function
  
    impure elemental module function multiply_FractionVector_IntFraction(this,that) result(output) 
      type(FractionVector), intent(in) :: this
      type(IntFraction),    intent(in) :: that
      type(FractionVector)             :: output
    end function
  
    impure elemental module function multiply_IntFraction_FractionVector(this,that) result(output) 
      type(IntFraction),    intent(in) :: this
      type(FractionVector), intent(in) :: that
      type(FractionVector)             :: output
    end function
  
    impure elemental module function multiply_FractionMatrix_integer(this, &
       & that) result(output) 
      type(FractionMatrix), intent(in) :: this
      integer,              intent(in) :: that
      type(FractionMatrix)             :: output
    end function
  
    impure elemental module function multiply_integer_FractionMatrix(this, &
       & that) result(output) 
      integer,              intent(in) :: this
      type(FractionMatrix), intent(in) :: that
      type(FractionMatrix)             :: output
    end function
  
    impure elemental module function multiply_FractionMatrix_IntFraction(this,that) result(output) 
      type(FractionMatrix), intent(in) :: this
      type(IntFraction),    intent(in) :: that
      type(FractionMatrix)             :: output
    end function
  
    impure elemental module function multiply_IntFraction_FractionMatrix(this,that) result(output) 
      type(IntFraction),    intent(in) :: this
      type(FractionMatrix), intent(in) :: that
      type(FractionMatrix)             :: output
    end function
  
    impure elemental module function dot_FractionVector_FractionVector(this, &
       & that) result(output) 
      type(FractionVector), intent(in) :: this
      type(FractionVector), intent(in) :: that
      type(IntFraction)                :: output
    end function
  
    impure elemental module function dot_FractionVector_IntVector(this,that) &
       & result(output) 
      type(FractionVector), intent(in) :: this
      type(IntVector),      intent(in) :: that
      type(IntFraction)                :: output
    end function
  
    impure elemental module function dot_IntVector_FractionVector(this,that) &
       & result(output) 
      type(IntVector),      intent(in) :: this
      type(FractionVector), intent(in) :: that
      type(IntFraction)                :: output
    end function
  
    impure elemental module function dot_FractionMatrix_FractionVector(this, &
       & that) result(output) 
      type(FractionMatrix), intent(in) :: this
      type(FractionVector), intent(in) :: that
      type(FractionVector)             :: output
    end function
  
    impure elemental module function dot_FractionMatrix_IntVector(this,that) &
       & result(output) 
      type(FractionMatrix), intent(in) :: this
      type(IntVector),      intent(in) :: that
      type(FractionVector)             :: output
    end function
  
    impure elemental module function dot_IntMatrix_FractionVector(this,that) &
       & result(output) 
      type(IntMatrix),      intent(in) :: this
      type(FractionVector), intent(in) :: that
      type(FractionVector)             :: output
    end function
  
    impure elemental module function dot_FractionVector_FractionMatrix(this, &
       & that) result(output) 
      type(FractionVector), intent(in) :: this
      type(FractionMatrix), intent(in) :: that
      type(FractionVector)             :: output
    end function
  
    impure elemental module function dot_FractionVector_IntMatrix(this,that) &
       & result(output) 
      type(FractionVector), intent(in) :: this
      type(IntMatrix),      intent(in) :: that
      type(FractionVector)             :: output
    end function
  
    impure elemental module function dot_IntVector_FractionMatrix(this,that) &
       & result(output) 
      type(IntVector),      intent(in) :: this
      type(FractionMatrix), intent(in) :: that
      type(FractionVector)             :: output
    end function
  
    impure elemental module function dot_FractionMatrix_FractionMatrix(this, &
       & that) result(output) 
      type(FractionMatrix), intent(in) :: this
      type(FractionMatrix), intent(in) :: that
      type(FractionMatrix)             :: output
    end function
  
    impure elemental module function dot_FractionMatrix_IntMatrix(this,that) &
       & result(output) 
      type(FractionMatrix), intent(in) :: this
      type(IntMatrix),      intent(in) :: that
      type(FractionMatrix)             :: output
    end function
  
    impure elemental module function dot_IntMatrix_FractionMatrix(this,that) &
       & result(output) 
      type(IntMatrix),      intent(in) :: this
      type(FractionMatrix), intent(in) :: that
      type(FractionMatrix)             :: output
    end function
  end interface
  
  interface operator(/)
    impure elemental module function divide_FractionVector_integer(this,that) &
       & result(output) 
      type(FractionVector), intent(in) :: this
      integer,              intent(in) :: that
      type(FractionVector)             :: output
    end function
  
    impure elemental module function divide_FractionVector_IntFraction(this, &
       & that) result(output) 
      type(FractionVector), intent(in) :: this
      type(IntFraction),    intent(in) :: that
      type(FractionVector)             :: output
    end function
  
    impure elemental module function divide_FractionMatrix_integer(this,that) &
       & result(output) 
      type(FractionMatrix), intent(in) :: this
      integer,              intent(in) :: that
      type(FractionMatrix)             :: output
    end function
  
    impure elemental module function divide_FractionMatrix_IntFraction(this, &
       & that) result(output) 
      type(FractionMatrix), intent(in) :: this
      type(IntFraction),    intent(in) :: that
      type(FractionMatrix)             :: output
    end function
  end interface
  
  interface sum
    ! ----------------------------------------------------------------------
    ! Sum.
    ! ----------------------------------------------------------------------
    module function sum_FractionVector(input) result(output) 
      type(FractionVector), intent(in) :: input(:)
      type(FractionVector)             :: output
    end function
  
    module function sum_FractionMatrix(input) result(output) 
      type(FractionMatrix), intent(in) :: input(:)
      type(FractionMatrix)             :: output
    end function
  end interface
  
  interface exp_2pii
    ! ----------------------------------------------------------------------
    ! exp(2*pi*i*input), cos(2*pi*input) and sin(2*pi*input).
    ! ----------------------------------------------------------------------
    impure elemental module function exp_2pii_IntFraction(input) &
       & result(output) 
      type(IntFraction), intent(in) :: input
      complex(dp)                   :: output
    end function
  end interface
  
  interface cos_2pi
    impure elemental module function cos_2pi_IntFraction(input) result(output) 
      type(IntFraction), intent(in) :: input
      real(dp)                      :: output
    end function
  end interface
  
  interface sin_2pi
    impure elemental module function sin_2pi_IntFraction(input) result(output) 
      type(IntFraction), intent(in) :: input
      real(dp)                      :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_FractionVector(this,input) 
      class(FractionVector), intent(out) :: this
      type(String),          intent(in)  :: input
    end subroutine
  end interface
  
  interface
    module function write_FractionVector(this) result(output) 
      class(FractionVector), intent(in) :: this
      type(String)                      :: output
    end function
  end interface
  
  interface FractionVector
    impure elemental module function new_FractionVector_String(input) &
       & result(this) 
      type(String), intent(in) :: input
      type(FractionVector)     :: this
    end function
  end interface
  
  interface
    module subroutine read_FractionMatrix(this,input) 
      class(FractionMatrix), intent(out) :: this
      type(String),          intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_FractionMatrix(this) result(output) 
      class(FractionMatrix), intent(in) :: this
      type(String), allocatable         :: output(:)
    end function
  end interface
  
  interface FractionMatrix
    module function new_FractionMatrix_Strings(input) result(this) 
      type(String), intent(in) :: input(:)
      type(FractionMatrix)     :: this
    end function
  
    impure elemental module function new_FractionMatrix_StringArray(input) &
       & result(this) 
      type(StringArray), intent(in) :: input
      type(FractionMatrix)          :: this
    end function
  end interface
end module
