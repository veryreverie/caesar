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
  
  interface FractionVector
    module procedure new_FractionVector_String
  end interface
  
  type, extends(Stringsable) :: FractionMatrix
    type(IntFraction), allocatable, private :: contents_(:,:)
  contains
    procedure, public :: read  => read_FractionMatrix
    procedure, public :: write => write_FractionMatrix
  end type
  
  interface FractionMatrix
    module procedure new_FractionMatrix_Strings
    module procedure new_FractionMatrix_StringArray
  end interface
  
  ! Conversions to and from vector and matrix types.
  interface vec
    module procedure vec_IntFractions
  end interface
  
  interface mat
    module procedure mat_IntFractions
    module procedure mat_IntFractions_shape
  end interface
  
  interface frac
    module procedure frac_FractionVector
    module procedure frac_FractionMatrix
  end interface
  
  interface fracvec
    module procedure fracvec_IntVector
  end interface
  
  interface fracmat
    module procedure fracmat_IntMatrix
  end interface
  
  interface intvec
    module procedure intvec_FractionVector
  end interface
  
  interface intmat
    module procedure intmat_FractionMatrix
  end interface
  
  interface dblevec
    module procedure dblevec_FractionVector
  end interface
  
  interface dblemat
    module procedure dblemat_FractionMatrix
  end interface
  
  ! Properties of the vectors and matrices.
  interface size
    module procedure size_FractionVector
    module procedure size_FractionMatrix
  end interface
  
  interface is_int
    module procedure is_int_FractionVector
    module procedure is_int_FractionMatrix
  end interface
  
  ! Comparison between vectors and matrices.
  interface operator(==)
    module procedure equality_FractionVector_FractionVector
    module procedure equality_FractionVector_IntVector
    module procedure equality_IntVector_FractionVector
    
    module procedure equality_FractionMatrix_FractionMatrix
    module procedure equality_FractionMatrix_IntMatrix
    module procedure equality_IntMatrix_FractionMatrix
  end interface
  
  interface operator(/=)
    module procedure non_equality_FractionVector_FractionVector
    module procedure non_equality_FractionVector_IntVector
    module procedure non_equality_IntVector_FractionVector
    
    module procedure non_equality_FractionMatrix_FractionMatrix
    module procedure non_equality_FractionMatrix_IntMatrix
    module procedure non_equality_IntMatrix_FractionMatrix
  end interface
  
  ! Matrix transpose.
  interface transpose
    module procedure transpose_FractionMatrix
  end interface
  
  ! Exact inversion of an integer matrix.
  interface invert
    module procedure invert_integers
    module procedure invert_IntMatrix
  end interface
  
  ! Linear algebra.
  interface operator(+)
    module procedure add_FractionVector_FractionVector
    module procedure add_FractionVector_IntVector
    module procedure add_IntVector_FractionVector
    
    module procedure add_FractionMatrix_FractionMatrix
    module procedure add_FractionMatrix_IntMatrix
    module procedure add_IntMatrix_FractionMatrix
  end interface
  
  interface operator(-)
    module procedure negative_FractionVector
    module procedure negative_FractionMatrix
    
    module procedure subtract_FractionVector_FractionVector
    module procedure subtract_FractionVector_IntVector
    module procedure subtract_IntVector_FractionVector
    
    module procedure subtract_FractionMatrix_FractionMatrix
    module procedure subtract_FractionMatrix_IntMatrix
    module procedure subtract_IntMatrix_FractionMatrix
  end interface
  
  interface operator(*)
    module procedure multiply_FractionVector_integer
    module procedure multiply_integer_FractionVector
    module procedure multiply_FractionVector_IntFraction
    module procedure multiply_IntFraction_FractionVector
    
    module procedure multiply_FractionMatrix_integer
    module procedure multiply_integer_FractionMatrix
    module procedure multiply_FractionMatrix_IntFraction
    module procedure multiply_IntFraction_FractionMatrix
    
    module procedure dot_FractionVector_FractionVector
    module procedure dot_FractionVector_IntVector
    module procedure dot_IntVector_FractionVector
    
    module procedure dot_FractionMatrix_FractionVector
    module procedure dot_FractionMatrix_IntVector
    module procedure dot_IntMatrix_FractionVector
    
    module procedure dot_FractionVector_FractionMatrix
    module procedure dot_FractionVector_IntMatrix
    module procedure dot_IntVector_FractionMatrix
    
    module procedure dot_FractionMatrix_FractionMatrix
    module procedure dot_FractionMatrix_IntMatrix
    module procedure dot_IntMatrix_FractionMatrix
  end interface
  
  interface operator(/)
    module procedure divide_FractionVector_integer
    module procedure divide_FractionVector_IntFraction
    module procedure divide_FractionMatrix_integer
    module procedure divide_FractionMatrix_IntFraction
  end interface
  
  ! Sum.
  interface sum
    module procedure sum_FractionVector
    module procedure sum_FractionMatrix
  end interface
  
  ! exp(2*pi*i*input), cos(2*pi*input) and sin(2*pi*input)
  interface exp_2pii
    module procedure exp_2pii_IntFraction
  end interface
  
  interface cos_2pi
    module procedure cos_2pi_IntFraction
  end interface
  
  interface sin_2pi
    module procedure sin_2pi_IntFraction
  end interface
contains

! ----------------------------------------------------------------------
! Procedures involving contents_
! ----------------------------------------------------------------------

! Conversion to and from vector and matrix types.
function vec_IntFractions(input) result(output)
  implicit none
  
  type(IntFraction), intent(in) :: input(:)
  type(FractionVector)          :: output
  
  output%contents_ = input
end function

function mat_IntFractions(input) result(output)
  implicit none
  
  type(IntFraction), intent(in) :: input(:,:)
  type(FractionMatrix)          :: output
  
  output%contents_ = input
end function

function mat_IntFractions_shape(input,m,n) result(output)
  implicit none
  
  type(IntFraction), intent(in) :: input(:)
  integer,           intent(in) :: m
  integer,           intent(in) :: n
  type(FractionMatrix)          :: output
  
  type(IntFraction), allocatable :: contents(:,:)
  
  !output%contents_ = transpose(reshape(input, [m,n]))
  ! WORKAROUND to avoid internal compiler error in ifort 19.1.0.166.
  contents = reshape(input, [m,n])
  output = mat(transpose(contents))
end function

! Conversion to fraction(:). Effectively getters for contents_.
function frac_FractionVector(input) result(output)
  implicit none
  
  type(FractionVector), intent(in) :: input
  type(IntFraction), allocatable   :: output(:)
  
  if (allocated(input%contents_)) then
    output = input%contents_
  else
    call print_line(CODE_ERROR//': Trying to use the contents of a vector &
       &before it has been allocated.')
    call err()
  endif
end function

function frac_FractionMatrix(input) result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: input
  type(IntFraction), allocatable   :: output(:,:)
  
  if (allocated(input%contents_)) then
    output = input%contents_
  else
    call print_line(CODE_ERROR//': Trying to use the contents of a matrix &
       &before it has been allocated.')
    call err()
  endif
end function

! ----------------------------------------------------------------------
! Procedures not involving contents_
! ----------------------------------------------------------------------
! N.B. the number of procedures accessing contents_ directly is intentionally
!    limited for stability reasons.
! The above procedures behave well if contents_ has not been allocated,
!    and this good behaviour is automatically passed to the procedures below.

! Conversion to and from vector and matrix types.

impure elemental function fracvec_IntVector(input) result(output)
  implicit none
  
  type(IntVector), intent(in) :: input
  type(FractionVector)        :: output
  
  output = vec(frac(int(input)))
end function

impure elemental function fracmat_IntMatrix(input) result(output)
  implicit none
  
  type(IntMatrix), intent(in)    :: input
  type(FractionMatrix)           :: output
  
  output = mat(frac(int(input)))
end function

impure elemental function intvec_FractionVector(input) result(output)
  implicit none
  
  type(FractionVector), intent(in) :: input
  type(IntVector)                  :: output
  
  output = vec(int(frac(input)))
end function

impure elemental function intmat_FractionMatrix(input) result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: input
  type(IntMatrix)                  :: output
  
  output = mat(int(frac(input)))
end function

impure elemental function dblevec_FractionVector(input) result(output)
  implicit none
  
  type(FractionVector), intent(in) :: input
  type(RealVector)                 :: output
  
  output = vec(dble(frac(input)))
end function

impure elemental function dblemat_FractionMatrix(input) result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: input
  type(RealMatrix)                 :: output
  
  output = mat(dble(frac(input)))
end function

! Properties of the vectors and matrices.
function size_FractionVector(this) result(output)
  implicit none
  
  type(FractionVector), intent(in) :: this
  integer                          :: output
  
  output = size(frac(this))
end function

function size_FractionMatrix(this,dim) result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  integer,              intent(in) :: dim
  integer                          :: output
  
  output = size(frac(this), dim)
end function

impure elemental function is_int_FractionVector(this) result(output)
  implicit none
  
  type(FractionVector), intent(in) :: this
  logical                          :: output
  
  output = all(is_int(frac(this)))
end function

impure elemental function is_int_FractionMatrix(this) result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  logical                          :: output
  
  output = all(is_int(frac(this)))
end function

! Comparisons.
impure elemental function equality_FractionVector_FractionVector(this,that) &
   & result(output)
  implicit none
  
  type(FractionVector), intent(in) :: this
  type(FractionVector), intent(in) :: that
  logical                          :: output
  
  output = all(frac(this)==frac(that))
end function

impure elemental function equality_FractionVector_IntVector(this,that) result(output)
  implicit none
  
  type(FractionVector), intent(in) :: this
  type(IntVector),      intent(in) :: that
  logical                          :: output
  
  output = all(frac(this)==int(that))
end function

impure elemental function equality_IntVector_FractionVector(this,that) result(output)
  implicit none
  
  type(IntVector),      intent(in) :: this
  type(FractionVector), intent(in) :: that
  logical                          :: output
  
  output = all(int(this)==frac(that))
end function

impure elemental function equality_FractionMatrix_FractionMatrix(this,that) &
   & result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  type(FractionMatrix), intent(in) :: that
  logical                          :: output
  
  output = all(frac(this)==frac(that))
end function

impure elemental function equality_FractionMatrix_IntMatrix(this,that) result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  type(IntMatrix),      intent(in) :: that
  logical                          :: output
  
  output = all(frac(this)==int(that))
end function

impure elemental function equality_IntMatrix_FractionMatrix(this,that) result(output)
  implicit none
  
  type(IntMatrix),      intent(in) :: this
  type(FractionMatrix), intent(in) :: that
  logical                          :: output
  
  output = all(int(this)==frac(that))
end function

impure elemental function non_equality_FractionVector_FractionVector(this,that) &
   & result(output)
  implicit none
  
  type(FractionVector), intent(in) :: this
  type(FractionVector), intent(in) :: that
  logical                          :: output
  
  output = .not. this==that
end function

impure elemental function non_equality_FractionVector_IntVector(this,that) &
   & result(output)
  implicit none
  
  type(FractionVector), intent(in) :: this
  type(IntVector),      intent(in) :: that
  logical                          :: output
  
  output = .not. this==that
end function

impure elemental function non_equality_IntVector_FractionVector(this,that) &
   & result(output)
  implicit none
  
  type(IntVector),      intent(in) :: this
  type(FractionVector), intent(in) :: that
  logical                          :: output
  
  output = .not. this==that
end function

impure elemental function non_equality_FractionMatrix_FractionMatrix(this,that) &
   & result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  type(FractionMatrix), intent(in) :: that
  logical                          :: output
  
  output = .not. this==that
end function

impure elemental function non_equality_FractionMatrix_IntMatrix(this,that) &
   & result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  type(IntMatrix),      intent(in) :: that
  logical                          :: output
  
  output = .not. this==that
end function

impure elemental function non_equality_IntMatrix_FractionMatrix(this,that) &
   & result(output)
  implicit none
  
  type(IntMatrix),      intent(in) :: this
  type(FractionMatrix), intent(in) :: that
  logical                          :: output
  
  output = .not. this==that
end function

! Matrix transpose.
function transpose_FractionMatrix(this) result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  type(FractionMatrix)             :: output
  
  output = mat(transpose(frac(this)))
end function

! Inversion of a 3x3 integer matrix.
function invert_integers(input) result(output)
  implicit none
  
  integer, intent(in) :: input(:,:)
  type(IntFraction)   :: output(3,3)
  
  integer :: det
  integer :: i,j
  
  ! Check that the input is a 3x3 matrix.
  if (size(input,1)/=3 .or. size(input,2)/=3) then
    call print_line(CODE_ERROR//': Trying to invert matrix which is not 3x3.')
    call err()
  endif
  
  det = determinant(input)
  if (det==0) then
    call print_line(ERROR//': Trying to invert a Matrix with determinant=0.')
    call err()
  endif
 
  output(1,1) = frac(input(2,2)*input(3,3)-input(3,2)*input(2,3))
  output(1,2) = frac(input(3,2)*input(1,3)-input(1,2)*input(3,3))
  output(1,3) = frac(input(1,2)*input(2,3)-input(2,2)*input(1,3))
  output(2,1) = frac(input(2,3)*input(3,1)-input(3,3)*input(2,1))
  output(2,2) = frac(input(3,3)*input(1,1)-input(1,3)*input(3,1))
  output(2,3) = frac(input(1,3)*input(2,1)-input(2,3)*input(1,1))
  output(3,1) = frac(input(2,1)*input(3,2)-input(3,1)*input(2,2))
  output(3,2) = frac(input(3,1)*input(1,2)-input(1,1)*input(3,2))
  output(3,3) = frac(input(1,1)*input(2,2)-input(2,1)*input(1,2))
  
  do i=1,3
    do j=1,3
      output(j,i) = output(j,i)/det
    enddo
  enddo
end function

function invert_IntMatrix(input) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: input
  type(FractionMatrix)        :: output
  
  output = mat(invert(int(input)))
end function

! Linear algebra.
impure elemental function add_FractionVector_FractionVector(this,that) &
   & result(output)
  implicit none
  
  type(FractionVector), intent(in) :: this
  type(FractionVector), intent(in) :: that
  type(FractionVector)             :: output
  
  output = vec(frac(this) + frac(that))
end function

impure elemental function add_FractionVector_IntVector(this,that) &
   & result(output)
  implicit none
  
  type(IntVector),      intent(in) :: this
  type(FractionVector), intent(in) :: that
  type(FractionVector)             :: output
  
  output = vec(int(this) + frac(that))
end function

impure elemental function add_IntVector_FractionVector(this,that) &
   & result(output)
  implicit none
  
  type(FractionVector), intent(in) :: this
  type(IntVector),      intent(in) :: that
  type(FractionVector)             :: output
  
  output = vec(frac(this) + int(that))
end function

impure elemental function add_FractionMatrix_FractionMatrix(this,that) &
   & result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  type(FractionMatrix), intent(in) :: that
  type(FractionMatrix)             :: output
  
  output = mat(frac(this) + frac(that))
end function

impure elemental function add_FractionMatrix_IntMatrix(this,that) &
   & result(output)
  implicit none
  
  type(IntMatrix),      intent(in) :: this
  type(FractionMatrix), intent(in) :: that
  type(FractionMatrix)             :: output
  
  output = mat(int(this) + frac(that))
end function

impure elemental function add_IntMatrix_FractionMatrix(this,that) &
   & result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  type(IntMatrix),      intent(in) :: that
  type(FractionMatrix)             :: output
  
  output = mat(frac(this) + int(that))
end function

impure elemental function negative_FractionVector(this) result(output)
  implicit none
  
  type(FractionVector), intent(in) :: this
  type(FractionVector)             :: output
  
  output = vec(-frac(this))
end function

impure elemental function negative_FractionMatrix(this) result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  type(FractionMatrix)             :: output
  
  output = mat(-frac(this))
end function

impure elemental function subtract_FractionVector_FractionVector(this,that) &
   & result(output)
  implicit none
  
  type(FractionVector), intent(in) :: this
  type(FractionVector), intent(in) :: that
  type(FractionVector)             :: output
  
  output = vec(frac(this) - frac(that))
end function

impure elemental function subtract_FractionVector_IntVector(this,that) &
   & result(output)
  implicit none
  
  type(FractionVector), intent(in) :: this
  type(IntVector),      intent(in) :: that
  type(FractionVector)             :: output
  
  output = vec(frac(this) - int(that))
end function

impure elemental function subtract_IntVector_FractionVector(this,that) &
   & result(output)
  implicit none
  
  type(IntVector),      intent(in) :: this
  type(FractionVector), intent(in) :: that
  type(FractionVector)             :: output
  
  output = vec(int(this) - frac(that))
end function

impure elemental function subtract_FractionMatrix_FractionMatrix(this,that) &
   & result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  type(FractionMatrix), intent(in) :: that
  type(FractionMatrix)             :: output
  
  output = mat(frac(this) - frac(that))
end function

impure elemental function subtract_FractionMatrix_IntMatrix(this,that) &
   & result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  type(IntMatrix),      intent(in) :: that
  type(FractionMatrix)             :: output
  
  output = mat(frac(this) - int(that))
end function

impure elemental function subtract_IntMatrix_FractionMatrix(this,that) &
   & result(output)
  implicit none
  
  type(IntMatrix),      intent(in) :: this
  type(FractionMatrix), intent(in) :: that
  type(FractionMatrix)             :: output
  
  output = mat(int(this) - frac(that))
end function

impure elemental function multiply_FractionVector_integer(this,that) &
   & result(output)
  implicit none
  
  type(FractionVector), intent(in) :: this
  integer,              intent(in) :: that
  type(FractionVector)             :: output
  
  output = vec(frac(this) * that)
end function

impure elemental function multiply_integer_FractionVector(this,that) &
   & result(output)
  implicit none
  
  integer,              intent(in) :: this
  type(FractionVector), intent(in) :: that
  type(FractionVector)             :: output
  
  output = vec(this * frac(that))
end function

impure elemental function multiply_FractionVector_IntFraction(this,that) &
   & result(output)
  implicit none
  
  type(FractionVector), intent(in) :: this
  type(IntFraction),    intent(in) :: that
  type(FractionVector)             :: output
  
  output = vec(frac(this) * that)
end function

impure elemental function multiply_IntFraction_FractionVector(this,that) &
   & result(output)
  implicit none
  
  type(IntFraction),    intent(in) :: this
  type(FractionVector), intent(in) :: that
  type(FractionVector)             :: output
  
  output = vec(this * frac(that))
end function

impure elemental function multiply_FractionMatrix_integer(this,that) &
   & result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  integer,              intent(in) :: that
  type(FractionMatrix)             :: output
  
  output = mat(frac(this) * that)
end function

impure elemental function multiply_integer_FractionMatrix(this,that) &
   & result(output)
  implicit none
  
  integer,              intent(in) :: this
  type(FractionMatrix), intent(in) :: that
  type(FractionMatrix)             :: output
  
  output = mat(this * frac(that))
end function

impure elemental function multiply_FractionMatrix_IntFraction(this,that) &
   & result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  type(IntFraction),    intent(in) :: that
  type(FractionMatrix)             :: output
  
  output = mat(frac(this) * that)
end function

impure elemental function multiply_IntFraction_FractionMatrix(this,that) &
   & result(output)
  implicit none
  
  type(IntFraction),    intent(in) :: this
  type(FractionMatrix), intent(in) :: that
  type(FractionMatrix)             :: output
  
  output = mat(this * frac(that))
end function

impure elemental function dot_FractionVector_FractionVector(this,that) &
   & result(output)
  implicit none
  
  type(FractionVector), intent(in) :: this
  type(FractionVector), intent(in) :: that
  type(IntFraction)                :: output
  
  type(IntFraction), allocatable :: a(:)
  type(IntFraction), allocatable :: b(:)
  
  integer :: i
  
  if (size(this)/=size(that)) then
    call print_line( CODE_ERROR//': Dot product of vectors of different &
       &lengths.')
    call err()
  endif
  
  a = frac(this)
  b = frac(that)
  output = frac(0)
  do i=1,size(this)
    output = output + a(i)*b(i)
  enddo
end function

impure elemental function dot_FractionVector_IntVector(this,that) &
   & result(output)
  implicit none
  
  type(FractionVector), intent(in) :: this
  type(IntVector),      intent(in) :: that
  type(IntFraction)                :: output
  
  type(IntFraction), allocatable :: a(:)
  integer,           allocatable :: b(:)
  
  integer :: i
  
  if (size(this)/=size(that)) then
    call print_line( CODE_ERROR//': Dot product of vectors of different &
       &lengths.')
    call err()
  endif
  
  a = frac(this)
  b = int(that)
  output = frac(0)
  do i=1,size(this)
    output = output + a(i)*b(i)
  enddo
end function

impure elemental function dot_IntVector_FractionVector(this,that) &
   & result(output)
  implicit none
  
  type(IntVector),      intent(in) :: this
  type(FractionVector), intent(in) :: that
  type(IntFraction)                :: output
  
  integer,           allocatable :: a(:)
  type(IntFraction), allocatable :: b(:)
  
  integer :: i
  
  if (size(this)/=size(that)) then
    call print_line( CODE_ERROR//': Dot product of vectors of different &
       &lengths.')
    call err()
  endif
  
  a = int(this)
  b = frac(that)
  output = frac(0)
  do i=1,size(this)
    output = output + a(i)*b(i)
  enddo
end function

impure elemental function dot_FractionMatrix_FractionVector(this,that) &
   & result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  type(FractionVector), intent(in) :: that
  type(FractionVector)             :: output
  
  type(IntFraction), allocatable :: a(:,:)
  type(IntFraction), allocatable :: b(:)
  type(IntFraction), allocatable :: contents(:)
  
  integer :: i,ialloc
  
  if (size(this,2)/=size(that)) then
    call print_line(CODE_ERROR//': Dot product of matrix and vector of &
       &incompatible dimensions.')
    call err()
  endif
  
  a = frac(this)
  b = frac(that)
  allocate(contents(size(this,1)), stat=ialloc); call err(ialloc)
  contents = frac(0)
  do i=1,size(that)
    contents = contents + a(:,i)*b(i)
  enddo
  output = vec(contents)
end function

impure elemental function dot_FractionMatrix_IntVector(this,that) &
   & result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  type(IntVector),      intent(in) :: that
  type(FractionVector)             :: output
  
  type(IntFraction), allocatable :: a(:,:)
  integer,           allocatable :: b(:)
  type(IntFraction), allocatable :: contents(:)
  
  integer :: i,ialloc
  
  if (size(this,2)/=size(that)) then
    call print_line(CODE_ERROR//': Dot product of matrix and vector of &
       &incompatible dimensions.')
    call err()
  endif
  
  a = frac(this)
  b = int(that)
  allocate(contents(size(this,1)), stat=ialloc); call err(ialloc)
  contents = frac(0)
  do i=1,size(that)
    contents = contents + a(:,i)*b(i)
  enddo
  output = vec(contents)
end function

impure elemental function dot_IntMatrix_FractionVector(this,that) &
   & result(output)
  implicit none
  
  type(IntMatrix),      intent(in) :: this
  type(FractionVector), intent(in) :: that
  type(FractionVector)             :: output
  
  integer,           allocatable :: a(:,:)
  type(IntFraction), allocatable :: b(:)
  type(IntFraction), allocatable :: contents(:)
  
  integer :: i,ialloc
  
  if (size(this,2)/=size(that)) then
    call print_line(CODE_ERROR//': Dot product of matrix and vector of &
       &incompatible dimensions.')
    call err()
  endif
  
  a = int(this)
  b = frac(that)
  allocate(contents(size(this,1)), stat=ialloc); call err(ialloc)
  contents = frac(0)
  do i=1,size(that)
    contents = contents + a(:,i)*b(i)
  enddo
  output = vec(contents)
end function

impure elemental function dot_FractionVector_FractionMatrix(this,that) &
   & result(output)
  implicit none
  
  type(FractionVector), intent(in) :: this
  type(FractionMatrix), intent(in) :: that
  type(FractionVector)             :: output
  
  type(IntFraction), allocatable :: a(:)
  type(IntFraction), allocatable :: b(:,:)
  type(IntFraction), allocatable :: contents(:)
  
  integer :: i,ialloc
  
  if (size(this)/=size(that,1)) then
    call print_line(CODE_ERROR//': Dot product of vector and matrix of &
       &incompatible dimensions.')
    call err()
  endif
  
  a = frac(this)
  b = frac(that)
  allocate(contents(size(that,2)), stat=ialloc); call err(ialloc)
  contents = frac(0)
  do i=1,size(this)
    contents = contents + a(i)*b(i,:)
  enddo
  output = vec(contents)
end function

impure elemental function dot_FractionVector_IntMatrix(this,that) &
   & result(output)
  implicit none
  
  type(FractionVector), intent(in) :: this
  type(IntMatrix),      intent(in) :: that
  type(FractionVector)             :: output
  
  type(IntFraction), allocatable :: a(:)
  integer,           allocatable :: b(:,:)
  type(IntFraction), allocatable :: contents(:)
  
  integer :: i,ialloc
  
  if (size(this)/=size(that,1)) then
    call print_line(CODE_ERROR//': Dot product of vector and matrix of &
       &incompatible dimensions.')
    call err()
  endif
  
  a = frac(this)
  b = int(that)
  allocate(contents(size(that,2)), stat=ialloc); call err(ialloc)
  contents = frac(0)
  do i=1,size(this)
    contents = contents + a(i)*b(i,:)
  enddo
  output = vec(contents)
end function

impure elemental function dot_IntVector_FractionMatrix(this,that) &
   & result(output)
  implicit none
  
  type(IntVector),      intent(in) :: this
  type(FractionMatrix), intent(in) :: that
  type(FractionVector)             :: output
  
  integer,           allocatable :: a(:)
  type(IntFraction), allocatable :: b(:,:)
  type(IntFraction), allocatable :: contents(:)
  
  integer :: i,ialloc
  
  if (size(this)/=size(that,1)) then
    call print_line(CODE_ERROR//': Dot product of vector and matrix of &
       &incompatible dimensions.')
    call err()
  endif
  
  a = int(this)
  b = frac(that)
  allocate(contents(size(that,2)), stat=ialloc); call err(ialloc)
  contents = frac(0)
  do i=1,size(this)
    contents = contents + a(i)*b(i,:)
  enddo
  output = vec(contents)
end function

impure elemental function dot_FractionMatrix_FractionMatrix(this,that) &
   & result(output) 
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  type(FractionMatrix), intent(in) :: that
  type(FractionMatrix)             :: output
  
  type(IntFraction), allocatable :: a(:,:)
  type(IntFraction), allocatable :: b(:,:)
  type(IntFraction), allocatable :: contents(:,:)
  
  integer :: i,j,k,ialloc
  
  if (size(this,2)/=size(that,1)) then
    call print_line(CODE_ERROR//': Dot product of two matrices of &
       &incompatible dimensions.')
    call err()
  endif
  
  a = frac(this)
  b = frac(that)
  allocate(contents(size(this,1),size(that,2)), stat=ialloc); call err(ialloc)
  contents = frac(0)
  do i=1,size(that,2)
    do j=1,size(that,1)
      do k=1,size(this,1)
        contents(k,i) = contents(k,i) + a(k,j)*b(j,i)
      enddo
    enddo
  enddo
  output = mat(contents)
end function

impure elemental function dot_FractionMatrix_IntMatrix(this,that) &
   & result(output) 
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  type(IntMatrix),      intent(in) :: that
  type(FractionMatrix)             :: output
  
  type(IntFraction), allocatable :: a(:,:)
  integer,           allocatable :: b(:,:)
  type(IntFraction), allocatable :: contents(:,:)
  
  integer :: i,j,k,ialloc
  
  if (size(this,2)/=size(that,1)) then
    call print_line(CODE_ERROR//': Dot product of two matrices of &
       &incompatible dimensions.')
    call err()
  endif
  
  a = frac(this)
  b = int(that)
  allocate(contents(size(this,1),size(that,2)), stat=ialloc); call err(ialloc)
  contents = frac(0)
  do i=1,size(that,2)
    do j=1,size(that,1)
      do k=1,size(this,1)
        contents(k,i) = contents(k,i) + a(k,j)*b(j,i)
      enddo
    enddo
  enddo
  output = mat(contents)
end function

impure elemental function dot_IntMatrix_FractionMatrix(this,that) &
   & result(output) 
  implicit none
  
  type(IntMatrix),      intent(in) :: this
  type(FractionMatrix), intent(in) :: that
  type(FractionMatrix)             :: output
  
  integer,           allocatable :: a(:,:)
  type(IntFraction), allocatable :: b(:,:)
  type(IntFraction), allocatable :: contents(:,:)
  
  integer :: i,j,k,ialloc
  
  if (size(this,2)/=size(that,1)) then
    call print_line(CODE_ERROR//': Dot product of two matrices of &
       &incompatible dimensions.')
    call err()
  endif
  
  a = int(this)
  b = frac(that)
  allocate(contents(size(this,1),size(that,2)), stat=ialloc); call err(ialloc)
  contents = frac(0)
  do i=1,size(that,2)
    do j=1,size(that,1)
      do k=1,size(this,1)
        contents(k,i) = contents(k,i) + a(k,j)*b(j,i)
      enddo
    enddo
  enddo
  output = mat(contents)
end function

impure elemental function divide_FractionVector_integer(this,that) &
   & result(output) 
  implicit none
  
  type(FractionVector), intent(in) :: this
  integer,              intent(in) :: that
  type(FractionVector)             :: output
  
  output = vec(frac(this) / that)
end function

impure elemental function divide_FractionVector_IntFraction(this,that) &
   & result(output) 
  implicit none
  
  type(FractionVector), intent(in) :: this
  type(IntFraction),    intent(in) :: that
  type(FractionVector)             :: output
  
  output = vec(frac(this) / that)
end function

impure elemental function divide_FractionMatrix_integer(this,that) &
   & result(output) 
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  integer,              intent(in) :: that
  type(FractionMatrix)             :: output
  
  output = mat(frac(this) / that)
end function

impure elemental function divide_FractionMatrix_IntFraction(this,that) &
   & result(output) 
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  type(IntFraction),    intent(in) :: that
  type(FractionMatrix)             :: output
  
  output = mat(frac(this) / that)
end function

! ----------------------------------------------------------------------
! Sum.
! ----------------------------------------------------------------------
function sum_FractionVector(input) result(output)
  implicit none
  
  type(FractionVector), intent(in) :: input(:)
  type(FractionVector)             :: output
  
  integer :: i
  
  if (size(input)==0) then
    call print_line(CODE_ERROR//': Trying to take the sum of an empty array.')
    call err()
  endif
  
  output = input(1)
  do i=2,size(input)
    output = output+input(i)
  enddo
end function

function sum_FractionMatrix(input) result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: input(:)
  type(FractionMatrix)             :: output
  
  integer :: i
  
  if (size(input)==0) then
    call print_line(CODE_ERROR//': Trying to take the sum of an empty array.')
    call err()
  endif
  
  output = input(1)
  do i=2,size(input)
    output = output+input(i)
  enddo
end function

! ----------------------------------------------------------------------
! exp(2*pi*i*input), cos(2*pi*input) and sin(2*pi*input).
! ----------------------------------------------------------------------
impure elemental function exp_2pii_IntFraction(input) result(output)
  implicit none
  
  type(IntFraction), intent(in) :: input
  complex(dp)                   :: output
  
  output = exp_2pii(dble(input))
end function

impure elemental function cos_2pi_IntFraction(input) result(output)
  implicit none
  
  type(IntFraction), intent(in) :: input
  real(dp)                      :: output
  
  output = cos_2pi(dble(input))
end function

impure elemental function sin_2pi_IntFraction(input) result(output)
  implicit none
  
  type(IntFraction), intent(in) :: input
  real(dp)                      :: output
  
  output = sin_2pi(dble(input))
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_FractionVector(this,input)
  implicit none
  
  class(FractionVector), intent(out) :: this
  type(String),          intent(in)  :: input
  
  select type(this); type is(FractionVector)
    this = vec(frac(split_line(input)))
  end select
end subroutine

function write_FractionVector(this) result(output)
  implicit none
  
  class(FractionVector), intent(in) :: this
  type(String)                      :: output
  
  select type(this); type is(FractionVector)
    output = join(frac(this))
  end select
end function

impure elemental function new_FractionVector_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(FractionVector)     :: this
  
  call this%read(input)
end function

subroutine read_FractionMatrix(this,input)
  implicit none
  
  class(FractionMatrix), intent(out) :: this
  type(String),          intent(in)  :: input(:)
  
  type(IntFraction), allocatable :: line(:)
  type(IntFraction), allocatable :: contents(:,:)
  
  integer :: i,ialloc
  
  select type(this); type is(FractionMatrix)
    if (size(input)==0) then
      allocate(contents(0,0), stat=ialloc); call err(ialloc)
    else
      line = frac(split_line(input(1)))
      allocate( contents(size(input),size(line)), &
              & stat=ialloc); call err(ialloc)
      contents(1,:) = line
      do i=2,size(input)
        line = frac(split_line(input(i)))
        if (size(line)/=size(contents,2)) then
          call print_line(ERROR//': Reading matrix: rows of different &
             &lengths.')
          call err()
        endif
        contents(i,:) = line
      enddo
    endif
    
    this = mat(contents)
  class default
    call err()
  end select
end subroutine

function write_FractionMatrix(this) result(output)
  implicit none
  
  class(FractionMatrix), intent(in) :: this
  type(String), allocatable         :: output(:)
  
  type(IntFraction), allocatable :: contents(:,:)
  type(IntFraction), allocatable :: row(:)
  
  integer :: i,ialloc
  
  select type(this); type is(FractionMatrix)
    contents = frac(this)
    allocate(output(size(this,1)), stat=ialloc); call err(ialloc)
    do i=1,size(this,1)
      row = contents(i,:)
      output(i) = join(row)
    enddo
  class default
    call err()
  end select
end function

function new_FractionMatrix_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(FractionMatrix)     :: this
  
  call this%read(input)
end function

impure elemental function new_FractionMatrix_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(FractionMatrix)          :: this
  
  this = FractionMatrix(str(input))
end function
end module
