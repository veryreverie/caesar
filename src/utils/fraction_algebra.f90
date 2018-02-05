! ======================================================================
! Vectors and Matrices of type IntFraction.
! ======================================================================
module fraction_algebra_module
  use constants_module, only : dp
  use string_module
  use io_module
  
  use fraction_module
  use linear_algebra_module
  implicit none
  
  type, extends(Stringable) :: FractionVector
    type(IntFraction), allocatable, private :: contents_(:)
  contains
    generic,   public  :: assignment(=) => assign_FractionVector_IntFractions
    procedure, private ::                  assign_FractionVector_IntFractions
    
    procedure :: str => str_FractionVector
  end type
  
  type, extends(Printable) :: FractionMatrix
    type(IntFraction), allocatable, private :: contents_(:,:)
  contains
    generic,   public  :: assignment(=) => assign_FractionMatrix_IntFractions
    procedure, private ::                  assign_FractionMatrix_IntFractions
    
    procedure, public  :: str => str_FractionMatrix
  end type
  
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
    module procedure frac_IntVector
    module procedure frac_IntMatrix
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
contains

! ----------------------------------------------------------------------
! Procedures involving contents_
! ----------------------------------------------------------------------

! Assignment.
subroutine assign_FractionVector_IntFractions(output,input)
  implicit none
  
  class(FractionVector), intent(inout) :: output
  type(IntFraction),     intent(in)    :: input(:)
  
  output%contents_ = input
end subroutine

subroutine assign_FractionMatrix_IntFractions(output,input)
  implicit none
  
  class(FractionMatrix), intent(inout) :: output
  type(IntFraction),     intent(in)    :: input(:,:)
  
  output%contents_ = input
end subroutine

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
function vec_IntFractions(input) result(output)
  implicit none
  
  type(IntFraction), intent(in) :: input(:)
  type(FractionVector)          :: output
  
  output = input
end function

function mat_IntFractions(input) result(output)
  implicit none
  
  type(IntFraction), intent(in) :: input(:,:)
  type(FractionMatrix)          :: output
  
  output = input
end function

function mat_IntFractions_shape(input,m,n) result(output)
  implicit none
  
  type(IntFraction), intent(in) :: input(:)
  integer,           intent(in) :: m
  integer,           intent(in) :: n
  type(FractionMatrix)          :: output
  
  output = transpose(reshape(input, [m,n]))
end function

function frac_IntVector(input) result(output)
  implicit none
  
  type(IntVector), intent(in)    :: input
  type(IntFraction), allocatable :: output(:)
  
  output = frac(int(input))
end function

function frac_IntMatrix(input) result(output)
  implicit none
  
  type(IntMatrix), intent(in)    :: input
  type(IntFraction), allocatable :: output(:,:)
  
  output = frac(int(input))
end function

function dblevec_FractionVector(input) result(output)
  implicit none
  
  type(FractionVector), intent(in) :: input
  type(RealVector)                 :: output
  
  output = dble(frac(input))
end function

function dblemat_FractionMatrix(input) result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: input
  type(RealMatrix)                 :: output
  
  output = dble(frac(input))
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
  
  output = transpose(frac(this))
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
 
  output(1,1) = input(2,2)*input(3,3)-input(3,2)*input(2,3)
  output(1,2) = input(3,2)*input(1,3)-input(1,2)*input(3,3)
  output(1,3) = input(1,2)*input(2,3)-input(2,2)*input(1,3)
  output(2,1) = input(2,3)*input(3,1)-input(3,3)*input(2,1)
  output(2,2) = input(3,3)*input(1,1)-input(1,3)*input(3,1)
  output(2,3) = input(1,3)*input(2,1)-input(2,3)*input(1,1)
  output(3,1) = input(2,1)*input(3,2)-input(3,1)*input(2,2)
  output(3,2) = input(3,1)*input(1,2)-input(1,1)*input(3,2)
  output(3,3) = input(1,1)*input(2,2)-input(2,1)*input(1,2)
  
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
  
  output = invert(int(input))
end function

! Linear algebra.
impure elemental function add_FractionVector_FractionVector(this,that) result(output)
  implicit none
  
  type(FractionVector), intent(in) :: this
  type(FractionVector), intent(in) :: that
  type(FractionVector)             :: output
  
  output = frac(this) + frac(that)
end function

impure elemental function add_FractionVector_IntVector(this,that) result(output)
  implicit none
  
  type(IntVector),      intent(in) :: this
  type(FractionVector), intent(in) :: that
  type(FractionVector)             :: output
  
  output = int(this) + frac(that)
end function

impure elemental function add_IntVector_FractionVector(this,that) result(output)
  implicit none
  
  type(FractionVector), intent(in) :: this
  type(IntVector),      intent(in) :: that
  type(FractionVector)             :: output
  
  output = frac(this) + int(that)
end function

impure elemental function add_FractionMatrix_FractionMatrix(this,that) result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  type(FractionMatrix), intent(in) :: that
  type(FractionMatrix)             :: output
  
  output = frac(this) + frac(that)
end function

impure elemental function add_FractionMatrix_IntMatrix(this,that) result(output)
  implicit none
  
  type(IntMatrix),      intent(in) :: this
  type(FractionMatrix), intent(in) :: that
  type(FractionMatrix)             :: output
  
  output = int(this) + frac(that)
end function

impure elemental function add_IntMatrix_FractionMatrix(this,that) result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  type(IntMatrix),      intent(in) :: that
  type(FractionMatrix)             :: output
  
  output = frac(this) + int(that)
end function

impure elemental function negative_FractionVector(this) result(output)
  implicit none
  
  type(FractionVector), intent(in) :: this
  type(FractionVector)             :: output
  
  output = -frac(this)
end function

impure elemental function negative_FractionMatrix(this) result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  type(FractionMatrix)             :: output
  
  output = -frac(this)
end function

impure elemental function subtract_FractionVector_FractionVector(this,that) result(output)
  implicit none
  
  type(FractionVector), intent(in) :: this
  type(FractionVector), intent(in) :: that
  type(FractionVector)             :: output
  
  output = frac(this) - frac(that)
end function

impure elemental function subtract_FractionVector_IntVector(this,that) result(output)
  implicit none
  
  type(FractionVector), intent(in) :: this
  type(IntVector),      intent(in) :: that
  type(FractionVector)             :: output
  
  output = frac(this) - int(that)
end function

impure elemental function subtract_IntVector_FractionVector(this,that) result(output)
  implicit none
  
  type(IntVector),      intent(in) :: this
  type(FractionVector), intent(in) :: that
  type(FractionVector)             :: output
  
  output = int(this) - frac(that)
end function

impure elemental function subtract_FractionMatrix_FractionMatrix(this,that) result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  type(FractionMatrix), intent(in) :: that
  type(FractionMatrix)             :: output
  
  output = frac(this) - frac(that)
end function

impure elemental function subtract_FractionMatrix_IntMatrix(this,that) result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  type(IntMatrix),      intent(in) :: that
  type(FractionMatrix)             :: output
  
  output = frac(this) - int(that)
end function

impure elemental function subtract_IntMatrix_FractionMatrix(this,that) result(output)
  implicit none
  
  type(IntMatrix),      intent(in) :: this
  type(FractionMatrix), intent(in) :: that
  type(FractionMatrix)             :: output
  
  output = int(this) - frac(that)
end function

impure elemental function multiply_FractionVector_integer(this,that) result(output)
  implicit none
  
  type(FractionVector), intent(in) :: this
  integer,              intent(in) :: that
  type(FractionVector)             :: output
  
  output = frac(this) * that
end function

impure elemental function multiply_integer_FractionVector(this,that) result(output)
  implicit none
  
  integer,              intent(in) :: this
  type(FractionVector), intent(in) :: that
  type(FractionVector)             :: output
  
  output = this * frac(that)
end function

impure elemental function multiply_FractionVector_IntFraction(this,that) result(output)
  implicit none
  
  type(FractionVector), intent(in) :: this
  type(IntFraction),    intent(in) :: that
  type(FractionVector)             :: output
  
  output = frac(this) * that
end function

impure elemental function multiply_IntFraction_FractionVector(this,that) result(output)
  implicit none
  
  type(IntFraction),    intent(in) :: this
  type(FractionVector), intent(in) :: that
  type(FractionVector)             :: output
  
  output = this * frac(that)
end function

impure elemental function multiply_FractionMatrix_integer(this,that) result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  integer,              intent(in) :: that
  type(FractionMatrix)             :: output
  
  output = frac(this) * that
end function

impure elemental function multiply_integer_FractionMatrix(this,that) result(output)
  implicit none
  
  integer,              intent(in) :: this
  type(FractionMatrix), intent(in) :: that
  type(FractionMatrix)             :: output
  
  output = this * frac(that)
end function

impure elemental function multiply_FractionMatrix_IntFraction(this,that) result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  type(IntFraction),    intent(in) :: that
  type(FractionMatrix)             :: output
  
  output = frac(this) * that
end function

impure elemental function multiply_IntFraction_FractionMatrix(this,that) result(output)
  implicit none
  
  type(IntFraction),    intent(in) :: this
  type(FractionMatrix), intent(in) :: that
  type(FractionMatrix)             :: output
  
  output = this * frac(that)
end function

impure elemental function dot_FractionVector_FractionVector(this,that) result(output)
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
  output = 0
  do i=1,size(this)
    output = output + a(i)*b(i)
  enddo
end function

impure elemental function dot_FractionVector_IntVector(this,that) result(output)
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
  output = 0
  do i=1,size(this)
    output = output + a(i)*b(i)
  enddo
end function

impure elemental function dot_IntVector_FractionVector(this,that) result(output)
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
  output = 0
  do i=1,size(this)
    output = output + a(i)*b(i)
  enddo
end function

impure elemental function dot_FractionMatrix_FractionVector(this,that) result(output)
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
  output = contents
end function

impure elemental function dot_FractionMatrix_IntVector(this,that) result(output)
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
  output = contents
end function

impure elemental function dot_IntMatrix_FractionVector(this,that) result(output)
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
  output = contents
end function

impure elemental function dot_FractionVector_FractionMatrix(this,that) result(output)
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
  output = contents
end function

impure elemental function dot_FractionVector_IntMatrix(this,that) result(output)
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
  output = contents
end function

impure elemental function dot_IntVector_FractionMatrix(this,that) result(output)
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
  output = contents
end function

impure elemental function dot_FractionMatrix_FractionMatrix(this,that) result(output)
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
  output = contents
end function

impure elemental function dot_FractionMatrix_IntMatrix(this,that) result(output)
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
  output = contents
end function

impure elemental function dot_IntMatrix_FractionMatrix(this,that) result(output)
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
  output = contents
end function

impure elemental function divide_FractionVector_integer(this,that) result(output)
  implicit none
  
  type(FractionVector), intent(in) :: this
  integer,              intent(in) :: that
  type(FractionVector)             :: output
  
  output = frac(this) / that
end function

impure elemental function divide_FractionVector_IntFraction(this,that) result(output)
  implicit none
  
  type(FractionVector), intent(in) :: this
  type(IntFraction),    intent(in) :: that
  type(FractionVector)             :: output
  
  output = frac(this) / that
end function

impure elemental function divide_FractionMatrix_integer(this,that) result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  integer,              intent(in) :: that
  type(FractionMatrix)             :: output
  
  output = frac(this) / that
end function

impure elemental function divide_FractionMatrix_IntFraction(this,that) result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  type(IntFraction),    intent(in) :: that
  type(FractionMatrix)             :: output
  
  output = frac(this) / that
end function

! I/O overloads.
function str_FractionVector(this) result(output)
  implicit none
  
  class(FractionVector), intent(in) :: this
  type(String)                      :: output
  
  type(IntFraction), allocatable :: contents(:)
  
  integer :: i
  
  contents = frac(this)
  output = ''
  do i=1,size(this)
    output = output//' '//contents(i)
  enddo
end function

function str_FractionMatrix(this) result(output)
  implicit none
  
  class(FractionMatrix), intent(in) :: this
  type(String), allocatable         :: output(:)
  
  type(IntFraction), allocatable :: contents(:,:)
  
  integer :: i,ialloc
  
  allocate(output(size(this,1)), stat=ialloc); call err(ialloc)
  
  do i=1,size(this,1)
    output(i) = join(str(contents(i,:)))
  enddo
end function
end module
