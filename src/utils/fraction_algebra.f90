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
! Assignment.
! ----------------------------------------------------------------------
pure subroutine assign_FractionVector_IntFractions(output,input)
  implicit none
  
  class(FractionVector), intent(inout) :: output
  type(IntFraction),     intent(in)    :: input(:)
  
  output%contents_ = input
end subroutine

pure subroutine assign_FractionMatrix_IntFractions(output,input)
  implicit none
  
  class(FractionMatrix), intent(inout) :: output
  type(IntFraction),     intent(in)    :: input(:,:)
  
  output%contents_ = input
end subroutine

! ----------------------------------------------------------------------
! Conversion to and from vector and matrix types.
! ----------------------------------------------------------------------
pure function vec_IntFractions(input) result(output)
  implicit none
  
  type(IntFraction), intent(in) :: input(:)
  type(FractionVector)          :: output
  
  output = input
end function

pure function mat_IntFractions(input) result(output)
  implicit none
  
  type(IntFraction), intent(in) :: input(:,:)
  type(FractionMatrix)          :: output
  
  output = input
end function

pure function mat_IntFractions_shape(input,m,n) result(output)
  implicit none
  
  type(IntFraction), intent(in) :: input(:)
  integer,           intent(in) :: m
  integer,           intent(in) :: n
  type(FractionMatrix)          :: output
  
  output = transpose(reshape(input, [m,n]))
end function

pure function frac_FractionVector(input) result(output)
  implicit none
  
  type(FractionVector), intent(in) :: input
  type(IntFraction), allocatable   :: output(:)
  
  output = input%contents_
end function

pure function frac_FractionMatrix(input) result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: input
  type(IntFraction), allocatable   :: output(:,:)
  
  output = input%contents_
end function

pure function frac_IntVector(input) result(output)
  implicit none
  
  type(IntVector), intent(in)    :: input
  type(IntFraction), allocatable :: output(:)
  
  output = frac(int(input))
end function

pure function frac_IntMatrix(input) result(output)
  implicit none
  
  type(IntMatrix), intent(in)    :: input
  type(IntFraction), allocatable :: output(:,:)
  
  output = frac(int(input))
end function

! ----------------------------------------------------------------------
! Properties of the vectors and matrices.
! ----------------------------------------------------------------------
pure function size_FractionVector(this) result(output)
  implicit none
  
  type(FractionVector), intent(in) :: this
  integer                          :: output
  
  output = size(this%contents_)
end function

pure function size_FractionMatrix(this,dim) result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  integer,              intent(in) :: dim
  integer                          :: output
  
  output = size(this%contents_, dim)
end function

elemental function is_int_FractionVector(this) result(output)
  implicit none
  
  type(FractionVector), intent(in) :: this
  logical                          :: output
  
  output = all(is_int(this%contents_))
end function

elemental function is_int_FractionMatrix(this) result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  logical                          :: output
  
  output = all(is_int(this%contents_))
end function

! ----------------------------------------------------------------------
! Comparisons.
! ----------------------------------------------------------------------
elemental function equality_FractionVector_FractionVector(this,that) &
   & result(output)
  implicit none
  
  type(FractionVector), intent(in) :: this
  type(FractionVector), intent(in) :: that
  logical                          :: output
  
  output = all(this%contents_==that%contents_)
end function

elemental function equality_FractionVector_IntVector(this,that) result(output)
  implicit none
  
  type(FractionVector), intent(in) :: this
  type(IntVector),      intent(in) :: that
  logical                          :: output
  
  output = all(this%contents_==int(that))
end function

elemental function equality_IntVector_FractionVector(this,that) result(output)
  implicit none
  
  type(IntVector),      intent(in) :: this
  type(FractionVector), intent(in) :: that
  logical                          :: output
  
  output = all(int(this)==that%contents_)
end function

elemental function equality_FractionMatrix_FractionMatrix(this,that) &
   & result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  type(FractionMatrix), intent(in) :: that
  logical                          :: output
  
  output = all(this%contents_==that%contents_)
end function

elemental function equality_FractionMatrix_IntMatrix(this,that) result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  type(IntMatrix),      intent(in) :: that
  logical                          :: output
  
  output = all(this%contents_==int(that))
end function

elemental function equality_IntMatrix_FractionMatrix(this,that) result(output)
  implicit none
  
  type(IntMatrix),      intent(in) :: this
  type(FractionMatrix), intent(in) :: that
  logical                          :: output
  
  output = all(int(this)==that%contents_)
end function

elemental function non_equality_FractionVector_FractionVector(this,that) &
   & result(output)
  implicit none
  
  type(FractionVector), intent(in) :: this
  type(FractionVector), intent(in) :: that
  logical                          :: output
  
  output = .not. this==that
end function

elemental function non_equality_FractionVector_IntVector(this,that) &
   & result(output)
  implicit none
  
  type(FractionVector), intent(in) :: this
  type(IntVector),      intent(in) :: that
  logical                          :: output
  
  output = .not. this==that
end function

elemental function non_equality_IntVector_FractionVector(this,that) &
   & result(output)
  implicit none
  
  type(IntVector),      intent(in) :: this
  type(FractionVector), intent(in) :: that
  logical                          :: output
  
  output = .not. this==that
end function

elemental function non_equality_FractionMatrix_FractionMatrix(this,that) &
   & result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  type(FractionMatrix), intent(in) :: that
  logical                          :: output
  
  output = .not. this==that
end function

elemental function non_equality_FractionMatrix_IntMatrix(this,that) &
   & result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  type(IntMatrix),      intent(in) :: that
  logical                          :: output
  
  output = .not. this==that
end function

elemental function non_equality_IntMatrix_FractionMatrix(this,that) &
   & result(output)
  implicit none
  
  type(IntMatrix),      intent(in) :: this
  type(FractionMatrix), intent(in) :: that
  logical                          :: output
  
  output = .not. this==that
end function

! ----------------------------------------------------------------------
! Linear algebra.
! ----------------------------------------------------------------------
pure function add_FractionVector_FractionVector(this,that) result(output)
  implicit none
  
  type(FractionVector), intent(in) :: this
  type(FractionVector), intent(in) :: that
  type(FractionVector)             :: output
  
  output = this%contents_ + that%contents_
end function

pure function add_FractionVector_IntVector(this,that) result(output)
  implicit none
  
  type(IntVector),      intent(in) :: this
  type(FractionVector), intent(in) :: that
  type(FractionVector)             :: output
  
  output = int(this) + that%contents_
end function

pure function add_IntVector_FractionVector(this,that) result(output)
  implicit none
  
  type(FractionVector), intent(in) :: this
  type(IntVector),      intent(in) :: that
  type(FractionVector)             :: output
  
  output = this%contents_ + int(that)
end function

pure function add_FractionMatrix_FractionMatrix(this,that) result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  type(FractionMatrix), intent(in) :: that
  type(FractionMatrix)             :: output
  
  output = this%contents_ + that%contents_
end function

pure function add_FractionMatrix_IntMatrix(this,that) result(output)
  implicit none
  
  type(IntMatrix),      intent(in) :: this
  type(FractionMatrix), intent(in) :: that
  type(FractionMatrix)             :: output
  
  output = int(this) + that%contents_
end function

pure function add_IntMatrix_FractionMatrix(this,that) result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  type(IntMatrix),      intent(in) :: that
  type(FractionMatrix)             :: output
  
  output = this%contents_ + int(that)
end function

pure function negative_FractionVector(this) result(output)
  implicit none
  
  type(FractionVector), intent(in) :: this
  type(FractionVector)             :: output
  
  output = -this%contents_
end function

pure function negative_FractionMatrix(this) result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  type(FractionMatrix)             :: output
  
  output = -this%contents_
end function

pure function subtract_FractionVector_FractionVector(this,that) result(output)
  implicit none
  
  type(FractionVector), intent(in) :: this
  type(FractionVector), intent(in) :: that
  type(FractionVector)             :: output
  
  output = this%contents_ - that%contents_
end function

pure function subtract_FractionVector_IntVector(this,that) result(output)
  implicit none
  
  type(FractionVector), intent(in) :: this
  type(IntVector),      intent(in) :: that
  type(FractionVector)             :: output
  
  output = this%contents_ - int(that)
end function

pure function subtract_IntVector_FractionVector(this,that) result(output)
  implicit none
  
  type(IntVector),      intent(in) :: this
  type(FractionVector), intent(in) :: that
  type(FractionVector)             :: output
  
  output = int(this) - that%contents_
end function

pure function subtract_FractionMatrix_FractionMatrix(this,that) result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  type(FractionMatrix), intent(in) :: that
  type(FractionMatrix)             :: output
  
  output = this%contents_ - that%contents_
end function

pure function subtract_FractionMatrix_IntMatrix(this,that) result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  type(IntMatrix),      intent(in) :: that
  type(FractionMatrix)             :: output
  
  output = this%contents_ - int(that)
end function

pure function subtract_IntMatrix_FractionMatrix(this,that) result(output)
  implicit none
  
  type(IntMatrix),      intent(in) :: this
  type(FractionMatrix), intent(in) :: that
  type(FractionMatrix)             :: output
  
  output = int(this) - that%contents_
end function

pure function multiply_FractionVector_integer(this,that) result(output)
  implicit none
  
  type(FractionVector), intent(in) :: this
  integer,              intent(in) :: that
  type(FractionVector)             :: output
  
  output = this%contents_ * that
end function

pure function multiply_integer_FractionVector(this,that) result(output)
  implicit none
  
  integer,              intent(in) :: this
  type(FractionVector), intent(in) :: that
  type(FractionVector)             :: output
  
  output = this * that%contents_
end function

pure function multiply_FractionVector_IntFraction(this,that) result(output)
  implicit none
  
  type(FractionVector), intent(in) :: this
  type(IntFraction),    intent(in) :: that
  type(FractionVector)             :: output
  
  output = this%contents_ * that
end function

pure function multiply_IntFraction_FractionVector(this,that) result(output)
  implicit none
  
  type(IntFraction),    intent(in) :: this
  type(FractionVector), intent(in) :: that
  type(FractionVector)             :: output
  
  output = this * that%contents_
end function

pure function multiply_FractionMatrix_integer(this,that) result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  integer,              intent(in) :: that
  type(FractionMatrix)             :: output
  
  output = this%contents_ * that
end function

pure function multiply_integer_FractionMatrix(this,that) result(output)
  implicit none
  
  integer,              intent(in) :: this
  type(FractionMatrix), intent(in) :: that
  type(FractionMatrix)             :: output
  
  output = this * that%contents_
end function

pure function multiply_FractionMatrix_IntFraction(this,that) result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  type(IntFraction),    intent(in) :: that
  type(FractionMatrix)             :: output
  
  output = this%contents_ * that
end function

pure function multiply_IntFraction_FractionMatrix(this,that) result(output)
  implicit none
  
  type(IntFraction),    intent(in) :: this
  type(FractionMatrix), intent(in) :: that
  type(FractionMatrix)             :: output
  
  output = this * that%contents_
end function

pure function dot_FractionVector_FractionVector(this,that) result(output)
  implicit none
  
  type(FractionVector), intent(in) :: this
  type(FractionVector), intent(in) :: that
  type(IntFraction)                :: output
  
  integer :: i
  
  output = 0
  do i=1,size(this)
    output = output + this%contents_(i)*that%contents_(i)
  enddo
end function

pure function dot_FractionVector_IntVector(this,that) result(output)
  implicit none
  
  type(FractionVector), intent(in) :: this
  type(IntVector),      intent(in) :: that
  type(IntFraction)                :: output
  
  integer, allocatable :: ints(:)
  
  integer :: i
  
  ints = int(that)
  
  output = 0
  do i=1,size(this)
    output = output + this%contents_(i)*ints(i)
  enddo
end function

pure function dot_IntVector_FractionVector(this,that) result(output)
  implicit none
  
  type(IntVector),      intent(in) :: this
  type(FractionVector), intent(in) :: that
  type(IntFraction)                :: output
  
  integer, allocatable :: ints(:)
  
  integer :: i
  
  ints = int(this)
  
  output = 0
  do i=1,size(this)
    output = output + ints(i)*that%contents_(i)
  enddo
end function

pure function dot_FractionMatrix_FractionVector(this,that) result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  type(FractionVector), intent(in) :: that
  type(FractionVector)             :: output
  
  integer :: i
  
  output = that
  output%contents_ = frac(0)
  do i=1,size(that)
    output%contents_ = output%contents_ + this%contents_(:,i)*that%contents_(i)
  enddo
end function

pure function dot_FractionMatrix_IntVector(this,that) result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  type(IntVector),      intent(in) :: that
  type(FractionVector)             :: output
  
  integer, allocatable :: ints(:)
  
  integer :: i
  
  ints = int(that)
  
  output = frac(that)
  output%contents_ = frac(0)
  do i=1,size(that)
    output%contents_ = output%contents_ + this%contents_(:,i)*ints(i)
  enddo
end function

pure function dot_IntMatrix_FractionVector(this,that) result(output)
  implicit none
  
  type(IntMatrix),      intent(in) :: this
  type(FractionVector), intent(in) :: that
  type(FractionVector)             :: output
  
  integer, allocatable :: ints(:,:)
  
  integer :: i
  
  ints = int(this)
  
  output = that
  output%contents_ = frac(0)
  do i=1,size(that)
    output%contents_ = output%contents_ + ints(:,i)*that%contents_(i)
  enddo
end function

pure function dot_FractionVector_FractionMatrix(this,that) result(output)
  implicit none
  
  type(FractionVector), intent(in) :: this
  type(FractionMatrix), intent(in) :: that
  type(FractionVector)             :: output
  
  integer :: i
  
  output = this
  output%contents_ = frac(0)
  do i=1,size(this)
    output%contents_ = output%contents_ + this%contents_(i)*that%contents_(i,:)
  enddo
end function

pure function dot_FractionVector_IntMatrix(this,that) result(output)
  implicit none
  
  type(FractionVector), intent(in) :: this
  type(IntMatrix),      intent(in) :: that
  type(FractionVector)             :: output
  
  integer, allocatable :: ints(:,:)
  
  integer :: i
  
  ints = int(that)
  
  output = this
  output%contents_ = frac(0)
  do i=1,size(this)
    output%contents_ = output%contents_ + this%contents_(i)*ints(i,:)
  enddo
end function

pure function dot_IntVector_FractionMatrix(this,that) result(output)
  implicit none
  
  type(IntVector),      intent(in) :: this
  type(FractionMatrix), intent(in) :: that
  type(FractionVector)             :: output
  
  integer, allocatable :: ints(:)
  
  integer :: i
  
  ints = int(this)
  
  output = frac(this)
  output%contents_ = frac(0)
  do i=1,size(this)
    output%contents_ = output%contents_ + ints(i)*that%contents_(i,:)
  enddo
end function

pure function dot_FractionMatrix_FractionMatrix(this,that) result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  type(FractionMatrix), intent(in) :: that
  type(FractionMatrix)             :: output
  
  integer :: i,j
  
  output = this
  output%contents_ = frac(0)
  do i=1,size(this,1)
    do j=1,size(this,2)
      output%contents_(i,:) = output%contents_(i,:) &
                          & + this%contents_(i,j)*that%contents_(j,:)
    enddo
  enddo
end function

pure function dot_FractionMatrix_IntMatrix(this,that) result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  type(IntMatrix),      intent(in) :: that
  type(FractionMatrix)             :: output
  
  integer, allocatable :: ints(:,:)
  
  integer :: i,j
  
  ints = int(that)
  
  output = this
  output%contents_ = frac(0)
  do i=1,size(this,1)
    do j=1,size(this,2)
      output%contents_(i,:) = output%contents_(i,:) &
                          & + this%contents_(i,j)*ints(j,:)
    enddo
  enddo
end function

pure function dot_IntMatrix_FractionMatrix(this,that) result(output)
  implicit none
  
  type(IntMatrix),      intent(in) :: this
  type(FractionMatrix), intent(in) :: that
  type(FractionMatrix)             :: output
  
  integer, allocatable :: ints(:,:)
  
  integer :: i,j
  
  ints = int(this)
  
  output = frac(that)
  output%contents_ = frac(0)
  do i=1,size(this,1)
    do j=1,size(this,2)
      output%contents_(i,:) = output%contents_(i,:) &
                          & + ints(i,j)*that%contents_(j,:)
    enddo
  enddo
end function

pure function divide_FractionVector_integer(this,that) result(output)
  implicit none
  
  type(FractionVector), intent(in) :: this
  integer,              intent(in) :: that
  type(FractionVector)             :: output
  
  output = this%contents_ / that
end function

pure function divide_FractionVector_IntFraction(this,that) result(output)
  implicit none
  
  type(FractionVector), intent(in) :: this
  type(IntFraction),    intent(in) :: that
  type(FractionVector)             :: output
  
  output = this%contents_ / that
end function

pure function divide_FractionMatrix_integer(this,that) result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  integer,              intent(in) :: that
  type(FractionMatrix)             :: output
  
  output = this%contents_ / that
end function

pure function divide_FractionMatrix_IntFraction(this,that) result(output)
  implicit none
  
  type(FractionMatrix), intent(in) :: this
  type(IntFraction),    intent(in) :: that
  type(FractionMatrix)             :: output
  
  output = this%contents_ / that
end function

! ----------------------------------------------------------------------
! I/O overloads.
! ----------------------------------------------------------------------
pure function str_FractionVector(this) result(output)
  implicit none
  
  class(FractionVector), intent(in) :: this
  type(String)                      :: output
  
  output = join(str(this%contents_))
end function

function str_FractionMatrix(this) result(output)
  implicit none
  
  class(FractionMatrix), intent(in) :: this
  type(String), allocatable         :: output(:)
  
  integer :: i,ialloc
  
  allocate(output(size(this,1)), stat=ialloc); call err(ialloc)
  
  do i=1,size(this,1)
    output(i) = join(str(this%contents_(i,:)))
  enddo
end function
end module
