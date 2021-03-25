! ======================================================================
! Assorted linear algebra / vector algebra subroutines.
! Includes interfaces for BLAS and LAPACK routines.
! ======================================================================
! N.B. to avoid repeating code with multiple types,
!    this file is mostly built by the preprocessor.
! See unary_algebra.fpp and binary_algebra.fpp for
!    interfaces and function definitions.
module caesar_linear_algebra_module
  use caesar_foundations_module
  use caesar_io_module
  
  use caesar_lapack_wrapper_module
  implicit none
  
  private
  
  public :: IntVector
  public :: RealVector
  public :: ComplexVector
  public :: IntMatrix
  public :: RealMatrix
  public :: ComplexMatrix
  public :: zeroes
  public :: make_identity_matrix
  public :: nint
  public :: dblevec
  public :: dblemat
  public :: cmplxvec
  public :: cmplxmat
  public :: l2_norm
  public :: abs
  public :: real
  public :: aimag
  public :: conjg
  public :: matrices_commute
  public :: matrices_anticommute
  public :: operator(==)
  public :: operator(/=)
  public :: transpose
  public :: hermitian
  public :: determinant
  public :: invert
  public :: linear_least_squares
  
  ! Functionality defined in unary_algebra.fpp.
  public :: int
  public :: dble
  public :: cmplx
  public :: vec
  public :: mat
  public :: row_matrix
  public :: column_matrix
  public :: size
  public :: trace
  public :: sum
  
  ! Functionality defined in binary_algebra.fpp.
  public :: operator(+)
  public :: operator(*)
  public :: operator(/)
  public :: cross_product
  public :: outer_product
  public :: commutator
  public :: anticommutator
  
  ! Functionality defined in both unary_algebra.fpp and binary_algebra.fpp.
  public :: operator(-)
  
  ! --------------------------------------------------
  ! Vector and Matrix classes
  ! --------------------------------------------------
  
  type, extends(Stringable) :: IntVector
    integer, allocatable, private :: contents_(:)
  contains
    procedure, public :: check => check_IntVector
    
    procedure, public :: element => element_IntVector
    
    procedure, public :: read  => read_IntVector
    procedure, public :: write => write_IntVector
  end type
  
  type, extends(Stringable) :: RealVector
    real(dp), allocatable, private :: contents_(:)
  contains
    procedure, public :: check => check_RealVector
    
    procedure, public :: element => element_RealVector
    
    procedure, public :: read  => read_RealVector
    procedure, public :: write => write_RealVector
  end type
  
  type, extends(Stringable) :: ComplexVector
    complex(dp), allocatable, private :: contents_(:)
  contains
    procedure, public :: check => check_ComplexVector
    
    procedure, public :: element => element_ComplexVector
    
    procedure, public :: read  => read_ComplexVector
    procedure, public :: write => write_ComplexVector
  end type
  
  type, extends(Stringsable) :: IntMatrix
    integer, allocatable, private :: contents_(:,:)
  contains
    procedure, public :: check => check_IntMatrix
    
    procedure, public :: element => element_IntMatrix
    procedure, public :: row     => row_IntMatrix
    procedure, public :: column  => column_IntMatrix
    
    procedure, public :: read  => read_IntMatrix
    procedure, public :: write => write_IntMatrix
  end type
  
  type, extends(Stringsable) :: RealMatrix
    real(dp), allocatable, private :: contents_(:,:)
  contains
    procedure, public :: check => check_RealMatrix
    
    procedure, public :: element => element_RealMatrix
    procedure, public :: row     => row_RealMatrix
    procedure, public :: column  => column_RealMatrix
    
    procedure, public :: read  => read_RealMatrix
    procedure, public :: write => write_RealMatrix
  end type
  
  type, extends(Stringsable) :: ComplexMatrix
    complex(dp), allocatable, private :: contents_(:,:)
  contains
    procedure, public :: check => check_ComplexMatrix
    
    procedure, public :: element => element_ComplexMatrix
    procedure, public :: row     => row_ComplexMatrix
    procedure, public :: column  => column_ComplexMatrix
    
    procedure, public :: read  => read_ComplexMatrix
    procedure, public :: write => write_ComplexMatrix
  end type
  
  ! Conversion between types.
  interface dblevec
    module procedure dblevec_IntVector
  end interface
  
  interface dblemat
    module procedure dblemat_IntMatrix
  end interface
  
  interface cmplxvec
    module procedure cmplxvec_IntVector
    module procedure cmplxvec_RealVector
    module procedure cmplxvec_RealVectors
  end interface
  
  interface cmplxmat
    module procedure cmplxmat_IntMatrix
    module procedure cmplxmat_RealMatrix
    module procedure cmplxmat_RealMatrices
  end interface
  
  ! Comparison and arithmetic.
  interface operator(==)
    module procedure equality_IntVector_IntVector
    module procedure equality_IntMatrix_IntMatrix
  end interface
  
  interface operator(/=)
    module procedure non_equality_IntVector_IntVector
    module procedure non_equality_IntMatrix_IntMatrix
  end interface
  
  ! Other operations.
  interface zeroes
    module procedure zeroes_IntVector
    module procedure zeroes_IntMatrix
  end interface
  
  interface nint
    module procedure nint_RealVector
    module procedure nint_RealMatrix
  end interface
  
  interface real
    module procedure real_ComplexVector
    module procedure real_ComplexMatrix
  end interface
  
  interface aimag
    module procedure aimag_ComplexVector
    module procedure aimag_ComplexMatrix
  end interface
  
  interface abs
    module procedure abs_ComplexVector
    module procedure abs_ComplexMatrix
  end interface
  
  interface conjg
    module procedure conjg_ComplexVector
    module procedure conjg_ComplexMatrix
  end interface
  
  interface matrices_commute
    module procedure matrices_commute_IntMatrix_IntMatrix
  end interface
  
  interface matrices_anticommute
    module procedure matrices_anticommute_IntMatrix_IntMatrix
  end interface
  
  interface l2_norm
    module procedure l2_norm_RealVector
    module procedure l2_norm_ComplexVector
  end interface
  
  interface transpose
    module procedure transpose_IntMatrix
    module procedure transpose_RealMatrix
    module procedure transpose_ComplexMatrix
  end interface
  
  interface hermitian
    module procedure hermitian_ComplexMatrix
  end interface
  
  interface determinant
    module procedure determinant_integer
    module procedure determinant_real
    module procedure determinant_IntMatrix
    module procedure determinant_RealMatrix
  end interface
  
  interface invert
    module procedure invert_reals
    module procedure invert_RealMatrix
  end interface
  
  interface linear_least_squares
    module procedure linear_least_squares_reals_reals
    module procedure linear_least_squares_reals_RealVector
    module procedure linear_least_squares_RealMatrix_reals
    module procedure linear_least_squares_RealMatrix_RealVector
  end interface

! Include preprocessed procedure headers.
#include "linear_algebra_includes.fpp"

contains

! Include preprocessed procedure bodies.
#define MACRO_BODY
#include "linear_algebra_includes.fpp"

! ----------------------------------------------------------------------
! Conversions between types.
! ----------------------------------------------------------------------
impure elemental function dblevec_IntVector(input) result(output)
  implicit none
  
  class(IntVector), intent(in) :: input
  type(RealVector)             :: output
  
  call input%check()
  
  output = vec(real(input%contents_,dp))
end function

impure elemental function dblemat_IntMatrix(input) result(output)
  implicit none
  
  class(IntMatrix), intent(in) :: input
  type(RealMatrix)             :: output
  
  call input%check()
  
  output = mat(real(input%contents_,dp))
end function

impure elemental function cmplxvec_IntVector(input) result(output)
  implicit none
  
  class(IntVector), intent(in) :: input
  type(ComplexVector)          :: output
  
  call input%check()
  
  output = vec(cmplx(input%contents_,0.0_dp,dp))
end function

impure elemental function cmplxvec_RealVector(input) result(output)
  implicit none
  
  class(RealVector), intent(in) :: input
  type(ComplexVector)           :: output
  
  call input%check()
  
  output = vec(cmplx(input%contents_,0.0_dp,dp))
end function

impure elemental function cmplxvec_RealVectors(real,imag) &
   & result(output)
  implicit none
  
  class(RealVector), intent(in) :: real
  class(RealVector), intent(in) :: imag
  type(ComplexVector)           :: output
  
  if (size(real)/=size(imag)) then
    call print_line(CODE_ERROR//': Constructing complex vector from &
       &incompatible real vectors.')
    call err()
  endif
  
  output = vec(cmplx(real%contents_,imag%contents_,dp))
end function

impure elemental function cmplxmat_IntMatrix(input) result(output)
  implicit none
  
  class(IntMatrix), intent(in) :: input
  type(ComplexMatrix)          :: output
  
  call input%check()
  
  output = mat(cmplx(input%contents_,0.0_dp,dp))
end function

impure elemental function cmplxmat_RealMatrix(input) result(output)
  implicit none
  
  class(RealMatrix), intent(in) :: input
  type(ComplexMatrix)           :: output
  
  call input%check()
  
  output = mat(cmplx(input%contents_,0.0_dp,dp))
end function

impure elemental function cmplxmat_RealMatrices(real,imag) result(output)
  implicit none
  
  class(RealMatrix), intent(in) :: real
  class(RealMatrix), intent(in) :: imag
  type(ComplexMatrix)           :: output
  
  if (size(real,1)/=size(imag,1) .or. size(real,2)/=size(imag,2)) then
    call print_line(CODE_ERROR//': Constructing complex vector from &
       &incompatible real vectors.')
    call err()
  endif
  
  output = mat(cmplx(real%contents_,imag%contents_,dp))
end function

! ----------------------------------------------------------------------
! All other Vector and Matrix operations.
! ----------------------------------------------------------------------

! Makes a length n vector full of zeroes.
function zeroes_IntVector(n) result(output)
  implicit none
  
  integer, intent(in) :: n
  type(IntVector)     :: output
  
  integer :: i
  
  output = vec([(0,i=1,n)])
end function

! Makes a mxn matrix full of zeroes.
function zeroes_IntMatrix(m,n) result(output)
  implicit none
  
  integer, intent(in) :: m
  integer, intent(in) :: n
  type(IntMatrix)     :: output
  
  integer :: ialloc
  
  allocate(output%contents_(m,n), stat=ialloc); call err(ialloc)
  output%contents_ = 0
end function

! Makes an nxn identity matrix.
function make_identity_matrix(n) result(output)
  implicit none
  
  integer, intent(in) :: n
  type(IntMatrix)     :: output
  
  integer :: i,ialloc
  
  allocate(output%contents_(n,n), stat=ialloc); call err(ialloc)
  output%contents_ = 0
  do i=1,n
    output%contents_(i,i) = 1
  enddo
end function

! Conversion to nearest integer.
impure elemental function nint_RealVector(input) result(output)
  implicit none
  
  type(RealVector), intent(in) :: input
  type(IntVector)              :: output
  
  call input%check()
  
  output = vec(nint(input%contents_))
end function

impure elemental function nint_RealMatrix(input) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: input
  type(IntMatrix)              :: output
  
  call input%check()
  
  output = mat(nint(input%contents_))
end function

! Real part of a complex object.
impure elemental function real_ComplexVector(input) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: input
  type(RealVector)                :: output
  
  call input%check()
  
  output = vec(real(input%contents_,dp))
end function

impure elemental function real_ComplexMatrix(input) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: input
  type(RealMatrix)                :: output
  
  call input%check()
  
  output = mat(real(input%contents_,dp))
end function

! Imaginary part of a complex object.
impure elemental function aimag_ComplexVector(input) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: input
  type(RealVector)                :: output
  
  call input%check()
  
  output = vec(aimag(input%contents_))
end function

impure elemental function aimag_ComplexMatrix(input) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: input
  type(RealMatrix)                :: output
  
  call input%check()
  
  output = mat(aimag(input%contents_))
end function

! Element-wise abs of a complex object.
impure elemental function abs_ComplexVector(input) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: input
  type(RealVector)                :: output
  
  call input%check()
  
  output = vec(abs(input%contents_))
end function

impure elemental function abs_ComplexMatrix(input) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: input
  type(RealMatrix)                :: output
  
  call input%check()
  
  output = mat(abs(input%contents_))
end function

! Conjugate of a complex object.
impure elemental function conjg_ComplexVector(input) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: input
  type(ComplexVector)             :: output
  
  call input%check()
  
  output = vec(conjg(input%contents_))
end function

impure elemental function conjg_ComplexMatrix(input) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: input
  type(ComplexMatrix)             :: output
  
  call input%check()
  
  output = mat(conjg(input%contents_))
end function

! Check if two matrices commute or anti-commute.
function matrices_commute_IntMatrix_IntMatrix(this,that) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: this
  type(IntMatrix), intent(in) :: that
  logical                     :: output
  
  output = commutator(this,that)==zeroes(size(this,1),size(this,1))
end function

function matrices_anticommute_IntMatrix_IntMatrix(this,that) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: this
  type(IntMatrix), intent(in) :: that
  logical                     :: output
  
  output = anticommutator(this,that)==zeroes(size(this,1),size(this,1))
end function

! Equality.
impure elemental function equality_IntVector_IntVector(a,b) result(output)
  implicit none
  
  type(IntVector), intent(in) :: a
  type(IntVector), intent(in) :: b
  logical                     :: output
  
  if (size(a)/=size(b)) then
    output = .false.
  else
    output = all(int(a)==int(b))
  endif
end function

impure elemental function equality_IntMatrix_IntMatrix(a,b) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: a
  type(IntMatrix), intent(in) :: b
  logical                     :: output
  
  if (size(a,1)/=size(b,1) .or. size(a,2)/=size(b,2)) then
    output = .false.
  else
    output = all(int(a)==int(b))
  endif
end function

! Non-equality.
impure elemental function non_equality_IntVector_IntVector(a,b) result(output)
  implicit none
  
  type(IntVector), intent(in) :: a
  type(IntVector), intent(in) :: b
  logical                     :: output
  
  output = .not. a==b
end function

impure elemental function non_equality_IntMatrix_IntMatrix(a,b) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: a
  type(IntMatrix), intent(in) :: b
  logical                     :: output
  
  output = .not. a==b
end function

! L2 norm.
impure elemental function l2_norm_RealVector(input) result(output)
  implicit none
  
  type(RealVector), intent(in) :: input
  real(dp)                     :: output
  
  output = sqrt(input*input)
end function

impure elemental function l2_norm_ComplexVector(input) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: input
  real(dp)                        :: output
  
  output = sqrt(real(input*conjg(input)))
end function

! Transpose.
function transpose_IntMatrix(input) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: input
  type(IntMatrix)             :: output
  
  call input%check()
  
  output = mat(transpose(input%contents_))
end function

function transpose_RealMatrix(input) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: input
  type(RealMatrix)             :: output
  
  call input%check()
  
  output = mat(transpose(input%contents_))
end function

function transpose_ComplexMatrix(input) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: input
  type(ComplexMatrix)             :: output
  
  call input%check()
  
  output = mat(transpose(input%contents_))
end function

! Hermitian conjugate.
function hermitian_ComplexMatrix(input) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: input
  type(ComplexMatrix)             :: output
  
  call input%check()
  
  output = mat(transpose(conjg(input%contents_)))
end function

! Determinant. (3x3 only)
! N.B. In cases where the determinant of a real matrix is small,
!    it is also very error-prone.
function determinant_integer(input) result(output)
  implicit none
  
  integer, intent(in) :: input(:,:)
  integer             :: output
  
  ! Check that input is a 3x3 matrix.
  if (size(input,1)/=3 .or. size(input,2)/=3) then
    call print_line(CODE_ERROR//': Trying to take the determinant of a &
       &matrix which is not 3x3.')
    call err()
  endif
  
  output = input(1,1)*(input(2,2)*input(3,3)-input(2,3)*input(3,2)) &
       & + input(1,2)*(input(2,3)*input(3,1)-input(2,1)*input(3,3)) &
       & + input(1,3)*(input(2,1)*input(3,2)-input(2,2)*input(3,1))
end function

function determinant_real(input) result(output)
  implicit none
  
  real(dp), intent(in) :: input(:,:)
  real(dp)             :: output
  
  ! Check that input is a 3x3 matrix.
  if (size(input,1)/=3 .or. size(input,2)/=3) then
    call print_line(CODE_ERROR//': Trying to take the determinant of a &
       &matrix which is not 3x3.')
    call err()
  endif
  
  output = input(1,1)*(input(2,2)*input(3,3)-input(2,3)*input(3,2)) &
       & + input(1,2)*(input(2,3)*input(3,1)-input(2,1)*input(3,3)) &
       & + input(1,3)*(input(2,1)*input(3,2)-input(2,2)*input(3,1))
end function

function determinant_IntMatrix(input) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: input
  integer                     :: output
  
  call input%check()
  
  output = determinant(input%contents_)
end function

function determinant_RealMatrix(input) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: input
  real(dp)                     :: output
  
  call input%check()
  
  output = determinant(input%contents_)
end function

! Calculates the inverse of a matrix.
function invert_reals(input) result(output)
  implicit none
  
  real(dp), intent(in)  :: input(:,:)
  real(dp), allocatable :: output(:,:)
  
  ! LAPACK variables.
  integer               :: n       ! The size of the matrix.
  integer, allocatable  :: ipiv(:) ! Pivot indices.
  integer               :: info    ! 0 on success.
  integer               :: lwork
  real(dp), allocatable :: work(:)
  
  integer, allocatable :: pivots(:)
  
  integer :: i,ialloc
  
  n = size(input,1)
  if (size(input,2)/=n) then
    call print_line(CODE_ERROR//': Trying to invert non-square matrix.')
    call err()
  endif
  
  allocate( ipiv(n),     &
          & output(n,n), &
          & stat=ialloc); call err(ialloc)
  output = input
  
  ! Run LU factorisation.
  call dgetrf(n,n,output,n,ipiv,info)
  if (info/=0) then
    if (info<0) then
      call print_line(ERROR//' in matrix inversion: argument '//-info//&
         &' had an illegal value.')
    else
      pivots = [(i,i=1,size(ipiv))]
      do i=1,size(ipiv)
        pivots([i,ipiv(i)]) = pivots([ipiv(i),i])
      enddo
    endif
    call print_line(ERROR//' in matrix inversion: row '//pivots(info)//' is &
       &linearly dependent on previous rows.')
    call print_line('Matrix:')
    call print_lines(mat(input))
    call err()
  endif
  
  ! Get the size of the required workspace.
  allocate(work(n),stat=ialloc); call err(ialloc)
  call dgetri(n,output,n,ipiv,work,-1,info)
  if (info/=0) then
    call print_line(ERROR//' in matrix inversion: dgetri error code: '//info)
    call err()
  endif
  lwork = nint(work(1))
  deallocate(work, stat=ialloc); call err(ialloc)
  allocate(work(lwork), stat=ialloc); call err(ialloc)
  
  ! Run matrix inversion.
  call dgetri(n,output,n,ipiv,work,lwork,info)
  if (info/=0) then
    call print_line(ERROR//'in matrix inversion: dgetri error code: '//info)
    call err()
  endif
end function

function invert_RealMatrix(input) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: input
  type(RealMatrix)             :: output
  
  call input%check()
  
  output = mat(invert(input%contents_))
end function

! --------------------------------------------------
! Finds x which minimises the least-squares fit l=(a.x-b)**2.
! --------------------------------------------------
! If a has less than full rank, the routine will fail. In this case:
!    - If return_empty_on_error is not present or .false.,
!         an error will be thrown.
!    - If return_empty_on_error is .true., an empty vector will be returned.
function linear_least_squares_reals_reals(a,b,return_empty_on_error) &
   & result(output)
  implicit none
  
  real(dp), intent(in)           :: a(:,:)
  real(dp), intent(in)           :: b(:)
  logical,  intent(in), optional :: return_empty_on_error
  type(RealVector)               :: output
  
  ! LAPACK dgels variables.
  real(dp), allocatable :: a2(:,:)
  real(dp), allocatable :: b2(:,:)
  integer               :: m,n
  real(dp), allocatable :: work(:)
  integer               :: lwork
  integer               :: info
  
  ! Temporary variables.
  integer :: ialloc
  
  ! Check that input dimensions are consistent.
  m = size(a,1)
  n = size(a,2)
  if (size(b)/=m) then
    call err()
  elseif (m<n) then
    call print_line('Error in least-squares optimisation: there are more &
       &variables to be fit than input data-points.')
    call err()
  endif
  
  ! Copy a and b so that they are not changed by dgels.
  a2 = a
  b2 = reshape(b, [size(b),1])
  
  ! Calculate optimal workspace size.
  allocate(work(2*m*n), stat=ialloc); call err(ialloc)
  call dgels('N',m,n,1,a2,m,b2,m,work,-1,info)
  if (info/=0) then
    call print_line(ERROR//' in linear regression: dgels error code: '//info)
    call err()
  endif
  lwork = nint(work(1))
  deallocate(work, stat=ialloc); call err(ialloc)
  allocate(work(lwork), stat=ialloc); call err(ialloc)
  
  ! Run linear least-squares optimisation.
  call dgels('N',m,n,1,a2,m,b2,m,work,lwork,info)
  if (info/=0) then
    if (info>0) then
      if (present(return_empty_on_error)) then
        if (return_empty_on_error) then
          output = vec([real(dp)::])
          return
        endif
      endif
    endif
    
    call print_line(ERROR//' in linear regression: dgels error code: '//info)
    call err()
  endif
  
  output = vec(b2(:n,1))
end function

function linear_least_squares_reals_RealVector(a,b) result(x)
  implicit none
  
  real(dp),         intent(in) :: a(:,:)
  type(RealVector), intent(in) :: b
  type(RealVector)             :: x
  
  x = linear_least_squares(a,dble(b))
end function

function linear_least_squares_RealMatrix_reals(a,b) result(x)
  implicit none
  
  type(RealMatrix), intent(in) :: a
  real(dp),         intent(in) :: b(:)
  type(RealVector)             :: x
  
  x = linear_least_squares(dble(a),b)
end function

function linear_least_squares_RealMatrix_RealVector(a,b) result(x)
  implicit none
  
  type(RealMatrix), intent(in) :: a
  type(RealVector), intent(in) :: b
  type(RealVector)             :: x
  
  x = linear_least_squares(dble(a),dble(b))
end function
end module
