! ======================================================================
! Assorted linear algebra / vector algebra subroutines.
! Includes interfaces for BLAS and LAPACK routines.
! ======================================================================
module linear_algebra_module
  use constants_module, only : dp
  use string_module
  use io_module
  
  use printable_module
  implicit none
  
  ! --------------------------------------------------
  ! Vector and Matrix classes
  ! --------------------------------------------------
  type, extends(Stringable) :: IntVector
    integer, allocatable, private :: contents_(:)
  contains
    generic,   public  :: assignment(=) => assign_IntVector_integers
    procedure, private ::                  assign_IntVector_integers
    
    procedure, public  :: str => str_IntVector
  end type
  
  type, extends(Stringable) :: RealVector
    real(dp), allocatable, private :: contents_(:)
  contains
    generic,   public  :: assignment(=) => assign_RealVector_reals
    procedure, private ::                  assign_RealVector_reals
    
    procedure, public  :: str => str_RealVector
  end type
  
  type, extends(Stringable) :: ComplexVector
    complex(dp), allocatable, private :: contents_(:)
  contains
    generic,   public  :: assignment(=) => assign_ComplexVector_complexes
    procedure, private ::                  assign_ComplexVector_complexes
    
    procedure, public  :: str => str_ComplexVector
  end type
  
  type, extends(Printable) :: IntMatrix
    integer, allocatable, private :: contents_(:,:)
  contains
    generic,   public  :: assignment(=) => assign_IntMatrix_integers
    procedure, private ::                  assign_IntMatrix_integers
    
    procedure, public  :: str => str_IntMatrix
  end type
  
  type, extends(Printable) :: RealMatrix
    real(dp), allocatable, private :: contents_(:,:)
  contains
    generic,   public  :: assignment(=) => assign_RealMatrix_reals
    procedure, private ::                  assign_RealMatrix_reals
    
    procedure, public  :: str => str_RealMatrix
  end type
  
  type, extends(Printable) :: ComplexMatrix
    complex(dp), allocatable, private :: contents_(:,:)
  contains
    generic,   public  :: assignment(=) => assign_ComplexMatrix_complexes
    procedure, private ::                  assign_ComplexMatrix_complexes
    
    procedure, public  :: str => str_ComplexMatrix
  end type
  
  ! --------------------------------------------------
  ! Overloads for vector / matrix operations.
  ! --------------------------------------------------
  interface vec
    module procedure vec_integers
    module procedure vec_reals
    module procedure vec_complexes
  end interface
  
  interface mat
    module procedure mat_integers
    module procedure mat_reals
    module procedure mat_complexes
    module procedure mat_integers_shape
    module procedure mat_reals_shape
    module procedure mat_complexes_shape
  end interface
  
  interface zeroes
    module procedure zeroes_IntVector
    module procedure zeroes_IntMatrix
  end interface
  
  interface int
    module procedure int_IntVector
    module procedure int_IntMatrix
  end interface
  
  interface nint
    module procedure nint_RealVector
    module procedure nint_RealMatrix
  end interface
  
  interface dble
    module procedure dble_IntVector
    module procedure dble_IntMatrix
    module procedure dble_RealVector
    module procedure dble_RealMatrix
  end interface
  
  interface cmplx
    module procedure cmplx_IntVector
    module procedure cmplx_IntMatrix
    module procedure cmplx_RealVector
    module procedure cmplx_RealMatrix
    module procedure cmplx_ComplexVector
    module procedure cmplx_ComplexMatrix
    module procedure cmplx_RealVectors
    module procedure cmplx_RealMatrices
  end interface
  
  interface real
    module procedure real_ComplexVector
    module procedure real_ComplexMatrix
  end interface
  
  interface aimag
    module procedure aimag_ComplexVector
    module procedure aimag_ComplexMatrix
  end interface
  
  interface conjg
    module procedure conjg_ComplexVector
    module procedure conjg_ComplexMatrix
  end interface
  
  interface size
    module procedure size_IntVector
    module procedure size_RealVector
    module procedure size_ComplexVector
    module procedure size_IntMatrix
    module procedure size_RealMatrix
    module procedure size_ComplexMatrix
  end interface
  
  interface trace
    module procedure trace_IntMatrix
    module procedure trace_RealMatrix
    module procedure trace_ComplexMatrix
  end interface
  
  interface operator(==)
    module procedure equality_IntVector_IntVector
    module procedure equality_IntMatrix_IntMatrix
  end interface
  
  interface operator(/=)
    module procedure non_equality_IntVector_IntVector
    module procedure non_equality_IntMatrix_IntMatrix
  end interface
  
  interface operator(+)
    module procedure add_IntVector_IntVector
    module procedure add_IntVector_RealVector
    module procedure add_RealVector_IntVector
    module procedure add_RealVector_RealVector
    module procedure add_RealVector_ComplexVector
    module procedure add_ComplexVector_RealVector
    module procedure add_ComplexVector_ComplexVector
    
    module procedure add_IntMatrix_IntMatrix
    module procedure add_IntMatrix_RealMatrix
    module procedure add_RealMatrix_IntMatrix
    module procedure add_RealMatrix_RealMatrix
    module procedure add_RealMatrix_ComplexMatrix
    module procedure add_ComplexMatrix_RealMatrix
    module procedure add_ComplexMatrix_ComplexMatrix
  end interface
  
  interface operator(-)
    module procedure negative_IntVector
    module procedure negative_RealVector
    module procedure negative_ComplexVector
    module procedure negative_IntMatrix
    module procedure negative_RealMatrix
    module procedure negative_ComplexMatrix
    
    module procedure subtract_IntVector_IntVector
    module procedure subtract_IntVector_RealVector
    module procedure subtract_RealVector_IntVector
    module procedure subtract_RealVector_RealVector
    module procedure subtract_RealVector_ComplexVector
    module procedure subtract_ComplexVector_RealVector
    module procedure subtract_ComplexVector_ComplexVector
    
    module procedure subtract_IntMatrix_IntMatrix
    module procedure subtract_IntMatrix_RealMatrix
    module procedure subtract_RealMatrix_IntMatrix
    module procedure subtract_RealMatrix_RealMatrix
    module procedure subtract_RealMatrix_ComplexMatrix
    module procedure subtract_ComplexMatrix_RealMatrix
    module procedure subtract_ComplexMatrix_ComplexMatrix
  end interface
  
  interface operator(*)
    module procedure multiply_IntVector_integer
    module procedure multiply_integer_IntVector
    module procedure multiply_IntVector_real
    module procedure multiply_real_IntVector
    module procedure multiply_Intvector_complex
    module procedure multiply_complex_Intvector
    
    module procedure multiply_RealVector_integer
    module procedure multiply_integer_RealVector
    module procedure multiply_RealVector_real
    module procedure multiply_real_RealVector
    module procedure multiply_RealVector_complex
    module procedure multiply_complex_RealVector
    
    module procedure multiply_ComplexVector_integer
    module procedure multiply_integer_ComplexVector
    module procedure multiply_ComplexVector_real
    module procedure multiply_real_ComplexVector
    module procedure multiply_ComplexVector_complex
    module procedure multiply_complex_ComplexVector
    
    module procedure multiply_IntMatrix_integer
    module procedure multiply_integer_IntMatrix
    module procedure multiply_IntMatrix_real
    module procedure multiply_real_IntMatrix
    module procedure multiply_IntMatrix_complex
    module procedure multiply_complex_IntMatrix
    
    module procedure multiply_RealMatrix_integer
    module procedure multiply_integer_RealMatrix
    module procedure multiply_RealMatrix_real
    module procedure multiply_real_RealMatrix
    module procedure multiply_RealMatrix_complex
    module procedure multiply_complex_RealMatrix
    
    module procedure multiply_ComplexMatrix_integer
    module procedure multiply_integer_ComplexMatrix
    module procedure multiply_ComplexMatrix_real
    module procedure multiply_real_ComplexMatrix
    module procedure multiply_ComplexMatrix_complex
    module procedure multiply_complex_ComplexMatrix
    
    module procedure dot_IntVector_IntVector
    module procedure dot_IntVector_RealVector
    module procedure dot_RealVector_IntVector
    module procedure dot_RealVector_RealVector
    module procedure dot_RealVector_ComplexVector
    module procedure dot_ComplexVector_RealVector
    module procedure dot_ComplexVector_ComplexVector
    
    module procedure dot_IntVector_IntMatrix
    module procedure dot_IntVector_RealMatrix
    module procedure dot_RealVector_IntMatrix
    module procedure dot_RealVector_RealMatrix
    module procedure dot_RealVector_ComplexMatrix
    module procedure dot_ComplexVector_RealMatrix
    module procedure dot_ComplexVector_ComplexMatrix
    
    module procedure dot_IntMatrix_IntVector
    module procedure dot_IntMatrix_RealVector
    module procedure dot_RealMatrix_IntVector
    module procedure dot_RealMatrix_RealVector
    module procedure dot_RealMatrix_ComplexVector
    module procedure dot_ComplexMatrix_RealVector
    module procedure dot_ComplexMatrix_ComplexVector
    
    module procedure dot_IntMatrix_IntMatrix
    module procedure dot_IntMatrix_RealMatrix
    module procedure dot_RealMatrix_IntMatrix
    module procedure dot_RealMatrix_RealMatrix
    module procedure dot_RealMatrix_ComplexMatrix
    module procedure dot_ComplexMatrix_RealMatrix
    module procedure dot_ComplexMatrix_ComplexMatrix
  end interface
  
  interface operator(/)
    module procedure divide_IntVector_integer
    module procedure divide_IntVector_real
    module procedure divide_IntVector_complex
    module procedure divide_RealVector_integer
    module procedure divide_RealVector_real
    module procedure divide_RealVector_complex
    module procedure divide_ComplexVector_integer
    module procedure divide_ComplexVector_real
    module procedure divide_ComplexVector_complex
    
    module procedure divide_IntMatrix_integer
    module procedure divide_IntMatrix_real
    module procedure divide_IntMatrix_complex
    module procedure divide_RealMatrix_integer
    module procedure divide_RealMatrix_real
    module procedure divide_RealMatrix_complex
    module procedure divide_ComplexMatrix_integer
    module procedure divide_ComplexMatrix_real
    module procedure divide_ComplexMatrix_complex
  end interface
  
  interface l2_norm
    module procedure l2_norm_RealVector
    module procedure l2_norm_ComplexVector
  end interface
  
  interface outer_product
    module procedure outer_product_RealVector_RealVector
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
  
  ! --------------------------------------------------
  ! The eigen(values and vectors) of a matrix.
  ! --------------------------------------------------
  ! degeneracy(i) = 0 if eigenvector(i) is non-degenerate
  !               = j if eigenvector(i) is degenerate,
  ! where j is an arbitrary but unique id.
  type RealEigenstuff
    real(dp), allocatable :: evals(:)
    real(dp), allocatable :: evecs(:,:)
    integer,  allocatable :: degeneracy(:)
  end type
  
  type ComplexEigenstuff
    real(dp),    allocatable :: evals(:)
    complex(dp), allocatable :: evecs(:,:)
    integer,     allocatable :: degeneracy(:)
  end type

  interface size
    module procedure size_RealEigenstuff
    module procedure size_ComplexEigenstuff
  end interface
  
  interface calculate_eigenstuff
    module procedure calculate_RealEigenstuff_reals
    module procedure calculate_RealEigenstuff_RealMatrix
    module procedure calculate_ComplexEigenstuff_complexes
    module procedure calculate_ComplexEigenstuff_ComplexMatrix
  end interface
  
  ! --------------------------------------------------
  ! Linear least-squares optimisation.
  ! --------------------------------------------------
  interface linear_least_squares
    module procedure linear_least_squares_reals_reals
    module procedure linear_least_squares_reals_RealVector
    module procedure linear_least_squares_RealMatrix_reals
    module procedure linear_least_squares_RealMatrix_RealVector
  end interface
  
  ! --------------------------------------------------
  ! BLAS / LAPACK interface.
  ! --------------------------------------------------
  interface
    ! Copies a real vector. Equivalent to DY = DX.
    pure subroutine dcopy(N,DX,INCX,DY,INCY)
      import :: dp
      implicit none
      
      integer,  intent(in)  :: N     ! length of vectors
      real(dp), intent(in)  :: DX(*) ! input vector
      integer,  intent(in)  :: INCX  ! increment along DX
      real(dp), intent(out) :: DY(*) ! output vector
      integer,  intent(in)  :: INCY  ! increment along DY
    end subroutine
    
    ! Copies complex vector. Equivalent to ZY = ZX.
    pure subroutine zcopy(N,ZX,INCX,ZY,INCY)
      import :: dp
      implicit none
      
      integer,     intent(in)  :: N     ! length of vectors
      complex(dp), intent(in)  :: ZX(*) ! input vector
      integer,     intent(in)  :: INCX  ! increment along ZX
      complex(dp), intent(out) :: ZY(*) ! output vector
      integer,     intent(in)  :: INCY  ! increment along ZY
    end subroutine
    
    ! Real dot product. Returns DX.DY.
    pure function ddot(N,DX,INCX,DY,INCY) result(output)
      import :: dp
      implicit none
      
      integer,  intent(in) :: N     ! length of vectors
      real(dp), intent(in) :: DX(*) ! first vector
      integer,  intent(in) :: INCX  ! increment along DX
      real(dp), intent(in) :: DY(*) ! second vector
      integer,  intent(in) :: INCY  ! increment along DY
      real(dp)             :: output
    end function
    
    ! Multiplies real vector by real scalar. Equivalent to DX *= DA.
    pure subroutine dscal(N,DA,DX,INCX)
      import :: dp
      implicit none
      
      integer,  intent(in)    :: N     ! length of vector
      real(dp), intent(in)    :: DA    ! scalar
      real(dp), intent(inout) :: DX(*) ! vector
      integer,  intent(in)    :: INCX  ! increment along DX
    end subroutine
    
    ! Multiplies complex vector by complex scalar. Equivalent to ZX *= ZA.
    pure subroutine zscal(N,ZA,ZX,INCX)
      import :: dp
      implicit none
      
      integer,     intent(in)    :: N     ! length of vector
      complex(dp), intent(in)    :: ZA    ! scalar
      complex(dp), intent(inout) :: ZX(*) ! vector
      integer,     intent(in)    :: INCX  ! increment along ZX
    end subroutine
    
    ! Complex norm. Returns sqrt(X.X).
    pure function dznrm2(N,X,INCX) result(output)
      import :: dp
      implicit none
      
      integer,     intent(in) :: N      ! length of vector
      complex(dp), intent(in) :: X(*)   ! vector
      integer,     intent(in) :: INCX   ! increment along X
      real(dp)                :: output ! result
    end function
    
    ! Finds the eigenvalues of a hermitian matrix.
    pure subroutine zheev(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,RWORK,INFO)
      import :: dp
      implicit none
      
      character(1), intent(in)    :: JOBZ     ! N/V: if V, calculate eigenvecs.
      character(1), intent(in)    :: UPLO     ! U/L: upper/lower triangle.
      integer,      intent(in)    :: N        ! the order of A.
      integer,      intent(in)    :: LDA      ! the dimension of A.
      complex(dp),  intent(inout) :: A(LDA,*) ! Hermitian matrix.
      real(dp),     intent(out)   :: W(*)     ! eigenvalues of A.
      complex(dp),  intent(out)   :: WORK(*)  ! WORK(1) = optimal LWORK.
      integer,      intent(in)    :: LWORK    ! the length of WORK.
      real(dp),     intent(out)   :: RWORK(*) ! working array.
      integer,      intent(out)   :: INFO     ! 0 on success.
    end subroutine
    
    ! Finds the eigenvalues of a symmetric matrix.
    pure subroutine dsyev(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,INFO)
      import :: dp
      implicit none
      
      character(1), intent(in)    :: JOBZ     ! N/V: if V, calculate eigenvecs.
      character(1), intent(in)    :: UPLO     ! U/L: upper/lower triangle.
      integer,      intent(in)    :: N        ! the order of A.
      integer,      intent(in)    :: LDA      ! the dimension of A.
      real(dp),     intent(inout) :: A(LDA,*) ! symmetric matrix.
      real(dp),     intent(out)   :: W(*)     ! eigenvalues of A.
      real(dp),     intent(out)   :: WORK(*)  ! WORK(1) = optimal LWORK.
      integer,      intent(in)    :: LWORK    ! the length of WORK.
      integer,      intent(out)   :: INFO     ! 0 on success.
    end subroutine
    
    ! Minimises the least-squares fit L=(A.X-B)**2.
    pure subroutine dgels(TRANS,M,N,NRHS,A,LDA,B,LDB,WORK,LWORK,INFO)
      import :: dp
      implicit none
      
      character(1), intent(in)    :: TRANS    ! N/T: if T, A is transposed.
      integer,      intent(in)    :: M        ! size(A,1)=size(B,1).
      integer,      intent(in)    :: N        ! size(A,2)=size(X,1).
      integer,      intent(in)    :: NRHS     ! size(B,2)=size(X,2).
      integer,      intent(in)    :: LDA      ! The leading dimension of A.
      real(dp),     intent(inout) :: A(LDA,*) ! The matrix A.
      integer,      intent(in)    :: LDB      ! The leading dimension of B.
      real(dp),     intent(inout) :: B(LDB,*) ! The matrix B.
      real(dp),     intent(out)   :: WORK(*)  ! WORK(1) = optimal LWORK.
      integer,      intent(in)    :: LWORK    ! size(WORK).
      integer,      intent(out)   :: INFO     ! 0 on success.
    end subroutine
  end interface

contains

! ----------------------------------------------------------------------
! Vector and Matrix operations.
! ----------------------------------------------------------------------
! Assignment.
pure subroutine assign_IntVector_integers(output,input)
  implicit none
  
  class(IntVector), intent(inout) :: output
  integer,          intent(in)    :: input(:)
  
  output%contents_ = input
end subroutine

pure subroutine assign_RealVector_reals(output,input)
  implicit none
  
  class(RealVector), intent(inout) :: output
  real(dp),          intent(in)    :: input(:)
  
  output%contents_ = input
end subroutine

pure subroutine assign_ComplexVector_complexes(output,input)
  implicit none
  
  class(ComplexVector), intent(inout) :: output
  complex(dp),          intent(in)    :: input(:)
  
  output%contents_ = input
end subroutine

pure subroutine assign_IntMatrix_integers(output,input)
  implicit none
  
  class(IntMatrix), intent(inout) :: output
  integer,          intent(in)    :: input(:,:)
  
  output%contents_ = input
end subroutine

pure subroutine assign_RealMatrix_reals(output,input)
  implicit none
  
  class(RealMatrix), intent(inout) :: output
  real(dp),          intent(in)    :: input(:,:)
  
  output%contents_ = input
end subroutine

pure subroutine assign_ComplexMatrix_complexes(output,input)
  implicit none
  
  class(ComplexMatrix), intent(inout) :: output
  complex(dp),          intent(in)    :: input(:,:)
  
  output%contents_ = input
end subroutine

! Conversion to Vector and Matrix.
pure function vec_integers(input) result(output)
  implicit none
  
  integer, intent(in) :: input(:)
  type(IntVector)     :: output
  
  output = input
end function

pure function vec_reals(input) result(output)
  implicit none
  
  real(dp), intent(in) :: input(:)
  type(RealVector)     :: output
  
  output = input
end function

pure function vec_complexes(input) result(output)
  implicit none
  
  complex(dp), intent(in) :: input(:)
  type(ComplexVector)     :: output
  
  output = input
end function

pure function mat_integers(input) result(output)
  implicit none
  
  integer, intent(in) :: input(:,:)
  type(IntMatrix)     :: output
  
  output = input
end function

pure function mat_reals(input) result(output)
  implicit none
  
  real(dp), intent(in) :: input(:,:)
  type(RealMatrix)     :: output
  
  output = input
end function

pure function mat_complexes(input) result(output)
  implicit none
  
  complex(dp), intent(in) :: input(:,:)
  type(ComplexMatrix)     :: output
  
  output = input
end function

pure function mat_integers_shape(input,m,n) result(output)
  implicit none
  
  integer, intent(in) :: input(:)
  integer, intent(in) :: m
  integer, intent(in) :: n
  type(IntMatrix)     :: output
  
  output = transpose(reshape(input, [m,n]))
end function

pure function mat_reals_shape(input,m,n) result(output)
  implicit none
  
  real(dp), intent(in) :: input(:)
  integer,  intent(in) :: m
  integer,  intent(in) :: n
  type(RealMatrix)     :: output
  
  output = transpose(reshape(input, [m,n]))
end function

pure function mat_complexes_shape(input,m,n) result(output)
  implicit none
  
  complex(dp), intent(in) :: input(:)
  integer,     intent(in) :: m
  integer,     intent(in) :: n
  type(ComplexMatrix)     :: output
  
  output = transpose(reshape(input, [m,n]))
end function

! Makes a length n vector full of zeroes.
function zeroes_IntVector(n) result(output)
  implicit none
  
  integer, intent(in) :: n
  type(IntVector)     :: output
  
  integer :: ialloc
  
  allocate(output%contents_(n), stat=ialloc); call err(ialloc)
  output%contents_ = 0
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
  
  integer :: i
  
  output = zeroes(n,n)
  do i=1,n
    output%contents_(i,i) = 1
  enddo
end function

! Conversion to fundamental types.
pure function int_IntVector(input) result(output)
  implicit none
  
  type(IntVector), intent(in) :: input
  integer, allocatable        :: output(:)
  
  output = input%contents_
end function

pure function int_IntMatrix(input) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: input
  integer, allocatable        :: output(:,:)
  
  output = input%contents_
end function

pure function nint_RealVector(input) result(output)
  implicit none
  
  type(RealVector), intent(in) :: input
  integer, allocatable         :: output(:)
  
  output = nint(input%contents_)
end function

pure function nint_RealMatrix(input) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: input
  integer, allocatable         :: output(:,:)
  
  output = nint(input%contents_)
end function

pure function dble_IntVector(input) result(output)
  implicit none
  
  type(IntVector), intent(in) :: input
  real(dp), allocatable       :: output(:)
  
  output = int(input)
end function

pure function dble_IntMatrix(input) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: input
  real(dp), allocatable       :: output(:,:)
  
  output = int(input)
end function

pure function dble_RealVector(input) result(output)
  implicit none
  
  type(RealVector), intent(in) :: input
  real(dp), allocatable        :: output(:)
  
  output = input%contents_
end function

pure function dble_RealMatrix(input) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: input
  real(dp), allocatable        :: output(:,:)
  
  output = input%contents_
end function

pure function cmplx_IntVector(input) result(output)
  implicit none
  
  type(IntVector), intent(in) :: input
  complex(dp), allocatable    :: output(:)
  
  output = int(input)
end function

pure function cmplx_IntMatrix(input) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: input
  complex(dp), allocatable    :: output(:,:)
  
  output = int(input)
end function

pure function cmplx_RealVector(input) result(output)
  implicit none
  
  type(RealVector), intent(in) :: input
  complex(dp), allocatable     :: output(:)
  
  output = dble(input)
end function

pure function cmplx_RealMatrix(input) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: input
  complex(dp), allocatable     :: output(:,:)
  
  output = dble(input)
end function

pure function cmplx_ComplexVector(input) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: input
  complex(dp), allocatable        :: output(:)
  
  output = input%contents_
end function

pure function cmplx_ComplexMatrix(input) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: input
  complex(dp), allocatable        :: output(:,:)
  
  output = input%contents_
end function

pure function cmplx_RealVectors(real,imag) result(output)
  implicit none
  
  type(RealVector), intent(in) :: real
  type(RealVector), intent(in) :: imag
  complex(dp), allocatable     :: output(:)
  
  output = cmplx(dble(real),dble(imag),dp)
end function

pure function cmplx_RealMatrices(real,imag) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: real
  type(RealMatrix), intent(in) :: imag
  complex(dp), allocatable     :: output(:,:)
  
  output = cmplx(dble(real),dble(imag),dp)
end function

! Real part of a complex object.
pure function real_ComplexVector(input) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: input
  type(RealVector)                :: output
  
  output = real(input%contents_,dp)
end function

pure function real_ComplexMatrix(input) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: input
  type(RealMatrix)                :: output
  
  output = real(input%contents_,dp)
end function

! Imaginary part of a complex object.
pure function aimag_ComplexVector(input) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: input
  type(RealVector)                :: output
  
  output = aimag(input%contents_)
end function

pure function aimag_ComplexMatrix(input) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: input
  type(RealMatrix)                :: output
  
  output = aimag(input%contents_)
end function

! Conjugate of a complex object.
pure function conjg_ComplexVector(input) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: input
  type(ComplexVector)             :: output
  
  output = conjg(input%contents_)
end function

pure function conjg_ComplexMatrix(input) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: input
  type(ComplexMatrix)             :: output
  
  output = conjg(input%contents_)
end function

! size().
pure function size_IntVector(input) result(output)
  implicit none
  
  type(IntVector), intent(in) :: input
  integer                     :: output
  
  output = size(input%contents_)
end function

pure function size_RealVector(input) result(output)
  implicit none
  
  type(RealVector), intent(in) :: input
  integer                      :: output
  
  output = size(input%contents_)
end function

pure function size_ComplexVector(input) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: input
  integer                         :: output
  
  output = size(input%contents_)
end function

pure function size_IntMatrix(input,dim) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: input
  integer,         intent(in) :: dim
  integer                     :: output
  
  output = size(input%contents_, dim)
end function

pure function size_RealMatrix(input,dim) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: input
  integer,          intent(in) :: dim
  integer                      :: output
  
  output = size(input%contents_, dim)
end function

pure function size_ComplexMatrix(input,dim) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: input
  integer,             intent(in) :: dim
  integer                         :: output
  
  output = size(input%contents_, dim)
end function

! Trace
pure function trace_IntMatrix(input) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: input
  integer                     :: output
  
  integer :: i
  
  output = 0
  do i=1,size(input,1)
    output = output + input%contents_(i,i)
  enddo
end function

pure function trace_RealMatrix(input) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: input
  real(dp)                     :: output
  
  integer :: i
  
  output = 0
  do i=1,size(input,1)
    output = output + input%contents_(i,i)
  enddo
end function

pure function trace_ComplexMatrix(input) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: input
  complex(dp)                     :: output
  
  integer :: i
  
  output = 0
  do i=1,size(input,1)
    output = output + input%contents_(i,i)
  enddo
end function

! Equality.
elemental function equality_IntVector_IntVector(a,b) result(output)
  implicit none
  
  type(IntVector), intent(in) :: a
  type(IntVector), intent(in) :: b
  logical                     :: output
  
  output = all(a%contents_==b%contents_)
end function

elemental function equality_IntMatrix_IntMatrix(a,b) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: a
  type(IntMatrix), intent(in) :: b
  logical                     :: output
  
  output = all(a%contents_==b%contents_)
end function

! Non-equality.
elemental function non_equality_IntVector_IntVector(a,b) result(output)
  implicit none
  
  type(IntVector), intent(in) :: a
  type(IntVector), intent(in) :: b
  logical                     :: output
  
  output = .not. a==b
end function

elemental function non_equality_IntMatrix_IntMatrix(a,b) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: a
  type(IntMatrix), intent(in) :: b
  logical                     :: output
  
  output = .not. a==b
end function

! Addition.
elemental function add_IntVector_IntVector(a,b) result(output)
  implicit none
  
  type(IntVector), intent(in) :: a
  type(IntVector), intent(in) :: b
  type(IntVector)             :: output
  
  output = a%contents_ + b%contents_
end function

elemental function add_IntVector_RealVector(a,b) result(output)
  implicit none
  
  type(IntVector),  intent(in) :: a
  type(RealVector), intent(in) :: b
  type(RealVector)             :: output
  
  output = a%contents_ + b%contents_
end function

elemental function add_RealVector_IntVector(a,b) result(output)
  implicit none
  
  type(RealVector), intent(in) :: a
  type(IntVector),  intent(in) :: b
  type(RealVector)             :: output
  
  output = a%contents_ + b%contents_
end function

elemental function add_RealVector_RealVector(a,b) result(output)
  implicit none
  
  type(RealVector), intent(in) :: a
  type(RealVector), intent(in) :: b
  type(RealVector)             :: output
  
  output = a%contents_ + b%contents_
end function

elemental function add_RealVector_ComplexVector(a,b) result(output)
  implicit none
  
  type(RealVector),    intent(in) :: a
  type(ComplexVector), intent(in) :: b
  type(ComplexVector)             :: output
  
  output = a%contents_ + b%contents_
end function

elemental function add_ComplexVector_RealVector(a,b) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: a
  type(RealVector),    intent(in) :: b
  type(ComplexVector)             :: output
  
  output = a%contents_ + b%contents_
end function

elemental function add_ComplexVector_ComplexVector(a,b) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: a
  type(ComplexVector), intent(in) :: b
  type(ComplexVector)             :: output
  
  output = a%contents_ + b%contents_
end function

elemental function add_IntMatrix_IntMatrix(a,b) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: a
  type(IntMatrix), intent(in) :: b
  type(IntMatrix)             :: output
  
  output = a%contents_ + b%contents_
end function

elemental function add_IntMatrix_RealMatrix(a,b) result(output)
  implicit none
  
  type(IntMatrix),  intent(in) :: a
  type(RealMatrix), intent(in) :: b
  type(RealMatrix)             :: output
  
  output = a%contents_ + b%contents_
end function

elemental function add_RealMatrix_IntMatrix(a,b) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: a
  type(IntMatrix),  intent(in) :: b
  type(RealMatrix)             :: output
  
  output = a%contents_ + b%contents_
end function

elemental function add_RealMatrix_RealMatrix(a,b) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: a
  type(RealMatrix), intent(in) :: b
  type(RealMatrix)             :: output
  
  output = a%contents_ + b%contents_
end function

elemental function add_RealMatrix_ComplexMatrix(a,b) result(output)
  implicit none
  
  type(RealMatrix),    intent(in) :: a
  type(ComplexMatrix), intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = a%contents_ + b%contents_
end function

elemental function add_ComplexMatrix_RealMatrix(a,b) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: a
  type(RealMatrix),    intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = a%contents_ + b%contents_
end function

elemental function add_ComplexMatrix_ComplexMatrix(a,b) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: a
  type(ComplexMatrix), intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = a%contents_ + b%contents_
end function

! Negative.
elemental function negative_IntVector(input) result(output)
  implicit none
  
  type(IntVector), intent(in) :: input
  type(IntVector)             :: output
  
  output = -input%contents_
end function

elemental function negative_RealVector(input) result(output)
  implicit none
  
  type(RealVector), intent(in) :: input
  type(RealVector)             :: output
  
  output = -input%contents_
end function

elemental function negative_ComplexVector(input) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: input
  type(ComplexVector)             :: output
  
  output = -input%contents_
end function

elemental function negative_IntMatrix(input) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: input
  type(IntMatrix)             :: output
  
  output = -input%contents_
end function

elemental function negative_RealMatrix(input) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: input
  type(RealMatrix)             :: output
  
  output = -input%contents_
end function

elemental function negative_ComplexMatrix(input) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: input
  type(ComplexMatrix)             :: output
  
  output = -input%contents_
end function

! Subtraction.
elemental function subtract_IntVector_IntVector(a,b) result(output)
  implicit none
  
  type(IntVector), intent(in) :: a
  type(IntVector), intent(in) :: b
  type(IntVector)             :: output
  
  output = a%contents_ - b%contents_
end function

elemental function subtract_IntVector_RealVector(a,b) result(output)
  implicit none
  
  type(IntVector),  intent(in) :: a
  type(RealVector), intent(in) :: b
  type(RealVector)             :: output
  
  output = a%contents_ - b%contents_
end function

elemental function subtract_RealVector_IntVector(a,b) result(output)
  implicit none
  
  type(RealVector), intent(in) :: a
  type(IntVector),  intent(in) :: b
  type(RealVector)             :: output
  
  output = a%contents_ - b%contents_
end function

elemental function subtract_RealVector_RealVector(a,b) result(output)
  implicit none
  
  type(RealVector), intent(in) :: a
  type(RealVector), intent(in) :: b
  type(RealVector)             :: output
  
  output = a%contents_ - b%contents_
end function

elemental function subtract_RealVector_ComplexVector(a,b) result(output)
  implicit none
  
  type(RealVector),    intent(in) :: a
  type(ComplexVector), intent(in) :: b
  type(ComplexVector)             :: output
  
  output = a%contents_ - b%contents_
end function

elemental function subtract_ComplexVector_RealVector(a,b) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: a
  type(RealVector),    intent(in) :: b
  type(ComplexVector)             :: output
  
  output = a%contents_ - b%contents_
end function

elemental function subtract_ComplexVector_ComplexVector(a,b) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: a
  type(ComplexVector), intent(in) :: b
  type(ComplexVector)             :: output
  
  output = a%contents_ - b%contents_
end function

elemental function subtract_IntMatrix_IntMatrix(a,b) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: a
  type(IntMatrix), intent(in) :: b
  type(IntMatrix)             :: output
  
  output = a%contents_ - b%contents_
end function

elemental function subtract_IntMatrix_RealMatrix(a,b) result(output)
  implicit none
  
  type(IntMatrix),  intent(in) :: a
  type(RealMatrix), intent(in) :: b
  type(RealMatrix)             :: output
  
  output = a%contents_ - b%contents_
end function

elemental function subtract_RealMatrix_IntMatrix(a,b) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: a
  type(IntMatrix),  intent(in) :: b
  type(RealMatrix)             :: output
  
  output = a%contents_ - b%contents_
end function

elemental function subtract_RealMatrix_RealMatrix(a,b) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: a
  type(RealMatrix), intent(in) :: b
  type(RealMatrix)             :: output
  
  output = a%contents_ - b%contents_
end function

elemental function subtract_RealMatrix_ComplexMatrix(a,b) result(output)
  implicit none
  
  type(RealMatrix),    intent(in) :: a
  type(ComplexMatrix), intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = a%contents_ - b%contents_
end function

elemental function subtract_ComplexMatrix_RealMatrix(a,b) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: a
  type(RealMatrix),    intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = a%contents_ - b%contents_
end function

elemental function subtract_ComplexMatrix_ComplexMatrix(a,b) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: a
  type(ComplexMatrix), intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = a%contents_ - b%contents_
end function

! Multiplication by a scalar.
pure function multiply_IntVector_integer(a,b) result(output)
  implicit none
  
  type(IntVector), intent(in) :: a
  integer,         intent(in) :: b
  type(IntVector)             :: output
  
  output = a%contents_*b
end function

pure function multiply_integer_IntVector(a,b) result(output)
  implicit none
  
  integer,         intent(in) :: a
  type(IntVector), intent(in) :: b
  type(IntVector)             :: output
  
  output = a*b%contents_
end function

pure function multiply_IntVector_real(a,b) result(output)
  implicit none
  
  type(IntVector), intent(in) :: a
  real(dp),        intent(in) :: b
  type(RealVector)            :: output
  
  output = a%contents_*b
end function

pure function multiply_real_IntVector(a,b) result(output)
  implicit none
  
  real(dp),        intent(in) :: a
  type(IntVector), intent(in) :: b
  type(RealVector)            :: output
  
  output = a*b%contents_
end function

pure function multiply_IntVector_complex(a,b) result(output)
  implicit none
  
  type(IntVector), intent(in) :: a
  complex(dp),     intent(in) :: b
  type(ComplexVector)         :: output
  
  output = a%contents_*b
end function

pure function multiply_complex_IntVector(a,b) result(output)
  implicit none
  
  complex(dp),     intent(in) :: a
  type(IntVector), intent(in) :: b
  type(ComplexVector)         :: output
  
  output = a*b%contents_
end function

pure function multiply_RealVector_integer(a,b) result(output)
  implicit none
  
  type(RealVector), intent(in) :: a
  integer,          intent(in) :: b
  type(RealVector)             :: output
  
  output = a%contents_*b
end function

pure function multiply_integer_RealVector(a,b) result(output)
  implicit none
  
  integer,          intent(in) :: a
  type(RealVector), intent(in) :: b
  type(RealVector)             :: output
  
  output = a*b%contents_
end function

pure function multiply_RealVector_real(a,b) result(output)
  implicit none
  
  type(RealVector), intent(in) :: a
  real(dp),         intent(in) :: b
  type(RealVector)             :: output
  
  output = a%contents_*b
end function

pure function multiply_real_RealVector(a,b) result(output)
  implicit none
  
  real(dp),         intent(in) :: a
  type(RealVector), intent(in) :: b
  type(RealVector)             :: output
  
  output = a*b%contents_
end function

pure function multiply_RealVector_complex(a,b) result(output)
  implicit none
  
  type(RealVector), intent(in) :: a
  complex(dp),      intent(in) :: b
  type(ComplexVector)          :: output
  
  output = a%contents_*b
end function

pure function multiply_complex_RealVector(a,b) result(output)
  implicit none
  
  complex(dp),      intent(in) :: a
  type(RealVector), intent(in) :: b
  type(ComplexVector)          :: output
  
  output = a*b%contents_
end function

pure function multiply_ComplexVector_integer(a,b) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: a
  integer,             intent(in) :: b
  type(ComplexVector)             :: output
  
  output = a%contents_*b
end function

pure function multiply_integer_ComplexVector(a,b) result(output)
  implicit none
  
  integer,             intent(in) :: a
  type(ComplexVector), intent(in) :: b
  type(ComplexVector)             :: output
  
  output = a*b%contents_
end function

pure function multiply_ComplexVector_real(a,b) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: a
  real(dp),            intent(in) :: b
  type(ComplexVector)             :: output
  
  output = a%contents_*b
end function

pure function multiply_real_ComplexVector(a,b) result(output)
  implicit none
  
  real(dp),            intent(in) :: a
  type(ComplexVector), intent(in) :: b
  type(ComplexVector)             :: output
  
  output = a*b%contents_
end function

pure function multiply_ComplexVector_complex(a,b) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: a
  complex(dp),         intent(in) :: b
  type(ComplexVector)             :: output
  
  output = a%contents_*b
end function

pure function multiply_complex_ComplexVector(a,b) result(output)
  implicit none
  
  complex(dp),         intent(in) :: a
  type(ComplexVector), intent(in) :: b
  type(ComplexVector)             :: output
  
  output = a*b%contents_
end function

pure function multiply_IntMatrix_integer(a,b) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: a
  integer,         intent(in) :: b
  type(IntMatrix)             :: output
  
  output = a%contents_*b
end function

pure function multiply_integer_IntMatrix(a,b) result(output)
  implicit none
  
  integer,         intent(in) :: a
  type(IntMatrix), intent(in) :: b
  type(IntMatrix)             :: output
  
  output = a*b%contents_
end function

pure function multiply_IntMatrix_real(a,b) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: a
  real(dp),        intent(in) :: b
  type(RealMatrix)            :: output
  
  output = a%contents_*b
end function

pure function multiply_real_IntMatrix(a,b) result(output)
  implicit none
  
  real(dp),        intent(in) :: a
  type(IntMatrix), intent(in) :: b
  type(RealMatrix)            :: output
  
  output = a*b%contents_
end function

pure function multiply_IntMatrix_complex(a,b) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: a
  complex(dp),     intent(in) :: b
  type(ComplexMatrix)         :: output
  
  output = a%contents_*b
end function

pure function multiply_complex_IntMatrix(a,b) result(output)
  implicit none
  
  complex(dp),     intent(in) :: a
  type(IntMatrix), intent(in) :: b
  type(ComplexMatrix)         :: output
  
  output = a*b%contents_
end function

pure function multiply_RealMatrix_integer(a,b) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: a
  integer,          intent(in) :: b
  type(RealMatrix)             :: output
  
  output = a%contents_*b
end function

pure function multiply_integer_RealMatrix(a,b) result(output)
  implicit none
  
  integer,          intent(in) :: a
  type(RealMatrix), intent(in) :: b
  type(RealMatrix)             :: output
  
  output = a*b%contents_
end function

pure function multiply_RealMatrix_real(a,b) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: a
  real(dp),         intent(in) :: b
  type(RealMatrix)             :: output
  
  output = a%contents_*b
end function

pure function multiply_real_RealMatrix(a,b) result(output)
  implicit none
  
  real(dp),         intent(in) :: a
  type(RealMatrix), intent(in) :: b
  type(RealMatrix)             :: output
  
  output = a*b%contents_
end function

pure function multiply_RealMatrix_complex(a,b) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: a
  complex(dp),      intent(in) :: b
  type(ComplexMatrix)          :: output
  
  output = a%contents_*b
end function

pure function multiply_complex_RealMatrix(a,b) result(output)
  implicit none
  
  complex(dp),      intent(in) :: a
  type(RealMatrix), intent(in) :: b
  type(ComplexMatrix)          :: output
  
  output = a*b%contents_
end function

pure function multiply_ComplexMatrix_integer(a,b) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: a
  integer,             intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = a%contents_*b
end function

pure function multiply_integer_ComplexMatrix(a,b) result(output)
  implicit none
  
  integer,             intent(in) :: a
  type(ComplexMatrix), intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = a*b%contents_
end function

pure function multiply_ComplexMatrix_real(a,b) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: a
  real(dp),            intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = a%contents_*b
end function

pure function multiply_real_ComplexMatrix(a,b) result(output)
  implicit none
  
  real(dp),            intent(in) :: a
  type(ComplexMatrix), intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = a*b%contents_
end function

pure function multiply_ComplexMatrix_complex(a,b) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: a
  complex(dp),         intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = a%contents_*b
end function

pure function multiply_complex_ComplexMatrix(a,b) result(output)
  implicit none
  
  complex(dp),         intent(in) :: a
  type(ComplexMatrix), intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = a*b%contents_
end function

! Dot products and matrix multiplication.
pure function dot_IntVector_IntVector(a,b) result(output)
  implicit none
  
  type(IntVector), intent(in) :: a
  type(IntVector), intent(in) :: b
  integer                     :: output
  
  output = dot_product(a%contents_, b%contents_)
end function

pure function dot_IntVector_RealVector(a,b) result(output)
  implicit none
  
  type(IntVector),  intent(in) :: a
  type(RealVector), intent(in) :: b
  real(dp)                     :: output
  
  output = dot_product(a%contents_, b%contents_)
end function

pure function dot_RealVector_IntVector(a,b) result(output)
  implicit none
  
  type(RealVector), intent(in) :: a
  type(IntVector),  intent(in) :: b
  real(dp)                     :: output
  
  output = dot_product(a%contents_, b%contents_)
end function

pure function dot_RealVector_RealVector(a,b) result(output)
  implicit none
  
  type(RealVector), intent(in) :: a
  type(RealVector), intent(in) :: b
  real(dp)                     :: output
  
  output = dot_product(a%contents_, b%contents_)
end function

pure function dot_RealVector_ComplexVector(a,b) result(output)
  implicit none
  
  type(RealVector),    intent(in) :: a
  type(ComplexVector), intent(in) :: b
  complex(dp)                     :: output
  
  output = dot_product(a%contents_, b%contents_)
end function

pure function dot_ComplexVector_RealVector(a,b) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: a
  type(RealVector),    intent(in) :: b
  complex(dp)                     :: output
  
  output = dot_product(a%contents_, b%contents_)
end function

pure function dot_ComplexVector_ComplexVector(a,b) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: a
  type(ComplexVector), intent(in) :: b
  complex(dp)                     :: output
  
  output = dot_product(a%contents_, b%contents_)
end function

pure function dot_IntVector_IntMatrix(a,b) result(output)
  implicit none
  
  type(IntVector), intent(in) :: a
  type(IntMatrix), intent(in) :: b
  type(IntVector)             :: output
  
  output = matmul(a%contents_, b%contents_)
end function

pure function dot_IntVector_RealMatrix(a,b) result(output)
  implicit none
  
  type(IntVector),  intent(in) :: a
  type(RealMatrix), intent(in) :: b
  type(RealVector)             :: output
  
  output = matmul(a%contents_, b%contents_)
end function

pure function dot_RealVector_IntMatrix(a,b) result(output)
  implicit none
  
  type(RealVector), intent(in) :: a
  type(IntMatrix),  intent(in) :: b
  type(RealVector)             :: output
  
  output = matmul(a%contents_, b%contents_)
end function

pure function dot_RealVector_RealMatrix(a,b) result(output)
  implicit none
  
  type(RealVector), intent(in) :: a
  type(RealMatrix), intent(in) :: b
  type(RealVector)             :: output
  
  output = matmul(a%contents_, b%contents_)
end function

pure function dot_RealVector_ComplexMatrix(a,b) result(output)
  implicit none
  
  type(RealVector),    intent(in) :: a
  type(ComplexMatrix), intent(in) :: b
  type(ComplexVector)             :: output
  
  output = matmul(a%contents_, b%contents_)
end function

pure function dot_ComplexVector_RealMatrix(a,b) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: a
  type(RealMatrix),    intent(in) :: b
  type(ComplexVector)             :: output
  
  output = matmul(a%contents_, b%contents_)
end function

pure function dot_ComplexVector_ComplexMatrix(a,b) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: a
  type(ComplexMatrix), intent(in) :: b
  type(ComplexVector)             :: output
  
  output = matmul(a%contents_, b%contents_)
end function

pure function dot_IntMatrix_IntVector(a,b) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: a
  type(IntVector), intent(in) :: b
  type(IntVector)             :: output
  
  output = matmul(a%contents_, b%contents_)
end function

pure function dot_IntMatrix_RealVector(a,b) result(output)
  implicit none
  
  type(IntMatrix),  intent(in) :: a
  type(RealVector), intent(in) :: b
  type(RealVector)             :: output
  
  output = matmul(a%contents_, b%contents_)
end function

pure function dot_RealMatrix_IntVector(a,b) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: a
  type(IntVector),  intent(in) :: b
  type(RealVector)             :: output
  
  output = matmul(a%contents_, b%contents_)
end function

pure function dot_RealMatrix_RealVector(a,b) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: a
  type(RealVector), intent(in) :: b
  type(RealVector)             :: output
  
  output = matmul(a%contents_, b%contents_)
end function

pure function dot_RealMatrix_ComplexVector(a,b) result(output)
  implicit none
  
  type(RealMatrix),    intent(in) :: a
  type(ComplexVector), intent(in) :: b
  type(ComplexVector)             :: output
  
  output = matmul(a%contents_, b%contents_)
end function

pure function dot_ComplexMatrix_RealVector(a,b) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: a
  type(RealVector),    intent(in) :: b
  type(ComplexVector)             :: output
  
  output = matmul(a%contents_, b%contents_)
end function

pure function dot_ComplexMatrix_ComplexVector(a,b) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: a
  type(ComplexVector), intent(in) :: b
  type(ComplexVector)             :: output
  
  output = matmul(a%contents_, b%contents_)
end function

pure function dot_IntMatrix_IntMatrix(a,b) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: a
  type(IntMatrix), intent(in) :: b
  type(IntMatrix)             :: output
  
  output = matmul(a%contents_, b%contents_)
end function

pure function dot_IntMatrix_RealMatrix(a,b) result(output)
  implicit none
  
  type(IntMatrix),  intent(in) :: a
  type(RealMatrix), intent(in) :: b
  type(RealMatrix)             :: output
  
  output = matmul(a%contents_, b%contents_)
end function

pure function dot_RealMatrix_IntMatrix(a,b) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: a
  type(IntMatrix),  intent(in) :: b
  type(RealMatrix)             :: output
  
  output = matmul(a%contents_, b%contents_)
end function

pure function dot_RealMatrix_RealMatrix(a,b) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: a
  type(RealMatrix), intent(in) :: b
  type(RealMatrix)             :: output
  
  output = matmul(a%contents_, b%contents_)
end function

pure function dot_RealMatrix_ComplexMatrix(a,b) result(output)
  implicit none
  
  type(RealMatrix),    intent(in) :: a
  type(ComplexMatrix), intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = matmul(a%contents_, b%contents_)
end function

pure function dot_ComplexMatrix_RealMatrix(a,b) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: a
  type(RealMatrix),    intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = matmul(a%contents_, b%contents_)
end function

pure function dot_ComplexMatrix_ComplexMatrix(a,b) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: a
  type(ComplexMatrix), intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = matmul(a%contents_, b%contents_)
end function

! Division by scalar.
pure function divide_IntVector_integer(a,b) result(output)
  implicit none
  
  type(IntVector), intent(in) :: a
  integer,         intent(in) :: b
  type(IntVector)             :: output
  
  output = a%contents_/b
end function

pure function divide_IntVector_real(a,b) result(output)
  implicit none
  
  type(IntVector), intent(in) :: a
  real(dp),        intent(in) :: b
  type(RealVector)            :: output
  
  output = a%contents_/b
end function

pure function divide_IntVector_complex(a,b) result(output)
  implicit none
  
  type(IntVector), intent(in) :: a
  complex(dp),     intent(in) :: b
  type(ComplexVector)         :: output
  
  output = a%contents_/b
end function

pure function divide_RealVector_integer(a,b) result(output)
  implicit none
  
  type(RealVector), intent(in) :: a
  integer,          intent(in) :: b
  type(RealVector)             :: output
  
  output = a%contents_/b
end function

pure function divide_RealVector_real(a,b) result(output)
  implicit none
  
  type(RealVector), intent(in) :: a
  real(dp),         intent(in) :: b
  type(RealVector)             :: output
  
  output = a%contents_/b
end function

pure function divide_RealVector_complex(a,b) result(output)
  implicit none
  
  type(RealVector),    intent(in) :: a
  complex(dp),         intent(in) :: b
  type(ComplexVector)             :: output
  
  output = a%contents_/b
end function

pure function divide_ComplexVector_integer(a,b) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: a
  integer,             intent(in) :: b
  type(ComplexVector)             :: output
  
  output = a%contents_/b
end function

pure function divide_ComplexVector_real(a,b) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: a
  real(dp),            intent(in) :: b
  type(ComplexVector)             :: output
  
  output = a%contents_/b
end function

pure function divide_ComplexVector_complex(a,b) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: a
  complex(dp),         intent(in) :: b
  type(ComplexVector)             :: output
  
  output = a%contents_/b
end function

pure function divide_IntMatrix_integer(a,b) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: a
  integer,         intent(in) :: b
  type(IntMatrix)             :: output
  
  output = a%contents_/b
end function

pure function divide_IntMatrix_real(a,b) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: a
  real(dp),        intent(in) :: b
  type(RealMatrix)            :: output
  
  output = a%contents_/b
end function

pure function divide_IntMatrix_complex(a,b) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: a
  complex(dp),     intent(in) :: b
  type(ComplexMatrix)         :: output
  
  output = a%contents_/b
end function

pure function divide_RealMatrix_integer(a,b) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: a
  integer,          intent(in) :: b
  type(RealMatrix)             :: output
  
  output = a%contents_/b
end function

pure function divide_RealMatrix_real(a,b) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: a
  real(dp),         intent(in) :: b
  type(RealMatrix)             :: output
  
  output = a%contents_/b
end function

pure function divide_RealMatrix_complex(a,b) result(output)
  implicit none
  
  type(RealMatrix),    intent(in) :: a
  complex(dp),         intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = a%contents_/b
end function

pure function divide_ComplexMatrix_integer(a,b) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: a
  integer,             intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = a%contents_/b
end function

pure function divide_ComplexMatrix_real(a,b) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: a
  real(dp),            intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = a%contents_/b
end function

pure function divide_ComplexMatrix_complex(a,b) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: a
  complex(dp),         intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = a%contents_/b
end function

! L2 norm.
elemental function l2_norm_RealVector(input) result(output)
  implicit none
  
  type(RealVector), intent(in) :: input
  real(dp)                     :: output
  
  output = sqrt(input*input)
end function

elemental function l2_norm_ComplexVector(input) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: input
  real(dp)                        :: output
  
  output = sqrt(real(input*conjg(input)))
end function

! Outer product.
function outer_product_RealVector_RealVector(a,b) result(output)
  implicit none
  
  type(RealVector), intent(in) :: a
  type(RealVector), intent(in) :: b
  type(RealMatrix)             :: output
  
  integer               :: i,ialloc
  
  allocate(output%contents_(size(a),size(b)), stat=ialloc); call err(ialloc)
  do i=1,size(b)
    output%contents_(:,i) = a%contents_*b%contents_(i)
  enddo
end function

! Transpose.
pure function transpose_IntMatrix(input) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: input
  type(IntMatrix)             :: output
  
  output = transpose(input%contents_)
end function

pure function transpose_RealMatrix(input) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: input
  type(RealMatrix)             :: output
  
  output = transpose(input%contents_)
end function

pure function transpose_ComplexMatrix(input) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: input
  type(ComplexMatrix)             :: output
  
  output = transpose(input%contents_)
end function

! Hermitian conjugate.
pure function hermitian_ComplexMatrix(input) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: input
  type(ComplexMatrix)             :: output
  
  output = transpose(conjg(input%contents_))
end function

! Determinant. (3x3 only)
function determinant_integer(A) result(determinant)
  implicit none
  
  integer, intent(in) :: A(:,:)
  integer             :: determinant
  
  ! Check that A is a 3x3 matrix.
  if (size(A,1)/=3 .or. size(A,2)/=3) then
    call err()
  endif
  
  determinant = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))&
             &+ A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))&
             &+ A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
end function

function determinant_real(A) result(determinant)
  implicit none
  
  real(dp), intent(in) :: A(:,:)
  real(dp)             :: determinant
  
  ! Check that A is a 3x3 matrix.
  if (size(A,1)/=3 .or. size(A,2)/=3) then
    call err()
  endif
  
  determinant = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))&
            & + A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))&
            & + A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
end function

function determinant_IntMatrix(input) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: input
  integer                     :: output
  
  output = determinant(input%contents_)
end function

function determinant_RealMatrix(input) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: input
  real(dp)                     :: output
  
  output = determinant(input%contents_)
end function

! calculates the 3x3 matrix inverse.
function invert_reals(input) result(output)
  implicit none
  
  real(dp), intent(in)  :: input(:,:)
  real(dp)              :: output(3,3)
  
  real(dp) :: inverse_determinant
  
  ! Check that the input is a 3x3 matrix.
  if (size(input,1)/=3 .or. size(input,2)/=3) then
    call err()
  endif
  
  inverse_determinant = 1/determinant(input)
  
  ! check for inverse_det=infinity.
  if (abs(inverse_determinant)>huge(0.0_dp)) then
    call print_line('Error in invert: singular matrix.')
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
  
  output = output * inverse_determinant
end function

function invert_RealMatrix(input) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: input
  type(RealMatrix)             :: output
  
  output = invert(input%contents_)
end function

! I/O overloads.
pure function str_IntVector(this) result(output)
  implicit none
  
  class(IntVector), intent(in) :: this
  type(String)                 :: output
  
  output = join(this%contents_)
end function

pure function str_RealVector(this) result(output)
  implicit none
  
  class(RealVector), intent(in) :: this
  type(String)                  :: output
  
  output = join(this%contents_)
end function

pure function str_ComplexVector(this) result(output)
  implicit none
  
  class(ComplexVector), intent(in) :: this
  type(String)                     :: output
  
  output = join(this%contents_)
end function

function str_IntMatrix(this) result(output)
  implicit none
  
  Class(IntMatrix), intent(in) :: this
  type(String), allocatable    :: output(:)
  
  integer :: i,ialloc
  
  allocate(output(size(this,1)), stat=ialloc); call err(ialloc)
  
  do i=1,size(this,1)
    output(i) = join(this%contents_(i,:))
  enddo
end function

function str_RealMatrix(this) result(output)
  implicit none
  
  Class(RealMatrix), intent(in) :: this
  type(String), allocatable     :: output(:)
  
  integer :: i,ialloc
  
  allocate(output(size(this,1)), stat=ialloc); call err(ialloc)
  
  do i=1,size(this,1)
    output(i) = join(this%contents_(i,:))
  enddo
end function

function str_ComplexMatrix(this) result(output)
  implicit none
  
  Class(ComplexMatrix), intent(in) :: this
  type(String), allocatable        :: output(:)
  
  integer :: i,ialloc
  
  allocate(output(size(this,1)), stat=ialloc); call err(ialloc)
  
  do i=1,size(this,1)
    output(i) = join(this%contents_(i,:))
  enddo
end function

! ----------------------------------------------------------------------
! Eigenvalue and Eigenvector wrappers.
! ----------------------------------------------------------------------

! --------------------------------------------------
! Returns the number of states of an Eigenstuff.
! --------------------------------------------------
pure function size_RealEigenstuff(estuff) result(output)
  implicit none
  
  type(RealEigenstuff), intent(in) :: estuff
  integer                          :: output
  
  output = size(estuff%evals)
end function

pure function size_ComplexEigenstuff(estuff) result(output)
  implicit none
  
  type(ComplexEigenstuff), intent(in) :: estuff
  integer                             :: output
  
  output = size(estuff%evals)
end function

! --------------------------------------------------
! Calculates the eigenvalues and eigenvectors of a real, symmetric matrix.
! --------------------------------------------------
function calculate_RealEigenstuff_reals(input) result(output)
  implicit none
  
  real(dp), intent(in) :: input(:,:)  ! A real, symmetric matrix.
  type(RealEigenstuff) :: output      ! The eigenvalues and eigenstates.
  
  ! dsyev variables.
  integer               :: n
  real(dp), allocatable :: work(:)
  integer               :: lwork
  integer               :: info
  
  ! Temporary variables.
  integer :: ialloc
  
  n = size(input,1)
  allocate( output%evals(n), &
          & work(3*n-1),     &
          & stat=ialloc); call err(ialloc)
  output%evecs = input
  
  ! calculate optimal lwork
  call dsyev('V', 'U', n, output%evecs(1,1), n, output%evals, &
    & work(1), -1, info)
  if (info /= 0) then
    call print_line("dsyev failed, info= "//info)
    call err()
  endif
  lwork = nint(work(1))
  deallocate(work, stat=ialloc); call err(ialloc)
  allocate(work(lwork), stat=ialloc); call err(ialloc)
  
  ! calculate eigenstuff
  call dsyev('V', 'U', n, output%evecs(1,1), n, output%evals, &
    & work(1), lwork, info)
  if (info /= 0) then
    call print_line("dsyev failed, info= "//info)
    call err()
  endif
end function

function calculate_RealEigenstuff_RealMatrix(input) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: input
  type(RealEigenstuff)         :: output
  
  output = calculate_eigenstuff(dble(input))
end function

! --------------------------------------------------
! Calculates the eigenvalues and eigenvectors of a complex, Hermitian matrix.
! --------------------------------------------------
function calculate_ComplexEigenstuff_complexes(input) result(output)
  implicit none
  
  complex(dp), intent(in) :: input(:,:)  ! a complex, hermitian matrix
  type(ComplexEigenstuff) :: output      ! the eigenvalues and eigenstates of a
  
  ! zheev variables.
  integer                  :: n
  complex(dp), allocatable :: work(:)
  integer                  :: lwork
  real(dp),    allocatable :: rwork(:)
  integer                  :: info
  
  ! Temporary variables.
  integer :: ialloc
  
  n = size(input,1)
  allocate( output%evals(n), &
          & rwork(3*n-2),    &
          & work(3*n-1),     &
          & stat=ialloc); call err(ialloc)
  output%evecs = input
  
  ! calculate optimal lwork
  call zheev('V', 'U', n, output%evecs(1,1), n, output%evals, &
    & work(1), -1, rwork, info)
  if (info /= 0) then
    call print_line("zheev failed, info= "//info)
    call err()
  endif
  lwork = nint(real(work(1)))
  deallocate(work, stat=ialloc); call err(ialloc)
  allocate(work(lwork), stat=ialloc); call err(ialloc)
  
  ! calculate eigenstuff
  call zheev('V', 'U', n, output%evecs(1,1), n, output%evals, &
    & work(1), lwork, rwork, info)
  if (info /= 0) then
    call print_line("zheev failed, info= "//info)
    call err()
  endif
end function

function calculate_ComplexEigenstuff_ComplexMatrix(input) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: input
  type(ComplexEigenstuff)         :: output
  
  output = calculate_eigenstuff(cmplx(input))
end function

! --------------------------------------------------
! Minimises the least-squares fit L=(a.x-b)**2.
! --------------------------------------------------
function linear_least_squares_reals_reals(a,b) result(x)
  implicit none
  
  real(dp), intent(in) :: a(:,:)
  real(dp), intent(in) :: b(:)
  type(RealVector)     :: x
  
  ! Array storage.
  real(dp), allocatable :: a2(:,:)
  
  ! dgels variables.
  integer               :: m,n
  real(dp), allocatable :: work(:)
  integer               :: lwork
  integer               :: info
  
  ! Temporary variables
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
  x = b
  
  ! Calculate optimal workspace size.
  allocate(work(2*m*n), stat=ialloc); call err(ialloc)
  call dgels('N',m,n,1,a2(1,1),m,x%contents_(1),m,work(1),-1,info)
  if (info/=0) then
    call print_line('dgels failed, info= '//info)
    call err()
  endif
  lwork = nint(work(1))
  deallocate(work, stat=ialloc); call err(ialloc)
  allocate(work(lwork), stat=ialloc); call err(ialloc)
  
  ! Run linear least-squares optimisation.
  call dgels('N',m,n,1,a2(1,1),m,x%contents_(1),m,work(1),lwork,info)
  if (info/=0) then
    call print_line('dgels failed, info= '//info)
    call err()
  endif
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
