! ======================================================================
! Assorted linear algebra / vector algebra subroutines.
! Includes interfaces for BLAS and LAPACK routines.
! ======================================================================
module linear_algebra_module
  use constants_module, only : dp
  use string_module
  use io_module
  implicit none
  
  ! --------------------------------------------------
  ! Vector and Matrix classes
  ! --------------------------------------------------
  type IntVector
    integer, allocatable, private :: contents(:)
  end type
  
  type RealVector
    real(dp), allocatable, private :: contents(:)
  end type
  
  type ComplexVector
    complex(dp), allocatable, private :: contents(:)
  end type
  
  type IntMatrix
    integer, allocatable, private :: contents(:,:)
  end type
  
  type RealMatrix
    real(dp), allocatable, private :: contents(:,:)
  end type
  
  type ComplexMatrix
    complex(dp), allocatable, private :: contents(:,:)
  end type
  
  ! --------------------------------------------------
  ! Overloads for vector / matrix operations.
  ! --------------------------------------------------
  interface assignment(=)
    module procedure assign_IntVector_integers
    module procedure assign_RealVector_reals
    module procedure assign_ComplexVector_complexes
    module procedure assign_IntMatrix_integers
    module procedure assign_RealMatrix_reals
    module procedure assign_ComplexMatrix_complexes
  end interface
  
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
  
  interface identity
    module procedure identity_IntMatrix
  end interface
  
  interface int
    module procedure int_IntVector
    module procedure int_IntMatrix
  end interface
  
  interface dble
    module procedure dble_RealVector
    module procedure dble_RealMatrix
  end interface
  
  interface cmplx
    module procedure cmplx_ComplexVector
    module procedure cmplx_ComplexMatrix
  end interface
  
  interface real
    module procedure real_ComplexVector
    module procedure real_ComplexMatrix
  end interface
  
  interface imag
    module procedure imag_ComplexVector
    module procedure imag_ComplexMatrix
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
    
    module procedure multiply_RealVector_integer
    module procedure multiply_integer_RealVector
    module procedure multiply_RealVector_real
    module procedure multiply_real_RealVector
    module procedure multiply_RealVector_complex
    module procedure multiply_complex_RealVector
    
    module procedure multiply_ComplexVector_real
    module procedure multiply_real_ComplexVector
    module procedure multiply_ComplexVector_complex
    module procedure multiply_complex_ComplexVector
    
    module procedure multiply_IntMatrix_integer
    module procedure multiply_integer_IntMatrix
    module procedure multiply_IntMatrix_real
    module procedure multiply_real_IntMatrix
    
    module procedure multiply_RealMatrix_integer
    module procedure multiply_integer_RealMatrix
    module procedure multiply_RealMatrix_real
    module procedure multiply_real_RealMatrix
    module procedure multiply_RealMatrix_complex
    module procedure multiply_complex_RealMatrix
    
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
    module procedure divide_RealVector_integer
    module procedure divide_RealVector_real
    module procedure divide_RealVector_complex
    module procedure divide_ComplexVector_real
    module procedure divide_ComplexVector_complex
    
    module procedure divide_IntMatrix_integer
    module procedure divide_IntMatrix_real
    module procedure divide_RealMatrix_integer
    module procedure divide_RealMatrix_real
    module procedure divide_RealMatrix_complex
    module procedure divide_ComplexMatrix_real
    module procedure divide_ComplexMatrix_complex
  end interface
  
  interface l2_norm
    module procedure l2_norm_RealVector
  end interface
  
  interface outer_product
    module procedure outer_product_RealVector_RealVector
  end interface
  
  interface transpose
    module procedure transpose_IntMatrix
    module procedure transpose_RealMatrix
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
    module procedure invert_real
    module procedure invert_RealMatrix
  end interface
  
  interface invert_int
    module procedure invert_int_integer
    module procedure invert_int_IntMatrix
  end interface
  
  ! --------------------------------------------------
  ! Overloads for printing vectors and matrices.
  ! --------------------------------------------------
  interface operator(//)
    module procedure concatenate_character_IntVector
    module procedure concatenate_IntVector_character
    module procedure concatenate_character_RealVector
    module procedure concatenate_RealVector_character
    module procedure concatenate_String_IntVector
    module procedure concatenate_IntVector_String
    module procedure concatenate_String_RealVector
    module procedure concatenate_RealVector_String
  end interface
  
  interface print_line
    module procedure print_line_IntVector
    module procedure print_line_IntVector_file
    module procedure print_line_RealVector
    module procedure print_line_RealVector_file
    module procedure print_line_ComplexVector
    module procedure print_line_ComplexVector_file
    module procedure print_line_IntMatrix
    module procedure print_line_IntMatrix_file
    module procedure print_line_RealMatrix
    module procedure print_line_RealMatrix_file
    module procedure print_line_ComplexMatrix
    module procedure print_line_ComplexMatrix_file
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
    pure real(dp) function ddot(N,DX,INCX,DY,INCY)
      import :: dp
      implicit none
      
      integer,  intent(in) :: N     ! length of vectors
      real(dp), intent(in) :: DX(*) ! first vector
      integer,  intent(in) :: INCX  ! increment along DX
      real(dp), intent(in) :: DY(*) ! second vector
      integer,  intent(in) :: INCY  ! increment along DY
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
  
  integer,         intent(in)  :: input(:)
  type(IntVector), intent(out) :: output
  
  output%contents = input
end subroutine

pure subroutine assign_RealVector_reals(output,input)
  implicit none
  
  real(dp),         intent(in)  :: input(:)
  type(RealVector), intent(out) :: output
  
  output%contents = input
end subroutine

pure subroutine assign_ComplexVector_complexes(output,input)
  implicit none
  
  complex(dp),         intent(in)  :: input(:)
  type(ComplexVector), intent(out) :: output
  
  output%contents = input
end subroutine

pure subroutine assign_IntMatrix_integers(output,input)
  implicit none
  
  integer,         intent(in)  :: input(:,:)
  type(IntMatrix), intent(out) :: output
  
  output%contents = input
end subroutine

pure subroutine assign_RealMatrix_reals(output,input)
  implicit none
  
  real(dp),         intent(in)  :: input(:,:)
  type(RealMatrix), intent(out) :: output
  
  output%contents = input
end subroutine

pure subroutine assign_ComplexMatrix_complexes(output,input)
  implicit none
  
  complex(dp),         intent(in)  :: input(:,:)
  type(ComplexMatrix), intent(out) :: output
  
  output%contents = input
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
  
  allocate(output%contents(n), stat=ialloc); call err(ialloc)
  output%contents = 0
end function

! Makes a mxn matrix full of zeroes.
function zeroes_IntMatrix(m,n) result(output)
  implicit none
  
  integer, intent(in) :: m
  integer, intent(in) :: n
  type(IntMatrix)     :: output
  
  integer :: ialloc
  
  allocate(output%contents(m,n), stat=ialloc); call err(ialloc)
  output%contents = 0
end function

! Makes an nxn identity matrix.
function identity_IntMatrix(n) result(output)
  implicit none
  
  integer, intent(in) :: n
  type(IntMatrix)     :: output
  
  integer :: i
  
  output = zeroes(n,n)
  do i=1,n
    output%contents(i,i) = 1
  enddo
end function

! Conversion to fundamental types.
pure function int_IntVector(input) result(output)
  implicit none
  
  type(IntVector), intent(in) :: input
  integer, allocatable        :: output(:)
  
  output = input%contents
end function

pure function int_IntMatrix(input) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: input
  integer, allocatable        :: output(:,:)
  
  output = input%contents
end function

pure function dble_RealVector(input) result(output)
  implicit none
  
  type(RealVector), intent(in) :: input
  real(dp), allocatable        :: output(:)
  
  output = input%contents
end function

pure function dble_RealMatrix(input) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: input
  real(dp), allocatable        :: output(:,:)
  
  output = input%contents
end function

pure function cmplx_ComplexVector(input) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: input
  complex(dp), allocatable        :: output(:)
  
  output = input%contents
end function

pure function cmplx_ComplexMatrix(input) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: input
  complex(dp), allocatable        :: output(:,:)
  
  output = input%contents
end function

! Real part.
pure function real_ComplexVector(input) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: input
  type(RealVector)                :: output
  
  output%contents = real(input%contents)
end function

pure function real_ComplexMatrix(input) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: input
  type(RealMatrix)                :: output
  
  output%contents = real(input%contents)
end function

! Imaginary part.
pure function imag_ComplexVector(input) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: input
  type(RealVector)                :: output
  
  output%contents = imag(input%contents)
end function

pure function imag_ComplexMatrix(input) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: input
  type(RealMatrix)                :: output
  
  output%contents = imag(input%contents)
end function

! size().
pure function size_IntVector(input) result(output)
  implicit none
  
  type(IntVector), intent(in) :: input
  integer                     :: output
  
  output = size(input%contents)
end function

pure function size_RealVector(input) result(output)
  implicit none
  
  type(RealVector), intent(in) :: input
  integer                      :: output
  
  output = size(input%contents)
end function

pure function size_ComplexVector(input) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: input
  integer                         :: output
  
  output = size(input%contents)
end function

pure function size_IntMatrix(input,dim) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: input
  integer,         intent(in) :: dim
  integer                     :: output
  
  output = size(input%contents, dim)
end function

pure function size_RealMatrix(input,dim) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: input
  integer,          intent(in) :: dim
  integer                      :: output
  
  output = size(input%contents, dim)
end function

pure function size_ComplexMatrix(input,dim) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: input
  integer,             intent(in) :: dim
  integer                         :: output
  
  output = size(input%contents, dim)
end function

! Trace
pure function trace_IntMatrix(input) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: input
  integer                     :: output
  
  integer :: i
  
  output = 0
  do i=1,size(input,1)
    output = output + input%contents(i,i)
  enddo
end function

pure function trace_RealMatrix(input) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: input
  real(dp)                     :: output
  
  integer :: i
  
  output = 0
  do i=1,size(input,1)
    output = output + input%contents(i,i)
  enddo
end function

pure function trace_ComplexMatrix(input) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: input
  complex(dp)                     :: output
  
  integer :: i
  
  output = 0
  do i=1,size(input,1)
    output = output + input%contents(i,i)
  enddo
end function

! Equality.
elemental function equality_IntVector_IntVector(a,b) result(output)
  implicit none
  
  type(IntVector), intent(in) :: a
  type(IntVector), intent(in) :: b
  logical                     :: output
  
  output = all(a%contents==b%contents)
end function

elemental function equality_IntMatrix_IntMatrix(a,b) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: a
  type(IntMatrix), intent(in) :: b
  logical                     :: output
  
  output = all(a%contents==b%contents)
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
pure function add_IntVector_IntVector(a,b) result(output)
  implicit none
  
  type(IntVector), intent(in) :: a
  type(IntVector), intent(in) :: b
  type(IntVector)             :: output
  
  output = a%contents + b%contents
end function

pure function add_IntVector_RealVector(a,b) result(output)
  implicit none
  
  type(IntVector),  intent(in) :: a
  type(RealVector), intent(in) :: b
  type(RealVector)             :: output
  
  output = a%contents + b%contents
end function

pure function add_RealVector_IntVector(a,b) result(output)
  implicit none
  
  type(RealVector), intent(in) :: a
  type(IntVector),  intent(in) :: b
  type(RealVector)             :: output
  
  output = a%contents + b%contents
end function

pure function add_RealVector_RealVector(a,b) result(output)
  implicit none
  
  type(RealVector), intent(in) :: a
  type(RealVector), intent(in) :: b
  type(RealVector)             :: output
  
  output = a%contents + b%contents
end function

pure function add_RealVector_ComplexVector(a,b) result(output)
  implicit none
  
  type(RealVector),    intent(in) :: a
  type(ComplexVector), intent(in) :: b
  type(ComplexVector)             :: output
  
  output = a%contents + b%contents
end function

pure function add_ComplexVector_RealVector(a,b) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: a
  type(RealVector),    intent(in) :: b
  type(ComplexVector)             :: output
  
  output = a%contents + b%contents
end function

pure function add_ComplexVector_ComplexVector(a,b) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: a
  type(ComplexVector), intent(in) :: b
  type(ComplexVector)             :: output
  
  output = a%contents + b%contents
end function

pure function add_IntMatrix_IntMatrix(a,b) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: a
  type(IntMatrix), intent(in) :: b
  type(IntMatrix)             :: output
  
  output = a%contents + b%contents
end function

pure function add_IntMatrix_RealMatrix(a,b) result(output)
  implicit none
  
  type(IntMatrix),  intent(in) :: a
  type(RealMatrix), intent(in) :: b
  type(RealMatrix)             :: output
  
  output = a%contents + b%contents
end function

pure function add_RealMatrix_IntMatrix(a,b) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: a
  type(IntMatrix),  intent(in) :: b
  type(RealMatrix)             :: output
  
  output = a%contents + b%contents
end function

pure function add_RealMatrix_RealMatrix(a,b) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: a
  type(RealMatrix), intent(in) :: b
  type(RealMatrix)             :: output
  
  output = a%contents + b%contents
end function

pure function add_RealMatrix_ComplexMatrix(a,b) result(output)
  implicit none
  
  type(RealMatrix),    intent(in) :: a
  type(ComplexMatrix), intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = a%contents + b%contents
end function

pure function add_ComplexMatrix_RealMatrix(a,b) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: a
  type(RealMatrix),    intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = a%contents + b%contents
end function

pure function add_ComplexMatrix_ComplexMatrix(a,b) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: a
  type(ComplexMatrix), intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = a%contents + b%contents
end function

! Negative.
pure function negative_IntVector(input) result(output)
  implicit none
  
  type(IntVector), intent(in) :: input
  type(IntVector)             :: output
  
  output = -input%contents
end function

pure function negative_RealVector(input) result(output)
  implicit none
  
  type(RealVector), intent(in) :: input
  type(RealVector)             :: output
  
  output = -input%contents
end function

pure function negative_ComplexVector(input) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: input
  type(ComplexVector)             :: output
  
  output = -input%contents
end function

pure function negative_IntMatrix(input) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: input
  type(IntMatrix)             :: output
  
  output = -input%contents
end function

pure function negative_RealMatrix(input) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: input
  type(RealMatrix)             :: output
  
  output = -input%contents
end function

pure function negative_ComplexMatrix(input) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: input
  type(ComplexMatrix)             :: output
  
  output = -input%contents
end function

! Subtraction.
pure function subtract_IntVector_IntVector(a,b) result(output)
  implicit none
  
  type(IntVector), intent(in) :: a
  type(IntVector), intent(in) :: b
  type(IntVector)             :: output
  
  output = a%contents - b%contents
end function

pure function subtract_IntVector_RealVector(a,b) result(output)
  implicit none
  
  type(IntVector),  intent(in) :: a
  type(RealVector), intent(in) :: b
  type(RealVector)             :: output
  
  output = a%contents - b%contents
end function

pure function subtract_RealVector_IntVector(a,b) result(output)
  implicit none
  
  type(RealVector), intent(in) :: a
  type(IntVector),  intent(in) :: b
  type(RealVector)             :: output
  
  output = a%contents - b%contents
end function

pure function subtract_RealVector_RealVector(a,b) result(output)
  implicit none
  
  type(RealVector), intent(in) :: a
  type(RealVector), intent(in) :: b
  type(RealVector)             :: output
  
  output = a%contents - b%contents
end function

pure function subtract_RealVector_ComplexVector(a,b) result(output)
  implicit none
  
  type(RealVector),    intent(in) :: a
  type(ComplexVector), intent(in) :: b
  type(ComplexVector)             :: output
  
  output = a%contents - b%contents
end function

pure function subtract_ComplexVector_RealVector(a,b) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: a
  type(RealVector),    intent(in) :: b
  type(ComplexVector)             :: output
  
  output = a%contents - b%contents
end function

pure function subtract_ComplexVector_ComplexVector(a,b) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: a
  type(ComplexVector), intent(in) :: b
  type(ComplexVector)             :: output
  
  output = a%contents - b%contents
end function

pure function subtract_IntMatrix_IntMatrix(a,b) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: a
  type(IntMatrix), intent(in) :: b
  type(IntMatrix)             :: output
  
  output = a%contents - b%contents
end function

pure function subtract_IntMatrix_RealMatrix(a,b) result(output)
  implicit none
  
  type(IntMatrix),  intent(in) :: a
  type(RealMatrix), intent(in) :: b
  type(RealMatrix)             :: output
  
  output = a%contents - b%contents
end function

pure function subtract_RealMatrix_IntMatrix(a,b) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: a
  type(IntMatrix),  intent(in) :: b
  type(RealMatrix)             :: output
  
  output = a%contents - b%contents
end function

pure function subtract_RealMatrix_RealMatrix(a,b) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: a
  type(RealMatrix), intent(in) :: b
  type(RealMatrix)             :: output
  
  output = a%contents - b%contents
end function

pure function subtract_RealMatrix_ComplexMatrix(a,b) result(output)
  implicit none
  
  type(RealMatrix),    intent(in) :: a
  type(ComplexMatrix), intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = a%contents - b%contents
end function

pure function subtract_ComplexMatrix_RealMatrix(a,b) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: a
  type(RealMatrix),    intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = a%contents - b%contents
end function

pure function subtract_ComplexMatrix_ComplexMatrix(a,b) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: a
  type(ComplexMatrix), intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = a%contents - b%contents
end function

! Multiplication by a scalar.
pure function multiply_IntVector_integer(a,b) result(output)
  implicit none
  
  type(IntVector), intent(in) :: a
  integer,         intent(in) :: b
  type(IntVector)             :: output
  
  output = a%contents*b
end function

pure function multiply_integer_IntVector(a,b) result(output)
  implicit none
  
  integer,         intent(in) :: a
  type(IntVector), intent(in) :: b
  type(IntVector)             :: output
  
  output = a*b%contents
end function

pure function multiply_IntVector_real(a,b) result(output)
  implicit none
  
  type(IntVector), intent(in) :: a
  real(dp),        intent(in) :: b
  type(RealVector)            :: output
  
  output = a%contents*b
end function

pure function multiply_real_IntVector(a,b) result(output)
  implicit none
  
  real(dp),        intent(in) :: a
  type(IntVector), intent(in) :: b
  type(RealVector)            :: output
  
  output = a*b%contents
end function

pure function multiply_RealVector_integer(a,b) result(output)
  implicit none
  
  type(RealVector), intent(in) :: a
  integer,          intent(in) :: b
  type(RealVector)             :: output
  
  output = a%contents*b
end function

pure function multiply_integer_RealVector(a,b) result(output)
  implicit none
  
  integer,          intent(in) :: a
  type(RealVector), intent(in) :: b
  type(RealVector)             :: output
  
  output = a*b%contents
end function

pure function multiply_RealVector_real(a,b) result(output)
  implicit none
  
  type(RealVector), intent(in) :: a
  real(dp),         intent(in) :: b
  type(RealVector)             :: output
  
  output = a%contents*b
end function

pure function multiply_real_RealVector(a,b) result(output)
  implicit none
  
  real(dp),         intent(in) :: a
  type(RealVector), intent(in) :: b
  type(RealVector)             :: output
  
  output = a*b%contents
end function

pure function multiply_RealVector_complex(a,b) result(output)
  implicit none
  
  type(RealVector), intent(in) :: a
  complex(dp),      intent(in) :: b
  type(ComplexVector)          :: output
  
  output = a%contents*b
end function

pure function multiply_complex_RealVector(a,b) result(output)
  implicit none
  
  complex(dp),      intent(in) :: a
  type(RealVector), intent(in) :: b
  type(ComplexVector)          :: output
  
  output = a*b%contents
end function

pure function multiply_ComplexVector_real(a,b) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: a
  real(dp),            intent(in) :: b
  type(ComplexVector)          :: output
  
  output = a%contents*b
end function

pure function multiply_real_ComplexVector(a,b) result(output)
  implicit none
  
  real(dp),            intent(in) :: a
  type(ComplexVector), intent(in) :: b
  type(ComplexVector)          :: output
  
  output = a*b%contents
end function

pure function multiply_ComplexVector_complex(a,b) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: a
  complex(dp),         intent(in) :: b
  type(ComplexVector)             :: output
  
  output = a%contents*b
end function

pure function multiply_complex_ComplexVector(a,b) result(output)
  implicit none
  
  complex(dp),         intent(in) :: a
  type(ComplexVector), intent(in) :: b
  type(ComplexVector)          :: output
  
  output = a*b%contents
end function

pure function multiply_IntMatrix_integer(a,b) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: a
  integer,         intent(in) :: b
  type(IntMatrix)             :: output
  
  output = a%contents*b
end function

pure function multiply_integer_IntMatrix(a,b) result(output)
  implicit none
  
  integer,         intent(in) :: a
  type(IntMatrix), intent(in) :: b
  type(IntMatrix)             :: output
  
  output = a*b%contents
end function

pure function multiply_IntMatrix_real(a,b) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: a
  real(dp),        intent(in) :: b
  type(RealMatrix)            :: output
  
  output = a%contents*b
end function

pure function multiply_real_IntMatrix(a,b) result(output)
  implicit none
  
  real(dp),        intent(in) :: a
  type(IntMatrix), intent(in) :: b
  type(RealMatrix)            :: output
  
  output = a*b%contents
end function

pure function multiply_RealMatrix_integer(a,b) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: a
  integer,          intent(in) :: b
  type(RealMatrix)             :: output
  
  output = a%contents*b
end function

pure function multiply_integer_RealMatrix(a,b) result(output)
  implicit none
  
  integer,          intent(in) :: a
  type(RealMatrix), intent(in) :: b
  type(RealMatrix)             :: output
  
  output = a*b%contents
end function

pure function multiply_RealMatrix_real(a,b) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: a
  real(dp),         intent(in) :: b
  type(RealMatrix)             :: output
  
  output = a%contents*b
end function

pure function multiply_real_RealMatrix(a,b) result(output)
  implicit none
  
  real(dp),         intent(in) :: a
  type(RealMatrix), intent(in) :: b
  type(RealMatrix)             :: output
  
  output = a*b%contents
end function

pure function multiply_RealMatrix_complex(a,b) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: a
  complex(dp),      intent(in) :: b
  type(ComplexMatrix)          :: output
  
  output = a%contents*b
end function

pure function multiply_complex_RealMatrix(a,b) result(output)
  implicit none
  
  complex(dp),      intent(in) :: a
  type(RealMatrix), intent(in) :: b
  type(ComplexMatrix)          :: output
  
  output = a*b%contents
end function

pure function multiply_ComplexMatrix_real(a,b) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: a
  real(dp),            intent(in) :: b
  type(ComplexMatrix)          :: output
  
  output = a%contents*b
end function

pure function multiply_real_ComplexMatrix(a,b) result(output)
  implicit none
  
  real(dp),            intent(in) :: a
  type(ComplexMatrix), intent(in) :: b
  type(ComplexMatrix)          :: output
  
  output = a*b%contents
end function

pure function multiply_ComplexMatrix_complex(a,b) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: a
  complex(dp),         intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = a%contents*b
end function

pure function multiply_complex_ComplexMatrix(a,b) result(output)
  implicit none
  
  complex(dp),         intent(in) :: a
  type(ComplexMatrix), intent(in) :: b
  type(ComplexMatrix)          :: output
  
  output = a*b%contents
end function

! Dot products and matrix multiplication.
pure function dot_IntVector_IntVector(a,b) result(output)
  implicit none
  
  type(IntVector), intent(in) :: a
  type(IntVector), intent(in) :: b
  integer                     :: output
  
  output = dot_product(a%contents, b%contents)
end function

pure function dot_IntVector_RealVector(a,b) result(output)
  implicit none
  
  type(IntVector),  intent(in) :: a
  type(RealVector), intent(in) :: b
  real(dp)                     :: output
  
  output = dot_product(a%contents, b%contents)
end function

pure function dot_RealVector_IntVector(a,b) result(output)
  implicit none
  
  type(RealVector), intent(in) :: a
  type(IntVector),  intent(in) :: b
  real(dp)                     :: output
  
  output = dot_product(a%contents, b%contents)
end function

pure function dot_RealVector_RealVector(a,b) result(output)
  implicit none
  
  type(RealVector), intent(in) :: a
  type(RealVector), intent(in) :: b
  real(dp)                     :: output
  
  output = dot_product(a%contents, b%contents)
end function

pure function dot_RealVector_ComplexVector(a,b) result(output)
  implicit none
  
  type(RealVector),    intent(in) :: a
  type(ComplexVector), intent(in) :: b
  complex(dp)                     :: output
  
  output = dot_product(a%contents, b%contents)
end function

pure function dot_ComplexVector_RealVector(a,b) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: a
  type(RealVector),    intent(in) :: b
  complex(dp)                     :: output
  
  output = dot_product(a%contents, b%contents)
end function

pure function dot_ComplexVector_ComplexVector(a,b) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: a
  type(ComplexVector),    intent(in) :: b
  complex(dp)                     :: output
  
  output = dot_product(a%contents, b%contents)
end function

pure function dot_IntVector_IntMatrix(a,b) result(output)
  implicit none
  
  type(IntVector), intent(in) :: a
  type(IntMatrix), intent(in) :: b
  type(IntVector)             :: output
  
  output = matmul(a%contents, b%contents)
end function

pure function dot_IntVector_RealMatrix(a,b) result(output)
  implicit none
  
  type(IntVector),  intent(in) :: a
  type(RealMatrix), intent(in) :: b
  type(RealVector)             :: output
  
  output = matmul(a%contents, b%contents)
end function

pure function dot_RealVector_IntMatrix(a,b) result(output)
  implicit none
  
  type(RealVector), intent(in) :: a
  type(IntMatrix),  intent(in) :: b
  type(RealVector)             :: output
  
  output = matmul(a%contents, b%contents)
end function

pure function dot_RealVector_RealMatrix(a,b) result(output)
  implicit none
  
  type(RealVector), intent(in) :: a
  type(RealMatrix), intent(in) :: b
  type(RealVector)             :: output
  
  output = matmul(a%contents, b%contents)
end function

pure function dot_RealVector_ComplexMatrix(a,b) result(output)
  implicit none
  
  type(RealVector),    intent(in) :: a
  type(ComplexMatrix), intent(in) :: b
  type(ComplexVector)             :: output
  
  output = matmul(a%contents, b%contents)
end function

pure function dot_ComplexVector_RealMatrix(a,b) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: a
  type(RealMatrix),    intent(in) :: b
  type(ComplexVector)             :: output
  
  output = matmul(a%contents, b%contents)
end function

pure function dot_ComplexVector_ComplexMatrix(a,b) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: a
  type(ComplexMatrix), intent(in) :: b
  type(ComplexVector)             :: output
  
  output = matmul(a%contents, b%contents)
end function

pure function dot_IntMatrix_IntVector(a,b) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: a
  type(IntVector), intent(in) :: b
  type(IntVector)             :: output
  
  output = matmul(a%contents, b%contents)
end function

pure function dot_IntMatrix_RealVector(a,b) result(output)
  implicit none
  
  type(IntMatrix),  intent(in) :: a
  type(RealVector), intent(in) :: b
  type(RealVector)             :: output
  
  output = matmul(a%contents, b%contents)
end function

pure function dot_RealMatrix_IntVector(a,b) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: a
  type(IntVector),  intent(in) :: b
  type(RealVector)             :: output
  
  output = matmul(a%contents, b%contents)
end function

pure function dot_RealMatrix_RealVector(a,b) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: a
  type(RealVector), intent(in) :: b
  type(RealVector)             :: output
  
  output = matmul(a%contents, b%contents)
end function

pure function dot_RealMatrix_ComplexVector(a,b) result(output)
  implicit none
  
  type(RealMatrix),    intent(in) :: a
  type(ComplexVector), intent(in) :: b
  type(ComplexVector)             :: output
  
  output = matmul(a%contents, b%contents)
end function

pure function dot_ComplexMatrix_RealVector(a,b) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: a
  type(RealVector),    intent(in) :: b
  type(ComplexVector)             :: output
  
  output = matmul(a%contents, b%contents)
end function

pure function dot_ComplexMatrix_ComplexVector(a,b) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: a
  type(ComplexVector), intent(in) :: b
  type(ComplexVector)             :: output
  
  output = matmul(a%contents, b%contents)
end function

pure function dot_IntMatrix_IntMatrix(a,b) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: a
  type(IntMatrix), intent(in) :: b
  type(IntMatrix)             :: output
  
  output = matmul(a%contents, b%contents)
end function

pure function dot_IntMatrix_RealMatrix(a,b) result(output)
  implicit none
  
  type(IntMatrix),  intent(in) :: a
  type(RealMatrix), intent(in) :: b
  type(RealMatrix)             :: output
  
  output = matmul(a%contents, b%contents)
end function

pure function dot_RealMatrix_IntMatrix(a,b) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: a
  type(IntMatrix),  intent(in) :: b
  type(RealMatrix)             :: output
  
  output = matmul(a%contents, b%contents)
end function

pure function dot_RealMatrix_RealMatrix(a,b) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: a
  type(RealMatrix), intent(in) :: b
  type(RealMatrix)             :: output
  
  output = matmul(a%contents, b%contents)
end function

pure function dot_RealMatrix_ComplexMatrix(a,b) result(output)
  implicit none
  
  type(RealMatrix),    intent(in) :: a
  type(ComplexMatrix), intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = matmul(a%contents, b%contents)
end function

pure function dot_ComplexMatrix_RealMatrix(a,b) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: a
  type(RealMatrix),    intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = matmul(a%contents, b%contents)
end function

pure function dot_ComplexMatrix_ComplexMatrix(a,b) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: a
  type(ComplexMatrix), intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = matmul(a%contents, b%contents)
end function

! Division by scalar.
pure function divide_IntVector_integer(a,b) result(output)
  implicit none
  
  type(IntVector), intent(in) :: a
  integer,         intent(in) :: b
  type(IntVector)             :: output
  
  output = a%contents/b
end function

pure function divide_IntVector_real(a,b) result(output)
  implicit none
  
  type(IntVector), intent(in) :: a
  real(dp),        intent(in) :: b
  type(RealVector)            :: output
  
  output = a%contents/b
end function

pure function divide_RealVector_integer(a,b) result(output)
  implicit none
  
  type(RealVector), intent(in) :: a
  integer,          intent(in) :: b
  type(RealVector)             :: output
  
  output = a%contents/b
end function

pure function divide_RealVector_real(a,b) result(output)
  implicit none
  
  type(RealVector), intent(in) :: a
  real(dp),         intent(in) :: b
  type(RealVector)             :: output
  
  output = a%contents/b
end function

pure function divide_RealVector_complex(a,b) result(output)
  implicit none
  
  type(RealVector),    intent(in) :: a
  complex(dp),         intent(in) :: b
  type(ComplexVector)             :: output
  
  output = a%contents/b
end function

pure function divide_ComplexVector_real(a,b) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: a
  real(dp),            intent(in) :: b
  type(ComplexVector)             :: output
  
  output = a%contents/b
end function

pure function divide_ComplexVector_complex(a,b) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: a
  complex(dp),         intent(in) :: b
  type(ComplexVector)             :: output
  
  output = a%contents/b
end function

pure function divide_IntMatrix_integer(a,b) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: a
  integer,         intent(in) :: b
  type(IntMatrix)             :: output
  
  output = a%contents/b
end function

pure function divide_IntMatrix_real(a,b) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: a
  real(dp),        intent(in) :: b
  type(RealMatrix)            :: output
  
  output = a%contents/b
end function

pure function divide_RealMatrix_integer(a,b) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: a
  integer,          intent(in) :: b
  type(RealMatrix)             :: output
  
  output = a%contents/b
end function

pure function divide_RealMatrix_real(a,b) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: a
  real(dp),         intent(in) :: b
  type(RealMatrix)             :: output
  
  output = a%contents/b
end function

pure function divide_RealMatrix_complex(a,b) result(output)
  implicit none
  
  type(RealMatrix),    intent(in) :: a
  complex(dp),         intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = a%contents/b
end function

pure function divide_ComplexMatrix_real(a,b) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: a
  real(dp),            intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = a%contents/b
end function

pure function divide_ComplexMatrix_complex(a,b) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: a
  complex(dp),         intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = a%contents/b
end function

! L2 norm.
pure function l2_norm_RealVector(input) result(output)
  implicit none
  
  type(RealVector), intent(in) :: input
  real(dp)                     :: output
  
  output = sqrt(input*input)
end function

! Outer product.
function outer_product_RealVector_RealVector(a,b) result(output)
  implicit none
  
  type(RealVector), intent(in) :: a
  type(RealVector), intent(in) :: b
  type(RealMatrix)             :: output
  
  integer               :: i,ialloc
  
  allocate(output%contents(size(a),size(b)), stat=ialloc); call err(ialloc)
  do i=1,size(b)
    output%contents(:,i) = a%contents*b%contents(i)
  enddo
end function

! Transpose.
pure function transpose_IntMatrix(input) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: input
  type(IntMatrix)             :: output
  
  output = transpose(input%contents)
end function

pure function transpose_RealMatrix(input) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: input
  type(RealMatrix)             :: output
  
  output = transpose(input%contents)
end function

! Hermitian conjugate.
pure function hermitian_ComplexMatrix(input) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: input
  type(ComplexMatrix)             :: output
  
  output = transpose(conjg(input%contents))
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
  
  output = determinant(input%contents)
end function

function determinant_RealMatrix(input) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: input
  real(dp)                     :: output
  
  output = determinant(input%contents)
end function

! calculates the inverse, B, of matrix A. Both are 3x3 matrices.
function invert_real(A) result(B)
  implicit none
  
  real(dp), intent(in)  :: A(:,:)
  real(dp)              :: B(3,3)
  
  real(dp) :: d      ! 1/det(A)
  real(dp) :: C(3,3) ! transpose(A)
  
  ! Check that A is a 3x3 matrix.
  if (size(A,1)/=3 .or. size(A,2)/=3) then
    call err()
  endif
  
  d = 1.0_dp/determinant(A)
  
  ! check for d=infinity or d=NaN
  if (abs(d)>huge(0.0_dp) .or. d<d) then
    call print_line('Error in invert: singular matrix.')
    call err()
  endif
  
  C = transpose(A)
  
  B(1,1) = (C(2,2)*C(3,3)-C(2,3)*C(3,2))*d
  B(1,2) = (C(2,3)*C(3,1)-C(2,1)*C(3,3))*d
  B(1,3) = (C(2,1)*C(3,2)-C(2,2)*C(3,1))*d
  B(2,1) = (C(3,2)*C(1,3)-C(3,3)*C(1,2))*d
  B(2,2) = (C(3,3)*C(1,1)-C(3,1)*C(1,3))*d
  B(2,3) = (C(3,1)*C(1,2)-C(3,2)*C(1,1))*d
  B(3,1) = (C(1,2)*C(2,3)-C(1,3)*C(2,2))*d
  B(3,2) = (C(1,3)*C(2,1)-C(1,1)*C(2,3))*d
  B(3,3) = (C(1,1)*C(2,2)-C(1,2)*C(2,1))*d
end function

function invert_RealMatrix(input) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: input
  type(RealMatrix)             :: output
  
  output = invert(input%contents)
end function

! Calculates B=invert(A)*|det(A)|
function invert_int_integer(A) result(B)
  implicit none
  
  integer, intent(in)  :: A(:,:)
  integer              :: B(3,3)
  
  integer :: C(3,3) ! transpose(A)
  integer :: d      ! 1 if det(A)>=0, -1 otherwise
  
  ! Check that A is a 3x3 matrix.
  if (size(A,1)/=3 .or. size(A,2)/=3) then
    call err()
  endif
  
  d = sign(1,determinant(A))
  
  C = transpose(A)
  
  B(1,1) = (C(2,2)*C(3,3)-C(2,3)*C(3,2))*d
  B(1,2) = (C(2,3)*C(3,1)-C(2,1)*C(3,3))*d
  B(1,3) = (C(2,1)*C(3,2)-C(2,2)*C(3,1))*d
  B(2,1) = (C(3,2)*C(1,3)-C(3,3)*C(1,2))*d
  B(2,2) = (C(3,3)*C(1,1)-C(3,1)*C(1,3))*d
  B(2,3) = (C(3,1)*C(1,2)-C(3,2)*C(1,1))*d
  B(3,1) = (C(1,2)*C(2,3)-C(1,3)*C(2,2))*d
  B(3,2) = (C(1,3)*C(2,1)-C(1,1)*C(2,3))*d
  B(3,3) = (C(1,1)*C(2,2)-C(1,2)*C(2,1))*d
end function

function invert_int_IntMatrix(input) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: input
  type(IntMatrix)             :: output
  
  output = invert_int(input%contents)
end function

! String concatenation functions.
function concatenate_character_IntVector(a,b) result(output)
  implicit none
  
  character(*),    intent(in) :: a
  type(IntVector), intent(in) :: b
  type(String)                :: output
  
  output = a//b%contents
end function

function concatenate_IntVector_character(a,b) result(output)
  implicit none
  
  type(IntVector), intent(in) :: a
  character(*),    intent(in) :: b
  type(String)                :: output
  
  output = a%contents//b
end function

function concatenate_character_RealVector(a,b) result(output)
  implicit none
  
  character(*),     intent(in) :: a
  type(RealVector), intent(in) :: b
  type(String)                 :: output
  
  output = a//b%contents
end function

function concatenate_RealVector_character(a,b) result(output)
  implicit none
  
  type(RealVector), intent(in) :: a
  character(*),     intent(in) :: b
  type(String)                 :: output
  
  output = a%contents//b
end function

function concatenate_String_IntVector(a,b) result(output)
  implicit none
  
  type(String),    intent(in) :: a
  type(IntVector), intent(in) :: b
  type(String)                :: output
  
  output = a//b%contents
end function

function concatenate_IntVector_String(a,b) result(output)
  implicit none
  
  type(IntVector), intent(in) :: a
  type(String),    intent(in) :: b
  type(String)                :: output
  
  output = a%contents//b
end function

function concatenate_String_RealVector(a,b) result(output)
  implicit none
  
  type(String),     intent(in) :: a
  type(RealVector), intent(in) :: b
  type(String)                 :: output
  
  output = a//b%contents
end function

function concatenate_RealVector_String(a,b) result(output)
  implicit none
  
  type(RealVector), intent(in) :: a
  type(String),     intent(in) :: b
  type(String)                 :: output
  
  output = a%contents//b
end function

! Print functions.
subroutine print_line_IntVector(input)
  implicit none
  
  type(IntVector), intent(in) :: input
  
  call print_line(input%contents)
end subroutine

subroutine print_line_IntVector_file(file_unit,input)
  implicit none
  
  integer,         intent(in) :: file_unit
  type(IntVector), intent(in) :: input
  
  call print_line(file_unit, input%contents)
end subroutine

subroutine print_line_RealVector(input)
  implicit none
  
  type(RealVector), intent(in) :: input
  
  call print_line(input%contents)
end subroutine

subroutine print_line_RealVector_file(file_unit,input)
  implicit none
  
  integer,          intent(in) :: file_unit
  type(RealVector), intent(in) :: input
  
  call print_line(file_unit, input%contents)
end subroutine

subroutine print_line_ComplexVector(input)
  implicit none
  
  type(ComplexVector), intent(in) :: input
  
  call print_line(input%contents)
end subroutine

subroutine print_line_ComplexVector_file(file_unit,input)
  implicit none
  
  integer,          intent(in) :: file_unit
  type(ComplexVector), intent(in) :: input
  
  call print_line(file_unit, input%contents)
end subroutine

subroutine print_line_IntMatrix(input)
  implicit none
  
  type(IntMatrix), intent(in) :: input
  
  integer :: i
  
  do i=1,size(input,1)
    call print_line(input%contents(i,:))
  enddo
end subroutine

subroutine print_line_IntMatrix_file(file_unit,input)
  implicit none
  
  integer,         intent(in) :: file_unit
  type(IntMatrix), intent(in) :: input
  
  integer :: i
  
  do i=1,size(input,1)
    call print_line(file_unit, input%contents(i,:))
  enddo
end subroutine

subroutine print_line_RealMatrix(input)
  implicit none
  
  type(RealMatrix), intent(in) :: input
  
  integer :: i
  
  do i=1,size(input,1)
    call print_line(input%contents(i,:))
  enddo
end subroutine

subroutine print_line_RealMatrix_file(file_unit,input)
  implicit none
  
  integer,          intent(in) :: file_unit
  type(RealMatrix), intent(in) :: input
  
  integer :: i
  
  do i=1,size(input,1)
    call print_line(file_unit, input%contents(i,:))
  enddo
end subroutine

subroutine print_line_ComplexMatrix(input)
  implicit none
  
  type(ComplexMatrix), intent(in) :: input
  
  integer :: i
  
  do i=1,size(input,1)
    call print_line(input%contents(i,:))
  enddo
end subroutine

subroutine print_line_ComplexMatrix_file(file_unit,input)
  implicit none
  
  integer,             intent(in) :: file_unit
  type(ComplexMatrix), intent(in) :: input
  
  integer :: i
  
  do i=1,size(input,1)
    call print_line(file_unit, input%contents(i,:))
  enddo
end subroutine

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
  elseif (m<=n) then
    call print_line('Error in least-squares optimisation: there are more &
       &variables to be fit than input data-points.')
    call err()
  endif
  
  ! Copy a and b so that they are not changed by dgels.
  a2 = dble(a)
  x = b
  
  ! Calculate optimal workspace size.
  allocate(work(2*m*n), stat=ialloc); call err(ialloc)
  call dgels('N',m,n,1,a2(1,1),m,x%contents(1),m,work(1),-1,info)
  if (info/=0) then
    call print_line('dgels failed, info= '//info)
    call err()
  endif
  lwork = nint(work(1))
  deallocate(work, stat=ialloc); call err(ialloc)
  allocate(work(lwork), stat=ialloc); call err(ialloc)
  
  ! Run linear least-squares optimisation.
  call dgels('N',m,n,1,a2(1,1),m,x%contents(1),m,work(1),lwork,info)
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
