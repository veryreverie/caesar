! ======================================================================
! Assorted linear algebra / vector algebra subroutines.
! Includes interfaces for BLAS and LAPACK routines.
! ======================================================================
module linear_algebra_submodule
  use precision_module
  use io_module
  implicit none
  
  private
  
  public :: IntVector
  public :: RealVector
  public :: ComplexVector
  public :: IntMatrix
  public :: RealMatrix
  public :: ComplexMatrix
  public :: assignment(=)
  public :: vec
  public :: mat
  public :: zeroes
  public :: make_identity_matrix
  public :: int
  public :: intvec
  public :: intmat
  public :: nint
  public :: dble
  public :: dblevec
  public :: dblemat
  public :: cmplx
  public :: cmplxvec
  public :: cmplxmat
  public :: row_matrix
  public :: real
  public :: aimag
  public :: abs
  public :: conjg
  public :: size
  public :: trace
  public :: commutator
  public :: matrices_commute
  public :: operator(==)
  public :: operator(/=)
  public :: operator(+)
  public :: operator(-)
  public :: operator(*)
  public :: operator(/)
  public :: sum
  public :: l2_norm
  public :: outer_product
  public :: transpose
  public :: hermitian
  public :: determinant
  public :: invert
  public :: linear_least_squares
  
  ! --------------------------------------------------
  ! Vector and Matrix classes
  ! --------------------------------------------------
  ! The actual Vector an Matrix classes extend from Vectorable and Matrixable
  !    classes so that fewer overloads are needed.
  ! Any type which extends ***able can be converted to *** by calling
  !    %to_***().
  type, abstract, extends(Stringable) :: ComplexVectorable
  contains
    procedure(to_ComplexVector_ComplexVectorable), deferred :: &
       & to_ComplexVector
  end type
  
  type, abstract, extends(ComplexVectorable) :: RealVectorable
  contains
    procedure(to_RealVector_RealVectorable), deferred :: to_RealVector
  end type
  
  type, abstract, extends(RealVectorable) :: IntVectorable
  contains
    procedure(to_IntVector_IntVectorable), deferred :: to_IntVector
  end type
  
  type, extends(IntVectorable) :: IntVector
    integer, allocatable, private :: contents_(:)
  contains
    procedure, public :: to_IntVector     => to_IntVector_IntVector
    procedure, public :: to_RealVector    => to_RealVector_IntVector
    procedure, public :: to_ComplexVector => to_ComplexVector_IntVector
    
    procedure, public :: read  => read_IntVector
    procedure, public :: write => write_IntVector
  end type
  
  type, extends(RealVectorable) :: RealVector
    real(dp), allocatable, private :: contents_(:)
  contains
    procedure, public :: to_RealVector    => to_RealVector_RealVector
    procedure, public :: to_ComplexVector => to_ComplexVector_RealVector
    
    procedure, public :: read  => read_RealVector
    procedure, public :: write => write_RealVector
  end type
  
  type, extends(ComplexVectorable) :: ComplexVector
    complex(dp), allocatable, private :: contents_(:)
  contains
    procedure, public :: to_ComplexVector => to_ComplexVector_ComplexVector
    
    procedure, public :: read  => read_ComplexVector
    procedure, public :: write => write_ComplexVector
  end type
  
  abstract interface
    function to_ComplexVector_ComplexVectorable(this) result(output)
      import ComplexVector
      import ComplexVectorable
      
      class(ComplexVectorable), intent(in) :: this
      type(ComplexVector)                  :: output
    end function
    
    function to_RealVector_RealVectorable(this) result(output)
      import RealVector
      import RealVectorable
      
      class(RealVectorable), intent(in) :: this
      type(RealVector)                  :: output
    end function
    
    function to_IntVector_IntVectorable(this) result(output)
      import IntVector
      import IntVectorable
      
      class(IntVectorable), intent(in) :: this
      type(IntVector)                  :: output
    end function
  end interface
  
  type, abstract, extends(Stringsable) :: ComplexMatrixable
  contains
    procedure(to_ComplexMatrix_ComplexMatrixable), deferred :: &
       & to_ComplexMatrix
  end type
  
  type, abstract, extends(ComplexMatrixable) :: RealMatrixable
  contains
    procedure(to_RealMatrix_RealMatrixable), deferred :: to_RealMatrix
  end type
  
  type, abstract, extends(RealMatrixable) :: IntMatrixable
  contains
    procedure(to_IntMatrix_IntMatrixable), deferred :: to_IntMatrix
  end type
  
  type, extends(IntMatrixable) :: IntMatrix
    integer, allocatable, private :: contents_(:,:)
  contains
    procedure, public :: to_IntMatrix     => to_IntMatrix_IntMatrix
    procedure, public :: to_RealMatrix    => to_RealMatrix_IntMatrix
    procedure, public :: to_ComplexMatrix => to_ComplexMatrix_IntMatrix
    
    procedure, public :: read  => read_IntMatrix
    procedure, public :: write => write_IntMatrix
  end type
  
  type, extends(RealMatrixable) :: RealMatrix
    real(dp), allocatable, private :: contents_(:,:)
  contains
    procedure, public :: to_RealMatrix    => to_RealMatrix_RealMatrix
    procedure, public :: to_ComplexMatrix => to_ComplexMatrix_RealMatrix
    
    procedure, public :: read  => read_RealMatrix
    procedure, public :: write => write_RealMatrix
  end type
  
  type, extends(ComplexMatrixable) :: ComplexMatrix
    complex(dp), allocatable, private :: contents_(:,:)
  contains
    procedure, public :: to_ComplexMatrix => to_ComplexMatrix_ComplexMatrix
    
    procedure, public :: read  => read_ComplexMatrix
    procedure, public :: write => write_ComplexMatrix
  end type
  
  abstract interface
    function to_ComplexMatrix_ComplexMatrixable(this) result(output)
      import ComplexMatrix
      import ComplexMatrixable
      
      class(ComplexMatrixable), intent(in) :: this
      type(ComplexMatrix)                  :: output
    end function
    
    function to_RealMatrix_RealMatrixable(this) result(output)
      import RealMatrix
      import RealMatrixable
      
      class(RealMatrixable), intent(in) :: this
      type(RealMatrix)                  :: output
    end function
    
    function to_IntMatrix_IntMatrixable(this) result(output)
      import IntMatrix
      import IntMatrixable
      
      class(IntMatrixable), intent(in) :: this
      type(IntMatrix)                  :: output
    end function
  end interface
  
  ! --------------------------------------------------
  ! Assignment and retrieval of contents.
  ! --------------------------------------------------
  interface assignment(=)
    module procedure assign_IntVector_integers
    module procedure assign_RealVector_reals
    module procedure assign_ComplexVector_complexes
    module procedure assign_IntMatrix_integers
    module procedure assign_RealMatrix_reals
    module procedure assign_ComplexMatrix_complexes
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
  
  ! --------------------------------------------------
  ! Conversion between types.
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
  
  interface intvec
    module procedure intvec_IntVectorable
  end interface
  
  interface intmat
    module procedure intmat_IntMatrixable
  end interface
  
  interface dblevec
    module procedure dblevec_RealVectorable
  end interface
  
  interface dblemat
    module procedure dblemat_RealMatrixable
  end interface
  
  interface cmplxvec
    module procedure cmplxvec_ComplexVectorable
    module procedure cmplxvec_RealVectorables
  end interface
  
  interface cmplxmat
    module procedure cmplxmat_ComplexMatrixable
    module procedure cmplxmat_RealMatrixables
  end interface
  
  ! --------------------------------------------------
  ! Comparison and arithmetic.
  ! --------------------------------------------------
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
  
  ! --------------------------------------------------
  ! Other operations.
  ! --------------------------------------------------
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
  
  interface commutator
    module procedure commutator_IntMatrix_IntMatrix
  end interface
  
  interface matrices_commute
    module procedure matrices_commute_IntMatrix_IntMatrix
  end interface
  
  interface row_matrix
    module procedure row_matrix_IntVectors
    module procedure row_matrix_RealVectors
    module procedure row_matrix_ComplexVectors
  end interface
  
  interface sum
    module procedure sum_IntVectors
    module procedure sum_RealVectors
    module procedure sum_ComplexVectors
    module procedure sum_IntMatrices
    module procedure sum_RealMatrices
    module procedure sum_ComplexMatrices
  end interface
  
  interface l2_norm
    module procedure l2_norm_RealVector
    module procedure l2_norm_ComplexVector
  end interface
  
  interface outer_product
    module procedure outer_product_RealVector_RealVector
    module procedure outer_product_ComplexVector_ComplexVector
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
  
  ! --------------------------------------------------
  ! BLAS / LAPACK interface.
  ! --------------------------------------------------
  interface
    ! Copies a real vector. Equivalent to dy = dx.
    subroutine dcopy(n,dx,incx,dy,incy)
      import :: dp
      implicit none
      
      integer,  intent(in)  :: n     ! Length of vectors
      real(dp), intent(in)  :: dx(*) ! Input vector
      integer,  intent(in)  :: incx  ! Increment along dx
      real(dp), intent(out) :: dy(*) ! Output vector
      integer,  intent(in)  :: incy  ! Increment along dy
    end subroutine
    
    ! Copies complex vector. Equivalent to zy = zx.
    subroutine zcopy(n,zx,incx,zy,incy)
      import :: dp
      implicit none
      
      integer,     intent(in)  :: n     ! Length of vectors
      complex(dp), intent(in)  :: zx(*) ! Input vector
      integer,     intent(in)  :: incx  ! Increment along zx
      complex(dp), intent(out) :: zy(*) ! Output vector
      integer,     intent(in)  :: incy  ! Increment along zy
    end subroutine
    
    ! Real dot product. Returns dx.dy.
    function ddot(n,dx,incx,dy,incy) result(output)
      import :: dp
      implicit none
      
      integer,  intent(in) :: n     ! Length of vectors
      real(dp), intent(in) :: dx(*) ! First vector
      integer,  intent(in) :: incx  ! Increment along dx
      real(dp), intent(in) :: dy(*) ! Second vector
      integer,  intent(in) :: incy  ! Increment along dy
      real(dp)             :: output
    end function
    
    ! Multiplies real vector by real scalar. Equivalent to dx *= da.
    subroutine dscal(n,da,dx,incx)
      import :: dp
      implicit none
      
      integer,  intent(in)    :: n     ! Length of vector
      real(dp), intent(in)    :: da    ! Scalar
      real(dp), intent(inout) :: dx(*) ! Vector
      integer,  intent(in)    :: incx  ! Increment along dx
    end subroutine
    
    ! Multiplies complex vector by complex scalar. Equivalent to zx *= za.
    subroutine zscal(n,za,zx,incx)
      import :: dp
      implicit none
      
      integer,     intent(in)    :: n     ! Length of vector
      complex(dp), intent(in)    :: za    ! Scalar
      complex(dp), intent(inout) :: zx(*) ! Vector
      integer,     intent(in)    :: incx  ! Increment along zx
    end subroutine
    
    ! Complex norm. Returns sqrt(x.x).
    function dznrm2(n,x,incx) result(output)
      import :: dp
      implicit none
      
      integer,     intent(in) :: n      ! Length of vector
      complex(dp), intent(in) :: x(*)   ! Vector
      integer,     intent(in) :: incx   ! Increment along x
      real(dp)                :: output ! Result
    end function
    
    ! Minimises the least-squares fit l=(a.x-b)**2.
    subroutine dgels(trans,m,n,nrhs,a,lda,b,ldb,work,lwork,info)
      import :: dp
      implicit none
      
      character(1), intent(in)    :: trans    ! n/t: if t, a is transposed.
      integer,      intent(in)    :: m        ! size(a,1)=size(b,1).
      integer,      intent(in)    :: n        ! size(a,2)=size(x,1).
      integer,      intent(in)    :: nrhs     ! size(b,2)=size(x,2).
      integer,      intent(in)    :: lda      ! The leading dimension of a.
      real(dp),     intent(inout) :: a(lda,*) ! The matrix a.
      integer,      intent(in)    :: ldb      ! The leading dimension of b.
      real(dp),     intent(inout) :: b(ldb,*) ! The matrix b.
      real(dp),     intent(out)   :: work(*)  ! work(1) = optimal lwork.
      integer,      intent(in)    :: lwork    ! size(work).
      integer,      intent(out)   :: info     ! 0 on success.
    end subroutine
    
    ! In-place LU factorisation. Required for dgetri. LU factorises a.
    subroutine dgetrf(m,n,a,lda,ipiv,info)
      import :: dp
      implicit none
      
      integer,  intent(in)    :: m        ! size(a,1).
      integer,  intent(in)    :: n        ! size(a,2).
      integer,  intent(in)    :: lda      ! The leading dimension of a.
      real(dp), intent(inout) :: a(lda,*) ! The matrix a.
      integer,  intent(out)   :: ipiv(*)  ! Pivot indices. size=min(n,m).
      integer,  intent(out)   :: info     ! 0 on success.
    end subroutine
    
    ! In-place matrix inversion. Inverts a.
    subroutine dgetri(n,a,lda,ipiv,work,lwork,info)
      import :: dp
      implicit none
      
      integer,  intent(in)    :: n        ! size(a,1)=size(a,2).
      integer,  intent(in)    :: lda      ! The leading dimension of a.
      real(dp), intent(inout) :: a(lda,*) ! The matrix to be inverted.
      integer,  intent(in)    :: ipiv(*)  ! Pivot indices. size=n.
      real(dp), intent(out)   :: work(*)  ! work(1) = optimal lwork.
      integer,  intent(in)    :: lwork    ! size(work).
      integer,  intent(out)   :: info     ! 0 on success.
    end subroutine
  end interface
contains

! ----------------------------------------------------------------------
! Vector and Matrix operations involving the private contents_ variable.
! ----------------------------------------------------------------------
! Assignment.
subroutine assign_IntVector_integers(output,input)
  implicit none
  
  type(IntVector), intent(out) :: output
  integer,         intent(in)  :: input(:)
  
  output%contents_ = input
end subroutine

subroutine assign_RealVector_reals(output,input)
  implicit none
  
  type(RealVector), intent(out) :: output
  real(dp),         intent(in)  :: input(:)
  
  output%contents_ = input
end subroutine

subroutine assign_ComplexVector_complexes(output,input)
  implicit none
  
  type(ComplexVector), intent(out) :: output
  complex(dp),         intent(in)  :: input(:)
  
  output%contents_ = input
end subroutine

subroutine assign_IntMatrix_integers(output,input)
  implicit none
  
  type(IntMatrix), intent(out) :: output
  integer,         intent(in)  :: input(:,:)
  
  output%contents_ = input
end subroutine

subroutine assign_RealMatrix_reals(output,input)
  implicit none
  
  type(RealMatrix), intent(out) :: output
  real(dp),         intent(in)  :: input(:,:)
  
  output%contents_ = input
end subroutine

subroutine assign_ComplexMatrix_complexes(output,input)
  implicit none
  
  type(ComplexMatrix), intent(out) :: output
  complex(dp),         intent(in)  :: input(:,:)
  
  output%contents_ = input
end subroutine

! Conversion to fundamental types. Effectively getters for contents_.
function int_IntVector(input) result(output)
  implicit none
  
  type(IntVector), intent(in) :: input
  integer, allocatable        :: output(:)
  
  if (allocated(input%contents_)) then
    output = input%contents_
  else
    call print_line(CODE_ERROR//': Trying to use the contents of a vector &
       &before it has been allocated.')
    call err()
  endif
end function

function int_IntMatrix(input) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: input
  integer, allocatable        :: output(:,:)
  
  if (allocated(input%contents_)) then
    output = input%contents_
  else
    call print_line(CODE_ERROR//': Trying to use the contents of a matrix &
       &before it has been allocated.')
    call err()
  endif
end function

function dble_RealVector(input) result(output)
  implicit none
  
  type(RealVector), intent(in) :: input
  real(dp), allocatable        :: output(:)
  
  if (allocated(input%contents_)) then
    output = input%contents_
  else
    call print_line(CODE_ERROR//': Trying to use the contents of a vector &
       &before it has been allocated.')
    call err()
  endif
end function

function dble_RealMatrix(input) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: input
  real(dp), allocatable        :: output(:,:)
  
  if (allocated(input%contents_)) then
    output = input%contents_
  else
    call print_line(CODE_ERROR//': Trying to use the contents of a matrix &
       &before it has been allocated.')
    call err()
  endif
end function

function cmplx_ComplexVector(input) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: input
  complex(dp), allocatable        :: output(:)
  
  if (allocated(input%contents_)) then
    output = input%contents_
  else
    call print_line(CODE_ERROR//': Trying to use the contents of a vector &
       &before it has been allocated.')
    call err()
  endif
end function

function cmplx_ComplexMatrix(input) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: input
  complex(dp), allocatable        :: output(:,:)
  
  if (allocated(input%contents_)) then
    output = input%contents_
  else
    call print_line(CODE_ERROR//': Trying to use the contents of a matrix &
       &before it has been allocated.')
    call err()
  endif
end function

! ----------------------------------------------------------------------
! Conversions between types.
! ----------------------------------------------------------------------
! Private conversion to other vector/matrix types.
function to_IntVector_IntVector(this) result(output)
  implicit none
  
  class(IntVector), intent(in) :: this
  type(IntVector)              :: output
  
  output = this
end function

function to_RealVector_IntVector(this) result(output)
  implicit none
  
  class(IntVector), intent(in) :: this
  type(RealVector)             :: output
  
  output = real(int(this),dp)
end function

function to_ComplexVector_IntVector(this) result(output)
  implicit none
  
  class(IntVector), intent(in) :: this
  type(ComplexVector)          :: output
  
  output = cmplx(int(this),0,dp)
end function

function to_RealVector_RealVector(this) result(output)
  implicit none
  
  class(RealVector), intent(in) :: this
  type(RealVector)              :: output
  
  output = this
end function

function to_ComplexVector_RealVector(this) result(output)
  implicit none
  
  class(RealVector), intent(in) :: this
  type(ComplexVector)           :: output
  
  output = cmplx(dble(this),0.0_dp,dp)
end function

function to_ComplexVector_ComplexVector(this) result(output)
  implicit none
  
  class(ComplexVector), intent(in) :: this
  type(ComplexVector)              :: output
  
  output = this
end function

function to_IntMatrix_IntMatrix(this) result(output)
  implicit none
  
  class(IntMatrix), intent(in) :: this
  type(IntMatrix)              :: output
  
  output = this
end function

function to_RealMatrix_IntMatrix(this) result(output)
  implicit none
  
  class(IntMatrix), intent(in) :: this
  type(RealMatrix)             :: output
  
  output = real(int(this),dp)
end function

function to_ComplexMatrix_IntMatrix(this) result(output)
  implicit none
  
  class(IntMatrix), intent(in) :: this
  type(ComplexMatrix)          :: output
  
  output = cmplx(int(this),0,dp)
end function

function to_RealMatrix_RealMatrix(this) result(output)
  implicit none
  
  class(RealMatrix), intent(in) :: this
  type(RealMatrix)              :: output
  
  output = this
end function

function to_ComplexMatrix_RealMatrix(this) result(output)
  implicit none
  
  class(RealMatrix), intent(in) :: this
  type(ComplexMatrix)           :: output
  
  output = cmplx(dble(this),0.0_dp,dp)
end function

function to_ComplexMatrix_ComplexMatrix(this) result(output)
  implicit none
  
  class(ComplexMatrix), intent(in) :: this
  type(ComplexMatrix)              :: output
  
  output = this
end function

! Public conversion to other vector/matrix types.
impure elemental function intvec_IntVectorable(input) result(output)
  implicit none
  
  class(IntVectorable), intent(in) :: input
  type(IntVector)                  :: output
  
  output = input%to_IntVector()
end function

impure elemental function intmat_IntMatrixable(input) result(output)
  implicit none
  
  class(IntMatrixable), intent(in) :: input
  type(IntMatrix)                  :: output
  
  output = input%to_IntMatrix()
end function

impure elemental function dblevec_RealVectorable(input) result(output)
  implicit none
  
  class(RealVectorable), intent(in) :: input
  type(RealVector)                  :: output
  
  output = input%to_RealVector()
end function

impure elemental function dblemat_RealMatrixable(input) result(output)
  implicit none
  
  class(RealMatrixable), intent(in) :: input
  type(RealMatrix)                  :: output
  
  output = input%to_RealMatrix()
end function

impure elemental function cmplxvec_ComplexVectorable(input) result(output)
  implicit none
  
  class(ComplexVectorable), intent(in) :: input
  type(ComplexVector)                  :: output
  
  output = input%to_ComplexVector()
end function

impure elemental function cmplxvec_RealVectorables(real,imag) result(output)
  implicit none
  
  class(RealVectorable), intent(in) :: real
  class(RealVectorable), intent(in) :: imag
  type(ComplexVector)               :: output
  
  output = cmplx(dble(real%to_RealVector()),dble(imag%to_RealVector()),dp)
end function

impure elemental function cmplxmat_ComplexMatrixable(input) result(output)
  implicit none
  
  class(ComplexMatrixable), intent(in) :: input
  type(ComplexMatrix)                  :: output
  
  output = input%to_ComplexMatrix()
end function

impure elemental function cmplxmat_RealMatrixables(real,imag) result(output)
  implicit none
  
  class(RealMatrixable), intent(in) :: real
  class(RealMatrixable), intent(in) :: imag
  type(ComplexMatrix)               :: output
  
  output = cmplx(dble(real%to_RealMatrix()),dble(imag%to_RealMatrix()),dp)
end function

! Conversion to Vector and Matrix from intrinsic types.
function vec_integers(input) result(output)
  implicit none
  
  integer, intent(in) :: input(:)
  type(IntVector)     :: output
  
  output = input
end function

function vec_reals(input) result(output)
  implicit none
  
  real(dp), intent(in) :: input(:)
  type(RealVector)     :: output
  
  output = input
end function

function vec_complexes(input) result(output)
  implicit none
  
  complex(dp), intent(in) :: input(:)
  type(ComplexVector)     :: output
  
  output = input
end function

function mat_integers(input) result(output)
  implicit none
  
  integer, intent(in) :: input(:,:)
  type(IntMatrix)     :: output
  
  output = input
end function

function mat_reals(input) result(output)
  implicit none
  
  real(dp), intent(in) :: input(:,:)
  type(RealMatrix)     :: output
  
  output = input
end function

function mat_complexes(input) result(output)
  implicit none
  
  complex(dp), intent(in) :: input(:,:)
  type(ComplexMatrix)     :: output
  
  output = input
end function

function mat_integers_shape(input,m,n) result(output)
  implicit none
  
  integer, intent(in) :: input(:)
  integer, intent(in) :: m
  integer, intent(in) :: n
  type(IntMatrix)     :: output
  
  output = transpose(reshape(input, [m,n]))
end function

function mat_reals_shape(input,m,n) result(output)
  implicit none
  
  real(dp), intent(in) :: input(:)
  integer,  intent(in) :: m
  integer,  intent(in) :: n
  type(RealMatrix)     :: output
  
  output = transpose(reshape(input, [m,n]))
end function

function mat_complexes_shape(input,m,n) result(output)
  implicit none
  
  complex(dp), intent(in) :: input(:)
  integer,     intent(in) :: m
  integer,     intent(in) :: n
  type(ComplexMatrix)     :: output
  
  output = transpose(reshape(input, [m,n]))
end function

! ----------------------------------------------------------------------
! All other Vector and Matrix operations.
! ----------------------------------------------------------------------
! N.B. the number of procedures accessing contents_ directly is intentionally
!    limited for stability reasons.
! The above procedures behave well if contents_ has not been allocated,
!    and this good behaviour is automatically passed to the procedures below.

! Makes a length n vector full of zeroes.
function zeroes_IntVector(n) result(output)
  implicit none
  
  integer, intent(in) :: n
  type(IntVector)     :: output
  
  integer :: i
  
  output = [(0,i=1,n)]
end function

! Makes a mxn matrix full of zeroes.
function zeroes_IntMatrix(m,n) result(output)
  implicit none
  
  integer, intent(in) :: m
  integer, intent(in) :: n
  type(IntMatrix)     :: output
  
  integer, allocatable :: contents(:,:)
  
  integer :: ialloc
  
  allocate(contents(m,n), stat=ialloc); call err(ialloc)
  contents = 0
  output = contents
end function

! Makes an nxn identity matrix.
function make_identity_matrix(n) result(output)
  implicit none
  
  integer, intent(in) :: n
  type(IntMatrix)     :: output
  
  integer, allocatable :: contents(:,:)
  
  integer :: i,ialloc
  
  allocate(contents(n,n), stat=ialloc); call err(ialloc)
  contents = 0
  do i=1,n
    contents(i,i) = 1
  enddo
  output = contents
end function

! Conversion to nearest integer.
function nint_RealVector(input) result(output)
  implicit none
  
  type(RealVector), intent(in) :: input
  integer, allocatable         :: output(:)
  
  output = nint(dble(input))
end function

function nint_RealMatrix(input) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: input
  integer, allocatable         :: output(:,:)
  
  output = nint(dble(input))
end function

! Makes a matrix whose rows are the input vectors.
function row_matrix_IntVectors(input) result(output)
  implicit none
  
  type(IntVector), intent(in) :: input(:)
  integer, allocatable        :: output(:,:)
  
  integer :: i,ialloc
  
  if (size(input)==0) then
    call print_line(CODE_ERROR//': Trying to make row matrix from empty &
       &array.')
    call err()
  endif
  
  do i=2,size(input)
    if (size(input(i))/=size(input(1))) then
      call print_line(CODE_ERROR//': Trying to make row matrix from &
         &inconsistent vectors.')
      call err()
    endif
  enddo
  
  allocate(output(size(input), size(input,1)), stat=ialloc); call err(ialloc)
  do i=1,size(input)
    output(i,:) = int(input(i))
  enddo
end function

function row_matrix_RealVectors(input) result(output)
  implicit none
  
  type(RealVector), intent(in) :: input(:)
  real(dp), allocatable        :: output(:,:)
  
  integer :: i,ialloc
  
  if (size(input)==0) then
    call print_line(CODE_ERROR//': Trying to make row matrix from empty &
       &array.')
    call err()
  endif
  
  do i=2,size(input)
    if (size(input(i))/=size(input(1))) then
      call print_line(CODE_ERROR//': Trying to make row matrix from &
         &inconsistent vectors.')
      call err()
    endif
  enddo
  
  allocate(output(size(input), size(input,1)), stat=ialloc); call err(ialloc)
  do i=1,size(input)
    output(i,:) = dble(input(i))
  enddo
end function

function row_matrix_ComplexVectors(input) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: input(:)
  complex(dp), allocatable        :: output(:,:)
  
  integer :: i,ialloc
  
  if (size(input)==0) then
    call print_line(CODE_ERROR//': Trying to make row matrix from empty &
       &array.')
    call err()
  endif
  
  do i=2,size(input)
    if (size(input(i))/=size(input(1))) then
      call print_line(CODE_ERROR//': Trying to make row matrix from &
         &inconsistent vectors.')
      call err()
    endif
  enddo
  
  allocate(output(size(input), size(input(1))), stat=ialloc); call err(ialloc)
  do i=1,size(input)
    output(i,:) = cmplx(input(i))
  enddo
end function

! Real part of a complex object.
impure elemental function real_ComplexVector(input) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: input
  type(RealVector)                :: output
  
  output = real(cmplx(input),dp)
end function

impure elemental function real_ComplexMatrix(input) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: input
  type(RealMatrix)                :: output
  
  output = real(cmplx(input),dp)
end function

! Imaginary part of a complex object.
impure elemental function aimag_ComplexVector(input) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: input
  type(RealVector)                :: output
  
  output = aimag(cmplx(input))
end function

impure elemental function aimag_ComplexMatrix(input) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: input
  type(RealMatrix)                :: output
  
  output = aimag(cmplx(input))
end function

! Element-wise abs of a complex object.
impure elemental function abs_ComplexVector(input) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: input
  type(RealVector)                :: output
  
  output = abs(cmplx(input))
end function

impure elemental function abs_ComplexMatrix(input) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: input
  type(RealMatrix)                :: output
  
  output = abs(cmplx(input))
end function

! Conjugate of a complex object.
impure elemental function conjg_ComplexVector(input) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: input
  type(ComplexVector)             :: output
  
  output = conjg(cmplx(input))
end function

impure elemental function conjg_ComplexMatrix(input) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: input
  type(ComplexMatrix)             :: output
  
  output = conjg(cmplx(input))
end function

! size().
function size_IntVector(input) result(output)
  implicit none
  
  type(IntVector), intent(in) :: input
  integer                     :: output
  
  output = size(int(input))
end function

function size_RealVector(input) result(output)
  implicit none
  
  type(RealVector), intent(in) :: input
  integer                      :: output
  
  output = size(dble(input))
end function

function size_ComplexVector(input) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: input
  integer                         :: output
  
  output = size(cmplx(input))
end function

function size_IntMatrix(input,dim) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: input
  integer,         intent(in) :: dim
  integer                     :: output
  
  output = size(int(input), dim)
end function

function size_RealMatrix(input,dim) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: input
  integer,          intent(in) :: dim
  integer                      :: output
  
  output = size(dble(input), dim)
end function

function size_ComplexMatrix(input,dim) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: input
  integer,             intent(in) :: dim
  integer                         :: output
  
  output = size(cmplx(input), dim)
end function

! Trace
function trace_IntMatrix(input) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: input
  integer                     :: output
  
  integer, allocatable :: contents(:,:)
  
  integer :: i
  
  contents = int(input)
  if (size(contents,1)/=size(contents,2)) then
    call print_line(CODE_ERROR//': Trying to take the trace of a non-square &
       &matrix.')
    call err()
  endif
  output = 0
  do i=1,size(contents,1)
    output = output + contents(i,i)
  enddo
end function

function trace_RealMatrix(input) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: input
  real(dp)                     :: output
  
  real(dp), allocatable :: contents(:,:)
  
  integer :: i
  
  contents = dble(input)
  if (size(contents,1)/=size(contents,2)) then
    call print_line(CODE_ERROR//': Trying to take the trace of a non-square &
       &matrix.')
    call err()
  endif
  output = 0
  do i=1,size(contents,1)
    output = output + contents(i,i)
  enddo
end function

function trace_ComplexMatrix(input) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: input
  complex(dp)                     :: output
  
  complex(dp), allocatable :: contents(:,:)
  
  integer :: i
  
  contents = cmplx(input)
  if (size(contents,1)/=size(contents,2)) then
    call print_line(CODE_ERROR//': Trying to take the trace of a non-square &
       &matrix.')
    call err()
  endif
  output = 0
  do i=1,size(contents,1)
    output = output + contents(i,i)
  enddo
end function

! Find the commutator of two matrices.
function commutator_IntMatrix_IntMatrix(this,that) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: this
  type(IntMatrix), intent(in) :: that
  type(IntMatrix)             :: output
  
  if (size(this,1)/=size(this,2)) then
    call print_line(CODE_ERROR//': Trying to find a commutator involving a &
       &non-square matrix.')
    call err()
  elseif (size(that,1)/=size(that,2)) then
    call print_line(CODE_ERROR//': Trying to find a commutator involving a &
       &non-square matrix.')
    call err()
  elseif (size(this,1)/=size(that,1)) then
    call print_line(CODE_ERROR//': Trying to find commutator of matrices of &
       &different sizes.')
    call err()
  endif
  
  output = this*that - that*this
end function

! Check if two matrices commute.
function matrices_commute_IntMatrix_IntMatrix(this,that) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: this
  type(IntMatrix), intent(in) :: that
  logical                     :: output
  
  output = commutator(this,that)==zeroes(size(this,1),size(this,1))
end function

! Equality.
impure elemental function equality_IntVector_IntVector(a,b) result(output)
  implicit none
  
  type(IntVector), intent(in) :: a
  type(IntVector), intent(in) :: b
  logical                     :: output
  
  output = all(int(a)==int(b))
end function

impure elemental function equality_IntMatrix_IntMatrix(a,b) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: a
  type(IntMatrix), intent(in) :: b
  logical                     :: output
  
  output = all(int(a)==int(b))
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

! Addition.
impure elemental function add_IntVector_IntVector(a,b) result(output)
  implicit none
  
  type(IntVector), intent(in) :: a
  type(IntVector), intent(in) :: b
  type(IntVector)             :: output
  
  output = int(a) + int(b)
end function

impure elemental function add_IntVector_RealVector(a,b) result(output)
  implicit none
  
  type(IntVector),  intent(in) :: a
  type(RealVector), intent(in) :: b
  type(RealVector)             :: output
  
  output = int(a) + dble(b)
end function

impure elemental function add_RealVector_IntVector(a,b) result(output)
  implicit none
  
  type(RealVector), intent(in) :: a
  type(IntVector),  intent(in) :: b
  type(RealVector)             :: output
  
  output = dble(a) + int(b)
end function

impure elemental function add_RealVector_RealVector(a,b) result(output)
  implicit none
  
  type(RealVector), intent(in) :: a
  type(RealVector), intent(in) :: b
  type(RealVector)             :: output
  
  output = dble(a) + dble(b)
end function

impure elemental function add_RealVector_ComplexVector(a,b) result(output)
  implicit none
  
  type(RealVector),    intent(in) :: a
  type(ComplexVector), intent(in) :: b
  type(ComplexVector)             :: output
  
  output = dble(a) + cmplx(b)
end function

impure elemental function add_ComplexVector_RealVector(a,b) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: a
  type(RealVector),    intent(in) :: b
  type(ComplexVector)             :: output
  
  output = cmplx(a) + dble(b)
end function

impure elemental function add_ComplexVector_ComplexVector(a,b) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: a
  type(ComplexVector), intent(in) :: b
  type(ComplexVector)             :: output
  
  output = cmplx(a) + cmplx(b)
end function

impure elemental function add_IntMatrix_IntMatrix(a,b) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: a
  type(IntMatrix), intent(in) :: b
  type(IntMatrix)             :: output
  
  output = int(a) + int(b)
end function

impure elemental function add_IntMatrix_RealMatrix(a,b) result(output)
  implicit none
  
  type(IntMatrix),  intent(in) :: a
  type(RealMatrix), intent(in) :: b
  type(RealMatrix)             :: output
  
  output = int(a) + dble(b)
end function

impure elemental function add_RealMatrix_IntMatrix(a,b) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: a
  type(IntMatrix),  intent(in) :: b
  type(RealMatrix)             :: output
  
  output = dble(a) + int(b)
end function

impure elemental function add_RealMatrix_RealMatrix(a,b) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: a
  type(RealMatrix), intent(in) :: b
  type(RealMatrix)             :: output
  
  output = dble(a) + dble(b)
end function

impure elemental function add_RealMatrix_ComplexMatrix(a,b) result(output)
  implicit none
  
  type(RealMatrix),    intent(in) :: a
  type(ComplexMatrix), intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = dble(a) + cmplx(b)
end function

impure elemental function add_ComplexMatrix_RealMatrix(a,b) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: a
  type(RealMatrix),    intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = cmplx(a) + dble(b)
end function

impure elemental function add_ComplexMatrix_ComplexMatrix(a,b) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: a
  type(ComplexMatrix), intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = cmplx(a) + cmplx(b)
end function

! Negative.
impure elemental function negative_IntVector(input) result(output)
  implicit none
  
  type(IntVector), intent(in) :: input
  type(IntVector)             :: output
  
  output = -int(input)
end function

impure elemental function negative_RealVector(input) result(output)
  implicit none
  
  type(RealVector), intent(in) :: input
  type(RealVector)             :: output
  
  output = -dble(input)
end function

impure elemental function negative_ComplexVector(input) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: input
  type(ComplexVector)             :: output
  
  output = -cmplx(input)
end function

impure elemental function negative_IntMatrix(input) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: input
  type(IntMatrix)             :: output
  
  output = -int(input)
end function

impure elemental function negative_RealMatrix(input) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: input
  type(RealMatrix)             :: output
  
  output = -dble(input)
end function

impure elemental function negative_ComplexMatrix(input) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: input
  type(ComplexMatrix)             :: output
  
  output = -cmplx(input)
end function

! Subtraction.
impure elemental function subtract_IntVector_IntVector(a,b) result(output)
  implicit none
  
  type(IntVector), intent(in) :: a
  type(IntVector), intent(in) :: b
  type(IntVector)             :: output
  
  output = int(a) - int(b)
end function

impure elemental function subtract_IntVector_RealVector(a,b) result(output)
  implicit none
  
  type(IntVector),  intent(in) :: a
  type(RealVector), intent(in) :: b
  type(RealVector)             :: output
  
  output = int(a) - dble(b)
end function

impure elemental function subtract_RealVector_IntVector(a,b) result(output)
  implicit none
  
  type(RealVector), intent(in) :: a
  type(IntVector),  intent(in) :: b
  type(RealVector)             :: output
  
  output = dble(a) - int(b)
end function

impure elemental function subtract_RealVector_RealVector(a,b) result(output)
  implicit none
  
  type(RealVector), intent(in) :: a
  type(RealVector), intent(in) :: b
  type(RealVector)             :: output
  
  output = dble(a) - dble(b)
end function

impure elemental function subtract_RealVector_ComplexVector(a,b) result(output)
  implicit none
  
  type(RealVector),    intent(in) :: a
  type(ComplexVector), intent(in) :: b
  type(ComplexVector)             :: output
  
  output = dble(a) - cmplx(b)
end function

impure elemental function subtract_ComplexVector_RealVector(a,b) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: a
  type(RealVector),    intent(in) :: b
  type(ComplexVector)             :: output
  
  output = cmplx(a) - dble(b)
end function

impure elemental function subtract_ComplexVector_ComplexVector(a,b) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: a
  type(ComplexVector), intent(in) :: b
  type(ComplexVector)             :: output
  
  output = cmplx(a) - cmplx(b)
end function

impure elemental function subtract_IntMatrix_IntMatrix(a,b) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: a
  type(IntMatrix), intent(in) :: b
  type(IntMatrix)             :: output
  
  output = int(a) - int(b)
end function

impure elemental function subtract_IntMatrix_RealMatrix(a,b) result(output)
  implicit none
  
  type(IntMatrix),  intent(in) :: a
  type(RealMatrix), intent(in) :: b
  type(RealMatrix)             :: output
  
  output = int(a) - dble(b)
end function

impure elemental function subtract_RealMatrix_IntMatrix(a,b) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: a
  type(IntMatrix),  intent(in) :: b
  type(RealMatrix)             :: output
  
  output = dble(a) - int(b)
end function

impure elemental function subtract_RealMatrix_RealMatrix(a,b) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: a
  type(RealMatrix), intent(in) :: b
  type(RealMatrix)             :: output
  
  output = dble(a) - dble(b)
end function

impure elemental function subtract_RealMatrix_ComplexMatrix(a,b) &
   & result(output)
  implicit none
  
  type(RealMatrix),    intent(in) :: a
  type(ComplexMatrix), intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = dble(a) - cmplx(b)
end function

impure elemental function subtract_ComplexMatrix_RealMatrix(a,b) &
   & result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: a
  type(RealMatrix),    intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = cmplx(a) - dble(b)
end function

impure elemental function subtract_ComplexMatrix_ComplexMatrix(a,b) &
   & result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: a
  type(ComplexMatrix), intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = cmplx(a) - cmplx(b)
end function

! Multiplication by a scalar.
impure elemental function multiply_IntVector_integer(a,b) result(output)
  implicit none
  
  type(IntVector), intent(in) :: a
  integer,         intent(in) :: b
  type(IntVector)             :: output
  
  output = int(a)*b
end function

impure elemental function multiply_integer_IntVector(a,b) result(output)
  implicit none
  
  integer,         intent(in) :: a
  type(IntVector), intent(in) :: b
  type(IntVector)             :: output
  
  output = a*int(b)
end function

impure elemental function multiply_IntVector_real(a,b) result(output)
  implicit none
  
  type(IntVector), intent(in) :: a
  real(dp),        intent(in) :: b
  type(RealVector)            :: output
  
  output = int(a)*b
end function

impure elemental function multiply_real_IntVector(a,b) result(output)
  implicit none
  
  real(dp),        intent(in) :: a
  type(IntVector), intent(in) :: b
  type(RealVector)            :: output
  
  output = a*int(b)
end function

impure elemental function multiply_IntVector_complex(a,b) result(output)
  implicit none
  
  type(IntVector), intent(in) :: a
  complex(dp),     intent(in) :: b
  type(ComplexVector)         :: output
  
  output = int(a)*b
end function

impure elemental function multiply_complex_IntVector(a,b) result(output)
  implicit none
  
  complex(dp),     intent(in) :: a
  type(IntVector), intent(in) :: b
  type(ComplexVector)         :: output
  
  output = a*int(b)
end function

impure elemental function multiply_RealVector_integer(a,b) result(output)
  implicit none
  
  type(RealVector), intent(in) :: a
  integer,          intent(in) :: b
  type(RealVector)             :: output
  
  output = dble(a)*b
end function

impure elemental function multiply_integer_RealVector(a,b) result(output)
  implicit none
  
  integer,          intent(in) :: a
  type(RealVector), intent(in) :: b
  type(RealVector)             :: output
  
  output = a*dble(b)
end function

impure elemental function multiply_RealVector_real(a,b) result(output)
  implicit none
  
  type(RealVector), intent(in) :: a
  real(dp),         intent(in) :: b
  type(RealVector)             :: output
  
  output = dble(a)*b
end function

impure elemental function multiply_real_RealVector(a,b) result(output)
  implicit none
  
  real(dp),         intent(in) :: a
  type(RealVector), intent(in) :: b
  type(RealVector)             :: output
  
  output = a*dble(b)
end function

impure elemental function multiply_RealVector_complex(a,b) result(output)
  implicit none
  
  type(RealVector), intent(in) :: a
  complex(dp),      intent(in) :: b
  type(ComplexVector)          :: output
  
  output = dble(a)*b
end function

impure elemental function multiply_complex_RealVector(a,b) result(output)
  implicit none
  
  complex(dp),      intent(in) :: a
  type(RealVector), intent(in) :: b
  type(ComplexVector)          :: output
  
  output = a*dble(b)
end function

impure elemental function multiply_ComplexVector_integer(a,b) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: a
  integer,             intent(in) :: b
  type(ComplexVector)             :: output
  
  output = cmplx(a)*b
end function

impure elemental function multiply_integer_ComplexVector(a,b) result(output)
  implicit none
  
  integer,             intent(in) :: a
  type(ComplexVector), intent(in) :: b
  type(ComplexVector)             :: output
  
  output = a*cmplx(b)
end function

impure elemental function multiply_ComplexVector_real(a,b) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: a
  real(dp),            intent(in) :: b
  type(ComplexVector)             :: output
  
  output = cmplx(a)*b
end function

impure elemental function multiply_real_ComplexVector(a,b) result(output)
  implicit none
  
  real(dp),            intent(in) :: a
  type(ComplexVector), intent(in) :: b
  type(ComplexVector)             :: output
  
  output = a*cmplx(b)
end function

impure elemental function multiply_ComplexVector_complex(a,b) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: a
  complex(dp),         intent(in) :: b
  type(ComplexVector)             :: output
  
  output = cmplx(a)*b
end function

impure elemental function multiply_complex_ComplexVector(a,b) result(output)
  implicit none
  
  complex(dp),         intent(in) :: a
  type(ComplexVector), intent(in) :: b
  type(ComplexVector)             :: output
  
  output = a*cmplx(b)
end function

impure elemental function multiply_IntMatrix_integer(a,b) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: a
  integer,         intent(in) :: b
  type(IntMatrix)             :: output
  
  output = int(a)*b
end function

impure elemental function multiply_integer_IntMatrix(a,b) result(output)
  implicit none
  
  integer,         intent(in) :: a
  type(IntMatrix), intent(in) :: b
  type(IntMatrix)             :: output
  
  output = a*int(b)
end function

impure elemental function multiply_IntMatrix_real(a,b) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: a
  real(dp),        intent(in) :: b
  type(RealMatrix)            :: output
  
  output = int(a)*b
end function

impure elemental function multiply_real_IntMatrix(a,b) result(output)
  implicit none
  
  real(dp),        intent(in) :: a
  type(IntMatrix), intent(in) :: b
  type(RealMatrix)            :: output
  
  output = a*int(b)
end function

impure elemental function multiply_IntMatrix_complex(a,b) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: a
  complex(dp),     intent(in) :: b
  type(ComplexMatrix)         :: output
  
  output = int(a)*b
end function

impure elemental function multiply_complex_IntMatrix(a,b) result(output)
  implicit none
  
  complex(dp),     intent(in) :: a
  type(IntMatrix), intent(in) :: b
  type(ComplexMatrix)         :: output
  
  output = a*int(b)
end function

impure elemental function multiply_RealMatrix_integer(a,b) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: a
  integer,          intent(in) :: b
  type(RealMatrix)             :: output
  
  output = dble(a)*b
end function

impure elemental function multiply_integer_RealMatrix(a,b) result(output)
  implicit none
  
  integer,          intent(in) :: a
  type(RealMatrix), intent(in) :: b
  type(RealMatrix)             :: output
  
  output = a*dble(b)
end function

impure elemental function multiply_RealMatrix_real(a,b) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: a
  real(dp),         intent(in) :: b
  type(RealMatrix)             :: output
  
  output = dble(a)*b
end function

impure elemental function multiply_real_RealMatrix(a,b) result(output)
  implicit none
  
  real(dp),         intent(in) :: a
  type(RealMatrix), intent(in) :: b
  type(RealMatrix)             :: output
  
  output = a*dble(b)
end function

impure elemental function multiply_RealMatrix_complex(a,b) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: a
  complex(dp),      intent(in) :: b
  type(ComplexMatrix)          :: output
  
  output = dble(a)*b
end function

impure elemental function multiply_complex_RealMatrix(a,b) result(output)
  implicit none
  
  complex(dp),      intent(in) :: a
  type(RealMatrix), intent(in) :: b
  type(ComplexMatrix)          :: output
  
  output = a*dble(b)
end function

impure elemental function multiply_ComplexMatrix_integer(a,b) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: a
  integer,             intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = cmplx(a)*b
end function

impure elemental function multiply_integer_ComplexMatrix(a,b) result(output)
  implicit none
  
  integer,             intent(in) :: a
  type(ComplexMatrix), intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = a*cmplx(b)
end function

impure elemental function multiply_ComplexMatrix_real(a,b) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: a
  real(dp),            intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = cmplx(a)*b
end function

impure elemental function multiply_real_ComplexMatrix(a,b) result(output)
  implicit none
  
  real(dp),            intent(in) :: a
  type(ComplexMatrix), intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = a*cmplx(b)
end function

impure elemental function multiply_ComplexMatrix_complex(a,b) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: a
  complex(dp),         intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = cmplx(a)*b
end function

impure elemental function multiply_complex_ComplexMatrix(a,b) result(output)
  implicit none
  
  complex(dp),         intent(in) :: a
  type(ComplexMatrix), intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = a*cmplx(b)
end function

! Dot products and matrix multiplication.
impure elemental function dot_IntVector_IntVector(a,b) result(output)
  implicit none
  
  type(IntVector), intent(in) :: a
  type(IntVector), intent(in) :: b
  integer                     :: output
  
  if (size(a)/=size(b)) then
    call print_line(CODE_ERROR//': Trying to take the dot product of two &
       &vectors of different lengths.')
    call err()
  endif
  
  output = dot_product(int(a), int(b))
end function

impure elemental function dot_IntVector_RealVector(a,b) result(output)
  implicit none
  
  type(IntVector),  intent(in) :: a
  type(RealVector), intent(in) :: b
  real(dp)                     :: output
  
  if (size(a)/=size(b)) then
    call print_line(CODE_ERROR//': Trying to take the dot product of two &
       &vectors of different lengths.')
    call err()
  endif
  
  output = dot_product(int(a), dble(b))
end function

impure elemental function dot_RealVector_IntVector(a,b) result(output)
  implicit none
  
  type(RealVector), intent(in) :: a
  type(IntVector),  intent(in) :: b
  real(dp)                     :: output
  
  if (size(a)/=size(b)) then
    call print_line(CODE_ERROR//': Trying to take the dot product of two &
       &vectors of different lengths.')
    call err()
  endif
  
  output = dot_product(dble(a), int(b))
end function

impure elemental function dot_RealVector_RealVector(a,b) result(output)
  implicit none
  
  type(RealVector), intent(in) :: a
  type(RealVector), intent(in) :: b
  real(dp)                     :: output
  
  if (size(a)/=size(b)) then
    call print_line(CODE_ERROR//': Trying to take the dot product of two &
       &vectors of different lengths.')
    call err()
  endif
  
  output = dot_product(dble(a), dble(b))
end function

impure elemental function dot_RealVector_ComplexVector(a,b) result(output)
  implicit none
  
  type(RealVector),    intent(in) :: a
  type(ComplexVector), intent(in) :: b
  complex(dp)                     :: output
  
  if (size(a)/=size(b)) then
    call print_line(CODE_ERROR//': Trying to take the dot product of two &
       &vectors of different lengths.')
    call err()
  endif
  
  output = dot_product(dble(a), cmplx(b))
end function

impure elemental function dot_ComplexVector_RealVector(a,b) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: a
  type(RealVector),    intent(in) :: b
  complex(dp)                     :: output
  
  if (size(a)/=size(b)) then
    call print_line(CODE_ERROR//': Trying to take the dot product of two &
       &vectors of different lengths.')
    call err()
  endif
  
  output = dot_product(cmplx(a), dble(b))
end function

impure elemental function dot_ComplexVector_ComplexVector(a,b) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: a
  type(ComplexVector), intent(in) :: b
  complex(dp)                     :: output
  
  if (size(a)/=size(b)) then
    call print_line(CODE_ERROR//': Trying to take the dot product of two &
       &vectors of different lengths.')
    call err()
  endif
  
  ! N.B. does not use dot_product since that takes the conjugate of the
  !    first argument.
  ! This behaviour would not be desirable, since (vec*mat)*vec should give the
  !    same result as vec*(mat*vec).
  output = sum(cmplx(a)*cmplx(b))
end function

impure elemental function dot_IntVector_IntMatrix(a,b) result(output)
  implicit none
  
  type(IntVector), intent(in) :: a
  type(IntMatrix), intent(in) :: b
  type(IntVector)             :: output
  
  if (size(a)/=size(b,1)) then
    call print_line(CODE_ERROR//': Trying to multiply a vector by a matrix &
       &of incompatible dimensions.')
    call err()
  endif
  
  output = matmul(int(a), int(b))
end function

impure elemental function dot_IntVector_RealMatrix(a,b) result(output)
  implicit none
  
  type(IntVector),  intent(in) :: a
  type(RealMatrix), intent(in) :: b
  type(RealVector)             :: output
  
  if (size(a)/=size(b,1)) then
    call print_line(CODE_ERROR//': Trying to multiply a vector by a matrix &
       &of incompatible dimensions.')
    call err()
  endif
  
  output = matmul(int(a), dble(b))
end function

impure elemental function dot_RealVector_IntMatrix(a,b) result(output)
  implicit none
  
  type(RealVector), intent(in) :: a
  type(IntMatrix),  intent(in) :: b
  type(RealVector)             :: output
  
  if (size(a)/=size(b,1)) then
    call print_line(CODE_ERROR//': Trying to multiply a vector by a matrix &
       &of incompatible dimensions.')
    call err()
  endif
  
  output = matmul(dble(a), int(b))
end function

impure elemental function dot_RealVector_RealMatrix(a,b) result(output)
  implicit none
  
  type(RealVector), intent(in) :: a
  type(RealMatrix), intent(in) :: b
  type(RealVector)             :: output
  
  if (size(a)/=size(b,1)) then
    call print_line(CODE_ERROR//': Trying to multiply a vector by a matrix &
       &of incompatible dimensions.')
    call err()
  endif
  
  output = matmul(dble(a), dble(b))
end function

impure elemental function dot_RealVector_ComplexMatrix(a,b) result(output)
  implicit none
  
  type(RealVector),    intent(in) :: a
  type(ComplexMatrix), intent(in) :: b
  type(ComplexVector)             :: output
  
  if (size(a)/=size(b,1)) then
    call print_line(CODE_ERROR//': Trying to multiply a vector by a matrix &
       &of incompatible dimensions.')
    call err()
  endif
  
  output = matmul(dble(a), cmplx(b))
end function

impure elemental function dot_ComplexVector_RealMatrix(a,b) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: a
  type(RealMatrix),    intent(in) :: b
  type(ComplexVector)             :: output
  
  if (size(a)/=size(b,1)) then
    call print_line(CODE_ERROR//': Trying to multiply a vector by a matrix &
       &of incompatible dimensions.')
    call err()
  endif
  
  output = matmul(cmplx(a), dble(b))
end function

impure elemental function dot_ComplexVector_ComplexMatrix(a,b) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: a
  type(ComplexMatrix), intent(in) :: b
  type(ComplexVector)             :: output
  
  if (size(a)/=size(b,1)) then
    call print_line(CODE_ERROR//': Trying to multiply a vector by a matrix &
       &of incompatible dimensions.')
    call err()
  endif
  
  output = matmul(cmplx(a), cmplx(b))
end function

impure elemental function dot_IntMatrix_IntVector(a,b) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: a
  type(IntVector), intent(in) :: b
  type(IntVector)             :: output
  
  if (size(a,2)/=size(b)) then
    call print_line(CODE_ERROR//': Trying to multiply a matrix by a vector &
       &of incompatible dimensions.')
    call err()
  endif
  
  output = matmul(int(a), int(b))
end function

impure elemental function dot_IntMatrix_RealVector(a,b) result(output)
  implicit none
  
  type(IntMatrix),  intent(in) :: a
  type(RealVector), intent(in) :: b
  type(RealVector)             :: output
  
  if (size(a,2)/=size(b)) then
    call print_line(CODE_ERROR//': Trying to multiply a matrix by a vector &
       &of incompatible dimensions.')
    call err()
  endif
  
  output = matmul(int(a), dble(b))
end function

impure elemental function dot_RealMatrix_IntVector(a,b) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: a
  type(IntVector),  intent(in) :: b
  type(RealVector)             :: output
  
  if (size(a,2)/=size(b)) then
    call print_line(CODE_ERROR//': Trying to multiply a matrix by a vector &
       &of incompatible dimensions.')
    call err()
  endif
  
  output = matmul(dble(a), int(b))
end function

impure elemental function dot_RealMatrix_RealVector(a,b) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: a
  type(RealVector), intent(in) :: b
  type(RealVector)             :: output
  
  if (size(a,2)/=size(b)) then
    call print_line(CODE_ERROR//': Trying to multiply a matrix by a vector &
       &of incompatible dimensions.')
    call err()
  endif
  
  output = matmul(dble(a), dble(b))
end function

impure elemental function dot_RealMatrix_ComplexVector(a,b) result(output)
  implicit none
  
  type(RealMatrix),    intent(in) :: a
  type(ComplexVector), intent(in) :: b
  type(ComplexVector)             :: output
  
  if (size(a,2)/=size(b)) then
    call print_line(CODE_ERROR//': Trying to multiply a matrix by a vector &
       &of incompatible dimensions.')
    call err()
  endif
  
  output = matmul(dble(a), cmplx(b))
end function

impure elemental function dot_ComplexMatrix_RealVector(a,b) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: a
  type(RealVector),    intent(in) :: b
  type(ComplexVector)             :: output
  
  if (size(a,2)/=size(b)) then
    call print_line(CODE_ERROR//': Trying to multiply a matrix by a vector &
       &of incompatible dimensions.')
    call err()
  endif
  
  output = matmul(cmplx(a), dble(b))
end function

impure elemental function dot_ComplexMatrix_ComplexVector(a,b) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: a
  type(ComplexVector), intent(in) :: b
  type(ComplexVector)             :: output
  
  if (size(a,2)/=size(b)) then
    call print_line(CODE_ERROR//': Trying to multiply a matrix by a vector &
       &of incompatible dimensions.')
    call err()
  endif
  
  output = matmul(cmplx(a), cmplx(b))
end function

impure elemental function dot_IntMatrix_IntMatrix(a,b) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: a
  type(IntMatrix), intent(in) :: b
  type(IntMatrix)             :: output
  
  if (size(a,2)/=size(b,1)) then
    call print_line(CODE_ERROR//': Trying to multiply two matrices of &
       & incompatible dimensions.')
    call err()
  endif
  
  output = matmul(int(a), int(b))
end function

impure elemental function dot_IntMatrix_RealMatrix(a,b) result(output)
  implicit none
  
  type(IntMatrix),  intent(in) :: a
  type(RealMatrix), intent(in) :: b
  type(RealMatrix)             :: output
  
  if (size(a,2)/=size(b,1)) then
    call print_line(CODE_ERROR//': Trying to multiply two matrices of &
       & incompatible dimensions.')
    call err()
  endif
  
  output = matmul(int(a), dble(b))
end function

impure elemental function dot_RealMatrix_IntMatrix(a,b) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: a
  type(IntMatrix),  intent(in) :: b
  type(RealMatrix)             :: output
  
  if (size(a,2)/=size(b,1)) then
    call print_line(CODE_ERROR//': Trying to multiply two matrices of &
       & incompatible dimensions.')
    call err()
  endif
  
  output = matmul(dble(a), int(b))
end function

impure elemental function dot_RealMatrix_RealMatrix(a,b) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: a
  type(RealMatrix), intent(in) :: b
  type(RealMatrix)             :: output
  
  if (size(a,2)/=size(b,1)) then
    call print_line(CODE_ERROR//': Trying to multiply two matrices of &
       & incompatible dimensions.')
    call err()
  endif
  
  output = matmul(dble(a), dble(b))
end function

impure elemental function dot_RealMatrix_ComplexMatrix(a,b) result(output)
  implicit none
  
  type(RealMatrix),    intent(in) :: a
  type(ComplexMatrix), intent(in) :: b
  type(ComplexMatrix)             :: output
  
  if (size(a,2)/=size(b,1)) then
    call print_line(CODE_ERROR//': Trying to multiply two matrices of &
       & incompatible dimensions.')
    call err()
  endif
  
  output = matmul(dble(a), cmplx(b))
end function

impure elemental function dot_ComplexMatrix_RealMatrix(a,b) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: a
  type(RealMatrix),    intent(in) :: b
  type(ComplexMatrix)             :: output
  
  if (size(a,2)/=size(b,1)) then
    call print_line(CODE_ERROR//': Trying to multiply two matrices of &
       & incompatible dimensions.')
    call err()
  endif
  
  output = matmul(cmplx(a), dble(b))
end function

impure elemental function dot_ComplexMatrix_ComplexMatrix(a,b) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: a
  type(ComplexMatrix), intent(in) :: b
  type(ComplexMatrix)             :: output
  
  if (size(a,2)/=size(b,1)) then
    call print_line(CODE_ERROR//': Trying to multiply two matrices of &
       & incompatible dimensions.')
    call err()
  endif
  
  output = matmul(cmplx(a), cmplx(b))
end function

! Division by scalar.
impure elemental function divide_IntVector_integer(a,b) result(output)
  implicit none
  
  type(IntVector), intent(in) :: a
  integer,         intent(in) :: b
  type(IntVector)             :: output
  
  output = int(a)/b
end function

impure elemental function divide_IntVector_real(a,b) result(output)
  implicit none
  
  type(IntVector), intent(in) :: a
  real(dp),        intent(in) :: b
  type(RealVector)            :: output
  
  output = int(a)/b
end function

impure elemental function divide_IntVector_complex(a,b) result(output)
  implicit none
  
  type(IntVector), intent(in) :: a
  complex(dp),     intent(in) :: b
  type(ComplexVector)         :: output
  
  output = int(a)/b
end function

impure elemental function divide_RealVector_integer(a,b) result(output)
  implicit none
  
  type(RealVector), intent(in) :: a
  integer,          intent(in) :: b
  type(RealVector)             :: output
  
  output = dble(a)/b
end function

impure elemental function divide_RealVector_real(a,b) result(output)
  implicit none
  
  type(RealVector), intent(in) :: a
  real(dp),         intent(in) :: b
  type(RealVector)             :: output
  
  output = dble(a)/b
end function

impure elemental function divide_RealVector_complex(a,b) result(output)
  implicit none
  
  type(RealVector),    intent(in) :: a
  complex(dp),         intent(in) :: b
  type(ComplexVector)             :: output
  
  output = dble(a)/b
end function

impure elemental function divide_ComplexVector_integer(a,b) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: a
  integer,             intent(in) :: b
  type(ComplexVector)             :: output
  
  output = cmplx(a)/b
end function

impure elemental function divide_ComplexVector_real(a,b) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: a
  real(dp),            intent(in) :: b
  type(ComplexVector)             :: output
  
  output = cmplx(a)/b
end function

impure elemental function divide_ComplexVector_complex(a,b) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: a
  complex(dp),         intent(in) :: b
  type(ComplexVector)             :: output
  
  output = cmplx(a)/b
end function

impure elemental function divide_IntMatrix_integer(a,b) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: a
  integer,         intent(in) :: b
  type(IntMatrix)             :: output
  
  output = int(a)/b
end function

impure elemental function divide_IntMatrix_real(a,b) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: a
  real(dp),        intent(in) :: b
  type(RealMatrix)            :: output
  
  output = int(a)/b
end function

impure elemental function divide_IntMatrix_complex(a,b) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: a
  complex(dp),     intent(in) :: b
  type(ComplexMatrix)         :: output
  
  output = int(a)/b
end function

impure elemental function divide_RealMatrix_integer(a,b) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: a
  integer,          intent(in) :: b
  type(RealMatrix)             :: output
  
  output = dble(a)/b
end function

impure elemental function divide_RealMatrix_real(a,b) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: a
  real(dp),         intent(in) :: b
  type(RealMatrix)             :: output
  
  output = dble(a)/b
end function

impure elemental function divide_RealMatrix_complex(a,b) result(output)
  implicit none
  
  type(RealMatrix),    intent(in) :: a
  complex(dp),         intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = dble(a)/b
end function

impure elemental function divide_ComplexMatrix_integer(a,b) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: a
  integer,             intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = cmplx(a)/b
end function

impure elemental function divide_ComplexMatrix_real(a,b) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: a
  real(dp),            intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = cmplx(a)/b
end function

impure elemental function divide_ComplexMatrix_complex(a,b) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: a
  complex(dp),         intent(in) :: b
  type(ComplexMatrix)             :: output
  
  output = cmplx(a)/b
end function

! Sum.
function sum_IntVectors(input) result(output)
  implicit none
  
  type(IntVector), intent(in) :: input(:)
  type(IntVector)             :: output
  
  integer :: i
  
  if (size(input)==0) then
    call print_line(CODE_ERROR//': Trying to take the sum of an empty array.')
    call err()
  endif
  
  output = input(1)
  do i=2,size(input)
    output = output + input(i)
  enddo
end function

function sum_RealVectors(input) result(output)
  implicit none
  
  type(RealVector), intent(in) :: input(:)
  type(RealVector)             :: output
  
  integer :: i
  
  if (size(input)==0) then
    call print_line(CODE_ERROR//': Trying to take the sum of an empty array.')
    call err()
  endif
  
  output = input(1)
  do i=2,size(input)
    output = output + input(i)
  enddo
end function

function sum_ComplexVectors(input) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: input(:)
  type(ComplexVector)             :: output
  
  integer :: i
  
  if (size(input)==0) then
    call print_line(CODE_ERROR//': Trying to take the sum of an empty array.')
    call err()
  endif
  
  output = input(1)
  do i=2,size(input)
    output = output + input(i)
  enddo
end function

function sum_IntMatrices(input) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: input(:)
  type(IntMatrix)             :: output
  
  integer :: i
  
  if (size(input)==0) then
    call print_line(CODE_ERROR//': Trying to take the sum of an empty array.')
    call err()
  endif
  
  output = input(1)
  do i=2,size(input)
    output = output + input(i)
  enddo
end function

function sum_RealMatrices(input) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: input(:)
  type(RealMatrix)             :: output
  
  integer :: i
  
  if (size(input)==0) then
    call print_line(CODE_ERROR//': Trying to take the sum of an empty array.')
    call err()
  endif
  
  output = input(1)
  do i=2,size(input)
    output = output + input(i)
  enddo
end function

function sum_ComplexMatrices(input) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: input(:)
  type(ComplexMatrix)             :: output
  
  integer :: i
  
  if (size(input)==0) then
    call print_line(CODE_ERROR//': Trying to take the sum of an empty array.')
    call err()
  endif
  
  output = input(1)
  do i=2,size(input)
    output = output + input(i)
  enddo
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

! Outer product.
! N.B. the conjugate is NOT taken in the complex case.
function outer_product_RealVector_RealVector(a,b) result(output)
  implicit none
  
  type(RealVector), intent(in) :: a
  type(RealVector), intent(in) :: b
  type(RealMatrix)             :: output
  
  real(dp), allocatable :: a_contents(:)
  real(dp), allocatable :: b_contents(:)
  real(dp), allocatable :: contents(:,:)
  
  integer :: i,j,ialloc
  
  a_contents = dble(a)
  b_contents = dble(b)
  allocate( contents(size(a_contents),size(b_contents)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(b_contents)
    do j=1,size(a_contents)
      contents(j,i) = a_contents(j)*b_contents(i)
    enddo
  enddo
  output = contents
end function

function outer_product_ComplexVector_ComplexVector(a,b) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: a
  type(ComplexVector), intent(in) :: b
  type(ComplexMatrix)             :: output
  
  complex(dp), allocatable :: a_contents(:)
  complex(dp), allocatable :: b_contents(:)
  complex(dp), allocatable :: contents(:,:)
  
  integer :: i,j,ialloc
  
  a_contents = cmplx(a)
  b_contents = cmplx(b)
  allocate( contents(size(a_contents),size(b_contents)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(b_contents)
    do j=1,size(a_contents)
      contents(j,i) = a_contents(j)*b_contents(i)
    enddo
  enddo
  output = contents
end function

! Transpose.
function transpose_IntMatrix(input) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: input
  type(IntMatrix)             :: output
  
  output = transpose(int(input))
end function

function transpose_RealMatrix(input) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: input
  type(RealMatrix)             :: output
  
  output = transpose(dble(input))
end function

function transpose_ComplexMatrix(input) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: input
  type(ComplexMatrix)             :: output
  
  output = transpose(cmplx(input))
end function

! Hermitian conjugate.
function hermitian_ComplexMatrix(input) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: input
  type(ComplexMatrix)             :: output
  
  output = transpose(conjg(cmplx(input)))
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
  
  output = input(1,1)*(input(2,2)*input(3,3)-input(2,3)*input(3,2))&
        &+ input(1,2)*(input(2,3)*input(3,1)-input(2,1)*input(3,3))&
        &+ input(1,3)*(input(2,1)*input(3,2)-input(2,2)*input(3,1))
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
  
  output = input(1,1)*(input(2,2)*input(3,3)-input(2,3)*input(3,2))&
       & + input(1,2)*(input(2,3)*input(3,1)-input(2,1)*input(3,3))&
       & + input(1,3)*(input(2,1)*input(3,2)-input(2,2)*input(3,1))
end function

function determinant_IntMatrix(input) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: input
  integer                     :: output
  
  output = determinant(int(input))
end function

function determinant_RealMatrix(input) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: input
  real(dp)                     :: output
  
  output = determinant(dble(input))
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
  
  integer :: ialloc
  
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
    call print_line(ERROR//' in LU factorisation: dgetrf error code: '//info)
    call err()
  endif
  
  ! Get the size of the required workspace.
  allocate(work(n),stat=ialloc); call err(ialloc)
  call dgetri(n,output,n,ipiv,work,-1,info)
  if (info/=0) then
    call print_line(ERROR//' in LU factorisation: dgetrf error code: '//info)
    call err()
  endif
  lwork = nint(work(1))
  deallocate(work, stat=ialloc); call err(ialloc)
  allocate(work(lwork), stat=ialloc); call err(ialloc)
  
  ! Run matrix inversion.
  call dgetri(n,output,n,ipiv,work,lwork,info)
  if (info/=0) then
    call print_line(ERROR//'in LU factorisation: dgetrf error code: '//info)
    call err()
  endif
end function

function invert_RealMatrix(input) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: input
  type(RealMatrix)             :: output
  
  output = invert(dble(input))
end function

! --------------------------------------------------
! Finds x which minimises the least-squares fit l=(a.x-b)**2.
! --------------------------------------------------
function linear_least_squares_reals_reals(a,b) result(output)
  implicit none
  
  real(dp), intent(in) :: a(:,:)
  real(dp), intent(in) :: b(:)
  type(RealVector)     :: output
  
  ! LAPACK dgels variables.
  real(dp), allocatable :: a2(:,:)
  real(dp), allocatable :: b2(:)
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
  b2 = b
  
  ! Calculate optimal workspace size.
  allocate(work(2*m*n), stat=ialloc); call err(ialloc)
  call dgels('N',m,n,1,a2(1,1),m,b2(1),m,work(1),-1,info)
  if (info/=0) then
    call print_line(ERROR//'in linear regression: dgels error code: '//info)
    call err()
  endif
  lwork = nint(work(1))
  deallocate(work, stat=ialloc); call err(ialloc)
  allocate(work(lwork), stat=ialloc); call err(ialloc)
  
  ! Run linear least-squares optimisation.
  call dgels('N',m,n,1,a2(1,1),m,b2(1),m,work(1),lwork,info)
  if (info/=0) then
    call print_line(ERROR//'in linear regression: dgels error code: '//info)
    call err()
  endif
  
  output = b2
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

! ----------------------------------------------------------------------
! I/O overloads.
! ----------------------------------------------------------------------
subroutine read_IntVector(this,input)
  implicit none
  
  class(IntVector), intent(out) :: this
  type(String),     intent(in)  :: input
  
  select type(this); type is(IntVector)
    this = int(split(input))
  end select
end subroutine

function write_IntVector(this) result(output)
  implicit none
  
  class(IntVector), intent(in) :: this
  type(String)                 :: output
  
  select type(this); type is(IntVector)
    output = join(int(this))
  end select
end function

subroutine read_RealVector(this,input)
  implicit none
  
  class(RealVector), intent(out) :: this
  type(String),      intent(in)  :: input
  
  select type(this); type is(RealVector)
    this = dble(split(input))
  end select
end subroutine

function write_RealVector(this) result(output)
  implicit none
  
  class(RealVector), intent(in) :: this
  type(String)                  :: output
  
  select type(this); type is(RealVector)
    output = join(dble(this))
  end select
end function

subroutine read_ComplexVector(this,input)
  implicit none
  
  class(ComplexVector), intent(out) :: this
  type(String),         intent(in)  :: input
  
  select type(this); type is(ComplexVector)
    this = cmplx(split(input))
  end select
end subroutine

function write_ComplexVector(this) result(output)
  implicit none
  
  class(ComplexVector), intent(in) :: this
  type(String)                     :: output
  
  select type(this); type is(ComplexVector)
    output = join(cmplx(this))
  end select
end function

subroutine read_IntMatrix(this,input)
  implicit none
  
  class(IntMatrix), intent(out) :: this
  type(String),     intent(in)  :: input(:)
  
  integer, allocatable :: line(:)
  integer, allocatable :: contents(:,:)
  
  integer :: i,ialloc
  
  select type(this); type is(IntMatrix)
    if (size(input)==0) then
      allocate(contents(0,0), stat=ialloc); call err(ialloc)
    else
      line = int(split(input(1)))
      allocate( contents(size(input),size(line)), &
              & stat=ialloc); call err(ialloc)
      contents(1,:) = line
      do i=2,size(input)
        line = int(split(input(i)))
        if (size(line)/=size(contents,2)) then
          call print_line(ERROR//': Reading matrix: rows of different &
             &lengths.')
          call err()
        endif
        contents(i,:) = line
      enddo
    endif
    
    this = contents
  end select
end subroutine

function write_IntMatrix(this) result(output)
  implicit none
  
  Class(IntMatrix), intent(in) :: this
  type(String), allocatable    :: output(:)
  
  integer, allocatable :: contents(:,:)
  
  integer :: i,ialloc
  
  select type(this); type is(IntMatrix)
    contents = int(this)
    allocate(output(size(contents,1)), stat=ialloc); call err(ialloc)
    do i=1,size(contents,1)
      output(i) = join(contents(i,:))
    enddo
  end select
end function

subroutine read_RealMatrix(this,input)
  implicit none
  
  class(RealMatrix), intent(out) :: this
  type(String),      intent(in)  :: input(:)
  
  real(dp), allocatable :: line(:)
  real(dp), allocatable :: contents(:,:)
  
  integer :: i,ialloc
  
  select type(this); type is(RealMatrix)
    if (size(input)==0) then
      allocate(contents(0,0), stat=ialloc); call err(ialloc)
    else
      line = dble(split(input(1)))
      allocate( contents(size(input),size(line)), &
              & stat=ialloc); call err(ialloc)
      contents(1,:) = line
      do i=2,size(input)
        line = dble(split(input(i)))
        if (size(line)/=size(contents,2)) then
          call print_line(ERROR//': Reading matrix: rows of different &
             &lengths.')
          call err()
        endif
        contents(i,:) = line
      enddo
    endif
    
    this = contents
  end select
end subroutine

function write_RealMatrix(this) result(output)
  implicit none
  
  Class(RealMatrix), intent(in) :: this
  type(String), allocatable     :: output(:)
  
  real(dp), allocatable :: contents(:,:)
  
  integer :: i,ialloc
  
  select type(this); type is(RealMatrix)
    contents = dble(this)
    allocate(output(size(contents,1)), stat=ialloc); call err(ialloc)
    do i=1,size(contents,1)
      output(i) = join(contents(i,:))
    enddo
  end select
end function

subroutine read_ComplexMatrix(this,input)
  implicit none
  
  class(ComplexMatrix), intent(out) :: this
  type(String),         intent(in)  :: input(:)
  
  complex(dp), allocatable :: line(:)
  complex(dp), allocatable :: contents(:,:)
  
  integer :: i,ialloc
  
  select type(this); type is(ComplexMatrix)
    if (size(input)==0) then
      allocate(contents(0,0), stat=ialloc); call err(ialloc)
    else
      line = cmplx(split(input(1)))
      allocate( contents(size(input),size(line)), &
              & stat=ialloc); call err(ialloc)
      contents(1,:) = line
      do i=2,size(input)
        line = cmplx(split(input(i)))
        if (size(line)/=size(contents,2)) then
          call print_line(ERROR//': Reading matrix: rows of different &
             &lengths.')
          call err()
        endif
        contents(i,:) = line
      enddo
    endif
    
    this = contents
  end select
end subroutine

function write_ComplexMatrix(this) result(output)
  implicit none
  
  Class(ComplexMatrix), intent(in) :: this
  type(String), allocatable        :: output(:)
  
  complex(dp), allocatable :: contents(:,:)
  
  integer :: i,ialloc
  
  select type(this); type is(ComplexMatrix)
    contents = cmplx(this)
    allocate(output(size(contents,1)), stat=ialloc); call err(ialloc)
    do i=1,size(contents,1)
      output(i) = join(contents(i,:))
    enddo
  end select
end function
end module
