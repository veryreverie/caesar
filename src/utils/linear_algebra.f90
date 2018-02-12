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
    procedure, public :: str => str_IntVector
  end type
  
  type, extends(Stringable) :: RealVector
    real(dp), allocatable, private :: contents_(:)
  contains
    procedure, public :: str => str_RealVector
  end type
  
  type, extends(Stringable) :: ComplexVector
    complex(dp), allocatable, private :: contents_(:)
  contains
    procedure, public :: str => str_ComplexVector
  end type
  
  type, extends(Printable) :: IntMatrix
    integer, allocatable, private :: contents_(:,:)
  contains
    procedure, public :: str => str_IntMatrix
  end type
  
  type, extends(Printable) :: RealMatrix
    real(dp), allocatable, private :: contents_(:,:)
  contains
    procedure, public :: str => str_RealMatrix
  end type
  
  type, extends(Printable) :: ComplexMatrix
    complex(dp), allocatable, private :: contents_(:,:)
  contains
    procedure, public :: str => str_ComplexMatrix
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
  
  interface int
    module procedure int_IntVector
    module procedure int_IntMatrix
  end interface
  
  interface nint
    module procedure nint_RealVector
    module procedure nint_RealMatrix
  end interface
  
  interface dble
    module procedure dble_RealVector
    module procedure dble_RealMatrix
  end interface
  
  interface dblevec
    module procedure dblevec_IntVector
  end interface
  
  interface dblemat
    module procedure dblemat_IntMatrix
  end interface
  
  interface cmplx
    module procedure cmplx_ComplexVector
    module procedure cmplx_ComplexMatrix
  end interface
  
  interface cmplxvec
    module procedure cmplxvec_IntVectors
    module procedure cmplxvec_RealVectors
  end interface
  
  interface cmplxmat
    module procedure cmplxmat_IntMatrices
    module procedure cmplxmat_RealMatrices
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
! All other Vector and Matrix operations.
! ----------------------------------------------------------------------
! N.B. the number of procedures accessing contents_ directly is intentionally
!    limited for stability reasons.
! The above procedures behave well if contents_ has not been allocated,
!    and this good behaviour is automatically passed to the procedures below.

! Conversion to Vector and Matrix.
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

! Makes a length n vector full of zeroes.
function zeroes_IntVector(n) result(output)
  implicit none
  
  integer, intent(in) :: n
  type(IntVector)     :: output
  
  integer, allocatable :: contents(:)
  
  integer :: ialloc
  
  allocate(contents(n), stat=ialloc); call err(ialloc)
  contents = 0
  output = contents
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

! Conversion to other fundamental types.
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

function dblevec_IntVector(input) result(output)
  implicit none
  
  type(IntVector), intent(in) :: input
  type(RealVector)            :: output
  
  output = real(int(input),dp)
end function

function dblemat_IntMatrix(input) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: input
  type(RealMatrix)            :: output
  
  output = real(int(input),dp)
end function

function cmplxvec_IntVectors(real,imag) result(output)
  implicit none
  
  type(IntVector), intent(in)           :: real
  type(IntVector), intent(in), optional :: imag
  type(ComplexVector)                   :: output
  
  if (present(imag)) then
    output = cmplx(int(real),int(imag),dp)
  else
    output = cmplx(int(real),0,dp)
  endif
end function

function cmplxmat_IntMatrices(real,imag) result(output)
  implicit none
  
  type(IntMatrix), intent(in)           :: real
  type(IntMatrix), intent(in), optional :: imag
  type(ComplexMatrix)                   :: output
  
  if (present(imag)) then
    output = cmplx(int(real),int(imag),dp)
  else
    output = cmplx(int(real),0,dp)
  endif
end function

function cmplxvec_RealVectors(real,imag) result(output)
  implicit none
  
  type(RealVector), intent(in)           :: real
  type(RealVector), intent(in), optional :: imag
  type(ComplexVector)                    :: output
  
  if (present(imag)) then
    output = cmplx(dble(real),dble(imag),dp)
  else
    output = cmplx(dble(real),0,dp)
  endif
end function

function cmplxmat_RealMatrices(real,imag) result(output)
  implicit none
  
  type(RealMatrix), intent(in)           :: real
  type(RealMatrix), intent(in), optional :: imag
  type(ComplexMatrix)                    :: output
  
  if (present(imag)) then
    output = cmplx(dble(real),dble(imag),dp)
  else
    output = cmplx(dble(real),0,dp)
  endif
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
    call print_line(ERROR//'in LU factorisation: dgetrf error code: '//info)
    call err()
  endif
  
  ! Get the size of the required workspace.
  allocate(work(n),stat=ialloc); call err(ialloc)
  call dgetri(n,output,n,ipiv,work,-1,info)
  if (info/=0) then
    call print_line(ERROR//'in LU factorisation: dgetrf error code: '//info)
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

! I/O overloads.
function str_IntVector(this) result(output)
  implicit none
  
  class(IntVector), intent(in) :: this
  type(String)                 :: output
  
  output = join(int(this))
end function

function str_RealVector(this) result(output)
  implicit none
  
  class(RealVector), intent(in) :: this
  type(String)                  :: output
  
  output = join(dble(this))
end function

function str_ComplexVector(this) result(output)
  implicit none
  
  class(ComplexVector), intent(in) :: this
  type(String)                     :: output
  
  output = join(cmplx(this))
end function

function str_IntMatrix(this) result(output)
  implicit none
  
  Class(IntMatrix), intent(in) :: this
  type(String), allocatable    :: output(:)
  
  integer, allocatable :: contents(:,:)
  
  integer :: i,ialloc
  
  contents = int(this)
  allocate(output(size(contents,1)), stat=ialloc); call err(ialloc)
  do i=1,size(contents,1)
    output(i) = join(contents(i,:))
  enddo
end function

function str_RealMatrix(this) result(output)
  implicit none
  
  Class(RealMatrix), intent(in) :: this
  type(String), allocatable     :: output(:)
  
  real(dp), allocatable :: contents(:,:)
  
  integer :: i,ialloc
  
  contents = dble(this)
  allocate(output(size(contents,1)), stat=ialloc); call err(ialloc)
  do i=1,size(contents,1)
    output(i) = join(contents(i,:))
  enddo
end function

function str_ComplexMatrix(this) result(output)
  implicit none
  
  Class(ComplexMatrix), intent(in) :: this
  type(String), allocatable        :: output(:)
  
  complex(dp), allocatable :: contents(:,:)
  
  integer :: i,ialloc
  
  contents = cmplx(this)
  allocate(output(size(contents,1)), stat=ialloc); call err(ialloc)
  do i=1,size(contents,1)
    output(i) = join(contents(i,:))
  enddo
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
end module
