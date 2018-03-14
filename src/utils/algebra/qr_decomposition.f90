! ======================================================================
! Routines for finding the QR decomposition of various matrices.
! ======================================================================
module qr_decomposition_submodule
  use precision_module
  use io_module
  
  use linear_algebra_submodule
  implicit none
  
  private
  
  public :: QRDecomposition
  public :: qr_decomposition
  public :: intersection_basis
  
  type QRDecomposition
    complex(dp), allocatable :: q(:,:)
    complex(dp), allocatable :: r(:,:)
  end type
  
  interface qr_decomposition
    module procedure qr_decomposition_complexes
    module procedure qr_decomposition_ComplexMatrix
  end interface
  
  interface intersection_basis
    module procedure intersection_basis_ComplexVectors
  end interface
  
  interface
    ! Find the QR decomposition of a complex matrix.
    ! N.B. the output is given in terms of reflectors.
    subroutine zgeqrf(m,n,a,lda,tau,work,lwork,info)
      import :: dp
      implicit none
      
      integer,     intent(in)    :: m        ! size(a,1).
      integer,     intent(in)    :: n        ! size(a,2).
      integer,     intent(in)    :: lda      ! The leading dimension of a.
      complex(dp), intent(inout) :: a(lda,*) ! Matrix a. Also R and reflectors.
      complex(dp), intent(out)   :: tau(*)   ! Factors of reflectors.
      complex(dp), intent(out)   :: work(*)  ! work(1) = optimal lwork.
      integer,     intent(in)    :: lwork    ! The length of work.
      integer,     intent(out)   :: info     ! 0 on success.
    end subroutine
    
    ! Convert reflectors into the Q matrix.
    subroutine zungqr(m,n,k,a,lda,tau,work,lwork,info)
      import :: dp
      implicit none
      
      integer,     intent(in)    :: m        ! size(a,1).
      integer,     intent(in)    :: n        ! size(a,2).
      integer,     intent(in)    :: k        ! No. reflectors.
      integer,     intent(in)    :: lda      ! The leading dimension of a.
      complex(dp), intent(inout) :: a(lda,*) ! The matrix a. Also the output.
      complex(dp), intent(in)    :: tau(*)   ! Factors of reflectors.
      complex(dp), intent(out)   :: work(*)  ! work(1) = optimal lwork.
      integer,     intent(in)    :: lwork    ! The length of work.
      integer,     intent(out)   :: info     ! 0 on success.
    end subroutine
  end interface
contains

! --------------------------------------------------
! Calculates the QR decomposition of a complex matrix.
! --------------------------------------------------
function qr_decomposition_complexes(input) result(output)
  implicit none
  
  complex(dp), intent(in) :: input(:,:)
  type(QRDecomposition)   :: output
  
  ! LAPACK variables.
  integer                  :: m       ! size(a,1).
  integer                  :: n       ! size(a,2).
  integer                  :: k       ! min(m,n).
  complex(dp), allocatable :: a(:,:)  ! The input matrix. Also the output.
  complex(dp), allocatable :: tau(:)  ! Further output coefficients.
  complex(dp), allocatable :: work(:) ! Workspace.
  integer                  :: lwork   ! Workspace.
  integer                  :: info    ! Error code.
  
  ! Output variables.
  complex(dp), allocatable :: q(:,:)
  complex(dp), allocatable :: r(:,:)
  
  integer :: i,ialloc
  
  m = size(input,1)
  n = size(input,2)
  k = min(m,n)
  a = input
  lwork = max(1,n)
  allocate( tau(max(1,k)), &
          & work(lwork),   &
          & stat=ialloc); call err(ialloc)
  
  ! calculate optimal lwork for zgeqrf.
  call zgeqrf( m     = m,    &
             & n     = n,    &
             & a     = a,    &
             & lda   = m,    &
             & tau   = tau,  &
             & work  = work, &
             & lwork = -1,   &
             & info  = info)
  if (info /= 0) then
    call print_line(ERROR//'in QR decomposition: zgeqrf error code: '//info)
    call err()
  endif
  lwork = nint(real(work(1)))
  deallocate(work, stat=ialloc); call err(ialloc)
  allocate(work(lwork), stat=ialloc); call err(ialloc)
  
  ! Run QR decomposition.
  call zgeqrf( m     = m,     &
             & n     = n,     &
             & a     = a,     &
             & lda   = m,     &
             & tau   = tau,   &
             & work  = work,  &
             & lwork = lwork, &
             & info  = info)
  if (info /= 0) then
    call print_line(ERROR//'in QR decomposition: zgeqrf error code: '//info)
    call err()
  endif
  
  ! Split the output of zgeqrf into R
  !    and the reflectors required to calculate Q.
  !
  ! A = (R R R)  or  A = (R R R R)
  !     (Q R R)          (Q R R R)
  !     (Q Q R)          (Q Q R R)
  !     (Q Q Q)
  allocate( q(m,m), &
          & r(m,n), &
          & stat=ialloc); call err(ialloc)
  r = cmplx(0.0_dp,0.0_dp,dp)
  do i=1,k
    r(i,i:) = a(i,i:)
  enddo
  
  q = cmplx(0.0_dp,0.0_dp,dp)
  do i=1,k
    q(i+1:,i) = a(i+1:,i)
  enddo
  
  ! Calculate optimal lwork for zungqr.
  lwork = max(1,m)
  deallocate(work, stat=ialloc); call err(ialloc)
  allocate(work(lwork), stat=ialloc); call err(ialloc)
  call zungqr( m     = m,    &
             & n     = m,    &
             & k     = k,    &
             & a     = q,    &
             & lda   = m,    &
             & tau   = tau,  &
             & work  = work, &
             & lwork = -1,   &
             & info  = info)
  if (info /= 0) then
    call print_line(ERROR//'in QR decomposition: zungqr error code: '//info)
    call err()
  endif
  lwork = nint(real(work(1)))
  deallocate(work, stat=ialloc); call err(ialloc)
  allocate(work(lwork), stat=ialloc); call err(ialloc)
  
  ! Assemble Q from reflectors.
  call zungqr( m     = m,     &
             & n     = m,     &
             & k     = k,     &
             & a     = q,     &
             & lda   = m,     &
             & tau   = tau,   &
             & work  = work,  &
             & lwork = lwork, &
             & info  = info)
  
  ! Construct output.
  output%q = q
  output%r = r
end function

function qr_decomposition_ComplexMatrix(input) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: input
  type(QRDecomposition)           :: output
  
  output = qr_decomposition(cmplx(input))
end function

! ----------------------------------------------------------------------
! Uses the Zassenhaus algorithm, based on, QR factorisation, to find
!    an orthonormal basis for the intersection of two vector subspaces,
!    given an orthonormal basis for each subspace.
! ----------------------------------------------------------------------
function intersection_basis_ComplexVectors(a,b) result(output)
  use logic_module
  implicit none
  
  type(ComplexVector), intent(in)  :: a(:)
  type(ComplexVector), intent(in)  :: b(:)
  type(ComplexVector), allocatable :: output(:)
  
  ! Working variables for QR decomposition.
  complex(dp), allocatable :: matrix(:,:)
  type(QRDecomposition)    :: qr
  
  ! The basis vectors, and their projections onto the input vectors.
  type(ComplexVector), allocatable :: left_vectors(:)
  real(dp),            allocatable :: left_projections(:)
  type(ComplexVector), allocatable :: right_vectors(:)
  real(dp),            allocatable :: right_projections(:)
  
  integer :: vector_length
  integer :: union_size
  integer :: i,j,ialloc
  
  ! Return early if nothing to process.
  if (size(a)==0 .or. size(b)==0) then
    output = [ComplexVector::]
    return
  endif
  
  ! Check vector dimensionalities are consistent, and record said lengths.
  vector_length = size(a(1))
  do i=2,size(a)
    if (size(a(i))/=vector_length) then
      call print_line(CODE_ERROR//': Trying to orhonormalise vectors of &
         &inconsistent lengths.')
      call err()
    endif
  enddo
  do i=1,size(b)
    if (size(b(i))/=vector_length) then
      call print_line(CODE_ERROR//': Trying to orhonormalise vectors of &
         &inconsistent lengths.')
      call err()
    endif
  enddo
  
  ! Arrange input vectors into matrix and perform QR decomposition.
  ! The Zassenhaus matrix is given by Z = (A A)
  !                                       (B 0)
  ! Where the rows of A and B are or vectors of a and b respectively.
  !
  ! After the QR factorisation of Z, R is given by R = (U *)
  !                                                    (0 I)
  ! Where:
  !  The rows of U are orthogonal vectors spanning the union of a and b.
  !  The rows of I are orthogonal vectors spanning the intersection of a and b.
  allocate( matrix(size(a)+size(b),2*vector_length), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(a)
    matrix(i,                :vector_length) = cmplx(a(i))
    matrix(i, vector_length+1:             ) = cmplx(a(i))
  enddo
  do i=1,size(b)
    matrix(size(a)+i,                :vector_length) = cmplx(b(i))
    matrix(size(a)+i, vector_length+1:             ) = cmplx(0.0_dp,0.0_dp,dp)
  enddo
  qr = qr_decomposition(matrix)
  
  ! Extract the vectors in [U,0] and [*,I], and calculate their lengths.
  allocate( left_vectors(size(a)+size(b)),      &
          & left_projections(size(a)+size(b)),  &
          & right_vectors(size(a)+size(b)),     &
          & right_projections(size(a)+size(b)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(a)+size(b)
    left_vectors(i) = qr%r(i,:vector_length)
    left_projections(i) = l2_norm(left_vectors(i))
    
    right_vectors(i) = qr%r(i,vector_length+1:)
    right_projections(i) = l2_norm(right_vectors(i))
  enddo
  
  ! Calculate the size of U by finding the first row of 0.
  union_size = last(left_projections>0.1_dp)
  
  ! Check that the U is larger than I, that the rows of 0 are zero,
  !    and that the rows of I are not small.
  if (union_size<size(matrix,1)-union_size) then
    call print_line(CODE_ERROR//': the union is smaller than the &
       &intersection. Check that a and b are orthonormal.')
    call err()
  elseif (any(left_projections(union_size+1:)>1e-10_dp)) then
    call print_line(CODE_ERROR//': non-zero vector below union in reduced &
       &matrix. Check a and b are orthonormal.')
    call err()
  elseif (any(right_projections(union_size+1:)<0.1_dp)) then
    call print_line(CODE_ERROR//': small vector in intersection. Check a and &
       &b are orthonormal.')
    call err()
  endif
  
  ! Construct the output, the normalised rows of I.
  output = right_vectors(union_size+1:) / right_projections(union_size+1:)
end function
end module
