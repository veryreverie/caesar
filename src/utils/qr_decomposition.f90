! ======================================================================
! Routines for finding the QR decomposition of various matrices.
! ======================================================================
module qr_decomposition_module
  use constants_module, only : dp
  use string_module
  use io_module
  
  use linear_algebra_module
  implicit none
  
  type QRDecomposition
    complex(dp), allocatable :: q(:,:)
    complex(dp), allocatable :: r(:,:)
  end type
  
  interface qr_decomposition
    module procedure qr_decomposition_complexes
    module procedure qr_decomposition_ComplexMatrix
  end interface
  
  interface orthonormalise
    module procedure orthonormalise_ComplexVectors
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

! --------------------------------------------------
! Construct an orthonormal basis from the given vectors using QR decomposition.
! If min_length is present then only vectors which are longer than min_length
!    before normalisation are returned.
! --------------------------------------------------
function orthonormalise_ComplexVectors(input,min_length) result(output)
  implicit none
  
  type(ComplexVector), intent(in)           :: input(:)
  real(dp),            intent(in), optional :: min_length
  type(ComplexVector), allocatable          :: output(:)
  
  integer :: no_vectors
  integer :: vector_length
  
  complex(dp), allocatable :: matrix(:,:)
  type(QRDecomposition)    :: qr
  
  integer :: no_independent_vectors
  
  integer :: i,ialloc
  
  no_vectors = size(input)
  if (no_vectors==0) then
    output = [ComplexVector::]
    return
  endif
  
  vector_length = size(input(1))
  do i=2,no_vectors
    if (size(input(i))/=vector_length) then
      call print_line(CODE_ERROR//': Trying to orhonormalise vectors of &
         &inconsistent lengths.')
      call err()
    endif
  enddo
  
  allocate(matrix(vector_length,no_vectors), stat=ialloc); call err(ialloc)
  do i=1,size(input)
    matrix(:,i) = cmplx(input(i))
  enddo
  qr = qr_decomposition(matrix)
  
  no_independent_vectors = min(no_vectors,vector_length)
  if (present(min_length)) then
    do i=1,no_independent_vectors
      if (abs(l2_norm(vec(qr%r(i,:))))<min_length) then
        no_independent_vectors = i-1
        exit
      endif
    enddo
  endif
  
  allocate(output(no_independent_vectors), stat=ialloc); call err(ialloc)
  do i=1,no_independent_vectors
    output(i) = vec(qr%q(:,i))
  enddo
end function
end module
