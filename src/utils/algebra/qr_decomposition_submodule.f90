submodule (caesar_qr_decomposition_module) caesar_qr_decomposition_submodule
  use caesar_algebra_module
contains

module procedure qr_decomposition_reals
  ! LAPACK variables.
  integer               :: m       ! size(a,1).
  integer               :: n       ! size(a,2).
  integer               :: k       ! min(m,n).
  real(dp), allocatable :: a(:,:)  ! The input matrix. Also the output.
  real(dp), allocatable :: tau(:)  ! Further output coefficients.
  real(dp), allocatable :: work(:) ! Workspace.
  integer               :: lwork   ! Workspace.
  integer               :: info    ! Error code.
  
  ! Output variables.
  real(dp), allocatable :: q(:,:)
  real(dp), allocatable :: r(:,:)
  
  integer :: i,ialloc
  
  m = size(input,1)
  n = size(input,2)
  
  if (m==0 .or. n==0) then
    output%q = int(make_identity_matrix(m))
    allocate(output%r(m,n), stat=ialloc); call err(ialloc)
    return
  endif
  
  k = min(m,n)
  a = input
  lwork = max(1,n)
  allocate( tau(max(1,k)), &
          & work(lwork),   &
          & stat=ialloc); call err(ialloc)
  
  ! calculate optimal lwork for zgeqrf.
  call dgeqrf( m     = m,    &
             & n     = n,    &
             & a     = a,    &
             & lda   = m,    &
             & tau   = tau,  &
             & work  = work, &
             & lwork = -1,   &
             & info  = info)
  if (info /= 0) then
    call print_line(ERROR//' in QR decomposition: dgeqrf error code: '//info)
    call err()
  endif
  lwork = nint(real(work(1)))
  deallocate(work, stat=ialloc); call err(ialloc)
  allocate(work(lwork), stat=ialloc); call err(ialloc)
  
  ! Run QR decomposition.
  call dgeqrf( m     = m,     &
             & n     = n,     &
             & a     = a,     &
             & lda   = m,     &
             & tau   = tau,   &
             & work  = work,  &
             & lwork = lwork, &
             & info  = info)
  if (info /= 0) then
    call print_line(ERROR//' in QR decomposition: dgeqrf error code: '//info)
    call err()
  endif
  
  ! Split the output of dgeqrf into R
  !    and the reflectors required to calculate Q.
  !
  ! A = (R R R)  or  A = (R R R R)
  !     (Q R R)          (Q R R R)
  !     (Q Q R)          (Q Q R R)
  !     (Q Q Q)
  allocate( q(m,m), &
          & r(m,n), &
          & stat=ialloc); call err(ialloc)
  r = 0.0_dp
  do i=1,k
    r(i,i:) = a(i,i:)
  enddo
  
  q = 0.0_dp
  do i=1,k
    q(i+1:,i) = a(i+1:,i)
  enddo
  
  ! Calculate optimal lwork for dorgqr.
  lwork = max(1,m)
  deallocate(work, stat=ialloc); call err(ialloc)
  allocate(work(lwork), stat=ialloc); call err(ialloc)
  call dorgqr( m     = m,    &
             & n     = m,    &
             & k     = k,    &
             & a     = q,    &
             & lda   = m,    &
             & tau   = tau,  &
             & work  = work, &
             & lwork = -1,   &
             & info  = info)
  if (info /= 0) then
    call print_line(ERROR//' in QR decomposition: dorgqr error code: '//info)
    call err()
  endif
  lwork = nint(real(work(1)))
  deallocate(work, stat=ialloc); call err(ialloc)
  allocate(work(lwork), stat=ialloc); call err(ialloc)
  
  ! Assemble Q from reflectors.
  call dorgqr( m     = m,     &
             & n     = m,     &
             & k     = k,     &
             & a     = q,     &
             & lda   = m,     &
             & tau   = tau,   &
             & work  = work,  &
             & lwork = lwork, &
             & info  = info)
  if (info /= 0) then
    call print_line(ERROR//' in QR decomposition: dorgqr error code: '//info)
    call err()
  endif
  
  ! Construct output.
  output%q = q
  output%r = r
end procedure

module procedure qr_decomposition_RealMatrix
  output = qr_decomposition(dble(input))
end procedure

module procedure qr_decomposition_complexes
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
  
  if (m==0 .or. n==0) then
    output%q = int(make_identity_matrix(m))
    allocate(output%r(m,n), stat=ialloc); call err(ialloc)
    return
  endif
  
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
    call print_line(ERROR//' in QR decomposition: zgeqrf error code: '//info)
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
    call print_line(ERROR//' in QR decomposition: zgeqrf error code: '//info)
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
    call print_line(ERROR//' in QR decomposition: zungqr error code: '//info)
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
  if (info /= 0) then
    call print_line(ERROR//' in QR decomposition: zungqr error code: '//info)
    call err()
  endif
  
  ! Construct output.
  output%q = q
  output%r = r
end procedure

module procedure qr_decomposition_ComplexMatrix
  output = qr_decomposition(cmplx(input))
end procedure
end submodule
