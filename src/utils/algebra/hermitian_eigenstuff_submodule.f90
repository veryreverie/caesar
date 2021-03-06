submodule (caesar_hermitian_eigenstuff_module) caesar_hermitian_eigenstuff_submodule
  use caesar_algebra_module
contains

module procedure diagonalise_symmetric_reals
  ! dsyev variables.
  integer               :: n       ! The size of the input.
  real(dp), allocatable :: a(:,:)  ! The input. Also eigenvectors.
  real(dp), allocatable :: w(:)    ! The eigenvalues.
  real(dp), allocatable :: work(:) ! Workspace.
  integer               :: lwork   ! Workspace.
  integer               :: info    ! Error code. 0 on success.
  
  ! Temporary variables.
  integer :: i,ialloc
  
  if (size(input,1)==0) then
    allocate(output(0), stat=ialloc); call err(ialloc)
    return
  endif
  
  if (size(input,1)/=size(input,2)) then
    call print_line(CODE_ERROR//': Trying to find the eigenvalues of a &
       &non-square matrix.')
    call err()
  endif
  
  if (present(basis)) then
    if (size(basis)==0) then
      allocate(output(0), stat=ialloc); call err(ialloc)
      return
    endif
    
    do i=1,size(basis)
      if (size(basis(i))/=size(input,1)) then
        call print_line(CODE_ERROR//': Trying to find the eigenavlues of a &
           &matrix in an incompatible basis.')
      endif
    enddo
  endif
  
  if (.not. present(basis)) then
    a = input
  else
    a = row_matrix(basis)
    a = dble(mat(a) * mat(input) * transpose(mat(a)))
  endif
  
  n = size(a,1)
  allocate( w(n),        &
          & work(3*n-1), &
          & stat=ialloc); call err(ialloc)
  
  ! calculate optimal lwork
  call dsyev( jobz  = 'V',  &
            & uplo  = 'U',  &
            & n     = n,    &
            & a     = a,    &
            & lda   = n,    &
            & w     = w,    &
            & work  = work, &
            & lwork = -1,   &
            & info  = info)
  if (info /= 0) then
    call print_line(ERROR//' in diagonalisation: dsyev error code: '//info)
    call err()
  endif
  lwork = nint(work(1))
  deallocate(work, stat=ialloc); call err(ialloc)
  allocate(work(lwork), stat=ialloc); call err(ialloc)
  
  ! calculate eigenstuff
  call dsyev( jobz  = 'V',   &
            & uplo  = 'U',   &
            & n     = n,     &
            & a     = a,     &
            & lda   = n,     &
            & w     = w,     &
            & work  = work,  &
            & lwork = lwork, &
            & info  = info)
  if (info /= 0) then
    call print_line(ERROR//' in diagonalisation: dsyev error code: '//info)
    call print_line('matrix:')
    call print_lines(mat(input))
    call err()
  endif
  
  allocate(output(n), stat=ialloc); call err(ialloc)
  do i=1,n
    output(i)%eval = w(i)
    if (.not. present(basis)) then
      output(i)%evec = a(:,i)
    else
      output(i)%evec = dble(sum(a(:,i)*basis))
    endif
  enddo
end procedure

module procedure diagonalise_symmetric_RealMatrix
  output = diagonalise_symmetric(dble(input), basis)
end procedure

module procedure diagonalise_hermitian_complexes
  ! zheev variables.
  integer                  :: n        ! The size of the input.
  complex(dp), allocatable :: a(:,:)   ! Input matrix. Also eigenvectors.
  real(dp),    allocatable :: w(:)     ! Eigenvalues.
  complex(dp), allocatable :: work(:)  ! Workspace.
  integer                  :: lwork    ! Workspace.
  real(dp),    allocatable :: rwork(:) ! Workspace.
  integer                  :: info     ! Error code. 0 on success.
  
  ! Temporary variables.
  integer :: i,ialloc
  
  if (size(input,1)==0) then
    allocate(output(0), stat=ialloc); call err(ialloc)
    return
  endif
  
  if (size(input,1)/=size(input,2)) then
    call print_line(CODE_ERROR//': Trying to find the eigenvalues of a &
       &non-square matrix.')
    call err()
  endif
  
  if (present(basis)) then
    if (size(basis)==0) then
      allocate(output(0), stat=ialloc); call err(ialloc)
      return
    endif
    
    do i=1,size(basis)
      if (size(basis(i))/=size(input,1)) then
        call print_line(CODE_ERROR//': Trying to find the eigenavlues of a &
           &matrix in an incompatible basis.')
      endif
    enddo
  endif
  
  if (.not. present(basis)) then
    a = input
  else
    a = row_matrix(basis)
    a = cmplx(mat(a) * mat(input) * hermitian(mat(a)))
  endif
  
  n = size(a,1)
  allocate( w(n),         &
          & rwork(3*n-2), &
          & work(3*n-1),  &
          & stat=ialloc); call err(ialloc)
  
  ! calculate optimal lwork
  call zheev( jobz  = 'V',   &
            & uplo  = 'U',   &
            & n     = n,     &
            & a     = a,     &
            & lda   = n,     &
            & w     = w,     &
            & work  = work,  &
            & lwork = -1,    &
            & rwork = rwork, &
            & info  = info)
  if (info /= 0) then
    call print_line(ERROR//' in diagonalisation: zheev error code: '//info)
    call err()
  endif
  lwork = nint(real(work(1)))
  deallocate(work, stat=ialloc); call err(ialloc)
  allocate(work(lwork), stat=ialloc); call err(ialloc)
  
  ! calculate eigenstuff
  call zheev( jobz  = 'V',   &
            & uplo  = 'U',   &
            & n     = n,     &
            & a     = a,     &
            & lda   = n,     &
            & w     = w,     &
            & work  = work,  &
            & lwork = lwork, &
            & rwork = rwork, &
            & info  = info)
  if (info /= 0) then
    call print_line(ERROR//' in diagonalisation: zheev error code: '//info)
    call err()
  endif
  
  allocate(output(n), stat=ialloc); call err(ialloc)
  do i=1,n
    output(i)%eval = w(i)
    
    if (.not. present(basis)) then
      output(i)%evec = a(:,i)
    else
      output(i)%evec = cmplx(sum(a(:,i)*basis))
    endif
  enddo
end procedure

module procedure diagonalise_hermitian_ComplexMatrix
  output = diagonalise_hermitian(cmplx(input), basis)
end procedure

module procedure diagonalise_complex_complexes
  ! zgeev variables.
  integer                  :: n        ! The size of the input.
  complex(dp), allocatable :: a(:,:)   ! The input.
  complex(dp), allocatable :: w(:)     ! Eigenvalues.
  complex(dp), allocatable :: vl(:,:)  ! Left eigenvectors.
  complex(dp), allocatable :: vr(:,:)  ! Right eigenvectors.
  complex(dp), allocatable :: work(:)  ! Workspace.
  integer                  :: lwork    ! Workspace.
  real(dp),    allocatable :: rwork(:) ! Workspace.
  integer                  :: info     ! Error code. 0 on success.
  
  ! Temporary variables.
  integer :: i,ialloc
  
  if (size(input,1)==0) then
    allocate(output(0), stat=ialloc); call err(ialloc)
    return
  endif
  
  if (size(input,1)/=size(input,2)) then
    call print_line(CODE_ERROR//': Trying to find the eigenvalues of a &
       &non-square matrix.')
    call err()
  endif
  
  if (present(basis)) then
    if (size(basis)==0) then
      allocate(output(0), stat=ialloc); call err(ialloc)
      return
    endif
    
    do i=1,size(basis)
      if (size(basis(i))/=size(input,1)) then
        call print_line(CODE_ERROR//': Trying to find the eigenavlues of a &
           &matrix in an incompatible basis.')
      endif
    enddo
  endif
  
  if (.not. present(basis)) then
    a = input
  else
    a = row_matrix(basis)
    a = cmplx(mat(a) * mat(input) * transpose(mat(a)))
  endif
  
  n = size(a,1)
  allocate( w(n),       &
          & vl(n,n),    &
          & vr(n,n),    &
          & rwork(2*n), &
          & work(2*n),  &
          & stat=ialloc); call err(ialloc)
  
  ! calculate optimal lwork.
  call zgeev( jobvl = 'V',   &
            & jobvr = 'V',   &
            & n     = n,     &
            & a     = a,     &
            & lda   = n,     &
            & w     = w,     &
            & vl    = vl,    &
            & ldvl  = n,     &
            & vr    = vr,    &
            & ldvr  = n,     &
            & work  = work,  &
            & lwork = -1,    &
            & rwork = rwork, &
            & info  = info)
  if (info /= 0) then
    call print_line(ERROR//' in diagonalisation: zgeev error code: '//info)
    call err()
  endif
  lwork = nint(real(work(1)))
  deallocate(work, stat=ialloc); call err(ialloc)
  allocate(work(lwork), stat=ialloc); call err(ialloc)
  
  ! calculate eigenstuff
  call zgeev( jobvl = 'V',   &
            & jobvr = 'V',   &
            & n     = n,     &
            & a     = a,     &
            & lda   = n,     &
            & w     = w,     &
            & vl    = vl,    &
            & ldvl  = n,     &
            & vr    = vr,    &
            & ldvr  = n,     &
            & work  = work,  &
            & lwork = lwork, &
            & rwork = rwork, &
            & info  = info)
  if (info /= 0) then
    call print_line(ERROR//' in diagonalisation: zheev error code: '//info)
    call err()
  endif
  
  allocate(output(n), stat=ialloc); call err(ialloc)
  do i=1,n
    output(i)%eval = w(i)
    if (.not. present(basis)) then
      output(i)%left_evec = vl(:,i)
      output(i)%right_evec = vr(:,i)
    else
      output(i)%left_evec = cmplx(sum(vl(:,i)*basis))
      output(i)%right_evec = cmplx(sum(vr(:,i)*conjg(basis)))
    endif
  enddo
end procedure

module procedure diagonalise_complex_ComplexMatrix
  output = diagonalise_complex(cmplx(input), basis)
end procedure
end submodule
