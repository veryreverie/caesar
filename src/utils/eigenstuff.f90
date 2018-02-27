! ======================================================================
! Routines for finding the eigenvalues and eigenvectors of various matrices.
! ======================================================================
module eigenstuff_module
  use constants_module, only : dp
  use string_module
  use io_module
  
  use linear_algebra_module
  use phase_module
  implicit none
  
  ! --------------------------------------------------
  ! The eigen(values and vectors) of a matrix.
  ! --------------------------------------------------
  type SymmetricEigenstuff
    real(dp)              :: eval
    real(dp), allocatable :: evec(:)
  end type
  
  type HermitianEigenstuff
    real(dp)                 :: eval
    complex(dp), allocatable :: evec(:)
  end type
  
  type ComplexEigenstuff
    complex(dp)              :: eval
    complex(dp), allocatable :: left_evec(:)
    complex(dp), allocatable :: right_evec(:)
  end type
  
  type UnitaryEigenstuff
    type(PhaseData)          :: eval
    complex(dp), allocatable :: evec(:)
  end type
  
  ! --------------------------------------------------
  ! Interfaces.
  ! --------------------------------------------------
  interface diagonalise_symmetric
    module procedure diagonalise_symmetric_reals
    module procedure diagonalise_symmetric_RealMatrix
  end interface
  
  interface diagonalise_hermitian
    module procedure diagonalise_hermitian_complexes
    module procedure diagonalise_hermitian_ComplexMatrix
  end interface
  
  interface diagonalise_complex
    module procedure diagonalise_complex_complexes
    module procedure diagonalise_complex_ComplexMatrix
  end interface
  
  interface diagonalise_unitary
    module procedure diagonalise_unitary_complexes
    module procedure diagonalise_unitary_ComplexMatrix
  end interface
  
  ! BLAS / LAPACK interface.
  interface
    ! Finds the eigenvalues of a real symmetric matrix.
    subroutine dsyev(jobz,uplo,n,a,lda,w,work,lwork,info)
      import :: dp
      implicit none
      
      character(1), intent(in)    :: jobz     ! N/V: if v, calculate eigenvecs.
      character(1), intent(in)    :: uplo     ! U/L: upper/lower triangle.
      integer,      intent(in)    :: n        ! The order of a.
      integer,      intent(in)    :: lda      ! The dimension of a.
      real(dp),     intent(inout) :: a(lda,*) ! Symmetric matrix.
      real(dp),     intent(out)   :: w(*)     ! Eigenvalues of a.
      real(dp),     intent(out)   :: work(*)  ! work(1) = optimal lwork.
      integer,      intent(in)    :: lwork    ! The length of work.
      integer,      intent(out)   :: info     ! 0 on success.
    end subroutine
    
    ! Finds the eigenvalues of a complex hermitian matrix.
    subroutine zheev(jobz,uplo,n,a,lda,w,work,lwork,rwork,info)
      import :: dp
      implicit none
      
      character(1), intent(in)    :: jobz     ! N/V: if v, calculate eigenvecs.
      character(1), intent(in)    :: uplo     ! U/L: upper/lower triangle.
      integer,      intent(in)    :: n        ! The order of a.
      integer,      intent(in)    :: lda      ! The dimension of a.
      complex(dp),  intent(inout) :: a(lda,*) ! Hermitian matrix.
      real(dp),     intent(out)   :: w(*)     ! Eigenvalues of a.
      complex(dp),  intent(out)   :: work(*)  ! work(1) = optimal lwork.
      integer,      intent(in)    :: lwork    ! The length of work.
      real(dp),     intent(out)   :: rwork(*) ! Working array.
      integer,      intent(out)   :: info     ! 0 on success.
    end subroutine
    
    ! Find the eigenvalues of a general complex matrix.
    subroutine zgeev(jobvl,jobvr,n,a,lda,w,vl,ldvl,vr,ldvr,work,lwork,rwork, &
       &info)
      import :: dp
      implicit none
      
      character(1), intent(in)    :: jobvl      ! N/V: if v, calc. left evecs.
      character(1), intent(in)    :: jobvr      ! N/V: if v, calc. right evecs.
      integer,      intent(in)    :: n          ! The order of a.
      integer,      intent(in)    :: lda        ! The dimension of a.
      complex(dp),  intent(inout) :: a(lda,*)   ! Complex matrix.
      complex(dp),  intent(out)   :: w(*)       ! Eigenvalues of a.
      integer,      intent(in)    :: ldvl       ! The dimension of vl.
      complex(dp),  intent(out)   :: vl(ldvl,*) ! Left eigenvectors of a.
      integer,      intent(in)    :: ldvr       ! The dimension of vr.
      complex(dp),  intent(out)   :: vr(ldvr,*) ! Right eigenvectors of a.
      complex(dp),  intent(out)   :: work(*)    ! work(1) = optimal lwork.
      integer,      intent(in)    :: lwork      ! The length of work.
      real(dp),     intent(out)   :: rwork(*)   ! Working array.
      integer,      intent(out)   :: info       ! 0 on success.
    end subroutine
  end interface
contains

! --------------------------------------------------
! Calculates the eigenvalues and eigenvectors of a real, symmetric matrix.
! --------------------------------------------------
function diagonalise_symmetric_reals(input) result(output)
  implicit none
  
  real(dp), intent(in)                   :: input(:,:)
  type(SymmetricEigenstuff), allocatable :: output(:)
  
  ! dsyev variables.
  integer               :: n       ! The size of the input.
  real(dp), allocatable :: a(:,:)  ! The input. Also eigenvectors.
  real(dp), allocatable :: w(:)    ! The eigenvalues.
  real(dp), allocatable :: work(:) ! Workspace.
  integer               :: lwork   ! Workspace.
  integer               :: info    ! Error code. 0 on success.
  
  ! Temporary variables.
  integer :: i,ialloc
  
  if (size(input,1)/=size(input,2)) then
    call print_line(CODE_ERROR//': Trying to find the eigenvalues of a &
       &non-square matrix.')
    call err()
  endif
  
  n = size(input,1)
  allocate( w(n),        &
          & work(3*n-1), &
          & stat=ialloc); call err(ialloc)
  a = input
  
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
    call print_line(ERROR//'in diagonalisation: dsyev error code: '//info)
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
    call print_line(ERROR//'in diagonalisation: dsyev error code: '//info)
    call err()
  endif
  
  allocate(output(n), stat=ialloc); call err(ialloc)
  do i=1,n
    output(i)%eval = w(i)
    output(i)%evec = a(:,i)
  enddo
end function

function diagonalise_symmetric_RealMatrix(input) result(output)
  implicit none
  
  type(RealMatrix), intent(in)           :: input
  type(SymmetricEigenstuff), allocatable :: output(:)
  
  output = diagonalise_symmetric(dble(input))
end function

! --------------------------------------------------
! Calculates the eigenvalues and eigenvectors of a complex, Hermitian matrix.
! --------------------------------------------------
function diagonalise_hermitian_complexes(input) result(output)
  implicit none
  
  complex(dp), intent(in)                :: input(:,:)
  type(HermitianEigenstuff), allocatable :: output(:)
  
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
  
  n = size(input,1)
  allocate( w(n),         &
          & rwork(3*n-2), &
          & work(3*n-1),  &
          & stat=ialloc); call err(ialloc)
  a = input
  
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
    call print_line(ERROR//'in diagonalisation: zheev error code: '//info)
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
    call print_line(ERROR//'in diagonalisation: zheev error code: '//info)
    call err()
  endif
  
  allocate(output(n), stat=ialloc); call err(ialloc)
  do i=1,n
    output(i)%eval = w(i)
    output(i)%evec = a(:,i)
  enddo
end function

function diagonalise_hermitian_ComplexMatrix(input) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in)        :: input
  type(HermitianEigenstuff), allocatable :: output(:)
  
  output = diagonalise_hermitian(cmplx(input))
end function

! --------------------------------------------------
! Calculates the eigenvalues and eigenvectors of
!    a general square complex matrix.
! --------------------------------------------------
function diagonalise_complex_complexes(input) result(output)
  implicit none
  
  complex(dp), intent(in)              :: input(:,:)
  type(ComplexEigenstuff), allocatable :: output(:)
  
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
  
  n = size(input,1)
  allocate( w(n),       &
          & vl(n,n),    &
          & vr(n,n),    &
          & rwork(2*n), &
          & work(2*n),  &
          & stat=ialloc); call err(ialloc)
  a = input
  
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
    call print_line(ERROR//'in diagonalisation: zgeev error code: '//info)
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
    call print_line(ERROR//'in diagonalisation: zheev error code: '//info)
    call err()
  endif
  
  allocate(output(n), stat=ialloc); call err(ialloc)
  do i=1,n
    output(i)%eval = w(i)
    output(i)%left_evec = vl(:,i)
    output(i)%right_evec = vr(:,i)
  enddo
end function

function diagonalise_complex_ComplexMatrix(input) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in)      :: input
  type(ComplexEigenstuff), allocatable :: output(:)
  
  output = diagonalise_complex(cmplx(input))
end function

! --------------------------------------------------
! Calculates the eigenvalues and eigenvectors of a unitary complex matrix.
! --------------------------------------------------
! Sorts outputs in order of the phase of the eigenvalue.
function diagonalise_unitary_complexes(input,order) result(output)
  use logic_module
  use qr_decomposition_module
  implicit none
  
  complex(dp), intent(in)              :: input(:,:)
  integer,     intent(in)              :: order
  type(UnitaryEigenstuff), allocatable :: output(:)
  
  type(ComplexEigenstuff), allocatable :: estuff(:)
  
  integer,             allocatable :: degenerate_modes(:)
  type(ComplexVector), allocatable :: degenerate_evecs(:)
  
  integer :: i,j,ialloc
  
  ! Diagonalise matrix.
  estuff = diagonalise_complex(input)
  
  ! Convert eigenvalues to exact representation.
  allocate(output(size(estuff)), stat=ialloc); call err(ialloc)
  do i=1,size(estuff)
    output(i)%eval = calculate_phase(estuff(i)%eval, order)
    output(i)%evec = estuff(i)%right_evec
  enddo
  
  ! Sort output in ascending order of phase.
  output = output(sort(output,compare_phases))
  
  ! Orthonormalise degenerate eigenvectors.
  do i=1,size(output)
    degenerate_modes = filter(output%eval==output(i)%eval)
    if (degenerate_modes(1)<i) then
      ! This mode has already been processed.
      cycle
    elseif (size(degenerate_modes)==1) then
      ! This mode is not degenerate with any other mode.
      cycle
    endif
    
    allocate( degenerate_evecs(size(degenerate_modes)), &
            & stat=ialloc); call err(ialloc)
    do j=1,size(degenerate_modes)
      degenerate_evecs(j) = vec(output(degenerate_modes(j))%evec)
    enddo
    degenerate_evecs = orthonormalise(degenerate_evecs)
    do j=1,size(degenerate_modes)
      output(degenerate_modes(j))%evec = cmplx(degenerate_evecs(j))
    enddo
    deallocate(degenerate_evecs, stat=ialloc); call err(ialloc)
  enddo
contains
  ! A lambda for ordering phases.
  ! Captures nothing.
  function compare_phases(this,that) result(output)
    implicit none
    
    class(*), intent(in) :: this
    class(*), intent(in) :: that
    logical              :: output
    
    select type(this); type is(UnitaryEigenstuff)
      select type(that); type is(UnitaryEigenstuff)
        output = this%eval%fraction < that%eval%fraction
      end select
    end select
  end function
end function

function diagonalise_unitary_ComplexMatrix(input,order) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in)      :: input
  integer,             intent(in)      :: order
  type(UnitaryEigenstuff), allocatable :: output(:)
  
  output = diagonalise_unitary(cmplx(input), order)
end function
end module
