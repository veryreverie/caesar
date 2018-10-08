! ======================================================================
! Provides the Lanczos algorithm, which calculates
!    the lowest n eigenvalues and eigenvectors of a real symmetric matrix.
! ======================================================================
module lanczos_module
  use precision_module
  use io_module
  
  use arpack_wrapper_module
  
  use linear_algebra_module
  use hermitian_eigenstuff_module
  implicit none
  
  private
  
  public :: lanczos
contains

function lanczos(input,no_eigenvalues,eigenvalue_convergence) result(output)
  implicit none
  
  type(RealMatrix), intent(in)           :: input
  integer,          intent(in)           :: no_eigenvalues
  real(dp),         intent(in)           :: eigenvalue_convergence
  type(SymmetricEigenstuff), allocatable :: output(:)
  
  ! ARPACK dsaupd variables.
  integer                   :: ido
  integer                   :: n
  real(dp),     allocatable :: resid(:)
  integer                   :: ncv
  real(dp),     allocatable :: v(:,:)
  integer,      allocatable :: iparam(:)
  integer,      allocatable :: ipntr(:)
  real(dp),     allocatable :: workd(:)
  real(dp),     allocatable :: workl(:)
  integer                   :: lworkl
  integer                   :: info
  
  ! ARPACK dseupd variables.
  logical,  allocatable :: select(:)
  real(dp), allocatable :: d(:)
  real(dp)              :: sigma
  
  ! Temporary variables.
  integer :: ialloc
  
  integer :: i
  
  if (size(input,1)/=size(input,2)) then
    call print_line(CODE_ERROR//': Trying to call Lanczos algorithm on &
       &non-square matrix.')
    call err()
  elseif (no_eigenvalues>=size(input,1)) then
    call print_line(CODE_ERROR//': Trying to call Lanczos algorithm and &
       &requesting at least as many eigenvalues as the dimension of the &
       &matrix.')
    call err()
  endif
  
  ido = 0
  n = size(input,1)
  ncv = min(2*no_eigenvalues,n) ! TODO: try varying this.
  lworkl = ncv*(ncv+8)
  info = 0
  allocate( resid(n),      &
          & v(n,ncv),      &
          & ipntr(11),     &
          & workd(3*n),    &
          & workl(lworkl), &
          & stat=ialloc); call err(ialloc)
  
  iparam = [ 1,       & ! 0 or 1.
           & 0,       & ! Not referenced.
           & huge(0), & ! input: max. iterations. output: iterations used.
           & 1,       & ! Must be 1.
           & 0,       & ! output: num converged evals.
           & 0,       & ! Not referenced.
           & 1,       & ! Mode. 1 for standard Lanczos.
           & 0,       & ! Output: num shifts to provide.
           & 0,       & ! Output: num OP*x operations used.
           & 0,       & ! Output: num B*x operations used.
           & 0        ] ! Output: num steps re-orthogonalisation used.
  
  do
    call dsaupd( ido    = ido,                    &
               & bmat   = 'I',                    &
               & n      = n,                      &
               & which  = 'SA',                   &
               & nev    = no_eigenvalues,         &
               & tol    = eigenvalue_convergence, &
               & resid  = resid,                  &
               & ncv    = ncv,                    &
               & v      = v,                      &
               & ldv    = n,                      &
               & iparam = iparam,                 &
               & ipntr  = ipntr,                  &
               & workd  = workd,                  &
               & workl  = workl,                  &
               & lworkl = lworkl,                 &
               & info   = info                    )
    if (info<0) then
      call print_line(ERROR//': ARPACK dsaupd error: info='//info)
      call err()
    elseif (ido==-1 .or. ido==1) then
      workd(ipntr(2):ipntr(2)+n-1) = dble( input                             &
                                       & * vec(workd(ipntr(1):ipntr(1)+n-1)) )
    else
      exit
    endif
  enddo
  
  allocate( select(ncv),       &
          & d(no_eigenvalues), &
          & stat=ialloc); call err(ialloc)
  call dseupd( rvec   = .true.,                 &
             & howmny = 'A',                    &
             & select = select,                 &
             & d      = d,                      &
             & z      = v,                      &
             & ldz    = n,                      &
             & sigma  = sigma,                  &
             & bmat   = 'I',                    &
             & n      = n,                      &
             & which  = 'SA',                   &
             & nev    = no_eigenvalues,         &
             & tol    = eigenvalue_convergence, &
             & resid  = resid,                  &
             & ncv    = ncv,                    &
             & ldv    = n,                      &
             & v      = v,                      &
             & iparam = iparam,                 &
             & ipntr  = ipntr,                  &
             & workd  = workd,                  &
             & workl  = workl,                  &
             & lworkl = lworkl,                 &
             & info   = info                    )
  if (info/=0) then
    call print_line(ERROR//': ARPACK dseupd error: info='//info)
    call err()
  elseif (iparam(5)<no_eigenvalues) then
    call print_line(ERROR//': the Lanczos algorithm has returned fewer than &
       &the requested number of eigenvalues.')
    call err()
  endif
  
  allocate(output(no_eigenvalues), stat=ialloc); call err(ialloc)
  do i=1,no_eigenvalues
    output(i) = SymmetricEigenstuff(d(i), v(:,i))
  enddo
end function
end module
