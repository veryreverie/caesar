! ======================================================================
! A dummy module for when ARPACK is not being linked.
! ======================================================================
module arpack_wrapper_module
  use precision_module
  use io_module
  implicit none
  
  private
  
  public :: dsaupd
  public :: dseupd
  
  public :: ARPACK_LINKED
  
  logical, parameter :: ARPACK_LINKED = .false.
contains
subroutine dsaupd(ido,bmat,n,which,nev,tol,resid,ncv,v,ldv,iparam,ipntr, &
   & workd,workl,lworkl,info)
  implicit none
  
  integer,      intent(inout) :: ido
  character(1), intent(in)    :: bmat
  integer,      intent(in)    :: n
  character(2), intent(in)    :: which
  integer,      intent(in)    :: nev
  real(dp),     intent(in)    :: tol
  real(dp),     intent(inout) :: resid(*)
  integer,      intent(in)    :: ncv
  integer,      intent(in)    :: ldv
  real(dp),     intent(out)   :: v(ldv,*)
  integer,      intent(inout) :: iparam(*)
  integer,      intent(out)   :: ipntr(*)
  real(dp),     intent(inout) :: workd(*)
  real(dp),     intent(inout) :: workl(*)
  integer,      intent(in)    :: lworkl
  integer,      intent(inout) :: info
  
  call print_line(ERROR//': Cannot perform Lanczos algorithm because Caesar &
     &has not been linked against ARPACK. Please use the CMake flag &
     &-DLINK_TO_ARPACK:LOGICAL=true to link against ARPACK.')
  call err()
  stop
end subroutine

subroutine dseupd(rvec,howmny,select,d,z,ldz,sigma,bmat,n,which,nev,tol, &
   & resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,info)
  implicit none
  
  logical,      intent(in)    :: rvec
  character(1), intent(in)    :: howmny
  logical,      intent(inout) :: select(*)
  real(dp),     intent(out)   :: d(*)
  integer,      intent(in)    :: ldz
  real(dp),     intent(out)   :: z(ldz,*)
  real(dp),     intent(in)    :: sigma
  character(1), intent(in)    :: bmat
  integer,      intent(in)    :: n
  character(2), intent(in)    :: which
  integer,      intent(in)    :: nev
  real(dp),     intent(in)    :: tol
  real(dp),     intent(inout) :: resid(*)
  integer,      intent(in)    :: ncv
  integer,      intent(in)    :: ldv
  real(dp),     intent(out)   :: v(ldv,*)
  integer,      intent(inout) :: iparam(*)
  integer,      intent(out)   :: ipntr(*)
  real(dp),     intent(inout) :: workd(*)
  real(dp),     intent(inout) :: workl(*)
  integer,      intent(in)    :: lworkl
  integer,      intent(inout) :: info
  
  call print_line(ERROR//': Cannot perform Lanczos algorithm because Caesar &
     &has not been linked against ARPACK. Please use the CMake flag &
     &-DLINK_TO_ARPACK:LOGICAL=true to link against ARPACK.')
  call err()
  stop
end subroutine
end module
