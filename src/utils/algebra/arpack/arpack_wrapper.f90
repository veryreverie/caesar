! ======================================================================
! An interface to various ARPACK routines.
! ======================================================================
module arpack_wrapper_module
  use precision_module
  use io_module
  implicit none
  
  private
  
  public :: ARPACK_LINKED
  
  logical, parameter :: ARPACK_LINKED = .true.
contains
subroutine ssaupd_interface()
  implicit none
  
  logical, allocatable :: random_var(:)
  
  call ssaupd()
end subroutine
end module
