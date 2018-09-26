! ======================================================================
! A dummy module for when ARPACK is not being linked.
! ======================================================================
module arpack_wrapper_module
  use precision_module
  use io_module
  implicit none
  
  private
  
  public :: ARPACK_LINKED
  
  logical, parameter :: ARPACK_LINKED = .false.
contains
end module
