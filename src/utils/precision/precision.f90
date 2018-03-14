! ======================================================================
! Specifies the precision of in-built types.
! ======================================================================
module precision_module
  implicit none
  
  private
  
  public :: dp
  
  ! Floating-point double precision.
  integer, parameter :: dp=selected_real_kind(15,300)
end module
