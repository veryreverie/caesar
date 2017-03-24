! Provides err(), which aborts to give a stacktrace, and
!    arr(a), which aborts if a=.false. or a/=0.
! This version only works with gfortran.
module err_module
  implicit none
  
  interface err
    module procedure err_none
    module procedure err_logical
    module procedure err_integer
  end interface
contains

subroutine err_none()
  implicit none
  
  call abort
end subroutine

subroutine err_logical(this)
  implicit none
  
  logical, intent(in) :: this
  
  if (.not. this) then
    call err()
  endif
end subroutine

subroutine err_integer(this)
  implicit none
  
  integer, intent(in) :: this
  
  if (this/=0) then
    call err()
  endif
end subroutine
end module
