! ======================================================================
! A dummy module for when QUIP is not being linked against.
! ======================================================================
module quip_wrapper_module
  use constants_module, only : dp
  use string_module
  use io_module
  implicit none
  
  type :: QuipResult
    real(dp)              :: energy
    real(dp), allocatable :: forces(:,:)
    real(dp)              :: virial(3,3)
  end type
contains

function call_quip(lattice,atomic_nos,positions,dir,seedname) result(output)
  implicit none
  
  real(dp),     intent(in) :: lattice(3,3)
  integer,      intent(in) :: atomic_nos(:)
  real(dp),     intent(in) :: positions(:,:)
  type(String), intent(in) :: dir
  type(String), intent(in) :: seedname
  type(QuipResult)         :: output
  
  call print_line(ERROR//': Caesar has not been linked against QUIP.')
  stop
end function
end module
