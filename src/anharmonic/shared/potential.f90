! ======================================================================
! The abstract potential type.
! ======================================================================
module potential_module
  use common_module
  implicit none
  
  private
  
  public :: Potential
  
  type, abstract, extends(Stringsable) :: Potential
  end type
contains
end module
