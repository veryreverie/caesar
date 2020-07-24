! ======================================================================
! An abstract type representing something which can be integrated.
! ======================================================================
module integratable_module
  use common_module
  
  use integrated_module
  implicit none
  
  public :: Integratable
  
  type, abstract, extends(NoDefaultConstructor) :: Integratable
  end type
contains
end module
