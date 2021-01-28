! ======================================================================
! An abstract type representing something which can be integrated.
! ======================================================================
module caesar_integratable_module
  use caesar_common_module
  
  use caesar_integrated_module
  implicit none
  
  public :: Integratable
  
  type, abstract, extends(NoDefaultConstructor) :: Integratable
  end type
contains
end module
