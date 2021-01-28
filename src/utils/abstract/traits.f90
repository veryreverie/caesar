! ======================================================================
! Provides base classes which can be extended to give useful traits.
! ======================================================================
module caesar_traits_module
  implicit none
  
  private
  
  public :: NoDefaultConstructor
  
  ! Provides a zero-sized private variable, which prevents any types which
  !    extend NoDefaultConstructor from having default constructors.
  type :: NoDefaultConstructor
    logical, private :: no_default_constructor_(0)
  end type
contains
end module
