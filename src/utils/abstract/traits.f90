!> Provides base classes which can be extended to give useful traits.
module caesar_traits_module
  implicit none
  
  private
  
  public :: NoDefaultConstructor
  
  !> Any type which extends [[NoDefaultConstructor(type)]] will not have the
  !>    default constructor.
  type :: NoDefaultConstructor
    logical, private :: no_default_constructor_(0)
  end type
contains
end module
