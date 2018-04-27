! ======================================================================
! An array of type String.
! Exists to allow heterogeneous 2-D arrays of type String.
! ======================================================================
module string_array_submodule
  use string_submodule
  use printable_submodule
  implicit none
  
  private
  
  public :: StringArray
  
  type, extends(Printable) :: StringArray
    type(String), allocatable :: strings(:)
  contains
    procedure, public :: to_String => to_String_StringArray
  end type
  
  interface size
    module procedure size_StringArray
  end interface
contains

function size_StringArray(this) result(output)
  implicit none
  
  type(StringArray), intent(in) :: this
  integer                       :: output
  
  output = size(this%strings)
end function

function to_String_StringArray(this) result(output)
  implicit none
  
  class(StringArray), intent(in) :: this
  type(String), allocatable      :: output(:)
  
  output = this%strings
end function
end module
