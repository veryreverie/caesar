! ======================================================================
! An abstract type, which allows extended types to be printed.
! ======================================================================
! Intended for types which require more than one line to print, compared to
!    stringable types which can be printed on one line.
! Any type which extends Printable can be:
!    - converted to String(:), using str(this).
!    - printed to stdout, using print_line(this).
!    - printed to file, using file%print_line(this).
! See example module below for how to use this module.
module printable_submodule
  use string_submodule
  implicit none
  
  private
  
  public :: Printable
  public :: str
  
  type, abstract :: Printable
  contains
    procedure(to_String_Printable), deferred :: to_String
  end type
  
  abstract interface
    recursive function to_String_Printable(this) result(output)
    import String
    import Printable
    implicit none
    
    class(Printable), intent(in) :: this
    type(String), allocatable    :: output(:)
    end function
  end interface
  
  interface str
    module procedure str_Printable
  end interface
contains

recursive function str_Printable(this) result(output)
  implicit none
  
  class(Printable), intent(in) :: this
  type(String), allocatable    :: output(:)
  
  output = this%to_String()
end function
end module

! ======================================================================
! An example module for how to extend the Printable type.
! ======================================================================
module printable_example_submodule
  use string_submodule
  use printable_submodule
  
  type, extends(Printable) :: PrintableExample
    type(String) :: row1
    type(String) :: row2
  contains
    procedure :: to_String => to_String_PrintableExample
  end type
contains

function to_String_PrintableExample(this) result(output)
  implicit none
  
  class(PrintableExample), intent(in) :: this
  type(String), allocatable           :: output(:)
  
  output = [ this%row1, &
           & this%row2 ]
end function
end module
