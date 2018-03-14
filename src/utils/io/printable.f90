! ======================================================================
! An abstract type, which allows extended types to be printed.
! ======================================================================
! Intended for types which require more than one line to print, compared to
!    stringable types which can be printed on one line.
! Any type which extends Printable can be:
!    - converted to an array of String(:), using this%str().
!    - printed to stdout, using print_line(this).
!    - printed to file, using file%print_line(this).
! See example module below for how to use this module.
module printable_submodule
  use string_submodule
  implicit none
  
  private
  
  public :: Printable
  
  type, abstract :: Printable
  contains
    procedure(str_Printable), deferred :: str
  end type
  
  abstract interface
    function str_Printable(this) result(output)
    import Printable, String
    implicit none
    
    class(Printable), intent(in) :: this
    type(String), allocatable    :: output(:)
    end function
  end interface
end module

! ======================================================================
! An example module for how to extend the Printable type.
! ======================================================================
module printable_example_submodule
  use string_submodule
  use printable_submodule
  
  type, extends(Printable) :: PrintableExample
    integer, allocatable :: contents(:)
  contains
    procedure :: str => str_PrintableExample
  end type
contains

function str_PrintableExample(this) result(output)
  implicit none
  
  class(PrintableExample), intent(in) :: this
  type(String), allocatable           :: output(:)
  
  output = str(this%contents)
end function
end module
