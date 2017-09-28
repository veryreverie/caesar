! ======================================================================
! An abstract type, which allows extended types to be turned into strings.
! ======================================================================
! Any type which extends Stringable can be:
!    - converted to String, using string=this or str(this).
!    - concatenated, using string//this or character//this.
!    - printed, using print_line(this)
! See example module below for how to use this module.
module stringable_module
  use string_module
  
  type, abstract :: Stringable
  contains
    generic :: assignment(=) => assign_String
    procedure(assign_String_Stringable), deferred, pass(that) :: assign_String
    
    generic :: operator(//) => concatenate_Stringable_character, &
                             & concatenate_character_Stringable, &
                             & concatenate_Stringable_String,    &
                             & concatenate_String_Stringable
    procedure             ::   concatenate_Stringable_character
    procedure, pass(that) ::   concatenate_character_Stringable
    procedure             ::   concatenate_Stringable_String
    procedure, pass(that) ::   concatenate_String_Stringable
  end type
  
  abstract interface
    subroutine assign_String_Stringable(this,that)
      import Stringable, String
      implicit none
      
      type(String),      intent(inout) :: this
      class(Stringable), intent(in)    :: that
    end subroutine
  end interface
  
  interface str
    module procedure str_Stringable
  end interface
contains

function str_Stringable(this) result(output)
  implicit none
  
  class(Stringable), intent(in) :: this
  type(String)                  :: output
  
  output = this
end function

function concatenate_Stringable_character(this,that) result(output)
  implicit none
  
  class(Stringable), intent(in) :: this
  character(*),      intent(in) :: that
  type(String)                  :: output
  
  output = str(this)//that
end function

function concatenate_character_Stringable(this,that) result(output)
  implicit none
  
  character(*),      intent(in) :: this
  class(Stringable), intent(in) :: that
  type(String)                  :: output
  
  output = this//str(that)
end function

function concatenate_Stringable_String(this,that) result(output)
  implicit none
  
  class(Stringable), intent(in) :: this
  type(String),      intent(in) :: that
  type(String)                  :: output
  
  output = str(this)//that
end function

function concatenate_String_Stringable(this,that) result(output)
  implicit none
  
  type(String),      intent(in) :: this
  class(Stringable), intent(in) :: that
  type(String)                  :: output
  
  output = this//str(that)
end function
end module

! ======================================================================
! An example module for how to extend the Stringable type.
! ======================================================================
module stringable_example_module
  use string_module
  use stringable_module
  
  type, extends(Stringable) :: StringableExample
    integer :: contents
  contains
    procedure, pass(that) :: assign_String => assign_String_StringableExample
  end type
contains

subroutine assign_String_StringableExample(this,that)
  implicit none
  
  type(String),             intent(inout) :: this
  class(StringableExample), intent(in)    :: that
  
  this = that%contents
end subroutine
end module
