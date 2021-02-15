!> Provides the [[StringReadable(type)]] abstract type.
module caesar_string_readable_module
  use caesar_foundations_module
  
  use caesar_string_module
  use caesar_string_array_module
  implicit none
  
  private
  
  public :: StringReadable
  
  !> An abstract type, which allows extended types to be read from a
  !>    [[String(type)]].
  !> Any type which extends [[StringReadable(type)]] must overload
  !>    `%read(input)`.
  type, abstract, extends(NoDefaultConstructor) :: StringReadable
  contains
    procedure(read_StringReadable), deferred :: read
  end type
  
  abstract interface
    !> Converts a [[String(type)]] to a type which extends
    !>    [[StringReadable(type)]].
    recursive subroutine read_StringReadable(this,input)
      import String
      import StringReadable
      implicit none
      
      class(StringReadable), intent(out) :: this
      type(String),          intent(in)  :: input
    end subroutine
  end interface
end module
