!> Provides the [[StringsReadable(type)]] abstract type.
module caesar_strings_readable_module
  use caesar_foundations_module
  
  use caesar_string_module
  use caesar_string_array_module
  implicit none
  
  private
  
  public :: StringsReadable
  
  !> An abstract type, which allows extended types to be read from a
  !>    [[String(type)]] array.
  !> Any type which extends [[StringsReadable(type)]] must overload
  !>    `%read(input)`.
  type, abstract, extends(NoDefaultConstructor) :: StringsReadable
  contains
    procedure(read_StringsReadable), deferred :: read
  end type
  
  abstract interface
    !> Converts a [[String(type)]] array to a type which extends
    !>    [[StringsReadable(type)]].
    recursive subroutine read_StringsReadable(this,input)
      import String
      import StringsReadable
      implicit none
      
      class(StringsReadable), intent(out) :: this
      type(String),           intent(in)  :: input(:)
    end subroutine
  end interface
end module
