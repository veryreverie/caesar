!> Provides the [[Stringable(type)]] abstract type.
module caesar_stringable_module
  use caesar_io_basic_module
  
  use caesar_string_array_module
  use caesar_string_readable_module
  use caesar_string_writeable_module
  implicit none
  
  private
  
  public :: Stringable
  
  !> An abstract type which allows extended types to be written to and read
  !>    from a [[String(type)]].
  !> Any type which extends [[StringWriteable(type)]] must overload
  !>    `%read(input)` and `%write()`
  type, abstract, extends(StringWriteable) :: Stringable
  contains
    procedure(read_Stringable), deferred :: read
  end type
  
  abstract interface
    !> Converts a [[String(type)]] to a type which extends
    !>    [[Stringable(type)]].
    recursive subroutine read_Stringable(this,input)
      import String
      import Stringable
      implicit none
      
      class(Stringable), intent(out) :: this
      type(String),      intent(in)  :: input
    end subroutine
  end interface
end module
