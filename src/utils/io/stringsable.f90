!> Provides the [[Stringsable(type)]] abstract type.
module caesar_stringsable_module
  use caesar_io_basic_module
  
  use caesar_string_array_module
  use caesar_strings_readable_module
  use caesar_strings_writeable_module
  implicit none
  
  private
  
  public :: Stringsable
  
  !> An abstract type, which allows extended types to be written to and read
  !>    from a [[String(type)]] array.
  !> Any type which extends [[Stringsable(type)]] must overload
  !>    `%read(input)` and `%write()`.
  type, abstract, extends(StringsWriteable) :: Stringsable
  contains
    procedure(read_Stringsable), deferred :: read
  end type
  
  abstract interface
    !> Converts a [[String(type)]] array to a type which extends
    !>    [[Stringsable(type)]].
    recursive subroutine read_Stringsable(this,input)
      import String
      import Stringsable
      implicit none
      
      class(Stringsable), intent(out) :: this
      type(String),       intent(in)  :: input(:)
    end subroutine
  end interface
end module
