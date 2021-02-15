!> Provides the [[StringsWriteable(type)]] abstract type.
module caesar_strings_writeable_module
  use caesar_foundations_module
  
  use caesar_string_module
  use caesar_print_settings_module
  use caesar_string_array_module
  implicit none
  
  private
  
  public :: StringsWriteable
  public :: str
  public :: print_lines
  
  !> An abstract type, which allows extended types to be written to a
  !>    [[String(type)]] array.
  !> Any type which extends [[StringsWriteable(type)]] must overload %write().
  type, abstract, extends(NoDefaultConstructor) :: StringsWriteable
  contains
    procedure(write_StringsWriteable), deferred :: write
  end type
  
  abstract interface
    !> Converts a type which extends [[StringWriteable(type)]] to
    !>    a [[String(type)]] array.
    recursive function write_StringsWriteable(this) result(output)
      import String
      import StringsWriteable
      implicit none
      
      class(StringsWriteable), intent(in) :: this
      type(String), allocatable           :: output(:)
    end function
  end interface
  
  interface str
    !> Converts a [[StringsWriteable(type)]] to a [[String(type)]] array.
    recursive module function str_StringsWriteable(this,settings) &
       & result(output) 
      class(StringsWriteable), intent(in)           :: this
      !> Options controlling formatting.
      type(PrintSettings),     intent(in), optional :: settings
      type(String), allocatable                     :: output(:)
    end function

    !> Converts a [[StringsWriteable(type)]] array to a [[String(type)]] array.
    recursive module function str_StringsWriteables_String(this, &
       & separating_line,settings) result(output) 
      class(StringsWriteable), intent(in)           :: this(:)
      !> If present, `separating_line` will be included between each
      !>    element of `this`.
      type(String),            intent(in), optional :: separating_line
      !> Options controlling formatting.
      type(PrintSettings),     intent(in), optional :: settings
      type(String), allocatable                     :: output(:)
    end function

    !> Converts a [[StringsWriteable(type)]] array to a [[String(type)]] array.
    recursive module function str_StringsWriteables_character(this, &
       & separating_line,settings) result(output) 
      class(StringsWriteable), intent(in)           :: this(:)
      !> If present, `separating_line` will be included between each
      !>    element of `this`.
      character(*),            intent(in)           :: separating_line
      !> Options controlling formatting.
      type(PrintSettings),     intent(in), optional :: settings
      type(String), allocatable                     :: output(:)
    end function
  end interface
  
  interface print_lines
    !> Converts a [[StringsWriteable]] to a [[String(type)]] array,
    !>    and prints it according to the settings of
    !>    [[caesar_intrinsics_module:print_lines(interface)]].
    module subroutine print_lines_StringsWriteable(this,settings) 
      class(StringsWriteable), intent(in)           :: this
      !> Options controlling formatting.
      type(PrintSettings),     intent(in), optional :: settings
    end subroutine

    !> Converts a [[StringsWriteable]] array to a [[String(type)]] array,
    !>    and prints it according to the settings of
    !>    [[caesar_intrinsics_module:print_lines(interface)]].
    module subroutine print_lines_StringsWriteables_String(this, &
       & separating_line,settings) 
      class(StringsWriteable), intent(in)           :: this(:)
      !> If present, `separating_line` will be included between each
      !>    element of `this`.
      type(String),            intent(in), optional :: separating_line
      !> Options controlling formatting.
      type(PrintSettings),     intent(in), optional :: settings
    end subroutine

    !> Converts a [[StringsWriteable]] array to a [[String(type)]] array,
    !>    and prints it according to the settings of
    !>    [[caesar_intrinsics_module:print_lines(interface)]].
    module subroutine print_lines_StringsWriteables_character(this, &
       & separating_line,settings) 
      class(StringsWriteable), intent(in)           :: this(:)
      !> If present, `separating_line` will be included between each
      !>    element of `this`.
      character(*),            intent(in)           :: separating_line
      !> Options controlling formatting.
      type(PrintSettings),     intent(in), optional :: settings
    end subroutine
  end interface
end module
