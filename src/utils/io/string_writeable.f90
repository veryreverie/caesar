!> Provides the [[StringWriteable(type)]] abstract type.
module caesar_string_writeable_module
  use caesar_foundations_module
  
  use caesar_string_module
  use caesar_print_settings_module
  use caesar_string_array_module
  implicit none
  
  private
  
  public :: StringWriteable
  public :: str
  public :: operator(//)
  public :: join
  public :: print_line
  public :: print_lines
  
  !> An abstract type, which allows extended types to be written to a
  !>    [[String(type)]].
  !> Any type which extends [[StringWriteable(type)]] must overload `%write()`.
  type, abstract, extends(NoDefaultConstructor) :: StringWriteable
  contains
    procedure(write_StringWriteable), deferred :: write
  end type
  
  abstract interface
    !> Converts a type which extends [[StringWriteable(type)]] to
    !>    [[String(type)]].
    recursive function write_StringWriteable(this) result(output)
      import String
      import StringWriteable
      implicit none
      
      class(StringWriteable), intent(in) :: this
      type(String)                       :: output
    end function
  end interface
  
  interface str
    !> Converts a [[StringWriteable(type)]] to [[String(type)]].
    recursive module function str_StringWriteable(this,settings) &
       & result(output) 
      class(StringWriteable), intent(in)           :: this
      !> Options controlling formatting.
      type(PrintSettings),    intent(in), optional :: settings
      type(String)                                 :: output
    end function

    !> Converts a [[StringWriteable(type)]] array to a [[String(type)]] array.
    recursive module function str_StringWriteables_String(this, &
       & separating_line,settings) result(output) 
      class(StringWriteable), intent(in)           :: this(:)
      !> If present, `separating_line` will be included between each
      !>    element of `this`.
      type(String),           intent(in), optional :: separating_line
      !> Options controlling formatting.
      type(PrintSettings),    intent(in), optional :: settings
      type(String), allocatable                    :: output(:)
    end function

    !> Converts a [[StringWriteable(type)]] array to a [[String(type)]] array.
    recursive module function str_StringWriteables_character(this, &
       & separating_line,settings) result(output) 
      class(StringWriteable), intent(in)           :: this(:)
      !> If present, `separating_line` will be included between each
      !>    element of `this`.
      character(*),           intent(in)           :: separating_line
      !> Options controlling formatting.
      type(PrintSettings),    intent(in), optional :: settings
      type(String), allocatable                    :: output(:)
    end function
  end interface
    
  interface operator(//)
    !> Concatenates [[StringWriteable(type)]] and `character(*)`.
    recursive module function concatenate_StringWriteable_character(this, &
       & that) result(output) 
      class(StringWriteable), intent(in) :: this
      character(*),           intent(in) :: that
      type(String)                       :: output
    end function

    !> Concatenates `character(*)` and [[StringWriteable(type)]].
    recursive module function concatenate_character_StringWriteable(this, &
       & that) result(output) 
      character(*),           intent(in) :: this
      class(StringWriteable), intent(in) :: that
      type(String)                       :: output
    end function

    !> Concatenates [[StringWriteable(type)]] and [[String(type)]].
    recursive module function concatenate_StringWriteable_String(this,that) &
       & result(output) 
      class(StringWriteable), intent(in) :: this
      type(String),           intent(in) :: that
      type(String)                       :: output
    end function

    !> Concatenates [[String(type)]] and [[StringWriteable(type)]].
    recursive module function concatenate_String_StringWriteable(this,that) &
       & result(output) 
      type(String),           intent(in) :: this
      class(StringWriteable), intent(in) :: that
      type(String)                       :: output
    end function
  end interface
  
  interface join
    !> Converts a [[StringWriteable(type)]] array to a [[String(type)]] array,
    !>    then concatenates this into a single [[String(type)]].
    recursive module function join_StringWriteable(this,delimiter,settings) &
       & result(output) 
      class(StringWriteable), intent(in)           :: this(:)
      !> `delimiter` will be inserted between each element of `this`.
      !> Defaults to a single space.
      character(*),           intent(in), optional :: delimiter
      !> Options controlling formatting.
      type(PrintSettings),    intent(in), optional :: settings
      type(String)                                 :: output
    end function
  end interface
  
  interface print_line
    !> Converts a [[StringWriteable]] to [[String(type)]], and prints it
    !>    according to the settings of
    !>    [[caesar_intrinsics_module:print_line(interface)]].
    module subroutine print_line_StringWriteable(this,settings) 
      class(StringWriteable), intent(in)           :: this
      !> Options controlling formatting.
      type(PrintSettings),    intent(in), optional :: settings
    end subroutine
  end interface
  
  interface print_lines
    !> Converts a [[StringWriteable]] array to a [[String(type)]] array,
    !>    and prints it according to the settings of
    !>    [[caesar_intrinsics_module:print_lines(interface)]].
    module subroutine print_lines_StringWriteables_String(this, &
       & separating_line,settings) 
      class(StringWriteable), intent(in)           :: this(:)
      !> If present, `separating_line` will be included between each
      !>    element of `this`.
      type(String),           intent(in), optional :: separating_line
      !> Options controlling formatting.
      type(PrintSettings),    intent(in), optional :: settings
    end subroutine
    
    !> Converts a [[StringWriteable]] array to a [[String(type)]] array,
    !>    and prints it according to the settings of
    !>    [[caesar_intrinsics_module:print_lines(interface)]].
    module subroutine print_lines_StringWriteables_character(this, &
       & separating_line,settings) 
      class(StringWriteable), intent(in)           :: this(:)
      !> If present, `separating_line` will be included between each
      !>    element of `this`.
      character(*),           intent(in)           :: separating_line
      !> Options controlling formatting.
      type(PrintSettings),    intent(in), optional :: settings
    end subroutine
  end interface
end module
