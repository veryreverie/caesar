!> Provides the [[StringArray(type)]] type.
module caesar_string_array_module
  use caesar_string_module
  implicit none
  
  private
  
  public :: StringArray
  public :: size
  public :: str
  public :: split_into_sections
  public :: operator(//)
  public :: join
  
  !> An array of type [[String(type)]].
  type :: StringArray
    type(String), allocatable :: strings(:)
  end type
  
  interface StringArray
    !> Constructor for [[StringArray(type)]]. Simply sets `this%strings`.
    module function new_StringArray_Strings(input) result(this) 
      type(String), intent(in) :: input(:)
      type(StringArray)        :: this
    end function
  end interface
  
  interface size
    !> Returns the size of `this%strings`.
    module function size_StringArray(this) result(output) 
      type(StringArray), intent(in) :: this
      integer                       :: output
    end function
  end interface
  
  interface str
    !> Converts a [[StringArray(type)]] to a [[String(type)]] array.
    !> Simply returns `this%strings`.
    module function str_StringArray(this) result(output) 
      type(StringArray), intent(in) :: this
      type(String), allocatable     :: output(:)
    end function

    !> Converts a [[StringArray(type)]] array to a [[String(type)]] array.
    !> Concatenates each `%string` array,
    !>    sepearating arrays by `separating_array` if present.
    module function str_StringArrays_String(this,separating_line) &
       & result(output) 
      type(StringArray), intent(in)           :: this(:)
      type(String),      intent(in), optional :: separating_line
      type(String), allocatable               :: output(:)
    end function

    !> Converts a [[StringArray(type)]] array to a [[String(type)]] array.
    !> Concatenates each `%string` array,
    !>    sepearating arrays by `separating_array` if present.
    module function str_StringArrays_character(this,separating_line) &
       & result(output) 
      type(StringArray), intent(in) :: this(:)
      character(*),      intent(in) :: separating_line
      type(String), allocatable     :: output(:)
    end function
  end interface
  
  interface split_into_sections
    !> Splits a [[String(type)]] array into a [[StringArray(type)]] array,
    !>    splitting by `separating_line`, which defaults to an empty string.
    module function split_into_sections_Strings_String(this,separating_line) &
       & result(output) 
      type(String), intent(in)           :: this(:)
      type(String), intent(in), optional :: separating_line
      type(StringArray), allocatable     :: output(:)
    end function

    !> Splits a [[String(type)]] array into a [[StringArray(type)]] array,
    !>    splitting by `separating_line`, which defaults to an empty string.
    module function split_into_sections_Strings_character(this, &
       & separating_line) result(output) 
      type(String), intent(in)       :: this(:)
      character(*), intent(in)       :: separating_line
      type(StringArray), allocatable :: output(:)
    end function

    !> Splits a [[StringArray(type)]] into a [[StringArray(type)]] array,
    !>    splitting by `separating_line`, which defaults to an empty string.
    module function split_into_sections_StringArray_String(this, &
       & separating_line) result(output) 
      type(StringArray), intent(in)           :: this
      type(String),      intent(in), optional :: separating_line
      type(StringArray), allocatable          :: output(:)
    end function

    !> Splits a [[StringArray(type)]] into a [[StringArray(type)]] array,
    !>    splitting by `separating_line`, which defaults to an empty string.
    module function split_into_sections_StringArray_character(this, &
       & separating_line) result(output) 
      type(StringArray), intent(in)  :: this
      character(*),      intent(in)  :: separating_line
      type(StringArray), allocatable :: output(:)
    end function
  end interface
  
  interface operator(//)
    !> Concatenate two [[StringArray(type)]]s.
    module function concatenate_StringArray_StringArray(this,that) &
       & result(output) 
      type(StringArray), intent(in) :: this
      type(StringArray), intent(in) :: that
      type(StringArray)             :: output
    end function

    !> Concatenate a [[StringArray(type)]] and a [[String(type)]].
    module function concatenate_StringArray_String(this,that) result(output) 
      type(StringArray), intent(in) :: this
      type(String),      intent(in) :: that
      type(StringArray)             :: output
    end function

    !> Concatenate a [[StringArray(type)]] and a `character(*)`.
    module function concatenate_StringArray_character(this,that) &
       & result(output) 
      type(StringArray), intent(in) :: this
      character(*),      intent(in) :: that
      type(StringArray)             :: output
    end function

    !> Concatenate a [[StringArray(type)]] and a [[String(type)]] array.
    module function concatenate_StringArray_Strings(this,that) result(output) 
      type(StringArray), intent(in) :: this
      type(String),      intent(in) :: that(:)
      type(StringArray)             :: output
    end function

    !> Concatenate a [[String(type)]] and a [[StringArray(type)]].
    module function concatenate_String_StringArray(this,that) result(output) 
      type(String),      intent(in) :: this
      type(StringArray), intent(in) :: that
      type(StringArray)             :: output
    end function

    !> Concatenate a `character(*)` and a [[StringArray(type)]].
    module function concatenate_character_StringArray(this,that) &
       & result(output) 
      character(*),      intent(in) :: this
      type(StringArray), intent(in) :: that
      type(StringArray)             :: output
    end function

    !> Concatenate a [[String(type)]] array and a [[StringArray(type)]].
    module function concatenate_Strings_StringArray(this,that) result(output) 
      type(String),      intent(in) :: this(:)
      type(StringArray), intent(in) :: that
      type(StringArray)             :: output
    end function
  end interface
  
  interface join
    !> Concatenate a [[StringArray(type)]] array to a single
    !>    [[StringArray(type)]], separating individual string arrays by
    !>    `separating_line` if present.
    module function join_StringArrays_String(input,separating_line) &
       & result(output) 
      type(StringArray), intent(in)           :: input(:)
      type(String),      intent(in), optional :: separating_line
      type(StringArray)                       :: output
    end function

    !> Concatenate a [[StringArray(type)]] array to a single
    !>    [[StringArray(type)]], separating individual string arrays by
    !>    `separating_line` if present.
    module function join_StringArrays_character(input,separating_line) &
       & result(output) 
      type(StringArray), intent(in) :: input(:)
      character(*),      intent(in) :: separating_line
      type(StringArray)             :: output
    end function
  end interface
end module
