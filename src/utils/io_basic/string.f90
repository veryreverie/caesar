!> Provides the [[String(type)]] class.
module caesar_string_module
  use caesar_precision_module
  use caesar_string_base_module
  use caesar_error_module
  implicit none
  
  private
  
  public :: String
  public :: str
  public :: len
  public :: join
  public :: lower_case
  public :: spaces
  public :: split_line
  public :: slice
  public :: trim
  public :: replace
  
  !> A simple heap-allocated String class.
  !> Allows for inhomogeneous character arrays for e.g. storing files.
  type, extends(StringBase) :: String
  contains
    generic,   public             :: operator(//) =>                &
                                   & concatenate_String_character_, &
                                   & concatenate_character_String_, &
                                   & concatenate_String_String_
    procedure, public             :: concatenate_String_character_
    procedure, public, pass(that) :: concatenate_character_String_
    procedure, public             :: concatenate_String_String_
  end type
  
  interface
    !> Concatenate a [[String(type)]] and a `character(*)`.
    module function concatenate_String_character_(this,that) result(output)
      class(String), intent(in) :: this
      character(*),  intent(in) :: that
      type(String)              :: output
    end function

    !> Concatenate a `character(*)` and a [[String(type)]].
    module function concatenate_character_String_(this,that) result(output)
      character(*),  intent(in) :: this
      class(String), intent(in) :: that
      type(String)              :: output
    end function

    !> Concatenate a [[String(type)]] and a [[String(type)]].
    module function concatenate_String_String_(this,that) result(output)
      class(String), intent(in) :: this
      class(String), intent(in) :: that
      type(String)              :: output
    end function
  end interface
  
  interface str
    !> Converts `character(*)` to [[String(type)]].
    impure elemental module function str_character(this) result(output)
      character(*), intent(in) :: this
      type(String)             :: output
    end function
  end interface
  
  interface len
    !> Returns the length of the [[String(type)]].
    !> Equivalent to the `character(*)` `len()` function.
    module function len_String(this) result(output)
      type(String), intent(in) :: this
      integer                  :: output
    end function
  end interface
  
  interface join
    !> Joins a [[String(type)]] array into one [[String(type)]].
    !> If a `delimiter` is given then this will be concatenated between
    !>    each string.
    !> If no `delimiter` is given then a single space will be used as a
    !>    delimiter.
    module function join_String(this,delimiter) result(output)
      type(String), intent(in)           :: this(:)
      character(*), intent(in), optional :: delimiter
      type(String)                       :: output
    end function
  end interface
  
  interface lower_case
    !> Converts a `character(*)` to a lower case [[String(type)]].
    impure elemental module function lower_case_character(input) result(output)
      character(*), intent(in) :: input
      character(len(input))    :: output
    end function

    !> Converts a [[String(type)]] to lower case.
    impure elemental module function lower_case_String(this) result(output)
      type(String), intent(in) :: this
      type(String)             :: output
    end function
  end interface
  
  interface spaces
    !> Returns a `character(*)` consisting of the requested number of spaces.
    module function spaces(no_spaces) result(output)
      integer, intent(in)       :: no_spaces
      character(:), allocatable :: output
    end function
  end interface
  
  interface split_line
    !> Split a `character(*)` into a [[String(type)]] array,
    !>    by a given `delimiter` or `delimiters`.
    !> If no delimiter or delimiters specified, splits by whitespace.
    module function split_line_character(input,delimiter,delimiters) &
       & result(output)
      character(*), intent(in)           :: input
      character(1), intent(in), optional :: delimiter
      character(1), intent(in), optional :: delimiters(:)
      type(String), allocatable          :: output(:)
    end function

    !> Split a [[String(type)]] into a [[String(type)]] array,
    !>    by a given `delimiter` or `delimiters`.
    !> If no delimiter or delimiters specified, splits by whitespace.
    module function split_line_String(input,delimiter,delimiters) &
       & result(output)
      type(String), intent(in)           :: input
      character(1), intent(in), optional :: delimiter
      character(1), intent(in), optional :: delimiters(:)
      type(String), allocatable          :: output(:)
    end function
  end interface
  
  interface first_index
    !> Behaves as the `index` intrinsic, but takes an array of delimiters.
    !> Returns 0 if no delimiter is found.
    module function first_index(input,delimiters) result(output)
      character(*), intent(in) :: input
      character(1), intent(in) :: delimiters(:)
      integer                  :: output
    end function
  end interface
  
  interface slice
    !> Takes a slice of a `character(*)`.
    !> `slice(input,first,last)` is equivalent to `input(first:last)`.
    module function slice_character(input,first,last) result(output)
      character(*), intent(in) :: input
      integer,      intent(in) :: first
      integer,      intent(in) :: last
      type(String)             :: output
    end function

    !> Takes a slice of a [[String(type)]].
    !> `slice(input,first,last)` is equivalent to slicing a `character(*)`
    !>    as `input(first:last)`.
    impure elemental module function slice_String(input,first,last) &
       & result(output)
      type(String), intent(in) :: input
      integer,      intent(in) :: first
      integer,      intent(in) :: last
      type(String)             :: output
    end function
  end interface
  
  interface trim
    !> Removes trailing spaces from a [[String(type)]].
    impure elemental module function trim_String(input) result(output)
      type(String), intent(in) :: input
      type(String)             :: output
    end function
  end interface
  
  interface replace
    !> Finds and replaces one character with another in a `character(*)`.
    module function replace_character(input,from,to) result(output)
      character(*), intent(in) :: input
      character(1), intent(in) :: from
      character(1), intent(in) :: to
      type(String)             :: output
    end function

    !> Finds and replaces one character with another in a [[String(type)]].
    module function replace_String(input,from,to) result(output)
      type(String), intent(in) :: input
      character(*), intent(in) :: from
      character(*), intent(in) :: to
      type(String)             :: output
    end function
  end interface
end module
