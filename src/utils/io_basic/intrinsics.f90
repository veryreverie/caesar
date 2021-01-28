!> Provides I/O functionality for intrinsic types.
module caesar_intrinsics_module
  use caesar_precision_module
  
  use caesar_string_base_module
  use caesar_string_module
  use caesar_error_module
  use caesar_print_settings_module
  use caesar_print_module
  implicit none
  
  private
  
  public :: str
  public :: left_pad
  public :: pad_int_to_str
  public :: join
  public :: lgcl
  public :: int
  public :: dble
  public :: cmplx
  public :: operator(//)
  public :: print_line
  public :: print_lines
  
  interface str
    !> Converts `logical` to [[String(type)]].
    impure elemental module function str_logical(input,settings) &
       & result(output) 
      implicit none
      
      logical, intent(in)                       :: input
      type(PrintSettings), intent(in), optional :: settings
      type(String)                              :: output
    end function

    !> Converts `integer` to [[String(type)]].
    impure elemental module function str_integer(input,settings) &
       & result(output) 
      implicit none
      
      integer, intent(in)                       :: input
      type(PrintSettings), intent(in), optional :: settings
      type(String)                              :: output
    end function

    !> Converts `real` to [[String(type)]].
    impure elemental module function str_real(input,settings) result(output) 
      implicit none
      
      real(dp), intent(in)                      :: input
      type(PrintSettings), intent(in), optional :: settings
      type(String)                              :: output
    end function

    !> Converts `complex` to [[String(type)]].
    impure elemental module function str_complex(input,settings) &
       & result(output) 
      implicit none
      
      complex(dp), intent(in)                   :: input
      type(PrintSettings), intent(in), optional :: settings
      type(String)                              :: output
    end function

    !> Converts a `logical` array to a [[String(type)]] array.
    !> If `separating_line` is given then it will be inserted between
    !>    each line of the output.
    module function str_logicals_character(input,separating_line,settings) &
       & result(output) 
      implicit none
      
      logical,             intent(in)           :: input(:)
      character(*),        intent(in), optional :: separating_line
      type(PrintSettings), intent(in), optional :: settings
      type(String), allocatable                 :: output(:)
    end function

    !> Converts a `logical` array to a [[String(type)]] array.
    !> If `separating_line` is given then it will be inserted between
    !>    each line of the output.
    module function str_logicals_String(input,separating_line,settings) &
       & result(output) 
      implicit none
      
      logical,             intent(in)           :: input(:)
      type(String),        intent(in)           :: separating_line
      type(PrintSettings), intent(in), optional :: settings
      type(String), allocatable                 :: output(:)
    end function

    !> Converts an `integer` array to a [[String(type)]] array.
    !> If `separating_line` is given then it will be inserted between
    !>    each line of the output.
    module function str_integers_character(input,separating_line,settings) &
       & result(output) 
      implicit none
      
      integer,             intent(in)           :: input(:)
      character(*),        intent(in), optional :: separating_line
      type(PrintSettings), intent(in), optional :: settings
      type(String), allocatable                 :: output(:)
    end function

    !> Converts an `integer` array to a [[String(type)]] array.
    !> If `separating_line` is given then it will be inserted between
    !>    each line of the output.
    module function str_integers_String(input,separating_line,settings) &
       & result(output) 
      implicit none
      
      integer,             intent(in)           :: input(:)
      type(String),        intent(in)           :: separating_line
      type(PrintSettings), intent(in), optional :: settings
      type(String), allocatable                 :: output(:)
    end function

    !> Converts a `real` array to a [[String(type)]] array.
    !> If `separating_line` is given then it will be inserted between
    !>    each line of the output.
    module function str_reals_character(input,separating_line,settings) &
       & result(output) 
      implicit none
      
      real(dp),            intent(in)           :: input(:)
      character(*),        intent(in), optional :: separating_line
      type(PrintSettings), intent(in), optional :: settings
      type(String), allocatable                 :: output(:)
    end function

    !> Converts a `real` array to a [[String(type)]] array.
    !> If `separating_line` is given then it will be inserted between
    !>    each line of the output.
    module function str_reals_String(input,separating_line,settings) &
       & result(output) 
      implicit none
      
      real(dp),            intent(in)           :: input(:)
      type(String),        intent(in)           :: separating_line
      type(PrintSettings), intent(in), optional :: settings
      type(String), allocatable                 :: output(:)
    end function

    !> Converts a `complex` array to a [[String(type)]] array.
    !> If `separating_line` is given then it will be inserted between
    !>    each line of the output.
    module function str_complexes_character(input,separating_line,settings) &
       & result(output) 
      implicit none
      
      complex(dp),         intent(in)           :: input(:)
      character(*),        intent(in), optional :: separating_line
      type(PrintSettings), intent(in), optional :: settings
      type(String), allocatable                 :: output(:)
    end function

    !> Converts a `complex` array to a [[String(type)]] array.
    !> If `separating_line` is given then it will be inserted between
    !>    each line of the output.
    module function str_complexes_String(input,separating_line,settings) &
       & result(output) 
      implicit none
      
      complex(dp),         intent(in)           :: input(:)
      type(String),        intent(in)           :: separating_line
      type(PrintSettings), intent(in), optional :: settings
      type(String), allocatable                 :: output(:)
    end function
  end interface
  
  interface left_pad
    !> Converts an integer to [[String(type)]], then left pads that string with
    !>    `pad_character` until it is as long as `match`.
    !> `pad_character` defaults to '0'.
    impure elemental module function left_pad_character(input,match, &
       & pad_character) result(output) 
      implicit none
      
      integer,      intent(in)           :: input
      character(*), intent(in)           :: match
      character(1), intent(in), optional :: pad_character
      type(String)                       :: output
    end function

    !> Converts an integer to [[String(type)]], then left pads that string with
    !>    `pad_character` until it is as long as `match`.
    !> `pad_character` defaults to '0'.
    impure elemental module function left_pad_String(input,match, &
       & pad_character) result(output) 
      implicit none
      
      integer,      intent(in)           :: input
      type(String), intent(in)           :: match
      character(1), intent(in), optional :: pad_character
      type(String)                       :: output
    end function
  end interface
  
  interface join
    !> Convert a `logical` array to a [[String(type)]] array,
    !>    then concatenate them into a single [[String(type)]].
    !> The strings will be separated by the `delimiter`, which defaults to a
    !>    single space.
    module function join_logicals(this,delimiter,settings) result(output) 
      implicit none
      
      logical,             intent(in)           :: this(:)
      character(*),        intent(in), optional :: delimiter
      type(PrintSettings), intent(in), optional :: settings
      type(String)                              :: output
    end function

    !> Convert an `integer` array to a [[String(type)]] array,
    !>    then concatenate them into a single [[String(type)]].
    !> The strings will be separated by the `delimiter`, which defaults to a
    !>    single space.
    module function join_integers(this,delimiter,pad_sign,settings) &
       & result(output) 
      implicit none
      
      integer,             intent(in)           :: this(:)
      character(*),        intent(in), optional :: delimiter
      !> If `pad_sign` is true then positive integers will be padded with
      !>    a single space, to align with negative integers.
      !> Defaults to true.
      logical,             intent(in), optional :: pad_sign
      type(PrintSettings), intent(in), optional :: settings
      type(String)                              :: output
    end function

    !> Convert a `real` array to a [[String(type)]] array,
    !>    then concatenate them into a single [[String(type)]].
    !> The strings will be separated by the `delimiter`, which defaults to a
    !>    single space.
    module function join_reals(this,delimiter,settings) result(output) 
      implicit none
      
      real(dp),            intent(in)           :: this(:)
      character(*),        intent(in), optional :: delimiter
      type(PrintSettings), intent(in), optional :: settings
      type(String)                              :: output
    end function

    !> Convert a `complex` array to a [[String(type)]] array,
    !>    then concatenate them into a single [[String(type)]].
    !> The strings will be separated by the `delimiter`, which defaults to a
    !>    single space.
    module function join_complexes(this,delimiter,settings) result(output) 
      implicit none
      
      complex(dp),         intent(in)           :: this(:)
      character(*),        intent(in), optional :: delimiter
      type(PrintSettings), intent(in), optional :: settings
      type(String)                              :: output
    end function
  end interface
  
  interface lgcl
    !> Converts from `character(*)` to `logical`.
    module function lgcl_character(this) result(output) 
      implicit none
      
      character(*), intent(in) :: this
      logical                  :: output
    end function

    !> Converts from [[String(type)]] to `logical`.
    impure elemental module function lgcl_String(this) result(output) 
      implicit none
      
      type(String), intent(in) :: this
      logical                  :: output
    end function
  end interface
  
  interface int
    !> Converts from `character(*)` to `integer`.
    module function int_character(this) result(output) 
      implicit none
      
      character(*), intent(in) :: this
      integer                  :: output
    end function

    !> Converts from [[String(type)]] to `integer`.
    impure elemental module function int_String(this) result(output) 
      implicit none
      
      type(String), intent(in) :: this
      integer                  :: output
    end function
  end interface
  
  interface dble
    !> Converts from `character(*)` to `real`.
    !> Can also convert fractions as strings to `real`,
    !>    e.g. `dble('4/5')` returns `0.8_dp`.
    recursive module function dble_character(this) result(output) 
      implicit none
      
      character(*), intent(in) :: this
      real(dp)                 :: output
    end function

    !> Converts from [[String(type)]] to `real`.
    !> Can also convert fractions as strings to real(dp), e.g. '4/5' -> 0.8_dp.
    impure elemental module function dble_String(this) result(output) 
      implicit none
      
      type(String), intent(in) :: this
      real(dp)                 :: output
    end function
  end interface
  
  interface cmplx
    !> Converts from `character(*)` to `complex`.
    module function cmplx_character(this) result(output) 
      implicit none
      
      character(*), intent(in) :: this
      complex(dp)              :: output
    end function

    !> Converts from [[String(type)]] to `complex`.
    impure elemental module function cmplx_String(this) result(output) 
      implicit none
      
      type(String), intent(in) :: this
      complex(dp)              :: output
    end function
  end interface

  interface operator(//)
    !> Concatenate a [[String(type)]] and a `logical`.
    module function concatenate_String_logical(this,that) result(output) 
      implicit none
      
      class(String), intent(in) :: this
      logical,       intent(in) :: that
      type(String)              :: output
    end function

    !> Concatenate a `logical` and a [[String(type)]].
    module function concatenate_logical_String(this,that) result(output) 
      implicit none
      
      logical,       intent(in) :: this
      class(String), intent(in) :: that
      type(String)              :: output
    end function

    !> Concatenate a [[String(type)]] and an `integer`.
    module function concatenate_String_integer(this,that) result(output) 
      implicit none
      
      class(String), intent(in) :: this
      integer,       intent(in) :: that
      type(String)              :: output
    end function

    !> Concatenate an `integer` and a [[String(type)]].
    module function concatenate_integer_String(this,that) result(output) 
      implicit none
      
      integer,       intent(in) :: this
      class(String), intent(in) :: that
      type(String)              :: output
    end function

    !> Concatenate a [[String(type)]] and a `real`.
    module function concatenate_String_real(this,that) result(output) 
      implicit none
      
      class(String), intent(in) :: this
      real(dp),      intent(in) :: that
      type(String)              :: output
    end function

    !> Concatenate a `real` and a [[String(type)]].
    module function concatenate_real_String(this,that) result(output) 
      implicit none
      
      real(dp),      intent(in) :: this
      class(String), intent(in) :: that
      type(String)              :: output
    end function

    !> Concatenate a [[String(type)]] and a `complex`.
    module function concatenate_String_complex(this,that) result(output) 
      implicit none
      
      class(String), intent(in) :: this
      complex(dp),   intent(in) :: that
      type(String)              :: output
    end function

    !> Concatenate a `complex` and a [[String(type)]].
    module function concatenate_complex_String(this,that) result(output) 
      implicit none
      
      complex(dp),   intent(in) :: this
      class(String), intent(in) :: that
      type(String)              :: output
    end function

    !> Concatenate a [[String(type)]] and a `logical` array.
    module function concatenate_String_logicals(this,that) result(output) 
      implicit none
      
      class(String), intent(in) :: this
      logical,       intent(in) :: that(:)
      type(String)              :: output
    end function

    !> Concatenate a `logical` array and a [[String(type)]].
    module function concatenate_logicals_String(this,that) result(output) 
      implicit none
      
      logical,       intent(in) :: this(:)
      class(String), intent(in) :: that
      type(String)              :: output
    end function

    !> Concatenate a [[String(type)]] and an `integer` array.
    module function concatenate_String_integers(this,that) result(output) 
      implicit none
      
      class(String), intent(in) :: this
      integer,       intent(in) :: that(:)
      type(String)              :: output
    end function

    !> Concatenate an `integer` array and a [[String(type)]].
    module function concatenate_integers_String(this,that) result(output) 
      implicit none
      
      integer,       intent(in) :: this(:)
      class(String), intent(in) :: that
      type(String)              :: output
    end function

    !> Concatenate a [[String(type)]] and a `real` array.
    module function concatenate_String_reals(this,that) result(output) 
      implicit none
      
      class(String), intent(in) :: this
      real(dp),      intent(in) :: that(:)
      type(String)              :: output
    end function

    !> Concatenate a `real` array and a [[String(type)]].
    module function concatenate_reals_String(this,that) result(output) 
      implicit none
      
      real(dp),      intent(in) :: this(:)
      class(String), intent(in) :: that
      type(String)              :: output
    end function

    !> Concatenate a [[String(type)]] and a `complex` array.
    module function concatenate_String_complexes(this,that) result(output) 
      implicit none
      
      class(String), intent(in) :: this
      complex(dp),   intent(in) :: that(:)
      type(String)              :: output
    end function

    !> Concatenate a `complex` array and a [[String(type)]].
    module function concatenate_complexes_String(this,that) result(output) 
      implicit none
      
      complex(dp),   intent(in) :: this(:)
      class(String), intent(in) :: that
      type(String)              :: output
    end function

    !> Concatenate a `character(*)` and a `logical`.
    module function concatenate_character_logical(this,that) result(output) 
      implicit none
      
      character(*), intent(in) :: this
      logical,      intent(in) :: that
      type(String)             :: output
    end function

    !> Concatenate a `logical` and a `character(*)`.
    module function concatenate_logical_character(this,that) result(output) 
      implicit none
      
      logical,      intent(in) :: this
      character(*), intent(in) :: that
      type(String)             :: output
    end function

    !> Concatenate a `character(*)` and an `integer`.
    module function concatenate_character_integer(this,that) result(output) 
      implicit none
      
      character(*), intent(in) :: this
      integer,      intent(in) :: that
      type(String)             :: output
    end function

    !> Concatenate an `integer` and a `character(*)`.
    module function concatenate_integer_character(this,that) result(output) 
      implicit none
      
      integer,      intent(in) :: this
      character(*), intent(in) :: that
      type(String)             :: output
    end function

    !> Concatenate a `character(*)` and a `real`.
    module function concatenate_character_real(this,that) result(output) 
      implicit none
      
      character(*), intent(in) :: this
      real(dp),     intent(in) :: that
      type(String)             :: output
    end function

    !> Concatenate a `real` and a `character(*)`.
    module function concatenate_real_character(this,that) result(output) 
      implicit none
      
      real(dp),     intent(in) :: this
      character(*), intent(in) :: that
      type(String)             :: output
    end function

    !> Concatenate a `character(*)` and a `complex`.
    module function concatenate_character_complex(this,that) result(output) 
      implicit none
      
      character(*), intent(in) :: this
      complex(dp),  intent(in) :: that
      type(String)             :: output
    end function

    !> Concatenate a `complex` and a `character(*)`.
    module function concatenate_complex_character(this,that) result(output) 
      implicit none
      
      complex(dp),  intent(in) :: this
      character(*), intent(in) :: that
      type(String)             :: output
    end function

    !> Concatenate a `character(*)` and a `logical` array.
    module function concatenate_character_logicals(this,that) result(output) 
      implicit none
      
      character(*), intent(in) :: this
      logical,      intent(in) :: that(:)
      type(String)             :: output
    end function

    !> Concatenate a `logical` array and a `character(*)`.
    module function concatenate_logicals_character(this,that) result(output) 
      implicit none
      
      logical,      intent(in) :: this(:)
      character(*), intent(in) :: that
      type(String)             :: output
    end function

    !> Concatenate a `character(*)` and an `integer` array.
    module function concatenate_character_integers(this,that) result(output) 
      implicit none
      
      character(*), intent(in) :: this
      integer,      intent(in) :: that(:)
      type(String)             :: output
    end function

    !> Concatenate an `integer` array and a `character(*)`.
    module function concatenate_integers_character(this,that) result(output) 
      implicit none
      
      integer,      intent(in) :: this(:)
      character(*), intent(in) :: that
      type(String)             :: output
    end function

    !> Concatenate a `character(*)` and a `real` array.
    module function concatenate_character_reals(this,that) result(output) 
      implicit none
      
      character(*), intent(in) :: this
      real(dp),     intent(in) :: that(:)
      type(String)             :: output
    end function

    !> Concatenate a `real` array and a `character(*)`.
    module function concatenate_reals_character(this,that) result(output) 
      implicit none
      
      real(dp),     intent(in) :: this(:)
      character(*), intent(in) :: that
      type(String)             :: output
    end function

    !> Concatenate a `character(*)` and a `complex` array.
    module function concatenate_character_complexes(this,that) result(output) 
      implicit none
      
      character(*), intent(in) :: this
      complex(dp),  intent(in) :: that(:)
      type(String)             :: output
    end function

    !> Concatenate a `complex` array and a `character(*)`.
    module function concatenate_complexes_character(this,that) result(output) 
      implicit none
      
      complex(dp),  intent(in) :: this(:)
      character(*), intent(in) :: that
      type(String)             :: output
    end function
  end interface
  
  interface print_line
    !> Prints a `logical` to the terminal.
    module subroutine print_line_logical(this,settings) 
      implicit none
      
      logical,             intent(in)           :: this
      type(PrintSettings), intent(in), optional :: settings
    end subroutine

    !> Prints an `integer` to the terminal.
    module subroutine print_line_integer(this,settings) 
      implicit none
      
      integer,             intent(in)           :: this
      type(PrintSettings), intent(in), optional :: settings
    end subroutine

    !> Prints a `real` to the terminal.
    module subroutine print_line_real(this,settings) 
      implicit none
      
      real(dp),            intent(in)           :: this
      type(PrintSettings), intent(in), optional :: settings
    end subroutine

    !> Prints a `complex` to the terminal.
    module subroutine print_line_complex(this,settings) 
      implicit none
      
      complex(dp),         intent(in)           :: this
      type(PrintSettings), intent(in), optional :: settings
    end subroutine

    !> Prints a `logical` array to the terminal, on one line.
    module subroutine print_line_logicals(this,settings) 
      implicit none
      
      logical,             intent(in)           :: this(:)
      type(PrintSettings), intent(in), optional :: settings
    end subroutine

    !> Prints an `integer` array to the terminal, on one line.
    module subroutine print_line_integers(this,settings) 
      implicit none
      
      integer,             intent(in)           :: this(:)
      type(PrintSettings), intent(in), optional :: settings
    end subroutine

    !> Prints a `real` array to the terminal, on one line.
    module subroutine print_line_reals(this,settings) 
      implicit none
      
      real(dp),            intent(in)           :: this(:)
      type(PrintSettings), intent(in), optional :: settings
    end subroutine

    !> Prints a `complex` array to the terminal, on one line.
    module subroutine print_line_complexes(this,settings) 
      implicit none
      
      complex(dp),         intent(in)           :: this(:)
      type(PrintSettings), intent(in), optional :: settings
    end subroutine
  end interface
  
  interface print_lines
    !> Prints a `logical` array to the terminal,
    !>    with one line per array element.
    module subroutine print_lines_logicals_character(this,separating_line, &
       & settings) 
      implicit none
      
      logical,             intent(in)           :: this(:)
      character(*),        intent(in), optional :: separating_line
      type(PrintSettings), intent(in), optional :: settings
    end subroutine

    !> Prints a `logical` array to the terminal,
    !>    with one line per array element.
    module subroutine print_lines_logicals_String(this,separating_line, &
       & settings) 
      implicit none
      
      logical,             intent(in)           :: this(:)
      type(String),        intent(in)           :: separating_line
      type(PrintSettings), intent(in), optional :: settings
    end subroutine

    !> Prints an `integer` array to the terminal,
    !>    with one line per array element.
    module subroutine print_lines_integers_character(this,separating_line, &
       & settings) 
      implicit none
      
      integer,             intent(in)           :: this(:)
      character(*),        intent(in), optional :: separating_line
      type(PrintSettings), intent(in), optional :: settings
    end subroutine

    !> Prints an `integer` array to the terminal,
    !>    with one line per array element.
    module subroutine print_lines_integers_String(this,separating_line, &
       & settings) 
      implicit none
      
      integer,             intent(in)           :: this(:)
      type(String),        intent(in)           :: separating_line
      type(PrintSettings), intent(in), optional :: settings
    end subroutine

    !> Prints a `real` array to the terminal,
    !>    with one line per array element.
    module subroutine print_lines_reals_character(this,separating_line, &
       & settings) 
      implicit none
      
      real(dp),            intent(in)           :: this(:)
      character(*),        intent(in), optional :: separating_line
      type(PrintSettings), intent(in), optional :: settings
    end subroutine

    !> Prints a `real` array to the terminal,
    !>    with one line per array element.
    module subroutine print_lines_reals_String(this,separating_line,settings) 
      implicit none
      
      real(dp),            intent(in)           :: this(:)
      type(String),        intent(in)           :: separating_line
      type(PrintSettings), intent(in), optional :: settings
    end subroutine

    !> Prints a `complex` array to the terminal,
    !>    with one line per array element.
    module subroutine print_lines_complexes_character(this,separating_line, &
       & settings) 
      implicit none
      
      complex(dp),         intent(in)           :: this(:)
      character(*),        intent(in), optional :: separating_line
      type(PrintSettings), intent(in), optional :: settings
    end subroutine

    !> Prints a `complex` array to the terminal,
    !>    with one line per array element.
    module subroutine print_lines_complexes_String(this,separating_line, &
       & settings) 
      implicit none
      
      complex(dp),         intent(in)           :: this(:)
      type(String),        intent(in)           :: separating_line
      type(PrintSettings), intent(in), optional :: settings
    end subroutine
  end interface
  
  interface pad_int_to_str
    !> Converts an integer to a String, and then pads the string so that
    !>   positive and negative integers end up being the same length.
    impure elemental module function pad_int_to_str(this,settings) &
       & result(output) 
      implicit none
      
      integer,             intent(in)           :: this
      type(PrintSettings), intent(in), optional :: settings
      type(String)                              :: output
    end function
  end interface
end module
