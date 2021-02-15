!> Provides settings for changing how string conversion and printing behave.
module caesar_print_settings_module
  use caesar_string_module
  implicit none
  
  private
  
  public :: PrintSettings
  public :: set_print_settings
  public :: unset_print_settings
  
  !> A list of settings for string conversion and printing.
  type :: PrintSettings
    !> The number of spaces prepended before each printed line.
    integer :: indent
    !> When a line is too long to fit in the terminal,
    !>    it is broken onto multiple lines.
    !> `overhang` defines the number of spaces (in addition to indent) which
    !>    are prepended before the lines beyond the first.
    integer :: overhang
    !> The number of decimal places stored when converting
    !>    a floating point number to a string.
    integer :: decimal_places
    !> The format used for converting floating point numbers to string.
    !> 'es' for scientific formatting, e.g. '-1.234E+002'.
    !> 'f' for fixed point formatting, e.g. '-123.4'.
    type(String) :: floating_point_format
    !> The number of characters before the decimal point
    !>    which are stored when converting floating point numbers to string.
    integer :: integer_digits
  end type
  
  interface PrintSettings
    !> Constructor for [[PrintSettings(type)]], taking string arguments as
    !>    `character(*)`.
    !> All arguments are optional. Those which are not given will default to
    !>    the those from the [[PrintSettings(type)]] currently in use.
    module function new_PrintSettings_character(indent,overhang, &
       & decimal_places,floating_point_format,integer_digits) result(this)
      integer,      intent(in), optional :: indent
      integer,      intent(in), optional :: overhang
      integer,      intent(in), optional :: decimal_places
      character(*), intent(in), optional :: floating_point_format
      integer,      intent(in), optional :: integer_digits
      type(PrintSettings)                :: this
    end function

    !> Constructor for [[PrintSettings(type)]], taking string arguments as
    !>    `[[String(type)]]`.
    !> All arguments are optional. Those which are not given will default to
    !>    the those from the [[PrintSettings(type)]] currently in use.
    module function new_PrintSettings_String(indent,overhang,decimal_places, &
       & floating_point_format,integer_digits) result(this)
      integer,      intent(in), optional :: indent
      integer,      intent(in), optional :: overhang
      integer,      intent(in), optional :: decimal_places
      type(String), intent(in)           :: floating_point_format
      integer,      intent(in), optional :: integer_digits
      type(PrintSettings)                :: this
    end function
  end interface
  
  interface set_print_settings
    !> Sets a [[PrintSettings(type)]] for use.
    module subroutine set_print_settings_PrintSettings(settings)
      type(PrintSettings), intent(in) :: settings
    end subroutine
    
    !> Constructs a [[PrintSettings(type)]], and sets it for use.
    module subroutine set_print_settings_arguments_character(indent,overhang, &
       & decimal_places,floating_point_format,integer_digits)
      integer,      intent(in), optional :: indent
      integer,      intent(in), optional :: overhang
      integer,      intent(in), optional :: decimal_places
      character(*), intent(in), optional :: floating_point_format
      integer,      intent(in), optional :: integer_digits
    end subroutine
    
    !> Constructs a [[PrintSettings(type)]], and sets it for use.
    module subroutine set_print_settings_arguments_String(indent,overhang, &
       & decimal_places,floating_point_format,integer_digits)
      integer,      intent(in), optional :: indent
      integer,      intent(in), optional :: overhang
      integer,      intent(in), optional :: decimal_places
      type(String), intent(in)           :: floating_point_format
      integer,      intent(in), optional :: integer_digits
    end subroutine
  end interface
  
  interface unset_print_settings
    !> Sets there to be no [[PrintSettings(type)]] in use.
    module subroutine unset_print_settings()
    end subroutine
  end interface
end module
