! ======================================================================
! Provides settings for changing how print_line and related functions work.
! ======================================================================
module print_settings_submodule
  use string_submodule
  implicit none
  
  private
  
  public :: PrintSettings
  public :: set_print_settings
  public :: unset_print_settings
  
  type :: PrintSettings
    ! The number of spaces before a line starts.
    integer :: indent
    ! The additional number of spaces before the continuation of a line which
    !    is too long to fit on the terminal.
    integer :: overhang
    ! The number of decimal places printed of a floating point number.
    integer :: decimal_places
    ! 'es' for scientific formatting, e.g. '-1.234E+002'.
    ! 'f' for fixed point formatting, e.g. '-123.4'.
    type(String) :: floating_point_format
    ! The number of characters before the decimal point
    !    in fixed-format printing.
    integer :: integer_digits
  end type
  
  interface PrintSettings
    module procedure new_PrintSettings
  end interface
  
  logical             :: USE_CURRENT_PRINT_SETTINGS = .false.
  type(PrintSettings) :: CURRENT_PRINT_SETTINGS
contains

! Constructor.
function new_PrintSettings(indent,overhang,decimal_places, &
   & floating_point_format,integer_digits) result(this)
  implicit none
  
  integer,      intent(in), optional :: indent
  integer,      intent(in), optional :: overhang
  integer,      intent(in), optional :: decimal_places
  type(String), intent(in), optional :: floating_point_format
  integer,      intent(in), optional :: integer_digits
  type(PrintSettings)                :: this
  
  if (.not. USE_CURRENT_PRINT_SETTINGS) then
    ! Default arguments.
    this%indent = 0
    this%overhang = 3
    this%decimal_places = 17
    this%floating_point_format = 'es'
    this%integer_digits = 1
  else
    ! Currently set settings.
    this = CURRENT_PRINT_SETTINGS
  endif
  
  ! Settings specific to this call.
  if (present(indent)) then
    this%indent = indent
  endif
  
  if (present(overhang)) then
    this%overhang = overhang
  endif
  
  if (present(decimal_places)) then
    this%decimal_places = decimal_places
  endif
  
  if (present(floating_point_format)) then
    this%floating_point_format = floating_point_format
  endif
  
  if (present(integer_digits)) then
    this%integer_digits = integer_digits
  endif
end function

subroutine set_print_settings(settings)
  implicit none
  
  type(PrintSettings), intent(in) :: settings
  
  CURRENT_PRINT_SETTINGS = settings
  USE_CURRENT_PRINT_SETTINGS = .true.
end subroutine

subroutine unset_print_settings()
  implicit none
  
  USE_CURRENT_PRINT_SETTINGS = .false.
end subroutine
end module
