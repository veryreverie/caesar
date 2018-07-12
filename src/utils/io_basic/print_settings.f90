! ======================================================================
! Provides settings for changing how print_line and related functions work.
! ======================================================================
module print_settings_submodule
  implicit none
  
  private
  
  public :: PrintSettings
  
  type :: PrintSettings
    ! The number of spaces before a line starts.
    integer :: indent
    ! The additional number of spaces before the continuation of a line which
    !    is too long to fit on the terminal.
    integer :: overhang
    ! The number of decimal places printed of a floating point number.
    integer :: decimal_places
  end type
  
  interface PrintSettings
    module procedure new_PrintSettings
  end interface
contains

! Constructor.
function new_PrintSettings(indent,overhang,decimal_places) result(this)
  implicit none
  
  integer, intent(in), optional :: indent
  integer, intent(in), optional :: overhang
  integer, intent(in), optional :: decimal_places
  type(PrintSettings)           :: this
  
  if (present(indent)) then
    this%indent = indent
  else
    this%indent = 0
  endif
  
  if (present(overhang)) then
    this%overhang = overhang
  else
    this%overhang = 3
  endif
  
  if (present(decimal_places)) then
    this%decimal_places = decimal_places
  else
    this%decimal_places = 17
  endif
end function
end module
