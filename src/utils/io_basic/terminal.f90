! ======================================================================
! A selection of terminal escape sequences.
! ======================================================================
module terminal_module
  implicit none
  
  private
  
  ! The terminal escape character.
  character(1), parameter, public :: ESC = achar(27)
  
  ! Terminal reset string. Removes colour and formatting.
  character(4), parameter, public :: RESET_ESC = ESC//'[0m'
  
  ! Terminal colour strings.
  character(5), parameter, public :: BLACK_ESC = ESC//'[30m'
  character(5), parameter, public :: RED_ESC = ESC//'[31m'
  character(5), parameter, public :: GREEN_ESC = ESC//'[32m'
  character(5), parameter, public :: YELLOW_ESC = ESC//'[33m'
  character(5), parameter, public :: BLUE_ESC = ESC//'[34m'
  character(5), parameter, public :: MAGENTA_ESC = ESC//'[35m'
  character(5), parameter, public :: CYAN_ESC = ESC//'[36m'
  character(5), parameter, public :: LIGHT_GRAY_ESC = ESC//'[37m'
  character(5), parameter, public :: DARK_GRAY_ESC = ESC//'[90m'
  character(5), parameter, public :: LIGHT_RED_ESC = ESC//'[91m'
  character(5), parameter, public :: LIGHT_GREEN_ESC = ESC//'[92m'
  character(5), parameter, public :: LIGHT_YELLOW_ESC = ESC//'[93m'
  character(5), parameter, public :: LIGHT_BLUE_ESC = ESC//'[94m'
  character(5), parameter, public :: LIGHT_MAGENTA_ESC = ESC//'[95m'
  character(5), parameter, public :: LIGHT_CYAN_ESC = ESC//'[96m'
  character(5), parameter, public :: WHITE_ESC = ESC//'[97m'
end module
