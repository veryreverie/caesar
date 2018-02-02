! ======================================================================
! A selection of terminal escape sequences.
! ======================================================================
module terminal_module
  implicit none
  
  ! The terminal escape character.
  character(1), parameter :: ESC = achar(27)
  
  ! Terminal reset string. Removes colour and formatting.
  character(4), parameter :: RESET_ESC = ESC//'[0m'
  
  ! Terminal colour strings.
  character(5), parameter :: BLACK_ESC = ESC//'[30m'
  character(5), parameter :: RED_ESC = ESC//'[31m'
  character(5), parameter :: GREEN_ESC = ESC//'[32m'
  character(5), parameter :: YELLOW_ESC = ESC//'[33m'
  character(5), parameter :: BLUE_ESC = ESC//'[34m'
  character(5), parameter :: MAGENTA_ESC = ESC//'[35m'
  character(5), parameter :: CYAN_ESC = ESC//'[36m'
  character(5), parameter :: LIGHT_GRAY_ESC = ESC//'[37m'
  character(5), parameter :: DARK_GRAY_ESC = ESC//'[90m'
  character(5), parameter :: LIGHT_RED_ESC = ESC//'[91m'
  character(5), parameter :: LIGHT_GREEN_ESC = ESC//'[92m'
  character(5), parameter :: LIGHT_YELLOW_ESC = ESC//'[93m'
  character(5), parameter :: LIGHT_BLUE_ESC = ESC//'[94m'
  character(5), parameter :: LIGHT_MAGENTA_ESC = ESC//'[95m'
  character(5), parameter :: LIGHT_CYAN_ESC = ESC//'[96m'
  character(5), parameter :: WHITE_ESC = ESC//'[97m'
end module
