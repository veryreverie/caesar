!> Provides a selection of terminal escape sequences.
module caesar_terminal_module
  implicit none
  
  private
  
  !> The terminal escape character.
  character(1), parameter, public :: ESC = achar(27)
  
  !> Terminal reset string. Removes colour and formatting.
  character(4), parameter, public :: RESET_ESC = ESC//'[0m'
  
  !> Terminal colour string. Colours subsequent output black.
  character(5), parameter, public :: BLACK_ESC = ESC//'[30m'
  !> Terminal colour string. Colours subsequent output red.
  character(5), parameter, public :: RED_ESC = ESC//'[31m'
  !> Terminal colour string. Colours subsequent output green.
  character(5), parameter, public :: GREEN_ESC = ESC//'[32m'
  !> Terminal colour string. Colours subsequent output yellow.
  character(5), parameter, public :: YELLOW_ESC = ESC//'[33m'
  !> Terminal colour string. Colours subsequent output blue.
  character(5), parameter, public :: BLUE_ESC = ESC//'[34m'
  !> Terminal colour string. Colours subsequent output magenta.
  character(5), parameter, public :: MAGENTA_ESC = ESC//'[35m'
  !> Terminal colour string. Colours subsequent output cyan.
  character(5), parameter, public :: CYAN_ESC = ESC//'[36m'
  !> Terminal colour string. Colours subsequent output light gray.
  character(5), parameter, public :: LIGHT_GRAY_ESC = ESC//'[37m'
  !> Terminal colour string. Colours subsequent output dark gray.
  character(5), parameter, public :: DARK_GRAY_ESC = ESC//'[90m'
  !> Terminal colour string. Colours subsequent output light red.
  character(5), parameter, public :: LIGHT_RED_ESC = ESC//'[91m'
  !> Terminal colour string. Colours subsequent output light green.
  character(5), parameter, public :: LIGHT_GREEN_ESC = ESC//'[92m'
  !> Terminal colour string. Colours subsequent output light yellow.
  character(5), parameter, public :: LIGHT_YELLOW_ESC = ESC//'[93m'
  !> Terminal colour string. Colours subsequent output light blue.
  character(5), parameter, public :: LIGHT_BLUE_ESC = ESC//'[94m'
  !> Terminal colour string. Colours subsequent output light magenta.
  character(5), parameter, public :: LIGHT_MAGENTA_ESC = ESC//'[95m'
  !> Terminal colour string. Colours subsequent output light cyan.
  character(5), parameter, public :: LIGHT_CYAN_ESC = ESC//'[96m'
  !> Terminal colour string. Colours subsequent output white.
  character(5), parameter, public :: WHITE_ESC = ESC//'[97m'
end module
