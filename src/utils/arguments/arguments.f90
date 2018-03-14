! ======================================================================
! Processes input arguments, and provides help functionality.
! ======================================================================
! This module is simply an interface for the various arguments submodules.
module arguments_module
  use keyword_submodule
  use help_submodule
  use dictionary_submodule
  use process_arguments_submodule
  use caesar_modes_submodule
  implicit none
end module
