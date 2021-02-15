! ======================================================================
! Provides low-level utilities, independent of the function of the code.
! Utilities include:
!    - Floating-point precision specification.
!    - I/O handling, including strings and file handling.
!    - Logical algorithms including sort, find, filter, map etc.
!    - Input argument processing.
!    - Algebraic routines, including linear algebra.
! ======================================================================
! This module is simply an interface to various subsidiary modules.
module caesar_utils_module
  use caesar_foundations_module
  use caesar_io_module
  use caesar_abstract_module
  use caesar_arguments_module
  use caesar_algebra_module
  use caesar_random_module
  implicit none
contains
subroutine startup_utils()
  implicit none
  
  call startup_io()
end subroutine
end module
