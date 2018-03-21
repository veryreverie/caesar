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
module utils_module
  use precision_module
  use traits_module
  use io_module
  use logic_module
  use arguments_module
  use algebra_module
  implicit none
end module
