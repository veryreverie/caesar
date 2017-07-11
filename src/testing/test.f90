! ======================================================================
! A test space, for temporary checking of misc. parts of the code.
! ======================================================================
module test_module
  use constants_module, only : dp
  use string_module
  use io_module
contains

! ----------------------------------------------------------------------
! Generates keywords and helptext.
! ----------------------------------------------------------------------
function test_keywords() result(keywords)
  use help_module
  implicit none
  
  type(KeywordData) :: keywords(0)
end function

subroutine test(arguments)
  use dictionary_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  type(String) :: wd
  type(String) :: filename
  
  complex(dp)               :: a
  complex(dp)               :: b
  integer                   :: output_file
  type(String), allocatable :: input_file(:)
  
  wd = item(arguments, 'working_directory')
  filename = wd//'/test.test'
  
  output_file = open_write_file(filename)
  
  a = cmplx( 1.0, 2.0,dp)
  call print_line(a)
  call print_line(output_file, a)
  a = cmplx(-1.0, 2.0,dp)
  call print_line(a)
  call print_line(output_file, a)
  a = cmplx( 1.0,-2.0,dp)
  call print_line(a)
  call print_line(output_file, a)
  a = cmplx(-1.0,-2.0,dp)
  call print_line(a)
  call print_line(output_file, a)
  
  close(output_file)
  
  input_file = read_lines(filename)
  b = cmplx(input_file(1))
  call print_line(b)
  b = cmplx(input_file(2))
  call print_line(b)
  b = cmplx(input_file(3))
  call print_line(b)
  b = cmplx(input_file(4))
  call print_line(b)
  
end subroutine
end module
