! ======================================================================
! A test space, for temporary checking of misc. parts of the code.
! ======================================================================
module test_module
  use common_module
  implicit none
  
  private
  
  public :: test
  
contains

! ----------------------------------------------------------------------
! Generates keywords and helptext.
! ----------------------------------------------------------------------
function test() result(output)
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'test'
  output%description = 'Runs temporary code for testing purposes.'
  output%keywords = [KeywordData::]
  output%main_subroutine => test_subroutine
  output%suppress_from_helptext = .true.
end function

! ----------------------------------------------------------------------
! Main function.
! ----------------------------------------------------------------------
subroutine test_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  type(String) :: thing
  
  call print_line(-1e111_dp)
  call print_line(-1e11_dp)
  call print_line(-1e1_dp)
  call print_line(-1e-1_dp)
  call print_line(-1e-11_dp)
  call print_line(-1e-111_dp)
  
  call print_line(cmplx(-1e111_dp,0.0_dp,dp))
  call print_line(cmplx(-1e11_dp,0.0_dp,dp))
  call print_line(cmplx(-1e1_dp,0.0_dp,dp))
  call print_line(cmplx(-1e-1_dp,0.0_dp,dp))
  call print_line(cmplx(-1e-11_dp,0.0_dp,dp))
  call print_line(cmplx(-1e-111_dp,0.0_dp,dp))
  
  call set_print_settings(decimal_places=3)
  call print_line(-1e111_dp)
  call print_line(-1e11_dp)
  call print_line(-1e1_dp)
  call print_line(-1e-1_dp)
  call print_line(-1e-11_dp)
  call print_line(-1e-111_dp)
  
  call print_line(cmplx(-1e111_dp,0.0_dp,dp))
  call print_line(cmplx(-1e11_dp,0.0_dp,dp))
  call print_line(cmplx(-1e1_dp,0.0_dp,dp))
  call print_line(cmplx(-1e-1_dp,0.0_dp,dp))
  call print_line(cmplx(-1e-11_dp,0.0_dp,dp))
  call print_line(cmplx(-1e-111_dp,0.0_dp,dp))
  
  call unset_print_settings()
  call set_print_settings(floating_point_format=str('f'))
  call print_line(-1e1_dp)
  call print_line(-1e-1_dp)
  call print_line(-1e-11_dp)
  call print_line(-1e-111_dp)
  
  call print_line(cmplx(-1e1_dp,0.0_dp,dp))
  call print_line(cmplx(-1e-1_dp,0.0_dp,dp))
  call print_line(cmplx(-1e-11_dp,0.0_dp,dp))
  call print_line(cmplx(-1e-111_dp,0.0_dp,dp))
end subroutine
end module
