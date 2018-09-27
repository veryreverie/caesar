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
  
  type(RealMatrix)      :: matrix
  integer               :: no_eigenvalues
  real(dp), allocatable :: eigenvalue_convergence
  
  type(SymmetricEigenstuff), allocatable :: estuff(:)
  
  integer :: i
  
  matrix = mat( [ 0.8_dp, 0.5_dp, 0.0_dp,   &
              &   0.5_dp, 0.8_dp, 0.0_dp,   &
              &   0.0_dp, 0.0_dp, 1.0_dp ], &
              & 3,3                         )
  
  no_eigenvalues = 2
  
  eigenvalue_convergence = 1e-8_dp
  
  estuff = diagonalise_symmetric(matrix)
  call print_line('')
  do i=1,no_eigenvalues
    call print_line(estuff(i)%eval)
    call print_line(estuff(i)%evec)
  enddo
  
  estuff = lanczos(matrix,no_eigenvalues,eigenvalue_convergence)
  call print_line('')
  do i=1,no_eigenvalues
    call print_line(estuff(i)%eval)
    call print_line(estuff(i)%evec)
  enddo
end subroutine
end module
