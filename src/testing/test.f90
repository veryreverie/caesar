! ======================================================================
! A test space, for temporary checking of misc. parts of the code.
! ======================================================================
module test_module
  use common_module
  implicit none
contains

! ----------------------------------------------------------------------
! Generates keywords and helptext.
! ----------------------------------------------------------------------
function test_keywords() result(keywords)
  implicit none
  
  type(KeywordData), allocatable :: keywords(:)
  
  keywords = [KeywordData::]
end function

function test_mode() result(output)
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'test'
  output%description = 'Runs temporary code for testing purposes.'
  output%keywords = test_keywords()
  output%main_subroutine => test
  output%suppress_from_helptext = .true.
end function

! ----------------------------------------------------------------------
! Main function.
! ----------------------------------------------------------------------
subroutine test(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  type(String) :: wd
  
  real(dp) :: one
  real(dp) :: zero
  real(dp) :: small
  
  type(RealVector), allocatable :: inputs(:)
  type(RealVector), allocatable :: outputs(:)
  
  integer :: i
  
  type(RealMatrix) :: matrix
  type(UnitaryEigenstuff), allocatable :: estuff(:)
  
  wd = arguments%value('working_directory')
  
  matrix = dblemat(mat([ 1,0,0, &
                       & 0,0,1, &
                       & 0,1,0],3,3))
  estuff = diagonalise_orthogonal(matrix,2)
  do i=1,size(estuff)
    call print_line(estuff(i)%eval)
    call print_line(estuff(i)%evec)
  enddo
  
  stop
  
  one = 1.0_dp
  zero = 0.0_dp
  small = 1e-10_dp
  
  inputs = [ vec([zero,zero,zero]),     &
           & vec([small,-small,small]), &
           & vec([zero,zero,small]),    &
           & vec([-small,-small,zero]), &
           & vec([zero,one,small]),     &
           & vec([small,one,zero]),     &
           & vec([small,small,zero]),   &
           & vec([one,one,small]) ]
  
  call print_line('Inputs')
  do i=1,size(inputs)
    call print_line(inputs(i))
  enddo
  
  outputs = orthonormal_basis(inputs,0.1_dp,1e-9_dp)
  
  call print_line('Outputs')
  do i=1,size(outputs)
    call print_line(outputs(i))
  enddo
  
  
end subroutine
end module
