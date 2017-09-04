! ======================================================================
! A test space, for temporary checking of misc. parts of the code.
! ======================================================================
module test_module
  use constants_module, only : dp
  use string_module
  use io_module
  
  interface cmplx2
    module procedure cmplx2_IntMatrix
  end interface
contains

! ----------------------------------------------------------------------
! Generates keywords and helptext.
! ----------------------------------------------------------------------
function test_keywords() result(keywords)
  use help_module
  implicit none
  
  type(KeywordData), allocatable :: keywords(:)
  
  keywords = [KeywordData::]
end function

subroutine test(arguments)
  use dictionary_module
  use linear_algebra_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  type(ComplexMatrix) :: test1
  complex(dp), allocatable :: test2(:,:)
  
  integer, allocatable :: int_mat(:,:)
  integer, allocatable :: int_mat_2(:,:)
  real(dp), allocatable :: real_mat(:,:)
  complex(dp), allocatable :: complex_mat(:,:)
  
  integer :: int_1
  integer :: int_2
  real(dp) :: real_1
  complex(dp) :: complex_1
  
  complex(dp) :: a
  
  a = cmplx(0_dp,0_dp,dp)
  
  int_1 = 0
  int_2 = int_1
  real_1 = int_1
  complex_1 = int_1
  
  call print_line('')
  call print_line(int_1)
  call print_line(int_2)
  call print_line(real_1)
  call print_line(complex_1)
  
  allocate(int_mat(2,2))
  int_mat = 0
  int_mat_2 = int_mat
  real_mat = int_mat
  complex_mat = int_mat
  complex_mat = reshape([a,a,a,a],[2,2])
  
  test2 = reshape([a,a,a,a],[2,2])
  test2 = reshape([0,0],[1,2])
  
  test2 = reshape([a,a,a,a,a,a,a,a,a],[3,3])
  
  call print_line('2x2')
  test2 = cmplx(zeroes(2,2))
  call print_line('3x3')
  test2 = cmplx(zeroes(3,3))
  call print_line('Done')
  call print_line(size(test2,1))
  call print_line(size(test2,2))
  test1 = test2
  
end subroutine

pure function cmplx2_IntMatrix(input) result(output)
  use linear_algebra_module
  implicit none
  
  type(IntMatrix), intent(in) :: input
  complex(dp), allocatable    :: output(:,:)
  
  output = int(input)
end function
end module
