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
  
  complex(dp) :: one = cmplx(1.0_dp,0.0_dp,dp)
  complex(dp) :: zero = cmplx(0.0_dp,0.0_dp,dp)
  complex(dp) :: small = cmplx(0.01_dp,0.0_dp,dp)
  complex(dp) :: tiny = cmplx(1e-6_dp,0.0_dp,dp)
  
  complex(dp), allocatable :: a(:,:)
  
  type(ComplexVector), allocatable :: invecs1(:)
  type(ComplexVector), allocatable :: invecs2(:)
  type(ComplexVector), allocatable :: outvecs1(:)
  type(ComplexVector), allocatable :: outvecs2(:)
  type(ComplexVector), allocatable :: intersection(:)
  
  type(QRDecomposition) :: qr
  
  integer :: i,j
  
  real(dp) :: l2_err
  
  wd = arguments%value('working_directory')
  
  call print_line('')
  call print_line('Invecs1:')
  invecs1 = [ vec([tiny,zero,zero]), &
            & vec([-tiny,tiny,zero]), &
            & vec([tiny,tiny,zero]), &
            & vec([tiny,zero,-tiny]), &
            & vec([one-tiny,one+tiny,one]),  &
            & vec([zero,tiny,zero]), &
            & vec([tiny,-tiny,tiny]), &
            & vec([tiny,zero,one])]
  do i=1,size(invecs1)
    call print_line(real(invecs1(i)))
  enddo
  
  call print_line('Outvecs1:')
  outvecs1 = orthonormal_basis( invecs1,               &
                              & shortest_valid=0.1_dp, &
                              & longest_invalid=10*abs(tiny))
  do i=1,size(outvecs1)
    call print_line(real(outvecs1(i)))
  enddo
  
  call print_line('')
  call print_line('Invecs2:')
  invecs2 = [ vec([tiny,-tiny,-tiny]),  &
            & vec([-one+tiny,tiny,-one]), &
            & vec([one-tiny,one-tiny,-one-tiny]), &
            & vec([one+tiny,tiny,one]) ]
  do i=1,size(invecs2)
    call print_line(real(invecs2(i)))
  enddo
  
  call print_line('')
  call print_line('Outvecs2:')
  outvecs2 = orthonormal_basis( invecs2,               &
                              & shortest_valid=0.1_dp, &
                              & longest_invalid=10*abs(tiny))
  do i=1,size(outvecs2)
    call print_line(real(outvecs2(i)))
  enddo
  
  call print_line('Intersection')
  intersection = intersection_basis(outvecs1,outvecs2)
  do i=1,size(intersection)
    call print_line(real(intersection(i)))
  enddo
  
end subroutine
end module
