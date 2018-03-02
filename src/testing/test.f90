! ======================================================================
! A test space, for temporary checking of misc. parts of the code.
! ======================================================================
module test_module
  use constants_module, only : dp
  use string_module
  use io_module
  implicit none
contains

! ----------------------------------------------------------------------
! Generates keywords and helptext.
! ----------------------------------------------------------------------
function test_keywords() result(keywords)
  use keyword_module
  implicit none
  
  type(KeywordData), allocatable :: keywords(:)
  
  keywords = [KeywordData::]
end function

function test_mode() result(output)
  use caesar_modes_module
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
  use dictionary_module
  use logic_module
  use logic_example_module
  use qr_decomposition_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  type(String) :: wd
  
  complex(dp) :: one = cmplx(1.0_dp,0.0_dp,dp)
  complex(dp) :: zero = cmplx(0.0_dp,0.0_dp,dp)
  complex(dp) :: small = cmplx(0.01_dp,0.0_dp,dp)
  
  complex(dp), allocatable :: a(:,:)
  
  type(ComplexVector), allocatable :: invecs(:)
  type(ComplexVector), allocatable :: outvecs(:)
  
  type(QRDecomposition) :: qr
  
  integer :: i
  
  integer, allocatable :: list1(:)
  integer, allocatable :: list2(:)
  
  wd = arguments%value('working_directory')
  
  call logic_example()
  
  stop
  
  a = cmplx(mat([one,one,one,zero],2,2))
  qr = qr_decomposition(a)
  call print_line('')
  call print_line('a:')
  call print_line(mat(real(a)))
  call print_line('q:')
  call print_line(mat(real(qr%q)))
  call print_line('r:')
  call print_line(mat(real(qr%r)))
  call print_line('qr:')
  call print_line(real(mat(qr%q)*mat(qr%r)))
  
  a = cmplx(mat([one,one,zero,zero,one,zero],2,3))
  qr = qr_decomposition(a)
  call print_line('')
  call print_line('a:')
  call print_line(mat(real(a)))
  call print_line('q:')
  call print_line(mat(real(qr%q)))
  call print_line('r:')
  call print_line(mat(real(qr%r)))
  call print_line('qr:')
  call print_line(real(mat(qr%q)*mat(qr%r)))
  
  a = cmplx(mat([one,one,zero,zero,one,zero],3,2))
  qr = qr_decomposition(a)
  call print_line('')
  call print_line('a:')
  call print_line(mat(real(a)))
  call print_line('q:')
  call print_line(mat(real(qr%q)))
  call print_line('r:')
  call print_line(mat(real(qr%r)))
  call print_line('qr:')
  call print_line(real(mat(qr%q)*mat(qr%r)))
  
  a = cmplx(mat([ small,small,zero,   &
                & zero,small,zero,  &
                & small,zero,zero,  &
                & zero,zero,zero, &
                & one,one,zero    &
                & ],3,5))
  qr = qr_decomposition(a)
  call print_line('')
  call print_line('a:')
  call print_line(mat(real(a)))
  call print_line('q:')
  call print_line(mat(real(qr%q)))
  call print_line('r:')
  call print_line(mat(real(qr%r)))
  call print_line('qr:')
  call print_line(real(mat(qr%q)*mat(qr%r)))
  
  invecs = [vec([one,zero,zero]), vec([one,one,zero])]
  outvecs = orthonormalise(invecs)
  call print_line('')
  call print_line('Invecs:')
  do i=1,size(invecs)
    call print_line(real(invecs(i)))
  enddo
  call print_line('Outvecs:')
  do i=1,size(outvecs)
    call print_line(real(outvecs(i)))
  enddo
  
  invecs = [vec([zero,one]), vec([one,zero]), vec([one,one])]
  outvecs = orthonormalise(invecs)
  call print_line('')
  call print_line('Invecs:')
  do i=1,size(invecs)
    call print_line(real(invecs(i)))
  enddo
  call print_line('Outvecs:')
  do i=1,size(outvecs)
    call print_line(real(outvecs(i)))
  enddo
  
  invecs = [vec([one,zero,zero]), vec([one,one,one]), vec([one,one,one])]
  outvecs = orthonormalise(invecs)
  call print_line('')
  call print_line('Invecs:')
  do i=1,size(invecs)
    call print_line(real(invecs(i)))
  enddo
  call print_line('Outvecs:')
  do i=1,size(outvecs)
    call print_line(real(outvecs(i)))
  enddo
  
  invecs = [vec([one,zero,zero]), vec([one,one,one]), vec([one,one,one])]
  outvecs = orthonormalise(invecs,0.1_dp)
  call print_line('')
  call print_line('Invecs:')
  do i=1,size(invecs)
    call print_line(real(invecs(i)))
  enddo
  call print_line('Outvecs:')
  do i=1,size(outvecs)
    call print_line(real(outvecs(i)))
  enddo
  
  invecs = [ vec([small,small,zero]),   &
           & vec([zero,small,zero]),  &
           & vec([small,zero,zero]),  &
           & vec([zero,zero,zero]), &
           & vec([one,one,zero]) ]
  outvecs = orthonormalise(invecs,1e-2_dp)
  call print_line('')
  call print_line('Invecs:')
  do i=1,size(invecs)
    call print_line(real(invecs(i)))
  enddo
  call print_line('Outvecs:')
  do i=1,size(outvecs)
    call print_line(real(outvecs(i)))
  enddo
  
end subroutine
end module
