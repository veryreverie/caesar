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
  use utils_module, only : gcd
  use dictionary_module
  
  use fraction_module
  use fraction_algebra_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  type(String) :: wd
  
  integer :: i1
  integer :: i2
  
  type(IntFraction) :: f1
  type(IntFraction) :: f2
  type(IntFraction) :: f3
  type(IntFraction) :: f4
  
  type(FractionVector) :: v1
  type(FractionVector) :: v2
  
  type(FractionMatrix) :: m1
  type(FractionMatrix) :: m2
  
  type(String), allocatable :: line(:)
  real(dp),     allocatable :: dbles(:)
  integer                   :: i
  
  complex(dp) :: thing
  
  wd = arguments%value('working_directory')
  
  thing = cmplx(1.2_dp*sqrt(2.0_dp), 3.4_dp*sqrt(2.0_dp), dp)
  
  call print_line(thing)
  write(*,'(es25.17,sp,es24.17,"i")') thing
  
  stop
  
  call print_line(dble('0.0'))
  call print_line(dble('0'))
  call print_line(dble('0/1'))
  call print_line(vec(dble(split('1/2 0 0'))))
  
  line = split('1/2 0 0')
  do i=1,size(line)
    call print_line('line('//i//'): "'//line(i)//'"')
  enddo
  
  dbles = dble(line)
  do i=1,size(dbles)
    call print_line('dbles('//i//'): "'//dbles(i)//'"')
  enddo
  
  call print_line(dble('1/2'))
  
  i1 = 6
  i2 = 3
  call print_line('6 3 3: '//i1//' '//i2//' '//gcd(i1,i2))
  
  i1 = 2
  i2 = 0
  call print_line('2 0 2: '//i1//' '//i2//' '//gcd(i1,i2))
  
  i1 = 5
  i2 = 7
  call print_line('5 7 1: '//i1//' '//i2//' '//gcd(i1,i2))
  
  f1 = 3
  f1 = f1 / 6
  call print_line('1/2  :'//f1)
  
  f2 = IntFraction(7,3)
  call print_line('7/3  :'//f2)
  
  f3 = frac(4)
  f3 = f3/8
  call print_line('1/2  :'//f3)
  
  f4 = IntFraction(1,-4)
  f4 = f4*2
  call print_line('-1/2 :'//f4)
  
  f4 = frac(str(f4))
  call print_line('-1/2 :'//f4)
  
  call print_line('-1/1 :'//f4-f3)
  call print_line('0/1  :'//f4+f3)
  
  call print_line('F    :'//is_int(f4))
  call print_line('T    :'//is_int(f4*2))
  call print_line('T    :'//is_int(f4+f3))
  call print_line('7/3  :'//dble(f2))
  call print_line('2    :'//int(f2))
  
  call print_line('')
  v1 = frac(vec([1,2,3]))
  call print_line('[1/1,2/1,3/1]: '//v1)
  
  v1 = v1/2
  call print_line('^/2          : '//v1)
  
  v2 = frac(vec([2,3,4]))
  call print_line('[2/1,3/1,4/1]: '//v2)
  
  v2 = v2/3
  call print_line('^/3          : '//v2)
  call print_line('10/3         : '//v1*v2)
  
  call print_line('')
  m1 = frac(mat( [ 1,2,3,4,5, &
               &   6,7,8,9,10 ], 5, 2))
  call print_line('m1 :')
  call print_line(m1)
  
  m1 = m1/IntFraction(4,3)
  call print_line('m1 / 4/3:')
  call print_line(m1)
  
  m2 = frac(mat( [ 3,5, &
               &   7,9  ], 2,2))
  call print_line('m2:')
  call print_line(m2)
  call print_line('m2*m1:')
  call print_line(m2*m1)
end subroutine
end module
