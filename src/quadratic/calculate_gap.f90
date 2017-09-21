module calculate_gap_module
  use constants_module, only : dp
  use string_module
  use io_module
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
function calculate_gap_keywords() result(keywords)
  use keyword_module
  implicit none
  
  type(KeywordData), allocatable :: keywords(:)
  
  keywords = [KeywordData::]
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine calculate_gap(arguments)
  use dictionary_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  integer  :: i
  integer  :: no_points
  real(dp) :: temp
  real(dp) :: top
  real(dp) :: bottom
  
  type(String) :: top_filename
  type(String) :: bottom_filename
  
  type(String), allocatable :: top_file(:)
  type(String), allocatable :: bottom_file(:)
  type(String), allocatable :: line(:)
  
  integer      :: gap_file
 
  call print_line('What is the file of the top band?')
  top_filename = read_line_from_user()
  call print_line('What is the file of the bottom band?')
  bottom_filename = read_line_from_user()
  call print_line('How many data points are there?')
  no_points = int(read_line_from_user())
  
  top_file = read_lines(top_filename)
  bottom_file = read_lines(bottom_filename)
  gap_file = open_write_file('gap.dat')
  
  do i=1,no_points
    line = split(top_file(i))
    temp = dble(line(1))
    top = dble(line(2))
    
    line = split(bottom_file(i))
    temp = dble(line(1))
    bottom = dble(line(2))
    
    call print_line(gap_file,temp//' '//top-bottom)
  enddo
  close(gap_file)
end subroutine
end module
