module calculate_gap_module
  use common_module
  implicit none
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
function calculate_gap_keywords() result(keywords)
  implicit none
  
  type(KeywordData), allocatable :: keywords(:)
  
  keywords = [KeywordData::]
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine calculate_gap(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  integer  :: i
  integer  :: no_points
  real(dp) :: temp
  real(dp) :: top
  real(dp) :: bottom
  
  type(String) :: top_filename
  type(String) :: bottom_filename
  
  type(IFile)               :: top_file
  type(IFile)               :: bottom_file
  type(String), allocatable :: line(:)
  
  type(OFile) :: gap_file
 
  call print_line('What is the file of the top band?')
  top_filename = read_line_from_user()
  call print_line('What is the file of the bottom band?')
  bottom_filename = read_line_from_user()
  call print_line('How many data points are there?')
  no_points = int(read_line_from_user())
  
  top_file = top_filename
  bottom_file = bottom_filename
  gap_file = 'gap.dat'
  
  do i=1,no_points
    line = split(top_file%line(i))
    temp = dble(line(1))
    top = dble(line(2))
    
    line = split(bottom_file%line(i))
    temp = dble(line(1))
    bottom = dble(line(2))
    
    call gap_file%print_line(temp//' '//top-bottom)
  enddo
end subroutine
end module
