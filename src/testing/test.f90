! ======================================================================
! A test space, for temporary checking of misc. parts of the code.
! ======================================================================
module test_module
  use common_module
  implicit none
  
  private
  
  public :: test_keywords
  public :: test_mode
  public :: test
contains

! ----------------------------------------------------------------------
! Generates keywords and helptext.
! ----------------------------------------------------------------------
function test_keywords() result(keywords)
  implicit none
  
  type(KeywordData), allocatable :: keywords(:)
  
  keywords = [                                                    &
     & KeywordData('ifile','ifile is a filename',is_path=.true.), &
     & KeywordData('ofile','ofile is a filename',is_path=.true.)]
end function

function test_mode() result(output)
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'test'
  output%description = 'Runs temporary code for testing purposes.'
  output%keywords = test_keywords()
  output%main_subroutine => test
  output%suppress_from_helptext = .true.
  output%suppress_settings_file = .true.
end function

subroutine temp(input)
  implicit none
  
  type(OFile), intent(in) :: input
  
  type(OFile) :: thing
  
  call print_line('=====')
  thing = input
  call print_line('-----')
end subroutine

! ----------------------------------------------------------------------
! Main function.
! ----------------------------------------------------------------------
subroutine test(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  type(String) :: wd
  
  type(String)                   :: in_filename
  type(IFile)                    :: in_file
  type(String)                   :: out_filename
  type(OFile)                    :: out_file
  type(OFile)                    :: out_file_2
  type(StringArray), allocatable :: strings(:)
  
  integer :: i
  
  wd = arguments%value('working_directory')
  
  in_filename = arguments%value('ifile')
  in_file = IFile(in_filename)
  strings = split_into_sections(in_file%lines())
  
  out_filename = arguments%value('ofile')
  call print_line('==================================================')
  !call print_line('========================================')
  out_file = OFile(out_filename)
  !call print_line('========================================')
  call print_line(out_file%counter%is_only_pointer())
  call print_line('==================================================')
  
  call temp(out_file)
  
  call print_line(out_file%counter%is_only_pointer())
  
  out_file_2 = out_file
  
  call print_line(out_file%counter%is_only_pointer())
  
  out_file_2 = OFile('paraqueet dispenser')
  
  call print_line(out_file%counter%is_only_pointer())
  
  stop
  
  
  call temp(out_file_2)
  
  do i=1,size(strings)
    call print_line('Section '//i//':')
    call out_file%print_line('Section '//i//':')
    call print_lines(strings(i)%strings)
    call out_file%print_lines(strings(i)%strings)
  enddo
end subroutine
end module
