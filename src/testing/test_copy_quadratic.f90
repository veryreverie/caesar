module test_copy_quadratic_module
  use constants_module, only : dp
  use string_module
  use io_module
contains

! ----------------------------------------------------------------------
! Generates keywords and helptext.
! ----------------------------------------------------------------------
function test_copy_quadratic_keywords() result(keywords)
  use keyword_module
  use setup_quadratic_module
  implicit none
  
  type(KeywordData), allocatable :: keywords(:)
  
  type(KeywordData), allocatable :: keywords_setup_quadratic(:)
  type(KeywordData)              :: keywords_test(0)
  
  integer :: ialloc
  
  keywords_setup_quadratic = setup_quadratic_keywords()
  
  allocate( keywords(size(keywords_setup_quadratic)+size(keywords_test)), &
          & stat=ialloc); call err(ialloc)
  
  keywords(:size(keywords_setup_quadratic)) = keywords_setup_quadratic
  keywords(size(keywords_setup_quadratic)+1:) = keywords_test
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine test_copy_quadratic(arguments)
  use dictionary_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Working directories.
  type(String) :: wd
  
  ! Terminal inputs.
  type(String) :: copy_dir
  
  ! Setup data
  type(String), allocatable :: user_input_file(:)
  type(String)              :: dft_code
  type(String)              :: seedname
  type(String)              :: harmonic_path
  type(String), allocatable :: no_supercells_file(:)
  integer                   :: no_supercells
  
  ! Directory and file names.
  type(String) :: sdir
  
  ! Temporary variables.
  integer :: i
  integer :: result_code
  
  ! --------------------------------------------------
  ! Read in arguments from user.
  ! --------------------------------------------------
  wd = arguments%value('working_directory')
  
  ! ----------------------------------------------------------------------
  ! Read in directory to compare against and copy from.
  ! ----------------------------------------------------------------------
  call print_line('')
  call print_line('This test will check dft input files against a previous &
     &calculations,')
  call print_line('   and copy over dft output files.')
  call print_line('')
  call print_line('Where is the quadratic directory for comparison?')
  copy_dir = format_path(read_line_from_user())
  
  ! ----------------------------------------------------------------------
  ! Read in previous settings.
  ! ----------------------------------------------------------------------
  user_input_file = read_lines(wd//'/user_input.txt')
  dft_code = user_input_file(1)
  seedname = user_input_file(2)
  harmonic_path = user_input_file(3)
  
  no_supercells_file = read_lines(harmonic_path//'/no_supercells.dat')
  no_supercells = int(no_supercells_file(1))
  
  do i=1,no_supercells
    sdir = wd//'/Supercell_'//i
    if (dft_code == 'castep') then
      result_code =  system_call( 'cp '// &
         & copy_dir//'/'//sdir//'/static/'//seedname//'.castep '// &
         & sdir//'/static/'//seedname//'.castep')
      if (result_code/=0) then
        call err()
      endif
    endif
  enddo
end subroutine
end module
