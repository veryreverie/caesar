module test_copy_quadratic_module
  use constants_module, only : dp
  use string_module
  use io_module
contains

subroutine test_copy_quadratic(wd)
  implicit none
  
  ! Working directories.
  type(String), intent(in) :: wd
  
  ! Terminal inputs.
  type(String) :: copy_dir
  
  ! Setup data
  type(String), allocatable :: user_input_file(:)
  type(String)              :: dft_code
  type(String)              :: seedname
  type(String)              :: harmonic_path
  type(String), allocatable :: no_sc_file(:)
  integer                   :: no_sc
  
  ! Directory and file names.
  type(String) :: sdir
  
  ! Temporary variables.
  integer :: i
  integer :: result_code
  
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
  
  no_sc_file = read_lines(harmonic_path//'/no_sc.dat')
  no_sc = int(no_sc_file(1))
  
  do i=1,no_sc
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
