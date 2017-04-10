module run_harmonic_module
  use constants_module, only : dp
  use string_module
  use io_module
contains

! ----------------------------------------------------------------------
! Runs DFT for harmonic calculations.
! ----------------------------------------------------------------------
subroutine run_harmonic(wd,cwd)
  use utils_module,     only : format_path
  use unique_directions_module
  implicit none
  
  ! Working directories.
  type(String), intent(in) :: wd
  type(String), intent(in) :: cwd
  
  ! Previous user inputs.
  type(String), allocatable :: no_sc_file(:)
  integer                   :: no_sc
  type(String), allocatable :: user_input_file(:)
  type(String)              :: dft_code
  type(String)              :: seedname
  
  ! Terminal inputs.
  integer      :: first_sc
  integer      :: last_sc
  integer      :: no_cores
  type(String) :: run_script
  
  ! Atom and direction data.
  type(UniqueDirections) :: unique_directions
  integer                :: atom
  character(1)           :: direction
  
  ! Temporary variables.
  integer                :: i,j
  type(String)           :: sdir
  type(String)           :: dir
  
  ! --------------------------------------------------
  ! Read in previous user inputs.
  ! --------------------------------------------------
  no_sc_file = read_lines(wd//'/no_sc.dat')
  no_sc = int(no_sc_file(1))
  
  user_input_file = read_lines(wd//'/user_input.txt')
  dft_code = user_input_file(1)
  seedname = user_input_file(2)
  
  ! --------------------------------------------------
  ! Get inputs from user.
  ! --------------------------------------------------
  call print_line('')
  call print_line('There are '//no_sc//' supercells in total.')
  call print_line('What is the first supercell to run?')
  first_sc = int(read_line_from_user())
  
  call print_line('')
  call print_line('What is the last supercell to run?')
  last_sc = int(read_line_from_user())
  
  call print_line('')
  call print_line('How many cores can be used?')
  no_cores = int(read_line_from_user())
  
  call print_line('')
  call print_line('What is the path to the script for running DFT?')
  call print_line('(An example script can be found in doc/input_files)')
  run_script = format_path(read_line_from_user(), cwd)
  
  ! --------------------------------------------------
  ! Check user inputs.
  ! --------------------------------------------------
  if (first_sc<=0) then
    call print_line('')
    call print_line('Error: first supercell must be > 0')
    call err()
  elseif (first_sc>no_sc) then
    call print_line('')
    call print_line('Error: first supercell must be <= '//no_sc)
    call err()
  elseif (last_sc<first_sc) then
    call print_line('')
    call print_line('Error: first supercell must be <= last supercell.')
    call err()
  elseif (last_sc>no_sc) then
    call print_line('')
    call print_line('Error: last supercell must be <= '//no_sc)
    call err()
  elseif (no_cores<=0) then
    call print_line('')
    call print_line('Error: no. cores must be >= 1')
    call err()
  elseif (.not. file_exists(run_script)) then
    call print_line('')
    call print_line('Error: '//run_script//' does not exist.')
  endif
  
  ! --------------------------------------------------
  ! Run calculations
  ! --------------------------------------------------
  do i=first_sc,last_sc
    sdir = wd//'/Supercell_'//i
    
    unique_directions = read_unique_directions_file( &
       & sdir//'/unique_directions.dat')
    
    do j=1,size(unique_directions)
      atom = unique_directions%atoms(j)
      direction = unique_directions%directions_char(j)
      
      dir = sdir//'/atom.'//atom//'.+d'//direction
      call print_line('Running calculation in directory '//dir)
      call system_call('cd '//wd//'; '//run_script//' '//dft_code //' '// &
                                                       & dir      //' '// &
                                                       & no_cores //' '// &
                                                       & seedname)
      
      dir = sdir//'/atom.'//atom//'.-d'//direction
      call print_line('Working in directory '//dir)
      call system_call('cd '//wd//'; '//run_script//' '//dft_code //' '// &
                                                       & dir      //' '// &
                                                       & no_cores //' '// &
                                                       & seedname)
    enddo
  enddo
end subroutine
end module
