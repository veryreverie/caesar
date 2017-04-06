module run_harmonic_module
  use constants_module, only : dp
  use string_module
  use io_module
contains

subroutine run_harmonic(wd,cwd)
  use constants_module, only : directions
  use utils_module,     only : format_path
  use unique_directions_module
  implicit none
  
  ! Working directories.
  type(String), intent(in) :: wd
  type(String), intent(in) :: cwd
  
  ! Terminal inputs.
  integer      :: first_sc
  integer      :: last_sc
  integer      :: no_cores
  type(String) :: run_script
  
  ! Previous user inputs.
  type(String), allocatable :: user_input_file(:)
  type(String)              :: dft_code
  type(String)              :: seedname
  
  ! Supercell data.
  type(String), allocatable :: no_sc_file(:)
  integer                   :: no_sc
  
  ! Atom and direction data.
  type(UniqueDirections) :: unique_directions
  integer                :: atom
  character(1)           :: direction
  logical                :: forces_calculated(3)
  
  integer                :: i,j,k
  type(String)           :: sdir
  type(String)           :: dir
  
  ! Read in supercell data.
  no_sc_file = read_lines(wd//'/no_sc.dat')
  no_sc = int(no_sc_file(1))
  
  ! Get inputs from user.
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
  
  ! Check user inputs.
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
  
  ! Read previous user inputs.
  user_input_file = read_lines(wd//'/user_input.txt')
  dft_code = user_input_file(1)
  seedname = user_input_file(2)
  
  ! Loop over supercells.
  do i=first_sc,last_sc
    sdir = wd//'/Supercell_'//i
    
    ! Read in unique directions.
    unique_directions = read_unique_directions_file( &
       & sdir//'/unique_directions.dat')
    
    do j=1,size(unique_directions)
      atom = unique_directions%unique_atoms(j)
      
      forces_calculated = (/ .true., .true., .true. /)
      if (unique_directions%xy_symmetry(j)/=0) then
        forces_calculated(2) = .false.
      endif
      
      if ( unique_directions%xz_symmetry(j)/=0 .or. &
         & unique_directions%yz_symmetry(j)/=0) then
        forces_calculated(3) = .false.
      endif
      
      do k=1,3
        if (.not. forces_calculated(k)) then
          cycle
        endif
        
        direction = directions(k)
        
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
  enddo
end subroutine
end module
