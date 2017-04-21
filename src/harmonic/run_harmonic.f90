! ======================================================================
! The second stage of Caesar.
! Runs DFT for harmonic calculations.
! ======================================================================
module run_harmonic_module
  use constants_module, only : dp
  use string_module
  use io_module
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
function run_harmonic_keywords() result(keywords)
  use help_module
  implicit none
  
  type(KeywordData) :: keywords(3)
  
  keywords = [ &
  & make_keyword('supercells_to_run', no_argument, 'supercells_to_run &
     &is first and last supercell to run. These should be specified as two &
     &integers separated by spaces.'),                                        &
  & make_keyword('no_cores', '1', 'no_cores is the number of cores on which &
     &DFT will be run. This is passed to the specified run script.'),         &
  & make_keyword('run_script', no_argument, 'run_script is the path to the &
     &script for running DFT. An example run script can be found in &
     &doc/input_files.')                                                      ]
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine run_harmonic(arguments,cwd)
  use utils_module,     only : format_path
  use unique_directions_module
  use dictionary_module
  implicit none
  
  ! Working directories.
  type(Dictionary), intent(in) :: arguments
  type(String),     intent(in) :: cwd
  
  ! Working directory
  type(String) :: wd
  
  ! Previous user inputs.
  type(Dictionary)          :: setup_harmonic_arguments
  type(String), allocatable :: no_sc_file(:)
  integer                   :: no_sc
  type(String)              :: dft_code
  type(String)              :: seedname
  
  ! Terminal inputs.
  integer      :: supercells_to_run(2)
  integer      :: no_cores
  type(String) :: run_script
  
  ! Atom and direction data.
  type(UniqueDirections) :: unique_directions
  integer                :: atom
  character(1)           :: direction
  
  ! Temporary variables.
  integer                   :: i,j
  type(String)              :: sdir
  type(String)              :: dir
  
  ! --------------------------------------------------
  ! Read in arguments to previous calculations.
  ! --------------------------------------------------
  no_sc_file = read_lines(wd//'/no_sc.dat')
  no_sc = int(no_sc_file(1))
  
  setup_harmonic_arguments = read_dictionary_file( &
     & wd//'/setup_harmonic.used_settings')
  dft_code = item(setup_harmonic_arguments, 'dft_code')
  seedname = item(setup_harmonic_arguments, 'seedname')
  
  ! --------------------------------------------------
  ! Get inputs from user.
  ! --------------------------------------------------
  wd = item(arguments, 'working_directory')
  supercells_to_run = int(split(item(arguments, 'supercells_to_run')))
  no_cores = int(item(arguments, 'no_cores'))
  run_script = format_path(item(arguments, 'run_script'), cwd)
  
  ! --------------------------------------------------
  ! Check user inputs.
  ! --------------------------------------------------
  if (supercells_to_run(1)<=0) then
    call print_line('')
    call print_line('Error: first supercell must be > 0')
    call err()
  elseif (supercells_to_run(1)>no_sc) then
    call print_line('')
    call print_line('Error: first supercell must be <= '//no_sc)
    call err()
  elseif (supercells_to_run(2)<supercells_to_run(1)) then
    call print_line('')
    call print_line('Error: first supercell must be <= last supercell.')
    call err()
  elseif (supercells_to_run(2)>no_sc) then
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
  do i=supercells_to_run(1),supercells_to_run(2)
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
