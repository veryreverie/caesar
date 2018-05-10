! ======================================================================
! The second stage of Caesar's anharmonic process.
! Runs DFT for anharmonic calculations.
! ======================================================================
module run_quadratic_module
  use common_module
  implicit none
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
function run_quadratic_keywords() result(keywords)
  implicit none
  
  type(KeywordData), allocatable :: keywords(:)
  
  keywords = [                                                                &
  & KeywordData( 'supercells_to_run',                                        &
  &               'supercells_to_run is the first and last supercell to run. &
  &These should be specified as two integers separated by spaces.'),          &
  & KeywordData( 'no_cores',                                                 &
  &               'no_cores is the number of cores on which DFT will be run. &
  &This is passed to the specified run script.',                              &
  &               default_value='1'),                                         &
  & KeywordData('run_script',                                                &
  &              'run_script is the path to the script for running DFT. An &
  &example run script can be found in doc/input_files.',                      &
  &              is_path=.true.) ]
end function

! ----------------------------------------------------------------------
! The main program.
! ----------------------------------------------------------------------
subroutine run_quadratic(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Working directory.
  type(String) :: wd
  
  ! Terminal inputs.
  integer      :: supercells_to_run(2)
  integer      :: no_cores
  type(String) :: run_script
  
  ! Previous user inputs.
  type(Dictionary)  :: setup_quadratic_arguments
  type(String)      :: dft_code
  type(String)      :: seedname
  type(String)      :: harmonic_path
  integer           :: no_samples
  
  ! Previously calculated data.
  type(StructureData)           :: structure
  type(IFile)                   :: no_supercells_file
  integer                       :: no_supercells
  type(QpointData), allocatable :: qpoints(:)
  
  ! Temporary variables.
  integer      :: i,j,k
  type(String) :: dir
  type(String) :: dft_input_filename
  integer      :: result_code
  
  ! --------------------------------------------------
  ! Get inputs from user.
  ! --------------------------------------------------
  wd = arguments%value('working_directory')
  supercells_to_run = int(split(arguments%value('supercells_to_run')))
  no_cores = int(arguments%value('no_cores'))
  run_script = arguments%value('run_script')
  
  ! Read in number of supercells.
  no_supercells_file = harmonic_path//'/no_supercells.dat'
  no_supercells = int(no_supercells_file%line(1))
  
  ! --------------------------------------------------
  ! Check user inputs.
  ! --------------------------------------------------
  if (supercells_to_run(1)<=0) then
    call print_line('')
    call print_line('Error: first supercell must be > 0')
    call err()
  elseif (supercells_to_run(1)>no_supercells) then
    call print_line('')
    call print_line('Error: first supercell must be <= '//no_supercells)
    call err()
  elseif (supercells_to_run(2)<supercells_to_run(1)) then
    call print_line('')
    call print_line('Error: first supercell must be <= last supercell.')
    call err()
  elseif (supercells_to_run(2)>no_supercells) then
    call print_line('')
    call print_line('Error: last supercell must be <= '//no_supercells)
    call err()
  elseif (no_cores<=0) then
    call print_line('')
    call print_line('Error: no. cores must be >= 1')
    call err()
  elseif (.not. file_exists(run_script)) then
    call print_line('')
    call print_line('Error: '//run_script//' does not exist.')
    call err()
  endif
  
  ! --------------------------------------------------
  ! Read in previous user inputs.
  ! --------------------------------------------------
  call setup_quadratic_arguments%read_file( &
     & wd//'/setup_quadratic.used_settings')
  dft_code = setup_quadratic_arguments%value('dft_code')
  seedname = setup_quadratic_arguments%value('seedname')
  harmonic_path = setup_quadratic_arguments%value('harmonic_path')
  no_samples = int(setup_quadratic_arguments%value('no_samples'))
  
  ! Read in structure.
  structure = read_structure_file(harmonic_path//'/structure.dat')
  
  ! Run static calculations.
  do i=supercells_to_run(1),supercells_to_run(2)
    dir=wd//'/Supercell_'//i
    call print_line('')
    call print_line('Running static calculation in directory '//dir)
    result_code = system_call( 'cd '//wd//'; '//run_script//' '// &
       & dft_code//' '//dir//' '//no_cores//' '//seedname)
    call print_line('Result code: '//result_code)
  enddo
  
  ! Read in IBZ q-points.
  qpoints = read_qpoints_file(harmonic_path//'/qpoints_ibz.dat')
  
  ! Loop over q-points.
  do i=1,size(qpoints)
    ! Ignore q-points outside of chosen supercells.
    if ( qpoints(i)%sc_id<supercells_to_run(1) .or. &
       & qpoints(i)%sc_id>supercells_to_run(2)) then
      cycle
    endif
    
    do j=1,structure%no_modes
      do k=-no_samples,no_samples
        dir = wd//'/qpoint_'//i//'/mode_'//j//'/amplitude_'//k
        dft_input_filename = make_dft_input_filename(dft_code,seedname)
        if (file_exists(dir//'/'//dft_input_filename)) then
          call print_line('')
          call print_line('Running calculation in directory '//dir)
          result_code = system_call( 'cd '//wd//'; '//run_script//' '// &
             & dft_code//' '//dir//' '//no_cores//' '//seedname)
          call print_line('Result code: '//result_code)
        endif
      enddo
    enddo
  enddo
end subroutine
end module