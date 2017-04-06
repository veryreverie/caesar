module run_quadratic_module
  use constants_module, only : dp
  use string_module
  use io_module
contains

subroutine run_quadratic(wd,cwd)
  use utils_module, only : format_path, make_dft_input_filename
  use structure_module
  use mapping_module
  use kpoints_module
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
  type(String)              :: harmonic_path
  type(MappingData)         :: mapping
  
  ! Previously calculated data.
  type(StructureData)       :: structure
  type(String), allocatable :: no_sc_file(:)
  integer                   :: no_sc
  type(KpointsIbz)          :: kpoints
  
  ! Temporary variables.
  integer      :: i,j,k
  type(String) :: dir
  type(String) :: dft_input_filename
  
  ! Read in number of supercells.
  no_sc_file = read_lines(harmonic_path//'/no_sc.dat')
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
  run_script = format_path(read_line_from_user(),cwd)
  
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
  
  ! Read in previous user inputs.
  user_input_file = read_lines(wd//'/user_input.txt')
  dft_code = user_input_file(1)
  seedname = user_input_file(2)
  harmonic_path = user_input_file(3)
  
  ! Read in maping data.
  mapping = read_mapping_file(wd//'/mapping.dat')
  
  ! Read in structure.
  structure = read_structure_file(harmonic_path//'/structure.dat')
  
  ! Run static calculations.
  do i=first_sc,last_sc
    dir=wd//'/Supercell_'//i
    call print_line('Running static calculation in directory '//dir)
    call system_call('cd '//wd//'; '//run_script//' '//dft_code //' '// &
                                                     & dir      //' '// &
                                                     & no_cores //' '// &
                                                     & seedname)
  enddo
  
  ! Read in IBZ K-points.
  kpoints = read_kpoints_ibz_file(harmonic_path//'/kpoints_ibz.dat')
  
  ! Loop over K-points.
  do i=1,size(kpoints)
    ! Ignore K-points outside of chosen supercells.
    if (kpoints%sc_ids(i)<first_sc .or. kpoints%sc_ids(i)>last_sc) then
      cycle
    endif
    
    do j=1,structure%no_modes
      do k=mapping%first,mapping%last
        dir = wd//'/kpoint_'//i//'/mode_'//j//'/amplitude_'//k
        dft_input_filename = make_dft_input_filename(dft_code,seedname)
        if (file_exists(dir//'/'//dft_input_filename)) then
          call print_line('Running calculation in directory '//dir)
          call system_call('cd '//wd//'; '//run_script//' '//dft_code //' '// &
                                                           & dir      //' '// &
                                                           & no_cores //' '// &
                                                           & seedname)
        endif
      enddo
    enddo
  enddo
end subroutine
end module
