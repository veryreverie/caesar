! ======================================================================
! Sets up anharmonic calculations.
! ======================================================================
module setup_anharmonic_module
  use constants_module, only : dp
  use string_module
  use io_module
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
function setup_anharmonic_keywords() result(keywords)
  use help_module
  implicit none
  
  type(KeywordData) :: keywords(6)
  
  keywords = [                                                                &
  & make_keyword( 'harmonic_path',                                            &
  &               'harmonic_path is the path to the directory where harmonic &
  &calculations were run.',                                                   &
  &               default_value='.',                                          &
  &               is_path=.true.),                                            &
  & make_keyword( 'temperature',                                              &
  &               'temperature is the temperature, in Kelvin, at which the &
  &simulation is run.',                                                       &
  &               default_value='0'),                                         &
  & make_keyword( 'grid_type',                                                &
  &               'grid_type specifies the sampling method. Options are &
  &"cubic", which uses a spline representation of the potential, and &
  &"spherical", which uses a polynomial representation of the potential.',    &
  &               default_value='cubic'),                                     &
  & make_keyword( 'max_energy',                                               &
  &               'max_energy is the maximum value of the potential up to &
  &which each normal mode will be evaluated.'),                               &
  & make_keyword( 'no_sampling_points',                                       &
  &               'no_sampling_points is the number of sampling points in &
  &each direction.'),                                                         &
  & make_keyword( 'coupling',                                                 &
  &               'coupling specifies the coupling between normal modes. &
  &Each set of coupled modes should be given as mode ids separated by spaces, &
  &and the sets should be separated by commas. Each q-point should be &
  &separated by semicolons. For example, if at q-point 1 modes 1,2 and 3 and &
  &modes 4 and 5 are coupled, and at q-point 2 modes 1 and 6 are coupled,     &
  &"coupling"="1 2 3, 4 5; 1 6".',                                            &
  &               default_value='') ]
end function

! ----------------------------------------------------------------------
! The main program.
! ----------------------------------------------------------------------
subroutine setup_anharmonic(arguments)
  use constants_module, only : kb_in_au
  use utils_module, only : mkdir
  use dictionary_module
  use dft_input_file_module
  use structure_module
  use qpoints_module
  use normal_mode_module
  use coupling_module
  use sampling_points_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! User inputs.
  type(String)                    :: harmonic_path
  real(dp)                        :: temperature
  type(String)                    :: grid_type
  real(dp)                        :: max_energy
  integer                         :: no_sampling_points
  type(String),       allocatable :: all_coupling(:)
  type(CoupledModes), allocatable :: coupling(:)
  
  ! Previous user inputs.
  type(Dictionary) :: setup_harmonic_arguments
  type(String)     :: seedname
  type(String)     :: dft_code
  
  ! File data.
  type(String) :: wd
  type(String) :: dft_input_filename
  
  type(String), allocatable :: no_supercells_file(:)
  
  ! Starting data.
  real(dp)                         :: thermal_energy
  type(StructureData)              :: structure
  integer                          :: no_supercells
  type(StructureData), allocatable :: supercells(:)
  type(QpointData),    allocatable :: qpoints(:)
  type(NormalMode),    allocatable :: modes(:)
  
  ! Normal mode data.
  real(dp), allocatable :: sample_spacing(:)
  
  ! Sampling points data.
  type(SamplingPoint), allocatable :: sampling_points(:)
  
  ! Supercell with displaced atoms.
  type(StructureData) :: supercell
  real(dp)            :: qr
  complex(dp)         :: exp_iqr
  type(ComplexVector) :: displacement
  
  ! Temporary variables.
  integer                   :: i,j,k,l,ialloc
  type(String), allocatable :: line(:)
  
  ! --------------------------------------------------
  ! Read inputs.
  ! --------------------------------------------------
  ! Read user inputs.
  wd = arguments%value('working_directory')
  harmonic_path = arguments%value('harmonic_path')
  temperature = dble(arguments%value('temperature'))
  grid_type = arguments%value('grid_type')
  max_energy = dble(arguments%value('max_energy'))
  no_sampling_points = int(arguments%value('no_sampling_points'))
  all_coupling = split(arguments%value('coupling'), ';')
  
  ! Read previous user inputs.
  call setup_harmonic_arguments%read_file( &
     & harmonic_path//'/setup_harmonic.used_settings')
  seedname = setup_harmonic_arguments%value('seedname')
  dft_code = setup_harmonic_arguments%value('dft_code')
  
  ! --------------------------------------------------
  ! Check inputs.
  ! --------------------------------------------------
  ! Check code is supported.
  if (dft_code/='castep') then
    call print_line('Error: the code '//dft_code//' is not yet supported.')
    stop
  endif
  
  ! Check input files exist.
  dft_input_filename = make_dft_input_filename(dft_code, seedname)
  dft_input_filename = wd//'/'//dft_input_filename
  if (.not. file_exists(dft_input_filename)) then
    call print_line('Error: the input file '//dft_input_filename// &
       & ' does not exist.')
    stop
  endif
  
  ! Check temperature is valid.
  if (temperature < 0) then
    call print_line('Error: temperature must be positive.')
    stop
  endif
  
  ! Check grid type.
  if (grid_type/='cubic') then
    call print_line('Error: the grid type '//grid_type//' is not yet &
       &supported.')
    stop
  endif
  
  ! --------------------------------------------------
  ! Read and calculcate starting data.
  ! --------------------------------------------------
  ! Calculate thermal energy.
  thermal_energy = temperature*kb_in_au
  
  ! Read in crystal structure.
  structure = read_structure_file(harmonic_path//'/structure.dat')
  
  ! Read in supercell structures.
  no_supercells_file = read_lines(harmonic_path//'/no_supercells.dat')
  no_supercells = int(no_supercells_file(1))
  allocate(supercells(no_supercells), stat=ialloc); call err(ialloc)
  do i=1,no_supercells
    supercells(i) = read_structure_file( &
       & harmonic_path//'/Supercell_'//i//'/structure.dat')
  enddo
  
  ! Read in q-points.
  qpoints = read_qpoints_file(harmonic_path//'/qpoints_ibz.dat')
  
  ! --------------------------------------------------
  ! Loop across q-points, running calculations at each point.
  ! --------------------------------------------------
  do i=1,size(qpoints)
    ! Make q-point directories.
    call mkdir(wd//'/qpoint_'//i)
    
    ! Read in normal modes.
    allocate(modes(structure%no_modes), stat=ialloc); call err(ialloc)
    do j=1,structure%no_modes
      modes(j) = read_normal_mode_file( &
         & harmonic_path//'/qpoint_'//i//'/mode_'//j//'.dat')
    enddo
    
    ! Parse coupling.
    line = split(all_coupling(i), ',')
    allocate(coupling(size(line)), stat=ialloc); call err(ialloc)
    do j=1,structure%no_modes
      coupling(j)%modes = int(split(line(j)))
    enddo
    coupling = calculate_all_coupling(coupling, structure%no_modes)
    
    ! Write out coupling.
    call write_coupling_file(coupling, &
       & wd//'/qpoint_'//i//'/coupling.dat')
    
    ! Calculate the sample spacing along each mode.
    ! Assumes the mode is harmonic.
    allocate(sample_spacing(structure%no_modes), stat=ialloc); call err(ialloc)
    do j=1,structure%no_modes
      sample_spacing(j) = sqrt(2.0_dp*max_energy) &
                      & / (modes(j)%frequency*no_sampling_points)
    enddo
    
    ! Generate sampling points.
    if (grid_type == 'cubic') then
      sampling_points = cubic_sampling_points( coupling,           &
                                             & structure%no_modes, &
                                             & no_sampling_points)
    endif
    
    ! Write out sampling points.
    call write_sampling_points_file(sampling_points, &
       & wd//'/qpoint_'//i//'/sampling_points.dat')
    
    ! Make dft a working directory containing a DFT input file at each
    !    sampling point.
    do j=1,size(sampling_points)
      call mkdir(wd//'/qpoint_'//i//'/sampling_point_'//j)
      
      supercell = supercells(qpoints(i)%sc_id)
      do k=1,supercell%no_atoms
        ! Calculate q.R and exp(i q.R)
        qr = qpoints(i)%qpoint * supercell%rvectors(supercell%atom_to_rvec(k))
        exp_iqr = cmplx(cos(qr), sin(qr), dp)
        
        ! Calculate displaced atomic co-ordinates.
        do l=1,structure%no_modes
          if (grid_type=='cubic') then
            displacement = modes(l)%displacements(supercell%atom_to_prim(k))
            supercell%atoms(k) = supercell%atoms(k)            &
                             & + sampling_points(j)%indices(l) &
                             & * real(displacement*exp_iqr)    &
                             & * sample_spacing(j)
          endif
        enddo
      enddo
      
      ! Write DFT input file.
      dft_input_filename = make_dft_input_filename(dft_code,seedname)
      call StructureData_to_dft_input_file( &
         & dft_code,                        &
         & supercell,                       &
         & wd//'/'//dft_input_filename,     &
         & wd//'/qpoint_'//i//'/sampling_pont_'//j//'/'//dft_input_filename)
    enddo
  enddo
end subroutine
end module
