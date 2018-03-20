! ======================================================================
! Sets up anharmonic calculations.
! ======================================================================
module setup_anharmonic_module
  use common_module
  
  use setup_harmonic_module
  
  use coupling_module
  use sampling_points_module
  use grid_types_module
  implicit none
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
function setup_anharmonic_keywords() result(keywords)
  implicit none
  
  type(KeywordData), allocatable :: keywords(:)
  
  keywords = [                                                                &
  & KeywordData( 'harmonic_path',                                            &
  &               'harmonic_path is the path to the directory where harmonic &
  &calculations were run.',                                                   &
  &               default_value='.',                                          &
  &               is_path=.true.),                                            &
  & KeywordData( 'temperature',                                              &
  &               'temperature is the temperature, in Kelvin, at which the &
  &simulation is run.',                                                       &
  &               default_value='0'),                                         &
  & KeywordData( 'grid_type',                                                &
  &               'grid_type specifies the sampling method. Options are &
  &"cubic", "octahedral", and "spherical".',                                  &
  &               default_value='cubic'),                                     &
  & KeywordData( 'max_energy',                                               &
  &               'max_energy is the maximum value of the potential up to &
  &which each normal mode will be evaluated.'),                               &
  & KeywordData( 'no_sampling_points',                                       &
  &               'no_sampling_points is the number of sampling points in &
  &each direction.'),                                                         &
  & KeywordData( 'coupling',                                                 &
  &               'coupling specifies the coupling between normal modes. &
  &Each set of coupled modes should be given as mode ids separated by spaces, &
  &and the sets should be separated by commas. Each q-point should be &
  &separated by semicolons. For example, if at q-point 1 modes 1, 2 and 3 are &
  &coupled and modes 4 and 5 are coupled, and at q-point 2 modes 1 and 6 are &
  &coupled, coupling="1 2 3, 4 5; 1 6". All couplings should be in ascending &
  &order, e.g. "1 2 3" rather than "1 3 2".',                                 &
  &               is_optional=.true.) ]
end function

function setup_anharmonic_mode() result(output)
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'setup_anharmonic'
  output%description = 'Sets up anharmonic calculations. Should be run after &
     &calculate_harmonic.'
  output%keywords = setup_anharmonic_keywords()
  output%main_subroutine => setup_anharmonic
end function

! ----------------------------------------------------------------------
! The main program.
! ----------------------------------------------------------------------
subroutine setup_anharmonic(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! User inputs.
  type(String)                    :: harmonic_path
  real(dp)                        :: temperature
  type(String)                    :: grid_type
  real(dp)                        :: max_energy
  integer                         :: no_sampling_points
  type(CoupledModes), allocatable :: coupling(:)
  
  ! Previous user inputs.
  type(Dictionary) :: setup_harmonic_arguments
  type(String)     :: seedname
  type(String)     :: file_type
  
  ! File data.
  type(String) :: wd
  type(String) :: input_filename
  type(IFile)  :: no_supercells_file
  
  ! Starting data.
  real(dp)                         :: thermal_energy
  type(StructureData)              :: structure
  integer                          :: no_supercells
  type(StructureData), allocatable :: supercells(:)
  type(QpointData),    allocatable :: qpoints(:)
  type(NormalMode),    allocatable :: modes(:)
  
  ! Main data structures.
  type(String),           allocatable :: all_coupling(:)
  type(CouplingSampling), allocatable :: setup_data(:)
  
  ! Normal mode data.
  real(dp), allocatable :: sample_spacing(:)
  
  ! Supercell with displaced atoms.
  type(StructureData)           :: supercell
  type(RealVector), allocatable :: displacement(:)
  
  ! Temporary variables.
  integer                   :: i,j,k,l,ialloc
  type(String), allocatable :: line(:)
  type(String)              :: qdir,cdir,sdir
  
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
  if (arguments%is_set('coupling')) then
    all_coupling = split(arguments%value('coupling'), ';')
  endif
  
  ! Read previous user inputs.
  setup_harmonic_arguments = setup_harmonic_keywords()
  call setup_harmonic_arguments%read_file( &
     & harmonic_path//'/setup_harmonic.used_settings')
  seedname = setup_harmonic_arguments%value('seedname')
  file_type = setup_harmonic_arguments%value('file_type')
  
  ! --------------------------------------------------
  ! Check inputs.
  ! --------------------------------------------------
  ! Check code is supported.
  if (file_type/='castep') then
    call print_line('Error: the code '//file_type//' is not yet supported.')
    stop
  endif
  
  ! Check temperature is valid.
  if (temperature < 0) then
    call print_line('Error: temperature must be positive.')
    stop
  endif
  
  ! Check grid type.
  if ( grid_type/='cubic'      .and. &
     & grid_type/='octahedral' .and. &
     & grid_type/='spherical') then
    call print_line('Error: the grid type '//grid_type//' is not yet &
       &supported.')
    stop
  endif
  
  ! Check dft input files exists.
  input_filename = make_input_filename(file_type,seedname)
  input_filename = wd//'/'//input_filename
  if (.not. file_exists(input_filename)) then
    call print_line('Error: The input file '//input_filename// &
       &' does not exist.')
    stop
  endif
  
  ! --------------------------------------------------
  ! Read and calculcate starting data.
  ! --------------------------------------------------
  ! Calculate thermal energy.
  thermal_energy = temperature*KB_IN_AU
  
  ! Read in crystal structure.
  structure = read_structure_file(harmonic_path//'/structure.dat')
  
  ! Read in supercell structures.
  no_supercells_file = harmonic_path//'/no_supercells.dat'
  no_supercells = int(no_supercells_file%line(1))
  allocate(supercells(no_supercells), stat=ialloc); call err(ialloc)
  do i=1,no_supercells
    supercells(i) = read_structure_file(                                &
       & harmonic_path//'/Supercell_'//left_pad(i,str(no_supercells))// &
       & '/structure.dat')
  enddo
  
  ! Read in q-points.
  qpoints = read_qpoints_file(harmonic_path//'/qpoints.dat')
  
  ! Check that the correct number of couplings have been specified.
  if (arguments%is_set('coupling')) then
    if (size(all_coupling)/=size(qpoints)) then
      call print_line('Error: the number of specified couplings does not &
         &match the number of q-points.')
      call err()
    endif
  else
    allocate(all_coupling(size(qpoints)), stat=ialloc); call err(ialloc)
    do i=1,size(all_coupling)
      all_coupling(i) = ''
    enddo
  endif
  
  ! --------------------------------------------------
  ! Loop across q-points, running calculations at each point.
  ! --------------------------------------------------
  do i=1,size(qpoints)
    ! Read in normal modes.
    allocate(modes(structure%no_modes), stat=ialloc); call err(ialloc)
    do j=1,structure%no_modes
      modes(j) = read_normal_mode_file( &
         & harmonic_path//'/qpoint_'//left_pad(i,str(size(modes)))// &
         & '/mode_'//left_pad(j,str(size(modes)))//'.dat')
    enddo
    
    ! Calculate the sample spacing along each mode.
    ! Assumes the mode is harmonic.
    allocate(sample_spacing(structure%no_modes), stat=ialloc); call err(ialloc)
    do j=1,structure%no_modes
      if (modes(j)%translational_mode) then
        sample_spacing(j) = 0
      else
        sample_spacing(j) = sqrt(2.0_dp*max_energy) &
                        & / (modes(j)%frequency*no_sampling_points)
      endif
    enddo
    
    ! Make q-point directories.
    qdir = wd//'/qpoint_'//left_pad(i,str(size(qpoints)))
    call mkdir(qdir)
    
    ! Parse coupling.
    line = split(all_coupling(i), ',')
    allocate(coupling(size(line)), stat=ialloc); call err(ialloc)
    do j=1,size(coupling)
      coupling(j)%modes = int(split(line(j)))
    enddo
    coupling = calculate_all_coupling(coupling, modes)
    
    ! Write out coupling.
    call write_coupling_file(coupling, qdir//'/coupling.dat')
    
    allocate(setup_data(size(coupling)), stat=ialloc); call err(ialloc)
    do j=1,size(coupling)
      setup_data(j)%coupling = coupling(j)
    enddo
    
    do j=1,size(setup_data)
      ! Make coupling directories.
      cdir = qdir//'/coupling_'//left_pad(j,str(size(setup_data)))
      call mkdir(cdir)
    
      ! Calculate indices of sampling points.
      setup_data(j)%sampling_points = generate_sampling_points(     &
                                    & grid_type,                    &
                                    & setup_data(j)%coupling%modes, &
                                    & structure%no_modes,           &
                                    & no_sampling_points,           &
                                    & sample_spacing)
      
      ! Write out sampling points.
      call write_sampling_points_file(setup_data(j)%sampling_points, &
         & cdir//'/sampling_points.dat')
      
      ! Make dft a working directory containing a DFT input file at each
      !    sampling point.
      do k=1,size(setup_data(j)%sampling_points)
        sdir = cdir//'/sampling_point_'// &
           & left_pad(k,size(setup_data(j)%sampling_points))
        call mkdir(sdir)
        
        ! Displace supercell atoms by a sum of normal modes.
        supercell = supercells(qpoints(i)%sc_id)
        displacement = normal_mode_to_cartesian(                       &
           & setup_data(j)%sampling_points(k)%displacement,            &
           & modes,                                                    &
           & qpoints(i),                                               &
           & supercell )
        do l=1,supercell%no_atoms
          call supercell%atoms(l)%set_cartesian_position( &
             & supercell%atoms(l)%cartesian_position() + displacement(l))
        enddo
        
        ! Write DFT input file.
        input_filename = make_input_filename(file_type,seedname)
        call StructureData_to_input_file( &
           & file_type,                   &
           & supercell,                   &
           & wd//'/'//input_filename,     &
           & sdir//'/'//input_filename)
      enddo
    enddo
    
    deallocate( modes,          &
              & sample_spacing, &
              & coupling,       &
              & setup_data,     &
              & stat=ialloc); call err(ialloc)
  enddo
end subroutine
end module
