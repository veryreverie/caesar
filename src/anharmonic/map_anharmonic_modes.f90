! ======================================================================
! Maps the anharmonic potential along each normal mode.
! ======================================================================
module map_anharmonic_modes_module
  use common_module
  
  use anharmonic_common_module
  use potentials_module

  use setup_harmonic_module
  
  use mode_map_module
  use setup_anharmonic_module
  implicit none
  
  private
  
  public :: map_anharmonic_modes
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
function map_anharmonic_modes() result(output)
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'map_anharmonic_modes'
  output%description = 'Maps the anharmonic potential along normal modes. &
     &Should be run after calculate_potential.'
  output%keywords = [                                                         &
     & KeywordData( 'no_single_mode_samples',                                 &
     &              'no_single_mode_samples is the number of points (either &
     &side of zero) along each mode at which the anharmonic potential will be &
     &sampled.',                                                              &
     &              default_value='100'),                                     &
     & KeywordData( 'validate_potential',                                     &
     &              'validate_potential specifies that the anharmonic &
     &potential should be verified against fresh electronic structure &
     &calculations when sampling. Depending on the electronic structure &
     &method, this is likely to be very computationally intensive.',          &
     &              default_value='false'),                                   &
     & KeywordData( 'run_script',                                             &
     &              'run_script is the path to the script for running DFT. An &
     &example run script can be found in doc/input_files.',                   &
     &               is_path=.true.),                                         &
     & KeywordData( 'no_cores',                                               &
     &              'no_cores is the number of cores on which the electronic &
     &structure calculation will be run. This is passed to the specified run &
     &script.',                                                               &
     &              default_value='1'),                                       &
     & KeywordData( 'no_nodes',                                               &
     &              'no_nodes is the number of nodes on which the electronic &
     &structure calculation will be run. This is passed to the specified run &
     &script.',                                                               &
     &              default_value='1'),                                       &
     & KeywordData( 'run_script_data',                                        &
     &              'run_script_data will be passed to the specified run &
     &script after all other arguments. This should be used to pass &
     &information not covered by the other arguments.',                       &
     &              default_value=''),                                        &
     & KeywordData( 'calculation_type',                                       &
     &              'calculation_type specifies whether any electronic &
     &structure calculations should be run in addition to the user-defined &
     &script. Settings are: "none" and "quip".',                              &
     &              default_value='none') ]
  output%main_subroutine => map_anharmonic_modes_subroutine
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine map_anharmonic_modes_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Input arguments.
  integer      :: no_single_mode_samples
  logical      :: validate_potential
  type(String) :: run_script
  integer      :: no_cores
  integer      :: no_nodes
  type(String) :: run_script_data
  type(String) :: calculation_type
  
  ! Arguments to setup_harmonic.
  type(Dictionary) :: setup_harmonic_arguments
  type(String)     :: seedname
  type(String)     :: file_type
  
  ! Anharmonic data.
  type(AnharmonicData)                  :: anharmonic_data
  type(StructureData)                   :: structure
  type(StructureData)                   :: anharmonic_supercell
  type(QpointData), allocatable         :: qpoints(:)
  type(RealMode),           allocatable :: real_modes(:)
  type(DegenerateSubspace), allocatable :: subspaces(:)
  real(dp)                              :: frequency_of_max_displacement
  
  ! Anharmonic potential.
  type(PotentialPointer) :: potential
  
  ! Electronic structure calculation handlers.
  type(CalculationWriter) :: calculation_writer
  type(CalculationRunner) :: calculation_runner
  type(CalculationReader) :: calculation_reader
  
  ! Variables for sampline potential.
  real(dp),      allocatable :: displacements(:)
  real(dp),      allocatable :: scaled_displacements(:)
  type(ModeMap), allocatable :: mode_maps(:)
  integer,       allocatable :: qpoint_modes(:)
  integer,       allocatable :: subspace_modes(:)
  type(ModeMap), allocatable :: subspace_mode_maps(:)
  
  ! Variables for validating potential.
  real(dp), allocatable       :: sampled_energies(:)
  type(RealModeForce)         :: sampled_force
  real(dp), allocatable       :: sampled_forces(:)
  type(IntMatrix)             :: supercell_matrix
  type(StructureData)         :: supercell
  type(RealMode)              :: mode
  type(RealModeDisplacement)  :: real_mode_displacement
  type(CartesianDisplacement) :: displacement
  type(StructureData)         :: displaced_structure
  type(ElectronicStructure)   :: electronic_structure
  
  ! Files and directories.
  type(IFile)  :: anharmonic_data_file
  type(IFile)  :: potential_file
  type(String) :: qpoint_dir
  type(String) :: subspace_dir
  type(OFile)  :: supercell_file
  type(OFile)  :: mode_maps_file
  type(String) :: mode_dir
  type(String) :: displacement_dir
  
  ! Temporary variables.
  integer :: i,j,k,ialloc
  
  ! --------------------------------------------------
  ! Read in inputs and previously calculated data.
  ! --------------------------------------------------
  no_single_mode_samples = int(arguments%value('no_single_mode_samples'))
  validate_potential = lgcl(arguments%value('validate_potential'))
  run_script = arguments%value('run_script')
  no_cores = int(arguments%value('no_cores'))
  no_nodes = int(arguments%value('no_nodes'))
  run_script_data = arguments%value('run_script_data')
  calculation_type = arguments%value('calculation_type')
  
  ! Read in setup_harmonic arguments.
  setup_harmonic_arguments = Dictionary(setup_harmonic())
  call setup_harmonic_arguments%read_file( &
          & 'setup_harmonic.used_settings' )
  seedname = setup_harmonic_arguments%value('seedname')
  file_type = setup_harmonic_arguments%value('file_type')
  
  ! Read in anharmonic data.
  anharmonic_data_file = IFile('anharmonic_data.dat')
  anharmonic_data = AnharmonicData(anharmonic_data_file%lines())
  
  structure = anharmonic_data%structure
  anharmonic_supercell = anharmonic_data%anharmonic_supercell
  qpoints = anharmonic_data%qpoints
  real_modes = anharmonic_data%real_modes
  subspaces = anharmonic_data%degenerate_subspaces
  frequency_of_max_displacement = anharmonic_data%frequency_of_max_displacement
  
  ! Read in anharmonic potential.
  potential_file = IFile('potential.dat')
  potential = PotentialPointer(potential_file%lines())
  
  ! --------------------------------------------------
  ! Initialise calculation handlers.
  ! --------------------------------------------------
  calculation_writer = CalculationWriter( file_type = file_type, &
                                        & seedname  = seedname   )
  
  calculation_runner = CalculationRunner(  &
     & file_type        = file_type,       &
     & seedname         = seedname,        &
     & run_script       = run_script,      &
     & no_cores         = no_cores,        &
     & no_nodes         = no_nodes,        &
     & run_script_data  = run_script_data, &
     & calculation_type = calculation_type )
  
  calculation_reader = CalculationReader()
  
  ! --------------------------------------------------
  ! Map the anharmonic potential.
  ! --------------------------------------------------
  ! Calculate displacements before scaling by 1/sqrt(frequency).
  displacements =                                               &
     &   [(i,i=-no_single_mode_samples,no_single_mode_samples)] &
     & * anharmonic_data%maximum_weighted_displacement          &
     & / no_single_mode_samples
  
  allocate(mode_maps(size(real_modes)), stat=ialloc); call err(ialloc)
  do i=1,size(real_modes)
    ! Scale displacement by 1/sqrt(frequency).
    scaled_displacements = displacements                              &
                       & * sqrt( frequency_of_max_displacement        &
                       &       / max( real_modes(i)%frequency,        &
                       &              frequency_of_max_displacement ) )
    
    ! Sample the model potential.
    mode_maps(i) = ModeMap( scaled_displacements, &
                          & real_modes(i),        &
                          & potential             )
  enddo
  
  ! --------------------------------------------------
  ! Validate the model potential by sampling the electronic structure at
  !    the same points as the model potential.
  ! --------------------------------------------------
  if (validate_potential) then
    do i=1,size(qpoints)
      if (qpoints(i)%paired_qpoint_id<qpoints(i)%id) then
        cycle
      endif
      qpoint_modes = filter(real_modes%qpoint_id_plus==qpoints(i)%id)
      if (qpoints(i)%paired_qpoint_id>qpoints(i)%id) then
        qpoint_modes = [ qpoint_modes,                                     &
                       & filter(real_modes%qpoint_id_minus==qpoints(i)%id) ]
      endif
      
      qpoint_dir = 'qpoint_'//left_pad(qpoints(i)%id, str(maxval(qpoints%id)))
      call mkdir(qpoint_dir)
      
      ! Construct and write out supercell.
      supercell_matrix = construct_supercell_matrix(qpoints(i), structure)
      supercell = construct_supercell( structure,       &
                                     & supercell_matrix )
      supercell_file = OFile(qpoint_dir//'/structure.dat')
      call supercell_file%print_lines(supercell)
      
      do j=1,size(qpoint_modes)
        mode = real_modes(qpoint_modes(j))
        mode_dir = qpoint_dir//'/real_mode_'// &
                 & left_pad(mode%id, str(maxval(real_modes%id)))
        call mkdir(mode_dir)
        allocate( sampled_energies(size(displacements)), &
                & sampled_forces(size(displacements)),   &
                & stat=ialloc); call err(ialloc)
        do k=1,size(sampled_energies)
          real_mode_displacement = RealModeDisplacement(     &
             & [mode],                                       &
             & [mode_maps(qpoint_modes(j))%displacements(k)] )
          displacement = CartesianDisplacement( real_mode_displacement, &
                                              & supercell,              &
                                              & real_modes,             &
                                              & qpoints                 )
          displaced_structure = displace_structure(supercell,displacement)
          
          displacement_dir = mode_dir//'/displacement_'// &
                           & left_pad(k,str(size(displacements)))
          call calculation_writer%write_calculation( &
                              & displaced_structure, &
                              & displacement_dir     )
          
          call calculation_runner%run_calculation(displacement_dir)
          
          electronic_structure = calculation_reader%read_calculation( &
                                                   & displacement_dir )
          
          sampled_energies(k) = electronic_structure%energy &
                            & / supercell%sc_size
          sampled_force = RealModeForce( electronic_structure%forces, &
                                       & supercell,                   &
                                       & real_modes,                  &
                                       & qpoints                      )
          sampled_forces(k) = sampled_force%force(mode)
        enddo
        sampled_energies = sampled_energies - potential%undisplaced_energy()
        mode_maps(qpoint_modes(j))%sampled_energies = sampled_energies
        mode_maps(qpoint_modes(j))%sampled_forces   = sampled_forces
        
        deallocate( sampled_energies, &
                  & sampled_forces,   &
                  & stat=ialloc); call err(ialloc)
      enddo
    enddo
  endif
  
  ! --------------------------------------------------
  ! Write out output.
  ! --------------------------------------------------
  do i=1,size(subspaces)
    subspace_modes = filter([(                         &
       & any(subspaces(i)%mode_ids==real_modes(j)%id), &
       & j=1,                                          &
       & size(real_modes)                              )])
    subspace_mode_maps = mode_maps(subspace_modes)
    subspace_dir = 'subspace_'//left_pad( subspaces(i)%id,          &
                                        & str(maxval(subspaces%id)) )
    call mkdir(subspace_dir)
    mode_maps_file = OFile(subspace_dir//'/anharmonic_mode_maps.dat')
    call mode_maps_file%print_line(                    &
       & 'Harmonic frequencies: '//subspaces%frequency )
    call mode_maps_file%print_line('')
    call mode_maps_file%print_lines(subspace_mode_maps, separating_line='')
  enddo
end subroutine
end module