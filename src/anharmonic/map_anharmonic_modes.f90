! ======================================================================
! Calculates anharmonic states using the anharmonic potential and VSCF.
! ======================================================================
module map_anharmonic_modes_module
  use common_module
  
  use anharmonic_common_module
  use polynomial_module
  
  use setup_harmonic_module
  use calculate_normal_modes_module
  
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
     &sampled when determining the effective frequency with which the &
     &harmonic basis along that mode will be constructed.',                   &
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
     &              'no_cores is the number of cores on which DFT will be &
     &run. This is passed to the specified run script.',                      &
     &               default_value='1'),                                      &
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
  
  ! Working directory.
  type(String) :: wd
  
  ! Input arguments.
  integer      :: no_single_mode_samples
  logical      :: validate_potential
  type(String) :: run_script
  integer      :: no_cores
  type(String) :: calculation_type
  
  ! Arguments to setup_anharmonic.
  type(Dictionary) :: setup_anharmonic_arguments
  type(String)     :: harmonic_path
  type(String)     :: potential_representation
  real(dp)         :: maximum_displacement
  real(dp)         :: frequency_of_max_displacement
  
  ! Arguments to calculate_normal_modes.
  type(Dictionary) :: calculate_normal_modes_arguments
  
  ! Arguments to setup_harmonic.
  type(Dictionary) :: setup_harmonic_arguments
  type(String)     :: seedname
  type(String)     :: file_type
  
  ! Maximum displacement in mass-weighted co-ordinates.
  real(dp) :: maximum_weighted_displacement
  
  ! Primitive structure.
  type(StructureData) :: structure
  
  ! Large anharmonic supercell and its q-points.
  type(StructureData)           :: anharmonic_supercell
  type(QpointData), allocatable :: qpoints(:)
  
  ! Normal modes.
  type(ComplexMode), allocatable :: complex_modes(:)
  type(RealMode),    allocatable :: real_modes(:)
  
  ! Anharmonic potential.
  type(PotentialPointer) :: potential
  
  ! Electronic structure calculation handlers.
  type(CalculationWriter) :: calculation_writer
  type(CalculationRunner) :: calculation_runner
  type(CalculationReader) :: calculation_reader
  
  ! Variables for generating effective mode frequencies.
  real(dp),                 allocatable :: displacements(:)
  real(dp),                 allocatable :: scaled_displacements(:)
  type(EffectiveFrequency), allocatable :: effective_frequencies(:)
  integer,                  allocatable :: qpoint_modes(:)
  type(EffectiveFrequency), allocatable :: qpoint_frequencies(:)
  
  ! Variables for validating potential.
  real(dp), allocatable       :: sampled_energies(:)
  type(RealModeForce)         :: sampled_force
  real(dp), allocatable       :: sampled_forces(:)
  type(IntMatrix)             :: supercell_matrix
  type(StructureData)         :: supercell
  type(ComplexMode)           :: complex_mode
  type(RealMode)              :: real_mode
  type(RealModeDisplacement)  :: real_mode_displacement
  type(CartesianDisplacement) :: displacement
  type(StructureData)         :: displaced_structure
  type(ElectronicStructure)   :: electronic_structure
  
  ! Files and directories.
  type(IFile)  :: structure_file
  type(IFile)  :: anharmonic_supercell_file
  type(IFile)  :: qpoints_file
  type(IFile)  :: complex_modes_file
  type(IFile)  :: real_modes_file
  type(IFile)  :: potential_file
  type(String) :: qpoint_dir
  type(OFile)  :: supercell_file
  type(OFile)  :: effective_frequencies_file
  type(String) :: mode_dir
  type(String) :: displacement_dir
  
  ! Temporary variables.
  integer :: i,j,k,ialloc
  
  ! --------------------------------------------------
  ! Read in inputs and previously calculated data.
  ! --------------------------------------------------
  
  wd = arguments%value('working_directory')
  
  no_single_mode_samples = int(arguments%value('no_single_mode_samples'))
  validate_potential = lgcl(arguments%value('validate_potential'))
  run_script = arguments%value('run_script')
  no_cores = int(arguments%value('no_cores'))
  calculation_type = arguments%value('calculation_type')
  
  ! Read in setup_anharmonic arguments.
  setup_anharmonic_arguments = Dictionary(setup_anharmonic())
  call setup_anharmonic_arguments%read_file( &
     & wd//'/setup_anharmonic.used_settings')
  harmonic_path = setup_anharmonic_arguments%value('harmonic_path')
  potential_representation = &
     & setup_anharmonic_arguments%value('potential_representation')
  maximum_displacement = &
     & dble(setup_anharmonic_arguments%value('maximum_displacement'))
  frequency_of_max_displacement = &
     & dble(setup_anharmonic_arguments%value('frequency_of_max_displacement'))
  
  ! Read in calculate_normal_modes arguments.
  calculate_normal_modes_arguments = Dictionary(calculate_normal_modes())
  call calculate_normal_modes_arguments%read_file( &
     & harmonic_path//'/calculate_normal_modes.used_settings')
  
  ! Read in setup_harmonic arguments.
  setup_harmonic_arguments = Dictionary(setup_harmonic())
  call setup_harmonic_arguments%read_file( &
     & harmonic_path//'/setup_harmonic.used_settings')
  seedname = setup_harmonic_arguments%value('seedname')
  file_type = setup_harmonic_arguments%value('file_type')
  
  ! Read in structure.
  structure_file = IFile(harmonic_path//'/structure.dat')
  structure = StructureData(structure_file%lines())
  
  ! Read in large anharmonic supercell and its q-points.
  anharmonic_supercell_file = IFile(wd//'/anharmonic_supercell.dat')
  anharmonic_supercell = StructureData(anharmonic_supercell_file%lines())
  
  qpoints_file = IFile(wd//'/qpoints.dat')
  qpoints = QpointData(qpoints_file%sections())
  
  ! Read in normal modes.
  complex_modes_file = IFile(wd//'/complex_modes.dat')
  complex_modes = ComplexMode(complex_modes_file%sections())
  
  real_modes_file = IFile(wd//'/real_modes.dat')
  real_modes = RealMode(real_modes_file%sections())
  
  ! Read in anharmonic potential.
  potential_file = IFile(wd//'/potential.dat')
  if (potential_representation=='polynomial') then
    potential = PolynomialPotential(potential_file%lines())
  else
    call print_line( ERROR//': Unrecognised potential representation : '// &
                   & potential_representation)
    call err()
  endif
  
  ! --------------------------------------------------
  ! Initialise calculation handlers.
  ! --------------------------------------------------
  calculation_writer = CalculationWriter( working_directory = wd,        &
                                        & file_type         = file_type, &
                                        & seedname          = seedname   )
  
  calculation_runner = CalculationRunner(   &
     & working_directory = wd,              &
     & file_type         = file_type,       &
     & seedname          = seedname,        &
     & run_script        = run_script,      &
     & no_cores          = no_cores,        &
     & calculation_type  = calculation_type )
  
  calculation_reader = CalculationReader()
  
  ! --------------------------------------------------
  ! Re-calculate maximum_weighted_displacement.
  ! --------------------------------------------------
  maximum_weighted_displacement = maximum_displacement &
                              & * sqrt(minval(structure%atoms%mass()))
  
  ! --------------------------------------------------
  ! Calculate effective harmonic potential, from which initial harmonic
  !    states are constructed.
  ! --------------------------------------------------
  ! Calculate displacements before scaling by 1/sqrt(frequency).
  displacements =                                               &
     &   [(j,j=-no_single_mode_samples,no_single_mode_samples)] &
     & * maximum_weighted_displacement                          &
     & / no_single_mode_samples
  
  allocate( effective_frequencies(size(complex_modes)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(complex_modes)
    ! Scale displacement by 1/sqrt(frequency).
    scaled_displacements = displacements                             &
                       & * sqrt( frequency_of_max_displacement       &
                       &       / max( complex_modes(i)%frequency,    &
                       &              frequency_of_max_displacement) )
    
    ! Sample the model potential to find the effective frequency.
    effective_frequencies(i) = EffectiveFrequency( scaled_displacements, &
                                                 & complex_modes(i),     &
                                                 & real_modes,           &
                                                 & potential             )
  enddo
  
  ! --------------------------------------------------
  ! Write effective frequencies to file, q-point by q-point.
  ! --------------------------------------------------
  do i=1,size(qpoints)
    qpoint_modes = filter(complex_modes%qpoint_id==qpoints(i)%id)
    
    qpoint_dir = wd//'/qpoint_'//left_pad( qpoints(i)%id,          &
                                         & str(maxval(qpoints%id)) )
    call mkdir(qpoint_dir)
    
    ! --------------------------------------------------
    ! Validate the model potential by sampling the electronic structure at
    !    the same points as the model potential.
    ! --------------------------------------------------
    if (validate_potential) then
      ! Construct and write out supercell.
      supercell_matrix = construct_supercell_matrix(qpoints(i), structure)
      supercell = construct_supercell( structure,       &
                                     & supercell_matrix )
      supercell_file = OFile(qpoint_dir//'/structure.dat')
      call supercell_file%print_lines(supercell)
      
      do j=1,size(qpoint_modes)
        complex_mode = complex_modes(qpoint_modes(j))
        real_mode = real_modes(                         &
           & first(real_modes%id==min( complex_mode%id, &
           &                           complex_mode%paired_id)))
        mode_dir = qpoint_dir//                           &
           & '/complex_mode_'//left_pad( complex_mode%id, &
           &                             str(size(complex_modes)))
        call mkdir(mode_dir)
        allocate( sampled_energies(size(displacements)), &
                & sampled_forces(size(displacements)),   &
                & stat=ialloc); call err(ialloc)
        do k=1,size(sampled_energies)
          real_mode_displacement = RealModeDisplacement(                 &
             & [real_mode],                                              &
             & [effective_frequencies(qpoint_modes(j))%displacements(k)] )
          displacement = CartesianDisplacement( real_mode_displacement, &
                                              & supercell,              &
                                              & real_modes,             &
                                              & qpoints                 )
          displaced_structure = displace_structure(supercell,displacement)
          
          displacement_dir = mode_dir// &
             & '/displacement_'//left_pad(k,str(size(displacements)))
          call calculation_writer%write_calculation( displaced_structure, &
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
          sampled_forces(k) = sampled_force%force(real_mode)
        enddo
        sampled_energies = sampled_energies                                   &
                       & - potential%energy(                                  &
                       &     RealModeDisplacement([RealSingleDisplacement::]) )
        effective_frequencies(qpoint_modes(j))%sampled_energies = &
           & sampled_energies
        effective_frequencies(qpoint_modes(j))%sampled_forces = sampled_forces
        
        deallocate( sampled_energies, &
                  & sampled_forces,   &
                  & stat=ialloc); call err(ialloc)
      enddo
    endif
    
    qpoint_frequencies = effective_frequencies(qpoint_modes)
    effective_frequencies_file = &
       & OFile(qpoint_dir//'/effective_frequencies.dat')
    call effective_frequencies_file%print_lines( qpoint_frequencies, &
                                               & separating_line='')
  enddo
end subroutine
end module
