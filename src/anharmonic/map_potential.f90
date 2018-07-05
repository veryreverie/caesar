! ======================================================================
! Maps out the potential at finite displacements along multiple modes.
! ======================================================================
module map_potential_module
  use common_module
  
  use anharmonic_common_module
  use polynomial_module
  
  use setup_harmonic_module
  use calculate_normal_modes_module
  
  use setup_anharmonic_module
  implicit none
  
  private
  
  public :: map_potential
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
function map_potential() result(output)
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'map_potential'
  output%description = 'Maps out the potential at finite displacements along &
     & multiple modes.'
  output%keywords = [                                                         &
  & KeywordData( 'no_single_mode_samples',                                    &
  &              'no_single_mode_samples is the number of points (either side &
  &of zero) along each mode at which the anharmonic potential will be sampled &
  &when determining the effective frequency with which the harmonic basis &
  &along that mode will be constructed.',                                     &
  &              default_value='100'),                                        &
  & KeywordData( 'validate_potential',                                        &
  &              'validate_potential specifies that the anharmonic potential &
  &should be verified against fresh electronic structure calculations when &
  &sampling. Depending on the electronic structure method, this is likely to &
  &be very computationally intensive.',                                       &
  &              default_value='false'),                                      &
  & KeywordData( 'run_script',                                                &
  &              'run_script is the path to the script for running DFT. An &
  &example run script can be found in doc/input_files.',                      &
  &               is_path=.true.),                                            &
  & KeywordData( 'no_cores',                                                  &
  &              'no_cores is the number of cores on which DFT will be run. &
  &This is passed to the specified run script.',                              &
  &               default_value='1') ]
  output%main_subroutine => map_potential_subroutine
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine map_potential_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Working directory.
  type(String) :: wd
  
  ! Input arguments.
  integer      :: no_single_mode_samples
  logical      :: validate_potential
  type(String) :: run_script
  integer      :: no_cores
  
  ! Arguments to setup_anharmonic.
  type(Dictionary) :: setup_anharmonic_arguments
  type(String)     :: harmonic_path
  type(String)     :: potential_representation
  real(dp)         :: maximum_displacement
  real(dp)         :: frequency_of_max_displacement
  
  ! Arguments to calculate_normal_modes.
  type(Dictionary) :: calculate_normal_modes_arguments
  type(String)     :: calculation_type
  
  ! Arguments to setup_harmonic.
  type(Dictionary) :: setup_harmonic_arguments
  type(String)     :: seedname
  type(String)     :: file_type
  real(dp)         :: symmetry_precision
  
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
  real(dp),                 allocatable :: scaled_displacements_i(:)
  real(dp),                 allocatable :: scaled_displacements_j(:)
  type(RealSingleDisplacement)          :: displacement_i
  type(RealSingleDisplacement)          :: displacement_j
  type(RealModeDisplacement)            :: displacement_ij
  type(QpointData),         allocatable :: qpoints_ij(:)
  type(EffectiveFrequency), allocatable :: effective_frequencies(:)
  integer,                  allocatable :: qpoint_modes(:)
  type(EffectiveFrequency), allocatable :: qpoint_frequencies(:)
  
  real(dp), allocatable :: harmonic_energy(:,:)
  real(dp), allocatable :: anharmonic_energy(:,:)
  real(dp), allocatable :: sampled_energy(:,:)
  
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
  type(IFile)  :: qpoints_file
  type(IFile)  :: complex_modes_file
  type(IFile)  :: real_modes_file
  type(IFile)  :: potential_file
  type(String) :: map_dir
  type(String) :: modes_dir
  type(String) :: displacement_dir
  type(OFile)  :: output_file
  
  ! Temporary variables.
  integer :: i,j,k,l,ialloc
  
  ! --------------------------------------------------
  ! Read in inputs and previously calculated data.
  ! --------------------------------------------------
  
  wd = arguments%value('working_directory')
  
  no_single_mode_samples = int(arguments%value('no_single_mode_samples'))
  validate_potential = lgcl(arguments%value('validate_potential'))
  run_script = arguments%value('run_script')
  no_cores = int(arguments%value('no_cores'))
  
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
  calculation_type = calculate_normal_modes_arguments%value('calculation_type')
  
  ! Read in setup_harmonic arguments.
  setup_harmonic_arguments = Dictionary(setup_harmonic())
  call setup_harmonic_arguments%read_file( &
     & harmonic_path//'/setup_harmonic.used_settings')
  seedname = setup_harmonic_arguments%value('seedname')
  file_type = setup_harmonic_arguments%value('file_type')
  symmetry_precision = &
     & dble(setup_harmonic_arguments%value('symmetry_precision'))
  
  ! Read in structure.
  structure = read_structure_file( harmonic_path//'/structure.dat', &
                                 & symmetry_precision)
  
  ! Read in large anharmonic supercell and its q-points.
  anharmonic_supercell = read_structure_file( &
           & wd//'/anharmonic_supercell.dat', &
           & symmetry_precision,              &
           & calculate_symmetry=.false.)
  
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
    potential = PolynomialPotential(StringArray(potential_file%lines()))
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
  
  calculation_runner = CalculationRunner( working_directory = wd,         &
                                        & file_type         = file_type,  &
                                        & seedname          = seedname,   &
                                        & run_script        = run_script, &
                                        & no_cores          = no_cores    )
  
  calculation_reader = CalculationReader(      &
     & working_directory  = wd,                &
     & file_type          = file_type,         &
     & seedname           = seedname,          &
     & calculation_type   = calculation_type,  &
     & symmetry_precision = symmetry_precision )
  
  ! --------------------------------------------------
  ! Re-calculate maximum_weighted_displacement.
  ! --------------------------------------------------
  maximum_weighted_displacement = maximum_displacement &
                              & * sqrt(minval(structure%atoms%mass()))
  
  ! --------------------------------------------------
  ! Calculate effective harmonic potential, from which initial harmonic
  !    states are constructed.
  ! --------------------------------------------------
  map_dir = wd//'/mapped_potential'
  call mkdir(map_dir)
  
  ! Calculate displacements before scaling by 1/sqrt(frequency).
  displacements =                                               &
     &   [(j,j=-no_single_mode_samples,no_single_mode_samples)] &
     & * maximum_weighted_displacement                          &
     & / no_single_mode_samples
  
  allocate( harmonic_energy(size(displacements), size(displacements)),   &
          & anharmonic_energy(size(displacements), size(displacements)), &
          & sampled_energy(size(displacements), size(displacements)),    &
          & stat=ialloc); call err(ialloc)
  do i=1,size(real_modes)
    do j=i+1,size(real_modes)
      modes_dir = map_dir//'/modes_'//                     &
                & left_pad(i,str(size(real_modes)))//'_'// &
                & left_pad(j,str(size(real_modes)))
      call mkdir(modes_dir)
      
      ! Scale displacement by 1/sqrt(frequency).
      scaled_displacements_i = displacements                             &
                           & * sqrt( frequency_of_max_displacement       &
                           &       / max( complex_modes(i)%frequency,    &
                           &              frequency_of_max_displacement) )
      scaled_displacements_j = displacements                             &
                           & * sqrt( frequency_of_max_displacement       &
                           &       / max( complex_modes(j)%frequency,    &
                           &              frequency_of_max_displacement) )
      
      if (validate_potential) then
        qpoints_ij = select_qpoints([real_modes(i),real_modes(j)],qpoints)
        supercell_matrix = construct_supercell_matrix( qpoints_ij, &
                                                     & structure   )
        supercell = construct_supercell( structure,                 &
                                       & supercell_matrix,          &
                                       & symmetry_precision,        &
                                       & calculate_symmetry=.false. )
        call write_structure_file(supercell, modes_dir//'/structure.dat')
      endif
      
      do k=1,size(scaled_displacements_i)
        do l=1,size(scaled_displacements_j)
          displacement_i = RealSingleDisplacement( real_modes(i),            &
                                                 & scaled_displacements_i(k) )
          displacement_j = RealSingleDisplacement( real_modes(j),            &
                                                 & scaled_displacements_j(l) )
          displacement_ij = RealModeDisplacement([ displacement_i, &
                                                 & displacement_j  ])
          anharmonic_energy(l,k) = potential%energy(displacement_ij)
          harmonic_energy(l,k) = 0.5_dp                        &
                             & * real_modes(i)%spring_constant &
                             & * scaled_displacements_i(k)     &
                             & * scaled_displacements_i(k)     &
                             & + 0.5_dp                        &
                             & * real_modes(j)%spring_constant &
                             & * scaled_displacements_j(l)     &
                             & * scaled_displacements_j(l)
          
          if (validate_potential) then
            displacement = CartesianDisplacement( displacement_ij, &
                                                & supercell,       &
                                                & real_modes,      &
                                                & qpoints          )
            displaced_structure = displace_structure(supercell,displacement)
            
            displacement_dir = modes_dir//'/displacement_'//               &
                             & left_pad(k,str(size(displacements)))//'_'// &
                             & left_pad(l,str(size(displacements)))
            call calculation_writer%write_calculation( displaced_structure, &
                                                     & displacement_dir     )
            call calculation_runner%run_calculation(displacement_dir)
            electronic_structure = calculation_reader%read_calculation( &
                                                     & displacement_dir )
            
            sampled_energy(l,k) = electronic_structure%energy &
                              & / supercell%sc_size
          endif
        enddo
      enddo
      
      anharmonic_energy = anharmonic_energy                            &
                      & - anharmonic_energy( no_single_mode_samples+1, &
                      &                      no_single_mode_samples+1  )
      sampled_energy = sampled_energy                            &
                   & - sampled_energy( no_single_mode_samples+1, &
                   &                   no_single_mode_samples+1  )
      
      output_file = OFile(modes_dir//'/potential.dat')
      
      call output_file%print_line('x-direction mode ID: '//real_modes(i)%id)
      call output_file%print_line('y-direction mode ID: '//real_modes(j)%id)
      call output_file%print_line('x-direction displacements:')
      call output_file%print_line(scaled_displacements_i)
      call output_file%print_line('y-direction displacements:')
      call output_file%print_line(scaled_displacements_j)
      call output_file%print_line('Harmonic energies:')
      call output_file%print_lines(mat(harmonic_energy))
      call output_file%print_line('Anharmonic energies:')
      call output_file%print_lines(mat(anharmonic_energy))
      if (validate_potential) then
        call output_file%print_line('Sampled energies:')
        call output_file%print_lines(mat(sampled_energy))
      endif
    enddo
  enddo
end subroutine
end module
