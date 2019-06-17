! ======================================================================
! Maps out the potential at finite displacements along multiple modes.
! ======================================================================
module map_potential_module
  use common_module
  
  use anharmonic_common_module
  use potentials_module
  implicit none
  
  private
  
  public :: startup_map_potential
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
subroutine startup_map_potential()
  implicit none
  
  type(CaesarMode) :: mode
  
  mode%mode_name = 'map_potential'
  mode%description = 'Maps out the potential at finite displacements along &
     & pairs of modes.'
  mode%keywords = [                                                           &
     & KeywordData( 'modes',                                                  &
     &              'modes is the IDs of the real modes along which &
     &displacements will be made. Modes should be given as a set of integer &
     &mode IDs separated by spaces. All pairs will be mapped, so e.g. &
     &modes="1 4 5" will result in the [1,4], [1,5] and [4,5] planes being &
     &mapped. If modes is not supplied then all mode pairs will be mapped.',  &
     &              is_optional=.true.),                                      &
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
  mode%main_subroutine => map_potential_subroutine
  
  call add_mode(mode)
end subroutine

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine map_potential_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Input arguments.
  integer, allocatable :: modes(:)
  integer              :: no_single_mode_samples
  logical              :: validate_potential
  type(String)         :: run_script
  integer              :: no_cores
  integer              :: no_nodes
  type(String)         :: run_script_data
  type(String)         :: calculation_type
  
  ! Arguments to setup_harmonic.
  type(Dictionary) :: setup_harmonic_arguments
  type(String)     :: seedname
  type(String)     :: file_type
  
  ! Anharmonic data.
  type(AnharmonicData)           :: anharmonic_data
  type(StructureData)            :: structure
  type(StructureData)            :: anharmonic_supercell
  type(QpointData),  allocatable :: qpoints(:)
  type(ComplexMode), allocatable :: complex_modes(:)
  type(RealMode),    allocatable :: real_modes(:)
  real(dp)                       :: frequency_of_max_displacement
  real(dp)                       :: maximum_weighted_displacement
  
  ! Selected real modes.
  type(RealMode), allocatable :: selected_modes(:)
  
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
  
  real(dp), allocatable :: harmonic_energy(:,:)
  real(dp), allocatable :: anharmonic_energy(:,:)
  real(dp), allocatable :: sampled_energy(:,:)
  
  ! Variables for validating potential.
  type(IntMatrix)             :: supercell_matrix
  type(StructureData)         :: supercell
  type(CartesianDisplacement) :: displacement
  type(StructureData)         :: displaced_structure
  type(ElectronicStructure)   :: electronic_structure
  
  ! Files and directories.
  type(IFile)  :: anharmonic_data_file
  type(IFile)  :: potential_file
  type(String) :: map_dir
  type(String) :: modes_dir
  type(OFile)  :: supercell_file
  type(String) :: displacement_dir
  type(OFile)  :: output_file
  
  ! Temporary variables.
  integer :: i,j,k,l,ialloc
  
  ! --------------------------------------------------
  ! Read in inputs and previously calculated data.
  ! --------------------------------------------------
  if (arguments%is_set('modes')) then
    modes = int(split_line(arguments%value('modes')))
  endif
  no_single_mode_samples = int(arguments%value('no_single_mode_samples'))
  validate_potential = lgcl(arguments%value('validate_potential'))
  run_script = arguments%value('run_script')
  no_cores = int(arguments%value('no_cores'))
  no_nodes = int(arguments%value('no_nodes'))
  run_script_data = arguments%value('run_script_data')
  calculation_type = arguments%value('calculation_type')
  
  ! Read in setup_harmonic arguments.
  setup_harmonic_arguments = Dictionary(CaesarMode('setup_harmonic'))
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
  complex_modes = anharmonic_data%complex_modes
  real_modes = anharmonic_data%real_modes
  frequency_of_max_displacement = anharmonic_data%frequency_of_max_displacement
  maximum_weighted_displacement = anharmonic_data%maximum_weighted_displacement
  
  ! Read in anharmonic potential.
  potential_file = IFile('potential.dat')
  potential = PotentialPointer(potential_file%lines())
  
  ! --------------------------------------------------
  ! Initialise calculation handlers.
  ! --------------------------------------------------
  calculation_writer = CalculationWriter( file_type = file_type, &
                                        & seedname  = seedname   )
  
  calculation_runner = CalculationRunner(      &
     & file_type           = file_type,        &
     & seedname            = seedname,         &
     & run_script          = run_script,       &
     & no_cores            = no_cores,         &
     & no_nodes            = no_nodes,         &
     & run_script_data     = run_script_data,  &
     & calculation_type    = calculation_type, &
     & use_forces          = .true.,           &
     & use_hessians        = .false.,          &
     & calculate_stress    = .false.,          &
     & exit_on_error       = .true.,           &
     & repeat_calculations = .true.            )

  
  calculation_reader = CalculationReader()
  
  ! --------------------------------------------------
  ! Calculate effective harmonic potential, from which initial harmonic
  !    states are constructed.
  ! --------------------------------------------------
  map_dir = 'mapped_potential'
  call mkdir(map_dir)
  
  ! Calculate displacements before scaling by 1/sqrt(frequency).
  displacements =                                               &
     &   [(j,j=-no_single_mode_samples,no_single_mode_samples)] &
     & * maximum_weighted_displacement                          &
     & / no_single_mode_samples
  
  ! Select the specified modes.
  if (arguments%is_set('modes')) then
    ! De-duplicate and sort modes.
    modes = modes(set(modes))
    modes = modes(sort(modes))
    
    allocate(selected_modes(size(modes)), stat=ialloc); call err(ialloc)
    do i=1,size(modes)
      selected_modes(i) = real_modes(first(real_modes%id==modes(i)))
    enddo
  else
    selected_modes = real_modes
  endif
  
  allocate( harmonic_energy(size(displacements), size(displacements)),   &
          & anharmonic_energy(size(displacements), size(displacements)), &
          & sampled_energy(size(displacements), size(displacements)),    &
          & stat=ialloc); call err(ialloc)
  do i=1,size(selected_modes)
    do j=i+1,size(selected_modes)
      modes_dir = map_dir//'/modes_'//                    &
         & left_pad( selected_modes(i)%id,                &
         &           str(maxval(real_modes%id)) ) //'_'// &
         & left_pad( selected_modes(j)%id,                &
         &           str(maxval(real_modes%id)) )
      call mkdir(modes_dir)
      
      ! Scale displacement by 1/sqrt(frequency).
      scaled_displacements_i = displacements                             &
                           & * sqrt( frequency_of_max_displacement       &
                           &       / max( selected_modes(i)%frequency,   &
                           &              frequency_of_max_displacement) )
      scaled_displacements_j = displacements                             &
                           & * sqrt( frequency_of_max_displacement       &
                           &       / max( selected_modes(j)%frequency,   &
                           &              frequency_of_max_displacement) )
      
      if (validate_potential) then
        qpoints_ij = select_qpoints( [selected_modes(i),selected_modes(j)], &
                                   & qpoints                                )
        supercell_matrix = construct_supercell_matrix( qpoints_ij, &
                                                     & structure   )
        supercell = construct_supercell( structure,       &
                                       & supercell_matrix )
        supercell_file = OFile(modes_dir//'/structure.dat')
        call supercell_file%print_lines(supercell)
      endif
      
      do k=1,size(scaled_displacements_i)
        do l=1,size(scaled_displacements_j)
          displacement_i = RealSingleDisplacement( selected_modes(i),        &
                                                 & scaled_displacements_i(k) )
          displacement_j = RealSingleDisplacement( selected_modes(j),        &
                                                 & scaled_displacements_j(l) )
          displacement_ij = RealModeDisplacement([ displacement_i, &
                                                 & displacement_j  ])
          anharmonic_energy(l,k) = potential%energy(displacement_ij)
          harmonic_energy(l,k) = 0.5_dp                            &
                             & * selected_modes(i)%spring_constant &
                             & * scaled_displacements_i(k)         &
                             & * scaled_displacements_i(k)         &
                             & + 0.5_dp                            &
                             & * selected_modes(j)%spring_constant &
                             & * scaled_displacements_j(l)         &
                             & * scaled_displacements_j(l)
          
          if (validate_potential) then
            displacement = CartesianDisplacement( displacement_ij, &
                                                & supercell,       &
                                                & selected_modes,  &
                                                & qpoints          )
            displaced_structure = displace_structure(supercell,displacement)
            
            displacement_dir = modes_dir                                 // &
                             & '/displacement_'                          // &
                             & left_pad(k,str(size(displacements)))//'_' // &
                             & left_pad(l,str(size(displacements)))
            call calculation_writer%write_calculation( displaced_structure, &
                                                     & displacement_dir     )
            call calculation_runner%run_calculation(displacement_dir)
            electronic_structure = calculation_reader%read_calculation( &
                                                     & displacement_dir )
            
            sampled_energy(l,k) = electronic_structure%energy() &
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
      
      call output_file%print_line('x-direction mode ID: '// &
                                     & selected_modes(i)%id )
      call output_file%print_line('y-direction mode ID: '// &
                                     & selected_modes(j)%id )
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
