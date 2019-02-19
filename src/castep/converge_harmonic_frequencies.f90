! ======================================================================
! Converges harmonic frequencies and free energies w/r/t
!    CASTEP cutoff energy and k-point spacing.
! Runs multiple harmonic caesar calculations with different
!    cutoff energies and k-point spacings.
! ======================================================================
module converge_harmonic_frequencies_module
  use common_module

  implicit none
  
  private
  
  public :: converge_harmonic_frequencies
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext
! ----------------------------------------------------------------------
function converge_harmonic_frequencies() result(output)
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'converge_harmonic_frequencies'
  output%description = 'Converges harmonic frequencies w/r/t cutoff energy &
     &and k-point spacing.'
  output%keywords = [                                                         &
     & KeywordData( 'file_type',                                              &
     &              'file_type is the file type which will be used for &
     &single-point energy calculations. Usual settings are: "castep", &
     &"caesar" and "xyz", but only CASTEP is supported by this convergence &
     &module at present',                                                     &
     &              default_value='castep'),                                  &
     & KeywordData( 'seedname',                                               &
     &              'seedname is the CASTEP seedname from which filenames &
     &are constructed.'),                                                     &
     & KeywordData( 'minimum_cutoff',                                         &
     &              'minimum_cutoff is the smallest cutoff energy which will &
     &be tested. minimum_cutoff must be an integer. This cut-off energy will &
     &be used when converging harmonic frequencies and free energies with &
     &respect to k-point sampling. minimum_cutoff should be given in &
     &Hartree.'),                                                             &
     & KeywordData( 'cutoff_step',                                            &
     &              'cutoff_step is the step between each cutoff energy &
     &which will be tested. cutoff_step must be an integer. cutoff_step &
     &should be given in Hartree.'),                                          &
     & KeywordData( 'maximum_cutoff',                                         &
     &              'maximum_cutoff is the cutoff energy at which calculation &
     &will be terminated if convergence has not been reached. maximum_cutoff &
     &must be an integer. maximum_cutoff should be given in Hartree.'),       &
     & KeywordData( 'maximum_k-point_spacing',                                &
     &              'maximum_k-point_spacing is the largest k-point spacing &
     &which will be tested. This is the k-point spacing that will be used &
     &when converging harmonic frequencies and free energies with respect to &
     &cut-off energy. maximum_k-point_spacing should be given in inverse &
     &Bohr.'),                                                                &
     & KeywordData( 'k-point_spacing_step',                                   &
     &              'k-point_spacing_step is the step between each k-point &
     &spacing which will be tested. k-point_spacing_step should be given in &
     &inverse Bohr.'),                                                        &
     & KeywordData( 'minimum_k-point_spacing',                                &
     &              'minimum_k-point_spacing is the k-point spacing at which &
     &calculation will be terminated if convergence has not been reached. &
     &minimum_k-point_spacing should be given in inverse Bohr.'),             &
     & KeywordData( 'converge_energies',                                      &
     &              'converge_energies determines whether convergence of the &
     &free energies is monitored as well as convergence of harmonic &
     &frequencies.',                                                          &
     &              default_value='false'),                                   &
     & KeywordData( 'convergence_mode',                                       &
     &              'convergence_mode determines which CASTEP parameters you &
     &wish to vary in your convergence testing. Options are "cutoff", &
     &"kpoints" and "both".',                                                 &
     &              default_value='both'),                                    &
     & KeywordData( 'freq_tolerance',                                         &
     &              'freq_tolerance is the tolerance used in the convergence &
     &of the harmonic frequencies.',                                          &
     &              default_value='1e-8'),                                    &
     & KeywordData( 'energy_tolerance',                                       &
     &              'energy_tolerance is the tolerance used in the &
     &convergence of the vibrational free energies. It is given in Hartree &
     &per primitive unit cell.',                                              &
     &              default_value='1e-6'),                                    &
     & KeywordData( 'convergence_count',                                      &
     &              'convergence_count is the number of consecutive &
     &frequencies and/or energies that must be within the defined tolerance &
     &of each other for convergence to be considered reached',                &
     &              default_value='3'),                                       &
     & KeywordData( 'q-point_grid',                                           &
     &              'q-point_grid is the number of q-points in each direction &
     &in a Monkhorst-Pack grid. This should be specified as three integers &
     &separated by spaces.'),                                                 &
     & KeywordData( 'symmetry_precision',                                     &
     &              'In order for a symmetry to be accepted, it must &
     &transform the position of every atom to within symmetry_precision of an &
     &atom of the same element. symmetry_precision should be given in Bohr.', &
     &              default_value='0.1'),                                     &
     & KeywordData( 'harmonic_displacement',                                  &
     &              'harmonic_displacement is the distance in bohr by which &
     &atoms will be displaced when mapping the harmonic Born-Oppenheimer &
     &surface.',                                                              &
     &              default_value='0.01'),                                    &
     & KeywordData( 'run_script',                                             &
     &              'run_script is the path to the script for running DFT. An &
     &example run script can be found in doc/input_files. This script should &
     &have executable file permissions and should run CASTEP',                &
     &              is_path=.true.),                                          &
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
     & KeywordData( 'acoustic_sum_rule',                                      &
     &              'acoustic_sum_rule specifies where the acoustic sum rule &
     &is applied. The options are "off", "forces", "matrices" and "both".',   &
     &              default_value='both'),                                    &
     & KeywordData( 'min_temperature',                                        &
     &              'min_temperature is the minimum temperature at which &
     &thermodynamic quantities are calculated. min_temperature should be &
     &given in Kelvin.',                                                      &
     &              default_value='0'),                                       &
     & KeywordData( 'max_temperature',                                        &
     &              'max_temperature is the maximum temperature at which &
     &thermodynamic quantities are calculated. min_temperature should be &
     &given in Kelvin.',                                                      &
     &              default_value='500'),                                     &
     & KeywordData( 'no_temperature_steps',                                   &
     &              'no_temperature_steps is the number of temperatures at &
     &which thermodynamic quantities are calculated.',                        &
     &              default_value='6'),                                       &
     & KeywordData( 'min_frequency',                                          &
     &              'min_frequency is the frequency below which modes will be &
     &ignored when calculating thermodynamic quantities. min_frequency should &
     &be given in Hartree.',                                                  &
     &              default_value='1e-8'),                                    &
     & KeywordData( 'path',                                                   &
     &              'path is the path through fractional reciprocal space &
     &which will be mapped by the phonon dispersion curve. The path should be &
     &specified as labels and q-points, separated by commas. The Gamma-point &
     &should be labelled G.',                                                 &
     &              default_value='G 0.0 0.0 0.0, R 0.5 0.5 0.5,              &
     &M 0.0 0.5 0.5, G 0.0 0.0 0.0, X 0.0 0.0 0.5'),                          &
     & KeywordData( 'no_dos_samples',                                         &
     &              'no_dos_samples is the number of points in reciprocal &
     &space at which the normal modes are calculated when calculating the &
     &vibrational density of states.',                                        &
     &              default_value='100000')                                   ]
  output%main_subroutine => converge_harmonic_frequencies_subroutine
end function

! ----------------------------------------------------------------------
! Main program
! ----------------------------------------------------------------------
subroutine converge_harmonic_frequencies_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! User inputs
  type(String) :: wd
  type(String) :: seedname
  type(String) :: file_type
  real(dp)     :: minimum_cutoff
  real(dp)     :: cutoff_step
  real(dp)     :: maximum_cutoff
  real(dp)     :: maximum_kpoint_spacing
  real(dp)     :: kpoint_spacing_step
  real(dp)     :: minimum_kpoint_spacing
  integer      :: grid(3)
  real(dp)     :: symmetry_precision
  real(dp)     :: harmonic_displacement
  integer      :: no_cores
  integer      :: no_nodes
  type(String) :: run_script
  logical      :: converge_energies
  logical      :: converge_cutoff
  logical      :: converge_kpoints
  type(String) :: acoustic_sum_rule
  real(dp)     :: min_temperature
  real(dp)     :: max_temperature
  integer      :: no_temperature_steps
  real(dp)     :: min_frequency
  type(String) :: path
  integer      :: no_dos_samples
  real(dp)     :: freq_tolerance
  real(dp)     :: energy_tolerance
  integer      :: convergence_count
  
  ! Structure.
  type(String)        :: input_filename
  type(StructureData) :: structure

  ! Energy cut-offs and kpoint spacings
  integer               :: no_cutoffs
  integer               :: no_kpoint_spacings
  real(dp), allocatable :: cutoffs(:)
  real(dp), allocatable :: kpoint_spacings(:)
  
  ! Variables to converge.
  type(RealVector), allocatable :: frequencies(:)
  type(RealVector), allocatable :: free_energies(:)
  
  ! Working variables
  integer :: no_qpoints
  logical :: frequencies_converged
  logical :: energies_converged
  
  ! Files and directories.
  type(String) :: dir
  type(OFile)  :: output_file
  
  ! Temporary variables
  integer :: i
  
  ! --------------------------------------------------
  ! Get settings from user, and check them.
  ! --------------------------------------------------
  wd = arguments%value('working_directory')
  file_type = arguments%value('file_type')
  seedname = arguments%value('seedname')
  minimum_cutoff = int(arguments%value('minimum_cutoff'))
  cutoff_step = int(arguments%value('cutoff_step'))
  maximum_cutoff = int(arguments%value('maximum_cutoff'))
  maximum_kpoint_spacing = dble(arguments%value('maximum_k-point_spacing'))
  kpoint_spacing_step = dble(arguments%value('k-point_spacing_step'))
  minimum_kpoint_spacing = dble(arguments%value('minimum_k-point_spacing'))
  grid = int(split_line(arguments%value('q-point_grid')))
  no_qpoints = grid(1)*grid(2)*grid(3)
  symmetry_precision = dble(arguments%value('symmetry_precision'))
  harmonic_displacement = dble(arguments%value('harmonic_displacement'))
  no_cores = int(arguments%value('no_cores'))
  no_nodes = int(arguments%value('no_nodes'))
  run_script = arguments%value('run_script')
  converge_energies = lgcl(arguments%value('converge_energies'))
  if (arguments%value('convergence_mode')=='cutoff') then
    converge_cutoff = .true.
    converge_kpoints = .false.
  elseif (arguments%value('convergence_mode')=='kpoints') then
    converge_cutoff = .false.
    converge_kpoints = .true.
  elseif (arguments%value('convergence_mode')=='both') then
    converge_cutoff = .true.
    converge_kpoints = .true.
  else
    call print_line(ERROR//': unexpected value for convergence_mode.')
    stop
  endif
  acoustic_sum_rule = arguments%value('acoustic_sum_rule')
  min_temperature = dble(arguments%value('min_temperature'))
  max_temperature = dble(arguments%value('max_temperature'))
  no_temperature_steps = int(arguments%value('no_temperature_steps'))
  min_frequency = dble(arguments%value('min_frequency'))
  path = arguments%value('path')
  no_dos_samples = int(arguments%value('no_dos_samples'))
  freq_tolerance = dble(arguments%value('freq_tolerance'))
  energy_tolerance = dble(arguments%value('energy_tolerance'))
  convergence_count = int(arguments%value('convergence_count'))
  
  if (maximum_cutoff<=minimum_cutoff) then
    call print_line(ERROR//': maximum_cutoff is smaller than minimum_cutoff.')
    stop
  elseif (maximum_kpoint_spacing<=minimum_kpoint_spacing) then
    call print_line(ERROR//': maximum_kpoints is smaller than &
       &minimum_kpoints.')
    stop
  elseif (file_type/='castep' .and. file_type/='quantum_espresso') then
    call print_line(ERROR//': castep and quantum_espresso are the only &
       &accepted file types for this mode.')
    stop
  elseif (max_temperature<=min_temperature) then
    call print_line(ERROR//': max_temperature is smaller than &
       &min_temperature.')
    stop
  elseif (convergence_count<1) then
    call print_line(ERROR//': convergence_count must be at least 1.')
    stop
  endif
  
  ! Read in structure.
  input_filename = make_input_filename(file_type, seedname)
  structure = input_file_to_StructureData(file_type, input_filename)
  
  ! Open output file.
  output_file = OFile('convergence.dat')
  
  ! Run cut-off energy convergence.
  if (converge_cutoff) then
    ! Initialise convergence checks and observable arrays.
    frequencies_converged = .false.
    energies_converged = .false.
    if (.not. converge_energies) then
     energies_converged = .true.
    endif
    
    frequencies = [RealVector::]
    free_energies = [RealVector::]
    
    ! Loop over cutoffs, running calculations at each.
    no_cutoffs = ceiling((maximum_cutoff-minimum_cutoff)/cutoff_step) &
             & + 1
    cutoffs = [( minimum_cutoff + (i-1)*cutoff_step, &
               & i=1,                                &
               & no_cutoffs                          )]
    do i=1,no_cutoffs
      dir = 'cutoff_'//                                                       &
         & left_pad(floor(cutoffs(i)),str(floor(cutoffs(no_cutoffs))))//'.'// &
         & left_pad(nint(1e2_dp*modulo(cutoffs(i),1.0_dp)), '  ')
      call mkdir(dir)
      
      frequencies = [ frequencies,                                   &
                    & calculate_frequencies( dir,                    &
                    &                        structure,              &
                    &                        maximum_kpoint_spacing, &
                    &                        cutoffs(i),             &
                    &                        seedname,               &
                    &                        no_qpoints,             &
                    &                        file_type,              &
                    &                        grid,                   &
                    &                        symmetry_precision,     &
                    &                        harmonic_displacement,  &
                    &                        run_script,             &
                    &                        no_cores,               &
                    &                        no_nodes,               &
                    &                        acoustic_sum_rule       ) ]
      
      if (i>convergence_count) then
        frequencies_converged = all(                                      &
           &   maximum_difference( frequencies(i),                        &
           &                       frequencies(i-convergence_count:i-1) ) &
           & < freq_tolerance                                             )
      endif

      if (converge_energies .and. .not. energies_converged) then
        free_energies = [ free_energies,                                 &
                        & calculate_free_energies( dir,                  &
                        &                          min_temperature,      &
                        &                          max_temperature,      &
                        &                          no_temperature_steps, &
                        &                          min_frequency,        &
                        &                          path,                 &
                        &                          no_dos_samples        ) ]

        if (i>convergence_count) then
          energies_converged = all(                                           &
             &   maximum_difference( free_energies(i),                        &
             &                       free_energies(i-convergence_count:i-1) ) &
             & < energy_tolerance                                             )
        endif
      endif
      
      call print_line('Cut-off energy: '//cutoffs(i)//' (Ha)')
      if (frequencies_converged .and. energies_converged) then
        call print_line('Convergence reached.')
        exit
      endif
    enddo
    
    if (.not. (frequencies_converged .or. energies_converged)) then
      call print_line('Convergence not reached.')
    endif
    
    ! Write output file.
    call output_file%print_line('Cutoff energy (Ha) | Mode frequencies (Ha)')
    do i=1,size(frequencies)
      call output_file%print_line(cutoffs(i)//' '//frequencies(i))
    enddo
    
    if (converge_energies) then
      call output_file%print_line('')
      call output_file%print_line('Cutoff energy (Ha) | Free energies &
         &(Ha per primitive cell)')
      do i=1,size(free_energies)
        call output_file%print_line(cutoffs(i)//' '//free_energies(i))
      enddo
    endif
  endif
  
  if (converge_cutoff .and. converge_kpoints) then
    call output_file%print_line('')
  endif
  
  ! Run k-points spacing convergence.
  if (converge_kpoints) then
    ! Initialise convergence checks and observable arrays.
    frequencies_converged = .false.
    energies_converged = .false.
    if (.not. converge_energies) then
     energies_converged = .true.
    endif
    
    frequencies = [RealVector::]
    free_energies = [RealVector::]
    
    ! Loop over spacings, running calculations at each.
    no_kpoint_spacings = ceiling(                            &
       &   (maximum_kpoint_spacing-minimum_kpoint_spacing)   &
       & / kpoint_spacing_step                             ) &
       & + 1
    kpoint_spacings = [( maximum_kpoint_spacing-(i-1)*kpoint_spacing_step, &
                       & i=1,                                              &
                       & no_kpoint_spacings                                )]
    do i=1,no_kpoint_spacings
      dir = 'kpoints_'//                                                      &
         & left_pad( floor(kpoint_spacings(i)),                               &
         &           str(floor(kpoint_spacings(no_kpoint_spacings))) )//'.'// &
         & left_pad( nint(1e4_dp*modulo(kpoint_spacings(i),1.0_dp)),          &
         &           '    '                                          )
      call mkdir(dir)
      
      frequencies = [ frequencies,                                  &
                    & calculate_frequencies( dir,                   &
                    &                        structure,             &
                    &                        kpoint_spacings(i),    &
                    &                        minimum_cutoff,        &
                    &                        seedname,              &
                    &                        no_qpoints,            &
                    &                        file_type,             &
                    &                        grid,                  &
                    &                        symmetry_precision,    &
                    &                        harmonic_displacement, &
                    &                        run_script,            &
                    &                        no_cores,              &
                    &                        no_nodes,              &
                    &                        acoustic_sum_rule      ) ]
      
      if (i>convergence_count) then
        frequencies_converged = all(                                      &
           &   maximum_difference( frequencies(i),                        &
           &                       frequencies(i-convergence_count:i-1) ) &
           & < freq_tolerance                                             )
      endif

      if (converge_energies .and. .not. energies_converged) then
        free_energies = [ free_energies,                                 &
                        & calculate_free_energies( dir,                  &
                        &                          min_temperature,      &
                        &                          max_temperature,      &
                        &                          no_temperature_steps, &
                        &                          min_frequency,        &
                        &                          path,                 &
                        &                          no_dos_samples        ) ]

        if (i>convergence_count) then
          energies_converged = all(                                           &
             &   maximum_difference( free_energies(i),                        &
             &                       free_energies(i-convergence_count:i-1) ) &
             & < energy_tolerance                                             )
        endif
      endif
      
      call print_line('K-point spacing: '//kpoint_spacings(i)//' (Bohr^-1)')
      if (frequencies_converged .and. energies_converged) then
        call print_line('Convergence reached.')
        exit
      endif
    enddo
    
    if (.not. (frequencies_converged .or. energies_converged)) then
      call print_line('Convergence not reached.')
    endif
    
    ! Write output file.
    call output_file%print_line('k-point spacing (Bohr^-1) | Mode frequencies &
       &(Ha)')
    do i=1,size(frequencies)
      call output_file%print_line(kpoint_spacings(i)//' '//frequencies(i))
    enddo
    
    if (converge_energies) then
      call output_file%print_line('')
      call output_file%print_line('k-point spacing (Bohr^-1) | Free energies &
         &(Ha per primitive cell)')
      do i=1,size(free_energies)
        call output_file%print_line(kpoint_spacings(i)//' '//free_energies(i))
      enddo
    endif
  endif
end subroutine

function calculate_frequencies(directory,structure,kpoint_spacing,cutoff,  &
   & seedname,no_qpoints,file_type,grid,symmetry_precision,                &
   & harmonic_displacement,run_script,no_cores,no_nodes,acoustic_sum_rule) &
   & result(output)
  implicit none
  
  type(String),        intent(in) :: directory
  type(StructureData), intent(in) :: structure
  real(dp),            intent(in) :: kpoint_spacing
  real(dp),            intent(in) :: cutoff
  type(String),        intent(in) :: seedname
  integer,             intent(in) :: no_qpoints
  type(String),        intent(in) :: file_type
  integer,             intent(in) :: grid(3)
  real(dp),            intent(in) :: symmetry_precision
  real(dp),            intent(in) :: harmonic_displacement
  type(String),        intent(in) :: run_script
  integer,             intent(in) :: no_cores
  integer,             intent(in) :: no_nodes
  type(String),        intent(in) :: acoustic_sum_rule
  type(RealVector)                :: output
  
  type(String) :: qpoint_dir
  type(IFile)  :: complex_modes_file
  
  type(ComplexMode), allocatable :: modes(:)
  
  integer :: i
  
  ! Write DFT files
  if (file_type=='castep') then
    call write_castep_files(directory, seedname, kpoint_spacing, cutoff)
  elseif (file_type=='quantum_espresso') then
    call write_qe_file(directory, seedname, structure, kpoint_spacing, cutoff)
  endif
  
  ! Call setup_harmonic.
  call call_caesar(                                              &
     & 'setup_harmonic -d '//directory                   //' '// &
     & '--file_type '//file_type                         //' '// &
     & '--seedname '//seedname                           //' '// &
     & '--q-point_grid '//grid                           //' '// &
     & '--symmetry_precision '//symmetry_precision       //' '// &
     & '--harmonic_displacement '//harmonic_displacement         )
  
  ! Call run_harmonic.
  call call_caesar( 'run_harmonic -d '//directory //' '// &
                  & '--run_script '//run_script   //' '// &
                  & '--no_cores '//no_cores       //' '// &
                  & '--no_nodes '//no_nodes       //' '// &
                  & '--exit_on_error true'        //' '// &
                  & '--repeat_calculations true'          )
  
  ! Call calculate_normal_modes.
  call call_caesar( 'calculate_normal_modes -d '//directory   //' '// &
                  & '--acoustic_sum_rule '//acoustic_sum_rule         )
  
  ! Read in normal modes.
  modes = [ComplexMode::]
  do i=1,no_qpoints
    qpoint_dir = directory//'/qpoint_'//left_pad(i, str(no_qpoints))
    complex_modes_file = IFile(qpoint_dir//'/complex_modes.dat')
    modes = [modes, ComplexMode(complex_modes_file%sections())]
  enddo
  
  output = vec(modes%frequency)
end function

function calculate_free_energies(directory,min_temperature,max_temperature, &
   & no_temperature_steps,min_frequency,path,no_dos_samples) result(output)
  implicit none
  
  type(String), intent(in) :: directory
  real(dp),     intent(in) :: min_temperature
  real(dp),     intent(in) :: max_temperature
  integer,      intent(in) :: no_temperature_steps
  real(dp),     intent(in) :: min_frequency
  type(String), intent(in) :: path
  integer,      intent(in) :: no_dos_samples
  type(RealVector)         :: output
  
  type(IFile)                          :: thermodynamics_file
  type(String),            allocatable :: lines(:)
  type(ThermodynamicData), allocatable :: thermodynamics(:)

  ! Call calculate_harmonic_observables.
  call call_caesar(                                            &
     & 'calculate_harmonic_observables -d '//directory //' '// &
     & '--min_temperature '//min_temperature           //' '// &
     & '--max_temperature '//max_temperature           //' '// &
     & '--no_temperature_steps '//no_temperature_steps //' '// &
     & '--min_frequency '//min_frequency               //' '// &
     & '--path '//path                                 //' '// &
     & '--no_dos_samples '//no_dos_samples                     )
  
  ! Read in thermodynamic data.
  thermodynamics_file = IFile(                                        &
     & directory//'/harmonic_observables/thermodynamic_variables.dat' )
  lines = thermodynamics_file%lines()
  thermodynamics = ThermodynamicData(lines(2:))
  
  output = vec(thermodynamics%free_energy)
end function

subroutine write_castep_files(directory,seedname,kpoint_spacing,cutoff)
  implicit none
  
  type(String), intent(in) :: directory
  type(String), intent(in) :: seedname
  real(dp),     intent(in) :: kpoint_spacing
  real(dp),     intent(in) :: cutoff
  
  type(IFile) :: input_cell_file
  type(OFile) :: output_cell_file
  type(IFile) :: input_param_file
  type(OFile) :: output_param_file
  
  logical :: line_replaced
  
  type(String), allocatable :: lines(:)
  type(String), allocatable :: line(:)
  
  integer :: i
  
  ! --------------------------------------------------
  ! Update cell file.
  ! --------------------------------------------------
  
  ! Read input cell file.
  input_cell_file = IFile(seedname//'.cell')
  lines = input_cell_file%lines()
  
  ! Update k-point spacing.
  line_replaced = .false.
  do i=1,size(lines)
    line = split_line(lower_case(lines(i)))
    if (size(line)>=1) then
      if ( line(1)=='kpoint_mp_spacing'  .or. &
         & line(1)=='kpoint_mp_grid'     .or. &
         & line(1)=='kpoints_mp_spacing' .or. &
         & line(1)=='kpoints_mp_grid'         ) then
        lines(i) = 'kpoints_mp_spacing '//kpoint_spacing/ANGSTROM_PER_BOHR
        line_replaced = .true.
        exit
      endif
    endif
  enddo
  
  if (.not. line_replaced) then
    lines = [lines, 'kpoints_mp_spacing '//kpoint_spacing/ANGSTROM_PER_BOHR]
  endif
  
  ! Write output cell file.
  output_cell_file = OFile(directory//'/'//seedname//'.cell')
  call output_cell_file%print_lines(lines)
  
  ! --------------------------------------------------
  ! Update param file.
  ! --------------------------------------------------
  
  ! Read input param file.
  input_param_file = IFile(seedname//'.param')
  lines = input_param_file%lines()
  
  ! Update electronic cut-off energy.
  line_replaced = .false.
  do i=1,size(lines)
    line = split_line(lower_case(lines(i)))
    if (line(1)=='cut_off_energy') then
      lines(i) = 'cut_off_energy '//cutoff*EV_PER_HARTREE
      line_replaced = .true.
      exit
    endif
  enddo
  
  if (.not. line_replaced) then
    lines = [lines, 'cut_off_energy '//cutoff*EV_PER_HARTREE]
  endif
  
  ! Write output param file.
  output_param_file = OFile(directory//'/'//seedname//'.param')
  call output_param_file%print_lines(lines)
end subroutine

subroutine write_qe_file(directory,seedname,structure,kpoint_spacing,cutoff)
  implicit none
  
  type(String),        intent(in) :: directory
  type(String),        intent(in) :: seedname
  type(StructureData), intent(in) :: structure
  real(dp),            intent(in) :: kpoint_spacing
  real(dp),            intent(in) :: cutoff
  
  type(IFile)       :: input_file
  type(QeInputFile) :: qe_input_file
  type(OFile)       :: output_file
  
  integer :: kpoint_grid(3)
  
  integer  :: ecutwfc_line
  integer  :: ecutrho_line
  real(dp) :: ecutwfc
  real(dp) :: ecutrho
  
  type(String), allocatable :: line(:)
  
  integer :: i
  
  ! Read input file.
  input_file = IFile(seedname//'.in')
  qe_input_file = QeInputFile(input_file%lines())
  
  ! Calculate k-point grid from spacing.
  kpoint_grid = calculate_kpoint_grid( kpoint_spacing,               &
                                     & dble(structure%recip_lattice) )
  
  line = split_line(qe_input_file%k_points(2))
  if (size(line)==3) then
    qe_input_file%k_points(2) = join(str(kpoint_grid))
  elseif (size(line)==6) then
    qe_input_file%k_points(2) = join([str(kpoint_grid), line(4:6)])
  else
    call print_line(ERROR//': k_points card has an unexpected number of &
       &entries.')
    stop
  endif
  
  ! Locate cutoff lines.
  ecutwfc_line = 0
  ecutrho_line = 0
  do i=1,size(qe_input_file%namelists)
    line = split_line(lower_case(qe_input_file%namelists(i)), delimiter='=')
    if (size(line)==2) then
      if (line(1)=='ecutwfc') then
        if (ecutwfc_line/=0) then
          call print_line(ERROR//': ecutwfc appears twice in input file.')
          stop
        endif
        ecutwfc_line = i
      elseif (line(1)=='ecutrho') then
        if (ecutrho_line/=0) then
          call print_line(ERROR//': ecutrho appears twice in input file.')
          stop
        endif
        ecutrho_line = i
      endif
    endif
  enddo
  
  if (ecutwfc_line==0) then
    call print_line(ERROR//': ecutwfc not present.')
  endif
  
  line = split_line(qe_input_file%namelists(ecutwfc_line), delimiter='=')
  ecutwfc = dble(line(2))
  qe_input_file%namelists(ecutwfc_line) = &
     & 'ecutwfc='//cutoff*EV_PER_HARTREE/EV_PER_RYDBERG
  
  if (ecutrho_line/=0) then
    line = split_line(qe_input_file%namelists(ecutrho_line), delimiter='=')
    ecutrho = dble(line(2))
    qe_input_file%namelists(ecutrho_line) = &
       & 'ecutrho='//ecutrho*cutoff/ecutwfc*EV_PER_HARTREE/EV_PER_RYDBERG
  endif
  
  ! Write output file.
  output_file = OFile(directory//'/'//seedname//'.in')
  call output_file%print_lines(qe_input_file)
end subroutine

! Find the maximum element-wise difference between two vectors.
impure elemental function maximum_difference(this,that) result(output)
  implicit none
  
  type(RealVector), intent(in) :: this
  type(RealVector), intent(in) :: that
  real(dp)                     :: output
  
  if (size(this)/=size(that)) then
    call print_line(CODE_ERROR//': Vectors of different lengths.')
    call err()
  endif
  
  output = maxval(abs(dble(this-that)))
end function
end module
