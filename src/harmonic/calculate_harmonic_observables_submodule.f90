submodule (caesar_calculate_harmonic_observables_module) caesar_calculate_harmonic_observables_submodule
  use caesar_harmonic_module
contains

module procedure calculate_harmonic_observables_mode
  output%mode_name = 'calculate_harmonic_observables'
  output%description = 'Calculates observables under the harmonic &
     &approximation: the phonon density of states and dispersion curve, &
     &and the energy, free energy and entropy per unit cell. Should be run &
     &after calculate_normal_modes.'
  output%keywords = [                                                         &
     & KeywordData( 'min_frequency',                                          &
     &              'min_frequency is the frequency below which modes will be &
     &ignored when calculating thermodynamic quantities. min_frequency should &
     &be given in Hartree.',                                                  &
     &              default_value='1e-8'),                                    &
     & KeywordData( 'calculate_dispersion',                                   &
     &              'calculate_dispersion specifies whether or not to &
     &calculate the phonon dispersion curve.',                                &
     &              default_value='true'),                                    &
     & KeywordData( 'path',                                                   &
     &              'path is the path through fractional reciprocal space &
     &which will be mapped by the phonon dispersion curve. The path should be &
     &specified as labels and q-points, separated by commas. Breaks in the &
     &path should be marked with pipes, "|". The Gamma-point should be &
     &labelled G.',                                                           &
     &              default_value='G 0.0 0.0 0.0, R 0.5 0.5 0.5, &
     &M 0.0 0.5 0.5, G 0.0 0.0 0.0, X 0.0 0.0 0.5 | R 0.5 0.5 0.5, &
     &X 0.0 0.0 0.5, M 0.0 0.5 0.5'),                                         &
     & KeywordData( 'no_path_points',                                         &
     &              'no_path_points is the number of q-points sampled along &
     &the entire dispersion curve path.',                                     &
     &              default_value='1000'),                                    &
     & KeywordData( 'calculate_dos_and_thermodynamics',                       &
     &              'calculate_dos_and_thermodynamics specifies whether or &
     &not to calculate the phonon density of states and thermodynamic &
     &variables such as the free energy and entropy.',                        &
     &              default_value='true'),                                    &
     & KeywordData( 'min_temperature',                                        &
     &              'min_temperature is the minimum temperature at which &
     &thermodynamic quantities are calculated. min_temperature should be &
     &given in Kelvin.'),                                                     &
     & KeywordData( 'max_temperature',                                        &
     &              'max_temperature is the maximum temperature at which &
     &thermodynamic quantities are calculated. min_temperature should be &
     &given in Kelvin.'),                                                     &
     & KeywordData( 'no_temperature_steps',                                   &
     &              'no_temperature_steps is the number of temperatures at &
     &which thermodynamic quantities are calculated.',                        &
     &              default_value='0'),                                       &
     & KeywordData( 'no_dos_samples',                                         &
     &              'no_dos_samples is the number of points in reciprocal &
     &space at which the normal modes are calculated when calculating the &
     &vibrational density of states.',                                        &
     &              default_value='10000')                                    ]
  output%main_subroutine => calculate_harmonic_observables_subroutine
end procedure

module procedure calculate_harmonic_observables_subroutine
  ! Inputs.
  type(RandomReal)      :: random_generator
  real(dp)              :: min_temperature
  real(dp)              :: max_temperature
  integer               :: no_temperature_steps
  real(dp)              :: min_frequency
  real(dp), allocatable :: thermal_energies(:)
  logical               :: calculate_dispersion
  type(String)          :: path
  integer               :: no_path_points
  logical               :: calculate_dos_and_thermodynamics
  integer               :: no_dos_samples
  
  ! Previous inputs.
  type(Dictionary) :: setup_harmonic_arguments
  type(String)     :: seedname
  
  ! Previously calculated data.
  type(StructureData)                :: structure
  type(StructureData)                :: large_supercell
  type(QpointData),      allocatable :: qpoints(:)
  type(DynamicalMatrix), allocatable :: dynamical_matrices(:)
  
  ! Hessian matrix and minimum image data.
  type(CartesianHessian)       :: hessian
  type(MinImages), allocatable :: min_images(:,:)
  
  ! Dispersion and density of states.
  type(PhononDispersion) :: phonon_dispersion
  type(PhononDos)        :: phonon_dos
  
  ! Files and directories.
  type(IFile)  :: structure_file
  type(IFile)  :: large_supercell_file
  type(IFile)  :: qpoints_file
  type(IFile)  :: dynamical_matrix_file
  type(String) :: output_dir
  type(OFile)  :: thermodynamic_file
  type(OFile)  :: logfile
  
  ! Temporary variables.
  integer :: i,ialloc
  
  ! --------------------------------------------------
  ! Read in arguments from user.
  ! --------------------------------------------------
  if (arguments%is_set('random_seed')) then
    random_generator = RandomReal(int(arguments%value('random_seed')))
  else
    random_generator = RandomReal()
  endif
  min_temperature = dble(arguments%value('min_temperature'))
  max_temperature = dble(arguments%value('max_temperature'))
  no_temperature_steps = int(arguments%value('no_temperature_steps'))
  min_frequency = dble(arguments%value('min_frequency'))
  calculate_dispersion = lgcl(arguments%value('calculate_dispersion'))
  if (calculate_dispersion) then
    path = arguments%value('path')
    no_path_points = int(arguments%value('no_path_points'))
  endif
  calculate_dos_and_thermodynamics = &
     & lgcl(arguments%value('calculate_dos_and_thermodynamics'))
  if (calculate_dos_and_thermodynamics) then
    no_dos_samples = int(arguments%value('no_dos_samples'))
  endif
  
  ! Check inputs.
  if ( (.not. calculate_dispersion) .and.       &
     & (.not. calculate_dos_and_thermodynamics) ) then
    call print_line(ERROR//': calculate_dispersion and &
       &calculate_dos_and_thermodynamics are both false. Nothing will be &
       &calculated.')
  endif
  
  if (calculate_dos_and_thermodynamics) then
    if (min_temperature<0) then
      call print_line(ERROR//': min_temperature must not be less than 0K.')
      call quit()
    elseif (max_temperature<min_temperature) then
      call print_line(ERROR//': max_temperature must not be less than &
         &min_temperature.')
      call quit()
    elseif (min_frequency<0) then
      call print_line(ERROR//': min_frequency must not be less than 0 &
         &Hartree.')
      call quit()
    endif
  endif
  
  ! --------------------------------------------------
  ! Read in previous arguments.
  ! --------------------------------------------------
  setup_harmonic_arguments = Dictionary(setup_harmonic_mode())
  call setup_harmonic_arguments%read_file('setup_harmonic.used_settings')
  seedname = setup_harmonic_arguments%value('seedname')
  
  ! --------------------------------------------------
  ! Read in previously calculated data.
  ! --------------------------------------------------
  structure_file = IFile('structure.dat')
  structure = StructureData(structure_file%lines())
  
  large_supercell_file = IFile('large_supercell.dat')
  large_supercell = StructureData(large_supercell_file%lines())
  
  qpoints_file = IFile('qpoints.dat')
  qpoints = QpointData(qpoints_file%sections())
  
  allocate( dynamical_matrices(size(qpoints)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(qpoints)
    dynamical_matrix_file = IFile(                   &
       & 'qpoint_'//left_pad(i,str(size(qpoints)))// &
       & '/dynamical_matrix.dat')
    dynamical_matrices(i) = DynamicalMatrix(dynamical_matrix_file%lines())
  enddo
  
  ! --------------------------------------------------
  ! Run calculations.
  ! --------------------------------------------------
  
  ! Make directory for harmonic observables.
  output_dir = 'harmonic_observables'
  call mkdir(output_dir)
  
  ! Open output files.
  logfile = OFile(output_dir//'/harmonic_observables_log.dat')
  
  ! Construct the Hessian matrix from dynamical matrices.
  hessian = reconstruct_hessian( large_supercell,    &
                               & qpoints,            &
                               & dynamical_matrices, &
                               & logfile             )
  
  ! Calculate minimum image distances.
  min_images = calculate_min_images(large_supercell)
  
  ! Calculate phonon dispersion curve.
  if (calculate_dispersion) then
    ! Generate harmonic phonon dispersion curve by interpolating between
    !    calculated q-points using Fourier interpolation.
    phonon_dispersion = PhononDispersion( supercell      = large_supercell, &
                                        & min_images     = min_images,      &
                                        & hessian        = hessian,         &
                                        & path_string    = path,            &
                                        & no_path_points = no_path_points,  &
                                        & logfile        = logfile          )
    
    call phonon_dispersion%write_files(output_dir, seedname, structure)
  endif
  
  ! Calculate phonon density of states and thermodynamic variables.
  if (calculate_dos_and_thermodynamics) then
    ! Generate thermal energies.
    ! thermal_energies(1)                    = min_temperature * kB.
    ! thermal_energies(no_temperature_steps) = max_temperature * kB.
    allocate( thermal_energies(no_temperature_steps), &
            & stat=ialloc); call err(ialloc)
    do i=1,no_temperature_steps
      thermal_energies(i) = KB_IN_AU                                   &
                        & * ( min_temperature*(no_temperature_steps-i) &
                        &   + max_temperature*(i-1) )                  &
                        & / (no_temperature_steps-1)
    enddo
    
    ! Generate harmonic phonon density of states, interpolating as above.
    phonon_dos = PhononDos( supercell        = large_supercell,  &
                          & min_images       = min_images,       &
                          & hessian          = hessian,          &
                          & thermal_energies = thermal_energies, &
                          & min_frequency    = min_frequency,    &
                          & no_dos_samples   = no_dos_samples,   &
                          & logfile          = logfile,          &
                          & random_generator = random_generator  )
    
    call phonon_dos%write_files(output_dir)
    
    thermodynamic_file = OFile(output_dir//'/thermodynamic_variables.dat')
    call thermodynamic_file%print_line( &
       &'kB * temperature (Hartree per cell) | &
       &Vibrational Energy per cell, U=<E>, (Hartree) | &
       &Vibrational Free Energy per cell, F=U-TS, (Hartree) | &
       &Vibrational Shannon Entropy per cell, S/k_B, (arb. units)')
    call thermodynamic_file%print_lines(phonon_dos%thermodynamic_data)
  endif
end procedure
end submodule
