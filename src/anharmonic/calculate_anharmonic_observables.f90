! ======================================================================
! Calculates, under the VSCF approximation:
!
! ======================================================================
! Should be run after calculate_states.
module calculate_anharmonic_observables_module
  use common_module
  
  use states_module
  use anharmonic_common_module
  use potentials_module
  
  use effective_frequency_module
  use vscf_thermodynamics_module
  implicit none
  
  private
  
  public :: calculate_anharmonic_observables
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
function calculate_anharmonic_observables() result(output)
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'calculate_anharmonic_observables'
  output%description = 'Calculates observables under the VSCF approximation. &
     &Should be run after calculate_vscf_potential.'
  output%keywords = [                                                         &
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
     & KeywordData( 'no_vscha_basis_states',                                  &
     &              'no_vscha_basis states is the number of states along each &
     &mode in the basis used for the VSCHA calculation.'),                    &
     & KeywordData( 'frequency_convergence',                                  &
     &              'frequency_convergence is the precision to which &
     &self-consistent frequencies will be converged when constructing the &
     &self-consistent anharmonic approximation to the VSCF potential. This &
     &should be given in Hartree.' ),                                         &
     & KeywordData( 'no_converged_calculations',                              &
     &              'no_converged_calculations is the number of consecutive &
     &calculations which must be converged to within frequency_convergence &
     &for the self-consistent anharmonic procedure to terminate.',            &
     &              default_value='5' ),                                      &
     & KeywordData( 'no_vscf_basis_states',                                   &
     &              'no_vscf_basis states is the number of states along each &
     &mode in the basis used for the VSCF calculation.'),                     &
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
     &              default_value='G 0.0 0.0 0.0, R 0.5 0.5 0.5, &
     &M 0.0 0.5 0.5, G 0.0 0.0 0.0, X 0.0 0.0 0.5'),                          &
     & KeywordData( 'no_dos_samples',                                         &
     &              'no_dos_samples is the number of points in reciprocal &
     &space at which the normal modes are calculated when calculating the &
     &vibrational density of states.',                                        &
     &              default_value='100000')                                   ]
  output%main_subroutine => calculate_anharmonic_observables_subroutine
end function

! ----------------------------------------------------------------------
! The main subroutine.
! ----------------------------------------------------------------------
subroutine calculate_anharmonic_observables_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Inputs.
  type(RandomReal)              :: random_generator
  integer                       :: random_generator_seed
  real(dp)                      :: min_temperature
  real(dp)                      :: max_temperature
  integer                       :: no_temperature_steps
  integer                       :: no_vscha_basis_states
  real(dp)                      :: frequency_convergence
  integer                       :: no_converged_calculations
  integer                       :: no_vscf_basis_states
  real(dp)                      :: min_frequency
  real(dp),         allocatable :: thermal_energies(:)
  type(String),     allocatable :: path_string(:)
  type(String),     allocatable :: path_labels(:)
  type(RealVector), allocatable :: path_qpoints(:)
  integer                       :: no_dos_samples
  
  ! Anharmonic data.
  type(AnharmonicData)                  :: anharmonic_data
  type(ComplexMode),        allocatable :: modes(:)
  type(QpointData),         allocatable :: qpoints(:)
  type(DegenerateSubspace), allocatable :: subspaces(:)
  type(StructureData)                   :: supercell
  
  ! Single-subspace potentials.
  type(PotentialPointer), allocatable :: subspace_potentials(:)
  
  ! Single-subspace basis. Only used to initalise frequency-finding.
  type(SubspaceBasis), allocatable :: subspace_bases(:)
  
  ! Finite-temperature effective harmonic frequencies.
  real(dp), allocatable :: initial_frequencies(:)
  real(dp), allocatable :: effective_frequencies(:,:)
  
  ! Modes and dynamical matrices at each q-point.
  type(DynamicalMatrix), allocatable :: dynamical_matrices(:)
  type(ComplexMode),     allocatable :: qpoint_modes(:)
  real(dp),              allocatable :: qpoint_frequencies(:)
  
  ! Hessian matrix and minimum image data.
  type(CartesianHessian)       :: hessian
  type(MinImages), allocatable :: min_images(:,:)
  
  ! Dispersion and density of states.
  type(PhononDispersion)               :: phonon_dispersion
  type(PhononDos)                      :: phonon_dos
  
  ! Thermodynamic data.
  type(VscfThermodynamics)             :: vscf_correction
  type(ThermodynamicData), allocatable :: vscha_thermodynamics(:)
  type(ThermodynamicData), allocatable :: vscha1_thermodynamics(:)
  type(ThermodynamicData), allocatable :: vscha2_thermodynamics(:)
  type(ThermodynamicData), allocatable :: uninterpolated_thermodynamics(:)
  type(ThermodynamicData), allocatable :: vscf1_thermodynamics(:)
  type(ThermodynamicData), allocatable :: vscf2_thermodynamics(:)
  type(ThermodynamicData), allocatable :: thermodynamics(:)
  
  ! VSCF correction data.
  real(dp)                :: vscf_frequency
  type(EnergySpectrum)    :: vscf_spectrum
  type(ThermodynamicData) :: vscha_thermo(2)
  
  ! Files and directories.
  type(IFile)  :: anharmonic_data_file
  type(IFile)  :: subspace_potentials_file
  type(IFile)  :: basis_file
  type(String) :: output_dir
  type(OFile)  :: logfile
  type(String) :: temperature_dir
  type(OFile)  :: temperature_file
  type(OFile)  :: dispersion_file
  type(OFile)  :: symmetry_points_file
  type(OFile)  :: sampled_qpoints_file
  type(OFile)  :: vscha_thermodynamic_file
  type(OFile)  :: vscha1_thermodynamic_file
  type(OFile)  :: vscha2_thermodynamic_file
  type(OFile)  :: uninterpolated_thermodynamic_file
  type(OFile)  :: vscf1_thermodynamic_file
  type(OFile)  :: vscf2_thermodynamic_file
  type(OFile)  :: thermodynamic_file
  type(OFile)  :: pdos_file
  
  ! Temporary variables.
  type(String), allocatable :: path_point(:)
  
  integer :: i,j,k,ialloc
  
  ! --------------------------------------------------
  ! Read in arguments from user.
  ! --------------------------------------------------
  if (arguments%is_set('random_seed')) then
    random_generator = RandomReal(int(arguments%value('random_seed')))
  else
    random_generator = RandomReal()
  endif
  random_generator_seed = random_generator%get_seed()
  min_temperature = dble(arguments%value('min_temperature'))
  max_temperature = dble(arguments%value('max_temperature'))
  no_temperature_steps = int(arguments%value('no_temperature_steps'))
  frequency_convergence = dble(arguments%value('frequency_convergence'))
  no_converged_calculations = int(arguments%value('no_converged_calculations'))
  min_frequency = dble(arguments%value('min_frequency'))
  no_vscha_basis_states = int(arguments%value('no_vscha_basis_states'))
  path_string = split_line(arguments%value('path'), ',')
  no_dos_samples = int(arguments%value('no_dos_samples'))
  no_vscf_basis_states = int(arguments%value('no_vscf_basis_states'))
  
  ! Check inputs.
  if (min_temperature<0) then
    call print_line(ERROR//': min_temperature must not be less than 0K.')
    stop
  elseif (max_temperature<min_temperature) then
    call print_line(ERROR//': max_temperature must not be less than &
       &min_temperature.')
    stop
  elseif (min_frequency<0) then
    call print_line(ERROR//': min_frequency must not be less than 0 Hartree.')
    stop
  endif
  
  ! Generate thermal energies.
  ! thermal_energies(1)                    = min_temperature * kB.
  ! thermal_energies(no_temperature_steps) = max_temperature * kB.
  if (no_temperature_steps==1) then
    thermal_energies = [KB_IN_AU*min_temperature]
  else
    allocate( thermal_energies(no_temperature_steps), &
            & stat=ialloc); call err(ialloc)
    do i=1,no_temperature_steps
      thermal_energies(i) = KB_IN_AU                                   &
                        & * ( min_temperature*(no_temperature_steps-i) &
                        &   + max_temperature*(i-1) )                  &
                        & / (no_temperature_steps-1.0_dp)
    enddo
  endif
  
  ! Generate path for dispersion calculation.
  allocate( path_qpoints(size(path_string)), &
          & path_labels(size(path_string)),  &
          & stat=ialloc); call err(ialloc)
  do i=1,size(path_string)
    path_point = split_line(path_string(i))
    path_labels(i) = path_point(1)
    path_qpoints(i) = dble(path_point(2:4))
  enddo
  
  ! Read in anharmonic data.
  anharmonic_data_file = IFile('anharmonic_data.dat')
  anharmonic_data = AnharmonicData(anharmonic_data_file%lines())
  
  modes = anharmonic_data%complex_modes
  qpoints = anharmonic_data%qpoints
  subspaces = anharmonic_data%degenerate_subspaces
  supercell = anharmonic_data%anharmonic_supercell
  
  ! Read in subspace potentials.
  subspace_potentials_file = IFile('subspace_potentials.dat')
  subspace_potentials = PotentialPointer(                                &
     & subspace_potentials_file%sections(separating_line=repeat('=',70)) )
  
  ! Read in subspace bases.
  basis_file = IFile('basis.dat')
  subspace_bases = SubspaceBasis(basis_file%sections())
  
  ! Make output directory, write out temperatures, and open logfile.
  output_dir = 'anharmonic_observables'
  call mkdir(output_dir)
  
  temperature_file = OFile(output_dir//'/temperatures.dat')
  call temperature_file%print_line('Thermal energies, KbT, (Ha):')
  call temperature_file%print_lines(thermal_energies)
  
  logfile = OFile(output_dir//'/anharmonic_observables_log.dat')
  
  ! Construct minimum image data.
  min_images = calculate_min_images(supercell)
  
  ! --------------------------------------------------
  ! Generate observables for each temperature in turn.
  ! --------------------------------------------------
  ! Initialise frequencies to the frequencies which minimise the energy
  !    of the harmonic ground state.
  initial_frequencies = subspace_bases%frequency
  allocate( effective_frequencies( size(subspace_potentials),      &
          &                        size(thermal_energies)),        &
          & dynamical_matrices(size(qpoints)),                     &
          & vscha_thermodynamics(size(thermal_energies)),          &
          & vscha1_thermodynamics(size(thermal_energies)),          &
          & vscha2_thermodynamics(size(thermal_energies)),          &
          & uninterpolated_thermodynamics(size(thermal_energies)), &
          & vscf1_thermodynamics(size(thermal_energies)),           &
          & vscf2_thermodynamics(size(thermal_energies)),           &
          & thermodynamics(size(thermal_energies)),                &
          & stat=ialloc); call err(ialloc)
  do i=1,size(thermal_energies)
    ! Make temperature directory.
    temperature_dir = output_dir//'/temperature_'// &
       & left_pad(i,str(size(thermal_energies)))
    call mkdir(temperature_dir)
    
    ! Calculate effective frequencies for each subspace.
    call print_line('')
    call print_line('Thermal energy '//i//' of '//size(thermal_energies)// &
       & ': KT = '//thermal_energies(i)//' Ha')
    do j=1,size(subspace_potentials)
      effective_frequencies(j,i) = calculate_effective_frequency( &
                                      & subspace_potentials(j),   &
                                      & subspaces(j),             &
                                      & anharmonic_data,          &
                                      & thermal_energies(i),      &
                                      & initial_frequencies(j),   &
                                      & no_vscha_basis_states,    &
                                      & frequency_convergence,    &
                                      & no_converged_calculations )
    enddo
    
    call print_line('Self-consistent harmonic frequencies calculated.')
    call print_line(effective_frequencies(:,i))
    
    ! The starting point for calculating the effective frequencies at the
    !    next temperature is the frequencies at this temperature.
    initial_frequencies = effective_frequencies(:,i)
    
    ! Assemble the dynamical matrices at each q-point from the complex modes
    !    and the effective frequencies.
    do j=1,size(qpoints)
      qpoint_modes = modes(filter(modes%qpoint_id==qpoints(j)%id))
      qpoint_frequencies = [(                                              &
         & effective_frequencies(                                          &
         &    first(subspaces%id==qpoint_modes(k)%subspace_id),            &
         &    i                                                 ),         &
         & k=1,                                                            &
         & size(qpoint_modes)                                              )]
      dynamical_matrices(j) = DynamicalMatrix(qpoint_modes,qpoint_frequencies)
    enddo
    
    ! Construct Hessian from dynamical matrices.
    hessian = reconstruct_hessian( supercell,          &
                                 & qpoints,            &
                                 & dynamical_matrices, &
                                 & logfile             )
    
    ! Generate self-consistent phonon dispersion curve by interpolating between
    !    calculated q-points using Fourier interpolation.
    phonon_dispersion = PhononDispersion( supercell,    &
                                        & min_images,   &
                                        & hessian,      &
                                        & path_labels,  &
                                        & path_qpoints, &
                                        & logfile       )
    
    symmetry_points_file = OFile(temperature_dir//'/high_symmetry_points.dat')
    call symmetry_points_file%print_lines( phonon_dispersion%path, &
                                         & separating_line=''      )
    
    dispersion_file = OFile(temperature_dir//'/phonon_dispersion_curve.dat')
    call dispersion_file%print_lines(phonon_dispersion%frequencies)
    
    ! Generate self-consistent phonon density of states,
    !    interpolating as above.
    ! Re-seed the random generator each time, so that the results at different
    !    temperatures can be compared free from random noise.
    random_generator = RandomReal(random_generator_seed)
    phonon_dos = PhononDos( supercell,             &
                          & min_images,            &
                          & hessian,               &
                          & [thermal_energies(i)], &
                          & min_frequency,         &
                          & no_dos_samples,        &
                          & logfile,               &
                          & random_generator       )
    
    ! Write out dos and dispersion at this temperature.
    sampled_qpoints_file = OFile(temperature_dir//'/sampled_qpoints.dat')
    call sampled_qpoints_file%print_line('q-point (x,y,z) | &
                                         &number of frequencies ignored')
    call sampled_qpoints_file%print_lines(phonon_dos%qpoints)
    
    pdos_file = OFile(temperature_dir//'/phonon_density_of_states.dat')
    call pdos_file%print_lines(phonon_dos%pdos)
    
    vscha_thermodynamics(i) = phonon_dos%thermodynamic_data(1)
  enddo
  
  call print_line('')
  call print_line('VSCHA calculation complete. Running VSCF calculation.')
  
  ! Calculate VSCF correction to thermodynamic quantities.
  vscf1_thermodynamics = vscha_thermodynamics * 0
  vscf2_thermodynamics = vscha_thermodynamics * 0
  vscha1_thermodynamics = vscha_thermodynamics * 0
  vscha2_thermodynamics = vscha_thermodynamics * 0
  uninterpolated_thermodynamics = vscha_thermodynamics * 0
  thermodynamics = vscha_thermodynamics * 0
  do i=1,size(subspaces)
    ! Calculate the VSCF spectrum,
    !    using the VSCHA basis at the lowest temperature.
    call print_line('Calculating VSCF spectrum in subspace '//subspaces(i)%id)
    vscf_frequency = effective_frequencies(i,1)
    vscf_spectrum = calculate_vscf_spectrum( vscf_frequency,         &
                                           & subspace_potentials(i), &
                                           & subspaces(i),           &
                                           & anharmonic_data,        &
                                           & no_vscf_basis_states    )
    call print_line('VSCF states span '//     &
       &  maxval(vscf_spectrum%vscf_energies) &
       & -minval(vscf_spectrum%vscf_energies) //' (Ha)')
    do j=1,size(thermal_energies)
      vscf1_thermodynamics(j) = vscf1_thermodynamics(j)                      &
                           & + ThermodynamicData( thermal_energies(j),     &
                           &                      vscf_spectrum%vscf_energies )
      vscf2_thermodynamics(j) = vscf2_thermodynamics(j)                      &
                           & + ThermodynamicData( thermal_energies(j),     &
                           &                      vscf_spectrum%vscha_energies )
      
      vscha_thermo = calculate_vscha_thermodynamics( &
                       & effective_frequencies(i,j), &
                       & subspace_potentials(i),     &
                       & subspaces(i),               &
                       & anharmonic_data,            &
                       & no_vscha_basis_states,      &
                       & thermal_energies(j)         )
      vscha1_thermodynamics(j) = vscha1_thermodynamics(j) + vscha_thermo(1)
      vscha2_thermodynamics(j) = vscha2_thermodynamics(j) + vscha_thermo(2)
      !uninterpolated_thermodynamics(j) =                     &
      !   &   uninterpolated_thermodynamics(j)                &
      !   & + ThermodynamicData( thermal_energies(j),         &
      !   &                      vscha_spectrum(1)%energies ) &
      !   & / supercell%sc_size
      !thermodynamics(j) =                                    &
      !   &   thermodynamics(j)                               &
      !   & + ThermodynamicData( thermal_energies(j),         &
      !   &                      vscha_spectrum(2)%energies ) &
      !   & / supercell%sc_size
    enddo
  enddo
  
  vscf1_thermodynamics = vscf1_thermodynamics / supercell%sc_size
  vscf2_thermodynamics = vscf2_thermodynamics / supercell%sc_size
  
  ! --------------------------------------------------
  ! Write out thermodynamic results.
  ! --------------------------------------------------
  vscha1_thermodynamic_file = OFile(                    &
     & output_dir//'/vscha_thermodynamic_variables.dat' )
  call vscha1_thermodynamic_file%print_line( &
     &'kB * temperature (Hartree per cell) | &
     &Vibrational Energy per cell, U=<E>, (Hartree) | &
     &Vibrational Free Energy per cell, F=U-TS, (Hartree) | &
     &Vibrational Shannon Entropy per cell, S/k_B, (arb. units)')
  call vscha1_thermodynamic_file%print_lines(vscha1_thermodynamics)
  
  vscha2_thermodynamic_file = OFile(                         &
     & output_dir//'/vscha_vscf_thermodynamic_variables.dat' )
  call vscha2_thermodynamic_file%print_line( &
     &'kB * temperature (Hartree per cell) | &
     &Vibrational Energy per cell, U=<E>, (Hartree) | &
     &Vibrational Free Energy per cell, F=U-TS, (Hartree) | &
     &Vibrational Shannon Entropy per cell, S/k_B, (arb. units)')
  call vscha2_thermodynamic_file%print_lines(vscha2_thermodynamics)
  
  ! Just VSCHA with Fourier interpolation.
  vscha_thermodynamic_file = OFile(                                  &
     & output_dir//'/interpolated_vscha_thermodynamic_variables.dat' )
  call vscha_thermodynamic_file%print_line( &
     &'kB * temperature (Hartree per cell) | &
     &Vibrational Energy per cell, U=<E>, (Hartree) | &
     &Vibrational Free Energy per cell, F=U-TS, (Hartree) | &
     &Vibrational Shannon Entropy per cell, S/k_B, (arb. units)')
  call vscha_thermodynamic_file%print_lines(vscha_thermodynamics)
  
  ! Just VSCF without Fourier interpolation.
  vscf1_thermodynamic_file = OFile(                     &
     & output_dir//'/vscf_thermodynamic_variables.dat' )
  call vscf1_thermodynamic_file%print_line( &
     &'kB * temperature (Hartree per cell) | &
     &Vibrational Energy per cell, U=<E>, (Hartree) | &
     &Vibrational Free Energy per cell, F=U-TS, (Hartree) | &
     &Vibrational Shannon Entropy per cell, S/k_B, (arb. units)')
  call vscf1_thermodynamic_file%print_lines(vscf1_thermodynamics)
  
  vscf2_thermodynamic_file = OFile(                                  &
     & output_dir//'/limited_vscha_vscf_thermodynamic_variables.dat' )
  call vscf2_thermodynamic_file%print_line( &
     &'kB * temperature (Hartree per cell) | &
     &Vibrational Energy per cell, U=<E>, (Hartree) | &
     &Vibrational Free Energy per cell, F=U-TS, (Hartree) | &
     &Vibrational Shannon Entropy per cell, S/k_B, (arb. units)')
  call vscf2_thermodynamic_file%print_lines(vscf2_thermodynamics)
  
  !! Interpolated VSCHA + uninterpolated (VSCF-VSCHA).
  !thermodynamic_file = OFile(output_dir//'/thermodynamic_variables.dat')
  !call thermodynamic_file%print_line( &
  !   &'kB * temperature (Hartree per cell) | &
  !   &Vibrational Energy per cell, U=<E>, (Hartree) | &
  !   &Vibrational Free Energy per cell, F=U-TS, (Hartree) | &
  !   &Vibrational Shannon Entropy per cell, S/k_B, (arb. units)')
  !call thermodynamic_file%print_lines(vscf_thermodynamics)
end subroutine
end module