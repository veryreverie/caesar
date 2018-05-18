! ======================================================================
! Calculates, under the harmonic approximation:
!    - The phonon dispersion curve along a specified path.
!    - The phonon density of states.
!    - The energy, free energy and entropy per unit cell.
! ======================================================================
! Should be run after calculate_normal_modes.
! Not required for anharmonic calculations.
module calculate_harmonic_observables_module
  use common_module
  
  use force_constants_module
  use dynamical_matrix_module
  use min_images_module
  use harmonic_properties_module
  use setup_harmonic_module
  implicit none
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
function calculate_harmonic_observables_keywords() result(keywords)
  implicit none
  
  type(KeywordData), allocatable :: keywords(:)
  
  keywords = [                                                                &
  & KeywordData( 'min_temperature',                                           &
  &              'min_temperature is the minimum temperature at which &
  &thermodynamic quantities are calculated. min_temperature should be given &
  &in Kelvin.'),                                                              &
  & KeywordData( 'max_temperature',                                           &
  &              'max_temperature is the maximum temperature at which &
  &thermodynamic quantities are calculated. min_temperature should be given &
  &in Kelvin.'),                                                              &
  & KeywordData( 'no_temperature_steps',                                      &
  &              'no_temperature_steps is the number of temperatures at which &
  &thermodynamic quantities are calculated.',                                 &
  &              default_value='0'),                                          &
  & KeywordData( 'min_frequency',                                             &
  &              'min_frequency is the frequency below which modes will be &
  &ignored when calculating thermodynamic quantities. min_frequency should be &
  &given in Hartree.',                                                        &
  &              default_value='1e-8'),                                       &
  & KeywordData( 'path',                                                      &
  &              'path is the path through fractional reciprocal space which &
  &will be mapped by the phonon dispersion curve. The path should be &
  &specified as labels and q-points, separated by commas. The Gamma-point &
  &should be labelled G.',                                                    &
  &              default_value='G 0.0 0.0 0.0, R 0.5 0.5 0.5, M 0.0 0.5 0.5, &
  &G 0.0 0.0 0.0, X 0.0 0.0 0.5'),                                            &
  & KeywordData( 'no_dos_samples',                                            &
  &              'no_dos_samples is the number of points in reciprocal space &
  &at which the normal modes are calculated when calculating the vibrational &
  &density of states.',                                                       &
  &              default_value='100000')]
end function

function calculate_harmonic_observables_mode() result(output)
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'calculate_harmonic_observables'
  output%description = 'Calculates several observables under the harmonic &
     &approximation: the phonon density of states and dispersion curve, &
     &and the energy, free energy and entropy per unit cell. Should be run &
     &after calculate_harmonic.'
  output%keywords = calculate_harmonic_observables_keywords()
  output%main_subroutine => calculate_harmonic_observables
end function

! ----------------------------------------------------------------------
! The main subroutine.
! ----------------------------------------------------------------------
subroutine calculate_harmonic_observables(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Inputs.
  type(String)                  :: wd
  type(RandomReal)              :: random_generator
  real(dp)                      :: min_temperature
  real(dp)                      :: max_temperature
  integer                       :: no_temperature_steps
  real(dp)                      :: min_frequency
  real(dp),         allocatable :: thermal_energies(:)
  type(String),     allocatable :: path_string(:)
  type(String),     allocatable :: path_labels(:)
  type(RealVector), allocatable :: path_qpoints(:)
  integer                       :: no_dos_samples
  
  ! Previous inputs.
  type(Dictionary) :: setup_harmonic_arguments
  real(dp)         :: symmetry_precision
  
  ! Previously calculated data.
  type(StructureData)                :: structure
  type(StructureData)                :: large_supercell
  type(QpointData),      allocatable :: qpoints(:)
  type(DynamicalMatrix), allocatable :: dynamical_matrices(:)
  
  ! Working variables.
  type(ForceConstants)          :: force_constants
  type(MinImages),  allocatable :: min_images(:,:)
  
  ! Dynamical matrix for checking.
  type(DynamicalMatrix) :: dyn_mat
  
  ! Files.
  type(IFile)                    :: qpoints_file
  type(IFile)                    :: dynamical_matrix_file
  type(OFile)                    :: logfile
  type(StringArray), allocatable :: file_sections(:)
  
  ! Temporary variables.
  type(String), allocatable :: path_point(:)
  integer                   :: i,ialloc
  
  ! --------------------------------------------------
  ! Read in arguments from user.
  ! --------------------------------------------------
  wd = arguments%value('working_directory')
  if (arguments%is_set('random_seed')) then
    random_generator = RandomReal(int(arguments%value('random_seed')))
  else
    random_generator = RandomReal()
  endif
  min_temperature = dble(arguments%value('min_temperature'))
  max_temperature = dble(arguments%value('max_temperature'))
  no_temperature_steps = int(arguments%value('no_temperature_steps'))
  min_frequency = dble(arguments%value('min_frequency'))
  path_string = split_line(arguments%value('path'), ',')
  no_dos_samples = int(arguments%value('no_dos_samples'))
  
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
  allocate( thermal_energies(no_temperature_steps), &
          & stat=ialloc); call err(ialloc)
  do i=1,no_temperature_steps
    thermal_energies(i) = KB_IN_AU                                   &
                      & * ( min_temperature*(no_temperature_steps-i) &
                      &   + max_temperature*(i-1) )                  &
                      & / (no_temperature_steps-1)
  enddo
  
  ! Generate path for dispersion calculation.
  allocate( path_qpoints(size(path_string)), &
          & path_labels(size(path_string)),  &
          & stat=ialloc); call err(ialloc)
  do i=1,size(path_string)
    path_point = split_line(path_string(i))
    path_labels(i) = path_point(1)
    path_qpoints(i) = dble(path_point(2:4))
  enddo
  
  ! --------------------------------------------------
  ! Read in previous arguments.
  ! --------------------------------------------------
  setup_harmonic_arguments = Dictionary(setup_harmonic_keywords())
  call setup_harmonic_arguments%read_file(wd//'/setup_harmonic.used_settings')
  symmetry_precision = &
     & dble(setup_harmonic_arguments%value('symmetry_precision'))
  
  ! --------------------------------------------------
  ! Read in previously calculated data.
  ! --------------------------------------------------
  structure = read_structure_file(wd//'/structure.dat', symmetry_precision)
  
  large_supercell = read_structure_file( wd//'/large_supercell.dat', &
                                       & symmetry_precision,         &
                                       & calculate_symmetry=.false.)
  
  qpoints_file = IFile(wd//'/qpoints.dat')
  file_sections = split_into_sections(qpoints_file%lines())
  allocate(qpoints(size(file_sections)), stat=ialloc); call err(ialloc)
  do i=1,size(qpoints)
    qpoints(i) = file_sections(i)
  enddo
  
  allocate( dynamical_matrices(size(qpoints)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(qpoints)
    dynamical_matrix_file = IFile(                        &
       & wd//'/qpoint_'//left_pad(i,str(size(qpoints)))// &
       & '/dynamical_matrix.dat')
    dynamical_matrices(i) = dynamical_matrix_file%lines()
  enddo
  
  ! --------------------------------------------------
  ! Run calculations.
  ! --------------------------------------------------
  
  logfile = OFile(wd//'/dos_and_dispersion_log.dat')
  
  ! Construct the matrix of force constants from dynamical matrices.
  force_constants = reconstruct_force_constants( large_supercell,    &
                                               & qpoints,            &
                                               & dynamical_matrices, &
                                               & logfile)
  
  ! Calculate minimum image distances.
  min_images = calculate_min_images(large_supercell)
  
  ! Reconstruct calculated dynamical matrices, to check for consistency.
  do i=1,size(qpoints)
    dyn_mat = DynamicalMatrix( dblevec(qpoints(i)%qpoint), &
                             & large_supercell,            &
                             & force_constants,            &
                             & min_images)
    call logfile%print_line('Comparing dynamical matrices before and after &
       &reconstruction of force constants.')
    call compare_dynamical_matrices(dynamical_matrices(i),dyn_mat,logfile)
  enddo
  
  ! Generate harmonic phonon dispersion curve by interpolating between
  !    calculated q-points using Fourier interpolation.
  call generate_dispersion( large_supercell,                    &
                          & min_images,                         &
                          & force_constants,                    &
                          & path_labels,                        &
                          & path_qpoints,                       &
                          & wd//'/phonon_dispersion_curve.dat', &
                          & wd//'/high_symmetry_points.dat',    &
                          & logfile)
  
  ! Generate harmonic phonon density of states, interpolating as above.
  call generate_dos( large_supercell,                     &
                   & min_images,                          &
                   & force_constants,                     &
                   & thermal_energies,                    &
                   & min_frequency,                       &
                   & no_dos_samples,                      &
                   & wd//'/thermodynamic_variables.dat',  &
                   & wd//'/phonon_density_of_states.dat', &
                   & logfile,                             &
                   & random_generator)
end subroutine
end module
