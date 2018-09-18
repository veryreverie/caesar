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
     &Should be run after calculate_states.'
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
     & KeywordData( 'no_basis_states',                                        &
     &              'no_basis states is the number of states along each mode &
     &in the basis.'),                                                        &
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
  type(String)                  :: wd
  type(RandomReal)              :: random_generator
  real(dp)                      :: min_temperature
  real(dp)                      :: max_temperature
  integer                       :: no_temperature_steps
  integer                       :: no_basis_states
  real(dp)                      :: frequency_convergence
  integer                       :: no_converged_calculations
  real(dp)                      :: min_frequency
  real(dp),         allocatable :: thermal_energies(:)
  type(String),     allocatable :: path_string(:)
  type(String),     allocatable :: path_labels(:)
  type(RealVector), allocatable :: path_qpoints(:)
  integer                       :: no_dos_samples
  
  ! Anharmonic data.
  type(AnharmonicData) :: anharmonic_data
  
  ! Single-subspace potentials.
  type(PotentialPointer), allocatable :: subspace_potentials(:)
  
  ! Single-subspace basis. Only used to initalise frequency-finding.
  type(SubspaceBasis), allocatable :: subspace_bases(:)
  
  ! Subspace.
  type(DegenerateSubspace) :: subspace
  
  ! Finite-temperature effective harmonic frequencies.
  real(dp), allocatable :: initial_frequencies(:)
  real(dp), allocatable :: effective_frequencies(:)
  
  ! Files and directories.
  type(IFile) :: anharmonic_data_file
  type(IFile) :: subspace_potentials_file
  type(IFile) :: basis_file
  
  integer :: i,j,ialloc
  
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
  frequency_convergence = dble(arguments%value('frequency_convergence'))
  no_converged_calculations = int(arguments%value('no_converged_calculations'))
  min_frequency = dble(arguments%value('min_frequency'))
  no_basis_states = int(arguments%value('no_basis_states'))
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
  
  ! Read in anharmonic data.
  anharmonic_data_file = IFile(wd//'/anharmonic_data.dat')
  anharmonic_data = AnharmonicData(anharmonic_data_file%lines())
  
  ! Read in subspace potentials.
  subspace_potentials_file = IFile(wd//'/subspace_potentials.dat')
  subspace_potentials = PotentialPointer(                                &
     & subspace_potentials_file%sections(separating_line=repeat('=',70)) )
  
  ! Read in subspace bases.
  basis_file = IFile(wd//'/basis.dat')
  subspace_bases = SubspaceBasis(basis_file%sections())
  
  ! --------------------------------------------------
  ! Generate observables for each temperature in turn.
  ! --------------------------------------------------
  ! Initialise frequencies to the frequencies which minimise the energy
  !    of the harmonic ground state.
  initial_frequencies = subspace_bases%frequency
  allocate( effective_frequencies(size(subspace_potentials)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(thermal_energies)
    ! Calculate effective frequencies for each subspace.
    call print_line('')
    call print_line('Thermal energy '//i//' of '//size(thermal_energies))
    do j=1,size(subspace_potentials)
      subspace = anharmonic_data%degenerate_subspaces(j)
      effective_frequencies(j) = calculate_effective_frequency( &
                                    & subspace_potentials(j),   &
                                    & subspace,                 &
                                    & anharmonic_data,          &
                                    & thermal_energies(i),      &
                                    & initial_frequencies(j),   &
                                    & no_basis_states,          &
                                    & frequency_convergence,    &
                                    & no_converged_calculations )
    enddo
    
    call print_line('Self-consistent harmonic frequencies calculated.')
    call print_line(effective_frequencies)
    
    ! The starting point for calculating the effective frequencies at the
    !    next temperature is the frequencies at this temperature.
    initial_frequencies = effective_frequencies
  enddo
  
  ! TODO
end subroutine
end module
