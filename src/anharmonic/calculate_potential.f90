! ======================================================================
! Calculates the anharmonic potential, using the results of run_anharmonic.
! ======================================================================
module calculate_potential_module
  use common_module
  
  use anharmonic_common_module
  use potentials_module
  
  use setup_harmonic_module
  
  use setup_anharmonic_module
  implicit none
  
  private
  
  public :: calculate_potential
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
function calculate_potential() result(output)
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'calculate_potential'
  output%description = 'Uses the results of run_anharmonic to calculate &
     &the anharmonic potential. Should be run after run_anharmonic.'
  output%keywords = [                                                      &
     & KeywordData( 'energy_to_force_ratio',                                  &
     &              'energy_to_force_ratio is the ratio of how penalised &
     &deviations in energy are compared to deviations in forces when the &
     &potential is being fitted. This should be given in units of (Hartree &
     &per primitive cell) divided by (Hartree per bohr). Due to the use of &
     &mass-weighted co-ordinates, in systems containing different elements &
     &forces along modes with higher contributions from heavier elements will &
     &be weighted less than this.'),                                          &
     & KeywordData( 'loto_direction',                                         &
     &              'loto_direction specifies the direction (in reciprocal &
     &co-ordinates from which the gamma point is approached when calculating &
     &LO/TO corrections. See setup_harmonic help for details. loto_direction &
     &may only be specified here if it was not specified in setup_harmonic, &
     &and if every symmetry of the system leaves it invariant, i.e. q.S=q for &
     &all S, where S is the symmetry matrix and q is loto_direction. See &
     &structure.dat for the list of symmetries.',                             &
     &              is_optional = .true.)                                     ]
  output%main_subroutine => calculate_potential_subroutine
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine calculate_potential_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Input arguments.
  real(dp)             :: energy_to_force_ratio
  real(dp)             :: weighted_energy_force_ratio
  type(Fractionvector) :: loto_direction
  logical              :: loto_direction_set
  
  ! Arguments to setup_harmonic.
  type(Dictionary) :: setup_harmonic_arguments
  
  ! Arguments to setup_anharmonic.
  type(Dictionary) :: setup_anharmonic_arguments
  type(String)     :: potential_representation
  integer          :: potential_expansion_order
  
  ! Electronic structure calculation reader.
  type(CalculationReader) :: calculation_reader
  
  ! Anharmonic data container.
  type(AnharmonicData) :: anharmonic_data
  
  ! Anharmonic potential.
  type(PotentialPointer) :: potential
  
  ! Files and directories.
  type(IFile)  :: anharmonic_data_file
  type(String) :: sampling_points_dir
  type(OFile)  :: logfile
  type(OFile)  :: potential_file
  
  ! Read in arguments.
  energy_to_force_ratio = dble(arguments%value('energy_to_force_ratio'))
  
  if (arguments%is_set('loto_direction')) then
    loto_direction = FractionVector(                      &
       & setup_harmonic_arguments%value('loto_direction') )
    loto_direction_set = .true.
  endif
  
  ! Read in setup_harmonic arguments.
  setup_harmonic_arguments = Dictionary(setup_harmonic())
  call setup_harmonic_arguments%read_file('setup_harmonic.used_settings')
  
  ! Read in setup_anharmonic arguments.
  setup_anharmonic_arguments = Dictionary(setup_anharmonic())
  call setup_anharmonic_arguments%read_file( &
          & 'setup_anharmonic.used_settings' )
  potential_representation = &
     & setup_anharmonic_arguments%value('potential_representation')
  potential_expansion_order = &
     & int(setup_anharmonic_arguments%value('potential_expansion_order'))
  
  ! Read in anharmonic data.
  anharmonic_data_file = IFile('anharmonic_data.dat')
  anharmonic_data = AnharmonicData(anharmonic_data_file%lines())
  
  ! Initialise LO/TO splitting if necessary.
  if (setup_harmonic_arguments%is_set('loto_direction')) then
    loto_direction = FractionVector(                      &
       & setup_harmonic_arguments%value('loto_direction') )
    loto_direction_set = .true.
  endif
  
  if (arguments%is_set('loto_direction')) then
    if (loto_direction_set) then
      call print_line(ERROR//': loto_direction may not be specified here, &
         &since it was already specified in setup_harmonic.')
      stop
    endif
    loto_direction = FractionVector(                      &
       & setup_harmonic_arguments%value('loto_direction') )
    loto_direction_set = .true.
    if (any(loto_breaks_symmetry(                     &
       & anharmonic_data%structure%symmetries%tensor, &
       & loto_direction                               ))) then
      call print_line(ERROR//': loto_direction has been specified in a &
         &direction which breaks symmetry. To specify this direction, please &
         &set loto_direction when running setup_harmonic.')
      stop
    endif
  endif
  
  ! Initialise calculation reader.
  if (loto_direction_set) then
    calculation_reader = CalculationReader(loto_direction)
  else
    calculation_reader = CalculationReader()
  endif
  
  ! Calculate weighted energy to force ratio.
  weighted_energy_force_ratio = &
     &   energy_to_force_ratio  &
     & * sqrt(maxval(anharmonic_data%structure%atoms%mass()))
  
  ! Initialise potential to the chosen representation
  if (potential_representation=='polynomial') then
    potential = PotentialPointer(                       &
       & PolynomialPotential(potential_expansion_order) )
  else
    call print_line( ERROR//': Unrecognised potential representation: '// &
                   & potential_representation)
    call err()
  endif
  
  ! Open logfile.
  logfile = OFile('setup_anharmonic_logfile.dat')
  
  ! Generate the potential itself.
  sampling_points_dir = 'sampling_points'
  call potential%generate_potential( anharmonic_data,             &
                                   & weighted_energy_force_ratio, &
                                   & sampling_points_dir,         &
                                   & calculation_reader,          &
                                   & logfile                      )
  
  ! Write the potential to file.
  potential_file = OFile('potential.dat')
  call potential_file%print_lines(potential)
end subroutine
end module
