! ======================================================================
! Calculates the anharmonic potential, using the results of run_anharmonic.
! ======================================================================
module calculate_potential_module
  use common_module
  
  use anharmonic_common_module
  use potentials_module
  
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
  &be weighted less than this.')]
  output%main_subroutine => calculate_potential_subroutine
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine calculate_potential_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Working directory.
  type(String) :: wd
  
  ! Input arguments.
  real(dp) :: energy_to_force_ratio
  real(dp) :: weighted_energy_force_ratio
  
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
  wd = arguments%value('working_directory')
  energy_to_force_ratio = dble(arguments%value('energy_to_force_ratio'))
  
  ! Read in setup_anharmonic arguments.
  setup_anharmonic_arguments = Dictionary(setup_anharmonic())
  call setup_anharmonic_arguments%read_file( &
     & wd//'/setup_anharmonic.used_settings')
  potential_representation = &
     & setup_anharmonic_arguments%value('potential_representation')
  potential_expansion_order = &
     & int(setup_anharmonic_arguments%value('potential_expansion_order'))
  
  ! Read in anharmonic data.
  anharmonic_data_file = IFile(wd//'/anharmonic_data.dat')
  anharmonic_data = AnharmonicData(anharmonic_data_file%lines())
  
  ! Initialise calculation reader.
  calculation_reader = CalculationReader()
  
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
  logfile = OFile(wd//'/setup_anharmonic_logfile.dat')
  
  ! Generate the potential itself.
  sampling_points_dir = wd//'/sampling_points'
  call potential%generate_potential( anharmonic_data,             &
                                   & weighted_energy_force_ratio, &
                                   & sampling_points_dir,         &
                                   & calculation_reader,          &
                                   & logfile                      )
  
  ! Write the potential to file.
  potential_file = OFile(wd//'/potential.dat')
  call potential_file%print_lines(potential)
end subroutine
end module
