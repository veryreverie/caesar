! ======================================================================
! Calculates anharmonic states using the anharmonic potential and VSCF.
! ======================================================================
module calculate_states_module
  use common_module
  
  use anharmonic_common_module
  use polynomial_module
  
  use setup_harmonic_module
  
  use setup_anharmonic_module
  implicit none
  
  private
  
  public :: calculate_states
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
function calculate_states() result(output)
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'calculate_states'
  output%description = 'Runs VSCF on the anharmonic potential to calculate &
     &states. Should be run after calculate_potential.'
  output%keywords = [                                                         &
  & KeywordData( 'no_single_mode_samples',                                    &
  &              'no_single_mode_samples is the number of points (either side &
  &of zero) along each mode at which the anharmonic potential will be sampled &
  &when determining the effective frequency with which the harmonic basis &
  &along that mode will be constructed.',                                     &
  &              default_value='100') ]
  output%main_subroutine => calculate_states_subroutine
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine calculate_states_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Working directory.
  type(String) :: wd
  
  ! Input arguments.
  integer :: no_single_mode_samples
  
  ! Arguments to setup_anharmonic.
  type(Dictionary) :: setup_anharmonic_arguments
  type(String)     :: harmonic_path
  type(String)     :: potential_representation
  real(dp)         :: maximum_displacement
  
  ! Arguments to setup_harmonic.
  type(Dictionary) :: setup_harmonic_arguments
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
  
  ! Files and directories.
  type(IFile) :: qpoints_file
  type(IFile) :: complex_modes_file
  type(IFile) :: real_modes_file
  type(IFile) :: potential_file
  
  ! --------------------------------------------------
  ! Read in inputs and previously calculated data.
  ! --------------------------------------------------
  
  wd = arguments%value('working_directory')
  
  no_single_mode_samples = int(arguments%value('no_single_mode_samples'))
  
  ! Read in setup_anharmonic arguments.
  setup_anharmonic_arguments = Dictionary(setup_anharmonic())
  call setup_anharmonic_arguments%read_file( &
     & wd//'/setup_anharmonic.used_settings')
  harmonic_path = setup_anharmonic_arguments%value('harmonic_path')
  potential_representation = &
     & setup_anharmonic_arguments%value('potential_representation')
  maximum_displacement = &
     & dble(setup_anharmonic_arguments%value('maximum_displacement'))
  
  ! Read in setup_harmonic arguments.
  setup_harmonic_arguments = Dictionary(setup_harmonic())
  call setup_harmonic_arguments%read_file( &
     & harmonic_path//'/setup_harmonic.used_settings')
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
  ! Re-calculate maximum_weighted_displacement.
  ! --------------------------------------------------
  maximum_weighted_displacement = maximum_displacement &
                              & * sqrt(minval(structure%atoms%mass()))
  
  ! --------------------------------------------------
  ! Calculate effective harmonic potential, from which initial harmonic
  !    states are constructed.
  ! --------------------------------------------------
  
end subroutine
end module
