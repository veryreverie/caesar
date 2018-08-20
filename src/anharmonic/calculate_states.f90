! ======================================================================
! Calculates anharmonic states, using the potential calculated by
!   calculate_potential.
! ======================================================================
module calculate_states_module
  use common_module
  
  use setup_harmonic_module
  use calculate_normal_modes_module
  
  use anharmonic_common_module
  use polynomial_module
  
  use setup_anharmonic_module
  implicit none
  
  private
  
  public :: calculate_states
contains

! ----------------------------------------------------------------------
! Generates keywords and helptext.
! ----------------------------------------------------------------------
function calculate_states() result(output)
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'calculate_states'
  output%description = 'Calculates anharmonic states using VSCF. Should be &
     &run after calculate_potential.'
  output%keywords = [                                                         &
     & KeywordData( 'frequency_convergence',                                  &
     &              'frequency_convergence is the precision to which &
     &frequencies will be converged when constructing the harmonic ground &
     &state.' ),                                                              &
     & KeywordData( 'max_pulay_iterations',                                   &
     &              'max_pulay_iterations is the maximum number of &
     &self-consistency iterations which will be passed into the Pulay &
     &scheme.',                                                               &
     &              default_value='20' ),                                     &
     & KeywordData( 'pre_pulay_iterations',                                   &
     &              'pre_pulay_iterations is the number of damped iterations &
     &which will be performed before the Pulay scheme is called.',            &
     &              default_value='2' ),                                      &
     & KeywordData( 'pre_pulay_damping',                                      &
     &              'pre_pulay_damping is the damping factor of the pre-Pulay &
     &iterations.',                                                           &
     &              default_value='0.1' )                                     ]
  output%main_subroutine => calculate_states_subroutine
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine calculate_states_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Input arguments.
  real(dp) :: frequency_convergence
  integer  :: max_pulay_iterations
  integer  :: pre_pulay_iterations
  real(dp) :: pre_pulay_damping
  
  ! Working directory,
  type(String) :: wd
  
  ! Arguments to setup_anharmonic.
  type(Dictionary) :: setup_anharmonic_arguments
  type(String)     :: harmonic_path
  type(String)     :: potential_representation
  integer          :: potential_expansion_order
  logical          :: vscf_basis_functions_only
  real(dp)         :: maximum_displacement
  real(dp)         :: frequency_of_max_displacement
  
  ! Arguments to setup_harmonic.
  type(Dictionary) :: setup_harmonic_arguments
  
  ! Arguments to calculate_normal_modes.
  type(Dictionary) :: calculate_normal_modes_arguments
  
  ! Primitive structure.
  type(StructureData) :: structure
  
  ! Large anharmonic supercell and its q-points.
  type(StructureData)           :: anharmonic_supercell
  type(QpointData), allocatable :: qpoints(:)
  
  ! Normal modes.
  type(ComplexMode), allocatable :: complex_modes(:)
  type(RealMode),    allocatable :: real_modes(:)
  
  ! maximum_displacement in mass-weighted co-ordinates.
  real(dp) :: maximum_weighted_displacement
  
  ! Degeneracy and coupling data.
  type(DegenerateSubspace), allocatable :: degenerate_subspaces(:)
  type(DegenerateSymmetry), allocatable :: degenerate_symmetries(:)
  type(SubspaceCoupling),   allocatable :: subspace_coupling(:)
  
  ! Anharmonic data container.
  type(AnharmonicData) :: anharmonic_data
  
  ! Anharmonic potential.
  type(PotentialPointer) :: potential
  
  ! Files and directories.
  type(IFile) :: structure_file
  type(IFile) :: anharmonic_supercell_file
  type(IFile) :: qpoints_file
  type(IFile) :: complex_modes_file
  type(IFile) :: real_modes_file
  type(IFile) :: subspaces_file
  type(IFile) :: subspace_coupling_file
  type(IFile) :: symmetry_file
  type(IFile) :: potential_file
  
  ! Temporary variables.
  integer :: i,ialloc
  
  ! --------------------------------------------------
  ! Read in inputs and previously calculated data.
  ! --------------------------------------------------
  
  wd = arguments%value('working_directory')
  frequency_convergence = dble(arguments%value('frequency_convergence'))
  max_pulay_iterations = int(arguments%value('max_pulay_iterations'))
  pre_pulay_iterations = int(arguments%value('pre_pulay_iterations'))
  pre_pulay_damping = dble(arguments%value('pre_pulay_damping'))
  
  ! Read in setup_anharmonic arguments.
  setup_anharmonic_arguments = Dictionary(setup_anharmonic())
  call setup_anharmonic_arguments%read_file( &
     & wd//'/setup_anharmonic.used_settings')
  harmonic_path = setup_anharmonic_arguments%value('harmonic_path')
  potential_representation = &
     & setup_anharmonic_arguments%value('potential_representation')
  potential_expansion_order = &
     & int(setup_anharmonic_arguments%value('potential_expansion_order'))
  vscf_basis_functions_only = &
     & lgcl(setup_anharmonic_arguments%value('vscf_basis_functions_only'))
  maximum_displacement = &
     & dble(setup_anharmonic_arguments%value('maximum_displacement'))
  frequency_of_max_displacement = &
     & dble(setup_anharmonic_arguments%value('frequency_of_max_displacement'))
  
  ! Read in setup_harmonic arguments.
  setup_harmonic_arguments = Dictionary(setup_harmonic())
  call setup_harmonic_arguments%read_file( &
     & harmonic_path//'/setup_harmonic.used_settings')
  
  ! Read in calculate_normal_modes arguments.
  calculate_normal_modes_arguments = Dictionary(calculate_normal_modes())
  call calculate_normal_modes_arguments%read_file( &
     & harmonic_path//'/calculate_normal_modes.used_settings')
  
  ! Read in structure.
  structure_file = IFile(harmonic_path//'/structure.dat')
  structure = StructureData(structure_file%lines())
  
  ! Read in large anharmonic supercell and its q-points.
  anharmonic_supercell_file = IFile(wd//'/anharmonic_supercell.dat')
  anharmonic_supercell = StructureData(anharmonic_supercell_file%lines())
  
  qpoints_file = IFile(wd//'/qpoints.dat')
  qpoints = QpointData(qpoints_file%sections())
  
  ! Read in normal modes.
  complex_modes_file = IFile(wd//'/complex_modes.dat')
  complex_modes = ComplexMode(complex_modes_file%sections())
  
  real_modes_file = IFile(wd//'/real_modes.dat')
  real_modes = RealMode(real_modes_file%sections())
  
  ! Read in degenerate subspaces and subspace coupling.
  subspaces_file = IFile(wd//'/degenerate_subspaces.dat')
  degenerate_subspaces = DegenerateSubspace(subspaces_file%sections())
  
  subspace_coupling_file = IFile(wd//'/subspace_coupling.dat')
  subspace_coupling = SubspaceCoupling(subspace_coupling_file%lines())
  
  ! Read in degenerate symmetries.
  symmetry_file = IFile(wd//'/symmetries.dat')
  degenerate_symmetries = DegenerateSymmetry(symmetry_file%sections())
  
  ! Read in anharmonic potential.
  potential_file = IFile(wd//'/potential.dat')
  if (potential_representation=='polynomial') then
    potential = PolynomialPotential(potential_file%lines())
  else
    call print_line( ERROR//': Unrecognised potential representation : '// &
                   & potential_representation)
    call err()
  endif
  
  ! Re-calculate maximum_weighted_displacement.
  maximum_weighted_displacement = maximum_displacement &
                              & * sqrt(minval(structure%atoms%mass()))
  
  ! --------------------------------------------------
  ! Load anharmonic data into container.
  ! --------------------------------------------------
  anharmonic_data = AnharmonicData( structure,                     &
                                  & anharmonic_supercell,          &
                                  & qpoints,                       &
                                  & complex_modes,                 &
                                  & real_modes,                    &
                                  & degenerate_subspaces,          &
                                  & degenerate_symmetries,         &
                                  & subspace_coupling,             &
                                  & vscf_basis_functions_only,     &
                                  & maximum_weighted_displacement, &
                                  & frequency_of_max_displacement )
  
  ! --------------------------------------------------
  ! Calculate frequencies from which to generate a harmonic basis.
  ! --------------------------------------------------
  call calculate_frequencies( potential,             &
                            & anharmonic_data,       &
                            & frequency_convergence, &
                            & max_pulay_iterations,  &
                            & pre_pulay_iterations,  &
                            & pre_pulay_damping      )
  
end subroutine
end module
