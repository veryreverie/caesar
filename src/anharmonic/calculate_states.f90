! ======================================================================
! Calculates anharmonic states, using the potential calculated by
!   calculate_potential.
! ======================================================================
module calculate_states_module
  use common_module
  
  use states_module
  use anharmonic_common_module
  use potentials_module
  
  use generate_basis_module
  use vscf_module
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
     & KeywordData( 'energy_convergence',                                     &
     &              'energy_convergence is the precision to which energies &
     &will be converged when constructing the VSCF ground state.' ),          &
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
     &              default_value='0.1' ),                                    &
     & KeywordData( 'no_basis_states',                                        &
     &              'no_basis states is the number of states along each mode &
     &in the basis.')                                                         ]
  output%main_subroutine => calculate_states_subroutine
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine calculate_states_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Working directory,
  type(String) :: wd
  
  ! Input arguments.
  real(dp) :: frequency_convergence
  real(dp) :: energy_convergence
  integer  :: max_pulay_iterations
  integer  :: pre_pulay_iterations
  real(dp) :: pre_pulay_damping
  integer  :: no_basis_states
  
  ! Anharmonic data container.
  type(AnharmonicData) :: anharmonic_data
  
  ! Anharmonic potential.
  type(PotentialPointer) :: potential
  
  ! Basis states.
  type(SubspaceBasis), allocatable :: basis(:)
  
  ! VSCF ground states.
  type(SubspacePotentialAndState), allocatable :: potentials_and_states(:)
  type(VscfGroundState),           allocatable :: ground_states(:)
  type(PotentialPointer),          allocatable :: subspace_potentials(:)
  
  ! Files and directories.
  type(IFile) :: anharmonic_data_file
  type(IFile) :: potential_file
  type(OFile) :: basis_file
  type(OFile) :: subspace_potentials_file
  type(OFile) :: ground_state_file
  
  ! Read in arguments.
  wd = arguments%value('working_directory')
  frequency_convergence = dble(arguments%value('frequency_convergence'))
  energy_convergence = dble(arguments%value('energy_convergence'))
  max_pulay_iterations = int(arguments%value('max_pulay_iterations'))
  pre_pulay_iterations = int(arguments%value('pre_pulay_iterations'))
  pre_pulay_damping = dble(arguments%value('pre_pulay_damping'))
  no_basis_states = int(arguments%value('no_basis_states'))
  
  ! Read in anharmonic data.
  anharmonic_data_file = IFile(wd//'/anharmonic_data.dat')
  anharmonic_data = AnharmonicData(anharmonic_data_file%lines())
  
  ! Read in anharmonic potential.
  potential_file = IFile(wd//'/potential.dat')
  potential = PotentialPointer(potential_file%lines())
  
  ! --------------------------------------------------
  ! Generate basis states by generating effective harmonic potential.
  ! --------------------------------------------------
  basis = generate_basis( potential,             &
                        & anharmonic_data,       &
                        & frequency_convergence, &
                        & max_pulay_iterations,  &
                        & pre_pulay_iterations,  &
                        & pre_pulay_damping,     &
                        & no_basis_states        )
  basis_file = OFile(wd//'/basis.dat')
  call basis_file%print_lines(basis, separating_line='')
  
  ! --------------------------------------------------
  ! Run VSCF to generate single-subspace potentials and ground states.
  ! --------------------------------------------------
  potentials_and_states = run_vscf( potential,            &
                                  & basis,                &
                                  & energy_convergence,   &
                                  & max_pulay_iterations, &
                                  & pre_pulay_iterations, &
                                  & pre_pulay_damping,    &
                                  & anharmonic_data       )
  
  subspace_potentials = potentials_and_states%potential
  subspace_potentials_file = OFile(wd//'/subspace_potentials.dat')
  call subspace_potentials_file%print_lines( subspace_potentials,           &
                                           & separating_line=repeat('=',70) )
  
  ground_states = potentials_and_states%state
  ground_state_file = OFile(wd//'/ground_state.dat')
  call ground_state_file%print_lines(ground_states, separating_line='')
end subroutine
end module
