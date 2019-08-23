! ======================================================================
! Run VSCF on a potential to generate VSCF ground states.
! ======================================================================
! Each step in the routine:
!    1) starts from a set of single-subspace potentials, P.
!    2) uses P to generate a set of single-subspace states, S(P),
!         such that the states diagonalise the potentials.
!    3) uses S(P) and P to generate the free energy, F(P,S(P)).
!    4) checks for convergence of F(P,S(P)), and returns if converged.
!    5) uses S(P) to generate a new set of single-subspace potentials, P(S(P)).
!    6) uses a Pulay scheme to generate the next P, which attempts to minimise
!          F(P,S(P)) subject to P(S(P))=P.
!
! The Pulay scheme sees P as a vector of real coefficients. It is down to the
!    implementation of the potential to define how P is generated.
module vscf_module
  use common_module
  
  use states_module
  
  use anharmonic_common_module
  use potentials_module
  
  use generate_subspace_potentials_module
  implicit none
  
  private
  
  public :: VscfOutput
  
  public :: run_vscf
  
  type, extends(NoDefaultConstructor) :: VscfOutput
    type(PotentialPointer)   :: potential
    type(BasisStatesPointer) :: states
  end type
  
  interface VscfOutput
    module procedure new_VscfOutput
  end interface
contains

! Constructor.
impure elemental function new_VscfOutput(potential,states) result(this)
  implicit none
  
  class(PotentialData), intent(in) :: potential
  Class(BasisStates),   intent(in) :: states
  type(VscfOutput)                 :: this
  
  this%potential = PotentialPointer(potential)
  this%states = BasisStatesPointer(states)
end function

! Vscf routine.
function run_vscf(potential,subspaces,subspace_bases,thermal_energy, &
   & energy_convergence,frequencies,no_converged_calculations,       &
   & max_pulay_iterations,pre_pulay_iterations,pre_pulay_damping,    &
   & anharmonic_data,random_generator,starting_configuration) result(output)
  implicit none
  
  class(PotentialData),     intent(in)           :: potential
  type(DegenerateSubspace), intent(in)           :: subspaces(:)
  class(SubspaceBasis),     intent(in)           :: subspace_bases(:)
  real(dp),                 intent(in)           :: thermal_energy
  real(dp),                 intent(in)           :: energy_convergence
  real(dp),                 intent(in)           :: frequencies(:)
  integer,                  intent(in)           :: no_converged_calculations
  integer,                  intent(in)           :: max_pulay_iterations
  integer,                  intent(in)           :: pre_pulay_iterations
  real(dp),                 intent(in)           :: pre_pulay_damping
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(RandomReal),         intent(in)           :: random_generator
  type(VscfOutput),         intent(in), optional :: starting_configuration(:)
  type(VscfOutput), allocatable                  :: output(:)
  
  ! Single-subspace potentials and states.
  type(PotentialPointer),   allocatable :: subspace_potentials(:)
  type(BasisStatesPointer), allocatable :: subspace_states(:)
  
  ! Potential coefficients.
  integer,  allocatable :: first_coefficient(:)
  integer,  allocatable :: last_coefficient(:)
  real(dp), allocatable :: coefficients(:)
  real(dp)              :: free_energy
  
  ! The Pulay scheme solver.
  type(PulaySolver) :: solver
  
  ! Thermodynamic variables.
  type(ThermodynamicData), allocatable :: thermodynamic_data(:)
  real(dp),                allocatable :: free_energies(:)
  
  ! Temporary variables.
  integer :: i,j,ialloc
  
  ! Generate initial single-subspace states,
  !    and use these states to generate initial single-subpsace potentials.
  call print_line('Generating initial states and potentials.')
  subspace_states = BasisStatesPointer(subspace_bases%initial_states( &
                                                    & subspaces,      &
                                                    & anharmonic_data ))
  subspace_potentials = generate_subspace_potentials( potential,       &
                                                    & subspaces,       &
                                                    & subspace_bases,  &
                                                    & subspace_states, &
                                                    & thermal_energy,  &
                                                    & anharmonic_data  )
  
  ! Construct a single vector of coefficients from the subspace potentials.
  allocate( coefficients(0),                    &
          & first_coefficient(size(subspaces)), &
          & last_coefficient(size(subspaces)),  &
          & stat=ialloc); call err(ialloc)
  do i=1,size(subspaces)
    first_coefficient(i) = size(coefficients)+1
    coefficients = [ coefficients,                                          &
                   & subspace_potentials(i)%coefficients( frequencies(i),   &
                   &                                      anharmonic_data ) ]
    last_coefficient(i) = size(coefficients)
  enddo
  
  ! Initialise Pulay solver.
  solver = PulaySolver( pre_pulay_iterations,      &
                      & pre_pulay_damping,         &
                      & max_pulay_iterations,      &
                      & energy_convergence,        &
                      & no_converged_calculations, &
                      & random_generator,          &
                      & coefficients               )
  
  ! Run Pulay scheme.
  i = 1
  free_energies = [real(dp)::]
  do
    ! Generate single-subspace potentials from Pulay scheme.
    coefficients = solver%get_x()
    do j=1,size(subspaces)
      call subspace_potentials(j)%set_coefficients(                &
         & coefficients(first_coefficient(j):last_coefficient(j)), &
         & frequencies(j),                                         &
         & anharmonic_data                                         )
    enddo
    
    ! Use single-subspace potentials to calculate new single-subspace states.
    call print_line('Generating single-subspace ground states.')
    subspace_states = BasisStatesPointer(                              &
       & subspace_bases%calculate_states( subspaces,                   &
       &                                  subspace_potentials,         &
       &                                  thermal_energy,              &
       &                                  energy_convergence,          &
       &                                  no_converged_calculations,   &
       &                                  max_pulay_iterations,        &
       &                                  pre_pulay_iterations,        &
       &                                  pre_pulay_damping,           &
       &                                  anharmonic_data            ) )
    
    ! Calculate the free energy from the potentials and states.
    thermodynamic_data = subspace_bases%thermodynamic_data( &
                          & thermal_energy,                 &
                          & subspace_states,                &
                          & subspaces,                      &
                          & subspace_potentials,            &
                          & anharmonic_data=anharmonic_data )
    
    ! Use single-subspace states to calculate new single-subspace potentials.
    call print_line('Generating single-subspace potentials.')
    subspace_potentials = generate_subspace_potentials( potential,       &
                                                      & subspaces,       &
                                                      & subspace_bases,  &
                                                      & subspace_states, &
                                                      & thermal_energy,  &
                                                      & anharmonic_data  )
    
    ! Update the Pulay scheme.
    do j=1,size(subspaces)
      coefficients(first_coefficient(j):last_coefficient(j)) = &
         & subspace_potentials(j)%coefficients(frequencies(j), anharmonic_data)
    enddo
    
    free_energy = sum(thermodynamic_data%free_energy) &
              & / anharmonic_data%anharmonic_supercell%sc_size
    
    call solver%set_f(coefficients, free_energy)
    
    ! Print progress.
    if (i>1) then
      call print_line('Self-consistency step '//i//'.')
      call print_line( 'Self-consistency error : '     // &
                     & solver%self_consistency_error() // &
                     & ' (Ha)'                            )
      call print_line( 'Change in free energy  : '  // &
                     & solver%free_energy_change()  // &
                     & ' (Ha)'                         )
    endif
    
    ! Check for convergence.
    if (solver%converged()) then
      output = VscfOutput(subspace_potentials, subspace_states)
      call print_line('Convergence reached.')
      return
    endif
    
    ! Increment the loop counter.
    i = i+1
  enddo
end function
end module
