! ======================================================================
! Run VSCF on a potential to generate VSCF ground states.
! ======================================================================
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

! ----------------------------------------------------------------------
! Main functions.
! ----------------------------------------------------------------------
function run_vscf(potential,subspaces,subspace_bases,thermal_energy,    &
   & energy_convergence,no_converged_calculations,max_pulay_iterations, &
   & pre_pulay_iterations,pre_pulay_damping,anharmonic_data)            &
   & result(output)
  implicit none
  
  class(PotentialData),     intent(in) :: potential
  type(DegenerateSubspace), intent(in) :: subspaces(:)
  class(SubspaceBasis),     intent(in) :: subspace_bases(:)
  real(dp),                 intent(in) :: thermal_energy
  real(dp),                 intent(in) :: energy_convergence
  integer,                  intent(in) :: no_converged_calculations
  integer,                  intent(in) :: max_pulay_iterations
  integer,                  intent(in) :: pre_pulay_iterations
  real(dp),                 intent(in) :: pre_pulay_damping
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(VscfOutput), allocatable        :: output(:)
  
  type(PulaySolver), allocatable :: solvers(:)
  
  type(PotentialPointer),   allocatable :: input_potentials(:)
  type(BasisStatesPointer), allocatable :: subspace_states(:)
  type(PotentialPointer),   allocatable :: output_potentials(:)
  
  type(PotentialPointer), allocatable :: in_potentials(:)
  type(PotentialPointer), allocatable :: out_potentials(:)
  
  type(ThermodynamicData), allocatable :: thermodynamic_data(:)
  real(dp),                allocatable :: free_energies(:)
  
  integer :: first_pulay_step
  
  integer :: i,j,k,ialloc
  
  ! Generate initial states,
  !    and use these states to generate initial potentials.
  call print_line('Generating initial states and potentials.')
  subspace_states = BasisStatesPointer(subspace_bases%initial_states( &
                                                    & subspaces,      &
                                                    & anharmonic_data ))
  input_potentials = generate_subspace_potentials( potential,       &
                                                 & subspaces,       &
                                                 & subspace_bases,  &
                                                 & subspace_states, &
                                                 & anharmonic_data  )
  ! Initialise Pulay solvers.
  solvers = [(                                            &
     & PulaySolver( pre_pulay_iterations,                 &
     &              pre_pulay_damping,                    &
     &              max_pulay_iterations,                 &
     &              input_potentials(i)%coefficients() ), &
     & i=1,                                               &
     & size(input_potentials)                             )]
  
  ! Run Pulay scheme.
  i = 1
  free_energies = [real(dp)::]
  do
    call print_line('Beginning VSCF self-consistency step '//i//'.')
    do j=1,size(input_potentials)
      call input_potentials(i)%set_coefficients(solvers(i)%get_x())
    enddo
    
    ! Use the single-subspace potentials to calculate the new states.
    call print_line('Generating single-subspace ground states.')
    subspace_states = BasisStatesPointer(                              &
       & subspace_bases%calculate_states( subspaces,                   &
       &                                  input_potentials,            &
       &                                  energy_convergence,          &
       &                                  no_converged_calculations,   &
       &                                  max_pulay_iterations,        &
       &                                  pre_pulay_iterations,        &
       &                                  pre_pulay_damping,           &
       &                                  anharmonic_data            ) )
    
    ! Generate the free energy from the states.
    thermodynamic_data = subspace_bases%thermodynamic_data( &
                          & thermal_energy,                 &
                          & subspace_states,                &
                          & subspaces,                      &
                          & input_potentials,               &
                          & anharmonic_data=anharmonic_data )
    free_energies = [ free_energies,                                 &
                    &   sum(thermodynamic_data%free_energy)          &
                    & / anharmonic_data%anharmonic_supercell%sc_size ]
    call print_line('Free energy: '//free_energies(i)//' (Ha).')
    
    ! Check whether the energies have converged by the normal convergence
    !    condition.
    j = i-no_converged_calculations
    if (j>0) then
      if (all( abs(free_energies(j:i-1)-free_energies(i)) &
           & < energy_convergence                         )) then
        output = VscfOutput(input_potentials, subspace_states)
        exit
      endif
    endif
    
    ! Check whether the last two sets of energies have converged to a much
    !    tighter convergence.
    ! This is needed to avoid numerical problems with the Pulay scheme caused
    !    by over-convergence.
    if (i>1) then
      if (abs(free_energies(i-1)-free_energies(i))<energy_convergence/100) then
        output = VscfOutput(input_potentials, subspace_states)
        exit
      endif
    endif
    
    ! Use the current single-subspace states to calculate the single-subspace
    !    potentials.
    call print_line('Generating single-subspace potentials.')
    output_potentials = generate_subspace_potentials( potential,       &
                                                    & subspaces,       &
                                                    & subspace_bases,  &
                                                    & subspace_states, &
                                                    & anharmonic_data  )
    
    ! If the energies have not converged, generate the next input potentials
    !    using either a damped iterative scheme or a Pulay scheme.
    do j=1,size(output_potentials)
      call solvers(j)%set_f(output_potentials(j)%coefficients())
    enddo
    
    ! Increment the loop counter.
    i = i+1
  enddo
end function
end module
