! ======================================================================
! Run VSCF on a potential to generate VSCF ground states.
! ======================================================================
module vscf_module
  use common_module
  
  use states_module
  
  use anharmonic_common_module
  use potentials_module
  implicit none
  
  private
  
  public :: VscfOutput
  
  public :: run_vscf
  
  type, extends(NoDefaultConstructor) :: VscfStep
    type(PotentialPointer),      allocatable :: input_potentials(:)
    type(SubspaceStatesPointer), allocatable :: states(:)
    type(EnergySpectra),         allocatable :: spectra(:)
    type(PotentialPointer),      allocatable :: output_potentials(:)
  end type
  
  interface VscfStep
    module procedure new_VscfStep
  end interface
  
  type, extends(NoDefaultConstructor) :: VscfOutput
    type(PotentialPointer)      :: potential
    type(SubspaceStatesPointer) :: states
  end type
  
  interface VscfOutput
    module procedure new_VscfOutput
  end interface
contains

! Constructors.
function new_VscfStep(input_potentials,states,spectra,output_potentials) &
   & result(this)
  implicit none
  
  class(PotentialData),  intent(in) :: input_potentials(:)
  class(SubspaceStates), intent(in) :: states(:)
  type(EnergySpectra),   intent(in) :: spectra(:)
  class(PotentialData),  intent(in) :: output_potentials(:)
  type(VscfStep)                    :: this
  
  if (size(input_potentials)/=size(states)) then
    call print_line(CODE_ERROR//': Input potentials and states do not match.')
    call err()
  elseif (size(input_potentials)/=size(spectra)) then
    call print_line(CODE_ERROR//': Input potentials and spectra do not match.')
    call err()
  elseif (size(input_potentials)/=size(output_potentials)) then
    call print_line(CODE_ERROR//': Input potentials and output potentials do &
       &not match.')
    call err()
  endif
  
  this%input_potentials = PotentialPointer(input_potentials)
  this%states = SubspaceStatesPointer(states)
  this%spectra = spectra
  this%output_potentials = PotentialPointer(output_potentials)
end function

impure elemental function new_VscfOutput(potential,states) result(this)
  implicit none
  
  class(PotentialData),  intent(in) :: potential
  Class(SubspaceStates), intent(in) :: states
  type(VscfOutput)                  :: this
  
  this%potential = PotentialPointer(potential)
  this%states = SubspaceStatesPointer(states)
end function

! ----------------------------------------------------------------------
! Main functions.
! ----------------------------------------------------------------------

function run_vscf(potential,subspaces,subspace_bases,energy_convergence,  &
   & no_converged_calculations,max_pulay_iterations,pre_pulay_iterations, &
   & pre_pulay_damping,anharmonic_data) result(output)
  implicit none
  
  class(PotentialData),     intent(in) :: potential
  type(DegenerateSubspace), intent(in) :: subspaces(:)
  class(SubspaceBasis),     intent(in) :: subspace_bases(:)
  real(dp),                 intent(in) :: energy_convergence
  integer,                  intent(in) :: no_converged_calculations
  integer,                  intent(in) :: max_pulay_iterations
  integer,                  intent(in) :: pre_pulay_iterations
  real(dp),                 intent(in) :: pre_pulay_damping
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(VscfOutput), allocatable        :: output(:)
  
  type(PotentialPointer),      allocatable :: input_potentials(:)
  type(SubspaceStatesPointer), allocatable :: subspace_states(:)
  type(EnergySpectra),         allocatable :: subspace_spectra(:)
  type(PotentialPointer),      allocatable :: output_potentials(:)
  
  type(PotentialPointer), allocatable :: in_potentials(:)
  type(PotentialPointer), allocatable :: out_potentials(:)
  
  type(VscfStep), allocatable :: vscf_steps(:)
  
  integer :: first_pulay_step
  
  integer :: i,j,k,ialloc
  
  vscf_steps = [VscfStep::]
  
  ! Generate initial states,
  !    and use these states to generate initial potentials.
  call print_line('Generating initial states and potentials.')
  subspace_states = SubspaceStatesPointer(subspace_bases%initial_states( &
                                                       & subspaces,      &
                                                       & anharmonic_data ))
  input_potentials = generate_subspace_potentials( potential,       &
                                                 & subspaces,       &
                                                 & subspace_bases,  &
                                                 & subspace_states, &
                                                 & anharmonic_data  )
    
  i = 1
  do
    call print_line('Beginning VSCF self-consistency step '//i//'.')
    
    ! Use the single-subspace potentials to calculate the new states.
    call print_line('Generating single-subspace ground states.')
    subspace_states = SubspaceStatesPointer(                &
       & subspace_bases%calculate_states( subspaces,        &
       &                                  input_potentials, &
       &                                  anharmonic_data   ) )
    
    ! Generate the energy spectra from the states.
    subspace_spectra = subspace_states%spectra( subspaces,      &
                                              & subspace_bases, &
                                              & anharmonic_data )
    call print_line( 'Ground-state energy: '//                    &
                   & sum(subspace_spectra%min_energy())//' (Ha).' )
    
    ! Use the current single-subspace states to calculate the single-subspace
    !    potentials.
    call print_line('Generating single-subspace potentials.')
    output_potentials = generate_subspace_potentials( potential,       &
                                                    & subspaces,       &
                                                    & subspace_bases,  &
                                                    & subspace_states, &
                                                    & anharmonic_data  )
    
    ! Store the input and output potentials, states and spectra as the next
    !    VSCF step.
    vscf_steps = [ vscf_steps,                   &
                 & VscfStep( input_potentials,   &
                 &           subspace_states,    &
                 &           subspace_spectra,   &
                 &           output_potentials ) ]
    
    ! Check whether the energies have converged by the normal convergence
    !    condition.
    if (i>no_converged_calculations) then
      if (all(steps_converged(                          &
         & vscf_steps(i-no_converged_calculations:i-1), &
         & vscf_steps(i),                               &
         & energy_convergence                           ))) then
        output = VscfOutput(input_potentials, subspace_states)
        exit
      endif
    endif
    
    ! Check whether the last two sets of energies have converged to a much
    !    tighter convergence.
    ! This is needed to avoid numerical problems with the Pulay scheme caused
    !    by over-convergence.
    if (i>1) then
      if (steps_converged( vscf_steps(i),         &
                         & vscf_steps(i-1),       &
                         & energy_convergence/100 )) then
        output = VscfOutput(input_potentials, subspace_states)
        exit
      endif
    endif
    
    ! If the energies have not converged, generate the next input potentials
    !    using either a damped iterative scheme or a Pulay scheme.
    if (i<=pre_pulay_iterations) then
      input_potentials = PotentialPointer(                       &
         & input_potentials%iterate_damped( output_potentials,   &
         &                                  pre_pulay_damping,   &
         &                                  anharmonic_data    ) )
    else
      first_pulay_step = max(1, i-max_pulay_iterations+1)
      do j=1,size(subspaces)
        in_potentials = [( vscf_steps(k)%input_potentials(j), &
                        &  k=first_pulay_step,                &
                        &  i                                  )]
        out_potentials = [( vscf_steps(k)%output_potentials(j), &
                         &  k=first_pulay_step,                 &
                         &  i                                   )]
        input_potentials(j) = PotentialPointer(                   &
           & input_potentials(j)%iterate_pulay( in_potentials,    &
           &                                    out_potentials,   &
           &                                    anharmonic_data ) )
      enddo
    endif
    
    ! Increment the loop counter.
    i = i+1
  enddo
end function

! Check if two vscf steps have converged w/r/t one another.
impure elemental function steps_converged(this,that,energy_convergence) &
   & result(output)
  implicit none
  
  type(VscfStep), intent(in) :: this
  type(VscfStep), intent(in) :: that
  real(dp),       intent(in) :: energy_convergence
  logical                    :: output
  
  output = all(converged(this%spectra, that%spectra, energy_convergence))
end function

! Appends a column to a spectra array.
! output(i,:) = [old_spectra(i,:), new_spectra(i)].
function append_spectra(old_spectra, new_spectra) result(output)
  implicit none
  
  type(EnergySpectra), intent(in)  :: old_spectra(:,:)
  type(EnergySpectra), intent(in)  :: new_spectra(:)
  type(EnergySpectra), allocatable :: output(:,:)
  
  integer :: n,ialloc
  
  if (size(old_spectra,1)/=size(new_spectra)) then
    call print_line(CODE_ERROR//': Incompatible sizes.')
    call err()
  endif
  
  n = size(old_spectra,2)
  allocate(output(size(new_spectra), n+1), stat=ialloc); call err(ialloc)
  output(:,:n) = old_spectra
  output(:,n+1) = new_spectra
end function

! Appends a column to a potential array.
! output(i,:) = [old_potentials(i,:), new_spectra(i)].
function append_potentials(old_potentials, new_potentials) result(output)
  implicit none
  
  type(PotentialPointer), intent(in)  :: old_potentials(:,:)
  type(PotentialPointer), intent(in)  :: new_potentials(:)
  type(PotentialPointer), allocatable :: output(:,:)
  
  integer :: n,ialloc
  
  if (size(old_potentials,1)/=size(new_potentials)) then
    call print_line(CODE_ERROR//': Incompatible sizes.')
    call err()
  endif
  
  n = size(old_potentials,2)
  allocate(output(size(new_potentials), n+1), stat=ialloc); call err(ialloc)
  output(:,:n) = old_potentials
  output(:,n+1) = new_potentials
end function
end module
