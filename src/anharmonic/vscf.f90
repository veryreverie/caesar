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
  
  public :: run_vscf
  public :: SubspacePotentialAndState
  
  type, extends(NoDefaultConstructor) :: SubspacePotentialAndState
    type(PotentialPointer) :: potential
    type(VscfState)        :: state
  end type
  
  interface SubspacePotentialAndState
    module procedure new_SubspacePotentialAndState
  end interface
contains

impure elemental function new_SubspacePotentialAndState(potential,state) &
   & result(this)
  implicit none
  
  type(PotentialPointer), intent(in) :: potential
  type(VscfState),        intent(in) :: state
  type(SubspacePotentialAndState)    :: this
  
  this%potential = potential
  this%state = state
end function

function run_vscf(potential,basis,energy_convergence,                     &
   & no_converged_calculations,max_pulay_iterations,pre_pulay_iterations, &
   & pre_pulay_damping,anharmonic_data) result(output)
  implicit none
  
  class(PotentialData), intent(in)             :: potential
  type(SubspaceBasis),  intent(in)             :: basis(:)
  real(dp),             intent(in)             :: energy_convergence
  integer,              intent(in)             :: no_converged_calculations
  integer,              intent(in)             :: max_pulay_iterations
  integer,              intent(in)             :: pre_pulay_iterations
  real(dp),             intent(in)             :: pre_pulay_damping
  type(AnharmonicData), intent(in)             :: anharmonic_data
  type(SubspacePotentialAndState), allocatable :: output(:)
  
  type(VscfState),                 allocatable :: state(:)
  type(SubspacePotentialAndState), allocatable :: potentials_and_states(:)
  real(dp),                        allocatable :: energies(:)
  
  type(RealVector), allocatable :: old_coefficients(:)
  type(RealVector), allocatable :: new_coefficients(:)
  
  type(RealVector) :: next_old_coefficients
  
  type(FractionVector), allocatable :: state_wavevectors(:)
  integer                           :: wavevector_changed
  
  integer :: i,j
  
  state = initial_ground_state(basis)
  energies = [0.0_dp]
  old_coefficients = [states_to_coefficients(state)]
  new_coefficients = [RealVector::]
  
  i = 1
  do
    call print_line('VSCF self-consistency step '//i//'.')
    ! Use the current state to calculate single-subspace potentials,
    !    and calculate the new ground states of these potentials.
    potentials_and_states = update(state,basis,potential,anharmonic_data)
    state = potentials_and_states%state
    energies = [energies, sum(potentials_and_states%state%energy)]
    new_coefficients = [new_coefficients, states_to_coefficients(state)]
    
    if (i==1) then
      wavevector_changed = i
    else
      if (any(state%wavevector/=state_wavevectors)) then
        wavevector_changed = i
      endif
    endif
    
    state_wavevectors = state%wavevector
    
    ! Use a damped iterative scheme or a Pulay scheme to converge towards the
    !    self-consistent solution where new coefficients = old coefficients.
    if (i-wavevector_changed<=pre_pulay_iterations) then
      next_old_coefficients = (1-pre_pulay_damping) * old_coefficients(i) &
                          & + pre_pulay_damping * new_coefficients(i)
    else
      j = max(2, i-max_pulay_iterations+1)
      next_old_coefficients = pulay(old_coefficients(j:), new_coefficients(j:))
    endif
    
    ! Normalise coefficients, and append them to the list.
    next_old_coefficients = normalise_coefficients( next_old_coefficients, &
                                                  & state                  )
    old_coefficients = [old_coefficients, next_old_coefficients]
    
    ! Use the next iteration of old_coefficients to update the current state.
    ! N.B. the energy of this state is unknown.
    state = coefficients_to_states(next_old_coefficients,state)
    
    ! Increment the loop counter.
    i = i+1
    
    ! Output the current energy to the terminal.
    call print_line('Ground-state energy: '//energies(i))
    call print_line('')
    
    ! Check whether the energies have converged.
    if (i>=no_converged_calculations) then
      j = i-no_converged_calculations+1
      if (all( abs(energies(j:i-1)-energies(i)) < energy_convergence )) then
        output = potentials_and_states
        exit
      endif
    endif
  enddo
end function

! For each subspace, integrate the potential across the ground states of all
!    other subspaces to give a single-subspace potential, and then diagonalise
!    this potential to give a new ground state.
function update(states,basis,potential,anharmonic_data) result(output)
  implicit none
  
  type(VscfState),      intent(in)             :: states(:)
  type(SubspaceBasis),  intent(in)             :: basis(:)
  class(PotentialData), intent(in)             :: potential
  type(AnharmonicData), intent(in)             :: anharmonic_data
  type(SubspacePotentialAndState), allocatable :: output(:)
  
  type(StructureData) :: supercell
  
  type(DegenerateSubspace) :: subspace
  
  type(PolynomialState),  allocatable :: subspace_states(:)
  type(PotentialPointer), allocatable :: subspace_potentials(:)
  
  type(SymmetricEigenstuff), allocatable :: wavevector_ground_states(:)
  
  real(dp), allocatable :: hamiltonian(:,:)
  
  type(HarmonicState) :: bra
  type(HarmonicState) :: ket
  
  type(SymmetricEigenstuff), allocatable :: estuff(:)
  
  type(VscfState), allocatable :: ground_states(:)
  
  integer :: i,j,k,l,m,ialloc
  
  supercell = anharmonic_data%anharmonic_supercell
  
  ! Generate the single-subspace potentials {V_i}, defined as
  !    V_i = (prod_{j/=i}<j|)V(prod_{j/=i}|j>).
  subspace_states = PolynomialState(states,basis)
  call print_line('Generating single-subspace potentials.')
  subspace_potentials = generate_subspace_potentials( potential,       &
                                                    & subspace_states, &
                                                    & anharmonic_data  )
  
  call print_line('Generating single-subspace ground states.')
  allocate( ground_states(size(states)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(states)
    subspace = anharmonic_data%degenerate_subspaces(i)
    
    ! Calculate the ground state at each wavevector.
    allocate( wavevector_ground_states(size(basis(i))), &
            & stat=ialloc); call err(ialloc)
    do j=1,size(basis(i))
      ! Calculate the Hamiltonian in the harmonic basis.
      allocate( hamiltonian( size(basis(i)%wavevectors(j)),    &
              &              size(basis(i)%wavevectors(j))  ), &
              & stat=ialloc); call err(ialloc)
      hamiltonian = 0.0_dp
      do k=1,size(basis(i)%wavevectors(j))
        bra = basis(i)%wavevectors(j)%harmonic_states(k)
        do l=1,size(basis(i)%wavevectors(j)%harmonic_couplings(k))
          m = basis(i)%wavevectors(j)%harmonic_couplings(k)%id(l)
          ket = basis(i)%wavevectors(j)%harmonic_states(m)
          hamiltonian(k,m) = subspace_potentials(i)%potential_energy( &
                         &                           bra,             &
                         &                           ket,             &
                         &                           anharmonic_data) &
                         & + kinetic_energy(bra, ket, subspace, supercell)
        enddo
      enddo
      
      ! Diagonalise Hamiltonian to get new ground state.
      estuff = diagonalise_symmetric(hamiltonian)
      
      ! Store the ground state at this wavevector in the
      !    wavevector-by-wavevector array.
      wavevector_ground_states(j) = estuff(1)
      
      deallocate(hamiltonian, stat=ialloc); call err(ialloc)
    enddo
    
    ! Find the wavevector with the lowest energy ground-state.
    j = minloc(wavevector_ground_states%eval,1)
    
    ground_states(i) = VscfState(                           &
       & subspace_id  = basis(i)%subspace_id,               &
       & wavevector   = basis(i)%wavevectors(j)%wavevector, &
       & degeneracy   = basis(i)%wavevectors(j)%degeneracy, &
       & energy       = wavevector_ground_states(j)%eval,   &
       & coefficients = wavevector_ground_states(j)%evec    )
    
    deallocate(wavevector_ground_states, stat=ialloc); call err(ialloc)
  enddo
  
  output = SubspacePotentialAndState( subspace_potentials, &
                                    & ground_states        )
end function

! Concatenates the coefficients of all the states together,
!    to allow the Pulay scheme to be called.
function states_to_coefficients(states) result(output)
  implicit none
  
  type(VscfState), intent(in) :: states(:)
  type(RealVector)            :: output
  
  integer :: i
  
  output = vec([( states(i)%coefficients, i=1, size(states) )])
end function

! Reverses coefficients(), generating a set of states from a vector of
!    coefficients.
! Takes an array of states to act as a template for the output states.
function coefficients_to_states(input,states) result(output)
  implicit none
  
  type(RealVector), intent(in)  :: input
  type(VscfState),  intent(in)  :: states(:)
  type(VscfState),  allocatable :: output(:)
  
  real(dp), allocatable :: coefficients(:)
  real(dp), allocatable :: state_coefficients(:)
  
  integer :: i,j,ialloc
  
  coefficients = dble(input)
  
  j = 0
  allocate(output(size(states)), stat=ialloc); call err(ialloc)
  do i=1,size(states)
    state_coefficients = coefficients(j+1:j+size(states(i)%coefficients))
    j = j+size(state_coefficients)
    output(i) = VscfState( subspace_id  = states(i)%subspace_id, &
                         & wavevector   = states(i)%wavevector,  &
                         & degeneracy   = states(i)%degeneracy,  &
                         & energy       = 0.0_dp,                &
                         & coefficients = state_coefficients     )
  enddo
end function

! Normalises a set of coefficients, using a set of states to identify which
!    coefficients correspond to which state.
function normalise_coefficients(input,states) result(output)
  implicit none
  
  type(RealVector), intent(in) :: input
  type(VscfState),  intent(in) :: states(:)
  type(RealVector)             :: output
  
  real(dp), allocatable :: coefficients(:)
  real(dp), allocatable :: state_coefficients(:)
  real(dp), allocatable :: new_coefficients(:)
  
  integer :: i,j
  
  coefficients = dble(input)
  new_coefficients = [real::]
  do i=1,size(states)
    j = size(new_coefficients)
    state_coefficients = coefficients(j+1:j+size(states(i)%coefficients))
    state_coefficients = state_coefficients / l2_norm(state_coefficients)
    new_coefficients = [new_coefficients, state_coefficients]
  enddo
  output = new_coefficients
end function
end module
