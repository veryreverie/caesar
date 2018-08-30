! ======================================================================
! Run VSCF on a potential to generate VSCF ground states.
! ======================================================================
module vscf_module
  use common_module
  
  use states_module
  
  use anharmonic_data_module
  use potential_module
  use potential_pointer_module
  implicit none
  
  private
  
  type StateAndEnergy
    type(VscfGroundState), allocatable :: state(:)
    real(dp)                           :: energy
  end type
  
  public :: run_vscf
contains

function run_vscf(potential,basis,energy_convergence,max_pulay_iterations, &
   & pre_pulay_iterations,pre_pulay_damping,anharmonic_data) result(output)
  implicit none
  
  class(PotentialData), intent(in)   :: potential
  type(SubspaceBasis),  intent(in)   :: basis(:)
  real(dp),             intent(in)   :: energy_convergence
  integer,              intent(in)   :: max_pulay_iterations
  integer,              intent(in)   :: pre_pulay_iterations
  real(dp),             intent(in)   :: pre_pulay_damping
  type(AnharmonicData), intent(in)   :: anharmonic_data
  type(VscfGroundState), allocatable :: output(:)
  
  type(VscfGroundState), allocatable :: old_guess(:)
  type(VscfGroundState), allocatable :: new_guess(:)
  
  type(VscfGroundState), allocatable :: state(:)
  type(StateAndEnergy)               :: state_and_energy
  real(dp),              allocatable :: energies(:)
  
  type(RealVector), allocatable :: old_coefficients(:)
  type(RealVector), allocatable :: new_coefficients(:)
  
  type(RealVector) :: next_old_coefficients
  
  integer :: i,j
  
  state = initial_ground_state(basis)
  energies = [0.0_dp]
  old_coefficients = [states_to_coefficients(state)]
  new_coefficients = [RealVector::]
  
  i = 1
  iter: do
    ! Use the current state to calculate single-subspace potentials,
    !    and calculate the new ground states of these potentials.
    state_and_energy = update(state,basis,potential,anharmonic_data)
    state = state_and_energy%state
    energies = [energies, state_and_energy%energy]
    new_coefficients = [new_coefficients, states_to_coefficients(state)]
    
    ! Use a damped iterative scheme or a Pulay scheme to converge towards the
    !    self-consistent solution where new coefficients = old coefficients.
    if (i<=pre_pulay_iterations) then
      next_old_coefficients = (1-pre_pulay_damping) * old_coefficients(i) &
                          & + pre_pulay_damping * new_coefficients(i)
    else
      j = max(2, i-max_pulay_iterations+1)
      next_old_coefficients = pulay(old_coefficients(j:), new_coefficients(j:))
    endif
    
    old_coefficients = [old_coefficients, next_old_coefficients]
    
    ! Increment the loop counter.
    i = i+1
    
    ! Output the current energy to the terminal.
    call print_line('Ground-state energy at self-consistency step '//i//' :&
                   & '//energies(i))
    
    ! Check whether the energies have converged.
    ! If the energies haven't converged, return to the top of the loop.
    if (i<=5) then
      cycle iter
    else
      do j=1,5
        if (energies(i)-energies(i-j)>energy_convergence) then
          cycle iter
        endif
      enddo
    endif
    
    ! If the frequencies are converged, break out of the loop.
    output = state
    exit iter
  enddo iter
end function

! For each subspace, integrate the potential across the ground states of all
!    other subspaces to give a single-subspace potential, and then diagonalise
!    this potential to give a new ground state.
function update(states,basis,potential,anharmonic_data) result(output)
  implicit none
  
  type(VscfGroundState), intent(in)  :: states(:)
  type(SubspaceBasis),   intent(in)  :: basis(:)
  class(PotentialData),  intent(in)  :: potential
  type(AnharmonicData),  intent(in)  :: anharmonic_data
  type(StateAndEnergy)               :: output
  
  type(StructureData) :: supercell
  
  type(DegenerateSubspace) :: subspace
  
  type(PotentialPointer) :: new_potential
  
  type(SymmetricEigenstuff), allocatable :: wavevector_ground_states(:)
  
  type(PotentialPointer) :: new_new_potential
  
  real(dp), allocatable :: hamiltonian(:,:)
  
  type(SubspaceState) :: bra
  type(SubspaceState) :: ket
  
  type(SymmetricEigenstuff), allocatable :: estuff(:)
  
  integer :: i,j,k,l,ialloc
  
  supercell = anharmonic_data%anharmonic_supercell
  
  allocate(output%state(size(states)), stat=ialloc); call err(ialloc)
  output%energy = 0.0_dp
  do i=1,size(states)
    subspace = anharmonic_data%degenerate_subspaces(i)
    
    ! Integrate the potential across all other subspaces, to give the
    !    single-subspace potential.
    new_potential = potential
    do j=1,size(states)
      if (j/=i) then
        call new_potential%braket( SumState(states(j),basis(j)), &
                                 & SumState(states(j),basis(j)), &
                                 & anharmonic_data               )
      endif
    enddo
    call new_potential%zero_energy()
    
    ! Calculate the ground state at each wavevector.
    allocate( wavevector_ground_states(size(basis(i))), &
            & stat=ialloc); call err(ialloc)
    do j=1,size(basis(i))
      ! Calculate the Hamiltonian in terms of the (non-orthonormal) states.
      allocate( hamiltonian( size(basis(i)%wavevectors(j)),    &
              &              size(basis(i)%wavevectors(j))  ), &
              & stat=ialloc); call err(ialloc)
      do k=1,size(basis(i)%wavevectors(j))
        do l=1,size(basis(i)%wavevectors(j))
          bra = basis(i)%wavevectors(j)%states(k)
          ket = basis(i)%wavevectors(j)%states(l)
          new_new_potential = new_potential
          call new_new_potential%braket( bra,            &
                                       & ket,            &
                                       & anharmonic_data )
          hamiltonian(l,k) = new_new_potential%undisplaced_energy() &
                         & + kinetic_energy( bra,                   &
                         &                   ket,                   &
                         &                   subspace,              &
                         &                   supercell )
        enddo
      enddo
      
      ! Convert the Hamiltonian into an orthonormal basis.
      hamiltonian = dble( basis(i)%wavevectors(j)%states_to_basis            &
                      & * mat(hamiltonian)                                   &
                      & * transpose(basis(i)%wavevectors(j)%states_to_basis) )
      
      ! Diagonalise Hamiltonian to get new ground state.
      estuff = diagonalise_symmetric(hamiltonian)
      
      ! Store the ground state at this wavevector in the
      !    wavevector-by-wavevector array.
      wavevector_ground_states(j) = estuff(1)
      
      deallocate(hamiltonian, stat=ialloc); call err(ialloc)
    enddo
    
    ! Find the wavevector with the lowest energy ground-state.
    j = minloc(wavevector_ground_states%eval,1)
    
    if (.not. is_int(basis(i)%wavevectors(j)%wavevector)) then
      call print_line(WARNING//': Ground state has gained non-zero &
         &wavevector. This will be ignored in favour of the zero-wavevector &
         &ground state.')
      j = first(is_int(basis(i)%wavevectors%wavevector))
    endif
    
    output%state(i) = VscfGroundState(                      &
       & subspace_id  = basis(i)%subspace_id,               &
       & wavevector   = basis(i)%wavevectors(j)%wavevector, &
       & coefficients = wavevector_ground_states(j)%evec    )
    output%energy = output%energy + wavevector_ground_states(j)%eval
    
    deallocate(wavevector_ground_states, stat=ialloc); call err(ialloc)
  enddo
end function

! Concatenates the coefficients of all the states together,
!    to allow the Pulay scheme to be called.
function states_to_coefficients(states) result(output)
  implicit none
  
  type(VscfGroundState), intent(in) :: states(:)
  type(RealVector)                  :: output
  
  integer :: i
  
  if (any(.not.is_int(states%wavevector))) then
    call print_line(ERROR//': The ground state has gained a finite &
       &wavevector. Caesar is at present unable to handle this.')
    call err()
  endif
  
  output = vec([( states(i)%coefficients, i=1, size(states) )])
end function

! Reverses coefficients(), generating a set of states from a vector of
!    coefficients.
! Takes an array of states to act as a template for the output states.
function coefficients_to_states(input,states) result(output)
  implicit none
  
  type(RealVector),      intent(in)  :: input
  type(VscfGroundState), intent(in)  :: states(:)
  type(VscfGroundState), allocatable :: output(:)
  
  real(dp), allocatable :: coefficients(:)
  real(dp), allocatable :: state_coefficients(:)
  
  integer :: i,j,ialloc
  
  coefficients = dble(input)
  
  j = 0
  allocate(output(size(states)), stat=ialloc); call err(ialloc)
  do i=1,size(states)
    state_coefficients = coefficients(j+1:j+size(states(i)%coefficients))
    j = j+size(state_coefficients)
    output(i) = VscfGroundState( subspace_id  = states(i)%subspace_id, &
                               & wavevector   = states(i)%wavevector,  &
                               & coefficients = state_coefficients     )
  enddo
end function
end module
