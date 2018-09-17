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
    type(VscfGroundState)  :: state
    real(dp)               :: energy
  end type
  
  interface SubspacePotentialAndState
    module procedure new_SubspacePotentialAndState
  end interface
contains

impure elemental function new_SubspacePotentialAndState(potential,state, &
   & energy) result(this)
  implicit none
  
  type(PotentialPointer), intent(in) :: potential
  type(VscfGroundState),  intent(in) :: state
  real(dp),               intent(in) :: energy
  type(SubspacePotentialAndState)    :: this
  
  this%potential = potential
  this%state = state
  this%energy = energy
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
  
  type(VscfGroundState), allocatable :: old_guess(:)
  type(VscfGroundState), allocatable :: new_guess(:)
  
  type(VscfGroundState),           allocatable :: state(:)
  type(SubspacePotentialAndState), allocatable :: potentials_and_states(:)
  real(dp),                        allocatable :: energies(:)
  
  type(RealVector), allocatable :: old_coefficients(:)
  type(RealVector), allocatable :: new_coefficients(:)
  
  type(RealVector) :: next_old_coefficients
  
  integer :: i,j
  
  state = initial_ground_state(basis)
  energies = [0.0_dp]
  old_coefficients = [states_to_coefficients(state)]
  new_coefficients = [RealVector::]
  
  i = 1
  do
    ! Use the current state to calculate single-subspace potentials,
    !    and calculate the new ground states of these potentials.
    potentials_and_states = update(state,basis,potential,anharmonic_data)
    state = potentials_and_states%state
    energies = [energies, sum(potentials_and_states%energy)]
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
  
  type(VscfGroundState), intent(in)            :: states(:)
  type(SubspaceBasis),   intent(in)            :: basis(:)
  class(PotentialData),  intent(in)            :: potential
  type(AnharmonicData),  intent(in)            :: anharmonic_data
  type(SubspacePotentialAndState), allocatable :: output(:)
  
  type(StructureData) :: supercell
  
  type(DegenerateSubspace) :: subspace
  
  type(PolynomialState),  allocatable :: subspace_states(:)
  type(PotentialPointer), allocatable :: subspace_potentials(:)
  
  type(SymmetricEigenstuff), allocatable :: wavevector_ground_states(:)
  
  real(dp), allocatable :: hamiltonian(:,:)
  
  type(MonomialState) :: bra
  type(MonomialState) :: ket
  
  type(SymmetricEigenstuff), allocatable :: estuff(:)
  
  type(VscfGroundState), allocatable :: ground_states(:)
  real(dp),              allocatable :: energies(:)
  
  integer :: i,j,k,l,m,ialloc
  
  supercell = anharmonic_data%anharmonic_supercell
  
  ! Generate the single-subspace potentials {V_i}, defined as
  !    V_i = (prod_{j/=i}<j|)V(prod_{j/=i}|j>).
  subspace_states = PolynomialState(states,basis)
  subspace_potentials = generate_subspace_potentials( potential,       &
                                                    & subspace_states, &
                                                    & anharmonic_data  )
  
  allocate( ground_states(size(states)), &
          & energies(size(states)),      &
          & stat=ialloc); call err(ialloc)
  do i=1,size(states)
    subspace = anharmonic_data%degenerate_subspaces(i)
    
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
          hamiltonian(l,k) = subspace_potentials(i)%potential_energy(   &
                         &                            bra,              &
                         &                            ket,              &
                         &                            anharmonic_data ) &
                         & + kinetic_energy( bra,                       &
                         &                   ket,                       &
                         &                   subspace,                  &
                         &                   supercell )
        enddo
      enddo
      
      ! Convert the Hamiltonian into an orthonormal basis.
      hamiltonian = &
         & basis(i)%wavevectors(j)%operator_states_to_basis(hamiltonian)
      
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
    
    ground_states(i) = VscfGroundState(                     &
       & subspace_id  = basis(i)%subspace_id,               &
       & wavevector   = basis(i)%wavevectors(j)%wavevector, &
       & coefficients = wavevector_ground_states(j)%evec    )
    energies(i) = wavevector_ground_states(j)%eval
    
    deallocate(wavevector_ground_states, stat=ialloc); call err(ialloc)
  enddo
  
  output = SubspacePotentialAndState( subspace_potentials, &
                                    & ground_states,       &
                                    & energies             )
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
