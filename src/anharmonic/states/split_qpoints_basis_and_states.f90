! ======================================================================
! A basis of states which treats q-points separately.
! ======================================================================
! N.B. only states at a single q-point are calculated and stored.
! All other q-points are included by symmetry.
module split_qpoints_basis_and_states_module
  use common_module
  
  use anharmonic_common_module
  
  use state_helper_module
  use monomial_state_module
  use harmonic_state_module
  use polynomial_state_module
  use state_conversion_module
  use wavevector_basis_module
  use split_qpoints_wavefunctions_module
  implicit none
  
  private
  
  public :: startup_split_qpoints_basis_and_states
  
  public :: SplitQpointsBasis
  public :: SplitQpointsState
  public :: SplitQpointsStates
  
  public :: size
  public :: initial_ground_state
  
  type, extends(NoDefaultConstructor) :: QpointModes
    type(QpointData)               :: qpoint
    type(QpointData)               :: paired_qpoint
    type(ComplexMode), allocatable :: modes(:)
    type(ComplexMode), allocatable :: paired_modes(:)
  end type
  
  interface QpointModes
    module procedure new_QpointModes
  end interface
  
  ! All states spanning the subspace.
  type, extends(SubspaceBasis) :: SplitQpointsBasis
    ! The maximum power of the monomial states.
    ! This is also the maximum occupation of the harmonic basis states.
    integer  :: maximum_power
    ! The expansion order of the potential.
    ! This is also the limit on coupling between states.
    integer :: expansion_order
    ! The ID and frequency of the subspace.
    integer  :: subspace_id
    real(dp) :: frequency
    ! The states, wavevector by wavevector.
    ! N.B. this only includes one wavevector from each symmetry-related set.
    type(WavevectorBasis), allocatable :: wavevectors(:)
    ! The q-points, and the modes at each q-point.
    type(QpointModes), allocatable :: qpoint_modes(:)
  contains
    procedure, public, nopass :: representation => &
                               & representation_SplitQpointsBasis
    
    ! Set the frequency of the basis.
    procedure, public :: set_frequency => set_frequency_SplitQpointsBasis
    
    ! Print the harmonic ground-state wavefunction of the basis.
    procedure, public :: ground_state_wavefunction
    
    ! Generate the first guess at states.
    procedure, public :: initial_states => initial_states_SplitQpointsBasis
    
    ! Generate the eigenstates of a single-subspace potential.
    procedure, public :: calculate_states => calculate_states_SplitQpointsBasis
    
    ! Generate the eigenstates of a single-subspace potential,
    !    at a single q-point.
    procedure, private :: calculate_split_states => &
                        & calculate_split_states_SplitQpointsBasis
    
    ! I/O.
    procedure, public :: read  => read_SplitQpointsBasis
    procedure, public :: write => write_SplitQpointsBasis
  end type
  
  interface SplitQpointsBasis
    module procedure new_SplitQpointsBasis
    module procedure new_SplitQpointsBasis_SubspaceBasis
    module procedure new_SplitQpointsBasis_subspace
    module procedure new_SplitQpointsBasis_Strings
    module procedure new_SplitQpointsBasis_StringArray
  end interface
  
  interface size
    module procedure size_SplitQpointsBasis
  end interface
  
  type, extends(SubspaceState) :: SplitQpointsState
    type(FractionVector)  :: wavevector
    integer               :: degeneracy
    real(dp)              :: energy
    real(dp), allocatable :: coefficients(:)
  contains
    procedure, public, nopass :: representation => &
                               & representation_SplitQpointsState
    
    procedure, public :: braket_SubspaceState => &
                       & braket_SubspaceState_SplitQpointsState
    procedure, public :: braket_ComplexUnivariate => &
                       & braket_ComplexUnivariate_SplitQpointsState
    procedure, public :: braket_ComplexMonomial => &
                       & braket_ComplexMonomial_SplitQpointsState
    procedure, public :: kinetic_energy => &
                       & kinetic_energy_SplitQpointsState
    procedure, public :: harmonic_potential_energy => &
                       & harmonic_potential_energy_SplitQpointsState
    procedure, public :: kinetic_stress => &
                       & kinetic_stress_SplitQpointsState
    
    procedure, public :: wavefunction => wavefunction_SplitQpointsState
    
    ! I/O.
    procedure, public :: read  => read_SplitQpointsState
    procedure, public :: write => write_SplitQpointsState
  end type
  
  interface SplitQpointsState
    module procedure new_SplitQpointsState
    module procedure new_SplitQpointsState_SubspaceState
    module procedure new_SplitQpointsState_Strings
    module procedure new_SplitQpointsState_StringArray
  end interface
  
  interface PolynomialState
    module procedure new_PolynomialState_SplitQpointsState
  end interface
  
  type, extends(SubspaceStates) :: SplitQpointsStates
    type(SplitQpointsState), allocatable :: vscf_states(:)
  contains
    procedure, public, nopass :: representation => &
                               & representation_SplitQpointsStates
    
    procedure, public :: spectra => spectra_SplitQpointsStates
    procedure, public :: wavefunctions => wavefunctions_SplitQpointsStates
    procedure, public :: integrate_potential => &
                       & integrate_potential_SplitQpointsStates
    procedure, public :: integrate_stress => &
                       & integrate_stress_SplitQpointsStates
    
    ! I/O.
    procedure, public :: read  => read_SplitQpointsStates
    procedure, public :: write => write_SplitQpointsStates
  end type
  
  interface SplitQpointsStates
    module procedure new_SplitQpointsStates
    module procedure new_SplitQpointsStates_SubspaceStates
    module procedure new_SplitQpointsStates_Strings
    module procedure new_SplitQpointsStates_StringArray
  end interface
contains

! Startup procedure.
subroutine startup_split_qpoints_basis_and_states()
  implicit none
  
  type(SplitQpointsBasis)  :: basis
  type(SplitQpointsState)  :: state
  type(SplitQpointsStates) :: states
  
  call basis%startup()
  call state%startup()
  call states%startup()
end subroutine

! QpointModes methods.
function new_QpointModes(qpoint,paired_qpoint,modes,paired_modes) &
   & result(this)
  implicit none
  
  type(QpointData),  intent(in) :: qpoint
  type(QpointData),  intent(in) :: paired_qpoint
  type(ComplexMode), intent(in) :: modes(:)
  type(ComplexMode), intent(in) :: paired_modes(:)
  type(QpointModes)             :: this
  
  this%qpoint = qpoint
  this%paired_qpoint = paired_qpoint
  this%modes = modes(sort(modes%id))
  this%paired_modes = paired_modes(sort(paired_modes%paired_id))
end function

! ----------------------------------------------------------------------
! SplitQpointsBasis methods.
! ----------------------------------------------------------------------
! Constructors and size functions.
function new_SplitQpointsBasis(maximum_power,expansion_order,subspace_id, &
   & frequency,wavevectors,qpoint_modes) result(this)
  implicit none
  
  integer,               intent(in) :: maximum_power
  integer,               intent(in) :: expansion_order
  integer,               intent(in) :: subspace_id
  real(dp),              intent(in) :: frequency
  type(WavevectorBasis), intent(in) :: wavevectors(:)
  type(QpointModes),     intent(in) :: qpoint_modes(:)
  type(SplitQpointsBasis)           :: this
  
  this%maximum_power   = maximum_power
  this%expansion_order = expansion_order
  this%subspace_id     = subspace_id
  this%frequency       = frequency
  this%wavevectors     = wavevectors
  this%qpoint_modes    = qpoint_modes
end function

recursive function new_SplitQpointsBasis_SubspaceBasis(input) result(this)
  implicit none
  
  class(SubspaceBasis), intent(in) :: input
  type(SplitQpointsBasis)          :: this
  
  select type(input); type is(SplitQpointsBasis)
    this = input
  type is(SubspaceBasisPointer)
    this = SplitQpointsBasis(input%basis())
  class default
    call err()
  end select
end function

function size_SplitQpointsBasis(this) result(output)
  implicit none
  
  type(SplitQpointsBasis), intent(in) :: this
  integer                             :: output
  
  output = size(this%wavevectors)
end function

! Type representation.
impure elemental function representation_SplitQpointsBasis() result(output)
  implicit none
  
  type(String) :: output
  
  output = 'split_qpoints'
end function

! Set the frequency of the basis.
impure elemental subroutine set_frequency_SplitQpointsBasis(this,frequency)
  implicit none
  
  class(SplitQpointsBasis), intent(inout) :: this
  real(dp),                 intent(in)    :: frequency
  
  this%frequency = frequency
  call this%wavevectors%set_frequency(frequency)
end subroutine

! Generates states up to a given power, spanning the whole subspace.
function new_SplitQpointsBasis_subspace(subspace,frequency,modes,qpoints, &
   & supercell,maximum_power,potential_expansion_order,symmetries)        &
   & result(output)
  implicit none
  
  type(DegenerateSubspace), intent(in) :: subspace
  real(dp),                 intent(in) :: frequency
  type(ComplexMode),        intent(in) :: modes(:)
  type(QpointData),         intent(in) :: qpoints(:)
  type(StructureData),      intent(in) :: supercell
  integer,                  intent(in) :: maximum_power
  integer,                  intent(in) :: potential_expansion_order
  type(SymmetryOperator),   intent(in) :: symmetries(:)
  type(SplitQpointsBasis)              :: output
  
  type(ComplexMode), allocatable :: subspace_modes(:)
  type(QpointData),  allocatable :: subspace_qpoints(:)
  type(QpointData),  allocatable :: qpoint_set(:)
  
  type(ComplexMode), allocatable :: basis_modes(:,:)
  
  type(QpointModes), allocatable :: qpoint_modes(:)
  
  type(WavevectorBasis), allocatable :: wavevectors(:)
  
  integer :: i,j,ialloc
  
  ! List the modes and corresponding q-points in the subspace.
  subspace_modes = subspace%modes(modes)
  subspace_qpoints = subspace%qpoints(modes, qpoints)
  
  qpoint_set = subspace_qpoints(set(subspace_qpoints%id))
  qpoint_set = qpoint_set(filter(qpoint_set%id<=qpoint_set%paired_qpoint_id))
  
  allocate(qpoint_modes(size(qpoint_set)), stat=ialloc); call err(ialloc)
  do i=1,size(qpoint_modes)
    qpoint_modes(i) = QpointModes(                                      &
       & qpoint_set(i),                                                 &
       & subspace_qpoints(first(                                        &
       &     subspace_qpoints%id==qpoint_set(i)%paired_qpoint_id )),    &
       & subspace_modes(filter(                                         &
       &    subspace_modes%qpoint_id==qpoint_set(i)%id )),              &
       & subspace_modes(filter(                                         &
       &    subspace_modes%qpoint_id==qpoint_set(i)%paired_qpoint_id )) )
  enddo
  
  ! Generate a basis at the first q-point only.
  wavevectors = WavevectorBasis( subspace,                  &
                               & frequency,                 &
                               & modes,                     &
                               & qpoints,                   &
                               & maximum_power,             &
                               & potential_expansion_order, &
                               & symmetries,                &
                               & qpoint_modes(1)%qpoint     )
  
  output = SplitQpointsBasis( maximum_power,             &
                            & potential_expansion_order, &
                            & subspace%id,               &
                            & frequency,                 &
                            & wavevectors,               &
                            & qpoint_modes               )
end function

! Returns the harmonic ground-state wavefunction for the basis.
function ground_state_wavefunction(this,subspace,supercell) result(output)
  implicit none
  
  class(SplitQpointsBasis), intent(in) :: this
  type(DegenerateSubspace), intent(in) :: subspace
  type(StructureData),      intent(in) :: supercell
  type(String)                         :: output
  
  real(dp) :: mass
  
  real(dp)                  :: coefficient
  type(String), allocatable :: terms(:)
  
  integer :: i
  
  ! Calculate the (geometric) average mass.
  mass = product(supercell%atoms%mass())
  mass = mass**(1.0_dp/size(supercell%atoms))
  
  ! Calculate the coefficient.
  coefficient = 1
  terms = [String::]
  do i=1,size(subspace)
    if (subspace%mode_ids(i)==subspace%paired_ids(i)) then
      ! |0_i> = sqrt(sqrt(m*w/pi)) exp(- 1/2 N w (u_i)^2 )
      coefficient = coefficient * (mass*this%frequency/PI)**0.25_dp
      terms = [terms, 'u'//subspace%mode_ids(i)//'^2']
    elseif (subspace%mode_ids(i)<subspace%paired_ids(i)) then
      ! |0_i,0_j> = sqrt(2*m*w/pi) exp(- N w |u_i|^2 )
      coefficient = coefficient * (2*mass*this%frequency/PI)**0.5_dp
      terms = [ terms,                                                    &
              & '2*u'//subspace%mode_ids(i)//'*u'//subspace%paired_ids(i) ]
    endif
  enddo
  
  output = coefficient//'*e^('               // &
     & (-supercell%sc_size*this%frequency/2) // &
     & '*('//join(terms,'+')//'))'
end function

! Generate an initial guess at states.
impure elemental function initial_states_SplitQpointsBasis(this,subspace, &
   & anharmonic_data) result(output)
  implicit none
  
  class(SplitQpointsBasis), intent(in) :: this
  type(DegenerateSubspace), intent(in) :: subspace
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(SubspaceStatesPointer)          :: output
  
  type(SplitQpointsState) :: ground_state
  
  integer :: ialloc
  
  ! Generate the state |0>.
  ground_state = initial_ground_state(this)
  
  ! Generate the set of states {|0>}.
  output = SubspaceStatesPointer(SplitQpointsStates([ground_state]))
end function

! Generate initial guess. This is simply the basis state |0>, i.e. the
!    ground state of the effective harmonic potential from which the basis
!    states were generated.
impure elemental function initial_ground_state(basis) result(output)
  implicit none
  
  type(SplitQpointsBasis), intent(in) :: basis
  type(SplitQpointsState)             :: output
  
  type(WavevectorBasis) :: wavevector_basis
  
  real(dp), allocatable :: coefficients(:)
  
  integer :: i
  
  ! Find the wavevector [0,0,0].
  wavevector_basis = basis%wavevectors(            &
     & first(is_int(basis%wavevectors%wavevector)) )
  
  ! Construct the coefficient vector in the basis of monomial states.
  ! All coefficients are zero, except for the coefficient of |0>, which is one.
  coefficients = [( 0.0_dp, i=1, size(wavevector_basis) )]
  coefficients(                                                      &
     & first(wavevector_basis%harmonic_states%total_occupation()==0) ) = 1
  
  ! Convert the coefficients into the orthonormal basis.
  coefficients = wavevector_basis%coefficients_states_to_basis(coefficients)
  
  ! Construct output.
  output = SplitQpointsState(                      &
     & subspace_id  = basis%subspace_id,           &
     & wavevector   = wavevector_basis%wavevector, &
     & degeneracy   = wavevector_basis%degeneracy, &
     & energy       = 0.0_dp,                      &
     & coefficients = coefficients                 )
end function

! Calculate the eigenstates of a single-subspace potential.
impure elemental function calculate_states_SplitQpointsBasis(this,subspace, &
   & subspace_potential,energy_convergence,no_converged_calculations,       &
   & max_pulay_iterations,pre_pulay_iterations,pre_pulay_damping,           &
   & anharmonic_data) result(output)
  implicit none
  
  class(SplitQpointsBasis), intent(in) :: this
  type(DegenerateSubspace), intent(in) :: subspace
  class(PotentialData),     intent(in) :: subspace_potential
  real(dp),                 intent(in) :: energy_convergence
  integer,                  intent(in) :: no_converged_calculations
  integer,                  intent(in) :: max_pulay_iterations
  integer,                  intent(in) :: pre_pulay_iterations
  real(dp),                 intent(in) :: pre_pulay_damping
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(SubspaceStatesPointer)          :: output
  
  output = SubspaceStatesPointer(split_vscf( subspace_potential,        &
                                           & subspace,                  &
                                           & this,                      &
                                           & energy_convergence,        &
                                           & no_converged_calculations, &
                                           & max_pulay_iterations,      &
                                           & pre_pulay_iterations,      &
                                           & pre_pulay_damping,         &
                                           & anharmonic_data            ))
end function

function split_vscf(potential,subspace,basis,energy_convergence,          &
   & no_converged_calculations,max_pulay_iterations,pre_pulay_iterations, &
   & pre_pulay_damping,anharmonic_data) result(output)
  implicit none
  
  class(PotentialData),     intent(in) :: potential
  type(DegenerateSubspace), intent(in) :: subspace
  type(SplitQpointsBasis),  intent(in) :: basis
  real(dp),                 intent(in) :: energy_convergence
  integer,                  intent(in) :: no_converged_calculations
  integer,                  intent(in) :: max_pulay_iterations
  integer,                  intent(in) :: pre_pulay_iterations
  real(dp),                 intent(in) :: pre_pulay_damping
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(SplitQpointsStates)             :: output
  
  type(SplitQPointsStates) :: initial_states
  type(SplitQpointsState)  :: initial_state
  type(PotentialPointer)   :: initial_potential
  
  type(PotentialPointer),   allocatable :: input_potentials(:)
  type(SplitQpointsStates), allocatable :: states(:)
  type(SplitQpointsState)               :: state
  type(PotentialPointer),   allocatable :: output_potentials(:)
  
  integer :: first_pulay_step
  
  integer :: i,j
  
  initial_states = SplitQpointsStates(basis%initial_states( subspace,       &
                                                          & anharmonic_data ))
  initial_state = initial_states%vscf_states(      &
     & minloc(initial_states%vscf_states%energy,1) )
  initial_potential = PotentialPointer(potential)
  do i=2,size(basis%qpoint_modes)
    initial_potential = braket( initial_state,               &
                              & initial_potential,           &
                              & subspace,                    &
                              & basis,                       &
                              & anharmonic_data,             &
                              & basis%qpoint_modes(i)%qpoint )
  enddo
  call initial_potential%zero_energy()
  
  input_potentials = [initial_potential]
  states = [SplitQpointsStates::]
  output_potentials = [PotentialPointer::]
  i = 1
  do
    ! Calculate new states.
    states = [ states,                                                    &
             & basis%calculate_split_states( subspace,                    &
             &                               input_potentials(i),         &
             &                               anharmonic_data            ) ]
    state = states(i)%vscf_states(minloc(states(i)%vscf_states%energy,1))
    
    ! Use states to calculate new potential.
    output_potentials = [output_potentials, PotentialPointer(potential)]
    do j=2,size(basis%qpoint_modes)
      output_potentials(i) = braket( state,                       &
                                   & output_potentials(i),        &
                                   & subspace,                    &
                                   & basis,                       &
                                   & anharmonic_data,             &
                                   & basis%qpoint_modes(i)%qpoint )
    enddo
    call output_potentials(i)%zero_energy()
    
    ! Check for convergence.
    if (i>no_converged_calculations) then
      if (all(steps_converged(                      &
         & states(i),                               &
         & states(i-no_converged_calculations:i-1), &
         & energy_convergence                       ))) then
      endif
      output = states(i)
      exit
    endif
    
    ! Check for over-convergence.
    ! This is needed to avoid numerical problems with the Pulay scheme.
    if (i>1) then
      if (steps_converged(states(i),states(i-1),energy_convergence/100)) then
        output = states(i)
        exit
      endif
    endif
    
    ! If convergence has not been reached, generate the next input potential.
    if (i<=pre_pulay_iterations) then
      input_potentials = [                                         &
         & input_potentials,                                       &
         & PotentialPointer(input_potentials(i)%iterate_damped(    &
         &                                output_potentials(i),    &
         &                                pre_pulay_damping,       &
         &                                anharmonic_data       )) ]
    else
      first_pulay_step = max(1,i-max_pulay_iterations+1)
      input_potentials = [                                        &
         & input_potentials,                                      &
         & PotentialPointer(input_potentials(i)%iterate_pulay(    &
         &              input_potentials(first_pulay_step:i),     &
         &              output_potentials(first_pulay_step:i),    &
         &              anharmonic_data                        )) ]
    endif
    
    ! Increment the loop counter.
    i = i+1
  enddo
end function

impure elemental function calculate_split_states_SplitQpointsBasis(this, &
   & subspace,subspace_potential,anharmonic_data) result(output)
  implicit none
  
  class(SplitQpointsBasis), intent(in) :: this
  type(DegenerateSubspace), intent(in) :: subspace
  class(PotentialData),     intent(in) :: subspace_potential
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(SplitQpointsStates)             :: output
  
  ! Variables for constructing the basis.
  type(StructureData) :: supercell
  
  ! Variables for calculating and diagonalising the VSCF Hamiltonian.
  integer                                :: no_states
  type(HarmonicState)                    :: bra
  type(HarmonicState)                    :: ket
  real(dp),                  allocatable :: hamiltonian(:,:)
  type(SymmetricEigenstuff), allocatable :: estuff(:)
  
  type(SplitQpointsState)              :: vscf_state
  type(SplitQpointsState), allocatable :: vscf_states(:)
  
  integer :: i,j,k,l,ialloc
  
  supercell = anharmonic_data%anharmonic_supercell
  
  vscf_states = [SplitQpointsState::]
  do i=1,size(this)
    no_states = size(this%wavevectors(i))
    allocate( hamiltonian(no_states,no_states), &
            & stat=ialloc); call err(ialloc)
    hamiltonian = 0
    do j=1,no_states
      bra = this%wavevectors(i)%harmonic_states(j)
      do k=1,size(this%wavevectors(i)%harmonic_couplings(j))
        l = this%wavevectors(i)%harmonic_couplings(j)%id(k)
        ket = this%wavevectors(i)%harmonic_states(l)
        
        hamiltonian(j,l) = kinetic_energy( bra,                  &
                       &                   ket,                  &
                       &                   subspace,             &
                       &                   this,                 &
                       &                   anharmonic_data )     &
                       & + potential_energy( bra,                &
                       &                     subspace_potential, &
                       &                     ket,                &
                       &                     subspace,           &
                       &                     this,               &
                       &                     anharmonic_data )
      enddo
    enddo
    
    estuff = diagonalise_symmetric(hamiltonian)
    do j=1,size(estuff)
      vscf_state = SplitQpointsState(                     &
         & subspace_id  = this%subspace_id,               &
         & wavevector   = this%wavevectors(i)%wavevector, &
         & degeneracy   = this%wavevectors(i)%degeneracy, &
         & energy       = estuff(j)%eval,                 &
         & coefficients = estuff(j)%evec                  )
      
      ! Make energies extensive rather than intensive.
      vscf_state%energy = vscf_state%energy * supercell%sc_size
      
      vscf_states = [vscf_states, vscf_state]
    enddo
    deallocate(hamiltonian, stat=ialloc); call err(ialloc)
  enddo
  
  output = SplitQpointsStates(vscf_states)
end function

impure elemental function steps_converged(this,that,energy_convergence) &
   & result(output)
  implicit none
  
  type(SplitQpointsStates), intent(in) :: this
  type(SplitQpointsStates), intent(in) :: that
  real(dp),                 intent(in) :: energy_convergence
  logical                              :: output
  
  output = all( abs(this%vscf_states%energy-that%vscf_states%energy) &
            & < energy_convergence                                   )
end function

! ----------------------------------------------------------------------
! SplitQpointsState methods.
! ----------------------------------------------------------------------
! Constructor.
function new_SplitQpointsState(subspace_id,wavevector,degeneracy,energy, &
   & coefficients) result(this)
  implicit none
  
  integer,              intent(in) :: subspace_id
  type(FractionVector), intent(in) :: wavevector
  integer,              intent(in) :: degeneracy
  real(dp),             intent(in) :: energy
  real(dp),             intent(in) :: coefficients(:)
  type(SplitQpointsState)          :: this
  
  this%subspace_id  = subspace_id
  this%wavevector   = wavevector
  this%degeneracy   = degeneracy
  this%energy       = energy
  this%coefficients = coefficients
end function

recursive function new_SplitQpointsState_SubspaceState(input) result(this)
  implicit none
  
  class(SubspaceState), intent(in) :: input
  type(SplitQpointsState)          :: this
  
  select type(input); type is(SplitQpointsState)
    this = input
  type is(SubspaceStatePointer)
    this = SplitQpointsState(input%state())
  class default
    call err()
  end select
end function

! Type representation.
impure elemental function representation_SplitQpointsState() result(output)
  implicit none
  
  type(String) :: output
  
  output = 'split_qpoints'
end function

! Construct a PolynomialState from a SplitQpointsState.
! N.B. this state does not span the entire subspace,
!    only the modes at one q-point.
impure elemental function new_PolynomialState_SplitQpointsState(state,basis, &
   & qpoint) result(output)
  implicit none
  
  type(SplitQpointsState), intent(in)           :: state
  class(SubspaceBasis),    intent(in)           :: basis
  type(QpointData),        intent(in), optional :: qpoint
  type(PolynomialState)                         :: output
  
  type(SplitQpointsBasis) :: split_basis
  
  real(dp), allocatable :: coefficients(:)
  
  integer :: i
  
  split_basis = SplitQpointsBasis(basis)
  
  i = first(split_basis%wavevectors%wavevector==state%wavevector)
  
  coefficients = split_basis%wavevectors(i)%coefficients_basis_to_states( &
                                                     & state%coefficients )
  
  output = PolynomialState( split_basis%subspace_id,                    &
                          & split_basis%wavevectors(i)%monomial_states, &
                          & coefficients                                )
  
  if (present(qpoint)) then
    output = change_qpoint(output,SplitQpointsBasis(basis),qpoint)
  endif
end function

impure elemental function change_qpoint(state,basis,qpoint) result(output)
  implicit none
  
  type(PolynomialState),   intent(in) :: state
  type(SplitQpointsBasis), intent(in) :: basis
  type(QpointData),        intent(in) :: qpoint
  type(PolynomialState)               :: output
  
  type(QpointModes)    :: qpoint_from
  type(QpointModes)    :: qpoint_to
  integer              :: max_mode_id
  integer, allocatable :: mode_map(:)
  type(Group)          :: mode_group
  
  integer :: i
  
  ! Construct the group which maps the modes at the first q-point to the
  !    modes at the specified q-point.
  qpoint_from = basis%qpoint_modes(1)
  qpoint_to   = basis%qpoint_modes(                               &
                 & first(basis%qpoint_modes%qpoint%id==qpoint%id) )
  max_mode_id = maxval([qpoint_from%modes%id, qpoint_from%paired_modes%id])
  mode_map = [(0,i=1,max_mode_id)]
  do i=1,size(qpoint_to%modes)
    mode_map(qpoint_from%modes(i)%id) = qpoint_to%modes(i)%id
    mode_map(qpoint_from%paired_modes(i)%id) = qpoint_to%paired_modes(i)%id
  enddo
  mode_group = Group(mode_map)
  
  ! Use the group to change the state.
  output = state%change_modes(mode_group)
end function

impure elemental function wavefunction_SplitQpointsState(this,basis, &
   & supercell) result(output)
  implicit none
  
  class(SplitQpointsState), intent(in) :: this
  type(SplitQpointsBasis),  intent(in) :: basis
  type(StructureData),      intent(in) :: supercell
  type(String)                         :: output
  
  type(WavevectorBasis) :: wavevector_basis
  
  real(dp),     allocatable :: coefficients(:)
  type(String), allocatable :: terms(:)
  
  type(String) :: state
  
  integer :: i,ialloc
  
  wavevector_basis = basis%wavevectors(                     &
     & first(basis%wavevectors%wavevector==this%wavevector) )
  
  coefficients = wavevector_basis%coefficients_basis_to_states( &
                                            & this%coefficients )
  allocate(terms(size(coefficients)), stat=ialloc); call err(ialloc)
  do i=1,size(coefficients)
    state = wavevector_basis%monomial_states(i)%wavefunction( &
                                           & basis%frequency, &
                                           & supercell        )
    ! Trim the '|0>' from the state.
    state = slice(state,1,len(state)-3)
    
    terms(i) = state
  enddo
  output = '('//join(terms,' + ')//')|0>'
end function

impure elemental function braket_SubspaceState_SplitQpointsState(this, &
   & ket,subspace,subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SplitQpointsState), intent(in)           :: this
  class(SubspaceState),     intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  real(dp)                                       :: output
  
  type(PolynomialState) :: bra_state
  type(PolynomialState) :: ket_state
  
  bra_state = PolynomialState(this,subspace_basis)
  if (present(ket)) then
    ket_state = PolynomialState(SplitQpointsState(ket), subspace_basis)
    output = bra_state%braket( ket_state,      &
                             & subspace,       &
                             & subspace_basis, &
                             & anharmonic_data )
  else
    output = bra_state%braket( subspace        = subspace,       &
                             & subspace_basis  = subspace_basis, &
                             & anharmonic_data = anharmonic_data )
  endif
end function

impure elemental function braket_ComplexUnivariate_SplitQpointsState(    &
   & this,univariate,ket,subspace,subspace_basis,anharmonic_data,qpoint) &
   & result(output)
  implicit none
  
  class(SplitQpointsState), intent(in)           :: this
  type(ComplexUnivariate),  intent(in)           :: univariate
  class(SubspaceState),     intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(QpointData),         intent(in), optional :: qpoint
  type(ComplexMonomial)                          :: output
  
  type(PolynomialState) :: bra_state
  type(PolynomialState) :: ket_state
  
  bra_state = PolynomialState(this,subspace_basis,qpoint)
  if (present(ket)) then
    ket_state = PolynomialState(SplitQpointsState(ket),subspace_basis,qpoint)
    output = bra_state%braket( univariate,      &
                             & ket_state,       &
                             & subspace,        &
                             & subspace_basis,  &
                             & anharmonic_data, &
                             & qpoint           )
  else
    output = bra_state%braket( univariate      = univariate,      &
                             & subspace        = subspace,        &
                             & subspace_basis  = subspace_basis,  &
                             & anharmonic_data = anharmonic_data, &
                             & qpoint          = qpoint           )
  endif
end function

impure elemental function braket_ComplexMonomial_SplitQpointsState(this, &
   & monomial,ket,subspace,subspace_basis,anharmonic_data,qpoint)        &
   & result(output)
  implicit none
  
  class(SplitQpointsState), intent(in)           :: this
  type(ComplexMonomial),    intent(in)           :: monomial
  class(SubspaceState),     intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(QpointData),         intent(in), optional :: qpoint
  type(ComplexMonomial)                          :: output
  
  type(PolynomialState) :: bra_state
  type(PolynomialState) :: ket_state
  
  bra_state = PolynomialState(this,subspace_basis,qpoint)
  if (present(ket)) then
    ket_state = PolynomialState(SplitQpointsState(ket),subspace_basis,qpoint)
    output = bra_state%braket( monomial,        &
                             & ket_state,       &
                             & subspace,        &
                             & subspace_basis,  &
                             & anharmonic_data, &
                             & qpoint           )
  else
    output = bra_state%braket( monomial        = monomial,        &
                             & subspace        = subspace,        &
                             & subspace_basis  = subspace_basis,  &
                             & anharmonic_data = anharmonic_data, &
                             & qpoint          = qpoint           )
  endif
end function

impure elemental function kinetic_energy_SplitQpointsState(this,ket, &
   & subspace,subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SplitQpointsState), intent(in)           :: this
  class(SubspaceState),     intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  real(dp)                                       :: output
  
  type(PolynomialState) :: bra_state
  type(PolynomialState) :: ket_state
  
  bra_state = PolynomialState(this,subspace_basis)
  if (present(ket)) then
    ket_state = PolynomialState(SplitQpointsState(ket), subspace_basis)
    output = bra_state%kinetic_energy( ket_state,      &
                                     & subspace,       &
                                     & subspace_basis, &
                                     & anharmonic_data )
  else
    output = bra_state%kinetic_energy( subspace        = subspace,       &
                                     & subspace_basis  = subspace_basis, &
                                     & anharmonic_data = anharmonic_data )
  endif
end function

impure elemental function harmonic_potential_energy_SplitQpointsState( &
   & this,ket,subspace,subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SplitQpointsState), intent(in)           :: this
  class(SubspaceState),     intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  real(dp)                                       :: output
  
  type(PolynomialState) :: bra_state
  type(PolynomialState) :: ket_state
  
  bra_state = PolynomialState(this,subspace_basis)
  if (present(ket)) then
    ket_state = PolynomialState(SplitQpointsState(ket), subspace_basis)
    output = bra_state%harmonic_potential_energy( ket_state,      &
                                                & subspace,       &
                                                & subspace_basis, &
                                                & anharmonic_data )
  else
    output = bra_state%harmonic_potential_energy( &
              & subspace        = subspace,       &
              & subspace_basis  = subspace_basis, &
              & anharmonic_data = anharmonic_data )
  endif
end function

impure elemental function kinetic_stress_SplitQpointsState(this,ket, &
   & subspace,subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SplitQpointsState), intent(in)           :: this
  class(SubspaceState),     intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(RealMatrix)                               :: output
  
  type(PolynomialState) :: bra_state
  type(PolynomialState) :: ket_state
  
  bra_state = PolynomialState(this,subspace_basis)
  if (present(ket)) then
    ket_state = PolynomialState(SplitQpointsState(ket), subspace_basis)
    output = bra_state%kinetic_stress( ket_state,      &
                                     & subspace,       &
                                     & subspace_basis, &
                                     & anharmonic_data )
  else
    output = bra_state%kinetic_stress( subspace        = subspace,       &
                                     & subspace_basis  = subspace_basis, &
                                     & anharmonic_data = anharmonic_data )
  endif
end function

! ----------------------------------------------------------------------
! SplitQpointsStates methods.
! ----------------------------------------------------------------------
! Constructor.
function new_SplitQpointsStates(vscf_states) result(this)
  implicit none
  
  type(SplitQpointsState), intent(in) :: vscf_states(:)
  type(SplitQpointsStates)            :: this
  
  this%vscf_states = vscf_states
end function

recursive function new_SplitQpointsStates_SubspaceStates(input) result(this)
  implicit none
  
  class(SubspaceStates), intent(in) :: input
  type(SplitQpointsStates)          :: this
  
  select type(input); type is(SplitQpointsStates)
    this = input
  type is(SubspaceStatesPointer)
    this = SplitQpointsStates(input%states())
  class default
    call err()
  end select
end function

! Type representation.
impure elemental function representation_SplitQpointsStates() result(output)
  implicit none
  
  type(String) :: output
  
  output = 'split_qpoints'
end function

! Energy spectra.
impure elemental function spectra_SplitQpointsStates(this,subspace,     &
   & subspace_potential,subspace_stress,subspace_basis,anharmonic_data) &
   & result(output)
  implicit none
  
  class(SplitQpointsStates), intent(in)           :: this
  type(DegenerateSubspace),  intent(in)           :: subspace
  class(PotentialData),      intent(in)           :: subspace_potential
  class(StressData),         intent(in), optional :: subspace_stress
  class(SubspaceBasis),      intent(in)           :: subspace_basis
  type(AnharmonicData),      intent(in)           :: anharmonic_data
  type(EnergySpectra)                             :: output
  
  type(SplitQpointsBasis)       :: split_basis
  type(RealMatrix), allocatable :: stress(:)
  type(EnergySpectrum)          :: energy_spectrum
  
  integer :: i,ialloc
  
  split_basis = SplitQpointsBasis(subspace_basis)
  
  if (present(subspace_stress)) then
    allocate(stress(size(this%vscf_states)), stat=ialloc); call err(ialloc)
    do i=1,size(this%vscf_states)
      stress(i) = potential_stress( this%vscf_states(i), &
              &                     subspace_stress,     &
              &                     subspace,            &
              &                     subspace_basis,      &
              &                     anharmonic_data   )  &
              & + kinetic_stress( this%vscf_states(i),   &
              &                   subspace,              &
              &                   subspace_basis,        &
              &                   anharmonic_data      )
    enddo
    energy_spectrum =  EnergySpectrum( this%vscf_states%energy,     &
                                     & this%vscf_states%degeneracy, &
                                     & stress                       )
  else
    energy_spectrum = EnergySpectrum( this%vscf_states%energy,    &
                                    & this%vscf_states%degeneracy )
  endif
  
  output = EnergySpectra([( energy_spectrum,               &
                          & i=1,                           &
                          & size(split_basis%qpoint_modes) )])
end function

! Wavefunctions.
impure elemental function wavefunctions_SplitQpointsStates(this,subspace, &
   & subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SplitQpointsStates), intent(in) :: this
  type(DegenerateSubspace),  intent(in) :: subspace
  class(SubspaceBasis),      intent(in) :: subspace_basis
  type(AnharmonicData),      intent(in) :: anharmonic_data
  type(SubspaceWavefunctionsPointer)    :: output
  
  type(SplitQpointsBasis)         :: split_basis
  type(String)                    :: ground_state
  type(String), allocatable       :: state_wavefunctions(:)
  type(SplitQpointsWavefunctions) :: wavefunctions
  
  integer :: ialloc
  
  split_basis = SplitQpointsBasis(subspace_basis)
  
  ! Construct the wavefunction of |0>.
  ground_state = split_basis%ground_state_wavefunction( &
                 & subspace,                            &
                 & anharmonic_data%anharmonic_supercell )
  
  ! Construct the wavefunctions of each state, in terms of |0>.
  state_wavefunctions = this%vscf_states%wavefunction( &
                & split_basis,                         &
                & anharmonic_data%anharmonic_supercell )
  
  ! Pack the wavefunctions into the output,
  !    of concrete type SplitQpointsWavefunctions.
  wavefunctions = SplitQpointsWavefunctions( subspace%id,                 &
                                           & subspace%mode_ids,           &
                                           & subspace%paired_ids,         &
                                           & ground_state,                &
                                           & this%vscf_states%energy,     &
                                           & this%vscf_states%degeneracy, &
                                           & state_wavefunctions          )
  
  ! Convert the output to abstract class SubspaceWavefunctions.
  output = SubspaceWavefunctionsPointer(wavefunctions)
end function

! Integrate a potential.
impure elemental function integrate_potential_SplitQpointsStates(this, &
   & potential,subspace,subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SplitQpointsStates), intent(in) :: this
  class(PotentialData),      intent(in) :: potential
  type(DegenerateSubspace),  intent(in) :: subspace
  class(SubspaceBasis),      intent(in) :: subspace_basis
  type(AnharmonicData),      intent(in) :: anharmonic_data
  type(PotentialPointer)                :: output
  
  type(SplitQpointsState) :: ground_state
  type(SplitQpointsBasis) :: split_basis
  
  integer :: i
  
  ! Identify the ground state.
  ground_state = this%vscf_states(minloc(this%vscf_states%energy,1))
  
  ! Braket the potential between the ground state.
  split_basis = SplitQpointsBasis(subspace_basis)
  output = PotentialPointer(potential)
  do i=1,size(split_basis%qpoint_modes)
    output = braket( ground_state,                      &
                   & output,                            &
                   & subspace,                          &
                   & subspace_basis,                    &
                   & anharmonic_data,                   &
                   & split_basis%qpoint_modes(i)%qpoint )
  enddo
end function

! Integrate a stress.
impure elemental function integrate_stress_SplitQpointsStates(this,stress, &
   & subspace,subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SplitQpointsStates), intent(in) :: this
  class(StressData),         intent(in) :: stress
  type(DegenerateSubspace),  intent(in) :: subspace
  class(SubspaceBasis),      intent(in) :: subspace_basis
  type(AnharmonicData),      intent(in) :: anharmonic_data
  type(StressPointer)                   :: output
  
  type(SplitQpointsState) :: ground_state
  
  ! TODO: include q-points.
  ! MULTIPLE STATES
  
  ! Identify the ground state.
  ground_state = this%vscf_states(minloc(this%vscf_states%energy,1))
  
  ! Braket the stress between the ground state.
  output = braket( ground_state,   &
                 & stress,         &
                 & subspace,       &
                 & subspace_basis, &
                 & anharmonic_data )
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_SplitQpointsBasis(this,input)
  implicit none
  
  class(SplitQpointsBasis), intent(out) :: this
  type(String),             intent(in)  :: input(:)
  
  integer                            :: maximum_power
  integer                            :: expansion_order
  integer                            :: subspace_id
  real(dp)                           :: frequency
  type(WavevectorBasis), allocatable :: wavevectors(:)
  
  integer, allocatable :: starting_lines(:)
  integer, allocatable :: ending_lines(:)
  
  type(String), allocatable :: line(:)
  type(String), allocatable :: wavevector_lines(:)
  
  integer :: i,ialloc
  
  ! TODO
  call err()
  
  !select type(this); type is(SplitQpointsBasis)
  !  line = split_line(input(1))
  !  maximum_power = int(line(4))
  !  
  !  line = split_line(input(2))
  !  expansion_order = int(line(4))
  !  
  !  line = split_line(input(3))
  !  subspace_id = int(line(3))
  !  
  !  line = split_line(input(4))
  !  frequency = dble(line(3))
  !  
  !  starting_lines = [integer::]
  !  do i=5,size(input)
  !    line = split_line(input(i))
  !    if (size(line)>0) then
  !      if (line(1)=='Wavevector') then
  !        starting_lines = [starting_lines, i]
  !      endif
  !    endif
  !  enddo
  !  
  !  ending_lines = [starting_lines(2:)-1, size(input)]
  !  
  !  allocate(wavevectors(size(starting_lines)), stat=ialloc); call err(ialloc)
  !  do i=1,size(starting_lines)
  !    wavevector_lines = [input(:4), input(starting_lines(i):ending_lines(i))]
  !    wavevectors(i) = WavevectorBasis(wavevector_lines)
  !  enddo
  !  
  !  this = SplitQpointsBasis( maximum_power,   &
  !                          & expansion_order, &
  !                          & subspace_id,     &
  !                          & frequency,       &
  !                          & wavevectors      )
  !class default
  !  call err()
  !end select
end subroutine

function write_SplitQpointsBasis(this) result(output)
  implicit none
  
  class(SplitQpointsBasis), intent(in) :: this
  type(String), allocatable            :: output(:)
  
  type(String), allocatable :: wavevector_strings(:)
  
  integer :: i
  
  ! TODO
  call err()
  
  !select type(this); type is(SplitQpointsBasis)
  !  output = [ 'Maximum power   : '//this%maximum_power,   &
  !           & 'Expansion order : '//this%expansion_order, &
  !           & 'Subspace        : '//this%subspace_id,     &
  !           & 'Frequency       : '//this%frequency        ]
  !  do i=1,size(this%wavevectors)
  !    wavevector_strings = str(this%wavevectors(i))
  !    output = [ output, wavevector_strings(5:) ]
  !  enddo
  !class default
  !  call err()
  !end select
end function

function new_SplitQpointsBasis_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(SplitQpointsBasis)  :: this
  
  call this%read(input)
end function

impure elemental function new_SplitQpointsBasis_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(SplitQpointsBasis)       :: this
  
  this = SplitQpointsBasis(str(input))
end function

subroutine read_SplitQpointsState(this,input)
  implicit none
  
  class(SplitQpointsState), intent(out) :: this
  type(String),             intent(in)  :: input(:)
  
  integer               :: subspace_id
  type(FractionVector)  :: wavevector
  integer               :: degeneracy
  real(dp)              :: energy
  real(dp), allocatable :: coefficients(:)
  
  type(String), allocatable :: line(:)
  
  select type(this); type is(SplitQpointsState)
    line = split_line(input(1))
    subspace_id = int(line(3))
    
    line = split_line(input(2))
    wavevector = FractionVector(join(line(3:5)))
    
    line = split_line(input(3))
    degeneracy = int(line(3))
    
    line = split_line(input(4))
    energy = dble(line(3))
    
    line = split_line(input(5))
    coefficients = dble(line(3:))
    
    this = SplitQpointsState( subspace_id, &
                            & wavevector,  &
                            & degeneracy,  &
                            & energy,      &
                            & coefficients )
  class default
    call err()
  end select
end subroutine

function write_SplitQpointsState(this) result(output)
  implicit none
  
  class(SplitQpointsState), intent(in) :: this
  type(String), allocatable            :: output(:)
  
  select type(this); type is(SplitQpointsState)
    output = [ 'Subspace     : '//this%subspace_id, &
             & 'Wavevector   : '//this%wavevector,  &
             & 'Degeneracy   : '//this%degeneracy,  &
             & 'Energy       : '//this%energy,      &
             & 'Coefficients : '//this%coefficients ]
  class default
    call err()
  end select
end function

function new_SplitQpointsState_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(SplitQpointsState)  :: this
  
  call this%read(input)
end function

impure elemental function new_SplitQpointsState_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(SplitQpointsState)       :: this
  
  this = SplitQpointsState(str(input))
end function

subroutine read_SplitQpointsStates(this,input)
  implicit none
  
  class(SplitQpointsStates), intent(out) :: this
  type(String),              intent(in)  :: input(:)
  
  type(SplitQpointsState), allocatable :: vscf_states(:)
  
  select type(this); type is(SplitQpointsStates)
    vscf_states = SplitQpointsState(                           &
       & split_into_sections( input,                           &
       &                      separating_line=repeat('=',50) ) )
    this = SplitQpointsStates(vscf_states)
  class default
    call err()
  end select
end subroutine

function write_SplitQpointsStates(this) result(output)
  implicit none
  
  class(SplitQpointsStates), intent(in) :: this
  type(String), allocatable            :: output(:)
  
  select type(this); type is(SplitQpointsStates)
    output = str(this%vscf_states, separating_line=repeat('=',50))
  class default
    call err()
  end select
end function

function new_SplitQpointsStates_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(SplitQpointsStates) :: this
  
  call this%read(input)
end function

impure elemental function new_SplitQpointsStates_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(SplitQpointsStates)      :: this
  
  this = SplitQpointsStates(str(input))
end function
end module
