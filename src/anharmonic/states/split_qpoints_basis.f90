! ======================================================================
! A basis of states which treats q-points separately.
! ======================================================================
! N.B. only states at a single q-point are calculated and stored.
! All other q-points are included by symmetry.
module split_qpoints_basis_module
  use common_module
  
  use anharmonic_common_module
  
  use wavevector_state_module
  use wavevector_states_module
  use wavevector_basis_module
  use split_qpoints_wavefunctions_module
  
  implicit none
  
  private
  
  public :: startup_split_qpoints_basis
  
  public :: SplitQpointsBasis
  
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
    ! The q-point which is currently being integrated.
    integer :: qpoint_mode
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
    
    ! Return the modes spanned by this basis.
    procedure, public :: modes => modes_SplitQpointsBasis
    
    ! Procedures involving individual states.
    procedure, public :: inner_product => &
                       & inner_product_WavevectorState
    procedure, public :: braket_ComplexMonomial => &
                       & braket_ComplexMonomial_WavevectorState
    procedure, public :: kinetic_energy => &
                       & kinetic_energy_WavevectorState
    procedure, public :: harmonic_potential_energy => &
                       & harmonic_potential_energy_WavevectorState
    procedure, public :: kinetic_stress => &
                       & kinetic_stress_WavevectorState
    procedure, public :: wavefunction => wavefunction_WavevectorState
    
    ! Procedures involving sets of states.
    procedure, public :: spectra => spectra_WavevectorStates
    procedure, public :: wavefunctions => wavefunctions_WavevectorStates
    procedure, public :: integrate_ComplexMonomial => &
                       & integrate_ComplexMonomial_WavevectorStates
    
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
contains

! Startup procedure.
subroutine startup_split_qpoints_basis()
  implicit none
  
  type(SplitQpointsBasis) :: basis
  
  call basis%startup()
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
! Constructors.
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
  this%qpoint_mode     = 0
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

function pick_qpoint(this,qpoint_mode) result(output)
  implicit none
  
  type(SplitQpointsBasis), intent(in) :: this
  integer,                 intent(in) :: qpoint_mode
  type(SplitQpointsBasis)             :: output
  
  output = this
  output%qpoint_mode = qpoint_mode
end function

! Type representation.
impure elemental function representation_SplitQpointsBasis() result(output)
  implicit none
  
  type(String) :: output
  
  output = 'split q-points basis'
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
  type(BasisStatesPointer)             :: output
  
  type(WavevectorState) :: ground_state
  
  integer :: ialloc
  
  ! Generate the state |0>.
  ground_state = initial_ground_state(this)
  
  ! Generate the set of states {|0>}.
  output = BasisStatesPointer(WavevectorStates([ground_state],[0.0_dp]))
end function

! Generate initial guess. This is simply the basis state |0>, i.e. the
!    ground state of the effective harmonic potential from which the basis
!    states were generated.
impure elemental function initial_ground_state(basis) result(output)
  implicit none
  
  type(SplitQpointsBasis), intent(in) :: basis
  type(WavevectorState)               :: output
  
  type(WavevectorBasis) :: wavevector_basis
  
  ! Find the basis at wavevector [0,0,0].
  wavevector_basis = basis%wavevectors(            &
     & first(is_int(basis%wavevectors%wavevector)) )
  
  ! Construct the ground state at [0,0,0].
  output = wavevector_basis%initial_ground_state()
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
  type(BasisStatesPointer)             :: output
  
  output = BasisStatesPointer(split_vscf( subspace_potential,        &
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
  type(WavevectorStates)               :: output
  
  type(WavevectorStates) :: initial_states
  type(WavevectorState)  :: initial_state
  type(PotentialPointer) :: initial_potential
  
  type(PotentialPointer), allocatable :: input_potentials(:)
  type(WavevectorStates), allocatable :: states(:)
  type(WavevectorState)               :: state
  type(PotentialPointer), allocatable :: output_potentials(:)
  
  integer :: first_pulay_step
  
  integer :: i,j
  
  call print_line( 'Running inter-subspace VSCF in subspace '// &
                 & subspace%id//'.')
  initial_states = WavevectorStates(basis%initial_states( subspace,       &
                                                        & anharmonic_data ))
  initial_state = initial_states%states(minloc(initial_states%energies,1) )
  initial_potential = PotentialPointer(potential)
  do i=2,size(basis%qpoint_modes)
    call initial_potential%braket( initial_state,                          &
                                 & subspace        = subspace,             &
                                 & subspace_basis  = pick_qpoint(basis,i), &
                                 & anharmonic_data = anharmonic_data       )
  enddo
  call initial_potential%zero_energy()
  
  input_potentials = [initial_potential]
  states = [WavevectorStates::]
  output_potentials = [PotentialPointer::]
  i = 1
  do
    ! Calculate new states.
    states = [ states,                                              &
             & basis%calculate_split_states( subspace,              &
             &                               input_potentials(i),   &
             &                               anharmonic_data      ) ]
    state = states(i)%states(minloc(states(i)%energies,1))
    
    ! Use states to calculate new potential.
    output_potentials = [output_potentials, PotentialPointer(potential)]
    do j=2,size(basis%qpoint_modes)
      call output_potentials(i)%braket(            &
         & state,                                  &
         & subspace        = subspace,             &
         & subspace_basis  = pick_qpoint(basis,j), &
         & anharmonic_data = anharmonic_data       )
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
  type(WavevectorStates)               :: output
  
  type(WavevectorState)              :: state
  type(WavevectorState), allocatable :: states(:)
  real(dp),              allocatable :: energies(:)
  
  type(WavevectorStates) :: wavevector_states
  
  integer :: i,j
  
  states = [WavevectorState::]
  energies = [real(dp)::]
  do i=1,size(this%wavevectors)
    wavevector_states = this%wavevectors(i)%calculate_states( &
                                        & subspace_potential, &
                                        & anharmonic_data     )
    states = [states, wavevector_states%states]
    energies = [energies, wavevector_states%energies]
  enddo
  
  output = WavevectorStates(states, energies)
end function

impure elemental function steps_converged(this,that,energy_convergence) &
   & result(output)
  implicit none
  
  type(WavevectorStates), intent(in) :: this
  type(WavevectorStates), intent(in) :: that
  real(dp),               intent(in) :: energy_convergence
  logical                            :: output
  
  output = all(abs(this%energies-that%energies) < energy_convergence)
end function

impure elemental function wavefunction_WavevectorState(this,state, &
   & supercell) result(output)
  implicit none
  
  class(SplitQpointsBasis), intent(in) :: this
  type(WavevectorState),    intent(in) :: state
  type(StructureData),      intent(in) :: supercell
  type(String)                         :: output
  
  ! TODO
  output = ''
end function

function modes_SplitQpointsBasis(this,subspace,anharmonic_data) result(output)
  implicit none
  
  class(SplitQpointsBasis), intent(in) :: this
  type(DegenerateSubspace), intent(in) :: subspace
  type(AnharmonicData),     intent(in) :: anharmonic_data
  integer, allocatable                 :: output(:)
  
  if (this%qpoint_mode==0) then
    output = subspace%mode_ids
  else
    if (    this%qpoint_modes(this%qpoint_mode)%modes(1)%id        &
       & == this%qpoint_modes(this%qpoint_mode)%modes(1)%paired_id ) then
      output = this%qpoint_modes(this%qpoint_mode)%modes%id
    else
      output = [ this%qpoint_modes(this%qpoint_mode)%modes%id,       &
               & this%qpoint_modes(this%qpoint_mode)%modes%paired_id ]
    endif
  endif
end function

impure elemental function inner_product_WavevectorState(this,bra,ket, &
   & subspace,anharmonic_data) result(output)
  implicit none
  
  class(SplitQpointsBasis), intent(in)           :: this
  class(BasisState),        intent(in)           :: bra
  class(BasisState),        intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  real(dp)                                       :: output
  
  type(WavevectorState) :: split_bra
  type(WavevectorState) :: split_ket
  
  integer :: i
  
  split_bra = WavevectorState(bra)
  
  i = first(this%wavevectors%wavevector == split_bra%wavevector)
  
  if (present(ket)) then
    split_ket = WavevectorState(ket)
    output = this%wavevectors(i)%inner_product( split_bra,       &
                                              & split_ket,       &
                                              & anharmonic_data  )
  else
    output = this%wavevectors(i)%inner_product( &
            & bra             = split_bra,      &
            & anharmonic_data = anharmonic_data )
  endif
end function

impure elemental function braket_ComplexMonomial_WavevectorState(this, &
   & bra,monomial,ket,subspace,anharmonic_data) result(output)
  implicit none
  
  class(SplitQpointsBasis), intent(in)           :: this
  class(BasisState),        intent(in)           :: bra
  type(ComplexMonomial),    intent(in)           :: monomial
  class(BasisState),        intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(ComplexMonomial)                          :: output
  
  type(ComplexUnivariate), allocatable :: qpoint_modes(:)
  type(ComplexUnivariate), allocatable :: non_qpoint_modes(:)
  type(ComplexUnivariate), allocatable :: transformed_qpoint_modes(:)
  type(ComplexMonomial)                :: qpoint_monomial
  
  type(WavevectorState) :: split_bra
  type(WavevectorState) :: split_ket
  
  integer :: i
  
  ! Check that qpoint_mode is set.
  if (this%qpoint_mode==0) then
    call print_line(CODE_ERROR//': split q-point basis must have a specified &
       &q-point in order to calculate <|X|>.')
    call err()
  endif
  
  ! Separate the modes at the given q-point from the other modes in the
  !    monomial.
  qpoint_modes = monomial%modes(                                        &
     & ids        = this%qpoint_modes(this%qpoint_mode)%modes%id,       &
     & paired_ids = this%qpoint_modes(this%qpoint_mode)%modes%paired_id )
  non_qpoint_modes = monomial%modes(                              &
     & exclude_ids = this%qpoint_modes(this%qpoint_mode)%modes%id )
  
  ! Construct a monomial with the terms of qpoint_modes, but transformed to
  !    the q-point at which the states are stored.
  ! Integrating this monomial at the first q-point is equivalent to
  !    integrating the input monomial in the basis of states at the
  !    selected q-point, but is much more efficient.
  transformed_qpoint_modes = ComplexUnivariate(   &
     & mode         = this%qpoint_modes(1)%modes, &
     & power        = qpoint_modes%power,         &
     & paired_power = qpoint_modes%paired_power   )
  qpoint_monomial = ComplexMonomial( coefficient = cmplx(1.0_dp,0.0_dp,dp), &
                                   & modes       = transformed_qpoint_modes )
  
  ! Integrate the separated modes.
  split_bra = WavevectorState(bra)
  
  i = first(this%wavevectors%wavevector == split_bra%wavevector)
  
  if (present(ket)) then
    split_ket = WavevectorState(ket)
    qpoint_monomial = this%wavevectors(i)%braket( split_bra,       &
                                                & qpoint_monomial, &
                                                & split_ket,       &
                                                & anharmonic_data  )
  else
    qpoint_monomial = this%wavevectors(i)%braket( &
             & bra             = split_bra,       &
             & monomial        = qpoint_monomial, &
             & anharmonic_data = anharmonic_data  )
  endif
  
  ! Reconstruct the output.
  output = ComplexMonomial(                                            &
     & coefficient = monomial%coefficient*qpoint_monomial%coefficient, &
     & modes       = non_qpoint_modes                                  )
end function

impure elemental function kinetic_energy_WavevectorState(this,bra,ket, &
   & subspace,anharmonic_data) result(output)
  implicit none
  
  class(SplitQpointsBasis), intent(in)           :: this
  class(BasisState),        intent(in)           :: bra
  class(BasisState),        intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  real(dp)                                       :: output
  
  type(WavevectorState) :: split_bra
  type(WavevectorState) :: split_ket
  
  integer :: i
  
  split_bra = WavevectorState(bra)
  
  i = first(this%wavevectors%wavevector == split_bra%wavevector)
  
  if (present(ket)) then
    split_ket = WavevectorState(ket)
    output = this%wavevectors(i)%kinetic_energy( split_bra,      &
                                               & split_ket,      &
                                               & anharmonic_data )
  else
    output = this%wavevectors(i)%kinetic_energy( &
             & bra             = split_bra,      &
             & anharmonic_data = anharmonic_data )
  endif
end function

impure elemental function harmonic_potential_energy_WavevectorState(this, &
   & bra,ket,subspace,anharmonic_data) result(output)
  implicit none
  
  class(SplitQpointsBasis), intent(in)           :: this
  class(BasisState),        intent(in)           :: bra
  class(BasisState),        intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  real(dp)                                       :: output
  
  type(WavevectorState) :: split_bra
  type(WavevectorState) :: split_ket
  
  integer :: i
  
  split_bra = WavevectorState(bra)
  
  i = first(this%wavevectors%wavevector == split_bra%wavevector)
  
  if (present(ket)) then
    split_ket = WavevectorState(ket)
    output = this%wavevectors(i)%harmonic_potential_energy( split_bra,      &
                                                          & split_ket,      &
                                                          & anharmonic_data )
  else
    output = this%wavevectors(i)%harmonic_potential_energy( &
                        & bra             = split_bra,      &
                        & anharmonic_data = anharmonic_data )
  endif
end function

impure elemental function kinetic_stress_WavevectorState(this,bra,ket, &
   & subspace,stress_prefactors,anharmonic_data) result(output)
  implicit none
  
  class(SplitQpointsBasis), intent(in)           :: this
  class(BasisState),        intent(in)           :: bra
  class(BasisState),        intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  type(StressPrefactors),   intent(in)           :: stress_prefactors
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(RealMatrix)                               :: output
  
  type(WavevectorState) :: split_bra
  type(WavevectorState) :: split_ket
  
  integer :: i
  
  split_bra = WavevectorState(bra)
  
  i = first(this%wavevectors%wavevector == split_bra%wavevector)
  
  if (present(ket)) then
    split_ket = WavevectorState(ket)
    output = this%wavevectors(i)%kinetic_stress( split_bra,         &
                                               & split_ket,         &
                                               & stress_prefactors, &
                                               & anharmonic_data    )
  else
    output = this%wavevectors(i)%kinetic_stress( &
        & bra               = split_bra,         &
        & stress_prefactors = stress_prefactors, &
        & anharmonic_data   = anharmonic_data    )
  endif
end function

! Energy spectra.
impure elemental function spectra_WavevectorStates(this,states,subspace,   &
   & subspace_potential,subspace_stress,stress_prefactors,anharmonic_data) &
   & result(output)
  implicit none
  
  class(SplitQpointsBasis), intent(in)           :: this
  class(BasisStates),       intent(in)           :: states
  type(DegenerateSubspace), intent(in)           :: subspace
  class(PotentialData),     intent(in)           :: subspace_potential
  class(StressData),        intent(in), optional :: subspace_stress
  type(StressPrefactors),   intent(in), optional :: stress_prefactors
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(EnergySpectra)                            :: output
  
  type(WavevectorStates)        :: split_states
  type(RealMatrix), allocatable :: stress(:)
  type(EnergySpectrum)          :: energy_spectrum
  
  integer :: i,ialloc
  
  if (present(subspace_stress) .neqv. present(stress_prefactors)) then
    call err()
  endif
  
  split_states = WavevectorStates(states)
  
  if (present(subspace_stress)) then
    allocate( stress(size(split_states%states)), &
            & stat=ialloc); call err(ialloc)
    do i=1,size(split_states%states)
      stress(i) = potential_stress( split_states%states(i),        &
              &                     subspace_stress,               &
              &                     subspace,                      &
              &                     this,                          &
              &                     anharmonic_data              ) &
              & + this%kinetic_stress(                             &
              &        bra               = split_states%states(i), &
              &        subspace          = subspace,               &
              &        stress_prefactors = stress_prefactors,      &
              &        anharmonic_data   = anharmonic_data         )
    enddo
    energy_spectrum =  EnergySpectrum(split_states%energies, stresses=stress)
  else
    energy_spectrum = EnergySpectrum(split_states%energies)
  endif
  
  output = EnergySpectra([( energy_spectrum,        &
                          & i=1,                    &
                          & size(this%qpoint_modes) )])
end function

! Wavefunctions.
impure elemental function wavefunctions_WavevectorStates(this,states, &
   & subspace,anharmonic_data) result(output)
  implicit none
  
  class(SplitQpointsBasis),  intent(in) :: this
  class(BasisStates),        intent(in) :: states
  type(DegenerateSubspace),  intent(in) :: subspace
  type(AnharmonicData),      intent(in) :: anharmonic_data
  type(SubspaceWavefunctionsPointer)    :: output
  
  type(WavevectorStates)          :: split_states
  type(String)                    :: ground_state
  type(String), allocatable       :: state_wavefunctions(:)
  type(SplitQpointsWavefunctions) :: wavefunctions
  
  integer :: ialloc
  
  split_states = WavevectorStates(states)
  
  ! Construct the wavefunction of |0>.
  ground_state = this%ground_state_wavefunction( &
          & subspace,                            &
          & anharmonic_data%anharmonic_supercell )
  
  ! Construct the wavefunctions of each state, in terms of |0>.
  state_wavefunctions = this%wavefunction(  &
     & split_states%states,                 &
     & anharmonic_data%anharmonic_supercell )
  
  ! Pack the wavefunctions into the output,
  !    of concrete type SplitQpointsWavefunctions.
  wavefunctions = SplitQpointsWavefunctions( subspace%id,           &
                                           & subspace%mode_ids,     &
                                           & subspace%paired_ids,   &
                                           & ground_state,          &
                                           & split_states%energies, &
                                           & state_wavefunctions    )
  
  ! Convert the output to abstract class SubspaceWavefunctions.
  output = SubspaceWavefunctionsPointer(wavefunctions)
end function

! Integrate a monomial.
impure elemental function integrate_ComplexMonomial_WavevectorStates(this, &
   & states,monomial,subspace,anharmonic_data) result(output)
  implicit none
  
  class(SplitQpointsBasis),  intent(in) :: this
  class(BasisStates),        intent(in) :: states
  type(ComplexMonomial),     intent(in) :: monomial
  type(DegenerateSubspace),  intent(in) :: subspace
  type(AnharmonicData),      intent(in) :: anharmonic_data
  type(ComplexMonomial)                 :: output
  
  type(WavevectorStates)  :: split_states
  type(WavevectorState)   :: ground_state
  type(SplitQpointsBasis) :: qpoint_basis
  
  integer :: i
  
  split_states = WavevectorStates(states)
  
  ! Identify the ground state.
  ground_state = split_states%states(minloc(split_states%energies,1))
  
  ! Braket the potential between the ground state.
  do i=1,size(this%qpoint_modes)
    qpoint_basis = pick_qpoint(this,i)
    output = qpoint_basis%braket( bra             = ground_state,   &
                                & monomial        = monomial,       &
                                & subspace        = subspace,       &
                                & anharmonic_data = anharmonic_data )
  enddo
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
end module
