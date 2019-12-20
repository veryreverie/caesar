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
  use core_shell_thermodynamics_module
  use calculate_weights_module
  implicit none
  
  private
  
  public :: startup_split_qpoints_basis
  
  public :: SplitQpointsBasis
  
  public :: initial_ground_state
  
  type, extends(Stringable) :: QpointModeIDs
    integer, allocatable :: mode_ids(:)
    integer, allocatable :: paired_mode_ids(:)
  contains
    procedure, public :: read  => read_QpointModeIDs
    procedure, public :: write => write_QpointModeIDs
  end type
  
  interface QpointModeIDs
    module procedure new_QpointModeIDs
    module procedure new_QpointModeIDs_String
  end interface
  
  ! All states spanning the subspace.
  type, extends(SubspaceBasis) :: SplitQpointsBasis
    ! The maximum power of the monomial states.
    ! This is also the maximum occupation of the harmonic basis states.
    integer  :: maximum_power
    ! The expansion order of the potential.
    ! This is also the limit on coupling between states.
    integer :: expansion_order
    ! The ID of the subspace.
    integer  :: subspace_id
    ! The states, wavevector by wavevector.
    ! N.B. this only includes one wavevector from each symmetry-related set.
    type(WavevectorBasis), allocatable :: wavevectors(:)
    ! The q-points, and the modes at each q-point.
    type(QpointModeIDs), allocatable :: qpoints(:)
    
    ! The q-points which should be integrated.
    integer, allocatable, private :: qpoints_to_integrate_(:)
  contains
    procedure, public, nopass :: representation => &
                               & representation_SplitQpointsBasis
    
    ! Print the harmonic ground-state wavefunction of the basis.
    procedure, public :: ground_state_wavefunction
    
    ! Generate the first guess at states.
    procedure, public :: initial_states => initial_states_SplitQpointsBasis
    
    ! Generate the eigenstates of a single-subspace potential.
    procedure, public :: calculate_states => calculate_states_SplitQpointsBasis
    
    ! Process a subspace potential or stress.
    ! For a SplitQpointsBasis, this converts the single-subspace object
    !    into a single-q-point object.
    procedure, public :: process_subspace_potential => &
                       & process_subspace_potential_SplitQpointsBasis
    procedure, public :: process_subspace_stress => &
                       & process_subspace_stress_SplitQpointsBasis
    
    ! Return the modes spanned by this basis.
    procedure, public :: mode_ids => mode_ids_SplitQpointsBasis
    procedure, public :: paired_mode_ids => paired_mode_ids_SplitQpointsBasis
    
    ! Procedures involving individual states.
    procedure, public :: inner_product => &
                       & inner_product_SplitQpointsBasis
    procedure, public :: integrate_BasisState => &
                       & integrate_BasisState_SplitQpointsBasis
    procedure, public :: kinetic_energy => &
                       & kinetic_energy_SplitQpointsBasis
    procedure, public :: harmonic_potential_energy => &
                       & harmonic_potential_energy_SplitQpointsBasis
    procedure, public :: kinetic_stress => &
                       & kinetic_stress_SplitQpointsBasis
    procedure, public :: wavefunction => wavefunction_SplitQpointsBasis
    
    ! Procedures involving sets of states.
    procedure, public :: thermodynamic_data => &
                       & thermodynamic_data_SplitQpointsBasis
    procedure, public :: wavefunctions => wavefunctions_SplitQpointsBasis
    procedure, public :: integrate_BasisStates => &
                       & integrate_BasisStates_SplitQpointsBasis
    
    ! I/O.
    procedure, public :: read  => read_SplitQpointsBasis
    procedure, public :: write => write_SplitQpointsBasis
    
    ! Private helper functions.
    procedure, private :: integrate_potential
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

! QpointModeIDs methods.
function new_QpointModeIDs(mode_ids,paired_mode_ids) result(this)
  implicit none
  
  integer, intent(in) :: mode_ids(:)
  integer, intent(in) :: paired_mode_ids(:)
  type(QpointModeIDs) :: this
  
  this%mode_ids = mode_ids(sort(mode_ids))
  this%paired_mode_ids = paired_mode_ids(sort(paired_mode_ids))
end function

! ----------------------------------------------------------------------
! SplitQpointsBasis methods.
! ----------------------------------------------------------------------
! Constructors.
function new_SplitQpointsBasis(maximum_power,expansion_order,subspace_id, &
   & frequency,wavevectors,qpoints) result(this)
  implicit none
  
  integer,               intent(in) :: maximum_power
  integer,               intent(in) :: expansion_order
  integer,               intent(in) :: subspace_id
  real(dp),              intent(in) :: frequency
  type(WavevectorBasis), intent(in) :: wavevectors(:)
  type(QpointModeIDs),   intent(in) :: qpoints(:)
  type(SplitQpointsBasis)           :: this
  
  integer :: i
  
  this%maximum_power         = maximum_power
  this%expansion_order       = expansion_order
  this%subspace_id           = subspace_id
  this%frequency             = frequency
  this%wavevectors           = wavevectors
  this%qpoints               = qpoints
  this%qpoints_to_integrate_ = [(i, i=1, size(qpoints))]
end function

recursive function new_SplitQpointsBasis_SubspaceBasis(input) result(this)
  implicit none
  
  class(SubspaceBasis), intent(in) :: input
  type(SplitQpointsBasis)          :: this
  
  select type(input); type is(SplitQpointsBasis)
    this = input
  type is(SubspaceBasisPointer)
    ! WORKAROUND: ifort doesn't recognise the interface to this function
    !    from within this function, so the full name is used instead.
    !this = SplitQpointsBasis(input%basis())
    this = new_SplitQpointsBasis_SubspaceBasis(input%basis())
  class default
    call err()
  end select
end function

! Type representation.
impure elemental function representation_SplitQpointsBasis() result(output)
  implicit none
  
  type(String) :: output
  
  output = 'split q-points basis'
end function

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
  
  type(ComplexMode),   allocatable :: qpoint_modes(:)
  type(QpointModeIDs), allocatable :: qpoint_mode_ids(:)
  
  type(WavevectorBasis), allocatable :: wavevectors(:)
  
  integer :: i,j,ialloc
  
  ! List the modes and corresponding q-points in the subspace.
  ! The q-points are de-duplicated, and then q-points with id>paired_id are
  !    removed.
  subspace_modes = subspace%modes(modes)
  subspace_qpoints = subspace%qpoints(modes, qpoints)
  
  subspace_qpoints = subspace_qpoints(set(subspace_qpoints%id))
  subspace_qpoints = subspace_qpoints(                                &
     & filter(subspace_qpoints%id<=subspace_qpoints%paired_qpoint_id) )
  
  allocate( qpoint_mode_ids(size(subspace_qpoints)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(qpoint_mode_ids)
    qpoint_modes = subspace_modes(                                &
       & filter(subspace_modes%qpoint_id==subspace_qpoints(i)%id) )
    qpoint_mode_ids(i) = QpointModeIDs( qpoint_modes%id,       &
                                      & qpoint_modes%paired_id )
  enddo
  
  ! Generate a basis at the first q-point only.
  wavevectors = WavevectorBasis( subspace,                  &
                               & frequency,                 &
                               & modes,                     &
                               & qpoints,                   &
                               & maximum_power,             &
                               & potential_expansion_order, &
                               & symmetries,                &
                               & subspace_qpoints(1)        )
  
  output = SplitQpointsBasis( maximum_power,             &
                            & potential_expansion_order, &
                            & subspace%id,               &
                            & frequency,                 &
                            & wavevectors,               &
                            & qpoint_mode_ids            )
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
  
  integer :: i,ialloc
  
  ! Calculate the (geometric) average mass.
  mass = product(supercell%atoms%mass())
  mass = mass**(1.0_dp/size(supercell%atoms))
  
  ! Calculate the coefficient.
  coefficient = 1
  allocate(terms(0), stat=ialloc); call err(ialloc)
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
   & thermal_energy,anharmonic_data) result(output)
  implicit none
  
  class(SplitQpointsBasis), intent(in) :: this
  type(DegenerateSubspace), intent(in) :: subspace
  real(dp),                 intent(in) :: thermal_energy
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(BasisStatesPointer)             :: output
  
  type(WavevectorState) :: ground_state
  
  integer :: ialloc
  
  ! Generate the state |0>.
  ground_state = initial_ground_state(this)
  
  ! Generate the set of states {|0>}.
  output = BasisStatesPointer(WavevectorStates( &
                & subspace_id = subspace%id,    &
                & states      = [ground_state], &
                & energies    = [0.0_dp],       &
                & weights     = [1.0_dp]        ))
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
! Uses a VSCF scheme between q-points in the subspace.
! See vscf.90 for details.
impure elemental function calculate_states_SplitQpointsBasis(this,subspace, &
   & subspace_potential,thermal_energy,convergence_data,anharmonic_data)    &
   & result(output)
  implicit none
  
  class(SplitQpointsBasis), intent(in) :: this
  type(DegenerateSubspace), intent(in) :: subspace
  class(PotentialData),     intent(in) :: subspace_potential
  real(dp),                 intent(in) :: thermal_energy
  type(ConvergenceData),    intent(in) :: convergence_data
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(BasisStatesPointer)             :: output
  
  type(WavevectorStates)             :: wavevector_states
  type(WavevectorState), allocatable :: states(:)
  real(dp),              allocatable :: energies(:)
  real(dp),              allocatable :: weights(:)
  
  integer :: i,ialloc
  
  allocate( states(0),   &
          & energies(0), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(this%wavevectors)
    wavevector_states = WavevectorStates(        &
       & this%wavevectors(i)%calculate_states(   &
       &            subspace,                    &
       &            subspace_potential,          &
       &            thermal_energy,              &
       &            convergence_data,            &
       &            anharmonic_data            ) )
    states = [states, wavevector_states%states]
    energies = [energies, wavevector_states%energies]
  enddo
  
  weights = calculate_weights(energies, thermal_energy)
  
  output = BasisStatesPointer(WavevectorStates( subspace%id, &
                                              & states,      &
                                              & energies,    &
                                              & weights      ))
end function

! Integrates the potential over all but the first q-point in the subspace.
impure elemental function process_subspace_potential_SplitQpointsBasis(this, &
   & potential,states,subspace,thermal_energy,anharmonic_data) result(output)
  implicit none
  
  class(SplitQpointsBasis), intent(in)    :: this
  class(PotentialData),     intent(in)    :: potential
  class(BasisStates),       intent(inout) :: states
  type(DegenerateSubspace), intent(in)    :: subspace
  real(dp),                 intent(in)    :: thermal_energy
  type(AnharmonicData),     intent(in)    :: anharmonic_data
  type(PotentialPointer)                  :: output
  
  type(SplitQpointsBasis) :: basis
  
  type(PotentialPointer) :: integrated_potential
  real(dp)               :: correction_energy
  
  integer :: i
  
  ! Calculate (prod_{j/=i}<j|)V(prod_{j/=i}|j>).
  output = PotentialPointer(potential)
  
  if (size(this%qpoints)==1) then
    return
  endif
  
  basis = this
  
  basis%qpoints_to_integrate_ = [(i,i=2,size(this%qpoints))]
  
  call output%braket( states,                           &
                    & subspace        = subspace,       &
                    & thermal_energy  = thermal_energy, &
                    & subspace_basis  = basis,          &
                    & whole_subspace  = .false.,        &
                    & anharmonic_data = anharmonic_data )
  
  ! Calculate (prod_i<i|)V(prod_i|i>), and use this to correct for
  !    over-counting. See main vscf method for details.
  integrated_potential = output
  
  basis%qpoints_to_integrate_ = [1]
  
  call integrated_potential%braket( states,                           &
                                  & subspace        = subspace,       &
                                  & thermal_energy  = thermal_energy, &
                                  & subspace_basis  = basis,          &
                                  & anharmonic_data = anharmonic_data )
  
  correction_energy = integrated_potential%undisplaced_energy() &
                  & * (1.0_dp-size(this%qpoints))/size(this%qpoints)
  
  call output%add_constant(correction_energy)
end function

impure elemental function process_subspace_stress_SplitQpointsBasis(this, &
   & stress,states,subspace,thermal_energy,anharmonic_data) result(output)
  implicit none
  
  class(SplitQpointsBasis), intent(in)    :: this
  class(StressData),        intent(in)    :: stress
  class(BasisStates),       intent(inout) :: states
  type(DegenerateSubspace), intent(in)    :: subspace
  real(dp),                 intent(in)    :: thermal_energy
  type(AnharmonicData),     intent(in)    :: anharmonic_data
  type(StressPointer)                     :: output
  
  type(SplitQpointsBasis) :: basis
  
  type(StressPointer) :: integrated_stress
  type(RealMatrix)    :: correction_stress
  
  integer :: i
  
  ! Calculate (prod_{j/=i}<j|)V(prod_{j/=i}|j>).
  output = StressPointer(stress)
  
  if (size(this%qpoints)==1) then
    return
  endif
  
  basis = this
  
  basis%qpoints_to_integrate_ = [(i,i=2,size(this%qpoints))]
  
  call output%braket( states,                           &
                    & subspace        = subspace,       &
                    & thermal_energy  = thermal_energy, &
                    & subspace_basis  = basis,          &
                    & whole_subspace  = .false.,        &
                    & anharmonic_data = anharmonic_data )
  
  ! Calculate (prod_i<i|)V(prod_i|i>), and use this to correct for
  !    over-counting. See main vscf method for details.
  integrated_stress = output
  
  basis%qpoints_to_integrate_ = [1]
  
  call integrated_stress%braket( states,                           &
                               & subspace        = subspace,       &
                               & thermal_energy  = thermal_energy, &
                               & subspace_basis  = basis,          &
                               & anharmonic_data = anharmonic_data )
  
  correction_stress = integrated_stress%undisplaced_stress() &
                  & * (1.0_dp-size(this%qpoints))/size(this%qpoints)
  
  call output%add_constant(correction_stress)
end function

function integrate_potential(this,potential,states,thermal_energy,subspace, &
   & anharmonic_data) result(output)
  implicit none
  
  class(SplitQpointsBasis), intent(in)    :: this
  class(PotentialData),     intent(in)    :: potential
  type(WavevectorStates),   intent(inout) :: states
  real(dp),                 intent(in)    :: thermal_energy
  type(DegenerateSubspace), intent(in)    :: subspace
  type(AnharmonicData),     intent(in)    :: anharmonic_data
  type(PotentialPointer)                  :: output
  
  type(SplitQpointsBasis) :: basis
  
  type(PotentialPointer) :: integrated_potential
  real(dp)               :: correction_energy
  
  integer :: i
  
  basis = this
  
  ! Calculate (prod_{j/=i}<j|)V(prod_{j/=i}|j>).
  output = PotentialPointer(potential)
  
  if (size(this%qpoints)==1) then
    return
  endif
  
  basis%qpoints_to_integrate_ = [(i,i=2,size(this%qpoints))]
  
  call output%braket( states,                           &
                    & subspace        = subspace,       &
                    & thermal_energy  = thermal_energy, &
                    & subspace_basis  = basis,          &
                    & whole_subspace  = .false.,        &
                    & anharmonic_data = anharmonic_data )
  
  ! Calculate (prod_i<i|)V(prod_i|i>), and use this to correct for
  !    over-counting. See main vscf method for details.
  integrated_potential = output
  
  basis%qpoints_to_integrate_ = [1]
  
  call integrated_potential%braket( states,                           &
                                  & subspace        = subspace,       &
                                  & thermal_energy  = thermal_energy, &
                                  & subspace_basis  = basis,          &
                                  & anharmonic_data = anharmonic_data )
  
  correction_energy = integrated_potential%undisplaced_energy() &
                  & * (1.0_dp-size(this%qpoints))/size(this%qpoints)
  
  call output%add_constant(correction_energy)
end function

impure elemental function wavefunction_SplitQpointsBasis(this,state, &
   & supercell) result(output)
  implicit none
  
  class(SplitQpointsBasis), intent(in) :: this
  type(WavevectorState),    intent(in) :: state
  type(StructureData),      intent(in) :: supercell
  type(String)                         :: output
  
  ! TODO
  output = ''
end function

function mode_ids_SplitQpointsBasis(this,subspace,anharmonic_data) &
   & result(output)
  implicit none
  
  class(SplitQpointsBasis), intent(in) :: this
  type(DegenerateSubspace), intent(in) :: subspace
  type(AnharmonicData),     intent(in) :: anharmonic_data
  integer, allocatable                 :: output(:)
  
  type(QpointModeIDs), allocatable :: qpoints(:)
  
  integer :: i
  
  output = [( this%qpoints(this%qpoints_to_integrate_(i))%mode_ids, &
            & i=1,                                                  &
            & size(this%qpoints_to_integrate_)                      )]
end function

function paired_mode_ids_SplitQpointsBasis(this,subspace,anharmonic_data) &
   & result(output)
  implicit none
  
  class(SplitQpointsBasis), intent(in) :: this
  type(DegenerateSubspace), intent(in) :: subspace
  type(AnharmonicData),     intent(in) :: anharmonic_data
  integer, allocatable                 :: output(:)
  
  type(QpointModeIDs), allocatable :: qpoints(:)
  
  integer :: i
  
  output = [( this%qpoints(this%qpoints_to_integrate_(i))%paired_mode_ids, &
            & i=1,                                                         &
            & size(this%qpoints_to_integrate_)                             )]
end function

impure elemental function inner_product_SplitQpointsBasis(this,bra,ket, &
   & subspace,anharmonic_data) result(output)
  implicit none
  
  class(SplitQpointsBasis), intent(in)           :: this
  class(BasisState),        intent(in)           :: bra
  class(BasisState),        intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  real(dp)                                       :: output
  
  type(WavevectorState)              :: split_bra
  type(WavevectorState), allocatable :: split_ket
  
  integer :: i
  
  split_bra = WavevectorState(bra)
  if (present(ket)) then
    split_ket = WavevectorState(ket)
  endif
  
  i = first(this%wavevectors%wavevector == split_bra%wavevector)
  
  output = this%wavevectors(i)%inner_product( split_bra,       &
                                            & split_ket,       &
                                            & subspace,        &
                                            & anharmonic_data  )
  
  output = output * size(this%qpoints_to_integrate_)
end function

impure elemental function integrate_BasisState_SplitQpointsBasis(this, &
   & bra,monomial,ket,subspace,anharmonic_data) result(output)
  implicit none
  
  class(SplitQpointsBasis), intent(in)           :: this
  class(BasisState),        intent(in)           :: bra
  type(SparseMonomial),     intent(in)           :: monomial
  class(BasisState),        intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  complex(dp)                                    :: output
  
  type(WavevectorState)              :: split_bra
  type(WavevectorState), allocatable :: split_ket
  
  type(SparseMonomial), allocatable :: monomials(:)
  
  integer :: i,j
  
  ! Convert the bra (and the ket if present) to type WavevectorState.
  split_bra = WavevectorState(bra)
  if (present(ket)) then
    split_ket = WavevectorState(ket)
  endif
  
  ! Find the basis at the correct wavevector for the states.
  i = first(this%wavevectors%wavevector == split_bra%wavevector)
  
  ! Split the monomial by q-point.
  monomials = split_monomial(this,monomial)
  
  ! Integrate across all q-points to be integrated.
  output = product(this%wavevectors(i)%integrate( split_bra,       &
                                                & monomials,       &
                                                & split_ket,       &
                                                & subspace,        &
                                                & anharmonic_data  ))
end function

! Helper function. Takes a monomial, and splits into an array of monomials,
!    where the output monomials correspond to each of this%qpoints in turn.
! Transforms each monomial so it appears to be at the first q-point.
! Integrating this monomial at the first q-point is equivalent to
!    integrating the input monomial in the basis of states at the
!    selected q-point, but is much more efficient.
function split_monomial(this,monomial) result(output)
  implicit none
  
  class(SplitQpointsBasis), intent(in) :: this
  type(SparseMonomial),     intent(in) :: monomial
  type(SparseMonomial), allocatable    :: output(:)
  
  type(QpointModeIDs), allocatable :: qpoints(:)
  
  integer :: modes_per_qpoint
  
  integer :: i,ialloc
  
  qpoints = this%qpoints(this%qpoints_to_integrate_)
  
  modes_per_qpoint = size(monomial%modes)/size(qpoints)
  
  allocate(output(size(qpoints)), stat=ialloc); call err(ialloc)
  do i=1,size(qpoints)
    ! Separate off the modes at each q-point.
    output(i) = SparseMonomial(                           &
       & modes=monomial%modes( (i-1)*modes_per_qpoint+1   &
       &                     :  i   *modes_per_qpoint   ) )
    ! Transform each mode to be at the first q-point.
    output(i)%modes%id = this%qpoints(1)%mode_ids
    output(i)%modes%paired_id = this%qpoints(1)%paired_mode_ids
  enddo
end function

impure elemental function kinetic_energy_SplitQpointsBasis(this,bra,ket, &
   & subspace,anharmonic_data) result(output)
  implicit none
  
  class(SplitQpointsBasis), intent(in)           :: this
  class(BasisState),        intent(in)           :: bra
  class(BasisState),        intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  real(dp)                                       :: output
  
  type(WavevectorState)              :: split_bra
  type(WavevectorState), allocatable :: split_ket
  
  integer :: i
  
  split_bra = WavevectorState(bra)
  if (present(ket)) then
    split_ket = WavevectorState(ket)
  endif
  
  i = first(this%wavevectors%wavevector == split_bra%wavevector)
  
  output = this%wavevectors(i)%kinetic_energy( split_bra,      &
                                             & split_ket,      &
                                             & subspace,       &
                                             & anharmonic_data )
  
  output = output * size(this%qpoints_to_integrate_)
end function

impure elemental function harmonic_potential_energy_SplitQpointsBasis(this, &
   & bra,ket,subspace,anharmonic_data) result(output)
  implicit none
  
  class(SplitQpointsBasis), intent(in)           :: this
  class(BasisState),        intent(in)           :: bra
  class(BasisState),        intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  real(dp)                                       :: output
  
  type(WavevectorState)              :: split_bra
  type(WavevectorState), allocatable :: split_ket
  
  integer :: i
  
  split_bra = WavevectorState(bra)
  if (present(ket)) then
    split_ket = WavevectorState(ket)
  endif
  
  i = first(this%wavevectors%wavevector == split_bra%wavevector)
  
  output = this%wavevectors(i)%harmonic_potential_energy( split_bra,      &
                                                        & split_ket,      &
                                                        & subspace,       &
                                                        & anharmonic_data )
  
  output = output * size(this%qpoints_to_integrate_)
end function

impure elemental function kinetic_stress_SplitQpointsBasis(this,bra,ket, &
   & subspace,stress_prefactors,anharmonic_data) result(output)
  implicit none
  
  class(SplitQpointsBasis), intent(in)           :: this
  class(BasisState),        intent(in)           :: bra
  class(BasisState),        intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  type(StressPrefactors),   intent(in)           :: stress_prefactors
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(RealMatrix)                               :: output
  
  type(WavevectorState)              :: split_bra
  type(WavevectorState), allocatable :: split_ket
  
  integer :: i
  
  split_bra = WavevectorState(bra)
  if (present(ket)) then
    split_ket = WavevectorState(ket)
  endif
  
  i = first(this%wavevectors%wavevector == split_bra%wavevector)
  
  output = this%wavevectors(i)%kinetic_stress( split_bra,         &
                                             & split_ket,         &
                                             & subspace,          &
                                             & stress_prefactors, &
                                             & anharmonic_data    )
  
  output = output * size(this%qpoints_to_integrate_)
end function

! Thermodynamic data. Energy, entropy, free energy etc.
impure elemental function thermodynamic_data_SplitQpointsBasis(this,    &
   & thermal_energy,states,subspace,subspace_potential,subspace_stress, &
   & stress_prefactors,anharmonic_data) result(output)
  implicit none
  
  class(SplitQpointsBasis), intent(in)           :: this
  real(dp),                 intent(in)           :: thermal_energy
  class(BasisStates),       intent(in)           :: states
  type(DegenerateSubspace), intent(in)           :: subspace
  class(PotentialData),     intent(in)           :: subspace_potential
  class(StressData),        intent(in), optional :: subspace_stress
  type(StressPrefactors),   intent(in), optional :: stress_prefactors
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(ThermodynamicData)                        :: output
  
  type(RealMatrix) :: stress
  
  integer :: i
  
  ! Calculate the thermodynamic properties for one q-point in the system.
  output = core_shell_thermodynamics(     &
     & thermal_energy,                    &
     & this%frequency,                    &
     & size(subspace)/size(this%qpoints), &
     & subspace,                          &
     & this%wavevectors,                  &
     & states,                            &
     & subspace_potential,                &
     & subspace_stress,                   &
     & stress_prefactors,                 &
     & anharmonic_data                    )
  
  ! Multiply by the number of q-points, and add in the constant term.
  output = output * size(this%qpoints)
  
  ! Average the stress across all symmetries.
  ! Each subspace obeys all symmetries, so the stress must too.
  ! This is needed because the properties are only calculated at one q-point,
  !    which will not in general obey the symmetries.
  if (allocated(output%stress)) then
    associate(symmetries=>anharmonic_data%structure%symmetries)
      stress = dblemat(zeroes(3,3))
      do i=1,size(symmetries)
        stress = stress                         &
             & + symmetries(i)%cartesian_tensor &
             & * output%stress                  &
             & * transpose(symmetries(i)%cartesian_tensor)
      enddo
      output%stress = stress / size(symmetries)
    end associate
  endif
end function

! Wavefunctions.
impure elemental function wavefunctions_SplitQpointsBasis(this,states, &
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
impure elemental function integrate_BasisStates_SplitQpointsBasis(this, &
   & states,thermal_energy,monomial,subspace,anharmonic_data) result(output)
  implicit none
  
  class(SplitQpointsBasis),  intent(in) :: this
  class(BasisStates),        intent(in) :: states
  real(dp),                  intent(in) :: thermal_energy
  type(SparseMonomial),      intent(in) :: monomial
  type(DegenerateSubspace),  intent(in) :: subspace
  type(AnharmonicData),      intent(in) :: anharmonic_data
  complex(dp)                           :: output
  
  type(WavevectorStates)  :: split_states
  
  type(SparseMonomial), allocatable :: monomials(:)
  
  integer :: i
  
  ! Convert the states to type WavevectorStates.
  split_states = WavevectorStates(states)
  
  ! Split the monomial by q-point.
  monomials = split_monomial(this,monomial)
  
  ! Integrate across all q-points to be integrated.
  output = 1.0_dp
  do i=1,size(this%qpoints_to_integrate_)
    output = output                       &
         & * integrate( this%wavevectors, &
         &              split_states,     &
         &              thermal_energy,   &
         &              monomials(i),     &
         &              anharmonic_data   )
  enddo
    
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_QpointModeIDs(this,input)
  implicit none
  
  class(QpointModeIDs), intent(out) :: this
  type(String),         intent(in)  :: input
  
  integer, allocatable :: mode_ids(:)
  integer, allocatable :: paired_mode_ids(:)
  
  type(String), allocatable :: line(:)
  
  select type(this); type is(QpointModeIDs)
    line = tokens(input, delimiter=',')
    
    mode_ids = int(tokens(line(1), first=2))
    paired_mode_ids = int(tokens(line(2), first=3))
    
    this = QpointModeIDs(mode_ids, paired_mode_ids)
  class default
    call err()
  end select
end subroutine

function write_QpointModeIDs(this) result(output)
  implicit none
  
  class(QpointModeIDs), intent(in) :: this
  type(String)                         :: output
  
  select type(this); type is(QpointModeIDs)
    output = 'Modes: '//this%mode_ids// &
           & ', Paired modes: '//this%paired_mode_ids
  class default
    call err()
  end select
end function

impure elemental function new_QpointModeIDs_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(QpointModeIDs)      :: this
  
  call this%read(input)
end function


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
  
  select type(this); type is(SplitQpointsBasis)
    output = [ 'Maximum power    : '//this%maximum_power,   &
             & 'Expansion order  : '//this%expansion_order, &
             & 'Subspace         : '//this%subspace_id,     &
             & 'Frequency        : '//this%frequency,       &
             & str('Wavevectors:'),                         &
             & str(repeat('-',50))                          ]
    do i=1,size(this%wavevectors)
      wavevector_strings = str(this%wavevectors(i))
      output = [ output, wavevector_strings(5:), str(repeat('-',50)) ]
    enddo
    output = [output, str('q-point mode ids:')]
    do i=1,size(this%qpoints)
      output = [output, str(this%qpoints(i))]
    enddo
  class default
    call err()
  end select
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
