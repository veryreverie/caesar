! ======================================================================
! A basis harmonic states for a given subspace at a given wavevector.
! N.B. the wavevector in this context is the wavevector of the state,
!    not the wavevector of the mode.
! e.g. the state |p> along a mode at q-point q has a wavevector of p*q.
! ======================================================================
! N.B. The WavevectorBasis has a polymorphic allocatable component,
!   harmonic_states_. Since compiler support for such things is not fantastic,
!   care should be taken when modifying the code referring to this component.
module caesar_wavevector_basis_module
  use caesar_common_module
  
  use caesar_anharmonic_common_module
  
  use caesar_harmonic_state_1d_module
  use caesar_harmonic_state_2d_module
  use caesar_harmonic_state_real_module
  use caesar_harmonic_state_complex_module
  use caesar_harmonic_braket_real_module
  use caesar_harmonic_braket_complex_module
  use caesar_coupled_states_module
  use caesar_density_matrix_module
  use caesar_wavevector_state_module
  use caesar_wavevector_states_module
  use caesar_calculate_weights_module
  implicit none
  
  private
  
  public :: startup_wavevector_basis
  
  public :: WavevectorBasis
  
  public :: calculate_states
  public :: core_harmonic_observables
  public :: core_effective_harmonic_observables
  public :: core_vci_observables
  public :: free_energy_gradient
  
  public :: HarmonicBraKetReal
  public :: HarmonicBraKetComplex
  
  type, extends(SubspaceBasis) :: WavevectorBasis
    integer                                    :: maximum_power
    integer                                    :: expansion_order
    integer                                    :: subspace_id
    type(FractionVector)                       :: wavevector
    class(SubspaceState), allocatable, private :: harmonic_states_(:)
    type(CoupledStates),  allocatable, private :: harmonic_couplings_(:)
  contains
    procedure, public, nopass :: representation => &
                               & representation_WavevectorBasis
    
    ! State procedures.
    procedure, public :: inner_product => &
                       & inner_product_WavevectorBasis
    procedure, public :: integrate_BasisState => &
                       & integrate_BasisState_WavevectorBasis
    procedure, public :: kinetic_energy => &
                       & kinetic_energy_WavevectorBasis
    procedure, public :: harmonic_potential_energy => &
                       & harmonic_potential_energy_WavevectorBasis
    procedure, public :: kinetic_stress => &
                       & kinetic_stress_WavevectorBasis
    
    ! Basis procedures.
    procedure, public :: initial_ground_state => &
                       & initial_ground_state_WavevectorBasis
    procedure, public :: initial_states => &
                       & initial_states_WavevectorBasis
    procedure, public :: calculate_states => &
                       & calculate_states_WavevectorBasis
    
    ! Mode IDs.
    procedure, public :: mode_ids => mode_ids_WavevectorBasis
    procedure, public :: paired_mode_ids => paired_mode_ids_WavevectorBasis
    
    ! I/O.
    procedure, public :: read  => read_WavevectorBasis
    procedure, public :: write => write_WavevectorBasis
    
    ! Inherited procedures which do not apply to this type.
    procedure, public :: integrate_BasisStates => &
                       & integrate_BasisStates_WavevectorBasis
    procedure, public :: thermodynamic_data => &
                       & thermodynamic_data_WavevectorBasis
    procedure, public :: wavefunctions => &
                       & wavefunctions_WavevectorBasis
    procedure, public :: free_energy_gradient => &
                       & free_energy_gradient_WavevectorBasis
  end type
  
  interface WavevectorBasis
    module procedure new_WavevectorBasis
    module procedure new_WavevectorBasis_subspace
    module procedure new_WavevectorBasis_Strings
    module procedure new_WavevectorBasis_StringArray
  end interface
  
  interface core_harmonic_observables
    module procedure core_harmonic_observables_WavevectorBasis
  end interface
  
  interface core_effective_harmonic_observables
    module procedure core_effective_harmonic_observables_WavevectorBasis
  end interface
  
  interface core_vci_observables
    module procedure core_vci_observables_WavevectorBasis
  end interface
  
  interface free_energy_gradient
    module procedure free_energy_gradient_WavevectorBases
  end interface
  
  interface HarmonicBraKetReal
    module procedure new_HarmonicBraKetReal_WavevectorBasis
  end interface
  
  interface HarmonicBraKetComplex
    module procedure new_HarmonicBraKetComplex_WavevectorBasis
  end interface
  
  type, extends(NoDefaultConstructor) :: SelectedStatesHamiltonian
    integer, allocatable   :: selected_states(:)
    type(SparseRealMatrix) :: hamiltonian
  end type
contains

! Startup procedure.
subroutine startup_wavevector_basis()
  implicit none
  
  type(WavevectorBasis) :: basis
  
  call basis%startup()
end subroutine

! Type representation.
impure elemental function representation_WavevectorBasis() result(output)
  implicit none
  
  type(String) :: output
  
  output = 'wavevector basis'
end function

! Constructor and size function.
function new_WavevectorBasis(maximum_power,expansion_order,subspace_id, &
   & frequency,wavevector,harmonic_states,harmonic_couplings) result(this)
  implicit none
  
  integer,              intent(in) :: maximum_power
  integer,              intent(in) :: expansion_order
  integer,              intent(in) :: subspace_id
  real(dp),             intent(in) :: frequency
  type(FractionVector), intent(in) :: wavevector
  class(SubspaceState), intent(in) :: harmonic_states(:)
  type(CoupledStates),  intent(in) :: harmonic_couplings(:)
  type(WavevectorBasis)            :: this
  
  integer :: ialloc
  
  if (size(harmonic_couplings)/=size(harmonic_states)) then
    call print_line(CODE_ERROR//': harmonic states and harmonic couplings do &
       &not match.')
  endif
  
  this%maximum_power       = maximum_power
  this%expansion_order     = expansion_order
  this%subspace_id         = subspace_id
  this%frequency           = frequency
  this%wavevector          = wavevector
  this%harmonic_couplings_ = harmonic_couplings
  
  allocate( this%harmonic_states_, source=harmonic_states, &
          & stat=ialloc); call err(ialloc)
end function

! Construct a BraKet from a basis.
function new_HarmonicBraKetReal_WavevectorBasis(basis,state,anharmonic_data) &
   & result(this) 
  implicit none
  
  type(WavevectorBasis),   intent(in) :: basis
  type(HarmonicStateReal), intent(in) :: state
  type(AnharmonicData),    intent(in) :: anharmonic_data
  type(HarmonicBraKetReal)            :: this
  
  this = HarmonicBraKetReal(                         &
     & basis%subspace_id,                            &
     & state%mode_ids(),                             &
     & state%paired_mode_ids(),                      &
     & basis%frequency,                              &
     & anharmonic_data%anharmonic_supercell%sc_size, &
     & basis%maximum_power,                          &
     & basis%expansion_order                         )
end function

function new_HarmonicBraKetComplex_WavevectorBasis(basis,state, &
   & anharmonic_data) result(this) 
  implicit none
  
  type(WavevectorBasis),      intent(in) :: basis
  type(HarmonicStateComplex), intent(in) :: state
  type(AnharmonicData),       intent(in) :: anharmonic_data
  type(HarmonicBraKetComplex)            :: this
  
  this = HarmonicBraKetComplex(                      &
     & basis%subspace_id,                            &
     & state%mode_ids(),                             &
     & state%paired_mode_ids(),                      &
     & basis%frequency,                              &
     & anharmonic_data%anharmonic_supercell%sc_size, &
     & basis%maximum_power,                          &
     & basis%expansion_order                         )
end function

! ----------------------------------------------------------------------
! Generates states in a given subspace, up to a given power.
! ----------------------------------------------------------------------
! If qpoint is specified, only modes at that q-point are included.
function new_WavevectorBasis_subspace(subspace,frequency,modes,qpoints, &
   & maximum_power,potential_expansion_order,supercell_size,symmetries, &
   & qpoint) result(output) 
  implicit none
  
  type(DegenerateSubspace), intent(in)           :: subspace
  real(dp),                 intent(in)           :: frequency
  type(ComplexMode),        intent(in)           :: modes(:)
  type(QpointData),         intent(in)           :: qpoints(:)
  integer,                  intent(in)           :: maximum_power
  integer,                  intent(in)           :: potential_expansion_order
  integer,                  intent(in)           :: supercell_size
  type(SymmetryOperator),   intent(in)           :: symmetries(:)
  type(QpointData),         intent(in), optional :: qpoint
  type(WavevectorBasis), allocatable             :: output(:)
  
  ! Variables for generating single-mode bases.
  type(ComplexMode), allocatable :: subspace_modes(:)
  type(WavevectorBasis)          :: mode_basis
  
  ! The basis, before splitting by wavevector.
  ! N.B. basis%wavevector will be zero, which is a placeholder.
  type(WavevectorBasis) :: basis
  
  integer :: i
  
  ! List modes in the subspace, but only including one of each conjugate pair.
  subspace_modes = subspace%modes(modes)
  subspace_modes = subspace_modes(filter(          &
     & subspace_modes%paired_id>=subspace_modes%id ))
  
  ! Filter for modes at the given q-point.
  if (present(qpoint)) then
    subspace_modes = subspace_modes(filter(  &
       & subspace_modes%qpoint_id==qpoint%id ))
  endif
  
  do i=1,size(subspace_modes)
    ! Generate the basis along a single mode (or a mode and its pair).
    mode_basis = generate_mode_basis( subspace,                  &
                                    & frequency,                 &
                                    & subspace_modes(i),         &
                                    & maximum_power,             &
                                    & potential_expansion_order, &
                                    & supercell_size             )
    
    if (i==1) then
      ! Seed the basis as the single-mode basis.
      basis = mode_basis
    else
      ! Expand the basis to include the single-mode basis.
      basis = tensor_product(basis, mode_basis)
    endif
  enddo
  
  output = split_by_wavevector(basis,modes,qpoints,symmetries)
end function

function generate_mode_basis(subspace,frequency,mode,maximum_power, &
   & potential_expansion_order,supercell_size) result(output)
  implicit none
  
  type(DegenerateSubspace), intent(in) :: subspace
  real(dp),                 intent(in) :: frequency
  type(ComplexMode),        intent(in) :: mode
  integer,                  intent(in) :: maximum_power
  integer,                  intent(in) :: potential_expansion_order
  integer,                  intent(in) :: supercell_size
  type(WavevectorBasis)                :: output
  
  if (mode%paired_id==mode%id) then
    output = generate_mode_basis_1d( subspace,                  &
                                   & frequency,                 &
                                   & mode,                      &
                                   & maximum_power,             &
                                   & potential_expansion_order, &
                                   & supercell_size             )
  else
    output = generate_mode_basis_2d( subspace,                  &
                                   & frequency,                 &
                                   & mode,                      &
                                   & maximum_power,             &
                                   & potential_expansion_order, &
                                   & supercell_size             )
  endif
end function

function generate_mode_basis_1d(subspace,frequency,mode,maximum_power, &
   & potential_expansion_order,supercell_size) result(output)
  implicit none
  
  type(DegenerateSubspace), intent(in) :: subspace
  real(dp),                 intent(in) :: frequency
  type(ComplexMode),        intent(in) :: mode
  integer,                  intent(in) :: maximum_power
  integer,                  intent(in) :: potential_expansion_order
  integer,                  intent(in) :: supercell_size
  type(WavevectorBasis)                :: output
  
  type(HarmonicState1D),   allocatable :: harmonic_states_1d(:)
  type(HarmonicStateReal), allocatable :: harmonic_states(:)
  type(CoupledStates),     allocatable :: harmonic_couplings(:)
  
  integer  :: power
  
  integer, allocatable :: coupling_ids(:)
  integer, allocatable :: separations(:)
  
  integer :: i,j,ialloc
  
  harmonic_states_1d = HarmonicState1D( id         = mode%id,           &
                                      & occupation = [( power,          &
                                      &                 power=0,        &
                                      &                 maximum_power)] )
  harmonic_states = [( HarmonicStateReal([harmonic_states_1d(i)]), &
                     & i=1,                                        &
                     & size(harmonic_states_1d)                    )]
  
  ! Calculate which states have non-zero <i|H|j> elements.
  ! |i> = a^n|0> and |j> = a^m|0>.
  ! If |n-m|>expansion_order then <i|H|j>=0.
  allocate( harmonic_couplings(size(harmonic_states_1d)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(harmonic_states_1d)
    coupling_ids = [( j,                                   &
                   & j=max( 1,                             &
                   &        i-potential_expansion_order ), &
                   & min( size(harmonic_states_1d),        &
                   &      i+potential_expansion_order )    )]
    separations = [(abs(coupling_ids(j)-i), j=1, size(coupling_ids))]
    harmonic_couplings(i) = CoupledStates(coupling_ids,separations)
  enddo
  
  ! Construct output.
  output = WavevectorBasis( maximum_power      = maximum_power,             &
                          & expansion_order    = potential_expansion_order, &
                          & subspace_id        = subspace%id,               &
                          & frequency          = frequency,                 &
                          & wavevector         = fracvec(zeroes(3)),        &
                          & harmonic_states    = harmonic_states,           &
                          & harmonic_couplings = harmonic_couplings         )
end function

function generate_mode_basis_2d(subspace,frequency,mode,maximum_power, &
   & potential_expansion_order,supercell_size) result(output)
  implicit none
  
  type(DegenerateSubspace), intent(in) :: subspace
  real(dp),                 intent(in) :: frequency
  type(ComplexMode),        intent(in) :: mode
  integer,                  intent(in) :: maximum_power
  integer,                  intent(in) :: potential_expansion_order
  integer,                  intent(in) :: supercell_size
  type(WavevectorBasis)                :: output
  
  type(HarmonicState2D),      allocatable :: harmonic_states_2d(:)
  type(HarmonicStateComplex), allocatable :: harmonic_states(:)
  type(CoupledStates),        allocatable :: harmonic_couplings(:)
  
  integer :: total_power
  
  integer  :: power
  integer  :: separation
  
  integer, allocatable :: coupling_ids(:)
  integer, allocatable :: separations(:)
  
  integer :: i,j,ialloc
  
  
  harmonic_states_2d =                                                     &
     & [( [( HarmonicState2D( id                = mode%id,                 &
     &                        paired_id         = mode%paired_id,          &
     &                        occupation        = power,                   &
     &                        paired_occupation = total_power-power ),     &
     &       power = total_power,                                          &
     &       0,                                                            &
     &       -1                                                        )], &
     & total_power=0,                                                      &
     & maximum_power                                                       )]
  harmonic_states = [( HarmonicStateComplex([harmonic_states_2d(i)]), &
                     & i=1,                                           &
                     & size(harmonic_states_2d)                       )]
  
  ! Calculate which states have non-zero <i|H|j> elements.
  ! |i> = (a+)^(m+).(a-)^(m-).|0>
  ! |j> = (a+)^(n+).(a-)^(n-).|0>
  ! If |m+ - n+| + |m- - n-| > potential_expansion_order then <i|H|j>=0.
  allocate( harmonic_couplings(size(harmonic_states_2d)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(harmonic_states_2d)
    coupling_ids = [integer::]
    separations = [integer::]
    do j=1,size(harmonic_states_2d)
      separation = abs( harmonic_states_2d(i)%occupation()        &
               &      - harmonic_states_2d(j)%occupation() )      &
               & + abs( harmonic_states_2d(i)%paired_occupation() &
               &      - harmonic_states_2d(j)%paired_occupation() )
      if (separation<=potential_expansion_order) then
        coupling_ids = [coupling_ids, j]
        separations = [separations, separation]
      endif
    enddo
    harmonic_couplings(i) = CoupledStates(coupling_ids, separations)
  enddo
  
  ! Generate output.
  output = WavevectorBasis( maximum_power      = maximum_power,             &
                          & expansion_order    = potential_expansion_order, &
                          & subspace_id        = subspace%id,               &
                          & frequency          = frequency,                 &
                          & wavevector         = fracvec(zeroes(3)),        &
                          & harmonic_states    = harmonic_states,           &
                          & harmonic_couplings = harmonic_couplings         )
end function

! Takes the basis over one set of modes in the subspace, and the basis over a
!    disjoint set of modes in the same subspace, and constructs the basis
!    over the union of the two sets, up to the maximum power.
! e.g. if 'this has states |0>, |u1> and |u1^2>,
!    and 'that' has states |0>, |u2> and |u2^2>, then the output will have
!    states |0>, |u1>, |u2>, |u1^2>, |u1u2> and |u2^2>.
! N.B. the states of 'that' must be sorted in ascending order of total power.
function tensor_product(this,that) result(output)
  implicit none
  
  type(WavevectorBasis), intent(in) :: this
  type(WavevectorBasis), intent(in) :: that
  type(WavevectorBasis)             :: output
  
  integer, allocatable :: this_occupations(:)
  integer, allocatable :: that_occupations(:)
  integer, allocatable :: no_states(:)
  integer, allocatable :: this_to_output(:)
  
  type(CoupledStates) :: mapped_coupling
  
  ! Output variables.
  integer                                 :: maximum_power
  integer                                 :: expansion_order
  integer                                 :: subspace_id
  real(dp)                                :: frequency
  type(FractionVector)                    :: wavevector
  type(HarmonicStateReal),    allocatable :: harmonic_states_real(:)
  type(HarmonicStateComplex), allocatable :: harmonic_states_complex(:)
  class(SubspaceState),       allocatable :: harmonic_states(:)
  type(CoupledStates),        allocatable :: harmonic_couplings(:)
  
  integer :: i,j,ialloc
  
  ! Check that 'this' and 'that' are consistent.
  if (this%expansion_order/=that%expansion_order) then
    call print_line(CODE_ERROR//': Expansion orders do not match.')
    call err()
  elseif (this%subspace_id/=that%subspace_id) then
    call print_line(CODE_ERROR//': Subspaces do not match.')
    call err()
  endif
  
  maximum_power   = this%maximum_power
  expansion_order = this%expansion_order
  subspace_id     = this%subspace_id
  frequency       = this%frequency
  wavevector      = this%wavevector
  
  ! Check that the states in 'that' are in ascending order
  !    of total occupation.
  this_occupations = this%harmonic_states_%occupation()
  that_occupations = that%harmonic_states_%occupation()
  if ( any(that_occupations(2:)                          &
   & < that_occupations(:size(that%harmonic_states_)-1)) ) then
    call print_line(CODE_ERROR//': Powers are not in ascending order.')
    call err()
  endif
  
  ! Count the number of output states corresponding to each state in 'this'.
  no_states = [( count( that_occupations+this_occupations(i)    &
               &    <= maximum_power                         ), &
               & i=1,                                           &
               & size(this%harmonic_states_)                    )]
  this_to_output = [( sum(no_states(:i-1)),        &
                    & i=1,                         &
                    & size(this%harmonic_states_) )]
  
  ! Generate monomial and harmonic states.
  associate( this_states=>this%harmonic_states_, &
           & that_states=>that%harmonic_states_  )
    select type(this_states); type is(HarmonicStateReal)
      select type(that_states); type is(HarmonicStateReal)
        harmonic_states_real = [( [( prod_real( this_states(i),        &
                                &               that_states(j)  ),     &
                                &    j=1,                              &
                                &    no_states(i)                  )], &
                                & i=1,                                 &
                                & size(this_states)                    )]
        allocate( harmonic_states, source=harmonic_states_real, &
                & stat=ialloc); call err(ialloc)
      end select
    type is(HarmonicStateComplex)
      select type(that_states); type is(HarmonicStateComplex)
        harmonic_states_complex = [( [( prod_complex( this_states(i),        &
                                   &                  that_states(j)  ),     &
                                   &    j=1,                                 &
                                   &    no_states(i)                     )], &
                                   & i=1,                                    &
                                   & size(this_states)                       )]
        allocate( harmonic_states, source=harmonic_states_complex, &
                & stat=ialloc); call err(ialloc)
      end select
    end select
  end associate
  
  ! Generate harmonic couplings.
  allocate( harmonic_couplings(size(harmonic_states)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(this%harmonic_states_)
    mapped_coupling = this%harmonic_couplings_(i)
    call mapped_coupling%map_ids(this_to_output)
    do j=1,no_states(i)
      harmonic_couplings(this_to_output(i)+j) = tensor_product_couplings( &
                                           & this%harmonic_couplings_(i), &
                                           & mapped_coupling,             &
                                           & no_states,                   &
                                           & that%harmonic_couplings_(j), &
                                           & expansion_order              )
    enddo
  enddo
  
  ! Update output.
  output = WavevectorBasis( maximum_power      = maximum_power,     &
                          & expansion_order    = expansion_order,   &
                          & subspace_id        = subspace_id,       &
                          & frequency          = frequency,         &
                          & wavevector         = wavevector,        &
                          & harmonic_states    = harmonic_states,   &
                          & harmonic_couplings = harmonic_couplings )
end function

! Splits up a WavevectorBasis by wavevector.
function split_by_wavevector(input,modes,qpoints,symmetries) result(output)
  implicit none
  
  type(WavevectorBasis),  intent(in) :: input
  type(ComplexMode),      intent(in) :: modes(:)
  type(QpointData),       intent(in) :: qpoints(:)
  type(SymmetryOperator), intent(in) :: symmetries(:)
  type(WavevectorBasis), allocatable :: output(:)
  
  type(QpointData),     allocatable :: wavevector_qpoints(:)
  type(QpointData),     allocatable :: basis_qpoints(:)
  
  integer, allocatable :: wavevector_states(:)
  integer, allocatable :: new_ids(:)
  integer, allocatable :: ids(:)
  integer, allocatable :: separations(:)
  integer, allocatable :: allowed_couplings(:)
  
  ! Output variables.
  type(HarmonicStateReal),    allocatable :: harmonic_states_real(:)
  type(HarmonicStateComplex), allocatable :: harmonic_states_complex(:)
  class(SubspaceState),       allocatable :: harmonic_states(:)
  type(CoupledStates),        allocatable :: couplings(:)
  
  integer :: i,j,ialloc
  
  wavevector_qpoints = [(                                                  &
     & QpointData(                                                         &
     &    qpoint = input%harmonic_states_(i)%wavevector(modes,qpoints),    &
     &    id               = 0,                                            &
     &    paired_qpoint_id = 0                                          ), &
     & i=1,                                                                &
     & size(input%harmonic_states_)                                        )]
  
  basis_qpoints = qpoints(filter(                 &
     & [( any(wavevector_qpoints==qpoints(i)),    &
     &    i=1,                                    &
     &    size(qpoints)                        )] ))
  
  allocate(output(size(basis_qpoints)), stat=ialloc); call err(ialloc)
  new_ids = [(0,i=1,size(input%harmonic_states_))]
  do i=1,size(output)
    ! Identify the states at the given wavevector.
    wavevector_states = filter(wavevector_qpoints==basis_qpoints(i))
    couplings         = input%harmonic_couplings_(wavevector_states)
    associate(states=>input%harmonic_states_)
      select type(states); type is(HarmonicStateReal)
        harmonic_states_real = states(wavevector_states)
        allocate( harmonic_states, source=harmonic_states_real, &
                & stat=ialloc); call err(ialloc)
      type is(HarmonicStateComplex)
        harmonic_states_complex = states(wavevector_states)
        allocate( harmonic_states, source=harmonic_states_complex, &
                & stat=ialloc); call err(ialloc)
      end select
    end associate
    
    ! Re-map IDs to correspond to states at the given wavevector only.
    new_ids = 0
    new_ids(wavevector_states) = [(j,j=1,size(wavevector_states))]
    do j=1,size(wavevector_states)
      ids = couplings(j)%ids()
      separations = couplings(j)%separations()
      ids = new_ids(ids)
      allowed_couplings = filter(ids>0)
      couplings(j) = CoupledStates( ids(allowed_couplings),        &
                                  & separations(allowed_couplings) )
    enddo
    
    ! Construct output.
    output(i) = WavevectorBasis( input%maximum_power,     &
                               & input%expansion_order,   &
                               & input%subspace_id,       &
                               & input%frequency,         &
                               & basis_qpoints(i)%qpoint, &
                               & harmonic_states,         &
                               & couplings                )
    deallocate(harmonic_states, stat=ialloc); call err(ialloc)
  enddo
end function

! ----------------------------------------------------------------------
! Operations involving WavevectorStates.
! ----------------------------------------------------------------------
impure elemental function inner_product_WavevectorBasis(this,bra,ket, &
   & subspace,anharmonic_data) result(output)
  implicit none
  
  class(WavevectorBasis),   intent(in)                   :: this
  class(BasisState),        intent(in),           target :: bra
  class(BasisState),        intent(in), optional, target :: ket
  type(DegenerateSubspace), intent(in)                   :: subspace
  type(AnharmonicData),     intent(in)                   :: anharmonic_data
  real(dp)                                               :: output
  
  type(WavevectorState) :: bra_
  type(WavevectorState) :: ket_
  
  bra_ = WavevectorState(bra)
  
  if (present(ket)) then
    ket_ = WavevectorState(ket)
    if (any(bra_%state_ids/=ket_%state_ids)) then
      call print_line(CODE_ERROR//': Bra and Ket do not match.')
      call err()
    endif
    ! The basis is orthonormal, so the dot product of the states is
    !    simply the dot product of the coefficients.
    output = dot_product(bra_%coefficients, ket_%coefficients)
  else
    ! States are orthonormal, so <p|p>=1.
    output = 1.0_dp
  endif
end function

impure elemental function integrate_BasisState_WavevectorBasis(this,bra, &
   & monomial,ket,subspace,anharmonic_data) result(output)
  implicit none
  
  class(WavevectorBasis),   intent(in)                   :: this
  class(BasisState),        intent(in),           target :: bra
  type(SparseMonomial),     intent(in)                   :: monomial
  class(BasisState),        intent(in), optional, target :: ket
  type(DegenerateSubspace), intent(in)                   :: subspace
  type(AnharmonicData),     intent(in)                   :: anharmonic_data
  complex(dp)                                            :: output
  
  type(WavevectorState) :: bra_
  type(WavevectorState) :: ket_
  
  bra_ = WavevectorState(bra)
  if (present(ket)) then
    ket_ = WavevectorState(ket)
  else
    ket_ = bra_
  endif
  
  if (any(bra_%state_ids/=ket_%state_ids)) then
    call print_line(CODE_ERROR//': Bra and Ket do not match.')
    call err()
  endif
  
  associate(harmonic_states=>this%harmonic_states_)
    select type(harmonic_states); type is(HarmonicStateReal)
      output = calculate_integral_real( this,            &
                                      & harmonic_states, &
                                      & bra_,            &
                                      & ket_,            &
                                      & monomial,        &
                                      & anharmonic_data  )
    type is(HarmonicStateComplex)
      output = calculate_integral_complex( this,            &
                                         & harmonic_states, &
                                         & bra_,            &
                                         & ket_,            &
                                         & monomial,        &
                                         & anharmonic_data  )
    end select
  end associate
end function

function calculate_integral_real(basis,harmonic_states,bra,ket,monomial, &
   & anharmonic_data) result(output) 
  implicit none
  
  type(WavevectorBasis),   intent(in)         :: basis
  type(HarmonicStateReal), intent(in), target :: harmonic_states(:)
  type(WavevectorState),   intent(in)         :: bra
  type(WavevectorState),   intent(in)         :: ket
  type(SparseMonomial),    intent(in)         :: monomial
  type(AnharmonicData),    intent(in)         :: anharmonic_data
  complex(dp)                                 :: output
  
  type(HarmonicBraKetReal) :: braket
  
  complex(dp) :: term
  
  integer :: i,j,k,l
  
  ! Initialise the braket. This can be done from any state.
  braket = HarmonicBraKetReal(basis,harmonic_states(1),anharmonic_data)
  
  ! Calculate the integral of the potential.
  output = 0.0_dp
  do i=1,size(bra%state_ids)
    braket%bra_ => harmonic_states(bra%state_ids(i))
    do j=1,size(basis%harmonic_couplings_(bra%state_ids(i)))
      k = basis%harmonic_couplings_(i)%id(j)
      
      ! Ignore terms with k<i. These are added when i and k are reversed.
      if (k<i) then
        cycle
      endif
      
      l = first(bra%state_ids==k, default=0)
      if (l==0) then
        cycle
      endif
      
      braket%ket_ => harmonic_states(k)
      
      term = bra%coefficients(i)                         &
         & * braket%integrate(monomial, anharmonic_data) &
         & * ket%coefficients(l)
      
      if (k==i) then
        output = output + term
      else
        ! Include both <k||i> and <i||k>.
        output = output + 2*term
      endif
    enddo
  enddo
end function

function calculate_integral_complex(basis,harmonic_states,bra,ket,monomial, &
   & anharmonic_data) result(output) 
  implicit none
  
  type(WavevectorBasis),      intent(in)         :: basis
  type(HarmonicStateComplex), intent(in), target :: harmonic_states(:)
  type(WavevectorState),      intent(in)         :: bra
  type(WavevectorState),      intent(in)         :: ket
  type(SparseMonomial),       intent(in)         :: monomial
  type(AnharmonicData),       intent(in)         :: anharmonic_data
  complex(dp)                                    :: output
  
  type(HarmonicBraKetComplex) :: braket
  
  complex(dp) :: term
  
  integer :: i,j,k,l
  
  ! Initialise the braket. This can be done from any state.
  braket = HarmonicBraKetComplex(basis, harmonic_states(1), anharmonic_data)
  
  ! Calculate the integral of the potential.
  output = 0.0_dp
  do i=1,size(harmonic_states)
    braket%bra_ => harmonic_states(i)
    do j=1,size(basis%harmonic_couplings_(i))
      k = basis%harmonic_couplings_(i)%id(j)
      
      ! Ignore terms with k<i. These are added when i and k are reversed.
      if (k<i) then
        cycle
      endif
      
      l = first(bra%state_ids==k, default=0)
      if (l==0) then
        cycle
      endif
      
      braket%ket_ => harmonic_states(k)
      
      term = bra%coefficients(i)                         &
         & * braket%integrate(monomial, anharmonic_data) &
         & * ket%coefficients(l)
      
      if (k==i) then
        output = output + term
      else
        ! Include both <k||i> and <i||k>.
        output = output + 2*term
      endif
    enddo
  enddo
end function

! Integrate a monomial between sets of states.
impure elemental function integrate_BasisStates_WavevectorBasis(this,states, &
   & monomial,subspace,anharmonic_data) result(output)
  implicit none
  
  class(WavevectorBasis),   intent(in)         :: this
  class(BasisStates),       intent(in), target :: states
  type(SparseMonomial),     intent(in)         :: monomial
  type(DegenerateSubspace), intent(in)         :: subspace
  type(AnharmonicData),     intent(in)         :: anharmonic_data
  complex(dp)                                  :: output
  
  type(WavevectorStates), pointer :: states_
  
  integer :: i
  
  if (size(this%harmonic_States_)==0) then
    output = 0.0_dp
    return
  endif
  
  ! Convert the states to type WavevectorStates.
  states_ => wavevector_states_pointer(states)
  
  ! Find the density matrix corresponding to this wavevector.
  i = first( states_%density_matrices%wavevector==this%wavevector, &
           & default=0                                             )
  if (i==0) then
    output = 0.0_dp
    return
  elseif (.not. allocated(states_%density_matrices(i)%values)) then
    output = 0.0_dp
    return
  endif
  
  ! Integrate the monomial w/r/t the density matrix.
  associate(harmonic_states=>this%harmonic_states_)
    select type(harmonic_states); type is(HarmonicStateReal)
      output = integrate_real( this,                        &
           &                   harmonic_states,             &
           &                   states_%density_matrices(i), &
           &                   monomial,                    &
           &                   anharmonic_data              )
    type is(HarmonicStateComplex)
      output = integrate_complex( this,                        &
           &                      harmonic_states,             &
           &                      states_%density_matrices(i), &
           &                      monomial,                    &
           &                      anharmonic_data              )
    end select
  end associate
end function

function integrate_real(basis,harmonic_states,density_matrix,monomial, &
   & anharmonic_data) result(output) 
  implicit none
  
  type(WavevectorBasis),   intent(in)         :: basis
  type(HarmonicStateReal), intent(in), target :: harmonic_states(:)
  type(DensityMatrix),     intent(in)         :: density_matrix
  type(SparseMonomial),    intent(in)         :: monomial
  type(AnharmonicData),    intent(in)         :: anharmonic_data
  complex(dp)                                 :: output
  
  type(HarmonicBraKetReal) :: braket
  
  complex(dp) :: term
  
  integer :: i
  
  ! Initialise the braket. This can be done from any state.
  braket = HarmonicBraKetReal(basis, harmonic_states(1), anharmonic_data)
  
  output = 0
  do i=1,size(density_matrix%values)
    ! Ignore terms with ket<bra. These are added when i and k are reversed.
    if (density_matrix%ket_ids(i)<density_matrix%bra_ids(i)) then
      cycle
    endif
    
    ! Calculate <bra|X|ket> in the harmonic basis,
    !    and thermally weight by the density matrix.
    braket%bra_ => harmonic_states(density_matrix%bra_ids(i))
    braket%ket_ => harmonic_states(density_matrix%ket_ids(i))
    term = braket%integrate(monomial, anharmonic_data) &
       & * density_matrix%values(i)
    
    ! If bra==ket, only include <bra|X|bra>,
    !    otherwise include <bra|X|ket> and <ket|X|bra>.
    if (density_matrix%bra_ids(i)==density_matrix%ket_ids(i)) then
      output = output + term
    else
      output = output + 2*term
    endif
  enddo
end function

function integrate_complex(basis,harmonic_states,density_matrix,monomial, &
   & anharmonic_data) result(output) 
  implicit none
  
  type(WavevectorBasis),      intent(in)         :: basis
  type(HarmonicStateComplex), intent(in), target :: harmonic_states(:)
  type(DensityMatrix),        intent(in)         :: density_matrix
  type(SparseMonomial),       intent(in)         :: monomial
  type(AnharmonicData),       intent(in)         :: anharmonic_data
  complex(dp)                                    :: output
  
  type(HarmonicBraKetComplex) :: braket
  
  complex(dp) :: term
  
  integer :: i
  
  ! Initialise the braket. This can be done from any state.
  braket = HarmonicBraKetComplex(basis, harmonic_states(1), anharmonic_data)
  
  output = 0
  do i=1,size(density_matrix%values)
    ! Ignore terms with ket<bra. These are added when i and k are reversed.
    if (density_matrix%ket_ids(i)<density_matrix%bra_ids(i)) then
      cycle
    endif
    
    ! Calculate <bra|X|ket> in the harmonic basis,
    !    and thermally weight by the density matrix.
    braket%bra_ => harmonic_states(density_matrix%bra_ids(i))
    braket%ket_ => harmonic_states(density_matrix%ket_ids(i))
    term = braket%integrate(monomial, anharmonic_data) &
       & * density_matrix%values(i)
    
    ! If bra==ket, only include <bra|X|bra>,
    !    otherwise include <bra|X|ket> and <ket|X|bra>.
    if (density_matrix%bra_ids(i)==density_matrix%ket_ids(i)) then
      output = output + term
    else
      output = output + 2*term
    endif
  enddo
end function

impure elemental function kinetic_energy_WavevectorBasis(this,bra,ket, &
   & subspace,anharmonic_data) result(output)
  implicit none
  
  class(WavevectorBasis),   intent(in)                   :: this
  class(BasisState),        intent(in),           target :: bra
  class(BasisState),        intent(in), optional, target :: ket
  type(DegenerateSubspace), intent(in)                   :: subspace
  type(AnharmonicData),     intent(in)                   :: anharmonic_data
  real(dp)                                               :: output
  
  type(WavevectorState) :: bra_
  type(WavevectorState) :: ket_
  
  bra_ = WavevectorState(bra)
  if (present(ket)) then
    ket_ = WavevectorState(ket)
  else
    ket_ = bra_
  endif
  
  associate(harmonic_states=>this%harmonic_states_)
    select type(harmonic_states); type is(HarmonicStateReal)
      output = kinetic_energy_real( this,            &
                                  & harmonic_states, &
                                  & bra_,            &
                                  & ket_,            &
                                  & anharmonic_data  )
    type is(HarmonicStateComplex)
      output = kinetic_energy_complex( this,            &
                                     & harmonic_states, &
                                     & bra_,            &
                                     & ket_,            &
                                     & anharmonic_data  )
    end select
  end associate
end function

function kinetic_energy_real(basis,harmonic_states,bra,ket,anharmonic_data) &
   & result(output) 
  implicit none
  
  type(WavevectorBasis),   intent(in)         :: basis
  type(HarmonicStateReal), intent(in), target :: harmonic_states(:)
  type(WavevectorState),   intent(in)         :: bra
  type(WavevectorState),   intent(in)         :: ket
  type(AnharmonicData),    intent(in)         :: anharmonic_data
  real(dp)                                    :: output
  
  type(HarmonicBraKetReal) :: braket
  
  real(dp) :: term
  
  integer :: i,j,k
  
  ! Initialise the braket. This can be done from any state.
  braket = HarmonicBraKetReal(basis, harmonic_states(1), anharmonic_data)
  
  ! Calculate the kinetic energy.
  output = 0.0_dp
  do i=1,size(harmonic_states)
    braket%bra_ => harmonic_states(i)
    do j=1,size(basis%harmonic_couplings_(i))
      k = basis%harmonic_couplings_(i)%id(j)
      
      ! Ignore terms with k<i. These are added when i and k are reversed.
      if (k<i) then
        cycle
      endif
      
      braket%ket_ => harmonic_states(k)
      
      term = bra%coefficients(i)                    &
         & * braket%kinetic_energy(anharmonic_data) &
         & * ket%coefficients(k)
      
      if (k==i) then
        output = output + term
      else
        ! Include both <k||i> and <i||k>.
        output = output + 2*term
      endif
    enddo
  enddo
end function

function kinetic_energy_complex(basis,harmonic_states,bra,ket, &
   & anharmonic_data) result(output) 
  implicit none
  
  type(WavevectorBasis),      intent(in)         :: basis
  type(HarmonicStateComplex), intent(in), target :: harmonic_states(:)
  type(WavevectorState),      intent(in)         :: bra
  type(WavevectorState),      intent(in)         :: ket
  type(AnharmonicData),       intent(in)         :: anharmonic_data
  real(dp)                                       :: output
  
  type(HarmonicBraKetComplex) :: braket
  
  real(dp) :: term
  
  integer :: i,j,k
  
  ! Initialise the braket. This can be done from any state.
  braket = HarmonicBraKetComplex(basis, harmonic_states(1), anharmonic_data)
  
  ! Calculate the kinetic energy.
  output = 0.0_dp
  do i=1,size(harmonic_states)
    braket%bra_ => harmonic_states(i)
    do j=1,size(basis%harmonic_couplings_(i))
      k = basis%harmonic_couplings_(i)%id(j)
      
      ! Ignore terms with k<i. These are added when i and k are reversed.
      if (k<i) then
        cycle
      endif
      
      braket%ket_ => harmonic_states(k)
      
      term = bra%coefficients(i)                    &
         & * braket%kinetic_energy(anharmonic_data) &
         & * ket%coefficients(k)
      
      if (k==i) then
        output = output + term
      else
        ! Include both <k||i> and <i||k>.
        output = output + 2*term
      endif
    enddo
  enddo
end function

impure elemental function harmonic_potential_energy_WavevectorBasis(this, &
   & bra,ket,subspace,anharmonic_data) result(output)
  implicit none
  
  class(WavevectorBasis),   intent(in)                   :: this
  class(BasisState),        intent(in),           target :: bra
  class(BasisState),        intent(in), optional, target :: ket
  type(DegenerateSubspace), intent(in)                   :: subspace
  type(AnharmonicData),     intent(in)                   :: anharmonic_data
  real(dp)                                               :: output
  
  type(WavevectorState) :: bra_
  type(WavevectorState) :: ket_
  
  bra_ = WavevectorState(bra)
  if (present(ket)) then
    ket_ = WavevectorState(ket)
  else
    ket_ = bra_
  endif
  
  associate(harmonic_states=>this%harmonic_states_)
    select type(harmonic_states); type is(HarmonicStateReal)
      output = harmonic_potential_energy_real( this,            &
                                             & harmonic_states, &
                                             & bra_,            &
                                             & ket_,            &
                                             & anharmonic_data  )
    type is(HarmonicStateComplex)
      output = harmonic_potential_energy_complex( this,            &
                                                & harmonic_states, &
                                                & bra_,            &
                                                & ket_,            &
                                                & anharmonic_data  )
    end select
  end associate
end function

function harmonic_potential_energy_real(basis,harmonic_states,bra,ket, &
   & anharmonic_data) result(output) 
  implicit none
  
  type(WavevectorBasis),   intent(in)         :: basis
  type(HarmonicStateReal), intent(in), target :: harmonic_states(:)
  type(WavevectorState),   intent(in)         :: bra
  type(WavevectorState),   intent(in)         :: ket
  type(AnharmonicData),    intent(in)         :: anharmonic_data
  real(dp)                                    :: output
  
  type(HarmonicBraKetReal) :: braket
  
  real(dp) :: term
  
  integer :: i,j,k
  
  ! Initialise the braket. This can be done from any state.
  braket = HarmonicBraKetReal(basis, harmonic_states(1), anharmonic_data)
  
  ! Calculate the integral of the potential.
  output = 0.0_dp
  do i=1,size(harmonic_states)
    braket%bra_ => harmonic_states(i)
    do j=1,size(basis%harmonic_couplings_(i))
      k = basis%harmonic_couplings_(i)%id(j)
      
      ! Ignore terms with k<i. These are added when i and k are reversed.
      if (k<i) then
        cycle
      endif
      
      braket%ket_ => harmonic_states(k)
      
      term = bra%coefficients(i)                               &
         & * braket%harmonic_potential_energy(anharmonic_data) &
         & * ket%coefficients(k)
      
      if (k==i) then
        output = output + term
      else
        ! Include both <k||i> and <i||k>.
        output = output + 2*term
      endif
    enddo
  enddo
end function

function harmonic_potential_energy_complex(basis,harmonic_states,bra,ket, &
   & anharmonic_data) result(output) 
  implicit none
  
  type(WavevectorBasis),      intent(in)         :: basis
  type(HarmonicStateComplex), intent(in), target :: harmonic_states(:)
  type(WavevectorState),      intent(in)         :: bra
  type(WavevectorState),      intent(in)         :: ket
  type(AnharmonicData),       intent(in)         :: anharmonic_data
  real(dp)                                       :: output
  
  type(HarmonicBraKetComplex) :: braket
  
  real(dp) :: term
  
  integer :: i,j,k
  
  ! Initialise the braket. This can be done from any state.
  braket = HarmonicBraKetComplex(basis, harmonic_states(1), anharmonic_data)
  
  ! Calculate the integral of the potential.
  output = 0.0_dp
  do i=1,size(harmonic_states)
    braket%bra_ => harmonic_states(i)
    do j=1,size(basis%harmonic_couplings_(i))
      k = basis%harmonic_couplings_(i)%id(j)
      
      ! Ignore terms with k<i. These are added when i and k are reversed.
      if (k<i) then
        cycle
      endif
      
      braket%ket_ => harmonic_states(k)
      
      term = bra%coefficients(i)                               &
         & * braket%harmonic_potential_energy(anharmonic_data) &
         & * ket%coefficients(k)
      
      if (k==i) then
        output = output + term
      else
        ! Include both <k||i> and <i||k>.
        output = output + 2*term
      endif
    enddo
  enddo
end function

impure elemental function kinetic_stress_WavevectorBasis(this,bra,ket, &
   & subspace,stress_prefactors,anharmonic_data) result(output)
  implicit none
  
  class(WavevectorBasis),   intent(in)                   :: this
  class(BasisState),        intent(in),           target :: bra
  class(BasisState),        intent(in), optional, target :: ket
  type(DegenerateSubspace), intent(in)                   :: subspace
  type(StressPrefactors),   intent(in)                   :: stress_prefactors
  type(AnharmonicData),     intent(in)                   :: anharmonic_data
  type(RealMatrix)                                       :: output
  
  type(WavevectorState) :: bra_
  type(WavevectorState) :: ket_
  
  bra_ = WavevectorState(bra)
  if (present(ket)) then
    ket_ = WavevectorState(ket)
  else
    ket_ = bra_
  endif
  
  associate(harmonic_states=>this%harmonic_states_)
    select type(harmonic_states); type is(HarmonicStateReal)
      output = kinetic_stress_real( this,              &
                                  & harmonic_states,   &
                                  & bra_,              &
                                  & ket_,              &
                                  & stress_prefactors, &
                                  & anharmonic_data    )
    type is(HarmonicStateComplex)
      output = kinetic_stress_complex( this,              &
                                     & harmonic_states,   &
                                     & bra_,              &
                                     & ket_,              &
                                     & stress_prefactors, &
                                     & anharmonic_data    )
    end select
  end associate
end function

function kinetic_stress_real(basis,harmonic_states,bra,ket, &
   & stress_prefactors,anharmonic_data) result(output) 
  implicit none
  
  type(WavevectorBasis),   intent(in)         :: basis
  type(HarmonicStateReal), intent(in), target :: harmonic_states(:)
  type(WavevectorState),   intent(in)         :: bra
  type(WavevectorState),   intent(in)         :: ket
  type(StressPrefactors),  intent(in)         :: stress_prefactors
  type(AnharmonicData),    intent(in)         :: anharmonic_data
  type(RealMatrix)                            :: output
  
  type(HarmonicBraKetReal) :: braket
  
  type(RealMatrix) :: term
  
  integer :: i,j,k
  
  ! Initialise the braket. This can be done from any state.
  braket = HarmonicBraKetReal(basis, harmonic_states(1), anharmonic_data)
  
  ! Calculate the kinetic stress.
  output = dblemat(zeroes(3,3))
  do i=1,size(harmonic_states)
    braket%bra_ => harmonic_states(i)
    do j=1,size(basis%harmonic_couplings_(i))
      k = basis%harmonic_couplings_(i)%id(j)
      
      ! Ignore terms with k<i. These are added when i and k are reversed.
      if (k<i) then
        cycle
      endif
      
      braket%ket_ => harmonic_states(k)
      
      term = bra%coefficients(i)                                       &
         & * braket%kinetic_stress(stress_prefactors, anharmonic_data) &
         & * ket%coefficients(k)
      
      if (k==i) then
        output = output + term
      else
        ! Include both <k||i> and <i||k>.
        output = output + 2*term
      endif
    enddo
  enddo
end function

function kinetic_stress_complex(basis,harmonic_states,bra,ket, &
   & stress_prefactors,anharmonic_data) result(output) 
  implicit none
  
  type(WavevectorBasis),      intent(in)         :: basis
  type(HarmonicStateComplex), intent(in), target :: harmonic_states(:)
  type(WavevectorState),      intent(in)         :: bra
  type(WavevectorState),      intent(in)         :: ket
  type(StressPrefactors),     intent(in)         :: stress_prefactors
  type(AnharmonicData),       intent(in)         :: anharmonic_data
  type(RealMatrix)                               :: output
  
  type(HarmonicBraKetComplex) :: braket
  
  type(RealMatrix) :: term
  
  integer :: i,j,k
  
  ! Initialise the braket. This can be done from any state.
  braket = HarmonicBraKetComplex(basis, harmonic_states(1), anharmonic_data)
  
  ! Calculate the kinetic stress.
  output = dblemat(zeroes(3,3))
  do i=1,size(harmonic_states)
    braket%bra_ => harmonic_states(i)
    do j=1,size(basis%harmonic_couplings_(i))
      k = basis%harmonic_couplings_(i)%id(j)
      
      ! Ignore terms with k<i. These are added when i and k are reversed.
      if (k<i) then
        cycle
      endif
      
      braket%ket_ => harmonic_states(k)
      
      term = bra%coefficients(i)                                       &
         & * braket%kinetic_stress(stress_prefactors, anharmonic_data) &
         & * ket%coefficients(k)
      
      if (k==i) then
        output = output + term
      else
        ! Include both <k||i> and <i||k>.
        output = output + 2*term
      endif
    enddo
  enddo
end function

! ----------------------------------------------------------------------
! Operations which generate WavevectorStates.
! ----------------------------------------------------------------------
! Generate initial guess. This is simply the basis state |0>, i.e. the
!    ground state of the effective harmonic potential from which the basis
!    states were generated.
impure elemental function initial_ground_state_WavevectorBasis(this) &
   & result(output)
  implicit none
  
  class(WavevectorBasis), intent(in) :: this
  type(WavevectorState)              :: output
  
  integer :: state_id
  
  state_id = first(this%harmonic_states_%occupation()==0)
  
  output = WavevectorState( this%subspace_id, &
                          & this%wavevector,  &
                          & [state_id],       &
                          & [1.0_dp]          )
end function

impure elemental function initial_states_WavevectorBasis(this,subspace, &
   & thermal_energy,anharmonic_data) result(output)
  implicit none
  
  class(WavevectorBasis),   intent(in) :: this
  type(DegenerateSubspace), intent(in) :: subspace
  real(dp),                 intent(in) :: thermal_energy
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(BasisStatesPointer)             :: output
  
  type(WavevectorState)  :: state
  type(WavevectorStates) :: states
  
  state = this%initial_ground_state()
  
  states = WavevectorStates( subspace_id     = this%subspace_id, &
                           & states          = [state],          &
                           & energies        = [0.0_dp]          )
  
  states%density_matrices = [DensityMatrix(   &
     & wavevector = this%wavevector,          &
     & couplings  = this%harmonic_couplings_, &
     & states     = [state],                  &
     & weights    = [1.0_dp]                  )]

  output = BasisStatesPointer(states)
end function

impure elemental function calculate_states_WavevectorBasis(this,subspace,    &
   & subspace_potential,thermal_energy,state_energy_cutoff,convergence_data, &
   & anharmonic_data) result(output) 
  implicit none
  
  class(WavevectorBasis),   intent(in) :: this
  type(DegenerateSubspace), intent(in) :: subspace
  class(PotentialBase),     intent(in) :: subspace_potential
  real(dp),                 intent(in) :: thermal_energy
  real(dp),                 intent(in) :: state_energy_cutoff
  type(ConvergenceData),    intent(in) :: convergence_data
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(BasisStatesPointer)             :: output
  
  type(SelectedStatesHamiltonian)        :: hamiltonian
  type(SymmetricEigenstuff), allocatable :: estuff(:)
  
  type(WavevectorState), allocatable :: wavevector_states(:)
  
  integer :: i,ialloc
  
  ! Calculate the Hamiltonian.
  associate(harmonic_states=>this%harmonic_states_)
    select type(harmonic_states); type is(HarmonicStateReal)
      hamiltonian = calculate_hamiltonian_real( &
            & this,                             &
            & harmonic_states,                  &
            & state_energy_cutoff,              &
            & subspace_potential,               &
            & anharmonic_data = anharmonic_data )
    type is(HarmonicStateComplex)
      hamiltonian = calculate_hamiltonian_complex( &
               & this,                             &
               & harmonic_states,                  &
               & state_energy_cutoff,              &
               & subspace_potential,               &
               & anharmonic_data = anharmonic_data )
    end select
  end associate
  
  estuff = diagonalise_symmetric(mat(hamiltonian%hamiltonian))
  
  allocate(wavevector_states(size(estuff)), stat=ialloc); call err(ialloc)
  do i=1,size(estuff)
    wavevector_states(i) = WavevectorState( this%subspace_id,            &
                                          & this%wavevector,             &
                                          & hamiltonian%selected_states, &
                                          & estuff(i)%evec               )
  enddo
  
  ! N.B. there is not enough information at this stage to calculate weights.
  ! These must be calculated once the states from each wavevector are collated.
  output = BasisStatesPointer(WavevectorStates(       &
     & subspace_id     = this%subspace_id,            &
     & states          = wavevector_states,           &
     & energies        = estuff%eval                  ))
end function

function calculate_hamiltonian_real(basis,harmonic_states,          &
   & state_energy_cutoff,subspace_potential,include_kinetic_energy, &
   & anharmonic_data) result(output) 
  implicit none
  
  type(WavevectorBasis),   intent(in)           :: basis
  type(HarmonicStateReal), intent(in), target   :: harmonic_states(:)
  real(dp),                intent(in)           :: state_energy_cutoff
  class(PotentialBase),    intent(in)           :: subspace_potential
  logical,                 intent(in), optional :: include_kinetic_energy
  type(AnharmonicData),    intent(in)           :: anharmonic_data
  type(SelectedStatesHamiltonian)               :: output
  
  logical :: include_kinetic_energy_
  
  type(HarmonicBraKetReal) :: braket
  
  real(dp), allocatable :: diagonal(:)
  real(dp)              :: min_element
  
  type(IntArray1d), allocatable :: selected_couplings(:)
  
  integer                        :: no_elements
  type(RealElement), allocatable :: elements(:)
  real(dp)                       :: element
  
  integer :: i,j,k,l,ialloc
  
  include_kinetic_energy_ = set_default(include_kinetic_energy, .false.)
  
  ! Initialise the braket. This can be done from any state.
  braket = HarmonicBraKetReal(basis, harmonic_states(1), anharmonic_data)
  
  ! Calculate the diagonal of the Hamiltonian for all states.
  allocate(diagonal(size(harmonic_states)), stat=ialloc); call err(ialloc)
  do i=1,size(harmonic_states)
    braket%bra_ => harmonic_states(i)
    braket%ket_ => harmonic_states(i)
    diagonal(i) = braket%kinetic_energy(anharmonic_data)               &
              & + subspace_potential%potential_energy( braket,         &
              &                                        anharmonic_data )
  enddo
  
  ! Select the states within state_energy_cutoff of the minimum energy.
  min_element = minval(diagonal)
  output%selected_states = filter(diagonal-min_element<state_energy_cutoff)
  selected_couplings = selected_states_couplings( basis%harmonic_couplings_, &
                                                & output%selected_states     )
  
  ! Calculate the Hamiltonian in the basis of selected states.
  no_elements = sum([( size(selected_couplings(i)%i), &
                     & i=1,                           &
                     & size(selected_couplings)       )])
  allocate( elements(no_elements), stat=ialloc); call err(ialloc)
  l = 0
  do i=1,size(selected_couplings)
    braket%bra_ => harmonic_states(output%selected_states(i))
    do j=1,size(selected_couplings(i))
      k = selected_couplings(i)%i(j)
      
      ! Ignore k<i. These are added when k and i are reversed.
      if (k<i) then
        cycle
      endif
      
      braket%ket_ => harmonic_states(output%selected_states(k))
      
      element = braket%kinetic_energy(anharmonic_data)               &
            & + subspace_potential%potential_energy( braket,         &
            &                                        anharmonic_data )
      
      if (abs(element)<1e-300_dp) then
        cycle
      endif
      
      l = l+1
      elements(l) = RealElement(i, k, element)
      if (k/=i) then
        l = l+1
        elements(l) = RealElement(k, i, element)
      endif
    enddo
  enddo
  
  output%hamiltonian = SparseRealMatrix( size(output%selected_states), &
                                       & size(output%selected_states), &
                                       & elements(:l)                  )
end function

function calculate_hamiltonian_complex(basis,harmonic_states,       &
   & state_energy_cutoff,subspace_potential,include_kinetic_energy, &
   & anharmonic_data) result(output) 
  implicit none
  
  type(WavevectorBasis),      intent(in)           :: basis
  type(HarmonicStateComplex), intent(in), target   :: harmonic_states(:)
  real(dp),                   intent(in)           :: state_energy_cutoff
  class(PotentialBase),       intent(in)           :: subspace_potential
  logical,                    intent(in), optional :: include_kinetic_energy
  type(AnharmonicData),       intent(in)           :: anharmonic_data
  type(SelectedStatesHamiltonian)                  :: output
  
  logical :: include_kinetic_energy_
  
  type(HarmonicBraKetComplex) :: braket
  
  real(dp), allocatable :: diagonal(:)
  real(dp)              :: min_element
  
  type(IntArray1d), allocatable :: selected_couplings(:)
  
  integer                        :: no_elements
  type(RealElement), allocatable :: elements(:)
  real(dp)                       :: element
  
  integer :: i,j,k,l,ialloc
  
  include_kinetic_energy_ = set_default(include_kinetic_energy, .false.)
  
  ! Initialise the braket. This can be done from any state.
  braket = HarmonicBraKetComplex(basis, harmonic_states(1), anharmonic_data)
  
  ! Calculate the diagonal of the Hamiltonian for all states.
  allocate(diagonal(size(harmonic_states)), stat=ialloc); call err(ialloc)
  do i=1,size(harmonic_states)
    braket%bra_ => harmonic_states(i)
    braket%ket_ => harmonic_states(i)
    diagonal(i) = braket%kinetic_energy(anharmonic_data)               &
              & + subspace_potential%potential_energy( braket,         &
              &                                        anharmonic_data )
  enddo
  
  ! Select the states within state_energy_cutoff of the minimum energy.
  min_element = minval(diagonal)
  output%selected_states = filter(diagonal-min_element<state_energy_cutoff)
  selected_couplings = selected_states_couplings( basis%harmonic_couplings_, &
                                                & output%selected_states     )
  
  ! Calculate the Hamiltonian in the basis of selected states.
  no_elements = sum([( size(selected_couplings(i)%i), &
                     & i=1,                           &
                     & size(selected_couplings)       )])
  allocate( elements(no_elements), stat=ialloc); call err(ialloc)
  l = 0
  do i=1,size(selected_couplings)
    braket%bra_ => harmonic_states(output%selected_states(i))
    do j=1,size(selected_couplings(i))
      k = selected_couplings(i)%i(j)
      
      ! Ignore k<i. These are added when k and i are reversed.
      if (k<i) then
        cycle
      endif
      
      braket%ket_ => harmonic_states(output%selected_states(k))
      
      element = braket%kinetic_energy(anharmonic_data)               &
            & + subspace_potential%potential_energy( braket,         &
            &                                        anharmonic_data )
      
      if (abs(element)<1e-300_dp) then
        cycle
      endif
      
      l = l+1
      elements(l) = RealElement(i, k, element)
      if (k/=i) then
        l = l+1
        elements(l) = RealElement(k, i, element)
      endif
    enddo
  enddo
  
  output%hamiltonian = SparseRealMatrix( size(output%selected_states), &
                                       & size(output%selected_states), &
                                       & elements(:l)                  )
end function

! Calculate the eigenstates of a wavevector basis.
function calculate_states(basis,subspace,subspace_potential,thermal_energy, &
   & state_energy_cutoff,convergence_data,anharmonic_data) result(output) 
  implicit none
  
  type(WavevectorBasis),    intent(in) :: basis(:)
  type(DegenerateSubspace), intent(in) :: subspace
  class(PotentialBase),     intent(in) :: subspace_potential
  real(dp),                 intent(in) :: thermal_energy
  real(dp),                 intent(in) :: state_energy_cutoff
  type(ConvergenceData),    intent(in) :: convergence_data
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(BasisStatesPointer)             :: output
  
  type(WavevectorStates)             :: wavevector_states
  type(WavevectorState), allocatable :: states(:)
  real(dp),              allocatable :: energies(:)
  real(dp),              allocatable :: weights(:)
  
  integer              :: key
  integer, allocatable :: keys(:,:)
  
  integer :: i,ialloc
  
  allocate( keys(2,size(basis)),    &
          & stat=ialloc); call err(ialloc)
  
  allocate( states(0),   &
          & energies(0), &
          & stat=ialloc); call err(ialloc)
  key = 0
  do i=1,size(basis)
    wavevector_states = WavevectorStates(      &
       & basis(i)%calculate_states(            &
       &            subspace,                  &
       &            subspace_potential,        &
       &            thermal_energy,            &
       &            state_energy_cutoff,       &
       &            convergence_data,          &
       &            anharmonic_data          ) )
    keys(1,i) = key+1
    keys(2,i) = key+size(wavevector_states%states)
    key = keys(2,i)
    states = [states, wavevector_states%states]
    energies = [energies, wavevector_states%energies]
  enddo
  
  weights = calculate_weights(energies,thermal_energy)
  
  if (thermal_energy>0) then
    if (count(weights>1e-100_dp)==1) then
      call print_line('')
      call print_line(WARNING//': State energy differences in subspace '// &
         &subspace%id//' are very large compared to the temperature. This may &
         &block convergence.')
      call print_line('Relative Energies / (kB T):')
      call print_line((energies-minval(energies))/thermal_energy)
    endif
  endif
  
  wavevector_states = WavevectorStates( subspace%id, &
                                      & states,      &
                                      & energies,    &
                                      & weights      )
  
  allocate( wavevector_states%density_matrices(size(basis)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(basis)
    wavevector_states%density_matrices(i) = DensityMatrix( &
                           & basis(i)%wavevector,          &
                           & basis(i)%harmonic_couplings_, &
                           & states(keys(1,i):keys(2,i)),  &
                           & weights(keys(1,i):keys(2,i))  )
  enddo
  
  output = BasisStatesPointer(wavevector_states)
end function

! ----------------------------------------------------------------------
! Operations which calculate thermodynamic data.
! ----------------------------------------------------------------------
! Calculate the harmonic expectation of the harmonic potential,
!    but only using the states in the basis.
! N.B. this normalises the sum of thermal weights of the included states
!    to be one,
!    rather than normalising to the sum across all states to be one.
function core_harmonic_observables_WavevectorBasis(bases,thermal_energy, &
   & stress,stress_prefactors,anharmonic_data) result(output)
  implicit none
  
  type(WavevectorBasis),  intent(in), target   :: bases(:)
  real(dp),               intent(in)           :: thermal_energy
  class(StressBase),      intent(in), optional :: stress
  type(StressPrefactors), intent(in), optional :: stress_prefactors
  type(AnharmonicData),   intent(in)           :: anharmonic_data
  type(ThermodynamicData)                      :: output
  
  type(HarmonicBraKetReal)    :: braket_real
  type(HarmonicBraKetComplex) :: braket_complex
  
  real(dp), allocatable :: kinetic_energies(:)
  
  type(RealMatrix), allocatable :: stresses(:)
  real(dp),         allocatable :: volume
  
  integer :: no_states
  
  integer :: i,j,k,ialloc
  
  if (present(stress) .neqv. present(stress_prefactors)) then
    call print_line(CODE_ERROR//': Either both or neither of stress and &
       &stress_prefactors must be present.')
    call err()
  endif
  
  no_states = sum([(size(bases(i)%harmonic_states_), i=1, size(bases))])
  
  allocate(kinetic_energies(no_states), stat=ialloc); call err(ialloc)
  if (present(stress)) then
    allocate(stresses(no_states), stat=ialloc); call err(ialloc)
    volume = anharmonic_data%structure%volume
  endif
  
  k = 0
  do i=1,size(bases)
    associate(harmonic_states=>bases(i)%harmonic_states_)
      select type(harmonic_states); type is(HarmonicStateReal)
        braket_real = HarmonicBraKetReal( bases(i),              &
                                        & harmonic_states(1), &
                                        & anharmonic_data     )
        do j=1,size(harmonic_states)
          braket_real%bra_ => harmonic_states(j)
          braket_real%ket_ => harmonic_states(j)
          kinetic_energies(j+k) = braket_real%kinetic_energy(anharmonic_data)
          if (present(stress)) then
            stresses(j+k) = braket_real%kinetic_stress( stress_prefactors,   &
                        &                               anharmonic_data    ) &
                        & + stress%potential_stress( braket_real,            &
                        &                            anharmonic_data )
          endif
        enddo
      type is(HarmonicStateComplex)
        braket_complex = HarmonicBraKetComplex( bases(i),           &
                                              & harmonic_states(1), &
                                              & anharmonic_data     )
        do j=1,size(harmonic_states)
          braket_complex%bra_ => harmonic_states(j)
          braket_complex%ket_ => harmonic_states(j)
          kinetic_energies(j+k) = braket_complex%kinetic_energy( &
                                               & anharmonic_data )
          if (present(stress)) then
            stresses(j+k) = braket_complex%kinetic_stress(           &
                        &               stress_prefactors,           &
                        &               anharmonic_data    )         &
                        & + stress%potential_stress( braket_complex, &
                        &                            anharmonic_data )
          endif
        enddo
      end select
      k = k+size(harmonic_states)
    end associate
  enddo
  
  ! N.B. U = <|T+V|> = 2<|T|> for the harmonic potential.
  output = ThermodynamicData( thermal_energy,     &
                            & 2*kinetic_energies, &
                            & stresses,           &
                            & volume              )
end function

! Calculate the harmonic expectation of a potential,
!    but only using the states in the basis.
! N.B. this normalises the sum of thermal weights of the included states
!    to be one,
!    rather than normalising to the sum across all states to be one.
function core_effective_harmonic_observables_WavevectorBasis(bases,     &
   & thermal_energy,potential,stress,stress_prefactors,anharmonic_data) &
   & result(output)
  implicit none
  
  type(WavevectorBasis),  intent(in), target   :: bases(:)
  real(dp),               intent(in)           :: thermal_energy
  class(PotentialBase),   intent(in)           :: potential
  class(StressBase),      intent(in), optional :: stress
  type(StressPrefactors), intent(in), optional :: stress_prefactors
  type(AnharmonicData),   intent(in)           :: anharmonic_data
  type(ThermodynamicData)                      :: output
  
  integer :: no_states
  
  type(HarmonicBraKetReal)    :: braket_real
  type(HarmonicBraKetComplex) :: braket_complex
  
  real(dp), allocatable :: harmonic_kinetic_energies(:)
  real(dp), allocatable :: harmonic_energy_differences(:)
  real(dp), allocatable :: anharmonic_potential_energies(:)
  real(dp), allocatable :: boltzmann_factors(:)
  
  integer :: ground_state
  
  real(dp) :: anharmonic_minus_harmonic
  
  integer :: i,j,k,ialloc
  
  ! Calculate properties with a harmonic potential.
  output = core_harmonic_observables( bases,             &
                                    & thermal_energy,    &
                                    & stress,            &
                                    & stress_prefactors, &
                                    & anharmonic_data    )
  
  no_states = sum([(size(bases(i)%harmonic_states_), i=1, size(bases))])
  
  ! Calculate <|T|> for the harmonic potential,
  !    and <|V|> for the input potential.
  allocate( harmonic_kinetic_energies(no_states),     &
          & anharmonic_potential_energies(no_states), &
          & stat=ialloc); call err(ialloc)
  k = 0
  do i=1,size(bases)
    associate(harmonic_states=>bases(i)%harmonic_states_)
      select type(harmonic_states); type is(HarmonicStateReal)
        braket_real = HarmonicBraKetReal( bases(i),           &
                                        & harmonic_states(1), &
                                        & anharmonic_data     )
        do j=1,size(harmonic_states)
          braket_real%bra_ => harmonic_states(j)
          braket_real%ket_ => harmonic_states(j)
          harmonic_kinetic_energies(j+k) = braket_real%kinetic_energy( &
                                                     & anharmonic_data )
          anharmonic_potential_energies(j+k) = potential%potential_energy( &
                                                         & braket_real,    &
                                                         & anharmonic_data )
        enddo
      type is(HarmonicStateComplex)
        braket_complex = HarmonicBraKetComplex( bases(i),           &
                                              & harmonic_states(1), &
                                              & anharmonic_data     )
        do j=1,size(harmonic_states)
          braket_complex%bra_ => harmonic_states(j)
          braket_complex%ket_ => harmonic_states(j)
          harmonic_kinetic_energies(j+k) = braket_complex%kinetic_energy( &
                                                        & anharmonic_data )
          anharmonic_potential_energies(j+k) = potential%potential_energy( &
                                                         & braket_complex, &
                                                         & anharmonic_data )
        enddo
      end select
      k = k + size(harmonic_states)
    end associate
  enddo
  
  ! Identify the ground state.
  ground_state = minloc(harmonic_kinetic_energies, 1)
  
  ! Calculate <|T+V|> - min(<|T+V|>) for the harmonic potential.
  ! N.B. <|V|> = <|T|> for the harmonic potential.
  harmonic_energy_differences = 2*( harmonic_kinetic_energies               &
                                & - harmonic_kinetic_energies(ground_state) )
  
  ! Calculate the difference in <|V|> between the harmonic
  !    and anharmonic potentials.
  if ( maxval(harmonic_energy_differences) < 1e20_dp*thermal_energy) then
    ! Normal temperature range.
    boltzmann_factors = [exp( -harmonic_energy_differences &
                          & / thermal_energy               )]
    anharmonic_minus_harmonic = sum( boltzmann_factors                   &
                            &      * ( anharmonic_potential_energies     &
                            &        - harmonic_kinetic_energies     ) ) &
                            & / sum(boltzmann_factors)
  else
    ! Very low temperature limit.
    anharmonic_minus_harmonic = anharmonic_potential_energies(ground_state) &
                            & - harmonic_kinetic_energies(ground_state)
  endif
  
  ! Add the change in <|V|> to U, F, H and G.
  output%energy = output%energy + anharmonic_minus_harmonic
  output%free_energy = output%free_energy + anharmonic_minus_harmonic
  if (allocated(output%enthalpy)) then
    output%enthalpy = output%enthalpy + anharmonic_minus_harmonic
    output%gibbs = output%gibbs + anharmonic_minus_harmonic
  endif
end function

! Calculate thermodynamic data for the core VCI states.
function core_vci_observables_WavevectorBasis(bases,thermal_energy,states, &
   & subspace,subspace_potential,subspace_stress,stress_prefactors,        &
   & anharmonic_data) result(output)
  implicit none
  
  type(WavevectorBasis),    intent(in)                  :: bases(:)
  real(dp),                 intent(in)                  :: thermal_energy
  class(BasisStates),       intent(in),          target :: states
  type(DegenerateSubspace), intent(in)                  :: subspace
  class(PotentialBase),     intent(in)                  :: subspace_potential
  class(StressBase),        intent(in), optional        :: subspace_stress
  type(StressPrefactors),   intent(in), optional        :: stress_prefactors
  type(AnharmonicData),     intent(in)                  :: anharmonic_data
  type(ThermodynamicData)                               :: output
  
  type(WavevectorStates) :: wavevector_states
  
  type(RealMatrix) :: stress
  real(dp)         :: volume
  
  real(dp) :: internal_energy
  real(dp) :: energy_difference
  
  integer :: i,j
  
  if (present(subspace_stress) .neqv. present(stress_prefactors)) then
    call print_line(CODE_ERROR//': Only one of subspace_stress and &
       &stress_prefactors passed.')
    call err()
  endif
  
  wavevector_states = WavevectorStates(states)
  
  output = ThermodynamicData( thermal_energy,            &
                            & wavevector_states%energies )
  
  ! Calculate stress.
  internal_energy = 0.0_dp
  if (present(subspace_stress)) then
    stress = dblemat(zeroes(3,3))
  endif
  
  do i=1,size(bases)
    j = first(    wavevector_states%density_matrices%wavevector &
             & == bases(i)%wavevector,                          &
             & default=0                                        )
    if (j==0) then
      cycle
    elseif (size(wavevector_states%density_matrices(j)%values)==0) then
      cycle
    endif
    
    associate(harmonic_states=>bases(i)%harmonic_states_)
      select type(harmonic_states); type is(HarmonicStateReal)
        internal_energy = internal_energy                           &
                      & + internal_energy_real(                     &
                      &      bases(i),                              &
                      &      harmonic_states,                       &
                      &      wavevector_states%density_matrices(j), &
                      &      subspace_potential,                    &
                      &      anharmonic_data                        )
        if (present(subspace_stress)) then
          stress = stress                                              &
               & + stress_real( bases(i),                              &
               &                harmonic_states,                       &
               &                wavevector_states%density_matrices(j), &
               &                subspace_stress,                       &
               &                stress_prefactors,                     &
               &                anharmonic_data                        )
        endif
      type is(HarmonicStateComplex)
        internal_energy = internal_energy                           &
                      & + internal_energy_complex(                  &
                      &      bases(i),                              &
                      &      harmonic_states,                       &
                      &      wavevector_states%density_matrices(j), &
                      &      subspace_potential,                    &
                      &      anharmonic_data                        )
        if (present(subspace_stress)) then
          stress = stress                                                 &
               & + stress_complex( bases(i),                              &
               &                   harmonic_states,                       &
               &                   wavevector_states%density_matrices(j), &
               &                   subspace_stress,                       &
               &                   stress_prefactors,                     &
               &                   anharmonic_data                        )
        endif
      end select
    end associate
  enddo
  
  energy_difference = internal_energy-output%energy
  output%energy = output%energy+energy_difference
  output%free_energy = output%free_energy+energy_difference
  
  if (present(subspace_stress)) then
    volume = anharmonic_data%structure%volume
    call output%set_stress(stress,volume)
  endif
end function

function internal_energy_real(basis,harmonic_states,density_matrix,potential, &
   & anharmonic_data) result(output) 
  implicit none
  
  type(WavevectorBasis),   intent(in)         :: basis
  type(HarmonicStateReal), intent(in), target :: harmonic_states(:)
  type(DensityMatrix),     intent(in)         :: density_matrix
  class(PotentialBase),    intent(in)         :: potential
  type(AnharmonicData),    intent(in)         :: anharmonic_data
  real(dp)                                    :: output
  
  type(HarmonicBraKetReal) :: braket
  
  real(dp) :: term
  
  integer :: i
  
  ! Initialise the braket. This can be done from any state.
  braket = HarmonicBraKetReal(basis, harmonic_states(1), anharmonic_data)
  
  output = 0.0_dp
  do i=1,size(density_matrix%values)
    ! Ignore terms with ket<bra. These are added when i and k are reversed.
    if (density_matrix%ket_ids(i)<density_matrix%bra_ids(i)) then
      cycle
    endif
    
    ! Calculate <bra|X|ket> in the harmonic basis,
    !    and thermally weight by the density matrix.
    braket%bra_ => harmonic_states(density_matrix%bra_ids(i))
    braket%ket_ => harmonic_states(density_matrix%ket_ids(i))
    term = ( potential%potential_energy(braket, anharmonic_data)   &
       &   + braket%kinetic_energy(anharmonic_data)              ) &
       & * density_matrix%values(i)
    
    ! If bra==ket, only include <bra|X|bra>,
    !    otherwise include <bra|X|ket> and <ket|X|bra>.
    if (density_matrix%bra_ids(i)==density_matrix%ket_ids(i)) then
      output = output + term
    else
      output = output + 2*term
    endif
  enddo
end function

function internal_energy_complex(basis,harmonic_states,density_matrix, &
   & potential,anharmonic_data) result(output) 
  implicit none
  
  type(WavevectorBasis),      intent(in)         :: basis
  type(HarmonicStateComplex), intent(in), target :: harmonic_states(:)
  type(DensityMatrix),        intent(in)         :: density_matrix
  class(PotentialBase),       intent(in)         :: potential
  type(AnharmonicData),       intent(in)         :: anharmonic_data
  real(dp)                                       :: output
  
  type(HarmonicBraKetComplex) :: braket
  
  real(dp) :: term
  
  integer :: i
  
  ! Initialise the braket. This can be done from any state.
  braket = HarmonicBraKetComplex(basis, harmonic_states(1), anharmonic_data)
  
  output = 0.0_dp
  do i=1,size(density_matrix%values)
    ! Ignore terms with ket<bra. These are added when i and k are reversed.
    if (density_matrix%ket_ids(i)<density_matrix%bra_ids(i)) then
      cycle
    endif
    
    ! Calculate <bra|X|ket> in the harmonic basis,
    !    and thermally weight by the density matrix.
    braket%bra_ => harmonic_states(density_matrix%bra_ids(i))
    braket%ket_ => harmonic_states(density_matrix%ket_ids(i))
    term = ( potential%potential_energy(braket, anharmonic_data)   &
       &   + braket%kinetic_energy(anharmonic_data)              ) &
       & * density_matrix%values(i)
    
    ! If bra==ket, only include <bra|X|bra>,
    !    otherwise include <bra|X|ket> and <ket|X|bra>.
    if (density_matrix%bra_ids(i)==density_matrix%ket_ids(i)) then
      output = output + term
    else
      output = output + 2*term
    endif
  enddo
end function

function stress_real(basis,harmonic_states,density_matrix,stress, &
   & stress_prefactors,anharmonic_data) result(output) 
  implicit none
  
  type(WavevectorBasis),   intent(in)         :: basis
  type(HarmonicStateReal), intent(in), target :: harmonic_states(:)
  type(DensityMatrix),     intent(in)         :: density_matrix
  class(StressBase),       intent(in)         :: stress
  type(StressPrefactors),  intent(in)         :: stress_prefactors
  type(AnharmonicData),    intent(in)         :: anharmonic_data
  type(RealMatrix)                            :: output
  
  type(HarmonicBraKetReal) :: braket
  
  type(RealMatrix) :: term
  
  integer :: i
  
  ! Initialise the braket. This can be done from any state.
  braket = HarmonicBraKetReal(basis, harmonic_states(1), anharmonic_data)
  
  output = dblemat(zeroes(3,3))
  do i=1,size(density_matrix%values)
    ! Ignore terms with ket<bra. These are added when i and k are reversed.
    if (density_matrix%ket_ids(i)<density_matrix%bra_ids(i)) then
      cycle
    endif
    
    ! Calculate <bra|X|ket> in the harmonic basis,
    !    and thermally weight by the density matrix.
    braket%bra_ => harmonic_states(density_matrix%bra_ids(i))
    braket%ket_ => harmonic_states(density_matrix%ket_ids(i))
    term = ( stress%potential_stress(braket, anharmonic_data)            &
       &   + braket%kinetic_stress(stress_prefactors, anharmonic_data) ) &
       & * density_matrix%values(i)
    
    ! If bra==ket, only include <bra|X|bra>,
    !    otherwise include <bra|X|ket> and <ket|X|bra>.
    if (density_matrix%bra_ids(i)==density_matrix%ket_ids(i)) then
      output = output + term
    else
      output = output + 2*term
    endif
  enddo
end function

function stress_complex(basis,harmonic_states,density_matrix,stress, &
   & stress_prefactors,anharmonic_data) result(output) 
  implicit none
  
  type(WavevectorBasis),      intent(in)         :: basis
  type(HarmonicStateComplex), intent(in), target :: harmonic_states(:)
  type(DensityMatrix),        intent(in)         :: density_matrix
  class(StressBase),          intent(in)         :: stress
  type(StressPrefactors),     intent(in)         :: stress_prefactors
  type(AnharmonicData),       intent(in)         :: anharmonic_data
  type(RealMatrix)                               :: output
  
  type(HarmonicBraKetComplex) :: braket
  
  type(RealMatrix) :: term
  
  integer :: i
  
  ! Initialise the braket. This can be done from any state.
  braket = HarmonicBraKetComplex(basis, harmonic_states(1), anharmonic_data)
  
  output = dblemat(zeroes(3,3))
  do i=1,size(density_matrix%values)
    ! Ignore terms with ket<bra. These are added when i and k are reversed.
    if (density_matrix%ket_ids(i)<density_matrix%bra_ids(i)) then
      cycle
    endif
    
    ! Calculate <bra|X|ket> in the harmonic basis,
    !    and thermally weight by the density matrix.
    braket%bra_ => harmonic_states(density_matrix%bra_ids(i))
    braket%ket_ => harmonic_states(density_matrix%ket_ids(i))
    term = ( stress%potential_stress(braket, anharmonic_data)            &
       &   + braket%kinetic_stress(stress_prefactors, anharmonic_data) ) &
       & * density_matrix%values(i)
    
    ! If bra==ket, only include <bra|X|bra>,
    !    otherwise include <bra|X|ket> and <ket|X|bra>.
    if (density_matrix%bra_ids(i)==density_matrix%ket_ids(i)) then
      output = output + term
    else
      output = output + 2*term
    endif
  enddo
end function

! ----------------------------------------------------------------------
! Mode IDs.
! ----------------------------------------------------------------------
function mode_ids_WavevectorBasis(this,subspace,anharmonic_data) result(output)
  implicit none
  
  class(WavevectorBasis),   intent(in) :: this
  type(DegenerateSubspace), intent(in) :: subspace
  type(AnharmonicData),     intent(in) :: anharmonic_data
  integer, allocatable                 :: output(:)
  
  output = this%harmonic_states_(1)%mode_ids()
end function

function paired_mode_ids_WavevectorBasis(this,subspace,anharmonic_data) &
   & result(output)
  implicit none
  
  class(WavevectorBasis),   intent(in) :: this
  type(DegenerateSubspace), intent(in) :: subspace
  type(AnharmonicData),     intent(in) :: anharmonic_data
  integer, allocatable                 :: output(:)
  
  output = this%harmonic_states_(1)%paired_mode_ids()
end function

! Calculate the derivative of the free energy.
function free_energy_gradient_WavevectorBasis(this,subspace_potential,   &
   & basis_functions,subspace,states,thermal_energy,state_energy_cutoff, &
   & anharmonic_data) result(output) 
  implicit none
  
  class(WavevectorBasis),   intent(in) :: this
  class(PotentialBase),     intent(in) :: subspace_potential
  class(PotentialBase),     intent(in) :: basis_functions(:)
  type(DegenerateSubspace), intent(in) :: subspace
  class(BasisStates),       intent(in) :: states
  real(dp),                 intent(in) :: thermal_energy
  real(dp),                 intent(in) :: state_energy_cutoff
  type(AnharmonicData),     intent(in) :: anharmonic_data
  real(dp), allocatable                :: output(:)
  
  ! Calculating the free energy gradient is not simply additive across the
  !    wavevectors, so free_energy_gradient_WavevectorBases should be called
  !    instead.
  call err()
end function

function free_energy_gradient_WavevectorBases(bases,subspace_potential,  &
   & basis_functions,subspace,states,thermal_energy,state_energy_cutoff, &
   & anharmonic_data) result(output) 
  implicit none
  
  type(WavevectorBasis),    intent(in) :: bases(:)
  class(PotentialBase),     intent(in) :: subspace_potential
  class(PotentialBase),     intent(in) :: basis_functions(:)
  type(DegenerateSubspace), intent(in) :: subspace
  class(BasisStates),       intent(in) :: states
  real(dp),                 intent(in) :: thermal_energy
  real(dp),                 intent(in) :: state_energy_cutoff
  type(AnharmonicData),     intent(in) :: anharmonic_data
  real(dp), allocatable                :: output(:)
  
  type(WavevectorStates) :: states_
  
  logical,             allocatable :: wavevector_included(:)
  type(DensityMatrix), allocatable :: density_matrices(:)
  
  type(IntArray1D),    allocatable :: states_at_wavevector(:)
  
  type(WavevectorState), allocatable :: eigenstates(:)
  real(dp),              allocatable :: state_energies(:)
  real(dp),              allocatable :: state_weights(:)
  
  real(dp),  allocatable :: transformation(:,:)
  
  type(HarmonicStateReal),         allocatable :: real_states(:)
  type(HarmonicStateComplex),      allocatable :: complex_states(:)
  
  type(SelectedStatesHamiltonian), allocatable :: subspace_hamiltonians_(:)
  type(SelectedStatesHamiltonian), allocatable :: basis_hamiltonians_(:,:)
  
  type(RealMatrix), allocatable :: subspace_hamiltonians(:)
  type(RealMatrix), allocatable :: basis_hamiltonians(:,:)
  
  real(dp), allocatable :: basis_expectations(:)
  
  real(dp) :: dmixing_dc
  real(dp) :: dp_dc
  
  integer :: i,j,k,l,ialloc
  
  output = [(0.0_dp, i=1, size(basis_functions))]
  
  ! The gradient cannot be calculated at zero temperature.
  if (thermal_energy<1e-300_dp) then
    return
  endif
  
  states_ = WavevectorStates(states)
  
  allocate( wavevector_included(size(bases)),              &
          & density_matrices(size(bases)),                 &
          & states_at_wavevector(size(bases)),             &
          & subspace_hamiltonians_(size(bases)),           &
          & subspace_hamiltonians(size(bases)),            &
          & basis_hamiltonians_( size(basis_functions),    &
          &                      size(bases)            ), &
          & basis_hamiltonians( size(basis_functions),     &
          &                     size(bases)            ),  &
          & stat=ialloc); call err(ialloc)
  wavevector_included = .true.
  basis_expectations = [(0.0_dp, i=1, size(basis_functions))]
  do i=1,size(bases)
    j = first( states_%density_matrices%wavevector==bases(i)%wavevector, &
             & default=0                                                 )
    if (j==0) then
      wavevector_included(i) = .false.
      cycle
    elseif (.not. allocated(states_%density_matrices(j)%values)) then
      wavevector_included(i) = .false.
      cycle
    elseif (thermal_energy<1e-300_dp) then
      wavevector_included(i) = .false.
      cycle
    endif
    
    density_matrices(i) = states_%density_matrices(j)
    
    states_at_wavevector(i) = array(filter(             &
       & states_%states%wavevector==bases(i)%wavevector &
       & .and. states_%weights>1e-300_dp                ))
    if (size(states_at_wavevector(i)%i)==0) then
      wavevector_included(i) = .false.
      cycle
    endif
    
    eigenstates = states_%states(states_at_wavevector(i)%i)
    state_energies = states_%energies(states_at_wavevector(i)%i)
    state_weights = states_%weights(states_at_wavevector(i)%i)
    
    allocate( transformation( size(eigenstates),                   &
            &                 size(eigenstates(1)%coefficients) ), &
            & stat=ialloc); call err(ialloc)
    do j=1,size(eigenstates)
      transformation(j,:) = eigenstates(j)%coefficients
    enddo
    
    ! Calculate the Hamiltonian in the harmonic basis for subspace_potential
    !    and each basis_function.
    associate(harmonic_states=>bases(i)%harmonic_states_)
      select type(harmonic_states); type is(HarmonicStateReal)
        real_states = harmonic_states(eigenstates(1)%state_ids)
        subspace_hamiltonians_(i) = calculate_hamiltonian_real( &
                            & bases(i),                         &
                            & real_states,                      &
                            & huge(1.0_dp),                     &
                            & subspace_potential,               &
                            & anharmonic_data = anharmonic_data )
        do j=1,size(basis_functions)
          basis_hamiltonians_(j,i) = calculate_hamiltonian_real( &
                             & bases(i),                         &
                             & real_states,                      &
                             & huge(1.0_dp),                     &
                             & basis_functions(j),               &
                             & include_kinetic_energy = .false., &
                             & anharmonic_data = anharmonic_data )
        enddo
      type is(HarmonicStateComplex)
        complex_states = harmonic_states(eigenstates(1)%state_ids)
        subspace_hamiltonians_(i) = calculate_hamiltonian_complex( &
                               & bases(i),                         &
                               & complex_states,                   &
                               & huge(1.0_dp),                     &
                               & subspace_potential,               &
                               & anharmonic_data = anharmonic_data )
        do j=1,size(basis_functions)
          basis_hamiltonians_(j,i) = calculate_hamiltonian_complex( &
                                & bases(i),                         &
                                & complex_states,                   &
                                & huge(1.0_dp),                     &
                                & basis_functions(j),               &
                                & include_kinetic_energy = .false., &
                                & anharmonic_data = anharmonic_data )
        enddo
      end select
    end associate
    
    ! Transform the Hamiltonians into the eigenbasis.
    subspace_hamiltonians(i) = mat(transformation)                   &
                           & * subspace_hamiltonians_(i)%hamiltonian &
                           & * mat(transpose(transformation))
    do j=1,size(basis_functions)
      basis_hamiltonians(j,i) = mat(transformation)                  &
                            & * basis_hamiltonians_(j,i)%hamiltonian &
                            & * mat(transpose(transformation))
    enddo
    
    deallocate(transformation, stat=ialloc); call err(ialloc)
    
    ! Calculate the contribution to <V_j> from this wavevector.
    do j=1,size(basis_functions)
      basis_expectations(j) = basis_expectations(j)           &
         & + sum( state_weights                               &
         &      * [( basis_hamiltonians(j,i)%element(k,k),    &
         &           k=1,                                     &
         &           size(eigenstates)                     )] )
    enddo
  enddo
  
  if (.not. any(wavevector_included)) then
    return
  endif
  
  ! The subspace Hamiltonian is H.
  ! The auxilliary Hamiltonian is sum_j c_j V_j + constant terms,
  !    where each c_j is a coefficient and each V_j is a basis function.
  ! The free energy of the potential V with the states from the auxilliary
  !    Hamiltonian is F.
  ! dF/d(c_j) =
  ! sum_k P_k[ (<V_j>-<k|V_j|k>) * (<k|H|k>+T*ln(P_k)) / T
  !          + (sum_l/=k <k|H|l><l|V_j|k>/(E_k-E_l) + c.c.) ]
  do i=1,size(bases)
    if (.not. wavevector_included(i)) then
      cycle
    endif
    
    eigenstates = states_%states(states_at_wavevector(i)%i)
    state_energies = states_%energies(states_at_wavevector(i)%i)
    state_weights = states_%weights(states_at_wavevector(i)%i)
    
    do j=1,size(basis_functions)
      do k=1,size(eigenstates)
        ! Calculate sum_l P_k <k|H|l><l|V_j|k>/(E_k-E_l) + c.c.
        ! N.B. In the eigenbasis, these operators are real, so (+ c.c.) just
        !    introduces a factor of 2.
        do l=1,k-1
          ! Calculate P_k*<l|V_j|k>/(E_k-E_l) + P_l*<k|V_j|l>/(E_l-E_k)
          !         = <l|V_j|k>*(P_k-P_l)/(E_k-E_l).
          ! This is the mixing of state |k> with state |l> under an
          !    infinitesimal change in c_j, equal to [d(<k|)/d(c_i) |l>].
          if (   abs(state_energies(l)-state_energies(k)) &
             & > thermal_energy*1e-10_dp                  ) then
            ! |E_k-E_l|/T is large. Calculate (P_k-P_l)/(E_k-E_l) directly.
            dmixing_dc = 2*basis_hamiltonians(j,i)%element(l,k) &
                     & * (state_weights(k)-state_weights(l))    &
                     & / (state_energies(k)-state_energies(l))
          elseif (   abs(state_energies(k)-state_energies(l)) &
                 & > thermal_energy*1e-20_dp                  ) then
            ! |E_k-E_l|/T < 1e-20, so O((T/|E_k-E_l|)^2) can be neglected.
            ! P_k-P_l = P_k * ( 1-e^((E_k-E_l)/T) )
            !         = P_k ( -(E_k-E_l)/T - ((E_k-E_l)/T)^2/2 + higher orders)
            ! 2*(P_k-P_l)/(E_k-E_l) = P_k
            !                       * (-1/T - (E_k-E_l)/(2T^2) + higher orders)
            !                       + P_l
            !                       * (-1/T - (E_l-E_k)/(2T^2) + higher orders)
            dmixing_dc = -basis_hamiltonians(j,i)%element(l,k)   &
                     & * ( (state_weights(k)+state_weights(l))   &
                     &   / thermal_energy                        &
                     &   + (state_weights(k)-state_weights(l))   &
                     &   * (state_energies(k)-state_energies(l)) &
                     &   / (2*thermal_energy**2)                 )
          else
            ! |E_k-E_l|/T < 1e-20, so O(T/|E_k-E_l|) can be neglected.
            ! P_k-P_l = P_k * ( 1-e^((E_k-E_l)/T) )
            !         = P_k ( (E_l-E_k)/T + higher orders )
            ! 2*(P_k-P_l)/(E_k-E_l) = -(P_k+P_l)/T + higher orders
            dmixing_dc = -basis_hamiltonians(j,i)%element(l,k) &
                     & * (state_weights(k)+state_weights(l)) / thermal_energy
          endif
          
          output(j) = output(j) &
                  & + dmixing_dc * subspace_hamiltonians(i)%element(k,l)
        enddo
        
        ! Calculate d(P_k)/d(c_j) * dF/d(P_k)
        !         = P_k*(<V_j>-<k|V_j|k>)/T * (<k|H|k>+T*ln(P_k)).
        dp_dc = ( basis_expectations(j)                  &
            &   - basis_hamiltonians(j,i)%element(k,k) ) &
            & * state_weights(k)/thermal_energy
        output(j) = output(j)                          &
           & + dp_dc                                   &
           & * ( subspace_hamiltonians(i)%element(k,k) &
           &   + thermal_energy*log(state_weights(k))  )
      enddo
    enddo
  enddo
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_WavevectorBasis(this,input)
  implicit none
  
  class(WavevectorBasis), intent(out) :: this
  type(String),           intent(in)  :: input(:)
  
  integer                                 :: maximum_power
  integer                                 :: expansion_order
  integer                                 :: subspace_id
  real(dp)                                :: frequency
  type(FractionVector)                    :: wavevector
  type(SubspaceStatePointer), allocatable :: harmonic_states(:)
  class(SubspaceState),       allocatable :: states(:)
  type(CoupledStates),        allocatable :: couplings(:)
  
  type(String), allocatable :: line(:)
  
  integer :: no_states
  
  integer :: i,ialloc
  
  ! TODO
  
  select type(this); type is(WavevectorBasis)
    line = split_line(input(1))
    maximum_power = int(line(4))
    
    line = split_line(input(2))
    expansion_order = int(line(4))
    
    line = split_line(input(3))
    subspace_id = int(line(3))
    
    line = split_line(input(4))
    frequency = dble(line(3))
    
    line = split_line(input(5))
    wavevector = FractionVector(join(line(3:5)))
    
    no_states = (size(input)-10)/5
    
    allocate( harmonic_states(no_states), &
            & couplings(no_states),       &
            & stat=ialloc); call err(ialloc)
    do i=1,no_states
      line = split_line(input(6+no_states+i))
      harmonic_states(i) = SubspaceStatePointer([ input(3:4),   &
                                                & str('State'), &
                                                & line(3)       ])
      
      couplings(i) = CoupledStates(input(7+2*no_states+i))
    enddo
    
    states = [( harmonic_states(i)%state(), &
              & i=1,                        &
              & size(harmonic_states)       )]
    
    this = WavevectorBasis( maximum_power,   &
                          & expansion_order, &
                          & subspace_id,     &
                          & frequency,       &
                          & wavevector,      &
                          & states,          &
                          & couplings        )
  class default
    call err()
  end select
end subroutine

function write_WavevectorBasis(this) result(output)
  implicit none
  
  class(WavevectorBasis), intent(in) :: this
  type(String), allocatable          :: output(:)
  
  type(String), allocatable :: state_strings(:)
  
  integer :: i
  
  select type(this); type is(WavevectorBasis)
    output = [ 'Maximum power   : '//this%maximum_power,   &
             & 'Expansion order : '//this%expansion_order, &
             & 'Subspace        : '//this%subspace_id,     &
             & 'Frequency       : '//this%frequency,       &
             & 'Wavevector      : '//this%wavevector,      &
             & str('Harmonic states')                      ]
    do i=1,size(this%harmonic_states_)
      state_strings = str(SubspaceStatePointer(this%harmonic_states_(i)))
      output = [output, '|'//i//'> = '//state_strings(size(state_strings))]
    enddo
    output = [ output,                                              &
             & str('Non-zero <i|H|j> integrals in harmonic basis'), &
             & str(this%harmonic_couplings_)                        ]
  class default
    call err()
  end select
end function

function new_WavevectorBasis_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(WavevectorBasis)    :: this
  
  call this%read(input)
end function

impure elemental function new_WavevectorBasis_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(WavevectorBasis)         :: this
  
  this = WavevectorBasis(str(input))
end function

! ----------------------------------------------------------------------
! Inherited procedures which do not apply to this type.
! ----------------------------------------------------------------------
impure elemental function thermodynamic_data_WavevectorBasis(this,      &
   & thermal_energy,states,subspace,subspace_potential,subspace_stress, &
   & stress_prefactors,anharmonic_data) result(output)
  implicit none
  
  class(WavevectorBasis),   intent(in)                  :: this
  real(dp),                 intent(in)                  :: thermal_energy
  class(BasisStates),       intent(in),          target :: states
  type(DegenerateSubspace), intent(in)                  :: subspace
  class(PotentialBase),     intent(in)                  :: subspace_potential
  class(StressBase),        intent(in), optional        :: subspace_stress
  type(StressPrefactors),   intent(in), optional        :: stress_prefactors
  type(AnharmonicData),     intent(in)                  :: anharmonic_data
  type(ThermodynamicData)                               :: output
  
  call err()
end function

impure elemental function wavefunctions_WavevectorBasis(this,states,subspace, &
   & anharmonic_data) result(output)
  implicit none
  
  class(WavevectorBasis),   intent(in)         :: this
  class(BasisStates),       intent(in), target :: states
  type(DegenerateSubspace), intent(in)         :: subspace
  type(AnharmonicData),     intent(in)         :: anharmonic_data
  type(SubspaceWavefunctionsPointer)           :: output
  
  call err()
end function
end module
