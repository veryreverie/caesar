! ======================================================================
! A basis harmonic states for a given subspace at a given wavevector.
! N.B. the wavevector in this context is the wavevector of the state,
!    not the wavevector of the mode.
! e.g. the state |p> along a mode at q-point q has a wavevector of p*q.
! ======================================================================
module wavevector_basis_module
  use common_module
  
  use anharmonic_common_module
  
  use harmonic_state_1d_module
  use harmonic_state_2d_module
  use harmonic_state_real_module
  use harmonic_state_complex_module
  use coupled_states_module
  use wavevector_state_module
  use wavevector_states_module
  use anharmonic_data_module
  implicit none
  
  private
  
  public :: startup_wavevector_basis
  
  public :: WavevectorBasis
  
  public :: core_harmonic_observables
  public :: core_effective_harmonic_observables
  public :: core_vci_observables
  public :: integrate
  
  type, extends(SubspaceBasis) :: WavevectorBasis
    integer                                          :: maximum_power
    integer                                          :: expansion_order
    integer                                          :: subspace_id
    type(FractionVector)                             :: wavevector
    type(SubspaceStatePointer), allocatable, private :: harmonic_states_(:)
    type(CoupledStates),        allocatable, private :: harmonic_couplings_(:)
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
  
  interface integrate
    module procedure integrate_SparseMonomial_WavevectorBases
  end interface
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
  
  integer,                    intent(in) :: maximum_power
  integer,                    intent(in) :: expansion_order
  integer,                    intent(in) :: subspace_id
  real(dp),                   intent(in) :: frequency
  type(FractionVector),       intent(in) :: wavevector
  type(SubspaceStatePointer), intent(in) :: harmonic_states(:)
  type(CoupledStates),        intent(in) :: harmonic_couplings(:)
  type(WavevectorBasis)                  :: this
  
  if (size(harmonic_couplings)/=size(harmonic_states)) then
    call print_line(CODE_ERROR//': harmonic states and harmonic couplings do &
       &not match.')
  endif
  
  this%maximum_power       = maximum_power
  this%expansion_order     = expansion_order
  this%subspace_id         = subspace_id
  this%frequency           = frequency
  this%wavevector          = wavevector
  this%harmonic_states_    = harmonic_states
  this%harmonic_couplings_ = harmonic_couplings
end function

! ----------------------------------------------------------------------
! Generates states in a given subspace, up to a given power.
! ----------------------------------------------------------------------
! If qpoint is specified, only modes at that q-point are included.
function new_WavevectorBasis_subspace(subspace,frequency,modes,qpoints, &
   & maximum_power,potential_expansion_order,symmetries,qpoint) result(output)
  implicit none
  
  type(DegenerateSubspace), intent(in)           :: subspace
  real(dp),                 intent(in)           :: frequency
  type(ComplexMode),        intent(in)           :: modes(:)
  type(QpointData),         intent(in)           :: qpoints(:)
  integer,                  intent(in)           :: maximum_power
  integer,                  intent(in)           :: potential_expansion_order
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
    mode_basis = generate_mode_basis( subspace,                 &
                                    & frequency,                &
                                    & subspace_modes(i),        &
                                    & maximum_power,            &
                                    & potential_expansion_order )
    
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
   & potential_expansion_order) result(output)
  implicit none
  
  type(DegenerateSubspace), intent(in) :: subspace
  real(dp),                 intent(in) :: frequency
  type(ComplexMode),        intent(in) :: mode
  integer,                  intent(in) :: maximum_power
  integer,                  intent(in) :: potential_expansion_order
  type(WavevectorBasis)                :: output
  
  if (mode%paired_id==mode%id) then
    output = generate_mode_basis_1d( subspace,                 &
                                   & frequency,                &
                                   & mode,                     &
                                   & maximum_power,            &
                                   & potential_expansion_order )
  else
    output = generate_mode_basis_2d( subspace,                 &
                                   & frequency,                &
                                   & mode,                     &
                                   & maximum_power,            &
                                   & potential_expansion_order )
  endif
end function

function generate_mode_basis_1d(subspace,frequency,mode,maximum_power, &
   & potential_expansion_order) result(output)
  implicit none
  
  type(DegenerateSubspace), intent(in) :: subspace
  real(dp),                 intent(in) :: frequency
  type(ComplexMode),        intent(in) :: mode
  integer,                  intent(in) :: maximum_power
  integer,                  intent(in) :: potential_expansion_order
  type(WavevectorBasis)                :: output
  
  type(HarmonicState1D),      allocatable :: harmonic_states_1d(:)
  type(SubspaceStatePointer), allocatable :: harmonic_states(:)
  type(CoupledStates),        allocatable :: harmonic_couplings(:)
  
  integer  :: power
  integer  :: state
  integer  :: term
  real(dp) :: coefficient
  
  integer, allocatable :: coupling_ids(:)
  integer, allocatable :: separations(:)
  
  integer :: i,j,ialloc
  
  harmonic_states_1d = HarmonicState1D( id         = mode%id,           &
                                      & occupation = [( power,          &
                                      &                 power=0,        &
                                      &                 maximum_power)] )
  harmonic_states = [( SubspaceStatePointer(HarmonicStateReal(     &
                     &                 subspace%id,                &
                     &                 frequency,                  &
                     &                 [harmonic_states_1d(i)] )), &
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
   & potential_expansion_order) result(output)
  implicit none
  
  type(DegenerateSubspace), intent(in) :: subspace
  real(dp),                 intent(in) :: frequency
  type(ComplexMode),        intent(in) :: mode
  integer,                  intent(in) :: maximum_power
  integer,                  intent(in) :: potential_expansion_order
  type(WavevectorBasis)                :: output
  
  type(HarmonicState2D),      allocatable :: harmonic_states_2d(:)
  type(SubspaceStatePointer), allocatable :: harmonic_states(:)
  type(CoupledStates),        allocatable :: harmonic_couplings(:)
  
  integer :: total_power
  
  integer  :: no_states
  integer  :: power
  integer  :: state
  integer  :: term
  real(dp) :: coefficient
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
  harmonic_states = [( SubspaceStatePointer(HarmonicStateComplex(     &
                     &                    subspace%id,                &
                     &                    frequency,                  &
                     &                    [harmonic_states_2d(i)] )), &
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
      separation = abs( harmonic_states(i)%occupation() &
                    & - harmonic_states(j)%occupation() )
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
  
  integer  :: id
  real(dp) :: coefficient
  integer  :: separation
  
  ! Output variables.
  integer                                 :: maximum_power
  integer                                 :: expansion_order
  integer                                 :: subspace_id
  real(dp)                                :: frequency
  type(FractionVector)                    :: wavevector
  type(SubspaceStatePointer), allocatable :: harmonic_states(:)
  type(CoupledStates),        allocatable :: harmonic_couplings(:)
  
  integer :: i,j,k,l,ialloc
  
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
  if (any(that_occupations(2:)<that_occupations(:size(that%harmonic_states_)-1))) then
    call print_line(CODE_ERROR//': Powers are not in ascending order.')
    call err()
  endif
  
  ! Count the number of output states corresponding to each state in 'this'.
  no_states = [( count(that_occupations + this_occupations(i) <= maximum_power), &
               & i=1,                                                  &
               & size(this%harmonic_states_)                          )]
  this_to_output = [( sum(no_states(:i-1)),        &
                    & i=1,                         &
                    & size(this%harmonic_states_) )]
  
  ! Generate monomial and harmonic states.
  allocate(harmonic_states(sum(no_states)), stat=ialloc); call err(ialloc)
  do i=1,size(this%harmonic_states_)
    harmonic_states(this_to_output(i)+1:this_to_output(i)+no_states(i)) = &
       & prod(this%harmonic_states_(i),that%harmonic_states_(:no_states(i)))
  enddo
  
  ! Generate harmonic couplings.
  allocate( harmonic_couplings(size(harmonic_states)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(this%harmonic_states_)
    do j=1,no_states(i)
      harmonic_couplings(this_to_output(i)+j) = CoupledStates()
      do k=1,size(this%harmonic_couplings_(i))
        do l=1,size(that%harmonic_couplings_(j))
          if ( that%harmonic_couplings_(j)%id(l) >          &
             & no_states(this%harmonic_couplings_(i)%id(k)) ) then
            exit
          endif
          id = this_to_output(this%harmonic_couplings_(i)%id(k)) &
           & + that%harmonic_couplings_(j)%id(l)
          separation = this%harmonic_couplings_(i)%separation(k) &
                   & + that%harmonic_couplings_(j)%separation(l)
          if (separation<=expansion_order) then
            call harmonic_couplings(this_to_output(i)+j)%add_coupling( &
                                             & id         = id,        &
                                             & separation = separation )
          endif
        enddo
      enddo
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

impure elemental function prod(lhs,rhs) result(output)
  implicit none
  
  type(SubspaceStatePointer), intent(in) :: lhs
  type(SubspaceStatePointer), intent(in) :: rhs
  type(SubspaceStatePointer)             :: output
  
  class(SubspaceState), allocatable :: lhs2
  class(SubspaceState), allocatable :: rhs2
  
  lhs2 = lhs%state()
  rhs2 = rhs%state()
  
  select type(lhs2); type is(HarmonicStateReal)
    select type(rhs2); type is (HarmonicStateReal)
      output = SubspaceStatePointer(prod_real(lhs2,rhs2))
    end select
  type is(HarmonicStateComplex)
    select type(rhs2); type is(HarmonicStateComplex)
      output = SubspaceStatePointer(prod_complex(lhs2,rhs2))
    end select
  end select
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
  
  integer, allocatable :: qpoint_set(:)
  integer, allocatable :: wavevector_states(:)
  integer, allocatable :: new_ids(:)
  integer, allocatable :: ids(:)
  integer, allocatable :: separations(:)
  integer, allocatable :: allowed_couplings(:)
  
  ! Output variables.
  type(SubspaceStatePointer), allocatable :: harmonic_states(:)
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
    harmonic_states = input%harmonic_states_(wavevector_states)
    couplings       = input%harmonic_couplings_(wavevector_states)
    
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
  enddo
end function

! ----------------------------------------------------------------------
! Operations involving WavevectorStates.
! ----------------------------------------------------------------------
impure elemental function inner_product_WavevectorBasis(this,bra,ket, &
   & subspace,anharmonic_data) result(output)
  implicit none
  
  class(WavevectorBasis),   intent(in)           :: this
  class(BasisState),        intent(in)           :: bra
  class(BasisState),        intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  real(dp)                                       :: output
  
  type(WavevectorState) :: bra_
  type(WavevectorState) :: ket_
  
  bra_ = WavevectorState(bra)
  
  if (present(ket)) then
    ket_ = WavevectorState(ket)
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
  
  class(WavevectorBasis),   intent(in)           :: this
  class(BasisState),        intent(in)           :: bra
  type(SparseMonomial),     intent(in)           :: monomial
  class(BasisState),        intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  complex(dp)                                    :: output
  
  type(WavevectorState) :: bra_
  type(WavevectorState) :: ket_
  
  complex(dp) :: term
  
  integer :: i,j,k
  
  bra_ = WavevectorState(bra)
  if (present(ket)) then
    ket_ = WavevectorState(ket)
  else
    ket_ = bra_
  endif
  
  output = 0.0_dp
  do i=1,size(this%harmonic_states_)
    do j=1,size(this%harmonic_couplings_(i))
      k = this%harmonic_couplings_(i)%id(j)
      
      ! Ignore terms with k<i. These are added when i and k are reversed.
      if (k<i) then
        cycle
      endif
      
      term = this%harmonic_states_(i)%integrate( monomial,                 &
                                               & this%harmonic_states_(k), &
                                               & anharmonic_data           )
      
      term = bra_%coefficients(i) * term * ket_%coefficients(k)
      
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
function integrate_SparseMonomial_WavevectorBases(bases,states, &
   & thermal_energy,monomial,anharmonic_data) result(output)
  implicit none
  
  class(WavevectorBasis), intent(in) :: bases(:)
  class(BasisStates),     intent(in) :: states
  real(dp),               intent(in) :: thermal_energy
  type(SparseMonomial),   intent(in) :: monomial
  type(AnharmonicData),   intent(in) :: anharmonic_data
  complex(dp)                        :: output
  
  type(WavevectorStates) :: states_
  
  integer, allocatable :: at_wavevector(:)
  
  complex(dp) :: term
  
  integer :: i,j,k,l,m
  
  states_ = WavevectorStates(states)
  
  output = 0.0_dp
  ! Loop over wavevector bases, summing the contribution from each.
  do i=1,size(bases)
    ! Filter the states to get only the states at this wavevector.
    at_wavevector = filter(states_%states%wavevector==bases(i)%wavevector)
    
    ! Perform a double loop across harmonic states, including only those
    !   for which <j|X|l> can be non-zero.
    do j=1,size(bases(i)%harmonic_states_)
      do k=1,size(bases(i)%harmonic_couplings_(j))
        l = bases(i)%harmonic_couplings_(j)%id(k)
        
        ! Ignore terms with l<j. These are added when j and l are reversed.
        if (l<j) then
          cycle
        endif
        
        ! Calculate <j|X|l> in the harmonic basis.
        term = bases(i)%harmonic_states_(j)%integrate( &
                       & monomial,                     &
                       & bases(i)%harmonic_states_(l), &
                       & anharmonic_data               )
        
        ! Calculate the contribution to the thermal average of <X>
        !    from <j|X|l>.
        ! <X> = sum_m <m|X|m> * w(m), where {|m>} are the eigenstates,
        !    and w(m) is the thermal weight of eigenstate m.
        ! <X> = sum_m (sum_j U(m,j)<j|) X (sum_l U(m,l)|l>) w(m)
        !     = sum_j sum_l ( <j|X|l> * sum_m(U(m,j) U(m,l) w(m)) )
        ! U(m,j) is the j'th component of the m'th eigenvector.
        term = term                                                      &
           & * sum([(   states_%states(at_wavevector(m))%coefficients(j) &
           &          * states_%states(at_wavevector(m))%coefficients(l) &
           &          * states_%weights(at_wavevector(m)),               &
           &          m=1,                                               &
           &          size(at_wavevector)                                )])
        
        ! If j==l, only include <j|X|j>, otherwise include <l|X|j> and <j|X|l>.
        if (l==j) then
          output = output + term
        else
          output = output + 2*term
        endif
      enddo
    enddo
  enddo
end function

impure elemental function kinetic_energy_WavevectorBasis(this,bra,ket, &
   & subspace,anharmonic_data) result(output)
  implicit none
  
  class(WavevectorBasis),   intent(in)           :: this
  class(BasisState),        intent(in)           :: bra
  class(BasisState),        intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  real(dp)                                       :: output
  
  type(WavevectorState) :: bra_
  type(WavevectorState) :: ket_
  
  real(dp) :: term
  
  integer :: i,j,k
  
  bra_ = WavevectorState(bra)
  if (present(ket)) then
    ket_ = WavevectorState(ket)
  else
    ket_ = bra_
  endif
  
  output = 0
  do i=1,size(this%harmonic_states_)
    do j=1,size(this%harmonic_couplings_(i))
      k = this%harmonic_couplings_(i)%id(j)
      term = this%harmonic_states_(i)%kinetic_energy( &
                          & this%harmonic_states_(k), &
                          & anharmonic_data            )
      output = output               &
           & + bra_%coefficients(i) &
           & * term                 &
           & * ket_%coefficients(k)
    enddo
  enddo
end function

impure elemental function harmonic_potential_energy_WavevectorBasis(this, &
   & bra,ket,subspace,anharmonic_data) result(output)
  implicit none
  
  class(WavevectorBasis),   intent(in)           :: this
  class(BasisState),        intent(in)           :: bra
  class(BasisState),        intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  real(dp)                                       :: output
  
  type(WavevectorState) :: bra_
  type(WavevectorState) :: ket_
  
  real(dp) :: term
  
  integer :: i,j,k
  
  bra_ = WavevectorState(bra)
  if (present(ket)) then
    ket_ = WavevectorState(ket)
  else
    ket_ = bra_
  endif
  
  output = 0
  do i=1,size(this%harmonic_states_)
    do j=1,size(this%harmonic_couplings_(i))
      k = this%harmonic_couplings_(i)%id(j)
      term = this%harmonic_states_(i)%harmonic_potential_energy( &
                                     & this%harmonic_states_(k), &
                                     & anharmonic_data           )
      if (present(ket)) then
        output = output               &
             & + bra_%coefficients(i) &
             & * term                 &
             & * ket_%coefficients(k)
      else
        output = output               &
             & + bra_%coefficients(i) &
             & * term                 &
             & * bra_%coefficients(k)
      endif
    enddo
  enddo
end function

impure elemental function kinetic_stress_WavevectorBasis(this,bra,ket, &
   & subspace,stress_prefactors,anharmonic_data) result(output)
  implicit none
  
  class(WavevectorBasis),   intent(in)           :: this
  class(BasisState),        intent(in)           :: bra
  class(BasisState),        intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  type(StressPrefactors),   intent(in)           :: stress_prefactors
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(RealMatrix)                               :: output
  
  type(WavevectorState) :: bra_
  type(WavevectorState) :: ket_
  
  type(RealMatrix) :: term
  
  integer :: i,j,k
  
  bra_ = WavevectorState(bra)
  if (present(ket)) then
    ket_ = WavevectorState(ket)
  else
    ket_ = bra_
  endif
  
  output = dblemat(zeroes(3,3))
  do i=1,size(this%harmonic_states_)
    do j=1,size(this%harmonic_couplings_(i))
      k = this%harmonic_couplings_(i)%id(j)
      term = this%harmonic_states_(i)%kinetic_stress( &
                          & this%harmonic_states_(k), &
                          & stress_prefactors,        &
                          & anharmonic_data           )
      output = output               &
           & + bra_%coefficients(i) &
           & * term                 &
           & * ket_%coefficients(k)
    enddo
  enddo
end function

! ----------------------------------------------------------------------
! Operations which generate WavevectorStates.
! ----------------------------------------------------------------------
impure elemental function initial_ground_state_WavevectorBasis(this) &
   & result(output)
  implicit none
  
  class(WavevectorBasis), intent(in) :: this
  type(WavevectorState)              :: output
  
  real(dp), allocatable :: coefficients(:)
  
  integer :: i
  
  coefficients = [( 0.0_dp, i=1, size(this%harmonic_states_) )]
  coefficients(first(this%harmonic_states_%occupation()==0)) = 1
  
  output = WavevectorState(this%subspace_id,this%wavevector,coefficients)
end function

impure elemental function initial_states_WavevectorBasis(this,subspace, &
   & thermal_energy,anharmonic_data) result(output)
  implicit none
  
  class(WavevectorBasis),   intent(in) :: this
  type(DegenerateSubspace), intent(in) :: subspace
  real(dp),                 intent(in) :: thermal_energy
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(BasisStatesPointer)             :: output
  
  output = BasisStatesPointer(WavevectorStates(     &
     & subspace_id = this%subspace_id,              &
     & states      = [this%initial_ground_state()], &
     & energies    = [0.0_dp],                      &
     & weights     = [1.0_dp]                       ))
end function

impure elemental function calculate_states_WavevectorBasis(this,subspace, &
   & subspace_potential,thermal_energy,energy_convergence,                &
   & no_converged_calculations,max_pulay_iterations,pre_pulay_iterations, &
   & pre_pulay_damping,anharmonic_data) result(output)
  implicit none
  
  class(WavevectorBasis),   intent(in) :: this
  type(DegenerateSubspace), intent(in) :: subspace
  class(PotentialData),     intent(in) :: subspace_potential
  real(dp),                 intent(in) :: thermal_energy
  real(dp),                 intent(in) :: energy_convergence
  integer,                  intent(in) :: no_converged_calculations
  integer,                  intent(in) :: max_pulay_iterations
  integer,                  intent(in) :: pre_pulay_iterations
  real(dp),                 intent(in) :: pre_pulay_damping
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(BasisStatesPointer)             :: output
  
  type(SubspaceStatePointer)             :: bra
  type(SubspaceStatePointer)             :: ket
  real(dp),                  allocatable :: hamiltonian(:,:)
  type(SymmetricEigenstuff), allocatable :: estuff(:)
  
  type(WavevectorState), allocatable :: wavevector_states(:)
  
  integer :: i,j,k,ialloc
  
  allocate( hamiltonian( size(this%harmonic_states_),    &
          &              size(this%harmonic_states_)  ), &
          & stat=ialloc); call err(ialloc)
  hamiltonian = 0
  do i=1,size(this%harmonic_states_)
    bra = this%harmonic_states_(i)
    do j=1,size(this%harmonic_couplings_(i))
      k = this%harmonic_couplings_(i)%id(j)
      ket = this%harmonic_states_(k)
      
      hamiltonian(i,k) = bra%kinetic_energy( ket,              &
                     &                       anharmonic_data ) &
                     & + potential_energy( bra,                &
                     &                     subspace_potential, &
                     &                     ket,                &
                     &                     anharmonic_data )
    enddo
  enddo
  
  estuff = diagonalise_symmetric(hamiltonian)
  
  wavevector_states = [( WavevectorState( this%subspace_id,    &
                      &                   this%wavevector,     &
                      &                   estuff(i)%evec    ), &
                      & i=1,                                   &
                      & size(estuff)                           )]
  
  ! N.B. there is not enough information at this stage to calculate weights.
  ! These must be calculated once the states from each wavevector are collated.
  output = BasisStatesPointer(WavevectorStates(  &
     & subspace_id = this%subspace_id,           &
     & states      = wavevector_states,          &
     & energies    = estuff%eval,                &
     & weights     = [(0.0_dp,i=1,size(estuff))] ))
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
  
  type(WavevectorBasis),  intent(in)           :: bases(:)
  real(dp),               intent(in)           :: thermal_energy
  class(StressData),      intent(in), optional :: stress
  type(StressPrefactors), intent(in), optional :: stress_prefactors
  type(AnharmonicData),   intent(in)           :: anharmonic_data
  type(ThermodynamicData)                      :: output
  
  real(dp), allocatable :: kinetic_energies(:)
  
  type(RealMatrix), allocatable :: stresses(:)
  real(dp),         allocatable :: volume
  
  integer :: i,j,k,ialloc
  
  ! Calculate <|T|>. N.B. this is also <|V|> for the harmonic potential.
  kinetic_energies = [( bases(i)%harmonic_states_%kinetic_energy(    &
                      &           anharmonic_data=anharmonic_data ), &
                      & i=1,                                         &
                      & size(bases)                                  )]
  
  if (present(stress) .neqv. present(stress_prefactors)) then
    call print_line(CODE_ERROR//': Either both or neither of stress and &
       &stress_prefactors must be present.')
    call err()
  endif
  if (present(stress)) then
    allocate(stresses(size(kinetic_energies)), stat=ialloc); call err(ialloc)
    k = 0
    do i=1,size(bases)
      do j=1,size(bases(i)%harmonic_states_)
        k = k+1
        stresses(k) = bases(i)%harmonic_states_(j)%kinetic_stress(       &
                  &         stress_prefactors = stress_prefactors,       &
                  &         anharmonic_data   = anharmonic_data    )     &
                  & + potential_stress(                                  &
                  &      state           = bases(i)%harmonic_states_(j), &
                  &      stress          = stress,                       &
                  &      anharmonic_data = anharmonic_data               )
      enddo
    enddo
    
    volume = anharmonic_data%structure%volume
  endif
  
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
  
  type(WavevectorBasis),  intent(in)           :: bases(:)
  real(dp),               intent(in)           :: thermal_energy
  class(PotentialData),   intent(in)           :: potential
  class(StressData),      intent(in), optional :: stress
  type(StressPrefactors), intent(in), optional :: stress_prefactors
  type(AnharmonicData),   intent(in)           :: anharmonic_data
  type(ThermodynamicData)                      :: output
  
  real(dp), allocatable :: harmonic_kinetic_energies(:)
  real(dp), allocatable :: harmonic_energy_differences(:)
  real(dp), allocatable :: anharmonic_potential_energies(:)
  real(dp), allocatable :: boltzmann_factors(:)
  
  integer :: ground_state
  
  real(dp) :: anharmonic_minus_harmonic
  
  integer :: i,j
  
  ! Calculate properties with a harmonic potential.
  output = core_harmonic_observables( bases,             &
                                    & thermal_energy,    &
                                    & stress,            &
                                    & stress_prefactors, &
                                    & anharmonic_data    )
  
  ! Calculate <|T|> for the harmonic potential.
  harmonic_kinetic_energies = [(                    &
     & bases(i)%harmonic_states_%kinetic_energy(    &
     &           anharmonic_data=anharmonic_data ), &
     & i=1,                                         &
     & size(bases)                                  )]
  
  ! Identify the ground state.
  ground_state = minloc(harmonic_kinetic_energies, 1)
  
  ! Calculate <|T+V|> - min(<|T+V|>) for the harmonic potential.
  ! N.B. <|V|> = <|T|> for the harmonic potential.
  harmonic_energy_differences = 2*( harmonic_kinetic_energies               &
                                & - harmonic_kinetic_energies(ground_state) )
  
  ! Calculate <|V|> for the input potential.
  anharmonic_potential_energies = [(                          &
     & [( potential_energy( bases(i)%harmonic_states_(j),     &
     &                      potential,                        &
     &                      anharmonic_data           ),      &
     &    j=1,                                                &
     &    size(bases(i)%harmonic_states_)                 )], &
     & i=1,                                                   &
     & size(bases)                                            )]
  
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
  
  type(WavevectorBasis),    intent(in)           :: bases(:)
  real(dp),                 intent(in)           :: thermal_energy
  class(BasisStates),       intent(in)           :: states
  type(DegenerateSubspace), intent(in)           :: subspace
  class(PotentialData),     intent(in)           :: subspace_potential
  class(StressData),        intent(in), optional :: subspace_stress
  type(StressPrefactors),   intent(in), optional :: stress_prefactors
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(ThermodynamicData)                        :: output
  
  type(WavevectorStates) :: wavevector_states
  
  type(RealMatrix), allocatable :: stress(:)
  real(dp),         allocatable :: volume
  
  integer :: i,j,ialloc
  
  if (present(subspace_stress) .neqv. present(stress_prefactors)) then
    call print_line(CODE_ERROR//': Only one of subspace_stress and &
       &stress_prefactors passed.')
    call err()
  endif
  
  wavevector_states = WavevectorStates(states)
  
  ! Calculate stress.
  if (present(subspace_stress)) then
    allocate( stress(size(wavevector_states%states)), &
            & stat=ialloc); call err(ialloc)
    do i=1,size(wavevector_states%states)
      associate(state => wavevector_states%states(i))
        j = first(bases%wavevector==state%wavevector)
        stress(i) = potential_stress( state,                  &
                &                     subspace_stress,        &
                &                     subspace,               &
                &                     bases(j),               &
                &                     anharmonic_data  )      &
                & + bases(j)%kinetic_stress(                  &
                &      bra               = state,             &
                &      subspace          = subspace,          &
                &      stress_prefactors = stress_prefactors, &
                &      anharmonic_data   = anharmonic_data    )
      end associate
    enddo
    
    volume = anharmonic_data%structure%volume
  endif
  
  output = ThermodynamicData( thermal_energy,             &
                            & wavevector_states%energies, &
                            & stress,                     &
                            & volume                      )
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
    
    this = WavevectorBasis( maximum_power,   &
                          & expansion_order, &
                          & subspace_id,     &
                          & frequency,       &
                          & wavevector,      &
                          & harmonic_states, &
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
      state_strings = str(this%harmonic_states_(i))
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
impure elemental function integrate_BasisStates_WavevectorBasis(this,states, &
   & thermal_energy,monomial,subspace,anharmonic_data) result(output)
  implicit none
  
  class(WavevectorBasis),   intent(in) :: this
  class(BasisStates),       intent(in) :: states
  real(dp),                 intent(in) :: thermal_energy
  type(SparseMonomial),     intent(in) :: monomial
  type(DegenerateSubspace), intent(in) :: subspace
  type(AnharmonicData),     intent(in) :: anharmonic_data
  complex(dp)                          :: output
  
  call err()
end function

impure elemental function thermodynamic_data_WavevectorBasis(this,      &
   & thermal_energy,states,subspace,subspace_potential,subspace_stress, &
   & stress_prefactors,anharmonic_data) result(output)
  implicit none
  
  class(WavevectorBasis),   intent(in) :: this
  real(dp),                 intent(in)           :: thermal_energy
  class(BasisStates),       intent(in)           :: states
  type(DegenerateSubspace), intent(in)           :: subspace
  class(PotentialData),     intent(in)           :: subspace_potential
  class(StressData),        intent(in), optional :: subspace_stress
  type(StressPrefactors),   intent(in), optional :: stress_prefactors
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(ThermodynamicData)                        :: output
  
  call err()
end function

impure elemental function wavefunctions_WavevectorBasis(this,states,subspace, &
   & anharmonic_data) result(output)
  implicit none
  
  class(WavevectorBasis),   intent(in) :: this
  class(BasisStates),       intent(in) :: states
  type(DegenerateSubspace), intent(in) :: subspace
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(SubspaceWavefunctionsPointer)   :: output
  
  call err()
end function
end module
