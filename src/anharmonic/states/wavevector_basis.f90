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
  
  public :: WavevectorBasis
  
  public :: harmonic_thermodynamics
  public :: harmonic_expectation
  
  type, extends(Stringsable) :: WavevectorBasis
    integer                                          :: maximum_power
    integer                                          :: expansion_order
    integer                                          :: subspace_id
    real(dp)                                         :: frequency
    type(FractionVector)                             :: wavevector
    type(SubspaceStatePointer), allocatable, private :: harmonic_states_(:)
    type(CoupledStates),        allocatable, private :: harmonic_couplings_(:)
  contains
    ! Set the frequency of the basis.
    procedure, public :: set_frequency => set_frequency_WavevectorBasis
    
    ! State procedures.
    procedure, public :: inner_product => &
                       & inner_product_WavevectorBasis
    procedure, public :: braket => &
                       & braket_WavevectorBasis
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
    
    ! I/O.
    procedure, public :: read  => read_WavevectorBasis
    procedure, public :: write => write_WavevectorBasis
  end type
  
  interface WavevectorBasis
    module procedure new_WavevectorBasis
    module procedure new_WavevectorBasis_subspace
    module procedure new_WavevectorBasis_Strings
    module procedure new_WavevectorBasis_StringArray
  end interface
  
  interface harmonic_thermodynamics
    module procedure harmonic_thermodynamics_WavevectorBasis
  end interface
  
  interface harmonic_expectation
    module procedure harmonic_expectation_WavevectorBasis
  end interface
contains

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

! Set the frequency of the basis.
impure elemental subroutine set_frequency_WavevectorBasis(this,frequency)
  implicit none
  
  class(WavevectorBasis), intent(inout) :: this
  real(dp),               intent(in)    :: frequency
  
  integer :: i
  
  this%frequency = frequency
  
  do i=1,size(this%harmonic_states_)
    call this%harmonic_states_(i)%set_frequency(frequency)
  enddo
end subroutine

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
   & anharmonic_data) result(output)
  implicit none
  
  class(WavevectorBasis), intent(in)           :: this
  type(WavevectorState),  intent(in)           :: bra
  type(WavevectorState),  intent(in), optional :: ket
  type(AnharmonicData),   intent(in)           :: anharmonic_data
  real(dp)                                     :: output
  
  if (present(ket)) then
    ! The basis is orthonormal, so the dot product of the states is
    !    simply the dot product of the coefficients.
    output = dot_product(bra%coefficients, ket%coefficients)
  else
    ! States are orthonormal, so <p|p>=1.
    output = 1.0_dp
  endif
end function

impure elemental function braket_WavevectorBasis(this,bra,monomial,ket, &
   & anharmonic_data) result(output)
  implicit none
  
  class(WavevectorBasis), intent(in)           :: this
  type(WavevectorState),  intent(in)           :: bra
  type(ComplexMonomial),  intent(in)           :: monomial
  type(WavevectorState),  intent(in), optional :: ket
  type(AnharmonicData),   intent(in)           :: anharmonic_data
  type(ComplexMonomial)                        :: output
  
  type(ComplexMonomial) :: term
  
  integer :: i,j,k
  
  do i=1,size(this%harmonic_states_)
    do j=1,size(this%harmonic_couplings_(i))
      k = this%harmonic_couplings_(i)%id(j)
      term = this%harmonic_states_(i)%braket( monomial,                 &
                                            & this%harmonic_states_(k), &
                                            & anharmonic_data           )
      
      if (i==1 .and. j==1) then
        output = term
        output%coefficient = 0
      endif
      
      if (present(ket)) then
        output%coefficient = output%coefficient  &
                         & + bra%coefficients(i) &
                         & * term%coefficient    &
                         & * ket%coefficients(k)
      else
        output%coefficient = output%coefficient  &
                         & + bra%coefficients(i) &
                         & * term%coefficient    &
                         & * bra%coefficients(k)
      endif
    enddo
  enddo
end function

impure elemental function kinetic_energy_WavevectorBasis(this,bra,ket, &
   & anharmonic_data) result(output)
  implicit none
  
  class(WavevectorBasis), intent(in)           :: this
  type(WavevectorState),  intent(in)           :: bra
  type(WavevectorState),  intent(in), optional :: ket
  type(AnharmonicData),   intent(in)           :: anharmonic_data
  real(dp)                                     :: output
  
  real(dp) :: term
  
  integer :: i,j,k
  
  output = 0
  do i=1,size(this%harmonic_states_)
    do j=1,size(this%harmonic_couplings_(i))
      k = this%harmonic_couplings_(i)%id(j)
      term = this%harmonic_states_(i)%kinetic_energy( &
                          & this%harmonic_states_(k), &
                          & anharmonic_data            )
      if (present(ket)) then
        output = output              &
             & + bra%coefficients(i) &
             & * term                &
             & * ket%coefficients(k)
      else
        output = output              &
             & + bra%coefficients(i) &
             & * term                &
             & * bra%coefficients(k)
      endif
    enddo
  enddo
end function

impure elemental function harmonic_potential_energy_WavevectorBasis(this, &
   & bra,ket,anharmonic_data) result(output)
  implicit none
  
  class(WavevectorBasis), intent(in)           :: this
  type(WavevectorState),  intent(in)           :: bra
  type(WavevectorState),  intent(in), optional :: ket
  type(AnharmonicData),   intent(in)           :: anharmonic_data
  real(dp)                                     :: output
  
  real(dp) :: term
  
  integer :: i,j,k
  
  output = 0
  do i=1,size(this%harmonic_states_)
    do j=1,size(this%harmonic_couplings_(i))
      k = this%harmonic_couplings_(i)%id(j)
      term = this%harmonic_states_(i)%harmonic_potential_energy( &
                                     & this%harmonic_states_(k), &
                                     & anharmonic_data           )
      if (present(ket)) then
        output = output              &
             & + bra%coefficients(i) &
             & * term                &
             & * ket%coefficients(k)
      else
        output = output              &
             & + bra%coefficients(i) &
             & * term                &
             & * bra%coefficients(k)
      endif
    enddo
  enddo
end function

impure elemental function kinetic_stress_WavevectorBasis(this,bra,ket, &
   & stress_prefactors,anharmonic_data) result(output)
  implicit none
  
  class(WavevectorBasis), intent(in)           :: this
  type(WavevectorState),  intent(in)           :: bra
  type(WavevectorState),  intent(in), optional :: ket
  type(StressPrefactors), intent(in)           :: stress_prefactors
  type(AnharmonicData),   intent(in)           :: anharmonic_data
  type(RealMatrix)                             :: output
  
  type(RealMatrix) :: term
  
  integer :: i,j,k
  
  output = dblemat(zeroes(3,3))
  do i=1,size(this%harmonic_states_)
    do j=1,size(this%harmonic_couplings_(i))
      k = this%harmonic_couplings_(i)%id(j)
      term = this%harmonic_states_(i)%kinetic_stress( &
                          & this%harmonic_states_(k), &
                          & stress_prefactors,        &
                          & anharmonic_data           )
      if (present(ket)) then
        output = output              &
             & + bra%coefficients(i) &
             & * term                &
             & * ket%coefficients(k)
      else
        output = output              &
             & + bra%coefficients(i) &
             & * term                &
             & * bra%coefficients(k)
      endif
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

function initial_states_WavevectorBasis(this,anharmonic_data) &
   & result(output)
  implicit none
  
  class(WavevectorBasis), intent(in) :: this
  type(AnharmonicData),   intent(in) :: anharmonic_data
  type(WavevectorStates)             :: output
  
  output = WavevectorStates( this%subspace_id,              &
                           & [this%initial_ground_state()], &
                           & [0.0_dp]                       )
end function

function calculate_states_WavevectorBasis(this,potential,anharmonic_data) &
   & result(output)
  implicit none
  
  class(WavevectorBasis), intent(in) :: this
  class(PotentialData),   intent(in) :: potential
  type(AnharmonicData),   intent(in) :: anharmonic_data
  type(WavevectorStates)             :: output
  
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
                     &                     potential,          &
                     &                     ket,                &
                     &                     anharmonic_data )   &
                     & * anharmonic_data%anharmonic_supercell%sc_size
    enddo
  enddo
  
  estuff = diagonalise_symmetric(hamiltonian)
  
  wavevector_states = [( WavevectorState( this%subspace_id,    &
                      &                   this%wavevector,     &
                      &                   estuff(i)%evec    ), &
                      & i=1,                                   &
                      & size(estuff)                           )]
  
  output = WavevectorStates(this%subspace_id, wavevector_states, estuff%eval)
end function

! Calculate the harmonic expectation of the harmonic potential,
!    but only using the states in the basis.
! N.B. this normalises the sum of thermal weights of the included states
!    to be one,
!    rather than normalising to the sum across all states to be one.
function harmonic_thermodynamics_WavevectorBasis(bases,thermal_energy, &
   & anharmonic_data) result(output)
  implicit none
  
  type(WavevectorBasis), intent(in) :: bases(:)
  real(dp),              intent(in) :: thermal_energy
  type(AnharmonicData),  intent(in) :: anharmonic_data
  type(ThermodynamicData)           :: output
  
  real(dp), allocatable :: kinetic_energies(:)
  
  integer :: i
  
  ! Calculate <|T|>. N.B. this is also <|V|> for the harmonic potential.
  kinetic_energies = [( bases(i)%harmonic_states_%kinetic_energy(    &
                      &           anharmonic_data=anharmonic_data ), &
                      & i=1,                                         &
                      & size(bases)                                  )]
  
  output = ThermodynamicData(thermal_energy, 2*kinetic_energies)
end function

! Calculate the harmonic expectation of a potential,
!    but only using the states in the basis.
! N.B. this normalises the sum of thermal weights of the included states
!    to be one,
!    rather than normalising to the sum across all states to be one.
function harmonic_expectation_WavevectorBasis(bases,potential, &
   & thermal_energy,anharmonic_data) result(output)
  implicit none
  
  type(WavevectorBasis), intent(in) :: bases(:)
  class(PotentialData),  intent(in) :: potential
  real(dp),              intent(in) :: thermal_energy
  type(AnharmonicData),  intent(in) :: anharmonic_data
  type(ThermodynamicData)           :: output
  
  real(dp), allocatable :: kinetic_energies(:)
  real(dp), allocatable :: energies(:)
  real(dp), allocatable :: boltzmann_factors(:)
  real(dp)              :: partition_function
  
  integer               :: ground_state
  real(dp), allocatable :: energy_difference(:)
  
  real(dp) :: anharmonic_minus_harmonic
  
  real(dp) :: energy
  real(dp) :: free_energy
  real(dp) :: entropy
  
  integer :: i,j,ialloc
  
  ! Calculate <|T|>. N.B. this is also <|V|> for the harmonic potential.
  kinetic_energies = [( bases(i)%harmonic_states_%kinetic_energy(    &
                      &           anharmonic_data=anharmonic_data ), &
                      & i=1,                                         &
                      & size(bases)                                  )]
  
  ! Calculate <|V|> for the input potential.
  energies = [(                                               &
     & [( potential_energy( bases(i)%harmonic_states_(j),     &
     &                      potential,                        &
     &                      anharmonic_data           ),      &
     &    j=1,                                                &
     &    size(bases(i)%harmonic_states_)                 )], &
     & i=1,                                                   &
     & size(bases)                                            )]
  energies = energies * anharmonic_data%anharmonic_supercell%sc_size
  
  ! Identify the ground state.
  ground_state = minloc(kinetic_energies, 1)
  
  ! Calculate E-E' = <|T+V|>-<|T+V|>' = 2<T>-2<T>'.
  energy_difference = [( 2.0_dp*( kinetic_energies(i)               &
                       &        - kinetic_energies(ground_state) ), &
                       & i=1,                                       &
                       & size(kinetic_energies)                     )]
  
  if (maxval(energy_difference) < 1e20_dp*thermal_energy) then
    ! Normal temperature range.
    boltzmann_factors = [exp(-energy_difference/thermal_energy)]
    partition_function = sum(boltzmann_factors)
    
    ! Calculate observables w/r/t the harmonic potential.
    energy = sum(boltzmann_factors*kinetic_energies) * 2 / partition_function
    free_energy = 2*kinetic_energies(ground_state) &
              & - thermal_energy*log(partition_function)
    entropy = (energy-free_energy)/thermal_energy
    
    ! Update the energy and free energy to account for the input potential.
    ! N.B. again, the harmonic kinetic and potential energies are the same.
    anharmonic_minus_harmonic = sum( boltzmann_factors             &
                            &      * (energies-kinetic_energies) ) &
                            & / partition_function
    energy = energy + anharmonic_minus_harmonic
    free_energy = free_energy + anharmonic_minus_harmonic
  else
    ! Very low temperature limit.
    energy = energies(ground_state) + kinetic_energies(ground_state)
    free_energy = energy
    entropy = 0
  endif
  
  output = ThermodynamicData(thermal_energy, energy, free_energy, entropy)
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
end module
