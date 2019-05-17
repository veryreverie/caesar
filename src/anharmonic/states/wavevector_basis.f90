! ======================================================================
! A basis of monomial and harmonic states for a given subspace
!    at a given wavevector.
! N.B. the wavevector in this context is the wavevector of the state,
!    not the wavevector of the mode.
! e.g. the state |p> along a mode at q-point q has a wavevector of p*q.
! ======================================================================
module wavevector_basis_module
  use common_module
  
  use monomial_state_module
  use harmonic_state_module
  use state_conversion_module
  use coupled_states_module
  implicit none
  
  private
  
  public :: WavevectorBasis
  public :: size
  
  type, extends(Stringsable) :: WavevectorBasis
    integer                                     :: maximum_power
    integer                                     :: expansion_order
    integer                                     :: subspace_id
    real(dp)                                    :: frequency
    type(FractionVector)                        :: wavevector
    integer                                     :: degeneracy
    type(MonomialState),   allocatable          :: monomial_states(:)
    type(HarmonicState),   allocatable          :: harmonic_states(:)
    type(StateConversion), allocatable, private :: states_to_basis_(:)
    type(StateConversion), allocatable, private :: basis_to_states_(:)
    type(CoupledStates),   allocatable          :: harmonic_couplings(:)
  contains
    ! Set the frequency of the basis.
    procedure, public :: set_frequency => set_frequency_WavevectorBasis
    
    ! Transform vectors of coefficients.
    procedure, public :: coefficients_states_to_basis
    procedure, public :: coefficients_basis_to_states
    
    ! Transform operator matrices.
    procedure, public :: operator_states_to_basis
    procedure, public :: operator_basis_to_states
    
    ! Transform Hamiltonian matrices.
    procedure, public :: hamiltonian_states_to_basis
    procedure, public :: hamiltonian_basis_to_states
    
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
  
  interface size
    module procedure size_WavevectorBasis
  end interface
contains

! Constructor and size function.
function new_WavevectorBasis(maximum_power,expansion_order,subspace_id,    &
   & frequency,wavevector,degeneracy,monomial_states,harmonic_states, &
   & states_to_basis,basis_to_states,harmonic_couplings) result(this)
  implicit none
  
  integer,               intent(in) :: maximum_power
  integer,               intent(in) :: expansion_order
  integer,               intent(in) :: subspace_id
  real(dp),              intent(in) :: frequency
  type(FractionVector),  intent(in) :: wavevector
  integer,               intent(in) :: degeneracy
  type(MonomialState),   intent(in) :: monomial_states(:)
  type(HarmonicState),   intent(in) :: harmonic_states(:)
  type(StateConversion), intent(in) :: states_to_basis(:)
  type(StateConversion), intent(in) :: basis_to_states(:)
  type(CoupledStates),   intent(in) :: harmonic_couplings(:)
  type(WavevectorBasis)             :: this
  
  if (size(harmonic_states)/=size(monomial_states)) then
    call print_line(CODE_ERROR//': monomial states and harmonic states do not &
       &match.')
    call err()
  elseif (size(states_to_basis)/=size(monomial_states)) then
    call print_line(CODE_ERROR//': monomial states and states to basis do not &
       &match.')
    call err()
  elseif (size(basis_to_states)/=size(monomial_states)) then
    call print_line(CODE_ERROR//': monomial states and basis to states do not &
       &match.')
    call err()
  elseif (size(harmonic_couplings)/=size(monomial_states)) then
    call print_line(CODE_ERROR//': monomial states and harmonic couplings do &
       &not match.')
  endif
  
  this%maximum_power      = maximum_power
  this%expansion_order    = expansion_order
  this%subspace_id        = subspace_id
  this%frequency          = frequency
  this%wavevector         = wavevector
  this%degeneracy         = degeneracy
  this%monomial_states    = monomial_states
  this%harmonic_states    = harmonic_states
  this%states_to_basis_   = states_to_basis
  this%basis_to_states_   = basis_to_states
  this%harmonic_couplings = harmonic_couplings
end function

function size_WavevectorBasis(this) result(output)
  implicit none
  
  type(WavevectorBasis), intent(in) :: this
  integer                           :: output
  
  output = size(this%monomial_states)
end function

! Set the frequency of the basis.
impure elemental subroutine set_frequency_WavevectorBasis(this,frequency)
  implicit none
  
  class(WavevectorBasis), intent(inout) :: this
  real(dp),               intent(in)    :: frequency
  
  this%frequency = frequency
  this%monomial_states%frequency = frequency
  this%harmonic_states%frequency = frequency
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
  
  integer, allocatable :: ids(:)
  
  type(ComplexUnivariate), allocatable :: univariates(:)
  type(ComplexMonomial),   allocatable :: monomials(:)
  type(MonomialState),     allocatable :: monomial_states(:)
  type(HarmonicState),     allocatable :: harmonic_states(:)
  type(StateConversion),   allocatable :: states_to_basis(:)
  type(StateConversion),   allocatable :: basis_to_states(:)
  type(CoupledStates),     allocatable :: harmonic_couplings(:)
  
  integer               :: no_states
  integer               :: power
  integer               :: state
  integer               :: term
  real(dp)              :: coefficient
  type(StateConversion) :: old
  
  integer, allocatable :: coupling_ids(:)
  integer, allocatable :: separations(:)
  
  integer :: i,j,ialloc
  
  no_states = 1+maximum_power
  
  ! Generate the monomial and harmonic states.
  ! The state monomial_states(id) is |u^{id-1}>, so ID=j+1.
  ! Ditto for harmonic states.
  allocate(ids(0:maximum_power), stat=ialloc); call err(ialloc)
  do power=0,maximum_power
    ids(power) = power+1
  enddo
  univariates = [(ComplexUnivariate(mode,power), power=0, maximum_power)]
  monomials = [( ComplexMonomial( coefficient = cmplx(1.0_dp,0.0_dp,dp),    &
               &                  modes       = [univariates(i)]         ), &
               & i=1,                                                       &
               & no_states                                                  )]
  monomials(1) = ComplexMonomial( coefficient = cmplx(1.0_dp,0.0_dp,dp), &
                                & modes       = [ComplexUnivariate::]    )
  monomial_states = [( MonomialState(subspace%id, frequency, monomials(i)), &
                     & i=1,                                                 &
                     & no_states                                            )]
  harmonic_states = [( HarmonicState(subspace%id, frequency, monomials(i)), &
                     & i=1,                                                 &
                     & no_states                                            )]
  
  ! A harmonic basis state |i> is defined in terms of monomial states |j> as
  !    |i> = sum h[i,j]|j>.
  ! h[0,0] = 1.
  ! h[i,j] = sqrt((2j-1)/i)       h[i-1,j-1]
  !        - (j+1)/sqrt(i(2j+1))) h[i-1,j+1]
  allocate(states_to_basis(no_states), stat=ialloc); call err(ialloc)
  states_to_basis(1) = StateConversion(ids=[1], coefficients=[1.0_dp])
  do state=2,no_states
    ! Construct h[i,j] from old = h[i-1,{j}].
    old = states_to_basis(state-1)
    
    states_to_basis(state) = StateConversion()
    i  = state-1
    do term=1,size(old)
      ! Add sqrt((2j-1)/i) h[i-1,j-1]
      j = univariates(old%id(term))%power + 1
      coefficient = sqrt((2*j-1.0_dp)/i) * old%coefficient(term)
      call states_to_basis(state)%add_element( id          = ids(j),     &
                                             & coefficient = coefficient )
      
      ! Subtract (j+1)/sqrt(i(2j+1)) h[i-1,j+1]
      j = univariates(old%id(term))%power - 1
      if (j>=0) then
        coefficient = (-(j+1)/sqrt(i*(2*j+1.0_dp))) * old%coefficient(term)
        call states_to_basis(state)%add_element( id          = ids(j),     &
                                               & coefficient = coefficient )
      endif
    enddo
  enddo
  
  ! A monomial state |i> is defined in terms of harmonic basis states |j> as
  !    |i> = sum g[i,j]|j>, where g is the inverse of h above.
  ! g[0,0] = 1.
  ! g[i,j] = sqrt(j/(2i-1))     g[i-1,j-1]
  !        + sqrt((j+1)/(2i-1)) g[i-1,j+1]
  allocate(basis_to_states(no_states), stat=ialloc); call err(ialloc)
  basis_to_states(1) = StateConversion(ids=[1], coefficients=[1.0_dp])
  do state=2,no_states
    ! Construct g[i,j] from g[i-1,{j}].
    old = basis_to_states(state-1)
    
    basis_to_states(state) = StateConversion()
    i  = state-1
    do term=1,size(old)
      ! Add sqrt(j/(2i-1)) g[i-1,j-1]
      j = univariates(old%id(term))%power + 1
      coefficient = sqrt(j/(2*i-1.0_dp)) * old%coefficient(term)
      call basis_to_states(state)%add_element( id          = ids(j),     &
                                             & coefficient = coefficient )
      
      ! Add sqrt((j+1)/(2i-1)) g[i-1,j+1]
      j = univariates(old%id(term))%power - 1
      if (j>=0) then
        coefficient = sqrt((j+1)/(2*i-1.0_dp)) * old%coefficient(term)
        call basis_to_states(state)%add_element( id          = ids(j),     &
                                               & coefficient = coefficient )
      endif
    enddo
  enddo
  
  ! Calculate which states have non-zero <i|H|j> elements.
  ! |i> = a^n|0> and |j> = a^m|0>.
  ! If |n-m|>expansion_order then <i|H|j>=0.
  allocate(harmonic_couplings(no_states), stat=ialloc); call err(ialloc)
  do i=1,no_states
    coupling_ids = [( j,                                         &
                   & j=max(1,i-potential_expansion_order),      &
                   & min(no_states,i+potential_expansion_order) )]
    separations = [(abs(coupling_ids(j)-i), j=1, size(coupling_ids))]
    harmonic_couplings(i) = CoupledStates(coupling_ids,separations)
  enddo
  
  ! Construct output.
  output = WavevectorBasis( maximum_power      = maximum_power,             &
                          & expansion_order    = potential_expansion_order, &
                          & subspace_id        = subspace%id,               &
                          & frequency          = frequency,                 &
                          & wavevector         = fracvec(zeroes(3)),        &
                          & degeneracy         = 1,                         &
                          & monomial_states    = monomial_states,           &
                          & harmonic_states    = harmonic_states,           &
                          & states_to_basis    = states_to_basis,           &
                          & basis_to_states    = basis_to_states,           &
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
  
  integer, allocatable :: ids(:,:) ! N.B. ids is zero-indexed.
  
  type(ComplexUnivariate), allocatable :: univariates(:)
  type(ComplexMonomial),   allocatable :: monomials(:)
  type(MonomialState),     allocatable :: monomial_states(:)
  type(HarmonicState),     allocatable :: harmonic_states(:)
  type(StateConversion),   allocatable :: states_to_basis(:)
  type(StateConversion),   allocatable :: basis_to_states(:)
  type(CoupledStates),     allocatable :: harmonic_couplings(:)
  
  integer               :: no_states
  integer               :: power
  integer               :: state
  integer               :: term
  integer               :: id
  real(dp)              :: coefficient
  type(StateConversion) :: old
  integer               :: separation
  
  integer, allocatable :: coupling_ids(:)
  integer, allocatable :: separations(:)
  
  integer :: ip,im,jp,jm
  
  integer :: i,j,ialloc
  
  no_states = ((1+maximum_power)*(2+maximum_power))/2
  
  ! Generate the monomial and harmonic states.
  allocate( ids(0:maximum_power,0:maximum_power), &
          & univariates(no_states),               &
          & stat=ialloc); call err(ialloc)
  id = 0
  do i=0,maximum_power
    do j=0,i
      ip = j
      im = i-j
      id = id+1
      ! Generate the monomial state |ip,im> = |(u+)^(ip) * (u-)^(im)> at
      !    monomial_states(id), and record the (ip,im)->id mapping.
      ids(ip,im) = id
      univariates(id) = ComplexUnivariate(mode, power=ip, paired_power=im)
    enddo
  enddo
  monomials = [( ComplexMonomial( coefficient = cmplx(1.0_dp,0.0_dp,dp),    &
               &                  modes       = [univariates(i)]         ), &
               & i=1,                                                       &
               & no_states                                                  )]
  monomials(1) = ComplexMonomial( coefficient = cmplx(1.0_dp,0.0_dp,dp), &
                                & modes       = [ComplexUnivariate::]    )
  monomial_states = [( MonomialState(subspace%id, frequency, monomials(i)), &
                     & i=1,                                                 &
                     & no_states                                            )]
  harmonic_states = [( HarmonicState(subspace%id, frequency, monomials(i)), &
                     & i=1,                                                 &
                     & no_states                                            )]
  
  ! A harmonic basis state |ip,im> is defined in terms of monomial states
  !    |jp,jm> as |ip,im> = sum h[ip,im,jp,jm]|jp,jm>.
  ! h[0,0,0,0} = 1.
  ! h[ip,im,jp,jm] = sqrt((jp+jm)/ip)         h[ip-1, im, jp-1, jm  ]
  !                - (jm+1)/sqrt(ip(jp+jm+1)) h[ip-1, im, jp,   jm+1]
  ! also
  ! h[ip,im,jp,jm] = sqrt((jp+jm)/im)          h[ip, im-1, jp,   jm-1]
  !                - (jp+1)/sqrt(im(jp+jm+1))  h[ip, im-1, jp+1, jm  ]
  allocate(states_to_basis(no_states), stat=ialloc); call err(ialloc)
  states_to_basis(1) = StateConversion(ids=[1], coefficients=[1.0_dp])
  do state=2,no_states
    ip = univariates(state)%power
    im = univariates(state)%paired_power
    
    ! If ip>0, construct h[ip,im,jp,jm] from h[ip-1,im,{jp},{jm}].
    ! Otherwise, construct h[ip,im,jp,jm] from h[ip,im-1,{jp},{jm}].
    
    states_to_basis(state) = StateConversion()
    if (ip>0) then
      old = states_to_basis(ids(ip-1,im))
      do term=1,size(old)
        ! Add sqrt((jp+jm)/ip) h[ip-1,im,jp-1,jm]
        jp = univariates(old%id(term))%power + 1
        jm = univariates(old%id(term))%paired_power
        coefficient = sqrt((1.0_dp*(jp+jm))/ip)*old%coefficient(term)
        call states_to_basis(state)%add_element( id          = ids(jp,jm), &
                                               & coefficient = coefficient )
        
        ! Subtract (jm+1)/sqrt(ip(jp+jm+1)) h[ip-1,im,jp,jm+1]
        jp = univariates(old%id(term))%power
        jm = univariates(old%id(term))%paired_power - 1
        if (jm>=0) then
          coefficient = -((jm+1)/sqrt(ip*(jp+jm+1.0_dp)))*old%coefficient(term)
          call states_to_basis(state)%add_element( id          = ids(jp,jm), &
                                                 & coefficient = coefficient )
        endif
      enddo
    else
      old = states_to_basis(ids(ip,im-1))
      do term=1,size(old)
        ! Add sqrt((jp+jm)/im) h[ip,im-1,jp-1,jm-1]
        jp = univariates(old%id(term))%power
        jm = univariates(old%id(term))%paired_power + 1
        coefficient = sqrt((1.0_dp*(jp+jm))/im)*old%coefficient(term)
        call states_to_basis(state)%add_element( id          = ids(jp,jm), &
                                               & coefficient = coefficient )
        
        ! Subtract (jp+1)/sqrt(im(jp+jm+1)) h[ip,im-1,jp+1,jm]
        jp = univariates(old%id(term))%power - 1
        jm = univariates(old%id(term))%paired_power
        if (jp>=0) then
          coefficient = -((jp+1)/sqrt(im*(jp+jm+1.0_dp)))*old%coefficient(term)
          call states_to_basis(state)%add_element( id          = ids(jp,jm), &
                                                 & coefficient = coefficient )
        endif
      enddo
    endif
  enddo
  
  ! A monomial state |ip,im> is defined in terms of harmonic basis states
  !    |jp,jm> as |ip,im> = sum g[ip,im,jp,jm]|jp,jm>,
  !    where g is the inverse of h above.
  ! g[0,0,0,0] = 1.
  ! g[ip,im,jp,jm] = sqrt(jp/(ip+im))    g[ip-1,im,jp-1,jm  ]
  !                + sqrt((jm+1)/(ip+im) g[ip-1,im,jp,  jm+1]
  ! also
  ! g[ip,im,jp,jm] = sqrt(jm/(ip+im))    g[ip-1,im,jp,  jm-1]
  !                + sqrt((jp+1)/(ip+im) g[ip-1,im,jp+1,jm  ]
  allocate(basis_to_states(no_states), stat=ialloc); call err(ialloc)
  basis_to_states(1) = StateConversion(ids=[1], coefficients=[1.0_dp])
  do state=2,no_states
    ip = univariates(state)%power
    im = univariates(state)%paired_power
    
    ! If ip>0, construct g[ip,im,jp,jm] from g[ip-1,im,{jp},{jm}].
    ! Otherwise, construct g[ip,im,jp,jm] from f[ip,im-1,{jp},{jm}].
    
    basis_to_states(state) = StateConversion()
    if (ip>0) then
      old = basis_to_states(ids(ip-1,im))
      do term=1,size(old)
        ! Add sqrt(jp/(ip+im)) g[ip-1,im,jp-1,jm]
        jp = univariates(old%id(term))%power + 1
        jm = univariates(old%id(term))%paired_power
        coefficient = sqrt((1.0_dp*jp)/(ip+im)) * old%coefficient(term)
        call basis_to_states(state)%add_element( id          = ids(jp,jm), &
                                               & coefficient = coefficient )
        
        ! Add sqrt((jm-1)/(ip+im)) g[ip-1,im,jp,jm+1]
        jp = univariates(old%id(term))%power
        jm = univariates(old%id(term))%paired_power - 1
        if (jm>=0) then
          coefficient = sqrt((jm+1.0_dp)/(ip+im)) * old%coefficient(term)
          call basis_to_states(state)%add_element( id          = ids(jp,jm), &
                                                 & coefficient = coefficient )
        endif
      enddo
    else
      old = basis_to_states(ids(ip,im-1))
      do term=1,size(old)
        ! Add sqrt(jm/(ip+im)) g[ip,im-1,jp,jm-1]
        jp = univariates(old%id(term))%power
        jm = univariates(old%id(term))%paired_power + 1
        coefficient = sqrt((1.0_dp*jm)/(ip+im)) * old%coefficient(term)
        call basis_to_states(state)%add_element( id          = ids(jp,jm), &
                                               & coefficient = coefficient )
        
        ! Add sqrt((jp-1)/(ip+im)) g[ip,im-1,jp+1,jm]
        jp = univariates(old%id(term))%power - 1
        jm = univariates(old%id(term))%paired_power
        if (jp>=0) then
          coefficient = sqrt((jp+1.0_dp)/(ip+im)) * old%coefficient(term)
          call basis_to_states(state)%add_element( id          = ids(jp,jm), &
                                                 & coefficient = coefficient )
        endif
      enddo
    endif
  enddo
  
  ! Calculate which states have non-zero <i|H|j> elements.
  ! |i> = (a+)^(m+).(a-)^(m-).|0>
  ! |j> = (a+)^(n+).(a-)^(n-).|0>
  ! If |m+ - n+| + |m- - n-| > potential_expansion_order then <i|H|j>=0.
  allocate(harmonic_couplings(no_states), stat=ialloc); call err(ialloc)
  do i=1,no_states
    coupling_ids = [integer::]
    separations = [integer::]
    do j=1,no_states
      separation = abs(univariates(i)%power-univariates(j)%power) &
               & + abs(univariates(i)%paired_power-univariates(j)%paired_power)
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
                          & degeneracy         = 1,                         &
                          & monomial_states    = monomial_states,           &
                          & harmonic_states    = harmonic_states,           &
                          & states_to_basis    = states_to_basis,           &
                          & basis_to_states    = basis_to_states,           &
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
  
  integer, allocatable :: this_powers(:)
  integer, allocatable :: that_powers(:)
  integer, allocatable :: no_states(:)
  integer, allocatable :: this_to_output(:)
  
  integer  :: id
  real(dp) :: coefficient
  integer  :: separation
  
  ! Output variables.
  integer                            :: maximum_power
  integer                            :: expansion_order
  integer                            :: subspace_id
  real(dp)                           :: frequency
  type(FractionVector)               :: wavevector
  integer                            :: degeneracy
  type(MonomialState),   allocatable :: monomial_states(:)
  type(HarmonicState),   allocatable :: harmonic_states(:)
  type(StateConversion), allocatable :: states_to_basis(:)
  type(StateConversion), allocatable :: basis_to_states(:)
  type(CoupledStates),   allocatable :: harmonic_couplings(:)
  
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
  degeneracy      = this%degeneracy
  
  ! Check that the states in 'that' are in ascending order of total power.
  this_powers = this%monomial_states%total_power()
  that_powers = that%monomial_states%total_power()
  if (any(that_powers(2:)<that_powers(:size(that)-1))) then
    call print_line(CODE_ERROR//': Powers are not in ascending order.')
    call err()
  endif
  
  ! Count the number of output states corresponding to each state in 'this'.
  no_states = [( count(that_powers + this_powers(i) <= maximum_power), &
               & i=1,                                                  &
               & size(this)                                            )]
  this_to_output = [( sum(no_states(:i-1)), i=1, size(this) )]
  
  ! Generate monomial and harmonic states.
  allocate( monomial_states(sum(no_states)), &
          & harmonic_states(sum(no_states)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(this)
    monomial_states(this_to_output(i)+1:this_to_output(i)+no_states(i)) = &
       & this%monomial_states(i) * that%monomial_states(:no_states(i))
    harmonic_states(this_to_output(i)+1:this_to_output(i)+no_states(i)) = &
       & this%harmonic_states(i) * that%harmonic_states(:no_states(i))
  enddo
  
  ! Generate states-to-basis conversion.
  allocate( states_to_basis(size(monomial_states)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(this)
    do j=1,no_states(i)
      states_to_basis(this_to_output(i)+j) = StateConversion()
      do k=1,size(this%states_to_basis_(i))
        do l=1,size(that%states_to_basis_(j))
          id = this_to_output(this%states_to_basis_(i)%id(k)) &
           & + that%states_to_basis_(j)%id(l)
          coefficient = this%states_to_basis_(i)%coefficient(k) &
                    & * that%states_to_basis_(j)%coefficient(l)
          call states_to_basis(this_to_output(i)+j)%add_element( &
                                     & id          = id,         &
                                     & coefficient = coefficient )
        enddo
      enddo
    enddo
  enddo
  
  ! Generate basis-to-states conversion.
  allocate( basis_to_states(size(monomial_states)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(this)
    do j=1,no_states(i)
      basis_to_states(this_to_output(i)+j) = StateConversion()
      do k=1,size(this%basis_to_states_(i))
        do l=1,size(that%basis_to_states_(j))
          id = this_to_output(this%basis_to_states_(i)%id(k)) &
           & + that%basis_to_states_(j)%id(l)
          coefficient = this%basis_to_states_(i)%coefficient(k) &
                    & * that%basis_to_states_(j)%coefficient(l)
          call basis_to_states(this_to_output(i)+j)%add_element( &
                                     & id          = id,         &
                                     & coefficient = coefficient )
        enddo
      enddo
    enddo
  enddo
  
  ! Generate harmonic couplings.
  allocate( harmonic_couplings(size(monomial_states)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(this)
    do j=1,no_states(i)
      harmonic_couplings(this_to_output(i)+j) = CoupledStates()
      do k=1,size(this%harmonic_couplings(i))
        do l=1,size(that%harmonic_couplings(j))
          if ( that%harmonic_couplings(j)%id(l) >          &
             & no_states(this%harmonic_couplings(i)%id(k)) ) then
            exit
          endif
          id = this_to_output(this%harmonic_couplings(i)%id(k)) &
           & + that%harmonic_couplings(j)%id(l)
          separation = this%harmonic_couplings(i)%separation(k) &
                   & + that%harmonic_couplings(j)%separation(l)
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
                          & degeneracy         = degeneracy,        &
                          & monomial_states    = monomial_states,   &
                          & harmonic_states    = harmonic_states,   &
                          & states_to_basis    = states_to_basis,   &
                          & basis_to_states    = basis_to_states,   &
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
  
  type(FractionVector), allocatable :: wavevectors(:)
  
  integer, allocatable :: qpoint_set(:)
  integer, allocatable :: wavevector_states(:)
  integer, allocatable :: new_ids(:)
  integer, allocatable :: ids(:)
  integer, allocatable :: separations(:)
  integer, allocatable :: allowed_couplings(:)
  
  ! Output variables.
  integer                            :: maximum_power
  integer                            :: expansion_order
  integer                            :: subspace_id
  real(dp)                           :: frequency
  type(FractionVector)               :: wavevector
  integer                            :: degeneracy
  type(MonomialState),   allocatable :: monomial_states(:)
  type(HarmonicState),   allocatable :: harmonic_states(:)
  type(StateConversion), allocatable :: states_to_basis(:)
  type(StateConversion), allocatable :: basis_to_states(:)
  type(CoupledStates),   allocatable :: couplings(:)
  
  integer :: i,j,ialloc
  
  wavevectors = [( input%monomial_states(i)%wavevector(modes,qpoints), &
                 & i=1,                                                &
                 & size(input)                                         )]
  
  qpoint_set = filter([( any(wavevectors==qpoints(i)%qpoint), &
                       & i=1,                                 &
                       & size(qpoints)                        )])
  
  allocate(output(size(qpoint_set)), stat=ialloc); call err(ialloc)
  new_ids = [(0,i=1,size(input))]
  do i=1,size(output)
    ! Copy over the maximum power, subspace ID and frequency.
    maximum_power = input%maximum_power
    expansion_order = input%expansion_order
    subspace_id = input%subspace_id
    frequency = input%frequency
    
    ! Identify the states at the given wavevector.
    wavevector = qpoints(qpoint_set(i))%qpoint
    degeneracy = 1
    wavevector_states = filter(wavevectors==wavevector)
    monomial_states = input%monomial_states(wavevector_states)
    harmonic_states = input%harmonic_states(wavevector_states)
    states_to_basis = input%states_to_basis_(wavevector_states)
    basis_to_states = input%basis_to_states_(wavevector_states)
    couplings       = input%harmonic_couplings(wavevector_states)
    
    ! Re-map IDs to correspond to states at the given wavevector only.
    new_ids = 0
    new_ids(wavevector_states) = [(j,j=1,size(wavevector_states))]
    do j=1,size(wavevector_states)
      states_to_basis(j) = StateConversion(                  &
         & ids          = new_ids(states_to_basis(j)%ids()), &
         & coefficients = states_to_basis(j)%coefficients()  )
      basis_to_states(j) = StateConversion(                  &
         & ids          = new_ids(basis_to_states(j)%ids()), &
         & coefficients = basis_to_states(j)%coefficients()  )
      ids = couplings(j)%ids()
      separations = couplings(j)%separations()
      ids = new_ids(ids)
      allowed_couplings = filter(ids>0)
      couplings(j) = CoupledStates( ids(allowed_couplings),        &
                                  & separations(allowed_couplings) )
    enddo
    
    ! Construct output.
    output(i) = WavevectorBasis( maximum_power,   &
                               & expansion_order, &
                               & subspace_id,     &
                               & frequency,       &
                               & wavevector,      &
                               & degeneracy,      &
                               & monomial_states, &
                               & harmonic_states, &
                               & states_to_basis, &
                               & basis_to_states, &
                               & couplings        )
  enddo
end function

! ----------------------------------------------------------------------
! Transform coefficients between monomial states and the orthonormal basis.
! ----------------------------------------------------------------------
function coefficients_states_to_basis(this,coefficients) result(output)
  implicit none
  
  class(WavevectorBasis), intent(in) :: this
  real(dp),               intent(in) :: coefficients(:)
  real(dp), allocatable              :: output(:)
  
  integer :: i
  
  output = [( 0.0_dp, i=1, size(this) )]
  
  do i=1,size(this)
    output(this%basis_to_states_(i)%ids()) =      &
       &   output(this%basis_to_states_(i)%ids()) &
       & + coefficients(i)                        &
       & * this%basis_to_states_(i)%coefficients()
  enddo
end function

function coefficients_basis_to_states(this,coefficients) result(output)
  implicit none
  
  class(WavevectorBasis), intent(in) :: this
  real(dp),               intent(in) :: coefficients(:)
  real(dp), allocatable              :: output(:)
  
  integer :: i
  
  output = [( 0.0_dp, i=1, size(this) )]
  
  do i=1,size(this)
      output(this%states_to_basis_(i)%ids()) =      &
         &   output(this%states_to_basis_(i)%ids()) &
         & + coefficients(i)                        &
         & * this%states_to_basis_(i)%coefficients()
  enddo
end function

! ----------------------------------------------------------------------
! Transform operators between monomial states and the orthonormal basis.
! ----------------------------------------------------------------------
function operator_states_to_basis(this,operator) result(output)
  implicit none
  
  class(WavevectorBasis), intent(in) :: this
  real(dp),               intent(in) :: operator(:,:)
  real(dp), allocatable              :: output(:,:)
  
  integer :: i,j,ialloc
  
  allocate(output(size(this), size(this)), stat=ialloc); call err(ialloc)
  output = 0.0_dp
  do i=1,size(this)
    do j=1,size(this)
      output(j,i) = dble( vec(this%states_to_basis_(j)%coefficients())    &
                      & * mat(operator( this%states_to_basis_(j)%ids(),   &
                      &                 this%states_to_basis_(i)%ids() )) &
                      & * vec(this%states_to_basis_(i)%coefficients())    )
    enddo
  enddo
end function

function operator_basis_to_states(this,operator) result(output)
  implicit none
  
  class(WavevectorBasis), intent(in) :: this
  real(dp),               intent(in) :: operator(:,:)
  real(dp), allocatable              :: output(:,:)
  
  integer :: i,j,ialloc
  
  allocate(output(size(this), size(this)), stat=ialloc); call err(ialloc)
  output = 0.0_dp
  do i=1,size(this)
    do j=1,size(this)
      output(j,i) = dble( vec(this%basis_to_states_(j)%coefficients())    &
                      & * mat(operator( this%basis_to_states_(j)%ids(),   &
                      &                 this%basis_to_states_(i)%ids() )) &
                      & * vec(this%basis_to_states_(i)%coefficients())    )
    enddo
  enddo
end function

! ----------------------------------------------------------------------
! As operator transformations, but uses symmetry to accelerate calculation.
! ----------------------------------------------------------------------
function hamiltonian_states_to_basis(this,hamiltonian) result(output)
  implicit none
  
  class(WavevectorBasis), intent(in) :: this
  real(dp),               intent(in) :: hamiltonian(:,:)
  real(dp), allocatable              :: output(:,:)
  
  integer :: i,j,ialloc
  
  allocate(output(size(this), size(this)), stat=ialloc); call err(ialloc)
  output = 0.0_dp
  do i=1,size(this)
    do j=1,size(this)
      output(j,i) = dble( vec(this%states_to_basis_(j)%coefficients())       &
                      & * mat(hamiltonian( this%states_to_basis_(j)%ids(),   &
                      &                    this%states_to_basis_(i)%ids() )) &
                      & * vec(this%states_to_basis_(i)%coefficients())       )
    enddo
  enddo
end function

function hamiltonian_basis_to_states(this,hamiltonian) result(output)
  implicit none
  
  class(WavevectorBasis), intent(in) :: this
  real(dp),               intent(in) :: hamiltonian(:,:)
  real(dp), allocatable              :: output(:,:)
  
  integer :: i,j,ialloc
  
  allocate(output(size(this), size(this)), stat=ialloc); call err(ialloc)
  output = 0.0_dp
  do i=1,size(this)
    do j=1,size(this)
      output(j,i) = dble( vec(this%basis_to_states_(j)%coefficients())       &
                      & * mat(hamiltonian( this%basis_to_states_(j)%ids(),   &
                      &                    this%basis_to_states_(i)%ids() )) &
                      & * vec(this%basis_to_states_(i)%coefficients())       )
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
  
  integer                            :: maximum_power
  integer                            :: expansion_order
  integer                            :: subspace_id
  real(dp)                           :: frequency
  type(FractionVector)               :: wavevector
  integer                            :: degeneracy
  type(MonomialState),   allocatable :: monomial_states(:)
  type(HarmonicState),   allocatable :: harmonic_states(:)
  type(StateConversion), allocatable :: states_to_basis(:)
  type(StateConversion), allocatable :: basis_to_states(:)
  type(CoupledStates),   allocatable :: couplings(:)
  
  type(String), allocatable :: line(:)
  
  integer :: no_states
  
  integer :: i,ialloc
  
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
    
    line = split_line(input(6))
    degeneracy = int(line(3))
    
    no_states = (size(input)-11)/5
    
    allocate( monomial_states(no_states), &
            & harmonic_states(no_states), &
            & states_to_basis(no_states), &
            & basis_to_states(no_states), &
            & couplings(no_states),       &
            & stat=ialloc); call err(ialloc)
    do i=1,no_states
      line = split_line(input(7+i))
      monomial_states(i) = MonomialState([input(3:4), str('State'), line(3)])
      
      line = split_line(input(8+no_states+i))
      harmonic_states(i) = HarmonicState([input(3:4), str('State'), line(3)])
      
      line = split_line(input(9+2*no_states+i))
      states_to_basis(i) = StateConversion(join(line(3:)))
      
      line = split_line(input(10+3*no_states+i))
      basis_to_states(i) = StateConversion(join(line(3:)))
      
      couplings(i) = CoupledStates(input(11+4*no_states+i))
    enddo
    
    this = WavevectorBasis( maximum_power,   &
                          & expansion_order, &
                          & subspace_id,     &
                          & frequency,       &
                          & wavevector,      &
                          & degeneracy,      &
                          & monomial_states, &
                          & harmonic_states, &
                          & states_to_basis, &
                          & basis_to_states, &
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
             & 'Degeneracy      : '//this%degeneracy,      &
             & str('Monomial states')                      ]
    do i=1,size(this%monomial_states)
      state_strings = str(this%monomial_states(i))
      output = [output, '|'//i//'> = '//state_strings(size(state_strings))]
    enddo
    
    output = [ output,                &
             & str('Harmonic states') ]
    do i=1,size(this%harmonic_states)
      state_strings = str(this%harmonic_states(i))
      output = [output, '|'//i//'> = '//state_strings(size(state_strings))]
    enddo
    
    output = [ output,                                   &
             & str('Harmonic states in monomial basis.') ]
    do i=1,size(this%states_to_basis_)
      output = [ output,                                       &
               & '|'//i//'> = '//str(this%states_to_basis_(i)) ] 
    enddo
    
    output = [ output,                                  &
             & str('Monomial states in harmonic Basis') ]
    do i=1,size(this%basis_to_states_)
      output = [ output,                                       &
               & '|'//i//'> = '//str(this%basis_to_states_(i)) ] 
    enddo
    
    output = [ output,                                              &
             & str('Non-zero <i|H|j> integrals in harmonic basis'), &
             & str(this%harmonic_couplings)                         ]
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
