! ======================================================================
! A state along each complex mode in a degenerate subspace,
!    defined as propotional to a monomial times |0>, where |0> is the
!    harmonic ground state.
! ======================================================================
! N.B. compare with HarmonicState:
!    - monomial states are normalised but not orthogonal.
!    - monomial states are easier to evaluate than harmonic states.
module monomial_state_module
  use common_module
  
  use anharmonic_common_module
  
  use state_helper_module
  implicit none
  
  private
  
  public :: startup_monomial_state
  
  public :: MonomialState
  
  public :: operator(*)
  public :: generate_monomial_states
  
  ! A MonomialState is a product of single-mode and double-mode states.
  ! |n_1,n_2,n_3,...,n_> = |n_1>|n_2,n_3>...|n_>.
  !
  ! If a mode is its own conjugate, (u_i)* = u_i,
  !    then the single-mode states along mode u_i are:
  ! |n_i> = prod_{k=1}^{n_i}[ sqrt(2Nw/(2k-1)) ] (u_i)^(n_i) |0_i>
  ! |0_i> = sqrt(sqrt(m*w/pi)) exp(- 1/2 N w (u_i)^2 )
  !
  ! If a mode is not its own conjugate, (u_i)* = u_j,
  !    then the double-mode states along modes u_i and u_j are:
  ! |n_i,n_j> = prod_{k=1}^{n_i+n_j}[ sqrt(4Nw/k) ] ui^ni uj^nj |0_i,0_j>
  ! |0_i,0_j> = sqrt(2*m*w/pi) exp(- N w |u_i|^2 )
  !
  ! In both cases, states are normalised (<n_i|n_i>=1, <n_i,n_j|n_i,n_j>=1),
  !    but are in general not orthogonal (<p_i|q_i>/=0).
  !
  ! w is the effective frequency of all modes in the subspace
  !    (in general not the same as the harmonic frequency).
  !
  ! m is the geometric average mass, arising as a result of mass-reduction
  !    of co-ordinates. This cancels in integrals, and so is not stored.
  !
  ! N is the number of primitive cells in the anharmonic supercell.
  !
  ! N.B. the coefficient of state_ is not used. Instead, states are implicitly
  !    normalised such that <state|state>=1. This removes the need to keep
  !    track of the various factors of two, pi, m (the geometric mean of the
  !    atomic masses, arising as a result of mass reduction) and N (the number
  !    of primitive cells in the supercell) , and allows the frequency to be
  !    changed without having to re-calculate coefficients.
  type, extends(SubspaceState) :: MonomialState
    real(dp)                       :: frequency
    type(ComplexMonomial), private :: state_
  contains
    procedure, public, nopass :: representation => representation_MonomialState
    
    procedure, public :: inner_product => &
                       & inner_product_MonomialState
    
    procedure, public :: braket_ComplexMonomial => &
                       & braket_ComplexMonomial_MonomialState
    procedure, public :: kinetic_energy => &
                       & kinetic_energy_MonomialState
    procedure, public :: harmonic_potential_energy => &
                       & harmonic_potential_energy_MonomialState
    procedure, public :: kinetic_stress => &
                       & kinetic_stress_MonomialState
    
    procedure, public :: change_modes => change_modes_MonomialState
    
    procedure, public :: total_power => total_power_MonomialState
    procedure, public :: wavevector => wavevector_MonomialState
    
    procedure, public :: wavefunction => wavefunction_MonomialState
    ! I/O.
    procedure, public :: read  => read_MonomialState
    procedure, public :: write => write_MonomialState
  end type
  
  interface MonomialState
    module procedure new_MonomialState
    module procedure new_MonomialState_SubspaceState
    module procedure new_MonomialState_Strings
    module procedure new_MonomialState_StringArray
  end interface
  
  interface operator(*)
    module procedure multiply_MonomialState_MonomialState
  end interface
  
  interface finite_overlap
    module procedure finite_overlap_MonomialStates
  end interface
contains

! Startup procedure.
subroutine startup_monomial_state()
  implicit none
  
  type(MonomialState) :: state
  
  call state%startup()
end subroutine

! ----------------------------------------------------------------------
! Constructor.
! ----------------------------------------------------------------------
function new_MonomialState(subspace_id,frequency,state) result(this)
  implicit none
  
  integer,               intent(in) :: subspace_id
  real(dp),              intent(in) :: frequency
  type(ComplexMonomial), intent(in) :: state
  type(MonomialState)               :: this
  
  this%subspace_id = subspace_id
  this%frequency   = frequency
  this%state_      = state
end function

recursive function new_MonomialState_SubspaceState(input) result(this)
  implicit none
  
  class(SubspaceState), intent(in) :: input
  type(MonomialState)              :: this
  
  select type(input); type is(MonomialState)
    this = input
  type is(SubspaceStatePointer)
    this = MonomialState(input%state())
  class default
    call err()
  end select
end function

! ----------------------------------------------------------------------
! Type representation.
! ----------------------------------------------------------------------
impure elemental function representation_MonomialState() result(output)
  implicit none
  
  type(String) :: output
  
  output = 'monomial'
end function

! ----------------------------------------------------------------------
! Multiply two monomial states together.
! ----------------------------------------------------------------------
impure elemental function multiply_MonomialState_MonomialState(this,that) &
   & result(output)
  implicit none
  
  type(MonomialState), intent(in) :: this
  type(MonomialState), intent(in) :: that
  type(MonomialState)             :: output
  
  if (this%subspace_id/=that%subspace_id) then
    call print_line(CODE_ERROR//': subspaces do not match.')
    call err()
  endif
  
  output = MonomialState( this%subspace_id,       &
                        & this%frequency,         &
                        & this%state_*that%state_ )
end function

! ----------------------------------------------------------------------
! Generates all monomial states in a subspace up to a given power.
! ----------------------------------------------------------------------
function generate_monomial_states(subspace,frequency,modes,maximum_power) &
   & result(output)
  implicit none
  
  type(DegenerateSubspace), intent(in) :: subspace
  real(dp),                 intent(in) :: frequency
  type(ComplexMode),        intent(in) :: modes(:)
  integer,                  intent(in) :: maximum_power
  type(MonomialState), allocatable     :: output(:)
  
  type(ComplexMode), allocatable :: subspace_modes(:)
  type(ComplexMonomial)          :: monomial
  type(MonomialState)            :: state
  
  subspace_modes = subspace%modes(modes)
  monomial = ComplexMonomial( coefficient = cmplx(1.0_dp,0.0_dp,dp), &
                            & modes       = [ComplexUnivariate::]    )
  state = MonomialState( subspace_id = subspace%id, &
                       & frequency   = frequency,   &
                       & state       = monomial     )
  output = generate_monomial_states_helper(subspace_modes,maximum_power,state)
end function

recursive function generate_monomial_states_helper(modes,power,state) &
   & result(output)
  implicit none
  
  type(ComplexMode),   intent(in)  :: modes(:)
  integer,             intent(in)  :: power
  type(MonomialState), intent(in)  :: state
  type(MonomialState), allocatable :: output(:)
  
  type(MonomialState) :: output_state
  
  integer :: i
  
  if (size(modes)==0 .or. power==0) then
    output = [state]
  else
    output = [ generate_monomial_states_helper( modes(2:),   &
             &                                  power,       &
             &                                  state      ) ]
    do i=1,power
      output_state = state
      output_state%state_ = output_state%state_ &
                        & * ComplexUnivariate(mode=modes(1), power=i)
      output = [ output,                                         &
               & generate_monomial_states_helper( modes(2:),     &
               &                                  power-i,       &
               &                                  output_state ) ]
    enddo
  endif
end function

! ----------------------------------------------------------------------
! Returns the total power of a given state.
! ----------------------------------------------------------------------
! The total power of the state |product_{q,i} (u_{q,i})^(n_{q,i})> is equal to
!    sum_{q,i}} n_{q,i}.
impure elemental function total_power_MonomialState(this) result(output)
  implicit none
  
  class(MonomialState), intent(in) :: this
  integer                          :: output
  
  output = this%state_%total_power()
end function

! ----------------------------------------------------------------------
! Returns the wavevector of a given state.
! ----------------------------------------------------------------------
! The wavevector of the state |product_{q,i} (u_{q,i})^(n_{q,i})> is equal to
!    sum_{q,i}} n_{q,i}q.
function wavevector_MonomialState(this,modes,qpoints) result(output)
  implicit none
  
  class(MonomialState), intent(in) :: this
  type(ComplexMode),    intent(in) :: modes(:)
  type(QpointData),     intent(in) :: qpoints(:)
  type(QpointData)                 :: output
  
  output = this%state_%wavevector(modes,qpoints)
end function

! ----------------------------------------------------------------------
! Returns the wavefunction of the state,
!    with all coefficients accounted for.
! ----------------------------------------------------------------------
impure elemental function wavefunction_MonomialState(this,frequency, &
   & supercell) result(output)
  implicit none
  
  class(MonomialState), intent(in) :: this
  real(dp),             intent(in) :: frequency
  type(StructureData),  intent(in) :: supercell
  type(String)                     :: output
  
  real(dp) :: coefficient
  
  type(ComplexUnivariate) :: mode
  
  integer :: i,k
  
  coefficient = 1
  do i=1,size(this%state_)
    mode = this%state_%mode(i)
    if (mode%id==mode%paired_id) then
      ! |n_i> = prod_{k=1}^{n_i}[ sqrt(2Nw/(2k-1)) ] (u_i)^(n_i) |0_i>
      do k=1,mode%power
        coefficient = coefficient * sqrt(2*supercell%sc_size*frequency/(2*k-1))
      enddo
    else
      ! |n_i,n_j> = prod_{k=1}^{n_i+n_j}[ sqrt(4Nw/k) ] ui^ni uj^nj |0_i,0_j>
      do k=1,mode%power+mode%paired_power
        coefficient = coefficient * sqrt(4*supercell%sc_size*frequency/k)
      enddo
    endif
  enddo
  
  if (size(this%state_)==0) then
    output = str(coefficient)//'|0>'
  else
    output = coefficient//'*'//join(this%state_%modes(),delimiter='*')//'|0>'
  endif
end function

! ----------------------------------------------------------------------
! SubspaceState methods.
! ----------------------------------------------------------------------
impure elemental function inner_product_MonomialState(this, &
   & ket,subspace,subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(MonomialState),     intent(in)           :: this
  class(SubspaceState),     intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  real(dp)                                       :: output
  
  type(MonomialState) :: monomial_ket
  type(StateHelper)   :: helper
  
  if (present(ket)) then
    monomial_ket = MonomialState(ket)
    helper = StateHelper(this%state_, monomial_ket%state_)
  else
    helper = StateHelper(this%state_, this%state_)
  endif
  
  output = product(braket_mode(helper%bra, helper%ket))
end function

impure elemental function braket_ComplexMonomial_MonomialState(this,monomial, &
   & ket,subspace,subspace_basis,anharmonic_data,qpoint) result(output)
  implicit none
  
  class(MonomialState),     intent(in)           :: this
  type(ComplexMonomial),    intent(in)           :: monomial
  class(SubspaceState),     intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(QpointData),         intent(in), optional :: qpoint
  type(ComplexMonomial)                          :: output
  
  type(MonomialState) :: monomial_ket
  type(StateHelper)   :: helper
  
  type(ComplexUnivariate), allocatable :: monomial_modes(:)
  complex(dp)                          :: coefficient
  
  integer :: i
  
  if (present(ket)) then
    ! Process the bra, ket and monomial.
    monomial_ket = MonomialState(ket)
    helper = StateHelper( this%state_,         &
                        & monomial_ket%state_, &
                        & monomial,            &
                        & subspace,            &
                        & anharmonic_data,     &
                        & qpoint               )
  else
    ! Process the bra, ket and monomial.
    helper = StateHelper( this%state_,     &
                        & this%state_,     &
                        & monomial,        &
                        & subspace,        &
                        & anharmonic_data, &
                        & qpoint           )
  endif
  
  ! Extract the modes in the monomial which have not been integrated.
  monomial_modes = monomial%modes(filter(.not.helper%monomial_integrated))
  
  ! Calculate the coefficient of the output.
  ! This is the input coefficient times the integral over all modes in the
  !    subspace.
  ! The helper function calculates this integral for each single- or
  !    double- mode, up to a factor of 1/(2Nw)^(n/2), where
  !    - N is the number of primitive cells in the anharmonic supercell.
  !    - w is the frequency of the modes in the subspace.
  !    - n is the power of the modes in the monomial which are integrated.
  coefficient = monomial%coefficient                                 &
            & * product(braket_mode_potential( helper%bra,           &
            &                                  helper%ket,           &
            &                                  helper%monomial ))    &
            & / sqrt( 2.0_dp                                         &
            &       * anharmonic_data%anharmonic_supercell%sc_size   &
            &       * this%frequency                               ) &
            & **(monomial%total_power()-sum(monomial_modes%total_power()))
  
  ! Construct output.
  output = ComplexMonomial( coefficient = coefficient,   &
                          & modes       = monomial_modes )
end function

impure elemental function kinetic_energy_MonomialState(this,ket, &
   & subspace,subspace_basis,anharmonic_data,qpoint) result(output)
  implicit none
  
  class(MonomialState),     intent(in)           :: this
  class(SubspaceState),     intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(QpointData),         intent(in), optional :: qpoint
  real(dp)                                       :: output
  
  type(MonomialState) :: monomial_ket
  type(StateHelper)   :: helper
  
  real(dp) :: prefactor
  
  integer :: i
  
  if (present(ket)) then
    ! Process the bra and ket.
    monomial_ket = MonomialState(ket)
    helper = StateHelper( this%state_,                       &
                        & monomial_ket%state_,               &
                        & subspace        = subspace,        &
                        & anharmonic_data = anharmonic_data, &
                        & qpoint          = qpoint           )
  else
    ! Process the bra and ket.
    helper = StateHelper( this%state_,                       &
                        & this%state_,                       &
                        & subspace        = subspace,        &
                        & anharmonic_data = anharmonic_data, &
                        & qpoint          = qpoint           )
  endif
  
  if (all(finite_overlap_mode(helper%bra, helper%ket))) then
    ! Calculate the prefactor, s.t. <bra|T|ket> = prefactor * w * <bra|ket>.
    ! <p1|<p2|...T|q1>|q2>... = <p1|T1|p1><p2|q2>... + <p1|q1><p2|T|q2>... + ...
    !                         = (t1+t2+...)*w*<p1|q1><p2|q2>...
    prefactor = sum(kinetic_energy_prefactor(helper%bra, helper%ket))
    ! Calculate output, multiply by w,
    !    and normalise by the number of primitive cells.
    output = prefactor                                                   &
         & * inner_product(this,subspace,subspace_basis,anharmonic_data) &
         & * this%frequency                                              &
         & / anharmonic_data%anharmonic_supercell%sc_size
  else
    ! If more than one single- or double-mode state has zero overlap,
    !    then the kinetic energy is zero.
    output = 0.0_dp
  endif
end function

impure elemental function harmonic_potential_energy_MonomialState( &
   & this,ket,subspace,subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(MonomialState),     intent(in)           :: this
  class(SubspaceState),     intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  real(dp)                                       :: output
  
  type(MonomialState) :: monomial_ket
  type(StateHelper)   :: helper
  
  real(dp) :: prefactor
  
  integer :: i
  
  if (present(ket)) then
    ! Process the bra and ket.
    monomial_ket = MonomialState(ket)
    helper = StateHelper( this%state_,                      &
                        & monomial_ket%state_,              &
                        & subspace        = subspace,       &
                        & anharmonic_data = anharmonic_data )
  else
    ! Process the bra and ket.
    helper = StateHelper( this%state_,                      &
                        & this%state_,                      &
                        & subspace = subspace,              &
                        & anharmonic_data = anharmonic_data )
  endif
  
  ! Calculate the prefactor, s.t. <bra|T|ket> = prefactor * <bra|ket>.
  prefactor = sum(harmonic_potential_energy_prefactor( helper%bra,    &
                                                     & helper%ket,    &
                                                     & this%frequency ))
  
  ! Calculate output, and normalise by the number of primitive cells.
  output = prefactor                                                   &
       & * inner_product(this,subspace,subspace_basis,anharmonic_data) &
       & / anharmonic_data%anharmonic_supercell%sc_size
end function

impure elemental function kinetic_stress_MonomialState(this,ket, &
   & subspace,subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(MonomialState),     intent(in)           :: this
  class(SubspaceState),     intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(RealMatrix)                               :: output
  
  ! TODO
end function

! ----------------------------------------------------------------------
! Change the modes of the state by the specified group.
! ----------------------------------------------------------------------
impure elemental function change_modes_MonomialState(this,mode_group) &
   & result(output)
  implicit none
  
  class(MonomialState), intent(in) :: this
  type(Group),          intent(in) :: mode_group
  type(MonomialState)              :: output
  
  type(ComplexUnivariate), allocatable :: univariates(:)
  type(ComplexMonomial)                :: monomial
  
  integer :: i,ialloc
  
  allocate(univariates(size(this%state_)), stat=ialloc); call err(ialloc)
  do i=1,size(this%state_)
    univariates(i) = ComplexUnivariate(         &
       & mode_group * this%state_%id(i),        &
       & mode_group * this%state_%paired_id(i), &
       & this%state_%power(i),                  &
       & this%state_%paired_power(i)            )
  enddo
  monomial = ComplexMonomial(this%state_%coefficient, univariates)
  
  output = this
  output%state_ = monomial
end function

! ----------------------------------------------------------------------
! Helpers for braket, kinetic_energy and potential_energy.
! ----------------------------------------------------------------------
! Returns whether or not braket(bra,ket) is non-zero.
impure elemental function finite_overlap_MonomialStates(bra,ket) result(output)
  implicit none
  
  type(MonomialState),      intent(in) :: bra
  type(MonomialState),      intent(in) :: ket
  logical                              :: output
  
  type(StateHelper) :: helper
  
  integer :: i
  
  ! Check that the bra and the ket cover the same subspace.
  if (bra%subspace_id/=ket%subspace_id) then
    call print_line(ERROR//': bra and ket from different subspaces.')
    call err()
  endif
  
  helper = StateHelper(bra%state_, ket%state_)
  do i=1,size(helper)
    if (.not. finite_overlap_mode(helper%bra(i),helper%ket(i))) then
      output = .false.
      return
    endif
  enddo
  output = .true.
end function

impure elemental function finite_overlap_mode(bra,ket) result(output)
  implicit none
  
  type(ComplexUnivariate), intent(in) :: bra
  type(ComplexUnivariate), intent(in) :: ket
  logical                             :: output
  
  integer :: p_i,p_j,q_i,q_j
  
  p_i = bra%power
  p_j = bra%paired_power
  q_i = ket%power
  q_j = ket%paired_power
  
  if (bra%id==bra%paired_id) then
    ! <p_i|q_i>  = 0 if p+q odd.
    !           /= 0 otherwise.
    output = modulo(p_i+q_i,2)==0
  else
    ! <p_i,p_j|q_i,q_j>  = 0 if p_i-p_j-q_i+q_j /= 0.
    !                   /= 0 otherwise.
    output = p_i-p_j-q_i+q_j==0
  endif
end function

impure elemental function braket_mode(bra,ket) result(output)
  implicit none
  
  type(ComplexUnivariate), intent(in) :: bra
  type(ComplexUnivariate), intent(in) :: ket
  real(dp)                            :: output
  
  integer :: p_i,p_j,q_i,q_j
  
  p_i = bra%power
  p_j = bra%paired_power
  q_i = ket%power
  q_j = ket%paired_power
  
  if (bra%id==bra%paired_id) then
    ! <p_i|q_i> = 0                                      if p+q odd.
    !           = f((p_i+q_i)/2) / sqrt(f(p_i)*f(q_i))   otherwise,
    ! where f(x) is the odd factorial of x, f(x) = prod_{k=1}^x [ 2k-1 ].
    if (modulo(p_i+q_i,2)==1) then
      output = 0.0_dp
    else
      output = exp( log_odd_factorial((p_i+q_i)/2)                         &
                & - 0.5_dp*(log_odd_factorial(p_i)+log_odd_factorial(q_i)) )
    endif
  else
    ! <p_i,p_j|q_i,q_j> = 0                          if p_i-p_j-q_i+q_j /= 0.
    !                   = ((p_i+p_j+q_i+q_j)/2)!
    !                   / sqrt((p_i+p_j)!(q_i+q_j)!) otherwise.
    if (p_i-p_j-q_i+q_j/=0) then
      output = 0.0_dp
    else
      output = exp( log_factorial((p_i+p_j+q_i+q_j)/2)                     &
                & - 0.5_dp*(log_factorial(p_i+p_j)+log_factorial(q_i+q_j)) )
    endif
  endif
end function

impure elemental function braket_mode_potential(bra,ket,potential) &
   & result(output)
  implicit none
  
  type(ComplexUnivariate), intent(in) :: bra
  type(ComplexUnivariate), intent(in) :: ket
  type(ComplexUnivariate), intent(in) :: potential
  real(dp)                            :: output
  
  integer :: p_i,p_j,q_i,q_j,n_i,n_j
  
  integer :: k
  
  p_i = bra%power
  p_j = bra%paired_power
  q_i = ket%power
  q_j = ket%paired_power
  n_i = potential%power
  n_j = potential%paired_power
  
  if (bra%id==bra%paired_id) then
    ! <p_i|(u_i)^(n_i)|q_i> = 0                     if p+n+q odd.
    !                       = 1/sqrt(2Nw)^{n_i}
    !                       * f((p_i+n_i+q_i)/2)
    !                       / sqrt(f(p_i)*f(q_i))   otherwise,
    ! where f(x) is the odd factorial of x, f(x) = prod_{k=1}^x [ 2k-1 ].
    !
    ! N.B. the factor of 1/sqrt(2Nw)^{n_i} is neglected.
    if (modulo(p_i+n_i+q_i,2)==1) then
      output = 0
    else
      output = exp( log_odd_factorial((p_i+n_i+q_i)/2)                     &
                & - 0.5_dp*(log_odd_factorial(p_i)+log_odd_factorial(q_i)) )
    endif
  else
    ! <p_i,p_j|(u_i)^(n_i) (u_j)^(n_j)|q_i,q_j> =
    !    = 0                                 if p_i-p_j-n_i+n_j-q_i+q_j /= 0.
    !
    !    = 1/sqrt(2Nw)^{n_i+n_j}
    !    * ((p_i+p_j+n_i+n_j+q_i+q_j)/2)!
    !    / sqrt( (p_i+p_j)! * (q_i+q_j)! )   otherwise.
    !
    ! N.B. the factor of 1/sqrt(2Nw)^{n_i+n_j} is neglected.
    if (p_i-p_j-n_i+n_j-q_i+q_j/=0) then
      output = 0
    else
      output = exp( log_factorial((p_i+p_j+n_i+n_j+q_i+q_j)/2)             &
                & - 0.5_dp*(log_factorial(p_i+p_j)+log_factorial(q_i+q_j)) )
    endif
  endif
end function

impure elemental function kinetic_energy_prefactor(bra,ket) result(output)
  implicit none
  
  type(ComplexUnivariate), intent(in) :: bra
  type(ComplexUnivariate), intent(in) :: ket
  real(dp)                            :: output
  
  integer :: p_i,p_j,q_i,q_j
  
  p_i = bra%power
  p_j = bra%paired_power
  q_i = ket%power
  q_j = ket%paired_power
  
  if (bra%id==bra%paired_id) then
    ! <p_i|T|q_i> = 1/2 w_i ( (1-p_i-q_i)/2 - 2p_iq_i/(p_i+q_i+1) ) <p_i|q_i>
    output = 0.5_dp * ( (1-p_i-q_i)/2.0_dp + (2.0_dp*p_i*q_i)/(p_i+q_i-1) )
  else
    ! <p_i,p_j|T|q_i,q_j> = w_i ( 1/2
    !                           + ((p_i-p_j)^2-(p_i-q_i)^2)
    !                           / (p_i+p_j+q_i+q_j)         )
    !                     * <p_i,p_j|q_i,q_j>
    !
    ! N.B. if p_i=p_j=q_i=q_j=0 then the second term is zero.
    output = 0.5_dp
    if (p_i+p_j+q_i+q_j/=0) then
      output = output                        &
           & + ((p_i-p_j)**2 - (p_i-q_i)**2) &
           & / (1.0_dp*(p_i+p_j+q_i+q_j))
    endif
  endif
end function

impure elemental function harmonic_potential_energy_prefactor(bra,ket, &
   & frequency) result(output)
  implicit none
  
  type(ComplexUnivariate), intent(in) :: bra
  type(ComplexUnivariate), intent(in) :: ket
  real(dp),                intent(in) :: frequency
  real(dp)                            :: output
  
  integer :: p_i,p_j,q_i,q_j
  
  p_i = bra%power
  p_j = bra%paired_power
  q_i = ket%power
  q_j = ket%paired_power
  
  if (bra%id==bra%paired_id) then
    ! <p_i|V|q_i> = 1/4 w_i ( 1 + p_i + q_i ) <p_i|q_i>
    if (modulo(p_i+q_i,2)==0) then
      output = 0.25_dp * frequency * (1 + p_i + q_i)
    else
      output = 0.0_dp
    endif
  else
    ! <p_i,p_j|V|q_i,q_j> = 1/4 w_i ( 2 + (p_i+p_j+q_i+q_j) ) <p_i,p_j|q_i,q_j>
    if (p_i+q_j==p_j+q_i) then
      output = 0.25_dp * frequency * (2+p_i+p_j+q_i+q_j)
    else
      output = 0.0_dp
    endif
  endif
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_MonomialState(this,input)
  implicit none
  
  class(MonomialState), intent(out) :: this
  type(String),         intent(in)  :: input(:)
  
  integer               :: subspace_id
  real(dp)              :: frequency
  type(ComplexMonomial) :: state
  
  type(String), allocatable :: line(:)
  
  select type(this); type is(MonomialState)
    line = split_line(input(1))
    subspace_id = int(line(3))
    
    line = split_line(input(2))
    frequency = dble(line(3))
    
    if (input(4)=='|0>') then
      state = ComplexMonomial( coefficient = cmplx(1.0_dp,0.0_dp,dp), &
                             & modes       = [ComplexUnivariate::]    )
    else
      state = ComplexMonomial(1.0_dp//'*'//slice(input(4),2,len(input(4))-1))
    endif
    
    this = MonomialState(subspace_id,frequency,state)
  class default
    call err()
  end select
end subroutine

function write_MonomialState(this) result(output)
  implicit none
  
  class(MonomialState), intent(in) :: this
  type(String), allocatable        :: output(:)
  
  type(String) :: state_string
  
  select type(this); type is(MonomialState)
    if (this%state_%total_power()==0) then
      state_string = '|0>'
    else
      state_string = '|'//join(this%state_%modes(),delimiter='*')//'>'
    endif
    output = [ 'Subspace  : '//this%subspace_id, &
             & 'Frequency : '//this%frequency,   &
             & str('State'),                     &
             & state_string                      ]
  class default
    call err()
  end select
end function

function new_MonomialState_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(MonomialState)      :: this
  
  call this%read(input)
end function

impure elemental function new_MonomialState_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(MonomialState)           :: this
  
  this = MonomialState(str(input))
end function
end module
