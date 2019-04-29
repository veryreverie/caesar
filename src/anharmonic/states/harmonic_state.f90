! ======================================================================
! A state along each complex mode in a degenerate subspace,
!    defined as a (normalised) product of harmonic creation operators
!    acting on |0>, where |0> is the harmonic ground state.
! ======================================================================
! N.B. compare with MonomialState:
!    - harmonic states are orthonormal.
!    - harmonic states are harder to evaluate than monomial states.
module harmonic_state_module
  use common_module
  
  use anharmonic_common_module
  
  use state_helper_module
  implicit none
  
  private
  
  public :: startup_harmonic_state
  
  public :: HarmonicState
  
  public :: operator(*)
  public :: generate_harmonic_states
  
  ! A HarmonicState is a product of single-mode and double-mode states
  ! A state |n_1,n_2,n_3,...,n_> = |n_1>|n_2,n_3>...|n_>.
  !
  ! If a mode is its own conjugate, (u_i)* = u_i,
  !    then the single-mode states along mode u_i are:
  ! |n_i> = 1/sqrt(n_i!) (a'_i)^(n_i) |0_i>
  ! |0_i> = sqrt(sqrt(m*w/pi)) exp(- 1/2 N w (u_i)^2 )
  !
  ! If a mode is not its own conjugate, (u_i)* = u_j,
  !    then the double-mode states along modes u_i and u_j are:
  ! |n_i,n_j> = 1/sqrt(n_i!n_j!) (a'_i)^(n_i) (a'_j)^(n_j) |0_i,0_j>
  ! |0_i,0_j> = sqrt(2*m*w/pi) exp(- N w |u_i|^2 )
  !
  ! where a'_i is the creation operator along mode i.
  !
  ! In both cases, states are orthonormal.
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
  type, extends(SubspaceState) :: HarmonicState
    real(dp)                       :: frequency
    type(ComplexMonomial), private :: state_
  contains
    procedure, public, nopass :: representation => &
                               & representation_HarmonicState
    
    procedure, public :: inner_product => &
                       & inner_product_HarmonicState
    
    procedure, public :: braket_ComplexMonomial => &
                       & braket_ComplexMonomial_HarmonicState
    procedure, public :: kinetic_energy => &
                       & kinetic_energy_HarmonicState
    procedure, public :: harmonic_potential_energy => &
                       & harmonic_potential_energy_HarmonicState
    procedure, public :: kinetic_stress => &
                       & kinetic_stress_HarmonicState
    
    procedure, public :: total_occupation => total_occupation_HarmonicState
    procedure, public :: wavevector => wavevector_HarmonicState
    ! I/O.
    procedure, public :: read  => read_HarmonicState
    procedure, public :: write => write_HarmonicState
  end type
  
  interface HarmonicState
    module procedure new_HarmonicState
    module procedure new_HarmonicState_SubspaceState
    module procedure new_HarmonicState_Strings
    module procedure new_HarmonicState_StringArray
  end interface
  
  interface operator(*)
    module procedure multiply_HarmonicState_HarmonicState
  end interface
  
  interface finite_overlap
    module procedure finite_overlap_HarmonicStates
  end interface
contains

! Startup procedure.
subroutine startup_harmonic_state()
  implicit none
  
  type(HarmonicState) :: state
  
  call state%startup()
end subroutine

! ----------------------------------------------------------------------
! Constructor.
! ----------------------------------------------------------------------
function new_HarmonicState(subspace_id,frequency,state) result(this)
  implicit none
  
  integer,               intent(in) :: subspace_id
  real(dp),              intent(in) :: frequency
  type(ComplexMonomial), intent(in) :: state
  type(HarmonicState)               :: this
  
  this%subspace_id = subspace_id
  this%frequency   = frequency
  this%state_      = state
end function

recursive function new_HarmonicState_SubspaceState(input) result(this)
  implicit none
  
  class(SubspaceState), intent(in) :: input
  type(HarmonicState)              :: this
  
  select type(input); type is(HarmonicState)
    this = input
  type is(SubspaceStatePointer)
    this = HarmonicState(input%state())
  class default
    call err()
  end select
end function

! ----------------------------------------------------------------------
! Type representation.
! ----------------------------------------------------------------------
impure elemental function representation_HarmonicState() result(output)
  implicit none
  
  type(String) :: output
  
  output = 'harmonic'
end function

! ----------------------------------------------------------------------
! Multiply two harmonic states together.
! ----------------------------------------------------------------------
impure elemental function multiply_HarmonicState_HarmonicState(this,that) &
   & result(output)
  implicit none
  
  type(HarmonicState), intent(in) :: this
  type(HarmonicState), intent(in) :: that
  type(HarmonicState)             :: output
  
  if (this%subspace_id/=that%subspace_id) then
    call print_line(CODE_ERROR//': subspaces do not match.')
    call err()
  endif
  
  output = HarmonicState( this%subspace_id,       &
                        & this%frequency,         &
                        & this%state_*that%state_ )
end function

! ----------------------------------------------------------------------
! Generates all harmonic states in a subspace up to a given power.
! ----------------------------------------------------------------------
function generate_harmonic_states(subspace,frequency,modes,maximum_power) &
   & result(output)
  implicit none
  
  type(DegenerateSubspace), intent(in) :: subspace
  real(dp),                 intent(in) :: frequency
  type(ComplexMode),        intent(in) :: modes(:)
  integer,                  intent(in) :: maximum_power
  type(HarmonicState), allocatable     :: output(:)
  
  type(ComplexMode), allocatable :: subspace_modes(:)
  type(ComplexMonomial)          :: harmonic
  type(HarmonicState)            :: state
  
  subspace_modes = subspace%modes(modes)
  harmonic = ComplexMonomial( coefficient = cmplx(1.0_dp,0.0_dp,dp), &
                            & modes       = [ComplexUnivariate::]    )
  state = HarmonicState( subspace_id = subspace%id, &
                       & frequency   = frequency,   &
                       & state       = harmonic     )
  output = generate_harmonic_states_helper(subspace_modes,maximum_power,state)
end function

recursive function generate_harmonic_states_helper(modes,power,state) &
   & result(output)
  implicit none
  
  type(ComplexMode),   intent(in)  :: modes(:)
  integer,             intent(in)  :: power
  type(HarmonicState), intent(in)  :: state
  type(HarmonicState), allocatable :: output(:)
  
  type(HarmonicState) :: output_state
  
  integer :: i
  
  if (size(modes)==0 .or. power==0) then
    output = [state]
  else
    output = [ generate_harmonic_states_helper( modes(2:),   &
             &                                  power,       &
             &                                  state      ) ]
    do i=1,power
      output_state = state
      output_state%state_ = output_state%state_ &
                        & * ComplexUnivariate(mode=modes(1), power=i)
      output = [ output,                                         &
               & generate_harmonic_states_helper( modes(2:),     &
               &                                  power-i,       &
               &                                  output_state ) ]
    enddo
  endif
end function

! ----------------------------------------------------------------------
! Returns the total occupation of a given state.
! ----------------------------------------------------------------------
! The occupation of the state product_{q,i} (a'_{q,i})^(n_{q,i})|0> is equal to
!    sum_{q,i}} n_{q,i}.
impure elemental function total_occupation_HarmonicState(this) result(output)
  implicit none
  
  class(HarmonicState), intent(in) :: this
  integer                          :: output
  
  output = this%state_%total_power()
end function

! ----------------------------------------------------------------------
! Returns the wavevector of a given state.
! ----------------------------------------------------------------------
! The wavevector of the state product_{q,i} (a'_{q,i})^(n_{q,i})|0> is equal to
!    sum_{q,i}} n_{q,i}q.
function wavevector_HarmonicState(this,modes,qpoints) result(output)
  implicit none
  
  class(HarmonicState), intent(in) :: this
  type(ComplexMode),    intent(in) :: modes(:)
  type(QpointData),     intent(in) :: qpoints(:)
  type(QpointData)                 :: output
  
  output = this%state_%wavevector(modes,qpoints)
end function

! ----------------------------------------------------------------------
! SubspaceState methods.
! ----------------------------------------------------------------------
impure elemental function inner_product_HarmonicState(this, &
   & ket,subspace,subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(HarmonicState),     intent(in)           :: this
  class(SubspaceState),     intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  real(dp)                                       :: output
  
  ! <p|q> = 1 if p=q.
  !       = 0 otherwise.
  
  if (present(ket)) then
    if (finite_overlap(this,HarmonicState(ket))) then
      output = 1.0_dp
    else
      output = 0.0_dp
    endif
  else
    output = 1.0_dp
  endif
end function

impure elemental function braket_ComplexMonomial_HarmonicState(this,monomial, &
   & ket,subspace,subspace_basis,anharmonic_data,qpoint) result(output)
  implicit none
  
  class(HarmonicState),     intent(in)           :: this
  type(ComplexMonomial),    intent(in)           :: monomial
  class(SubspaceState),     intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(QpointData),         intent(in), optional :: qpoint
  type(ComplexMonomial)                          :: output
  
  type(HarmonicState) :: harmonic_ket
  type(StateHelper)   :: helper
  
  type(ComplexUnivariate), allocatable :: monomial_modes(:)
  complex(dp)                          :: coefficient
  
  integer :: i
  
  if (present(ket)) then
    ! Process the bra, ket and monomial.
    harmonic_ket = HarmonicState(ket)
    helper = StateHelper( this%state_,         &
                        & harmonic_ket%state_, &
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
  ! The helper function calculates this integral for each single- or double-
  !    mode, up to a factor of 1/(2Nw)^(n/2), where
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

impure elemental function kinetic_energy_HarmonicState(this,ket, &
   & subspace,subspace_basis,anharmonic_data,qpoint) result(output)
  implicit none
  
  class(HarmonicState),     intent(in)           :: this
  class(SubspaceState),     intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(QpointData),         intent(in), optional :: qpoint
  real(dp)                                       :: output
  
  type(HarmonicState) :: harmonic_ket
  type(StateHelper)   :: helper
  
  logical, allocatable :: finite_overlaps(:)
  
  integer :: i
  
  if (present(ket)) then
    ! Process the bra and ket.
    harmonic_ket = HarmonicState(ket)
    helper = StateHelper( this%state_,                       &
                        & harmonic_ket%state_,               &
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
  
  ! Calculate which single- and double-mode states have non-zero overlaps.
  finite_overlaps = finite_overlap_mode(helper%bra,helper%ket)
  
  ! Calculate the kinetic energy, up to a factor of the frequency, w.
  if (all(finite_overlaps)) then
    output = sum(kinetic_energy_mode(helper%bra, helper%ket))
  elseif (count(.not.finite_overlaps)==1) then
    i = first(.not.finite_overlaps)
    output = kinetic_energy_mode(helper%bra(i), helper%ket(i))
  else
    output = 0
    return
  endif
  
  ! Multiply by w, and normalise by the number of primitive cells in the
  !    anharmonic supercell.
  output = output*this%frequency/anharmonic_data%anharmonic_supercell%sc_size
end function

impure elemental function harmonic_potential_energy_HarmonicState(this,ket, &
   & subspace,subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(HarmonicState),     intent(in)           :: this
  class(SubspaceState),     intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  real(dp)                                       :: output
  
  type(HarmonicState) :: harmonic_ket
  type(StateHelper)   :: helper
  
  logical, allocatable :: finite_overlaps(:)
  
  integer :: i
  
  
  if (present(ket)) then
    ! Process the bra and ket.
    harmonic_ket = HarmonicState(ket)
    helper = StateHelper( this%state_,                      &
                        & harmonic_ket%state_,              &
                        & subspace        = subspace,       &
                        & anharmonic_data = anharmonic_data )
  else
    ! Process the bra and ket.
    helper = StateHelper( this%state_,                      &
                        & this%state_,                      &
                        & subspace        = subspace,       &
                        & anharmonic_data = anharmonic_data )
  endif
  
  ! Calculate which single- and double-mode states have non-zero overlaps.
  finite_overlaps = finite_overlap_mode(helper%bra,helper%ket)
  
  ! Calculate the harmonic potential energy,
  !    up to a factor of the frequency, w.
  if (all(finite_overlaps)) then
    output = sum(harmonic_potential_energy_mode(helper%bra, helper%ket))
  elseif (count(.not.finite_overlaps)==1) then
    i = first(.not.finite_overlaps)
    output = harmonic_potential_energy_mode(helper%bra(i), helper%ket(i))
  else
    output = 0
    return
  endif
  
  ! Multiply by w, and normalise by the number of primitive cells in the
  !    anharmonic supercell.
  output = output*this%frequency/anharmonic_data%anharmonic_supercell%sc_size
end function

impure elemental function kinetic_stress_HarmonicState(this,ket, &
   & subspace,subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(HarmonicState),     intent(in)           :: this
  class(SubspaceState),     intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(RealMatrix)                               :: output

  ! TODO
end function

! ----------------------------------------------------------------------
! Helpers for braket, kinetic_energy and potential_energy.
! ----------------------------------------------------------------------
! Returns whether or not braket(bra,ket) is non-zero.
impure elemental function finite_overlap_HarmonicStates(bra,ket) result(output)
  implicit none
  
  type(HarmonicState),      intent(in) :: bra
  type(HarmonicState),      intent(in) :: ket
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
  
  ! <p|q> = 1 if p=q.
  !       = 0 otherwise.
  output = bra%power==ket%power .and. bra%paired_power==ket%paired_power
end function

impure elemental function braket_mode_potential(bra,ket,potential) &
   & result(output)
  implicit none
  
  type(ComplexUnivariate), intent(in) :: bra
  type(ComplexUnivariate), intent(in) :: ket
  type(ComplexUnivariate), intent(in) :: potential
  real(dp)                            :: output
  
  integer :: p_i,p_j,q_i,q_j,n_i,n_j
  
  integer :: m,d,a
  
  integer :: k
  
  real(dp) :: log_term
  
  p_i = bra%power
  p_j = bra%paired_power
  q_i = ket%power
  q_j = ket%paired_power
  n_i = potential%power
  n_j = potential%paired_power
  
  ! N.B. the expressions below for the overlap integrals are factored so that
  !   factorials of p and q are avoided, replaced with binomial coefficients
  !   with small denomianators. This allows for states with very large p and q
  !   to be included.
  ! It is expected that n is much smaller than p and q.
  
  if (bra%id==bra%paired_id) then
    ! I = <p_i|(u_i)^(n_i)|q_i)> = 0 if p+q+n is odd, or |p-q|>n.
    ! Otherwise, defining m=min(p,q), d=|p-q|, and a=(n-d)/2,
    ! I = sum_{k=0}^{min(a,m)}[ J(k) ]
    ! J = n_i!
    !   * sqrt(binomial(m,k)*binomial(m+d,k+d))
    !   / ( sqrt(2Nw)^n_i * 2^(a-k) * (a-k)! * sqrt((d+k)!k!) )
    !
    ! N.B. the factor of 1/sqrt(2Nw)^{n_i} is neglected.
    !
    ! Uses logs to avoid overflowing part way through calculations.
    output = 0
    m = min(p_i,q_i)
    d = abs(p_i-q_i)
    a = (n_i-d)/2
    if (modulo(n_i-d,2)==0) then
      do k=0,min(a,m)
        log_term = log_factorial(n_i)                               &
               & + 0.5_dp*(log_binomial(m,k)+log_binomial(m+d,k+d)) &
               & - (a-k)*log(2.0_dp)                                &
               & - log_factorial(a-k)                               &
               & - 0.5_dp*(log_factorial(k)+log_factorial(d+k))
        output = output + exp(log_term)
      enddo
    endif
  else
    ! I = <p_i,p_j|(u_i)^(n_i)*(u_j)^(n_j)|q_i,q_j>.
    ! I = 0 if p_i-p_j-n_i+n_j-q_i+q_j/=0.
    ! Otherwise, I=sum_{k=max(0,p_i-q_i)}^{min(n_i,p_i,n_j+p_i-q_i)}[J]
    ! J = binomial(p_j+n_i-k)
    !   * sqrt( binomial(n_i,k)
    !         * binomial(n_j,k+q_i-p_i)
    !         * binomial(p_i,k)
    !         * binomial(q_i,k+q_i-p_i)
    !         * n_i!
    !         * n_j!                    )
    !  / sqrt(2Nw)^n_i
    !
    ! N.B. the factor of 1/sqrt(2Nq)^{n_i+n_j} is neglected.
    !
    ! Uses logs to avoid overflowing part way through calculations.
    output = 0
    if (p_i-p_j-n_i+n_j-q_i+q_j==0) then
      do k=max(0,p_i-q_i),min(n_i,p_i,n_j+p_i-q_i)
        log_term = log_binomial(p_j+n_i-k,n_i-k)         &
               & + 0.5_dp*( log_binomial(n_i, k)         &
               &          + log_binomial(n_j, k+q_i-p_i) &
               &          + log_binomial(p_i, k)         &
               &          + log_binomial(q_i, k+q_i-p_i) &
               &          + log_factorial(n_i)           &
               &          + log_factorial(n_j)           )
        output = output + exp(log_term)
      enddo
    endif
  endif
end function

impure elemental function kinetic_energy_mode(bra,ket) result(output)
  implicit none
  
  type(ComplexUnivariate), intent(in) :: bra
  type(ComplexUnivariate), intent(in) :: ket
  real(dp)                            :: output
  
  integer :: p_i,p_j,q_i,q_j
  
  output = 0.0_dp
  
  if (bra%id==bra%paired_id) then
    ! <p|T|q> = -w/4 sqrt(q(q-1)) if q=p+2.
    !         =  w/4 (2q+1)       if q=p.
    !         = -w/4 sqrt(p(p-1)) if q=p-2.
    !         = 0                 otherwise.
    p_i = bra%power
    q_i = ket%power
    if (q_i==p_i+2) then
      output = output - 0.25_dp * sqrt(q_i*(q_i-1.0_dp))
    elseif (q_i==p_i) then
      output = output + 0.25_dp * (2*q_i+1)
    elseif (q_i==p_i-2) then
      output = output - 0.25_dp * sqrt(p_i*(p_i-1.0_dp))
    endif
  else
    ! <p_i,p_j|T|q_i,q_j> = -w/2 sqrt(q_iq_j) if q_i=p_i+1 and q_j=p_j+1.
    !                     =  w/2 (q_i+q_j+1)  if q_i=p_i   and q_j=p_j.
    !                     = -w/2 sqrt(p_ip_j) if q_i=p_i-1 and q_j=p_j-1.
    !                     =  0                otherwise.
    p_i = bra%power
    p_j = bra%paired_power
    q_i = ket%power
    q_j = ket%paired_power
    if (q_i==p_i+1 .and. q_j==p_j+1) then
      output = output - 0.5_dp * sqrt(real(q_i*q_j,dp))
    elseif (q_i==p_i .and. q_j==p_j) then
      output = output + 0.5_dp * (q_i+q_j+1)
    elseif (q_i==p_i-1 .and. q_j==p_j-1) then
      output = output - 0.5_dp * sqrt(real(p_i*p_j,dp))
    endif
  endif
end function

impure elemental function harmonic_potential_energy_mode(bra,ket) &
   & result(output)
  implicit none
  
  type(ComplexUnivariate), intent(in) :: bra
  type(ComplexUnivariate), intent(in) :: ket
  real(dp)                            :: output
  
  integer :: p_i,p_j,q_i,q_j
  
  output = 0.0_dp
  
  if (bra%id==bra%paired_id) then
    ! <p|Vh|q> = w/4 sqrt(q(q-1)) if q=p+2.
    !          = w/4 (2q+1)       if q=p.
    !          = w/4 sqrt(p(p-1)) if q=p-2.
    !          = 0                 otherwise.
    p_i = bra%power
    q_i = ket%power
    if (q_i==p_i+2) then
      output = output + 0.25_dp * sqrt(q_i*(q_i-1.0_dp))
    elseif (q_i==p_i) then
      output = output + 0.25_dp * (2*q_i+1)
    elseif (q_i==p_i-2) then
      output = output + 0.25_dp * sqrt(p_i*(p_i-1.0_dp))
    endif
  else
    ! <p_i,p_j|Vh|q_i,q_j> = w/2 sqrt(q_iq_j) if q_i=p_i+1 and q_j=p_j+1.
    !                      = w/2 (q_i+q_j+1)  if q_i=p_i   and q_j=p_j.
    !                      = w/2 sqrt(p_ip_j) if q_i=p_i-1 and q_j=p_j-1.
    !                      = 0                otherwise.
    p_i = bra%power
    p_j = bra%paired_power
    q_i = ket%power
    q_j = ket%paired_power
    if (q_i==p_i+1 .and. q_j==p_j+1) then
      output = output + 0.5_dp * sqrt(real(q_i*q_j,dp))
    elseif (q_i==p_i .and. q_j==p_j) then
      output = output + 0.5_dp * (q_i+q_j+1)
    elseif (q_i==p_i-1 .and. q_j==p_j-1) then
      output = output + 0.5_dp * sqrt(real(p_i*p_j,dp))
    endif
  endif
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_HarmonicState(this,input)
  implicit none
  
  class(HarmonicState), intent(out) :: this
  type(String),         intent(in)  :: input(:)
  
  integer               :: subspace_id
  real(dp)              :: frequency
  type(ComplexMonomial) :: state
  
  type(String), allocatable :: line(:)
  type(String)              :: state_string
  
  select type(this); type is(HarmonicState)
    line = split_line(input(1))
    subspace_id = int(line(3))
    
    line = split_line(input(2))
    frequency = dble(line(3))
    
    if (input(4)=='|0>') then
      state = ComplexMonomial( coefficient = cmplx(1.0_dp,0.0_dp,dp), &
                             & modes       = [ComplexUnivariate::]    )
    else
      state_string = slice(input(4),2,len(input(4))-1)
      state_string = replace(state_string,'a','u')
      state = ComplexMonomial(1.0_dp//'*'//state_string)
    endif
    
    this = HarmonicState(subspace_id,frequency,state)
  class default
    call err()
  end select
end subroutine

function write_HarmonicState(this) result(output)
  implicit none
  
  class(HarmonicState), intent(in) :: this
  type(String), allocatable        :: output(:)
  
  type(String) :: state_string
  
  select type(this); type is(HarmonicState)
    if (this%state_%total_power()==0) then
      state_string = '|0>'
    else
      state_string = '|'//join(this%state_%modes(),delimiter='*')//'>'
      state_string = replace(state_string,'u','a')
    endif
    output = [ 'Subspace  : '//this%subspace_id, &
             & 'Frequency : '//this%frequency,   &
             & str('State'),                     &
             & state_string                      ]
  class default
    call err()
  end select
end function

function new_HarmonicState_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(HarmonicState)      :: this
  
  call this%read(input)
end function

impure elemental function new_HarmonicState_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(HarmonicState)           :: this
  
  this = HarmonicState(str(input))
end function
end module
