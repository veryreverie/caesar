! ======================================================================
! A state along each complex mode in a degenerate subspace,
!    at a set of q-points for which 2q=G.
! ======================================================================
module monomial_state_real_module
  use common_module
  
  use anharmonic_common_module
  
  use state_helper_module
  use monomial_state_1d_module
  implicit none
  
  private
  
  public :: startup_monomial_state
  
  public :: MonomialStateReal
  
  public :: generate_monomial_states
  
  type, extends(SubspaceState) :: MonomialStateReal
    real(dp)                                    :: frequency
    type(MonomialState1D), private, allocatable :: modes_(:)
  contains
    procedure, public, nopass :: representation => representation_MonomialStateReal
    
    procedure, public :: change_modes => change_modes_MonomialStateReal
    
    procedure, public :: total_power => total_power_MonomialStateReal
    procedure, public :: wavevector => wavevector_MonomialStateReal
    
    procedure, public :: wavefunction => wavefunction_MonomialStateReal
    
    procedure, public :: inner_product => &
                       & inner_product_MonomialStateReal
    procedure, public :: braket_ComplexMonomial => &
                       & braket_ComplexMonomial_MonomialStateReal
    procedure, public :: kinetic_energy => &
                       & kinetic_energy_MonomialStateReal
    procedure, public :: harmonic_potential_energy => &
                       & harmonic_potential_energy_MonomialStateReal
    procedure, public :: kinetic_stress => &
                       & kinetic_stress_MonomialStateReal
    
    ! I/O.
    procedure, public :: read  => read_MonomialStateReal
    procedure, public :: write => write_MonomialStateReal
  end type
  
  interface MonomialStateReal
    module procedure new_MonomialStateReal
    module procedure new_MonomialStateReal_SubspaceState
    module procedure new_MonomialStateReal_Strings
    module procedure new_MonomialStateReal_StringArray
  end interface
  
  interface finite_overlap
    module procedure finite_overlap_MonomialStateReals
  end interface
contains

! Startup procedure.
subroutine startup_monomial_state()
  implicit none
  
  type(MonomialStateReal) :: state
  
  call state%startup()
end subroutine

! ----------------------------------------------------------------------
! Constructor.
! ----------------------------------------------------------------------
function new_MonomialStateReal(subspace_id,frequency,modes) result(this)
  implicit none
  
  integer,               intent(in) :: subspace_id
  real(dp),              intent(in) :: frequency
  type(MonomialState1D), intent(in) :: modes(:)
  type(MonomialStateReal)           :: this
  
  this%subspace_id = subspace_id
  this%frequency   = frequency
  this%modes_      = modes
end function

recursive function new_MonomialStateReal_SubspaceState(input) result(this)
  implicit none
  
  class(SubspaceState), intent(in) :: input
  type(MonomialStateReal)          :: this
  
  select type(input); type is(MonomialStateReal)
    this = input
  type is(SubspaceStatePointer)
    this = MonomialStateReal(input%state())
  class default
    call err()
  end select
end function

! ----------------------------------------------------------------------
! Type representation.
! ----------------------------------------------------------------------
impure elemental function representation_MonomialStateReal() result(output)
  implicit none
  
  type(String) :: output
  
  output = 'monomial'
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
  type(MonomialStateReal), allocatable :: output(:)
  
  integer                 :: ids(:)
  type(MonomialStateReal) :: state
  
  ids = subspace%mode_ids(sort(subspace%mode_ids))
  state = MonomialStateReal( subspace_id = subspace%id,        &
                           & frequency   = frequency,          &
                           & modes       = [MonomialState1D::] )
  output = generate_monomial_states_helper(ids,maximum_power,state)
end function

recursive function generate_monomial_states_helper(ids,power,state) &
   & result(output)
  implicit none
  
  integer,                 intent(in)  :: ids(:)
  integer,                 intent(in)  :: power
  type(MonomialStateReal), intent(in)  :: state
  type(MonomialStateReal), allocatable :: output(:)
  
  type(MonomialStateReal) :: output_state
  
  integer :: i
  
  if (size(ids)==0) then
    output = [state]
  else
    output = [( generate_monomial_states_helper(                           &
              &      ids   = ids(2:),                                      &
              &      power = power-i,                                      &
              &      state = MonomialStateReal(                            &
              &           subspace_id = state%subspace_id,                 &
              &           frequency   = state%frequency,                   &
              &           modes       = [ state%modes_,                    &
              &                           MonomialState1D(ids(1),i) ] ) ), &
              & i=0,                                                       &
              & power                                                      )]
  endif
end function

! ----------------------------------------------------------------------
! Returns the total power of a given state.
! ----------------------------------------------------------------------
! The total power of the state |product_{q,i} (u_{q,i})^(n_{q,i})> is equal to
!    sum_{q,i}} n_{q,i}.
impure elemental function total_power_MonomialStateReal(this) result(output)
  implicit none
  
  class(MonomialStateReal), intent(in) :: this
  integer                              :: output
  
  output = sum(this%modes_%total_power())
end function

! ----------------------------------------------------------------------
! Returns the wavevector of a given state.
! ----------------------------------------------------------------------
! The wavevector of the state |product_{q,i} (u_{q,i})^(n_{q,i})> is equal to
!    sum_{q,i}} n_{q,i}q.
function wavevector_MonomialStateReal(this,modes,qpoints) result(output)
  implicit none
  
  class(MonomialStateReal), intent(in) :: this
  type(ComplexMode),        intent(in) :: modes(:)
  type(QpointData),         intent(in) :: qpoints(:)
  type(FractionVector)                 :: output
  
  output = sum(this%modes_%wavevector(modes,qpoints))
end function

! ----------------------------------------------------------------------
! Returns the wavefunction of the state,
!    with all coefficients accounted for.
! ----------------------------------------------------------------------
impure elemental function wavefunction_MonomialStateReal(this,frequency, &
   & supercell) result(output)
  implicit none
  
  class(MonomialStateReal), intent(in) :: this
  real(dp),                 intent(in) :: frequency
  type(StructureData),      intent(in) :: supercell
  type(String)                         :: output
  
  real(dp) :: log_2nw
  real(dp) :: coefficient
  
  integer :: i
  
  ! |n> = sqrt((2Nw)^n/f(n)) u^n |0>,
  !    where f(n) is the odd factorial of n, f(n) = (2n)!/(n!*2^n).
  
  log_2nw = log(2*supercell%sc_size*frequency)
  coefficient = product(exp(0.5_dp*(            &
     &   this%modes_%power()*log_2nw            &
     & - log_odd_factorial(this%modes_%power()) )))
  
  if (size(this%modes_)==0) then
    output = str(coefficient)//'|0>'
  else
    output = coefficient                                                   // &
       & '*'                                                               // &
       & join(                                                                &
       & [( '(u'//this%modes_(i)%id()//'^'//this%modes_(i)%power()//')',      &
       &    i=1,                                                              &
       &    size(this%modes_)                                            )],  &
       & delimiter='*'                                                   ) // &
       & '|0>'
  endif
end function

! ----------------------------------------------------------------------
! SubspaceState methods.
! ----------------------------------------------------------------------
! Returns whether or not braket(bra,ket) is non-zero.
impure elemental function finite_overlap_MonomialStateReals(bra,ket) &
   & result(output)
  implicit none
  
  type(MonomialStateReal), intent(in) :: bra
  type(MonomialStateReal), intent(in) :: ket
  logical                             :: output
  
  output = all(bra%modes_%finite_overlap(ket%modes_))
end function

impure elemental function inner_product_MonomialStateReal(this, &
   & ket,subspace,subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(MonomialStateReal), intent(in)           :: this
  class(SubspaceState),     intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  real(dp)                                       :: output
  
  type(MonomialStateReal) :: monomial_ket
  
  if (present(ket)) then
    monomial_ket = MonomialStateReal(ket)
    output = product(this%modes_%inner_product(monomial_ket%modes_))
  else
    ! Modes are normalised, so <p|p>=1.
    output = 1
  endif
end function

impure elemental function braket_ComplexMonomial_MonomialStateReal(this,monomial, &
   & ket,subspace,subspace_basis,anharmonic_data,qpoint) result(output)
  implicit none
  
  class(MonomialStateReal),     intent(in)           :: this
  type(ComplexMonomial),    intent(in)           :: monomial
  class(SubspaceState),     intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(QpointData),         intent(in), optional :: qpoint
  type(ComplexMonomial)                          :: output
  
  type(MonomialStateReal) :: monomial_ket
  type(StateHelper)   :: helper
  
  type(ComplexUnivariate), allocatable :: monomial_modes(:)
  complex(dp)                          :: coefficient
  
  integer :: i
  
  monomial_modes = monomial%modes( ids        = this%modes_%id(), &
                                 & paired_ids = this%modes_%id()  )
  
  ! Calculate the coefficient of <bra|X|ket>,
  !    up to the factor of 1/sqrt(2Nw)^n.
  !    - N is the number of primitive cells in the anharmonic supercell.
  !    - w is the frequency of the modes in the subspace.
  !    - n is the power of the modes in the monomial which are integrated.
  if (present(ket)) then
    monomial_ket = MonomialStateReal(ket)
    coefficient = monomial%coefficient                             &
              & * product(this%modes_%braket( monomial_ket%modes_, &
              &                               monomial_modes       ))
  else
    coefficient = monomial%coefficient                       &
              & * product(this%modes_%braket( this%modes_,   &
              &                               monomial_modes ))
  endif
  
  ! Include the factor of (2Nw)^(n/2).
  coefficient = coefficient                                          &
           &  / sqrt( 2.0_dp                                         &
           &        * anharmonic_data%anharmonic_supercell%sc_size   &
           &        * this%frequency                               ) &
           & ** sum(monomial_modes%total_power())
  
  ! Construct the output, from the coefficient and the un-integrated modes.
  output = ComplexMonomial(                                        &
     & coefficient = coefficient,                                  &
     & modes       = monomial%modes(exclude_ids(this%modes_%id())) )
end function

impure elemental function kinetic_energy_MonomialStateReal(this,ket, &
   & subspace,subspace_basis,anharmonic_data,qpoint) result(output)
  implicit none
  
  class(MonomialStateReal), intent(in)           :: this
  class(SubspaceState),     intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(QpointData),         intent(in), optional :: qpoint
  real(dp)                                       :: output
  
  type(MonomialStateReal)            :: monomial_ket
  type(MonomialState1D), allocatable :: bra_modes(:)
  type(MonomialState1D), allocatable :: ket_modes(:)
  real(dp),              allocatable :: overlap(:)
  
  ! |p> = product_i |p_i>
  ! The kinetic energy is given by T = -(1/2N) sum_i d^2/d(u_i^2).
  ! <p_i|d^2/d(u_i)^2|q_i> is calculated up to a factor of 2Nw, so
  !    <p|T|q> = w * sum_i <p_i|%second_derivative(|q_i>) * (<p'|q'>),
  !    where |p'> is |p> excluding |p_i>, so |p>=|p_i>|p'>.
  
  if (present(qpoint)) then
    ! TODO: qpoint
  else
    bra_modes = this%modes_
  endif
  
  if (present(ket)) then
    monomial_ket = MonomialStateReal(ket)
    
    if (present(qpoint)) then
      ! TODO: qpoint
    else
      ket_modes = monomial_ket%modes_
    endif
    
    if (all(bra_modes%finite_overlap(ket_modes))) then
      ! All <p_i|q_i>/=0, so all <p'|q'>=<p|q>/<p_i|q_i>,
      !    so <p|T|q> = -w<p|q>
      !               * sum_i <p_i|%second_derivative(|q_i>) / <p_i|q_i>
      overlap = bra_modes%inner_product(ket_modes)
      output = -this%frequency*product(overlap) &
           & * sum(bra_modes%second_derivative(ket_modes)/overlap)
    else
      ! At least one <p_i|q_i>=0, and for monomial states,
      !    if <p_i|q_i>=0 then <p_i|T|q_i>=0 as well.
      output = 0
    endif
  else
    ! |p>=|q>, so <p'|q'>=1, and so
    !    <p|T|q> = -w*sum_i <p_i|%second_derivative().
    output = -this%frequency*sum(bra_modes%second_derivative())
  endif
end function

impure elemental function harmonic_potential_energy_MonomialStateReal( &
   & this,ket,subspace,subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(MonomialStateReal), intent(in)           :: this
  class(SubspaceState),     intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  real(dp)                                       :: output
  
  type(MonomialStateReal)              :: monomial_ket
  type(MonomialState1D),   allocatable :: bra_modes(:)
  type(MonomialState1D),   allocatable :: ket_modes(:)
  type(ComplexUnivariate), allocatable :: harmonic_potential(:)
  real(dp),                allocatable :: overlap(:)
  
  integer :: i
  
  if (present(qpoint)) then
    ! TODO: qpoint
  else
    bra_modes = this%modes_
  endif
  
  harmonic_potential = ComplexUnivariate( id           = bra_modes%id(), &
                                        & paired_id    = bra_modes%id(), &
                                        & power        = 2,              &
                                        & paired_power = 2               )
  
  if (present(ket)) then
    monomial_ket = MonomialStateReal(ket)
    
    if (present(qpoint)) then
      ! TODO: qpoint
    else
      ket_modes = monomial_ket%modes_
    endif
    
    if (all(bra_modes%finite_overlap(ket_modes))) then
      ! All <p_i|q_i>/=0, so all <p'|q'>=<p|q>/<p_i|q_i>,
      !    so <p|T|q> = -w<p|q>
      !               * sum_i <p_i|%second_derivative(|q_i>) / <p_i|q_i>
      overlap = bra_modes%inner_product(ket_modes)
      output = -this%frequency*product(overlap) &
           & * sum(bra_modes%second_derivative(ket_modes)/overlap)
    else
      ! At least one <p_i|q_i>=0, and for monomial states,
      !    if <p_i|q_i>=0 then <p_i|T|q_i>=0 as well.
      output = 0
    endif
  else
    ! |p>=|q>, so <p'|q'>=1, and so
    !    <p|V|q> = (w^2/4)*sum_i <p_i|%second_derivative().
    output = (this%frequency**2/4) &
         & * sum(bra_modes%braket(potential=harmonic_potential))
  endif
end function

impure elemental function kinetic_stress_MonomialStateReal(this,ket, &
   & subspace,subspace_basis,stress_prefactors,anharmonic_data) result(output)
  implicit none
  
  class(MonomialStateReal),     intent(in)           :: this
  class(SubspaceState),     intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(StressPrefactors),   intent(in)           :: stress_prefactors
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(RealMatrix)                               :: output
  
  ! TODO
end function

! ----------------------------------------------------------------------
! Change the modes of the state by the specified group.
! ----------------------------------------------------------------------
impure elemental function change_modes_MonomialStateReal(this,mode_group) &
   & result(output)
  implicit none
  
  class(MonomialStateReal), intent(in) :: this
  type(Group),          intent(in) :: mode_group
  type(MonomialStateReal)              :: output
  
  integer, allocatable :: ids(:)
  integer, allocatable :: powers(:)
  integer, allocatable :: sort_key(:)
  
  type(ComplexUnivariate), allocatable :: univariates(:)
  type(ComplexMonomial)                :: monomial
  
  integer :: i,ialloc
  
  ! Get the ids and powers of the single-mode terms.
  ids = this%modes_%id()
  powers = this%modes_%power()
  
  ! Change the ids according to the given group.
  ids = mode_group*ids
  
  ! Sort the modes by id.
  sort_key = sort(ids)
  ids = ids(sort_key)
  powers = powers(sort_key)
  
  ! Construct output using the new ids.
  output = MonomialStateReal( subspace_id = this%subspace_id,           &
                            & frequency   = this%frequency,             &
                            & modes       = MonomialState1D(ids,powers) )
end function

! ----------------------------------------------------------------------
! Helpers for braket, kinetic_energy and potential_energy.
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_MonomialStateReal(this,input)
  implicit none
  
  class(MonomialStateReal), intent(out) :: this
  type(String),             intent(in)  :: input(:)
  
  integer                            :: subspace_id
  real(dp)                           :: frequency
  type(MonomialState1D), allocatable :: modes(:)
  
  type(String), allocatable :: line(:)
  
  select type(this); type is(MonomialStateReal)
    line = split_line(input(1))
    subspace_id = int(line(3))
    
    line = split_line(input(2))
    frequency = dble(line(3))
    
    line = split_line(input(4),delimiter='>')
    line = [(line(i)//'>',i=1,size(line))]
    modes = MonomialState1D(line)
    
    this = MonomialStateReal(subspace_id,frequency,modes)
  class default
    call err()
  end select
end subroutine

function write_MonomialStateReal(this) result(output)
  implicit none
  
  class(MonomialStateReal), intent(in) :: this
  type(String), allocatable        :: output(:)
  
  select type(this); type is(MonomialStateReal)
    output = [ 'Subspace  : '//this%subspace_id, &
             & 'Frequency : '//this%frequency,   &
             & str('State'),                     &
             & join(str(this%modes_))            ]
  class default
    call err()
  end select
end function

function new_MonomialStateReal_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(MonomialStateReal)      :: this
  
  call this%read(input)
end function

impure elemental function new_MonomialStateReal_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(MonomialStateReal)           :: this
  
  this = MonomialStateReal(str(input))
end function
end module
