! ======================================================================
! A state along each complex mode in a degenerate subspace,
!    at a set of q-points for which 2q/=G.
! ======================================================================
module monomial_state_complex_module
  use common_module
  
  use anharmonic_common_module
  
  use monomial_state_2d_module
  implicit none
  
  private
  
  public :: startup_monomial_state_complex
  
  public :: MonomialStateComplex
  
  type, extends(SubspaceState) :: MonomialStateComplex
    real(dp)                                    :: frequency
    type(MonomialState2D), private, allocatable :: modes_(:)
  contains
    procedure, public, nopass :: representation => &
                               & representation_MonomialStateComplex
    
    procedure, public :: mode_ids => mode_ids_MonomialStateComplex
    procedure, public :: paired_mode_ids => &
                       & paired_mode_ids_MonomialStateComplex
    
    procedure, public :: occupation => occupation_MonomialStateComplex
    
    procedure, public :: change_modes => change_modes_MonomialStateComplex
    
    procedure, public :: wavevector => wavevector_MonomialStateComplex
    
    procedure, public :: wavefunction => wavefunction_MonomialStateComplex
    
    procedure, public :: inner_product => &
                       & inner_product_MonomialStateComplex
    procedure, public :: integrate => integrate_MonomialStateComplex
    procedure, public :: kinetic_energy => &
                       & kinetic_energy_MonomialStateComplex
    procedure, public :: harmonic_potential_energy => &
                       & harmonic_potential_energy_MonomialStateComplex
    procedure, public :: kinetic_stress => &
                       & kinetic_stress_MonomialStateComplex
    
    ! I/O.
    procedure, public :: read  => read_MonomialStateComplex
    procedure, public :: write => write_MonomialStateComplex
  end type
  
  interface MonomialStateComplex
    module procedure new_MonomialStateComplex
    module procedure new_MonomialStateComplex_SubspaceState
    module procedure new_MonomialStateComplex_Strings
    module procedure new_MonomialStateComplex_StringArray
  end interface
  
  interface finite_overlap
    module procedure finite_overlap_MonomialStateComplexs
  end interface
contains

! Startup procedure.
subroutine startup_monomial_state_complex()
  implicit none
  
  type(MonomialStateComplex) :: state
  
  call state%startup()
end subroutine

! ----------------------------------------------------------------------
! Constructor.
! ----------------------------------------------------------------------
function new_MonomialStateComplex(subspace_id,frequency,modes) result(this)
  implicit none
  
  integer,               intent(in) :: subspace_id
  real(dp),              intent(in) :: frequency
  type(MonomialState2D), intent(in) :: modes(:)
  type(MonomialStateComplex)        :: this
  
  this%subspace_id = subspace_id
  this%frequency   = frequency
  this%modes_      = modes
end function

recursive function new_MonomialStateComplex_SubspaceState(input) result(this)
  implicit none
  
  class(SubspaceState), intent(in) :: input
  type(MonomialStateComplex)       :: this
  
  select type(input); type is(MonomialStateComplex)
    this = input
  type is(SubspaceStatePointer)
    ! WORKAROUND: ifort doesn't recognise the interface to this function
    !    from within this function, so the full name is used instead.
    !this = MonomialStateComplex(input%state())
    this = new_MonomialStateComplex_SubspaceState(input%state())
  class default
    call err()
  end select
end function

! ----------------------------------------------------------------------
! Type representation.
! ----------------------------------------------------------------------
impure elemental function representation_MonomialStateComplex() result(output)
  implicit none
  
  type(String) :: output
  
  output = 'monomial complex'
end function

! ----------------------------------------------------------------------
! Returns the modes spanned by the state.
! ----------------------------------------------------------------------
function mode_ids_MonomialStateComplex(this) result(output)
  implicit none
  
  class(MonomialStateComplex), intent(in) :: this
  integer, allocatable                    :: output(:)
  
  output = this%modes_%id()
end function

function paired_mode_ids_MonomialStateComplex(this) result(output)
  implicit none
  
  class(MonomialStateComplex), intent(in) :: this
  integer, allocatable                    :: output(:)
  
  output = this%modes_%paired_id()
end function

! ----------------------------------------------------------------------
! Returns the total power of a given state.
! ----------------------------------------------------------------------
! The total power of the state product_{q,i} |(u_{q,i})^(n_{q,i})> is equal to
!    sum_{q,i} n_{q,i}.
impure elemental function occupation_MonomialStateComplex(this) result(output)
  implicit none
  
  class(MonomialStateComplex), intent(in) :: this
  integer                                 :: output
  
  output = sum(this%modes_%total_power())
end function

! ----------------------------------------------------------------------
! Returns the wavevector of a given state.
! ----------------------------------------------------------------------
! The wavevector of the state product_{q,i} |(u_{q,i})^(n_{q,i})> is equal to
!    sum_{q,i} n_{q,i}q.
function wavevector_MonomialStateComplex(this,modes,qpoints) result(output)
  implicit none
  
  class(MonomialStateComplex), intent(in) :: this
  type(ComplexMode),           intent(in) :: modes(:)
  type(QpointData),            intent(in) :: qpoints(:)
  type(FractionVector)                    :: output
  
  integer :: i
  
  output = sum([( this%modes_(i)%wavevector(modes,qpoints), &
                & i=1,                                      &
                & size(this%modes_)                         )])
end function

! ----------------------------------------------------------------------
! Returns the wavefunction of the state,
!    with all coefficients accounted for.
! ----------------------------------------------------------------------
impure elemental function wavefunction_MonomialStateComplex(this,frequency, &
   & supercell) result(output)
  implicit none
  
  class(MonomialStateComplex), intent(in) :: this
  real(dp),                 intent(in) :: frequency
  type(StructureData),      intent(in) :: supercell
  type(String)                         :: output
  
  real(dp) :: log_2nw
  real(dp) :: coefficient
  
  integer :: i
  
  ! |n_i,n_j> = sqrt((2Nw)^(n_i+n_j)/(n_i+n_j)!) (u_i)^(n_i) (u_j)^(n_j) |0>,
  
  log_2nw = log(2*supercell%sc_size*frequency)
  coefficient = product(exp(0.5_dp*(                                       &
     &   (this%modes_%power()+this%modes_%paired_power())*log_2nw          &
     & - log_odd_factorial(this%modes_%power()+this%modes_%paired_power()) )))
  
  if (size(this%modes_)==0) then
    output = str(coefficient)//'|0,0>'
  else
    output = coefficient                                          // &
       & '*'                                                      // &
       & join( [( '(u'//this%modes_(i)%id()          //              &
       &           '^'//this%modes_(i)%power()       //              &
       &           'u'//this%modes_(i)%paired_id()   //              &
       &           '^'//this%modes_(i)%paired_power()//')',          &
       &          i=1,                                               &
       &          size(this%modes_)                         )],      &
       &       delimiter='*'                                    ) // &
       & '|0,0>'
  endif
end function

! ----------------------------------------------------------------------
! SubspaceState methods.
! ----------------------------------------------------------------------
! Returns whether or not braket(bra,ket) is non-zero.
impure elemental function finite_overlap_MonomialStateComplexs(bra,ket, &
   & anharmonic_data) result(output)
  implicit none
  
  type(MonomialStateComplex), intent(in) :: bra
  type(MonomialStateComplex), intent(in) :: ket
  type(AnharmonicData),       intent(in) :: anharmonic_data
  logical                                :: output
  
  output = all(bra%modes_%finite_overlap(ket%modes_))
end function

impure elemental function inner_product_MonomialStateComplex(this, &
   & ket,anharmonic_data) result(output)
  implicit none
  
  class(MonomialStateComplex), intent(in)           :: this
  class(SubspaceState),        intent(in), optional :: ket
  type(AnharmonicData),        intent(in)           :: anharmonic_data
  real(dp)                                          :: output
  
  type(MonomialStateComplex) :: monomial_ket
  
  if (present(ket)) then
    monomial_ket = MonomialStateComplex(ket)
    output = product(this%modes_%inner_product(monomial_ket%modes_))
  else
    ! Modes are normalised, so <p|p>=1.
    output = 1
  endif
end function

impure elemental function integrate_MonomialStateComplex(this, &
   & monomial,ket,anharmonic_data) result(output)
  implicit none
  
  class(MonomialStateComplex), intent(in)           :: this
  type(SparseMonomial),        intent(in)           :: monomial
  class(SubspaceState),        intent(in), optional :: ket
  type(AnharmonicData),        intent(in)           :: anharmonic_data
  complex(dp)                                       :: output
  
  type(MonomialStateComplex) :: monomial_ket
  
  integer :: i
  
  ! Calculate the coefficient of <bra|X|ket>,
  !    up to the factor of 1/sqrt(2Nw)^n.
  !    - N is the number of primitive cells in the anharmonic supercell.
  !    - w is the frequency of the modes in the subspace.
  !    - n is the power of the modes in the monomial which are integrated.
  if (present(ket)) then
    monomial_ket = MonomialStateComplex(ket)
    output = product(this%modes_%braket( monomial_ket%modes_, &
                                       & monomial%modes       ))
  else
    output = product(this%modes_%braket( this%modes_,   &
                                       & monomial%modes ))
  endif
  
  ! Include the factor of (2Nw)^(n/2).
  output = output                                          &
      &  / sqrt( 2.0_dp                                         &
      &        * anharmonic_data%anharmonic_supercell%sc_size   &
      &        * this%frequency                               ) &
      & ** sum(monomial%modes%total_power())
end function

impure elemental function kinetic_energy_MonomialStateComplex(this,ket, &
   & anharmonic_data) result(output)
  implicit none
  
  class(MonomialStateComplex), intent(in)           :: this
  class(SubspaceState),        intent(in), optional :: ket
  type(AnharmonicData),        intent(in)           :: anharmonic_data
  real(dp)                                          :: output
  
  type(MonomialStateComplex)         :: monomial_ket
  type(MonomialState2D), allocatable :: bra_modes(:)
  type(MonomialState2D), allocatable :: ket_modes(:)
  real(dp),              allocatable :: overlap(:)
  
  ! |p> = product_i |p_i>
  ! The kinetic energy is given by T = -(1/N) sum_i d^2/d(u_i^2).
  ! <p_i|d^2/d(u_i)^2|q_i> is calculated up to a factor of 2Nw, so
  !    <p|T|q> = -2w * sum_i <p_i|%second_derivative(|q_i>) * (<p'|q'>),
  !    where |p'> is |p> excluding |p_i>, so |p>=|p_i>|p'>.
  
  bra_modes = this%modes_
  
  if (present(ket)) then
    monomial_ket = MonomialStateComplex(ket)
    ket_modes = monomial_ket%modes_
    
    if (all(bra_modes%finite_overlap(ket_modes))) then
      ! All <p_i|q_i>/=0, so all <p'|q'>=<p|q>/<p_i|q_i>,
      !    so <p|T|q> = -2w<p|q>
      !               * sum_i <p_i|%second_derivative(|q_i>) / <p_i|q_i>
      overlap = bra_modes%inner_product(ket_modes)
      output = -2*this%frequency*product(overlap) &
           & * sum(bra_modes%second_derivative(ket_modes)/overlap)
    else
      ! At least one <p_i|q_i>=0, and for monomial states,
      !    if <p_i|q_i>=0 then <p_i|T|q_i>=0 as well.
      output = 0
    endif
  else
    ! |p>=|q>, so <p'|q'>=1, and so
    !    <p|T|q> = -2w*sum_i <p_i|%second_derivative().
    output = -2*this%frequency*sum(bra_modes%second_derivative(bra_modes))
  endif
end function

impure elemental function harmonic_potential_energy_MonomialStateComplex( &
   & this,ket,anharmonic_data) result(output)
  implicit none
  
  class(MonomialStateComplex), intent(in)           :: this
  class(SubspaceState),        intent(in), optional :: ket
  type(AnharmonicData),        intent(in)           :: anharmonic_data
  real(dp)                                          :: output
  
  type(MonomialStateComplex)           :: monomial_ket
  type(MonomialState2D),   allocatable :: bra_modes(:)
  type(MonomialState2D),   allocatable :: ket_modes(:)
  type(ComplexUnivariate), allocatable :: harmonic_potential(:)
  real(dp),                allocatable :: overlap(:)
  
  integer :: i
  
  ! |p> = product_i |p_i>
  ! The harmonic potential energy is given by V = Nw^2 sum_i (u_i^2).
  ! <p_i|(u_i)^2|q_i> is calculated up to a factor of 2Nw, so
  !    <p|V|q> = (w/2) * sum_i <p_i|%braket(|q_i>,V_i) * (<p'|q'>),
  !    where |p'> is |p> excluding |p_i>, so |p>=|p_i>|p'>.
  
  bra_modes = this%modes_
  
  harmonic_potential = ComplexUnivariate(    &
     & id           = bra_modes%id(),        &
     & paired_id    = bra_modes%paired_id(), &
     & power        = 1,                     &
     & paired_power = 1                      )
  
  if (present(ket)) then
    monomial_ket = MonomialStateComplex(ket)
    ket_modes = monomial_ket%modes_
    
    if (all(bra_modes%finite_overlap(ket_modes))) then
      ! All <p_i|q_i>/=0, so all <p'|q'>=<p|q>/<p_i|q_i>,
      !    so <p|T|q> = (w^2/2)<p|q>
      !               * sum_i <p_i|%braket(|q_i>,V_i) / <p_i|q_i>
      overlap = bra_modes%inner_product(ket_modes)
      output = (this%frequency**2/2)*product(overlap) &
           & * sum(bra_modes%braket(ket_modes,harmonic_potential)/overlap)
    else
      ! At least one <p_i|q_i>=0, and for monomial states,
      !    if <p_i|q_i>=0 then <p_i|V_i|q_i>=0 as well.
      output = 0
    endif
  else
    ! |p>=|q>, so <p'|q'>=1, and so
    !    <p|V|q> = (w^2/2)*sum_i <p_i|%braket(|p_i>,V_i).
    output = (this%frequency**2/2) &
         & * sum(bra_modes%braket(bra_modes,harmonic_potential))
  endif
end function

impure elemental function kinetic_stress_MonomialStateComplex(this,ket, &
   & stress_prefactors,anharmonic_data) result(output)
  implicit none
  
  class(MonomialStateComplex), intent(in)           :: this
  class(SubspaceState),        intent(in), optional :: ket
  type(StressPrefactors),      intent(in)           :: stress_prefactors
  type(AnharmonicData),        intent(in)           :: anharmonic_data
  type(RealMatrix)                                  :: output
  
  type(MonomialStateComplex) :: monomial_ket
  
  logical,  allocatable :: finite_overlap(:)
  real(dp), allocatable :: overlaps(:)
  
  integer :: i,j,ialloc
  
  ! |p> = product_i |p_i>
  ! The kinetic stress is given by
  !    S = -(1/NV) sum_i (I_{i,i}d^2/d(u_i^2) + sum_{j/=i}I_{i,j}d^2/du_idu_j).
  ! <p_i|d^2/d(u_i)^2|q_i> and <p_i|d^2/du_idu_j|q_i> are calculated up to
  !    a factor of 2Nw, so
  !    <p|S|q> = -2w * (
  !      sum_i prefactor_{i,i}*<p_i|%second_derivative(|q_i>) * (<p'|q'>)
  !    + sum_{i,j} prefactor_{i,j}*<p_i|%plus_derivative(|q_i>)
  !                               *<p_j|%minus_derivative(|q_j>)*(<p''|q''>) ),
  ! where |p'> is |p> excluding |p_i>, so |p>=|p_i>|p'>,
  ! and |p''> is |p> excluding |p_i> and |p_j>, so |p>=|p_i>|p_j>|p''>.
  
  output = dblemat(zeroes(3,3))
  
  if (present(ket)) then
    monomial_ket = MonomialStateComplex(ket)
    
    finite_overlap = this%modes_%finite_overlap(monomial_ket%modes_)
    if (count(.not.finite_overlap)==0) then
      ! All <p_i|q_i> are finite, so
      !    <p'|q'> = <p|q>/<p_i|q_i>
      !    <p''|q''> = <p|q>/(<p_i|q_i><p_j|q_j>)
      ! For monomial states, if <p_i|q_i>/=0 then
      !    <p_i|d/d(u_i)|q_i>=0 by parity.
      ! S = -2w<p|q>
      !   *  sum_i prefactor_{i,i}*<p_i|%second_derivative(|q_i>)/<p_i|q_i>
      
      overlaps = this%modes_%inner_product(monomial_ket%modes_)
      output = sum( stress_prefactors%prefactor( this%modes_%id(),       &
           &                                     this%modes_%id()  )     &
           &      * this%modes_%second_derivative(monomial_ket%modes_)   &
           &      / overlaps                                           ) &
           & * (-2*this%frequency)*product(overlaps)
    elseif (count(.not.finite_overlap)==2) then
      ! Exactly two <p_i|q_i> are zero. Label these i and j.
      ! <p|S|q> = -2w*
      !   (prefactor_{i,j}+prefactor_{j,i})
      !   * <p_i|%first_derivative(|q_i>)
      !   * <p_j|%first_derivative(|q_j>)
      !   * <p''|q''>
      i = first(.not. finite_overlap)
      j = i + first(.not. finite_overlap(i+1:))
      
      ! Calculate <p|S|q> up to the factor of <p''|q''>.
      output = ( stress_prefactors%prefactor( this%modes_(i)%id(),         &
           &                                  this%modes_(j)%id()  )       &
           &   * this%modes_(i)%plus_derivative(monomial_ket%modes_(i))    &
           &   * this%modes_(j)%minus_derivative(monomial_ket%modes_(j))   &
           &   + stress_prefactors%prefactor( this%modes_(j)%id(),         &
           &                                  this%modes_(i)%id()  )       &
           &   * this%modes_(j)%plus_derivative(monomial_ket%modes_(j))    &
           &   * this%modes_(i)%minus_derivative(monomial_ket%modes_(i)) ) &
           & * (-2*this%frequency)
      
      ! Calculate and include the factor of <p''|q''>.
      overlaps = this%modes_%inner_product(monomial_ket%modes_)
      output = output                     &
           & * product(overlaps(:i-1))    &
           & * product(overlaps(i+1:j-1)) &
           & * product(overlaps(j+1:))
    else
      ! More than two <p_i|q_i>=0, so the whole expression is zero.
      return
    endif
  else
    ! |p>=|q>, so all first derivative expectations are zero.
    ! Also, <p|p>=1, so <p'|p'>=1.
    ! -> <p|S|p> = -2w sum_i prefactor_{i,i}<p_i|%second_derivative(|p_i>)
    output = -2*this%frequency                                   &
         & * sum( stress_prefactors%prefactor( this%modes_%id(), &
         &                                     this%modes_%id()) &
         &      * this%modes_%second_derivative(this%modes_)     )
  endif
  
  ! Divide by the volume.
  output = output / anharmonic_data%structure%volume
end function

! ----------------------------------------------------------------------
! Change the modes of the state by the specified group.
! ----------------------------------------------------------------------
impure elemental function change_modes_MonomialStateComplex(this,mode_group) &
   & result(output)
  implicit none
  
  class(MonomialStateComplex), intent(in) :: this
  type(Group),                 intent(in) :: mode_group
  type(MonomialStateComplex)              :: output
  
  integer, allocatable :: ids(:)
  integer, allocatable :: paired_ids(:)
  integer, allocatable :: powers(:)
  integer, allocatable :: paired_powers(:)
  integer, allocatable :: sort_key(:)
  
  type(ComplexUnivariate), allocatable :: univariates(:)
  type(ComplexMonomial)                :: monomial
  
  integer :: i,ialloc
  
  ! Get the ids and powers of the single-mode terms.
  ids = this%modes_%id()
  paired_ids = this%modes_%paired_id()
  powers = this%modes_%power()
  paired_powers = this%modes_%paired_power()
  
  ! Change the ids according to the given group.
  ids = mode_group*ids
  
  ! Sort the modes by id.
  sort_key = sort(ids)
  ids = ids(sort_key)
  paired_ids = paired_ids(sort_key)
  powers = powers(sort_key)
  paired_powers = paired_powers(sort_key)
  
  ! Construct output using the new ids.
  output = MonomialStateComplex(                                     &
     & subspace_id = this%subspace_id,                               &
     & frequency   = this%frequency,                                 &
     & modes       = MonomialState2D( id           = ids,            &
     &                                paired_id    = paired_ids,     &
     &                                power        = powers,         &
     &                                paired_power = paired_powers ) )
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_MonomialStateComplex(this,input)
  implicit none
  
  class(MonomialStateComplex), intent(out) :: this
  type(String),                intent(in)  :: input(:)
  
  integer                            :: subspace_id
  real(dp)                           :: frequency
  type(MonomialState2D), allocatable :: modes(:)
  
  type(String), allocatable :: line(:)
  
  integer :: i
  
  select type(this); type is(MonomialStateComplex)
    line = split_line(input(1))
    subspace_id = int(line(3))
    
    line = split_line(input(2))
    frequency = dble(line(3))
    
    line = split_line(input(4),delimiter='>')
    line = [(line(i)//'>',i=1,size(line))]
    modes = MonomialState2D(line)
    
    this = MonomialStateComplex(subspace_id,frequency,modes)
  class default
    call err()
  end select
end subroutine

function write_MonomialStateComplex(this) result(output)
  implicit none
  
  class(MonomialStateComplex), intent(in) :: this
  type(String), allocatable        :: output(:)
  
  select type(this); type is(MonomialStateComplex)
    output = [ 'Subspace  : '//this%subspace_id,    &
             & 'Frequency : '//this%frequency,      &
             & str('State'),                        &
             & join(str(this%modes_), delimiter='') ]
  class default
    call err()
  end select
end function

function new_MonomialStateComplex_Strings(input) result(this)
  implicit none
  
  type(String), intent(in)   :: input(:)
  type(MonomialStateComplex) :: this
  
  call this%read(input)
end function

impure elemental function new_MonomialStateComplex_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(MonomialStateComplex)           :: this
  
  this = MonomialStateComplex(str(input))
end function
end module
