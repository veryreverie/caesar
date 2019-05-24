! ======================================================================
! A state along each complex mode in a degenerate subspace,
!    at a set of q-points for which 2q/=G.
! ======================================================================
module harmonic_state_complex_module
  use common_module
  
  use anharmonic_common_module
  
  use harmonic_state_2d_module
  implicit none
  
  private
  
  public :: startup_harmonic_state_complex
  
  public :: HarmonicStateComplex
  
  public :: prod_complex
  
  type, extends(SubspaceState) :: HarmonicStateComplex
    real(dp)                                    :: frequency
    type(HarmonicState2D), private, allocatable :: modes_(:)
  contains
    procedure, public, nopass :: representation => representation_HarmonicStateComplex
    
    procedure, public :: modes => modes_HarmonicStateComplex
    
    procedure, public :: set_frequency => set_frequency_HarmonicStateComplex
    
    procedure, public :: change_modes => change_modes_HarmonicStateComplex
    
    procedure, public :: total_occupation => total_occupation_HarmonicStateComplex
    procedure, public :: wavevector => wavevector_HarmonicStateComplex
    
    procedure, public :: wavefunction => wavefunction_HarmonicStateComplex
    
    procedure, public :: inner_product => &
                       & inner_product_HarmonicStateComplex
    procedure, public :: braket_ComplexMonomial => &
                       & braket_ComplexMonomial_HarmonicStateComplex
    procedure, public :: kinetic_energy => &
                       & kinetic_energy_HarmonicStateComplex
    procedure, public :: harmonic_potential_energy => &
                       & harmonic_potential_energy_HarmonicStateComplex
    procedure, public :: kinetic_stress => &
                       & kinetic_stress_HarmonicStateComplex
    
    ! I/O.
    procedure, public :: read  => read_HarmonicStateComplex
    procedure, public :: write => write_HarmonicStateComplex
  end type
  
  interface HarmonicStateComplex
    module procedure new_HarmonicStateComplex
    module procedure new_HarmonicStateComplex_SubspaceState
    module procedure new_HarmonicStateComplex_Strings
    module procedure new_HarmonicStateComplex_StringArray
  end interface
  
  interface finite_overlap
    module procedure finite_overlap_HarmonicStateComplexs
  end interface
contains

function prod_complex(lhs,rhs) result(output)
  implicit none
  
  type(HarmonicStateComplex), intent(in) :: lhs
  type(HarmonicStateComplex), intent(in) :: rhs
  type(HarmonicStateComplex)             :: output
  
  output = HarmonicStateComplex(lhs%subspace_id, lhs%frequency, [lhs%modes_,rhs%modes_])
end function

! Startup procedure.
subroutine startup_harmonic_state_complex()
  implicit none
  
  type(HarmonicStateComplex) :: state
  
  call state%startup()
end subroutine

! ----------------------------------------------------------------------
! Constructor.
! ----------------------------------------------------------------------
function new_HarmonicStateComplex(subspace_id,frequency,modes) result(this)
  implicit none
  
  integer,               intent(in) :: subspace_id
  real(dp),              intent(in) :: frequency
  type(HarmonicState2D), intent(in) :: modes(:)
  type(HarmonicStateComplex)        :: this
  
  this%subspace_id = subspace_id
  this%frequency   = frequency
  this%modes_      = modes
end function

recursive function new_HarmonicStateComplex_SubspaceState(input) result(this)
  implicit none
  
  class(SubspaceState), intent(in) :: input
  type(HarmonicStateComplex)       :: this
  
  select type(input); type is(HarmonicStateComplex)
    this = input
  type is(SubspaceStatePointer)
    this = HarmonicStateComplex(input%state())
  class default
    call err()
  end select
end function

! ----------------------------------------------------------------------
! Type representation.
! ----------------------------------------------------------------------
impure elemental function representation_HarmonicStateComplex() result(output)
  implicit none
  
  type(String) :: output
  
  output = 'harmonic complex'
end function

! ----------------------------------------------------------------------
! Returns the modes spanned by the state.
! ----------------------------------------------------------------------
function modes_HarmonicStateComplex(this) result(output)
  implicit none
  
  class(HarmonicStateComplex), intent(in) :: this
  integer, allocatable                   :: output(:)
  
  output = [this%modes_%id(), this%modes_%paired_id()]
end function

! ----------------------------------------------------------------------
! Set the frequency.
! ----------------------------------------------------------------------
impure elemental subroutine set_frequency_HarmonicStateComplex(this,frequency)
  implicit none
  
  class(HarmonicStateComplex), intent(inout) :: this
  real(dp),                    intent(in)    :: frequency
  
  this%frequency = frequency
end subroutine

! ----------------------------------------------------------------------
! Returns the total occupation of a given state.
! ----------------------------------------------------------------------
! The total occupation of the state product_{q,i}|n_{q,i}> is equal to
!    sum_{q,i} n_{q,i}.
impure elemental function total_occupation_HarmonicStateComplex(this) &
   & result(output)
  implicit none
  
  class(HarmonicStateComplex), intent(in) :: this
  integer                                 :: output
  
  output = sum(this%modes_%total_occupation())
end function

! ----------------------------------------------------------------------
! Returns the wavevector of a given state.
! ----------------------------------------------------------------------
! The wavevector of the state product_{q,i} |n_{q,i}> is equal to
!    sum_{q,i} n_{q,i}q.
function wavevector_HarmonicStateComplex(this,modes,qpoints) result(output)
  implicit none
  
  class(HarmonicStateComplex), intent(in) :: this
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
impure elemental function wavefunction_HarmonicStateComplex(this,frequency, &
   & supercell) result(output)
  implicit none
  
  class(HarmonicStateComplex), intent(in) :: this
  real(dp),                 intent(in) :: frequency
  type(StructureData),      intent(in) :: supercell
  type(String)                         :: output
  
  ! TODO
  call err()
end function

! ----------------------------------------------------------------------
! SubspaceState methods.
! ----------------------------------------------------------------------
! Returns whether or not braket(bra,ket) is non-zero.
impure elemental function finite_overlap_HarmonicStateComplexs(bra,ket, &
   & anharmonic_data) result(output)
  implicit none
  
  type(HarmonicStateComplex), intent(in) :: bra
  type(HarmonicStateComplex), intent(in) :: ket
  type(AnharmonicData),       intent(in) :: anharmonic_data
  logical                                :: output
  
  output = all(bra%modes_%finite_overlap(ket%modes_))
end function

impure elemental function inner_product_HarmonicStateComplex(this, &
   & ket,anharmonic_data) result(output)
  implicit none
  
  class(HarmonicStateComplex), intent(in)           :: this
  class(SubspaceState),        intent(in), optional :: ket
  type(AnharmonicData),        intent(in)           :: anharmonic_data
  real(dp)                                          :: output
  
  type(HarmonicStateComplex) :: harmonic_ket
  
  if (present(ket)) then
    harmonic_ket = HarmonicStateComplex(ket)
    output = product(this%modes_%inner_product(harmonic_ket%modes_))
  else
    ! Modes are normalised, so <p|p>=1.
    output = 1
  endif
end function

impure elemental function braket_ComplexMonomial_HarmonicStateComplex(this,monomial, &
   & ket,anharmonic_data) result(output)
  implicit none
  
  class(HarmonicStateComplex), intent(in)           :: this
  type(ComplexMonomial),       intent(in)           :: monomial
  class(SubspaceState),        intent(in), optional :: ket
  type(AnharmonicData),        intent(in)           :: anharmonic_data
  type(ComplexMonomial)                             :: output
  
  type(HarmonicStateComplex) :: harmonic_ket
  
  type(ComplexUnivariate), allocatable :: monomial_modes(:)
  complex(dp)                          :: coefficient
  
  integer :: i
  
  monomial_modes = monomial%modes( ids        = this%modes_%id(),       &
                                 & paired_ids = this%modes_%paired_id() )
  
  ! Calculate the coefficient of <bra|X|ket>,
  !    up to the factor of 1/sqrt(2Nw)^n.
  !    - N is the number of primitive cells in the anharmonic supercell.
  !    - w is the frequency of the modes in the subspace.
  !    - n is the occupation of the modes in the monomial which are integrated.
  if (present(ket)) then
    harmonic_ket = HarmonicStateComplex(ket)
    coefficient = monomial%coefficient                             &
              & * product(this%modes_%braket( harmonic_ket%modes_, &
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
  output = ComplexMonomial(                                       &
     & coefficient = coefficient,                                 &
     & modes       = monomial%modes(exclude_ids=this%modes_%id()) )
end function

impure elemental function kinetic_energy_HarmonicStateComplex(this,ket, &
   & anharmonic_data) result(output)
  implicit none
  
  class(HarmonicStateComplex), intent(in)           :: this
  class(SubspaceState),        intent(in), optional :: ket
  type(AnharmonicData),        intent(in)           :: anharmonic_data
  real(dp)                                          :: output
  
  type(HarmonicStateComplex)         :: harmonic_ket
  type(HarmonicState2D), allocatable :: bra_modes(:)
  type(HarmonicState2D), allocatable :: ket_modes(:)
  
  integer :: i
  
  ! |p> = product_i |p_i>
  ! The kinetic energy is given by T = -(1/N) sum_i d^2/d(u_i^2).
  ! <p_i|d^2/d(u_i)^2|q_i> is calculated up to a factor of 2Nw, so
  !    <p|T|q> = -2w * sum_i <p_i|%second_derivative(|q_i>) * (<p'|q'>),
  !    where |p'> is |p> excluding |p_i>, so |p>=|p_i>|p'>.
  ! Since harmonic states are orthonormal, <p'|q'> is either 0 or 1.
  
  bra_modes = this%modes_
  
  if (present(ket)) then
    harmonic_ket = HarmonicStateComplex(ket)
    ket_modes = harmonic_ket%modes_
    
    if (all(bra_modes%finite_overlap(ket_modes))) then
      ! All <p_i|q_i>/=0, so |q>=|p>,
      !    so <p|T|q> = -2w * sum_i <p_i|%second_derivative(|q_i>).
      output = -2*this%frequency * sum(bra_modes%second_derivative(ket_modes))
    elseif (count(.not.bra_modes%finite_overlap(ket_modes))==1) then
      ! <p_i|q_i>=0, but all other <p_j|q_j>=1.
      ! <p|T|q> = -2w<p_i%second_derivative(|q_i>).
      i = first(.not. bra_modes%finite_overlap(ket_modes))
      output = -2*this%frequency * bra_modes(i)%second_derivative(ket_modes(i))
    else
      ! More than one <p_i|q_i>=0, so <p|V|q>=0.
      output = 0
    endif
  else
    ! |p>=|q>, so <p'|q'>=1, and so
    !    <p|T|q> = -2w*sum_i <p_i|%second_derivative().
    output = -2*this%frequency*sum(bra_modes%second_derivative(bra_modes))
  endif
end function

impure elemental function harmonic_potential_energy_HarmonicStateComplex( &
   & this,ket,anharmonic_data) result(output)
  implicit none
  
  class(HarmonicStateComplex), intent(in)           :: this
  class(SubspaceState),        intent(in), optional :: ket
  type(AnharmonicData),        intent(in)           :: anharmonic_data
  real(dp)                                          :: output
  
  type(HarmonicStateComplex)           :: harmonic_ket
  type(HarmonicState2D),   allocatable :: bra_modes(:)
  type(HarmonicState2D),   allocatable :: ket_modes(:)
  type(ComplexUnivariate), allocatable :: harmonic_potential(:)
  
  integer :: i
  
  ! |p> = product_i |p_i>
  ! The harmonic potential energy is given by V = Nw^2 sum_i (u_i^2).
  ! <p_i|(u_i)^2|q_i> is calculated up to a factor of 2Nw, so
  !    <p|V|q> = (w/2) * sum_i <p_i|%braket(|q_i>,V_i) * (<p'|q'>),
  !    where |p'> is |p> excluding |p_i>, so |p>=|p_i>|p'>.
  ! Since harmonic states are orthonormal, <p'|q'> is either 0 or 1.
  
  bra_modes = this%modes_
  
  harmonic_potential = ComplexUnivariate(    &
     & id           = bra_modes%id(),        &
     & paired_id    = bra_modes%paired_id(), &
     & power        = 1,                     &
     & paired_power = 1                      )
  
  if (present(ket)) then
    harmonic_ket = HarmonicStateComplex(ket)
    ket_modes = harmonic_ket%modes_
    
    if (all(bra_modes%finite_overlap(ket_modes))) then
      ! All <p_i|q_i>/=0, so |q>=|p>,
      !    so <p|T|q> = (w^2/2)
      !               * sum_i <p_i|%braket(|q_i>,V_i)
      output = (this%frequency**2/2) &
           & * sum(bra_modes%braket(ket_modes,harmonic_potential))
    elseif (count(.not.bra_modes%finite_overlap(ket_modes))==1) then
      ! <p_i|q_i>=0, but all other <p_j|q_j>=1.
      ! <p|T|q> = (w^2/2)<p_i%second_derivative(|q_i>).
      i = first(.not. bra_modes%finite_overlap(ket_modes))
      output = (this%frequency**2/2) &
           & * bra_modes(i)%braket(ket_modes(i),harmonic_potential(i))
    else
      ! More than one <p_i|q_i>=0, so <p|V|q>=0.
      output = 0
    endif
  else
    ! |p>=|q>, so <p'|q'>=1, and so
    !    <p|V|q> = (w^2/2)*sum_i <p_i|%braket(|p_i>,V_i).
    output = (this%frequency**2/2) &
         & * sum(bra_modes%braket(bra_modes,harmonic_potential))
  endif
end function

impure elemental function kinetic_stress_HarmonicStateComplex(this,ket, &
   & stress_prefactors,anharmonic_data) result(output)
  implicit none
  
  class(HarmonicStateComplex), intent(in)           :: this
  class(SubspaceState),        intent(in), optional :: ket
  type(StressPrefactors),      intent(in)           :: stress_prefactors
  type(AnharmonicData),        intent(in)           :: anharmonic_data
  type(RealMatrix)                                  :: output
  
  type(HarmonicStateComplex) :: harmonic_ket
  
  logical,  allocatable :: finite_overlap(:)
  
  integer :: i,j,ialloc
  
  ! |p> = product_i |p_i>
  ! The kinetic energy is given by
  !    S = -(1/N) sum_i (I_{i,i}d^2/d(u_i^2) + sum_{j/=i}I_{i,j}d^2/du_idu_j).
  ! <p_i|d^2/d(u_i)^2|q_i> and <p_i|d^2/du_idu_j|q_i> are calculated up to
  !    a factor of 2Nw, so
  !    <p|S|q> = -2w * (
  !      sum_i prefactor_{i,i}*<p_i|%second_derivative(|q_i>) * (<p'|q'>)
  !    + sum_{i,j} prefactor_{i,j}*<p_i|%plus_derivative(|q_i>)
  !                               *<p_j|%minus_derivative(|q_j>)*(<p''|q''>) ),
  ! where |p'> is |p> excluding |p_i>, so |p>=|p_i>|p'>,
  ! and |p''> is |p> excluding |p_i> and |p_j>, so |p>=|p_i>|p_j>|p''>.
  ! Since harmonic states are orthonormal, <p'|q'> and <p''|q''> are
  !    either 0 or 1.
  
  output = dblemat(zeroes(3,3))
  
  if (present(ket)) then
    harmonic_ket = HarmonicStateComplex(ket)
    
    finite_overlap = this%modes_%finite_overlap(harmonic_ket%modes_)
    if (count(.not.finite_overlap)==0) then
      ! All <p_i|q_i> are finite, |q>=|p>.
      ! -> <p|S|p> = -2w sum_i prefactor_{i,i}<p_i|%second_derivative(|p_i>)
      output = -2*this%frequency                                   &
           & * sum( stress_prefactors%prefactor( this%modes_%id(), &
           &                                     this%modes_%id()) &
           &      * this%modes_%second_derivative(this%modes_)     )
    elseif (count(.not.finite_overlap)==2) then
      ! Exactly two <p_i|q_i> are zero. Label these i and j.
      ! <p|S|q> = -2w*
      !   (prefactor_{i,j}+prefactor_{j,i})
      !   * <p_i|%first_derivative(|q_i>)
      !   * <p_j|%first_derivative(|q_j>)
      i = first(.not. finite_overlap)
      j = i + first(.not. finite_overlap(i+1:))
      
      output = ( stress_prefactors%prefactor( this%modes_(i)%id(),         &
           &                                  this%modes_(j)%id()  )       &
           &   * this%modes_(i)%plus_derivative(harmonic_ket%modes_(i))    &
           &   * this%modes_(j)%minus_derivative(harmonic_ket%modes_(j))   &
           &   + stress_prefactors%prefactor( this%modes_(j)%id(),         &
           &                                  this%modes_(i)%id()  )       &
           &   * this%modes_(j)%plus_derivative(harmonic_ket%modes_(j))    &
           &   * this%modes_(i)%minus_derivative(harmonic_ket%modes_(i)) ) &
           & * (-2*this%frequency)
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
end function

! ----------------------------------------------------------------------
! Change the modes of the state by the specified group.
! ----------------------------------------------------------------------
impure elemental function change_modes_HarmonicStateComplex(this,mode_group) &
   & result(output)
  implicit none
  
  class(HarmonicStateComplex), intent(in) :: this
  type(Group),                 intent(in) :: mode_group
  type(HarmonicStateComplex)              :: output
  
  integer, allocatable :: ids(:)
  integer, allocatable :: paired_ids(:)
  integer, allocatable :: occupations(:)
  integer, allocatable :: paired_occupations(:)
  integer, allocatable :: sort_key(:)
  
  type(ComplexUnivariate), allocatable :: univariates(:)
  type(ComplexMonomial)                :: monomial
  
  integer :: i,ialloc
  
  ! Get the ids and occupations of the single-mode terms.
  ids = this%modes_%id()
  paired_ids = this%modes_%paired_id()
  occupations = this%modes_%occupation()
  paired_occupations = this%modes_%paired_occupation()
  
  ! Change the ids according to the given group.
  ids = mode_group*ids
  
  ! Sort the modes by id.
  sort_key = sort(ids)
  ids = ids(sort_key)
  paired_ids = paired_ids(sort_key)
  occupations = occupations(sort_key)
  paired_occupations = paired_occupations(sort_key)
  
  ! Construct output using the new ids.
  output = HarmonicStateComplex(                                     &
     & subspace_id = this%subspace_id,                               &
     & frequency   = this%frequency,                                 &
     & modes       = HarmonicState2D( id           = ids,            &
     &                                paired_id    = paired_ids,     &
     &                                occupation        = occupations,         &
     &                                paired_occupation = paired_occupations ) )
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_HarmonicStateComplex(this,input)
  implicit none
  
  class(HarmonicStateComplex), intent(out) :: this
  type(String),                intent(in)  :: input(:)
  
  integer                            :: subspace_id
  real(dp)                           :: frequency
  type(HarmonicState2D), allocatable :: modes(:)
  
  type(String), allocatable :: line(:)
  
  integer :: i
  
  select type(this); type is(HarmonicStateComplex)
    line = split_line(input(1))
    subspace_id = int(line(3))
    
    line = split_line(input(2))
    frequency = dble(line(3))
    
    line = split_line(input(4),delimiter='>')
    line = [(line(i)//'>',i=1,size(line))]
    modes = HarmonicState2D(line)
    
    this = HarmonicStateComplex(subspace_id,frequency,modes)
  class default
    call err()
  end select
end subroutine

function write_HarmonicStateComplex(this) result(output)
  implicit none
  
  class(HarmonicStateComplex), intent(in) :: this
  type(String), allocatable        :: output(:)
  
  select type(this); type is(HarmonicStateComplex)
    output = [ 'Subspace  : '//this%subspace_id, &
             & 'Frequency : '//this%frequency,   &
             & str('State'),                     &
             & join(str(this%modes_))            ]
  class default
    call err()
  end select
end function

function new_HarmonicStateComplex_Strings(input) result(this)
  implicit none
  
  type(String), intent(in)   :: input(:)
  type(HarmonicStateComplex) :: this
  
  call this%read(input)
end function

impure elemental function new_HarmonicStateComplex_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(HarmonicStateComplex)           :: this
  
  this = HarmonicStateComplex(str(input))
end function
end module
