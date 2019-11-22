! ======================================================================
! A state along each complex mode in a degenerate subspace,
!    at a set of q-points for which 2q=G.
! ======================================================================
module harmonic_state_real_module
  use common_module
  
  use anharmonic_common_module
  
  use harmonic_state_1d_module
  implicit none
  
  private
  
  public :: startup_harmonic_state_real
  
  public :: HarmonicStateReal
  
  public :: prod_real
  
  type, extends(SubspaceState) :: HarmonicStateReal
    real(dp)                                    :: frequency
    type(HarmonicState1D), private, allocatable :: modes_(:)
  contains
    procedure, public, nopass :: representation => &
                               & representation_HarmonicStateReal
    
    procedure, public :: mode_ids => mode_ids_HarmonicStateReal
    procedure, public :: paired_mode_ids => paired_mode_ids_HarmonicStateReal
    
    procedure, public :: change_modes => change_modes_HarmonicStateReal
    
    procedure, public :: occupation => occupation_HarmonicStateReal
    
    procedure, public :: wavevector => wavevector_HarmonicStateReal
    
    procedure, public :: wavefunction => wavefunction_HarmonicStateReal
    
    procedure, public :: inner_product => &
                       & inner_product_HarmonicStateReal
    procedure, public :: integrate => integrate_HarmonicStateReal
    procedure, public :: kinetic_energy => &
                       & kinetic_energy_HarmonicStateReal
    procedure, public :: harmonic_potential_energy => &
                       & harmonic_potential_energy_HarmonicStateReal
    procedure, public :: kinetic_stress => &
                       & kinetic_stress_HarmonicStateReal
    
    ! I/O.
    procedure, public :: read  => read_HarmonicStateReal
    procedure, public :: write => write_HarmonicStateReal
  end type
  
  interface HarmonicStateReal
    module procedure new_HarmonicStateReal
    module procedure new_HarmonicStateReal_SubspaceState
    module procedure new_HarmonicStateReal_Strings
    module procedure new_HarmonicStateReal_StringArray
  end interface
  
  interface finite_overlap
    module procedure finite_overlap_HarmonicStateReals
  end interface
contains

function prod_real(lhs,rhs) result(output)
  implicit none
  
  type(HarmonicStateReal), intent(in) :: lhs
  type(HarmonicStateReal), intent(in) :: rhs
  type(HarmonicStateReal)             :: output
  
  output = HarmonicStateReal( lhs%subspace_id,        &
                            & lhs%frequency,          &
                            & [lhs%modes_,rhs%modes_] )
end function

! Startup procedure.
subroutine startup_harmonic_state_real()
  implicit none
  
  type(HarmonicStateReal) :: state
  
  call state%startup()
end subroutine

! ----------------------------------------------------------------------
! Constructor.
! ----------------------------------------------------------------------
function new_HarmonicStateReal(subspace_id,frequency,modes) result(this)
  implicit none
  
  integer,               intent(in) :: subspace_id
  real(dp),              intent(in) :: frequency
  type(HarmonicState1D), intent(in) :: modes(:)
  type(HarmonicStateReal)           :: this
  
  this%subspace_id = subspace_id
  this%frequency   = frequency
  this%modes_      = modes
end function

recursive function new_HarmonicStateReal_SubspaceState(input) result(this)
  implicit none
  
  class(SubspaceState), intent(in) :: input
  type(HarmonicStateReal)          :: this
  
  select type(input); type is(HarmonicStateReal)
    this = input
  type is(SubspaceStatePointer)
    ! WORKAROUND: ifort doesn't recognise the interface to this function
    !    from within this function, so the full name is used instead.
    !this = HarmonicStateReal(input%state())
    this = new_HarmonicStateReal_SubspaceState(input%state())
  class default
    call err()
  end select
end function

! ----------------------------------------------------------------------
! Type representation.
! ----------------------------------------------------------------------
impure elemental function representation_HarmonicStateReal() result(output)
  implicit none
  
  type(String) :: output
  
  output = 'harmonic real'
end function

! ----------------------------------------------------------------------
! Returns the modes spanned by the state.
! ----------------------------------------------------------------------
function mode_ids_HarmonicStateReal(this) result(output)
  implicit none
  
  class(HarmonicStateReal), intent(in) :: this
  integer, allocatable                 :: output(:)
  
  output = this%modes_%id()
end function

function paired_mode_ids_HarmonicStateReal(this) result(output)
  implicit none
  
  class(HarmonicStateReal), intent(in) :: this
  integer, allocatable                 :: output(:)
  
  output = this%modes_%id()
end function

! ----------------------------------------------------------------------
! Returns the total occupation of a given state.
! ----------------------------------------------------------------------
! The total occupation of the state product_{q,i} |n_{q,i}> is equal to
!    sum_{q,i} n_{q,i}.
impure elemental function occupation_HarmonicStateReal(this) result(output)
  implicit none
  
  class(HarmonicStateReal), intent(in) :: this
  integer                              :: output
  
  output = sum(this%modes_%total_occupation())
end function

! ----------------------------------------------------------------------
! Returns the wavevector of a given state.
! ----------------------------------------------------------------------
! The wavevector of the state product_{q,i} |n_{q,i}> is equal to
!    sum_{q,i} n_{q,i}q.
function wavevector_HarmonicStateReal(this,modes,qpoints) result(output)
  implicit none
  
  class(HarmonicStateReal), intent(in) :: this
  type(ComplexMode),        intent(in) :: modes(:)
  type(QpointData),         intent(in) :: qpoints(:)
  type(FractionVector)                 :: output
  
  integer :: i
  
  output = sum([( this%modes_(i)%wavevector(modes,qpoints), &
                & i=1,                                      &
                & size(this%modes_)                         )])
end function

! ----------------------------------------------------------------------
! Returns the wavefunction of the state,
!    with all coefficients accounted for.
! ----------------------------------------------------------------------
impure elemental function wavefunction_HarmonicStateReal(this,frequency, &
   & supercell) result(output)
  implicit none
  
  class(HarmonicStateReal), intent(in) :: this
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
impure elemental function finite_overlap_HarmonicStateReals(bra,ket, &
   & anharmonic_data) result(output)
  implicit none
  
  type(HarmonicStateReal), intent(in) :: bra
  type(HarmonicStateReal), intent(in) :: ket
  type(AnharmonicData),    intent(in) :: anharmonic_data
  logical                             :: output
  
  output = all(bra%modes_%finite_overlap(ket%modes_))
end function

impure elemental function inner_product_HarmonicStateReal(this, &
   & ket,anharmonic_data) result(output)
  implicit none
  
  class(HarmonicStateReal), intent(in)           :: this
  class(SubspaceState),     intent(in), optional :: ket
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  real(dp)                                       :: output
  
  type(HarmonicStateReal) :: harmonic_ket
  
  if (present(ket)) then
    harmonic_ket = HarmonicStateReal(ket)
    output = product(this%modes_%inner_product(harmonic_ket%modes_))
  else
    ! Modes are normalised, so <p|p>=1.
    output = 1
  endif
end function

impure elemental function integrate_HarmonicStateReal(this, &
   & monomial,ket,anharmonic_data) result(output)
  implicit none
  
  class(HarmonicStateReal), intent(in)           :: this
  type(SparseMonomial),     intent(in)           :: monomial
  class(SubspaceState),     intent(in), optional :: ket
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  complex(dp)                                    :: output
  
  integer :: i
  
  ! Calculate the coefficient of <bra|X|ket>,
  !    up to the factor of 1/sqrt(2Nw)^n.
  !    - N is the number of primitive cells in the anharmonic supercell.
  !    - w is the frequency of the modes in the subspace.
  !    - n is the occupation of the modes in the monomial which are integrated.
  
  ! N.B. this function is called many times, and so uses
  !    SubspaceStatePointer%state_ directly rather than calling
  !    HarmonicStateReal(ket), in order to improve runtimes.
  if (present(ket)) then
    select type(ket); type is(SubspaceStatePointer)
      associate(ket2=>ket%state_)
        select type(ket2); type is(HarmonicStateReal)
          output = product(this%modes_%braket( ket2%modes_,   &
                                             & monomial%modes ))
        class default
          call err()
        end select
      end associate
    type is(HarmonicStateReal)
      output = product(this%modes_%braket( ket%modes_,    &
                                         & monomial%modes ))
    class default
      call err()
    end select
  else
    output = product(this%modes_%braket( this%modes_,   &
                                       & monomial%modes ))
  endif
  
  ! Include the factor of (2Nw)^(n/2).
  output = output                                               &
      &  / sqrt( 2.0_dp                                         &
      &        * anharmonic_data%anharmonic_supercell%sc_size   &
      &        * this%frequency                               ) &
      & ** sum(monomial%modes%total_power())
end function

impure elemental function kinetic_energy_HarmonicStateReal(this,ket, &
   & anharmonic_data) result(output)
  implicit none
  
  class(HarmonicStateReal), intent(in)           :: this
  class(SubspaceState),     intent(in), optional :: ket
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  real(dp)                                       :: output
  
  type(HarmonicStateReal)            :: harmonic_ket
  type(HarmonicState1D), allocatable :: bra_modes(:)
  type(HarmonicState1D), allocatable :: ket_modes(:)
  
  integer :: i
  
  ! |p> = product_i |p_i>
  ! The kinetic energy is given by T = -(1/2N) sum_i d^2/d(u_i^2).
  ! <p_i|d^2/d(u_i)^2|q_i> is calculated up to a factor of 2Nw, so
  !    <p|T|q> = -w * sum_i <p_i|%second_derivative(|q_i>) * (<p'|q'>),
  !    where |p'> is |p> excluding |p_i>, so |p>=|p_i>|p'>.
  ! Since harmonic states are orthonormal, <p'|q'> is either 0 or 1.
  
  bra_modes = this%modes_
  
  if (present(ket)) then
    harmonic_ket = HarmonicStateReal(ket)
    ket_modes = harmonic_ket%modes_
    
    if (all(bra_modes%finite_overlap(ket_modes))) then
      ! All <p_i|q_i>/=0, so |q>=|p>,
      !    so <p|T|q> = -w * sum_i <p_i|%second_derivative(|q_i>).
      output = -this%frequency * sum(bra_modes%second_derivative(ket_modes))
    elseif (count(.not.bra_modes%finite_overlap(ket_modes))==1) then
      ! <p_i|q_i>=0, but all other <p_j|q_j>=1.
      ! <p|T|q> = -w<p_i%second_derivative(|q_i>).
      i = first(.not. bra_modes%finite_overlap(ket_modes))
      output = -this%frequency * bra_modes(i)%second_derivative(ket_modes(i))
    else
      ! More than one <p_i|q_i>=0, so <p|T|q>=0.
      output = 0
    endif
  else
    ! |p>=|q>, so <p'|q'>=1, and so
    !    <p|T|q> = -w*sum_i <p_i|%second_derivative().
    output = -this%frequency*sum(bra_modes%second_derivative(bra_modes))
  endif
end function

impure elemental function harmonic_potential_energy_HarmonicStateReal( &
   & this,ket,anharmonic_data) result(output)
  implicit none
  
  class(HarmonicStateReal), intent(in)           :: this
  class(SubspaceState),     intent(in), optional :: ket
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  real(dp)                                       :: output
  
  type(HarmonicStateReal)              :: harmonic_ket
  type(HarmonicState1D),   allocatable :: bra_modes(:)
  type(HarmonicState1D),   allocatable :: ket_modes(:)
  type(ComplexUnivariate), allocatable :: harmonic_potential(:)
  
  integer :: i
  
  ! |p> = product_i |p_i>
  ! The harmonic potential energy is given by V = (Nw^2/2) sum_i (u_i^2).
  ! <p_i|(u_i)^2|q_i> is calculated up to a factor of 2Nw, so
  !    <p|V|q> = (w/4) * sum_i <p_i|%braket(|q_i>,V_i) * (<p'|q'>),
  !    where |p'> is |p> excluding |p_i>, so |p>=|p_i>|p'>.
  ! Since harmonic states are orthonormal, <p'|q'> is either 0 or 1.
  
  bra_modes = this%modes_
  
  harmonic_potential = ComplexUnivariate( id           = bra_modes%id(), &
                                        & paired_id    = bra_modes%id(), &
                                        & power        = 2,              &
                                        & paired_power = 2               )
  
  if (present(ket)) then
    harmonic_ket = HarmonicStateReal(ket)
    ket_modes = harmonic_ket%modes_
    
    if (all(bra_modes%finite_overlap(ket_modes))) then
      ! All <p_i|q_i>/=0, so all <p'|q'>=<p|q>/<p_i|q_i>,
      !    so <p|T|q> = (w^2/4)
      !               * sum_i <p_i|%braket(|q_i>,V_i)
      output = (this%frequency**2/4) &
           & * sum(bra_modes%braket(ket_modes,harmonic_potential))
    elseif (count(.not.bra_modes%finite_overlap(ket_modes))==1) then
      ! <p_i|q_i>=0, but all other <p_j|q_j>=1.
      ! <p|V|q> = (w^2/4)<p_i%second_derivative(|q_i>).
      i = first(.not. bra_modes%finite_overlap(ket_modes))
      output = (this%frequency**2/4) &
           & * bra_modes(i)%braket(ket_modes(i),harmonic_potential(i))
    else
      ! More than one <p_i|q_i>=0, so <p|V|q>=0.
      output = 0
    endif
  else
    ! |p>=|q>, so <p'|q'>=1, and so
    !    <p|V|q> = (w^2/4)*sum_i <p_i|%braket(|p_i>,V_i).
    output = (this%frequency**2/4) &
         & * sum(bra_modes%braket(bra_modes,harmonic_potential))
  endif
end function

impure elemental function kinetic_stress_HarmonicStateReal(this,ket, &
   & stress_prefactors,anharmonic_data) result(output)
  implicit none
  
  class(HarmonicStateReal), intent(in)           :: this
  class(SubspaceState),     intent(in), optional :: ket
  type(StressPrefactors),   intent(in)           :: stress_prefactors
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(RealMatrix)                               :: output
  
  type(HarmonicStateReal) :: harmonic_ket
  
  logical,  allocatable :: finite_overlap(:)
  real(dp), allocatable :: first_derivatives(:)
  
  integer :: i,j,ialloc
  
  ! |p> = product_i |p_i>
  ! The kinetic stress is given by
  !    S = -(1/2NV) sum_i (I_{i,i}d^2/d(u_i^2) + sum_{j/=i}I_{i,j}d^2/du_idu_j).
  ! <p_i|d^2/d(u_i)^2|q_i> and <p_i|d^2/du_idu_j|q_i> are calculated up to
  !    a factor of 2Nw, so
  !    <p|S|q> = -w * (
  !      sum_i prefactor_{i,i}*<p_i|%second_derivative(|q_i>) * (<p'|q'>)
  !    + sum_{i,j} prefactor_{i,j}*<p_i|%first_derivative(|q_i>)
  !                               *<p_j|%first_derivative(|q_j>)*(<p''|q''>) ),
  ! where |p'> is |p> excluding |p_i>, so |p>=|p_i>|p'>,
  ! and |p''> is |p> excluding |p_i> and |p_j>, so |p>=|p_i>|p_j>|p''>.
  ! Since harmonic states are orthonormal, <p'|q'> and <p''|q''> are
  !    either 0 or 1.
  
  output = dblemat(zeroes(3,3))
  
  if (present(ket)) then
    harmonic_ket = HarmonicStateReal(ket)
    
    finite_overlap = this%modes_%finite_overlap(harmonic_ket%modes_)
    if (count(.not.finite_overlap)==0) then
      ! All <p_i|q_i> are finite, so |q>=|p>.
      ! -> <p|S|p> = -w sum_i prefactor_{i,i}<p_i|%second_derivative(|p_i>)
      output = -this%frequency                                     &
           & * sum( stress_prefactors%prefactor( this%modes_%id(), &
           &                                     this%modes_%id()) &
           &      * this%modes_%second_derivative(this%modes_)     )
    elseif (count(.not.finite_overlap)==2) then
      ! Exactly two <p_i|q_i> are zero. Label these i and j.
      ! <p|S|q> = -w*
      !   (prefactor_{i,j}+prefactor_{j,i})
      !   * <p_i|%first_derivative(|q_i>)
      !   * <p_j|%first_derivative(|q_j>)
      !   * <p''|q''>
      i = first(.not. finite_overlap)
      j = i + first(.not. finite_overlap(i+1:))
      
      output = ( stress_prefactors%prefactor( this%modes_(i)%id(),     &
           &                                  this%modes_(j)%id()  )   &
           &   + stress_prefactors%prefactor( this%modes_(j)%id(),     &
           &                                  this%modes_(i)%id()  ) ) &
           & * this%modes_(i)%first_derivative(harmonic_ket%modes_(i)) &
           & * this%modes_(j)%first_derivative(harmonic_ket%modes_(j)) &
           & * (-this%frequency)
    else
      ! More than two <p_i|q_i>=0, so the whole expression is zero.
      return
    endif
  else
    ! |p>=|q>, so all first derivative expectations are zero.
    ! -> <p|S|p> = -w sum_i prefactor_{i,i}<p_i|%second_derivative(|p_i>)
    output = -this%frequency                                     &
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
impure elemental function change_modes_HarmonicStateReal(this,mode_group) &
   & result(output)
  implicit none
  
  class(HarmonicStateReal), intent(in) :: this
  type(Group),              intent(in) :: mode_group
  type(HarmonicStateReal)              :: output
  
  integer, allocatable :: ids(:)
  integer, allocatable :: occupations(:)
  integer, allocatable :: sort_key(:)
  
  type(ComplexUnivariate), allocatable :: univariates(:)
  type(ComplexMonomial)                :: monomial
  
  integer :: i,ialloc
  
  ! Get the ids and occupations of the single-mode terms.
  ids = this%modes_%id()
  occupations = this%modes_%occupation()
  
  ! Change the ids according to the given group.
  ids = mode_group*ids
  
  ! Sort the modes by id.
  sort_key = sort(ids)
  ids = ids(sort_key)
  occupations = occupations(sort_key)
  
  ! Construct output using the new ids.
  output = HarmonicStateReal( subspace_id = this%subspace_id,                &
                            & frequency   = this%frequency,                  &
                            & modes       = HarmonicState1D(ids,occupations) )
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_HarmonicStateReal(this,input)
  implicit none
  
  class(HarmonicStateReal), intent(out) :: this
  type(String),             intent(in)  :: input(:)
  
  integer                            :: subspace_id
  real(dp)                           :: frequency
  type(HarmonicState1D), allocatable :: modes(:)
  
  type(String), allocatable :: line(:)
  
  integer :: i
  
  select type(this); type is(HarmonicStateReal)
    line = split_line(input(1))
    subspace_id = int(line(3))
    
    line = split_line(input(2))
    frequency = dble(line(3))
    
    line = split_line(input(4),delimiter='>')
    line = [(line(i)//'>',i=1,size(line))]
    modes = HarmonicState1D(line)
    
    this = HarmonicStateReal(subspace_id,frequency,modes)
  class default
    call err()
  end select
end subroutine

function write_HarmonicStateReal(this) result(output)
  implicit none
  
  class(HarmonicStateReal), intent(in) :: this
  type(String), allocatable        :: output(:)
  
  select type(this); type is(HarmonicStateReal)
    output = [ 'Subspace  : '//this%subspace_id,    &
             & 'Frequency : '//this%frequency,      &
             & str('State'),                        &
             & join(str(this%modes_), delimiter='') ]
  class default
    call err()
  end select
end function

function new_HarmonicStateReal_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(HarmonicStateReal)  :: this
  
  call this%read(input)
end function

impure elemental function new_HarmonicStateReal_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(HarmonicStateReal)       :: this
  
  this = HarmonicStateReal(str(input))
end function
end module
