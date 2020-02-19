! ======================================================================
! A harmonic state along a single mode at 2q=G,
!    |p> = a^p|0>.
! ======================================================================
! N.B. throughout, f(p) is the odd factorial of p, f(p) = (2p)!/(p! 2^p).
! N.B. for speed, any factors of (2Nw) are neglected here,
!    and should be added elsewhere.
!
! If a mode is its own conjugate, (u)* = u,
!    then the single-mode states along mode u are:
! |p> = 1/sqrt(p!) (a)^(p) |0>,
!    where a is the harmonic creation operator.
! The wavefunction of |0> is sqrt(sqrt(m*w/pi)) exp(- 1/2 N w (u)^2 ).
!
! States are orthonormalised, <p|p>=1 and <p|q>=0 if p/=q.
!
! w is the effective frequency of the mode
!    (in general not the same as the harmonic frequency).
!
! m is the geometric average mass, arising as a result of mass-reduction
!    of co-ordinates. This cancels in integrals, and so is not stored.
!
! N is the number of primitive cells in the anharmonic supercell.
module harmonic_state_1d_module
  use common_module
  
  use anharmonic_common_module
  implicit none
  
  private
  
  public :: HarmonicState1D
  
  type, extends(Stringable) :: HarmonicState1D
    integer, private :: id_
    integer, private :: occupation_
  contains
    procedure, public :: id               => id_HarmonicState1D
    procedure, public :: occupation       => occupation_HarmonicState1D
    procedure, public :: total_occupation => total_occupation_HarmonicState1D
    procedure, public :: wavevector       => wavevector_HarmonicState1D
    
    ! ------------------------------
    ! Objects of the form <p|X|q>
    ! ------------------------------
    ! Whether <p|q> is finite.
    procedure, public :: finite_overlap => &
                       & finite_overlap_HarmonicState1D
    ! <p|q>.
    procedure, public :: inner_product => &
                       & inner_product_HarmonicState1D
    ! <p|u^n|q>.
    procedure, public :: braket => &
                       & braket_HarmonicState1D
    procedure, public :: log_braket => &
                       & log_braket_HarmonicState1D
    ! <p|d/du|q>.
    procedure, public :: first_derivative => &
                       & first_derivative_HarmonicState1D
    ! <p|d2/du2|q>.
    procedure, public :: second_derivative => &
                       & second_derivative_HarmonicState1D
    
    ! ------------------------------
    ! I/O.
    ! ------------------------------
    procedure, public :: read  => read_HarmonicState1D
    procedure, public :: write => write_HarmonicState1D
  end type
  
  interface HarmonicState1D
    module procedure new_HarmonicState1D
    module procedure new_HarmonicState1D_String
  end interface
contains

! Constructor.
impure elemental function new_HarmonicState1D(id,occupation) result(this)
  implicit none
  
  integer, intent(in)   :: id
  integer, intent(in)   :: occupation
  type(HarmonicState1D) :: this
  
  this%id_         = id
  this%occupation_ = occupation
end function

! Getters.
impure elemental function id_HarmonicState1D(this) result(output)
  implicit none
  
  class(HarmonicState1D), intent(in) :: this
  integer                            :: output
  
  output = this%id_
end function

impure elemental function occupation_HarmonicState1D(this) result(output)
  implicit none
  
  class(HarmonicState1D), intent(in) :: this
  integer                            :: output
  
  output = this%occupation_
end function

! Returns the total occupation of the state.
impure elemental function total_occupation_HarmonicState1D(this) result(output)
  implicit none
  
  class(HarmonicState1D), intent(in) :: this
  integer                            :: output
  
  output = this%occupation_
end function

! Returns the wavevector of the state.
function wavevector_HarmonicState1D(this,modes,qpoints) result(output)
  implicit none
  
  class(HarmonicState1D), intent(in) :: this
  type(ComplexMode),      intent(in) :: modes(:)
  type(QpointData),       intent(in) :: qpoints(:)
  type(FractionVector)               :: output
  
  integer :: mode
  integer :: qpoint
  
  mode   = first(modes%id==this%id_)
  qpoint = first(qpoints%id==modes(mode)%qpoint_id)
  output = qpoints(qpoint)%qpoint * this%occupation_
end function

! ----------------------------------------------------------------------
! Objects of the form <p|X|q>, for a variety of operators X.
! ----------------------------------------------------------------------

! Whether <p|q> is non-zero.
!    = true   if p=q.
!    = false  otherwise.
impure elemental function finite_overlap_HarmonicState1D(bra,ket) &
   & result(output)
  implicit none
  
  class(HarmonicState1D), intent(in) :: bra
  class(HarmonicState1D), intent(in) :: ket
  logical                            :: output
  
  integer :: p,q
  
  p = bra%occupation_
  q = ket%occupation_
  
  output = p==q
end function

! <p|q>.
!    = 1   if p=q.
!    = 0   otherwise.
impure elemental function inner_product_HarmonicState1D(bra,ket) &
   & result(output)
  implicit none
  
  class(HarmonicState1D), intent(in) :: bra
  class(HarmonicState1D), intent(in) :: ket
  integer                            :: output
  
  integer :: p,q
  
  p = bra%occupation_
  q = ket%occupation_
  
  if (p==q) then
    output = 1
  else
    output = 0
  endif
end function

! <p|(u)^(n)|q>.  (Defining d=|p-q|).
!     = 0                 if p+n+q odd.
!     = 0                 if n<d.
!     = 1/sqrt(2Nw)^{n}
!     * (n!/(((n+d)/2)!2^((n-d)/2))
!     * sqrt(max(p,q)!/min(p,q)!)
!     * sum_{k=0}^{min((n-d)/2,p,q)}[ 2^k 
!                                   * binom((n+d)/2,d+k)
!                                   * binom(min(p,q),k)  ] otherwise.
impure elemental function braket_HarmonicState1D(bra,ket,potential,log_2nw, &
   & maximum_power,expansion_order) result(output) 
  implicit none
  
  class(HarmonicState1D),  intent(in) :: bra
  class(HarmonicState1D),  intent(in) :: ket
  type(ComplexUnivariate), intent(in) :: potential
  real(dp),                intent(in) :: log_2nw
  integer,                 intent(in) :: maximum_power
  integer,                 intent(in) :: expansion_order
  real(dp)                            :: output
  
  output = exp(bra%log_braket( ket,            &
                             & potential,      &
                             & log_2nw,        &
                             & maximum_power,  &
                             & expansion_order ))
end function

impure elemental function log_braket_HarmonicState1D(bra,ket,potential, &
   & log_2nw,maximum_power,expansion_order) result(output) 
  implicit none
  
  class(HarmonicState1D),  intent(in) :: bra
  class(HarmonicState1D),  intent(in) :: ket
  type(ComplexUnivariate), intent(in) :: potential
  real(dp),                intent(in) :: log_2nw
  integer,                 intent(in) :: maximum_power
  integer,                 intent(in) :: expansion_order
  real(dp)                            :: output
  
  real(dp), allocatable, save :: cache(:,:,:)
  logical,  allocatable, save :: cached(:,:,:)
  
  integer :: min_p,max_p,n,d
  
  integer :: k,ialloc
  
  if (.not. allocated(cache)) then
    allocate( cache(0:expansion_order,0:expansion_order,0:maximum_power),  &
            & cached(0:expansion_order,0:expansion_order,0:maximum_power), &
            & stat=ialloc); call err(ialloc)
    cached = .false.
  endif
  
  min_p = min(bra%occupation_, ket%occupation_)
  max_p = max(bra%occupation_, ket%occupation_)
  n = potential%power
  d = max_p-min_p
  
  if (modulo(min_p+max_p+n,2)==1 .or. n<d) then
    output = -1e100_dp
  else
    ! Rather than calculating the result directly, the result is calculated
    !    and cached the first time around, and then the cached result is used
    !    thereafter.
    if (.not. cached(n,d,min_p)) then
      cache(n,d,min_p) = log_factorial(n)                              &
                     & - log_factorial((n+d)/2)                        &
                     & - ((n-d)/2) * log(2.0_dp)                       &
                     & + 0.5_dp*( log_factorial(max_p)                 &
                     &          - log_factorial(min_p) )               &
                     & + log(sum([( exp( k*log(2.0_dp)                 &
                     &                 + log_binomial((n+d)/2, d+k)    &
                     &                 + log_binomial(min_p, k)     ), &
                     &              k=0,                               &
                     &              min(min_p,(n-d)/2)                 )]))
      cached(n,d,min_p) = .true.
    endif
    output = cache(n,d,min_p)-0.5_dp*n*log_2nw
  endif
end function

! <p|d/du|q>.
!    = -sqrt(2Nw) * sqrt(p)/2   if p=q+1.
!    =  sqrt(2Nw) * sqrt(q)/2   if q=p+1.
!    =  0                       otherwise.
!
! N.B. the factor of sqrt(2Nw) is neglected.
impure elemental function first_derivative_HarmonicState1D(bra,ket) &
   & result(output)
  implicit none
  
  class(HarmonicState1D), intent(in) :: bra
  class(HarmonicState1D), intent(in) :: ket
  real(dp)                           :: output
  
  integer :: p,q
  
  p = bra%occupation_
  q = ket%occupation_
  
  if (p==q+1) then
    output = -sqrt(p/4.0_dp)
  elseif (q==p+1) then
    output =  sqrt(q/4.0_dp)
  else
    output = 0
  endif
end function

! <p|d2/du2|q>.
!    = -2Nw (p+0.5)/2                       if p=q.
!    =  2Nw sqrt(max(p,q)*(max(p,q)-1))/4   if |p-q|=2.
!    =  0                                   otherwise.
!
! N.B. the factor of 2Nw is neglected.
impure elemental function second_derivative_HarmonicState1D(bra,ket) &
   & result(output)
  implicit none
  
  class(HarmonicState1D), intent(in) :: bra
  class(HarmonicState1D), intent(in) :: ket
  real(dp)                           :: output
  
  integer :: p,q
  
  p = bra%occupation_
  q = ket%occupation_
  
  if (p==q) then
    output = -(p+0.5_dp)/2
  elseif (abs(p-q)==2) then
    output = sqrt(max(p,q)*(max(p,q)-1.0_dp))/4
  else
    output = 0
  endif
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_HarmonicState1D(this,input)
  implicit none
  
  class(HarmonicState1D), intent(out) :: this
  type(String),           intent(in)  :: input
  
  type(String), allocatable :: line(:)
  
  select type(this); type is(HarmonicState1D)
    ! input = '|ai^j>', where i and j are integers.
    ! First remove the '|' and '>',
    !    and then split by the '^', to give ['ai^j'] 
    line = split_line(slice(input,2,len(input)-1),delimiter='^')
    
    ! Construct the state. id=i, occupation=j.
    this = HarmonicState1D( id         = int(line(1)), &
                          & occupation = int(line(2))  )
  class default
    call err()
  end select
end subroutine

function write_HarmonicState1D(this) result(output)
  implicit none
  
  class(HarmonicState1D), intent(in) :: this
  type(String)                       :: output
  
  select type(this); type is(HarmonicState1D)
    output = '|a'//this%id_//'^'//this%occupation_//'>'
  class default
    call err()
  end select
end function

impure elemental function new_HarmonicState1D_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(HarmonicState1D)    :: this
  
  call this%read(input)
end function
end module
