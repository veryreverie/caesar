! ======================================================================
! A monomial state along a single mode at 2q=G,
!    |p> ~ u^p|0>.
! ======================================================================
! N.B. throughout, f(p) is the odd factorial of p, f(p) = (2p)!/(p! 2^p).
! N.B. for speed, any factors of (2Nw) are neglected here,
!    and should be added elsewhere.
!
! If a mode is its own conjugate, (u)* = u,
!    then the single-mode states along mode u are:
! |p> = sqrt((2Nw)^(p) / f(p)) (u)^(p) |0>
! The wavefunction of |0> is sqrt(sqrt(m*w/pi)) exp(- 1/2 N w (u)^2 )
!
! States are normalised <p|p>=1,
!    but are in general not orthogonal <p|q>/=0.
!
! w is the effective frequency of the mode
!    (in general not the same as the harmonic frequency).
!
! m is the geometric average mass, arising as a result of mass-reduction
!    of co-ordinates. This cancels in integrals, and so is not stored.
!
! N is the number of primitive cells in the anharmonic supercell.
module monomial_state_1d_module
  use common_module
  
  use anharmonic_common_module
  implicit none
  
  private
  
  public :: MonomialState1D
  
  type, extends(Stringable) :: MonomialState1D
    integer, private :: id_
    integer, private :: power_
  contains
    procedure, public :: id           => id_MonomialState1D
    procedure, public :: power        => power_MonomialState1D
    procedure, public :: total_power  => total_power_MonomialState1D
    procedure, public :: wavevector   => wavevector_MonomialState1D
    
    ! ------------------------------
    ! Objects of the form <p|X|q>
    ! ------------------------------
    ! Whether <p|q> is finite.
    procedure, public :: finite_overlap => &
                       & finite_overlap_MonomialState1D
    ! <p|q>.
    procedure, public :: inner_product => &
                       & inner_product_MonomialState1D
    ! <p|u^n|q>.
    procedure, public :: braket => &
                       & braket_MonomialState1D
    ! <p|d/du|q>.
    procedure, public :: first_derivative => &
                       & first_derivative_MonomialState1D
    ! <p|d2/du2|q>.
    procedure, public :: second_derivative => &
                       & second_derivative_MonomialState1D
    
    ! ------------------------------
    ! I/O.
    ! ------------------------------
    procedure, public :: read  => read_MonomialState1D
    procedure, public :: write => write_MonomialState1D
  end type
  
  interface MonomialState1D
    module procedure new_MonomialState1D
    module procedure new_MonomialState1D_String
  end interface
contains

! Constructor.
impure elemental function new_MonomialState1D(id,power) result(this)
  implicit none
  
  integer, intent(in)   :: id
  integer, intent(in)   :: power
  type(MonomialState1D) :: this
  
  this%id_    = id
  this%power_ = power
end function

! Getters.
impure elemental function id_MonomialState1D(this) result(output)
  implicit none
  
  class(MonomialState1D), intent(in) :: this
  integer                            :: output
  
  output = this%id_
end function

impure elemental function power_MonomialState1D(this) result(output)
  implicit none
  
  class(MonomialState1D), intent(in) :: this
  integer                            :: output
  
  output = this%power_
end function

! Returns the total power of the state.
impure elemental function total_power_MonomialState1D(this) result(output)
  implicit none
  
  class(MonomialState1D), intent(in) :: this
  integer                            :: output
  
  output = this%power_
end function

! Returns the wavevector of the state.
function wavevector_MonomialState1D(this,modes,qpoints) result(output)
  implicit none
  
  class(MonomialState1D), intent(in) :: this
  type(ComplexMode),      intent(in) :: modes(:)
  type(QpointData),       intent(in) :: qpoints(:)
  type(FractionVector)               :: output
  
  integer :: mode
  integer :: qpoint
  
  mode   = first(modes%id==this%id_)
  qpoint = first(qpoints%id==modes(mode)%qpoint_id)
  output = qpoints(qpoint)%qpoint * this%power_
end function

! ----------------------------------------------------------------------
! Objects of the form <p|X|q>, for a variety of operators X.
! ----------------------------------------------------------------------

! Whether <p|q> is non-zero.
!    = false  if p+q odd.
!    = true   otherwise.
impure elemental function finite_overlap_MonomialState1D(bra,ket) &
   & result(output)
  implicit none
  
  class(MonomialState1D), intent(in) :: bra
  class(MonomialState1D), intent(in) :: ket
  logical                            :: output
  
  integer :: p,q
  
  p = bra%power_
  q = ket%power_
  
  output = modulo(p+q, 2)==0
end function

! <p|q>.
!    = 0                              if p+q odd.
!    = f((p+q)/2) / sqrt(f(p)*f(q))   otherwise,
! where f(x) is the odd factorial of x, f(x) = prod_{k=1}^x [ 2k-1 ].
impure elemental function inner_product_MonomialState1D(bra,ket) &
   & result(output)
  implicit none
  
  class(MonomialState1D), intent(in) :: bra
  class(MonomialState1D), intent(in) :: ket
  real(dp)                           :: output
  
  integer :: p,q
  
  p = bra%power_
  q = ket%power_
  
  if (modulo(p+q,2)==1) then
    output = 0
  else
    output = exp( log_odd_factorial((p+q)/2)                         &
              & - 0.5_dp*(log_odd_factorial(p)+log_odd_factorial(q)) )
  endif
end function

! <p|(u)^(n)|q>.
!     = 0                 if p+n+q odd.
!     = 1/sqrt(2Nw)^{n}
!     * f((p+n+q)/2)
!     / sqrt(f(p)*f(q))   otherwise,
! where f(x) is the odd factorial of x, f(x) = prod_{k=1}^x [ 2k-1 ].
impure elemental function braket_MonomialState1D(bra,ket,potential,log_2nw) &
   & result(output)
  implicit none
  
  class(MonomialState1D),  intent(in) :: bra
  class(MonomialState1D),  intent(in) :: ket
  type(ComplexUnivariate), intent(in) :: potential
  real(dp),                intent(in) :: log_2nw
  real(dp)                            :: output
  
  integer :: p,q,n
  
  p = bra%power_
  q = ket%power_
  n = potential%power
  
  if (modulo(p+n+q,2)==1) then
    output = 0
  else
    output = exp( log_odd_factorial((p+n+q)/2)                       &
              & - 0.5_dp*(log_odd_factorial(p)+log_odd_factorial(q)) &
              & - 0.5_dp*n*log_2nw                                   )
  endif
end function

! <p|d/du|q>.
!    = 0                 if p+q even.
!    = sqrt(2Nw)
!    * (q-p)/2
!    * f((p+q-1)/2)
!    / sqrt(f(p)*f(q))   otherwise,
! where f(x) is the odd factorial of x, f(x) = prod_{k=1}^x [ 2k-1 ].
!
! N.B. the factor of sqrt(2Nw) is neglected.
impure elemental function first_derivative_MonomialState1D(bra,ket) &
   & result(output)
  implicit none
  
  class(MonomialState1D), intent(in) :: bra
  class(MonomialState1D), intent(in) :: ket
  real(dp)                           :: output
  
  integer :: p,q
  
  p = bra%power_
  q = ket%power_
  
  if (modulo(p+q,2)==0) then
    output = 0
  else
    output = ((q-p)/2)                                               &
         & * exp( log_odd_factorial((p+q-1)/2)                       &
         &      - 0.5_dp*(log_odd_factorial(p)+log_odd_factorial(q)) )
  endif
end function

! <p|d2/du2|q>.
!    = 0                              if p+q odd.
!    = 2Nw
!    * (((p-q)^2-1)/(2(p+q-1)) - 1)
!    * f((p+q)/2) / sqrt(f(p)*f(q))   otherwise,
! where f(x) is the odd factorial of x, f(x) = prod_{k=1}^x [ 2k-1 ].
!
! N.B. the factor of 2Nw is neglected.
impure elemental function second_derivative_MonomialState1D(bra,ket) &
   & result(output)
  implicit none
  
  class(MonomialState1D), intent(in) :: bra
  class(MonomialState1D), intent(in) :: ket
  real(dp)                           :: output
  
  integer :: p,q
  
  p = bra%power_
  q = ket%power_
  
  if (modulo(p+q,2)==1) then
    output = 0
  else
    output = (((p-q)*(p-q)-1)/(2.0_dp*(p+q-1)) - 1)                  &
         & * exp( log_odd_factorial((p+q)/2)                         &
         &      - 0.5_dp*(log_odd_factorial(p)+log_odd_factorial(q)) )
  endif
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_MonomialState1D(this,input)
  implicit none
  
  class(MonomialState1D), intent(out) :: this
  type(String),           intent(in)  :: input
  
  type(String), allocatable :: line(:)
  
  select type(this); type is(MonomialState1D)
    ! input = '|ui^j>', where i and j are integers.
    ! First remove the '|' and '>',
    !    and then split by the '^', to give ['ui^j'] 
    line = split_line(slice(input,2,len(input)-1),delimiter='^')
    
    ! Construct the state. id=i, power=j.
    this = MonomialState1D( id    = int(line(1)), &
                          & power = int(line(2))  )
  class default
    call err()
  end select
end subroutine

function write_MonomialState1D(this) result(output)
  implicit none
  
  class(MonomialState1D), intent(in) :: this
  type(String)                       :: output
  
  select type(this); type is(MonomialState1D)
    output = '|u'//this%id_//'^'//this%power_//'>'
  class default
    call err()
  end select
end function

impure elemental function new_MonomialState1D_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(MonomialState1D)    :: this
  
  call this%read(input)
end function
end module
