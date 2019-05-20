! ======================================================================
! A monomial state along a mode at q and its -q pair,
!    |p_i,p_j> ~ (u_i^p_i)(u_j^p_j)|0_i,0_j>.
! ======================================================================
! N.B. throughout, f(p) is the odd factorial of p, f(p) = (2p)!/(p! 2^p).
! N.B. for speed, any factors of (2Nw) are neglected here,
!    and should be added elsewhere.
!
! If a mode is not its own conjugate, (u_i)* = u_j,
!    then the double-mode states along modes u_i and u_j are:
! |p_i,p_j> = sqrt((2Nw)^(p_i+p_j) / (p_i+p_j)!) u_i^p_i u_j^p_j |0_i,0_j>
! |0_i,0_j> = sqrt(2*m*w/pi) exp(- N w |u_i|^2 )
!
! States are normalised <p_i,p_j|p_i,p_j>=1,
!    but are in general not orthogonal <p_i,p_j|q_i,q_j>/=0.
!
! w is the effective frequency the pair of modes
!    (in general not the same as the harmonic frequency).
!
! m is the geometric average mass, arising as a result of mass-reduction
!    of co-ordinates. This cancels in integrals, and so is not stored.
!
! N is the number of primitive cells in the anharmonic supercell.
module monomial_state_2d_module
  use common_module
  
  use anharmonic_common_module
  implicit none
  
  private
  
  public :: MonomialState2D
  
  type, extends(Stringable) :: MonomialState2D
    integer, private :: id_
    integer, private :: paired_id_
    integer, private :: power_
    integer, private :: paired_power_
  contains
    procedure, public :: id           => id_MonomialState2D
    procedure, public :: paired_id    => paired_id_MonomialState2D
    procedure, public :: power        => power_MonomialState2D
    procedure, public :: paired_power => paired_power_MonomialState2D
    procedure, public :: total_power  => total_power_MonomialState2D
    procedure, public :: wavevector   => wavevector_MonomialState2D
    
    ! ------------------------------
    ! Objects of the form <p_i,p_j|X|q_i,q_j>
    ! ------------------------------
    ! Whether <p_i,p_j|q_i,q_j> is finite.
    procedure, public :: finite_overlap => &
                       & finite_overlap_MonomialState2D
    ! <p_i,p_j|q_i,q_j>.
    procedure, public :: inner_product => &
                       & inner_product_MonomialState2D
    ! <p_i,p_j|(u_i)^(n_i)(u_j)^(n_j)|q_i,q_j>.
    procedure, public :: braket => &
                       & braket_MonomialState2D
    ! <p_i,p_j|d/d(u_i)|q_i,q_j>.
    procedure, public :: plus_derivative => &
                       & plus_derivative_MonomialState2D
    ! <p_i,p_j|d/d(u_j)|q_i,q_j>.
    procedure, public :: minus_derivative => &
                       & minus_derivative_MonomialState2D
    ! <p_i,p_j|d2/d(u_i)d(u_j)|q_i,q_j>.
    procedure, public :: second_derivative => &
                       & second_derivative_MonomialState2D
    
    ! ------------------------------
    ! I/O.
    ! ------------------------------
    procedure, public :: read  => read_MonomialState2D
    procedure, public :: write => write_MonomialState2D
  end type
  
  interface MonomialState2D
    module procedure new_MonomialState2D
    module procedure new_MonomialState2D_String
  end interface
contains

! Constructor.
impure elemental function new_MonomialState2D(id,paired_id,power, &
   & paired_power) result(this)
  implicit none
  
  integer, intent(in)   :: id
  integer, intent(in)   :: paired_id
  integer, intent(in)   :: power
  integer, intent(in)   :: paired_power
  type(MonomialState2D) :: this
  
  if (id==paired_id) then
    call print_line(CODE_ERROR//': Mode is its own pair.')
    call err()
  endif
  
  this%id_           = id
  this%paired_id_    = paired_id
  this%power_        = power
  this%paired_power_ = paired_power
end function

! Getters.
impure elemental function id_MonomialState2D(this) result(output)
  implicit none
  
  class(MonomialState2D), intent(in) :: this
  integer                            :: output
  
  output = this%id_
end function

impure elemental function paired_id_MonomialState2D(this) result(output)
  implicit none
  
  class(MonomialState2D), intent(in) :: this
  integer                            :: output
  
  output = this%paired_id_
end function

impure elemental function power_MonomialState2D(this) result(output)
  implicit none
  
  class(MonomialState2D), intent(in) :: this
  integer                            :: output
  
  output = this%power_
end function

impure elemental function paired_power_MonomialState2D(this) result(output)
  implicit none
  
  class(MonomialState2D), intent(in) :: this
  integer                            :: output
  
  output = this%paired_power_
end function

! Returns the total power of the state.
impure elemental function total_power_MonomialState2D(this) result(output)
  implicit none
  
  class(MonomialState2D), intent(in) :: this
  integer                            :: output
  
  output = this%power_ + this%paired_power_
end function

! Returns the wavevector of the state.
function wavevector_MonomialState2D(this,modes,qpoints) result(output)
  implicit none
  
  class(MonomialState2D), intent(in) :: this
  type(ComplexMode),      intent(in) :: modes(:)
  type(QpointData),       intent(in) :: qpoints(:)
  type(FractionVector)               :: output
  
  integer :: mode
  integer :: qpoint
  
  mode   = first(modes%id==this%id_)
  qpoint = first(qpoints%id==modes(mode)%qpoint_id)
  output = qpoints(qpoint)%qpoint * (this%power_-this%paired_power_)
end function

! ----------------------------------------------------------------------
! Objects of the form <p|X|q>, for a variety of operators X.
! ----------------------------------------------------------------------

! Whether <p_i,p_j|q_i,q_j> is non-zero.
!    = false    if p_i-p_j-q_i+q_j /= 0.
!    = true     otherwise.
impure elemental function finite_overlap_MonomialState2D(bra,ket) &
   & result(output)
  implicit none
  
  class(MonomialState2D), intent(in) :: bra
  class(MonomialState2D), intent(in) :: ket
  logical                            :: output
  
  integer :: p_i,p_j,q_i,q_j
  
  p_i = bra%power_
  p_j = bra%paired_power_
  q_i = ket%power_
  q_j = ket%paired_power_
  
  output = p_i-p_j-q_i+q_j==0
end function

! <p_i,p_j|q_i,q_j>.
!    = 0                          if p_i-p_j-q_i+q_j /= 0.
!    = ((p_i+p_j+q_i+q_j)/2)!
!    / sqrt((p_i+p_j)!(q_i+q_j)!) otherwise.
impure elemental function inner_product_MonomialState2D(bra,ket) &
   & result(output)
  implicit none
  
  class(MonomialState2D), intent(in) :: bra
  class(MonomialState2D), intent(in) :: ket
  real(dp)                           :: output
  
  integer :: p_i,p_j,q_i,q_j
  
  p_i = bra%power_
  p_j = bra%paired_power_
  q_i = ket%power_
  q_j = ket%paired_power_
  
  if (p_i-p_j-q_i+q_j/=0) then
    output = 0.0_dp
  else
    output = exp( log_factorial((p_i+p_j+q_i+q_j)/2)                     &
              & - 0.5_dp*(log_factorial(p_i+p_j)+log_factorial(q_i+q_j)) )
  endif
end function

! <p_i,p_j|(u_i)^(n_i)(u_j)^(n_j)|q_i,q_j>.
!    = 0                                 if p_i-p_j-n_i+n_j-q_i+q_j /= 0.
!    = 1/sqrt(2Nw)^{n_i+n_j}
!    * ((p_i+p_j+n_i+n_j+q_i+q_j)/2)!
!    / sqrt( (p_i+p_j)! * (q_i+q_j)! )   otherwise.
!
! N.B. the factor of 1/sqrt(2Nw)^{n_i+n_j} is neglected.
impure elemental function braket_MonomialState2D(bra,ket,potential) &
   & result(output)
  implicit none
  
  class(MonomialState2D),  intent(in) :: bra
  class(MonomialState2D),  intent(in) :: ket
  type(ComplexUnivariate), intent(in) :: potential
  real(dp)                            :: output
  
  integer :: p_i,p_j,q_i,q_j,n_i,n_j
  
  p_i = bra%power_
  p_j = bra%paired_power_
  q_i = ket%power_
  q_j = ket%paired_power_
  n_i = potential%power
  n_j = potential%paired_power
  
  if (p_i-p_j-n_i+n_j-q_i+q_j/=0) then
    output = 0
  else
    output = exp( log_factorial((p_i+p_j+n_i+n_j+q_i+q_j)/2)             &
              & - 0.5_dp*(log_factorial(p_i+p_j)+log_factorial(q_i+q_j)) )
  endif
end function

! <p_i,p_j|d/d(u_i)|q_i,q_j>.
!    = 0                                       if q_i-p_i-q_j+p_j /= 1
!    = sqrt(2Nw)
!    * (q_i - sqrt((q+1)/(p+1))*(p+q+1)/4)
!    * ((p+q-1)/2)!/sqrt(p!q!)                 otherwise,
! where p=p_i+p_j and q=q_i+q_j
!
! N.B. the factor of sqrt(2Nw) is neglected.
impure elemental function plus_derivative_MonomialState2D(bra,ket) &
   & result(output)
  implicit none
  
  class(MonomialState2D), intent(in) :: bra
  class(MonomialState2D), intent(in) :: ket
  real(dp)                           :: output
  
  integer :: p_i,p_j,q_i,q_j
  integer :: p,q
  
  p_i = bra%power_
  p_j = bra%paired_power_
  q_i = ket%power_
  q_j = ket%paired_power_
  
  p = p_i+p_j
  q = q_i+q_j
  
  if (q_i-p_i-q_j+p_j/=1) then
    output = 0
  else
    output = (q_i - sqrt((q+1)/(p+1.0_dp))*(p+q+1)/4.0_dp)   &
         & * exp( log_factorial(p+q-1)                       &
         &      - 0.5_dp*(log_factorial(p)+log_factorial(q)) )
  endif
end function

! <p_i,p_j|d/d(u_j)|q_i,q_j>.
!    = 0                                       if q_j-p_j-q_i+p_i /= 1
!    = sqrt(2Nw)
!    * (q_j - sqrt((q+1)/(p+1))*(p+q+1)/4)
!    * ((p+q-1)/2)!/sqrt(p!q!)                 otherwise,
! where p=p_i+p_j and q=q_i+q_j
!
! N.B. the factor of sqrt(2Nw) is neglected.
impure elemental function minus_derivative_MonomialState2D(bra,ket) &
   & result(output)
  implicit none
  
  class(MonomialState2D), intent(in) :: bra
  class(MonomialState2D), intent(in) :: ket
  real(dp)                           :: output
  
  integer :: p_i,p_j,q_i,q_j
  integer :: p,q
  
  p_i = bra%power_
  p_j = bra%paired_power_
  q_i = ket%power_
  q_j = ket%paired_power_
  
  p = p_i+p_j
  q = q_i+q_j
  
  if (q_j-p_j-q_i+p_i/=1) then
    output = 0
  else
    output = (q_j - sqrt((q+1)/(p+1.0_dp))*(p+q+1)/4.0_dp)   &
         & * exp( log_factorial(p+q-1)                       &
         &      - 0.5_dp*(log_factorial(p)+log_factorial(q)) )
  endif
end function

! <p_i,p_j|d2/d(u_i)d(u_j)|q_i,q_j>.
!    = 0  if p_i-p_j-q_i+q_j /= 0.
!    = 2Nw
!    * ((p_i-q_j)(p_j-q_i)/(2(p+q))-1/4)
!    * ((p+q)/2)!/sqrt(p!q!)
! where p=p_i+p_j and q=q_i+q_j.
!
! N.B. the factor of 2Nw is neglected.
impure elemental function second_derivative_MonomialState2D(bra,ket) &
   & result(output)
  implicit none
  
  class(MonomialState2D), intent(in) :: bra
  class(MonomialState2D), intent(in) :: ket
  real(dp)                           :: output
  
  integer :: p_i,p_j,q_i,q_j
  integer :: p,q
  
  p_i = bra%power_
  p_j = bra%paired_power_
  q_i = ket%power_
  q_j = ket%paired_power_
  
  p = p_i+p_j
  q = q_i+q_j
  
  if (p_i-p_j-q_i+q_j/=0) then
    output = 0
  else
    output = ((p_i-q_j)*(p_j-q_i)/(2.0_dp*(p+q))-0.25_dp)    &
         & * exp( log_factorial((p+q)/2)                     &
         &      - 0.5_dp*(log_factorial(p)+log_factorial(q)) )
  endif
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_MonomialState2D(this,input)
  implicit none
  
  class(MonomialState2D), intent(out) :: this
  type(String),           intent(in)  :: input
  
  type(String), allocatable :: line(:)
  
  select type(this); type is(MonomialState2D)
    ! input = '|ui^j,uk^l>', where i, j, k and l are integers.
    ! First remove the '|' and '>', and then split by the comma,
    !    to give ['ui^j','uk^l'].
    line = split_line(slice(input,2,len(input)-1),delimiter=',')
    
    ! Split ['ui^j','uk^l'] into ['i','j','k','l'].
    line = [ split_line(slice(line(1),2,len(line(1))),delimiter='^'), &
           & split_line(slice(line(2),2,len(line(2))),delimiter='^')  ]
    
    ! Construct the state. id=i, power=j, paired_id=k, paired_power=l.
    this = MonomialState2D( id           = int(line(1)), &
                          & power        = int(line(2)), &
                          & paired_id    = int(line(3)), &
                          & paired_power = int(line(4))  )
  class default
    call err()
  end select
end subroutine

function write_MonomialState2D(this) result(output)
  implicit none
  
  class(MonomialState2D), intent(in) :: this
  type(String)                       :: output
  
  select type(this); type is(MonomialState2D)
    output = '|u'                                     // &
           &  this%id_//'^'//this%power_              // &
           & ',u'                                     // &
           & this%paired_id_//'^'//this%paired_power_ // &
           & '>'
  class default
    call err()
  end select
end function

impure elemental function new_MonomialState2D_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(MonomialState2D)    :: this
  
  call this%read(input)
end function
end module
