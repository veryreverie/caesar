! ======================================================================
! A harmonic state along a mode at q and its -q pair.
!    |n_i,n_j> ~ (a_i^n_i)(a_j^n_j)|0_i,0_j>.
! ======================================================================
! N.B. throughout, f(p) is the odd factorial of p, f(p) = (2p)!/(p! 2^p).
! N.B. for speed, any factors of (2Nw) are neglected here,
!    and should be added elsewhere.
!
! If a mode is not its own conjugate, (u_i)* = u_j,
!    then the double-mode states along modes u_i and u_j are:
! |p_i,p_j> = 1/sqrt(p_i!p_j!) a_i^p_i a_j^p_j |0_i,0_j>,
!    where a_i is the harmonic creation operator along mode i.
! |0_i,0_j> = sqrt(2*m*w/pi) exp(- N w |u_i|^2 )
!
! States are orthonormalised, <p_i,p_j|p_i,p_j>=1 and
!    <p_i,p_j|q_i,q_j>=0 if p_i/=q_i or p_j/=q_j.
!
! w is the effective frequency the pair of modes
!    (in general not the same as the harmonic frequency).
!
! m is the geometric average mass, arising as a result of mass-reduction
!    of co-ordinates. This cancels in integrals, and so is not stored.
!
! N is the number of primitive cells in the anharmonic supercell.
module harmonic_state_2d_module
  use common_module
  
  use anharmonic_common_module
  implicit none
  
  private
  
  public :: HarmonicState2D
  
  type, extends(Stringable) :: HarmonicState2D
    integer, private :: id_
    integer, private :: paired_id_
    integer, private :: occupation_
    integer, private :: paired_occupation_
  contains
    procedure, public :: id                => id_HarmonicState2D
    procedure, public :: paired_id         => paired_id_HarmonicState2D
    procedure, public :: occupation        => occupation_HarmonicState2D
    procedure, public :: paired_occupation => paired_occupation_HarmonicState2D
    procedure, public :: total_occupation  => total_occupation_HarmonicState2D
    procedure, public :: wavevector        => wavevector_HarmonicState2D
    
    ! ------------------------------
    ! Objects of the form <p_i,p_j|X|q_i,q_j>
    ! ------------------------------
    ! Whether <p_i,p_j|q_i,q_j> is finite.
    procedure, public :: finite_overlap => &
                       & finite_overlap_HarmonicState2D
    ! <p_i,p_j|q_i,q_j>.
    procedure, public :: inner_product => &
                       & inner_product_HarmonicState2D
    ! <p_i,p_j|(u_i)^(n_i)(u_j)^(n_j)|q_i,q_j>.
    procedure, public :: braket => &
                       & braket_HarmonicState2D
    procedure, public :: log_braket => &
                       & log_braket_HarmonicState2D
    ! <p_i,p_j|d/d(u_i)|q_i,q_j>.
    procedure, public :: plus_derivative => &
                       & plus_derivative_HarmonicState2D
    ! <p_i,p_j|d/d(u_j)|q_i,q_j>.
    procedure, public :: minus_derivative => &
                       & minus_derivative_HarmonicState2D
    ! <p_i,p_j|d2/d(u_i)d(u_j)|q_i,q_j>.
    procedure, public :: second_derivative => &
                       & second_derivative_HarmonicState2D
    
    ! ------------------------------
    ! I/O.
    ! ------------------------------
    procedure, public :: read  => read_HarmonicState2D
    procedure, public :: write => write_HarmonicState2D
  end type
  
  interface HarmonicState2D
    module procedure new_HarmonicState2D
    module procedure new_HarmonicState2D_String
  end interface
contains

! Constructor.
impure elemental function new_HarmonicState2D(id,paired_id,occupation, &
   & paired_occupation) result(this)
  implicit none
  
  integer, intent(in)   :: id
  integer, intent(in)   :: paired_id
  integer, intent(in)   :: occupation
  integer, intent(in)   :: paired_occupation
  type(HarmonicState2D) :: this
  
  if (id==paired_id) then
    call print_line(CODE_ERROR//': Mode is its own pair.')
    call err()
  endif
  
  this%id_                = id
  this%paired_id_         = paired_id
  this%occupation_        = occupation
  this%paired_occupation_ = paired_occupation
end function

! Getters.
impure elemental function id_HarmonicState2D(this) result(output)
  implicit none
  
  class(HarmonicState2D), intent(in) :: this
  integer                            :: output
  
  output = this%id_
end function

impure elemental function paired_id_HarmonicState2D(this) result(output)
  implicit none
  
  class(HarmonicState2D), intent(in) :: this
  integer                            :: output
  
  output = this%paired_id_
end function

impure elemental function occupation_HarmonicState2D(this) result(output)
  implicit none
  
  class(HarmonicState2D), intent(in) :: this
  integer                            :: output
  
  output = this%occupation_
end function

impure elemental function paired_occupation_HarmonicState2D(this) &
   & result(output)
  implicit none
  
  class(HarmonicState2D), intent(in) :: this
  integer                            :: output
  
  output = this%paired_occupation_
end function

! Returns the total occupation of the state.
impure elemental function total_occupation_HarmonicState2D(this) result(output)
  implicit none
  
  class(HarmonicState2D), intent(in) :: this
  integer                            :: output
  
  output = this%occupation_ + this%paired_occupation_
end function

! Returns the wavevector of the state.
function wavevector_HarmonicState2D(this,modes,qpoints) result(output)
  implicit none
  
  class(HarmonicState2D), intent(in) :: this
  type(ComplexMode),      intent(in) :: modes(:)
  type(QpointData),       intent(in) :: qpoints(:)
  type(FractionVector)               :: output
  
  integer :: mode
  integer :: qpoint
  
  mode   = first(modes%id==this%id_)
  qpoint = first(qpoints%id==modes(mode)%qpoint_id)
  output = qpoints(qpoint)%qpoint * (this%occupation_-this%paired_occupation_)
end function

! ----------------------------------------------------------------------
! Objects of the form <p|X|q>, for a variety of operators X.
! ----------------------------------------------------------------------

! Whether <p_i,p_j|q_i,q_j> is non-zero.
!    = true    if p_i=q_i and p_j=q_j.
!    = false   otherwise.
impure elemental function finite_overlap_HarmonicState2D(bra,ket) &
   & result(output)
  implicit none
  
  class(HarmonicState2D), intent(in) :: bra
  class(HarmonicState2D), intent(in) :: ket
  logical                            :: output
  
  integer :: p_i,p_j,q_i,q_j
  
  p_i = bra%occupation_
  p_j = bra%paired_occupation_
  q_i = ket%occupation_
  q_j = ket%paired_occupation_
  
  output = p_i==q_i .and. p_j==q_j
end function

! <p_i,p_j|q_i,q_j>.
!    = 1   if p_i=q_i and p_j=q_j.
!    = 0   otherwise.
impure elemental function inner_product_HarmonicState2D(bra,ket) &
   & result(output)
  implicit none
  
  class(HarmonicState2D), intent(in) :: bra
  class(HarmonicState2D), intent(in) :: ket
  integer                            :: output
  
  integer :: p_i,p_j,q_i,q_j
  
  p_i = bra%occupation_
  p_j = bra%paired_occupation_
  q_i = ket%occupation_
  q_j = ket%paired_occupation_
  
  if (p_i==q_i .and. p_j==q_j) then
    output = 1
  else
    output = 0
  endif
end function

! <p_i,p_j|(u_i)^(n_i)(u_j)^(n_j)|q_i,q_j>.
!    = 0                                 if p_i-p_j-n_i+n_j-q_i+q_j /= 0.
!    = 0                                 if n_i < p_i-q_i.
!    = 0                                 if n_i < q_j-p_j.
!    = 0                                 if n_j < p_j-q_j.
!    = 0                                 if n_j < q_i-p_i.
!    = 1/sqrt(2Nw)^{n_i+n_j}
!    * sqrt(p_i!q_i!/(p_j!q_j!))
!    * sum_{k=max(0,p_i-n_i,q_i-n_j)}^{min(p_i,q_i)}[
!                                    binom(n_i,p_i-k)
!                                  * binom(n_j,q_i-k)
!                                  * (n_i+p_j-p_i+k)!
!                                  / k!               ] otherwise.
impure elemental function braket_HarmonicState2D(bra,ket,potential,log_2nw, &
   & maximum_power,expansion_order) result(output) 
  implicit none
  
  class(HarmonicState2D),  intent(in) :: bra
  class(HarmonicState2D),  intent(in) :: ket
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

impure elemental function log_braket_HarmonicState2D(bra,ket,potential, &
   & log_2nw,maximum_power,expansion_order) result(output) 
  implicit none
  
  class(HarmonicState2D),  intent(in) :: bra
  class(HarmonicState2D),  intent(in) :: ket
  type(ComplexUnivariate), intent(in) :: potential
  real(dp),                intent(in) :: log_2nw
  integer,                 intent(in) :: maximum_power
  integer,                 intent(in) :: expansion_order
  real(dp)                            :: output
  
  real(dp), allocatable, save :: cache(:,:,:,:,:,:)
  logical,  allocatable, save :: cached(:,:,:,:,:,:)
  
  integer :: p_i,p_j,q_i,q_j,n_i,n_j
  
  integer :: k,ialloc
  
  if (.not. allocated(cache)) then
    allocate( cache( 0               :expansion_order,     &
            &        0               :expansion_order,     &
            &        -expansion_order:expansion_order,     &
            &        -expansion_order:expansion_order,     &
            &        0               :maximum_power,       &
            &        0               :maximum_power    ),  &
            & cached( 0               :expansion_order,    &
            &         0               :expansion_order,    &
            &         -expansion_order:expansion_order,    &
            &         -expansion_order:expansion_order,    &
            &         0               :maximum_power,      &
            &         0               :maximum_power    ), &
            & stat=ialloc); call err(ialloc)
    cached = .false.
  endif
  
  p_i = bra%occupation_
  p_j = bra%paired_occupation_
  q_i = ket%occupation_
  q_j = ket%paired_occupation_
  n_i = potential%power
  n_j = potential%paired_power
  
  if (p_i-p_j-n_i+n_j-q_i+q_j/=0 .or. n_i<p_i-q_i .or. n_i<q_j-p_j) then
    output = -1e100_dp
  else
    ! Rather than calculating the result directly, the result is calculated
    !    and cached the first time around, and then the cached result is used
    !    thereafter.
    if (.not.cached(n_i,n_j,p_i-q_i,p_j-q_j,p_i,p_j)) then
      cache(n_i,n_j,p_i-q_i,p_j-q_j,p_i,p_j) =               &
         &   0.5_dp*( log_factorial(p_i)                     &
         &          + log_factorial(q_i)                     &
         &          - log_factorial(p_j)                     &
         &          - log_factorial(q_j) )                   &
         & + log(sum([( exp( log_binomial(n_i,p_i-k)         &
         &                 + log_binomial(n_j,q_i-k)         &
         &                 + log_factorial(n_i+p_j-p_i+k)    &
         &                 - log_factorial(k)             ), &
         &              k=max(0,p_i-n_i,q_i-n_j),            &
         &              min(p_i,q_i)                         )]))
      cached(n_i,n_j,p_i-q_i,p_j-q_j,p_i,p_j) = .true.
    endif
    output = cache(n_i,n_j,p_i-q_i,p_j-q_j,p_i,p_j) &
         & -0.5_dp*(n_i+n_j)*log_2nw
  endif
end function

! <p_i,p_j|d/d(u_i)|q_i,q_j>.
!    = -sqrt(2Nw) * sqrt(p_i)/2   if p_i-q_i=1 and p_j=q_j.
!    =  sqrt(2Nw) * sqrt(q_j)/2   if q_j-p_j=1 and p_i=q_i.
!    =  0                         otherwise.
!
! N.B. the factor of sqrt(2Nw) is neglected.
impure elemental function plus_derivative_HarmonicState2D(bra,ket) &
   & result(output)
  implicit none
  
  class(HarmonicState2D), intent(in) :: bra
  class(HarmonicState2D), intent(in) :: ket
  real(dp)                           :: output
  
  integer :: p_i,p_j,q_i,q_j
  
  p_i = bra%occupation_
  p_j = bra%paired_occupation_
  q_i = ket%occupation_
  q_j = ket%paired_occupation_
  
  if (p_i-q_i==1 .and. p_j==q_j) then
    output = -sqrt(p_i/4.0_dp)
  elseif (q_j-p_j==1 .and. p_i==q_i) then
    output = sqrt(q_j/4.0_dp)
  else
    output = 0
  endif
end function

! <p_i,p_j|d/d(u_j)|q_i,q_j>.
!    = -sqrt(2Nw) * sqrt(p_j)/2   if p_j-q_j=1 and p_i=q_i.
!    =  sqrt(2Nw) * sqrt(q_i)/2   if q_i-p_i=1 and p_j=q_j.
!    =  0                         otherwise.
!
! N.B. the factor of sqrt(2Nw) is neglected.
impure elemental function minus_derivative_HarmonicState2D(bra,ket) &
   & result(output)
  implicit none
  
  class(HarmonicState2D), intent(in) :: bra
  class(HarmonicState2D), intent(in) :: ket
  real(dp)                           :: output
  
  integer :: p_i,p_j,q_i,q_j
  
  p_i = bra%occupation_
  p_j = bra%paired_occupation_
  q_i = ket%occupation_
  q_j = ket%paired_occupation_
  
  if (p_j-q_j==1 .and. p_i==q_i) then
    output = -sqrt(p_j/4.0_dp)
  elseif (q_i-p_i==1 .and. p_j==q_j) then
    output = sqrt(q_i/4.0_dp)
  else
    output = 0
  endif
end function

! <p_i,p_j|d2/d(u_i)d(u_j)|q_i,q_j>.
!    = -2Nw * (p_i+p_j+1)/4     if p_i=q_i and p_j=q_j.
!    =  2Nw * sqrt(p_i*p_j)/4   if p_i-q_i=p_j-q_j=1.
!    =  2Nw * sqrt(q_i*q_j)/4   if q_i-p_i=q_j-p_j=1.
!    =  0                       otherwise.
!
! N.B. the factor of 2Nw is neglected.
impure elemental function second_derivative_HarmonicState2D(bra,ket) &
   & result(output)
  implicit none
  
  class(HarmonicState2D), intent(in) :: bra
  class(HarmonicState2D), intent(in) :: ket
  real(dp)                           :: output
  
  integer :: p_i,p_j,q_i,q_j
  
  p_i = bra%occupation_
  p_j = bra%paired_occupation_
  q_i = ket%occupation_
  q_j = ket%paired_occupation_
  
  if (p_i==q_i .and. p_j==q_j) then
    output = -(p_i+p_j+1)/4.0_dp
  elseif (p_i-q_i==1 .and. p_j-q_j==1) then
    output = sqrt(p_i*p_j/16.0_dp)
  elseif (q_i-p_i==1 .and. q_j-p_j==1) then
    output = sqrt(q_i*q_j/16.0_dp)
  else
    output = 0
  endif
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_HarmonicState2D(this,input)
  implicit none
  
  class(HarmonicState2D), intent(out) :: this
  type(String),           intent(in)  :: input
  
  type(String), allocatable :: line(:)
  
  select type(this); type is(HarmonicState2D)
    ! input = '|ai^j,ak^l>', where i, j, k and l are integers.
    ! First remove the '|' and '>', and then split by the comma,
    !    to give ['ai^j','ak^l'].
    line = split_line(slice(input,2,len(input)-1),delimiter=',')
    
    ! Split ['ai^j','ak^l'] into ['i','j','k','l'].
    line = [ split_line(slice(line(1),2,len(line(1))),delimiter='^'), &
           & split_line(slice(line(2),2,len(line(2))),delimiter='^')  ]
    
    ! Construct the state. id=i, occupation=j, paired_id=k, paired_occupation=l.
    this = HarmonicState2D( id                = int(line(1)), &
                          & occupation        = int(line(2)), &
                          & paired_id         = int(line(3)), &
                          & paired_occupation = int(line(4))  )
  class default
    call err()
  end select
end subroutine

function write_HarmonicState2D(this) result(output)
  implicit none
  
  class(HarmonicState2D), intent(in) :: this
  type(String)                       :: output
  
  select type(this); type is(HarmonicState2D)
    output = '|a'                                          // &
           &  this%id_//'^'//this%occupation_              // &
           & ',a'                                          // &
           & this%paired_id_//'^'//this%paired_occupation_ // &
           & '>'
  class default
    call err()
  end select
end function

impure elemental function new_HarmonicState2D_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(HarmonicState2D)    :: this
  
  call this%read(input)
end function
end module
