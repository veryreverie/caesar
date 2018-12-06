! ======================================================================
! The building blocks of basis functions in real co-ordinates.
! ======================================================================
module real_polynomial_module
  use utils_module
  
  use real_mode_module
  use real_single_mode_displacement_module
  use real_single_mode_force_module
  use real_mode_displacement_module
  use real_mode_force_module
  implicit none
  
  private
  
  public :: RealUnivariate
  public :: RealMonomial
  public :: RealPolynomial
  public :: size
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: operator(-)
  public :: sum
  public :: select_mode
  public :: select_modes
  public :: select_displacement
  public :: select_displacements
  public :: select_force
  public :: select_forces
  public :: compare_real_monomials
  
  ! --------------------------------------------------
  ! Types, and conversions between types.
  ! --------------------------------------------------
  
  ! It is desirable to be able to convert:
  !    univariates -> monomials -> polynomials
  ! This is acheived by extending from the classes ComplexPolynomialable and
  !    ComplexMonomialable, representing types which can be converted to
  !    ComplexPolynomial and ComplexMonomial respectively.
  
  type, abstract, extends(Stringable) :: RealPolynomialable
  contains
    procedure(to_RealPolynomial_RealPolynomialable), deferred, public :: &
       & to_RealPolynomial
  end type
  
  type, abstract, extends(RealPolynomialable) :: RealMonomialable
  contains
    procedure(to_RealMonomial_RealMonomialable), deferred, public :: &
       & to_RealMonomial
  end type
  
  type, extends(RealMonomialable) :: RealUnivariate
    integer :: id
    integer :: paired_id
    integer :: power
    integer :: paired_power
  contains
    procedure, public :: to_RealMonomial   => to_RealMonomial_RealUnivariate
    procedure, public :: to_RealPolynomial => to_RealPolynomial_RealUnivariate
    
    procedure, public :: total_power => total_power_RealUnivariate
    
    procedure, public :: energy => energy_RealUnivariate
    procedure, public :: force  => force_RealUnivariate
    
    procedure, public :: read  => read_RealUnivariate
    procedure, public :: write => write_RealUnivariate
  end type
  
  interface RealUnivariate
    module procedure new_RealUnivariate
    module procedure new_RealUnivariate_RealMode
    module procedure new_RealUnivariate_String
  end interface
  
  type, extends(RealMonomialable) :: RealMonomial
    real(dp)                                   :: coefficient
    ! Modes is private so that it can be guarenteed to be sorted.
    type(RealUnivariate), allocatable, private :: modes_(:)
  contains
    procedure, public :: to_RealMonomial   => to_RealMonomial_RealMonomial
    procedure, public :: to_RealPolynomial => to_RealPolynomial_RealMonomial
    
    procedure, public :: mode => mode_RealMonomial
    procedure, public :: modes => modes_RealMonomial
    procedure, public :: id  => id_RealMonomial
    procedure, public :: ids => ids_RealMonomial
    procedure, public :: paired_id  => paired_id_RealMonomial
    procedure, public :: paired_ids => paired_ids_RealMonomial
    procedure, public :: power  => power_RealMonomial
    procedure, public :: powers => powers_RealMonomial
    procedure, public :: paired_power  => paired_power_RealMonomial
    procedure, public :: paired_powers => paired_powers_RealMonomial
    
    procedure, public :: simplify => simplify_RealMonomial
    
    procedure, public :: total_power => total_power_RealMonomial
    
    procedure, public :: energy => energy_RealMonomial
    procedure, public :: force  => force_RealMonomial
    
    procedure, public :: read  => read_RealMonomial
    procedure, public :: write => write_RealMonomial
  end type
  
  interface RealMonomial
    module procedure new_RealMonomial
    module procedure new_RealMonomial_RealMonomialable
    module procedure new_RealMonomial_String
  end interface
  
  type, extends(RealPolynomialable) :: RealPolynomial
    type(RealMonomial), allocatable :: terms(:)
  contains
    procedure, public :: to_RealPolynomial => to_RealPolynomial_RealPolynomial
    
    procedure, public :: simplify => simplify_RealPolynomial
    
    procedure, public :: energy => energy_RealPolynomial
    procedure, public :: force  => force_RealPolynomial
    
    procedure, public :: read  => read_RealPolynomial
    procedure, public :: write => write_RealPolynomial
  end type
  
  interface RealPolynomial
    module procedure new_RealPolynomial
    module procedure new_RealPolynomial_RealPolynomialable
    module procedure new_RealPolynomial_String
  end interface
  
  abstract interface
    function to_RealPolynomial_RealPolynomialable(this) result(output)
      import RealPolynomial
      import RealPolynomialable
      implicit none
      
      class(RealPolynomialable), intent(in) :: this
      type(RealPolynomial)                  :: output
    end function
    
    function to_RealMonomial_RealMonomialable(this) result(output)
      import RealMonomial
      import RealMonomialable
      implicit none
      
      class(RealMonomialable), intent(in) :: this
      type(RealMonomial)                  :: output
    end function
  end interface
  
  ! --------------------------------------------------
  ! Operations involving types.
  ! --------------------------------------------------
  
  interface size
    module procedure size_RealMonomial
    module procedure size_RealPolynomial
  end interface
  
  interface operator(*)
    module procedure multiply_RealMonomial_real
    module procedure multiply_real_RealMonomial
    
    module procedure multiply_RealPolynomial_real
    module procedure multiply_real_RealPolynomial
    
    module procedure multiply_RealMonomialable_RealMonomialable
    
    module procedure multiply_RealPolynomial_RealMonomialable
    module procedure multiply_RealMonomialable_RealPolynomial
  end interface
  
  interface operator(/)
    module procedure divide_RealMonomial_real
    
    module procedure divide_RealPolynomial_real
  end interface
  
  interface operator(+)
    module procedure add_RealPolynomialable_RealPolynomialable
  end interface
  
  interface operator(-)
    module procedure negative_RealPolynomialable
    
    module procedure subtract_RealPolynomialable_RealPolynomialable
  end interface
  
  interface sum
    module procedure sum_RealPolynomialables
  end interface
  
  interface select_mode
    module procedure select_mode_RealUnivariate
  end interface
  
  interface select_modes
    module procedure select_modes_RealUnivariates
  end interface
  
  interface select_displacement
    module procedure select_displacement_RealUnivariate
  end interface
  
  interface select_displacements
    module procedure select_displacements_RealUnivariates
  end interface
  
  interface select_force
    module procedure select_force_RealUnivariate
  end interface
  
  interface select_forces
    module procedure select_forces_RealUnivariates
  end interface
contains

! ----------------------------------------------------------------------
! Constructors.
! ----------------------------------------------------------------------
impure elemental function new_RealUnivariate(id,paired_id,power, &
   & paired_power) result(this)
  implicit none
  
  integer, intent(in)  :: id
  integer, intent(in)  :: paired_id
  integer, intent(in)  :: power
  integer, intent(in)  :: paired_power
  type(RealUnivariate) :: this
  
  if (id<=paired_id) then
    if (id==paired_id .and. power/=paired_power) then
      call print_line(CODE_ERROR//': modes the same, but powers differ.')
      call err()
    endif
    this%id           = id
    this%paired_id    = paired_id
    this%power        = power
    this%paired_power = paired_power
  else
    this%id           = paired_id
    this%paired_id    = id
    this%power        = paired_power
    this%paired_power = power
  endif
end function

function new_RealUnivariate_RealMode(mode,power,paired_power) result(this)
  implicit none
  
  type(RealMode), intent(in)           :: mode
  integer,        intent(in)           :: power
  integer,        intent(in), optional :: paired_power
  type(RealUnivariate)                 :: this
  
  if (present(paired_power)) then
    if (mode%id==mode%paired_id .and. power/=paired_power) then
      call print_line(ERROR//': Mode is its own pair, but power does not &
         & match paired_power.')
      call err()
    endif
    this = RealUnivariate( id           = mode%id,        &
                         & paired_id    = mode%paired_id, &
                         & power        = power,          &
                         & paired_power = paired_power    )
  else
    if (mode%id==mode%paired_id) then
      this = RealUnivariate( id           = mode%id, &
                           & paired_id    = mode%id, &
                           & power        = power,   &
                           & paired_power = power    )
    else
      this = RealUnivariate( id           = mode%id,        &
                           & paired_id    = mode%paired_id, &
                           & power        = power,          &
                           & paired_power = 0               )
    endif
  endif
end function

function new_RealMonomial(coefficient,modes) result(this)
  implicit none
  
  real(dp),             intent(in) :: coefficient
  type(RealUnivariate), intent(in) :: modes(:)
  type(RealMonomial)               :: this
  
  this%coefficient = coefficient
  this%modes_      = modes(sort(modes%id))
end function

function new_RealMonomial_RealMonomialable(input) result(this)
  implicit none
  
  class(RealMonomialable), intent(in) :: input
  type(RealMonomial)                  :: this
  
  this = input%to_RealMonomial()
end function

function new_RealPolynomial(terms) result(this)
  implicit none
  
  type(RealMonomial), intent(in) :: terms(:)
  type(RealPolynomial)           :: this
  
  this%terms = terms
end function

function new_RealPolynomial_RealPolynomialable(input) result(this)
  implicit none
  
  class(RealPolynomialable), intent(in) :: input
  type(RealPolynomial)                  :: this
  
  this = input%to_RealPolynomial()
end function

! ----------------------------------------------------------------------
! Conversions between types.
! ----------------------------------------------------------------------
function to_RealMonomial_RealUnivariate(this) result(output)
  implicit none
  
  class(RealUnivariate), intent(in) :: this
  type(RealMonomial)                :: output
  
  select type(this); type is(RealUnivariate)
    output = RealMonomial( coefficient = 1.0_dp, &
                         & modes       = [this])
  class default
    call err()
  end select
end function

function to_RealPolynomial_RealUnivariate(this) result(output)
  implicit none
  
  class(RealUnivariate), intent(in) :: this
  type(RealPolynomial)              :: output
  
  select type(this); type is(RealUnivariate)
    output = RealPolynomial([this%to_RealMonomial()])
  class default
    call err()
  end select
end function

function to_RealMonomial_RealMonomial(this) result(output)
  implicit none
  
  class(RealMonomial), intent(in) :: this
  type(RealMonomial)              :: output
  
  select type(this); type is(RealMonomial)
    output = this
  class default
    call err()
  end select
end function

function to_RealPolynomial_RealMonomial(this) result(output)
  implicit none
  
  class(RealMonomial), intent(in) :: this
  type(RealPolynomial)            :: output
  
  select type(this); type is(RealMonomial)
    output = RealPolynomial([this])
  class default
    call err()
  end select
end function

function to_RealPolynomial_RealPolynomial(this) result(output)
  implicit none
  
  class(RealPolynomial), intent(in) :: this
  type(RealPolynomial)              :: output
  
  select type(this); type is(RealPolynomial)
    output = this
  class default
    call err()
  end select
end function

! ----------------------------------------------------------------------
! Operations involving types.
! ----------------------------------------------------------------------

! The number of modes in a monomial, or terms in a polynomial.
function size_RealMonomial(this) result(output)
  implicit none
  
  type(RealMonomial), intent(in) :: this
  integer                        :: output
  
  output = size(this%modes_)
end function

function size_RealPolynomial(this) result(output)
  implicit none
  
  class(RealPolynomial), intent(in) :: this
  integer                           :: output
  
  output = size(this%terms)
end function

! Getters for monomials.
impure elemental function mode_RealMonomial(this,index) result(output)
  implicit none
  
  class(RealMonomial), intent(in) :: this
  integer,             intent(in) :: index
  type(RealUnivariate)            :: output
  
  output = this%modes_(index)
end function

function modes_RealMonomial(this,indices) result(output)
  implicit none
  
  class(RealMonomial), intent(in)           :: this
  integer,             intent(in), optional :: indices(:)
  type(RealUnivariate), allocatable         :: output(:)
  
  if (present(indices)) then
    output = this%modes_(indices)
  else
    output = this%modes_
  endif
end function

impure elemental function id_RealMonomial(this,index) result(output)
  implicit none
  
  class(RealMonomial), intent(in) :: this
  integer,             intent(in) :: index
  integer                         :: output
  
  output = this%modes_(index)%id
end function

function ids_RealMonomial(this,indices) result(output)
  implicit none
  
  class(RealMonomial), intent(in)           :: this
  integer,             intent(in), optional :: indices(:)
  integer, allocatable                      :: output(:)
  
  if (present(indices)) then
    output = this%modes_(indices)%id
  else
    output = this%modes_%id
  endif
end function

impure elemental function paired_id_RealMonomial(this,index) result(output)
  implicit none
  
  class(RealMonomial), intent(in) :: this
  integer,             intent(in) :: index
  integer                         :: output
  
  output = this%modes_(index)%paired_id
end function

function paired_ids_RealMonomial(this,indices) result(output)
  implicit none
  
  class(RealMonomial), intent(in)           :: this
  integer,             intent(in), optional :: indices(:)
  integer, allocatable                      :: output(:)
  
  if (present(indices)) then
    output = this%modes_(indices)%paired_id
  else
    output = this%modes_%paired_id
  endif
end function

impure elemental function power_RealMonomial(this,index) result(output)
  implicit none
  
  class(RealMonomial), intent(in) :: this
  integer,             intent(in) :: index
  integer                         :: output
  
  output = this%modes_(index)%power
end function

function powers_RealMonomial(this,indices) result(output)
  implicit none
  
  class(RealMonomial), intent(in)           :: this
  integer,             intent(in), optional :: indices(:)
  integer, allocatable                      :: output(:)
  
  if (present(indices)) then
    output = this%modes_(indices)%power
  else
    output = this%modes_%power
  endif
end function

impure elemental function paired_power_RealMonomial(this,index) &
   & result(output)
  implicit none
  
  class(RealMonomial), intent(in) :: this
  integer,             intent(in) :: index
  integer                         :: output
  
  output = this%modes_(index)%paired_power
end function

function paired_powers_RealMonomial(this,indices) result(output)
  implicit none
  
  class(RealMonomial), intent(in)           :: this
  integer,             intent(in), optional :: indices(:)
  integer, allocatable                      :: output(:)
  
  if (present(indices)) then
    output = this%modes_(indices)%paired_power
  else
    output = this%modes_%paired_power
  endif
end function

! Simplify a monomial or polynomial.
impure elemental subroutine simplify_RealMonomial(this)
  implicit none
  
  class(RealMonomial), intent(inout) :: this
  
  integer :: i
  
  ! Combine modes with the same ID, remove modes with power=paired_power=0.
  i = 1
  do while(i<=size(this))
    if (this%modes_(i)%power<0 .or. this%modes_(i)%paired_power<0) then
      call err()
    elseif (this%modes_(i)%power==0 .and. this%modes_(i)%paired_power==0) then
      this%modes_ = [this%modes_(:i-1), this%modes_(i+1:)]
      cycle
    endif
    
    if (i>1) then
      if (this%modes_(i)%id==this%modes_(i-1)%id) then
        this%modes_(i-1)%power = this%modes_(i-1)%power + this%modes_(i)%power
        this%modes_(i-1)%paired_power = this%modes_(i-1)%paired_power &
                                    & + this%modes_(i)%paired_power
        this%modes_ = [this%modes_(:i-1), this%modes_(i+1:)]
        cycle
      endif
    endif
    
    i = i+1
  enddo
end subroutine

impure elemental subroutine simplify_RealPolynomial(this)
  implicit none
  
  class(RealPolynomial), intent(inout) :: this
  
  type(RealMonomial), allocatable :: equivalent_monomials(:)
  type(RealMonomial), allocatable :: monomials(:)
  
  real(dp) :: max_coefficient
  
  integer :: i,ialloc
  
  call this%terms%simplify()
  
  ! Add together monomials with the same powers.
  monomials = this%terms(set(this%terms, compare_real_monomials))
  do i=1,size(monomials)
    equivalent_monomials = this%terms(                          &
       & filter(this%terms,compare_real_monomials,monomials(i)) )
    monomials(i)%coefficient = sum(equivalent_monomials%coefficient)
  enddo
  
  ! Remove any monomials with coefficient less than 1e10 times smaller than
  !    the largest coefficient.
  max_coefficient = maxval(abs(monomials%coefficient))
  monomials = monomials(                                          &
     & filter(abs(monomials%coefficient)>=max_coefficient/1e10_dp) )
  
  this%terms = monomials
end subroutine

! The total power of a univariate or monomial.
impure elemental function total_power_RealUnivariate(this) result(output)
  implicit none
  
  class(RealUnivariate), intent(in) :: this
  integer                           :: output
  
  if (this%id==this%paired_id) then
    output = this%power
  else
    output = this%power + this%paired_power
  endif
end function

impure elemental function total_power_RealMonomial(this) result(output)
  implicit none
  
  class(RealMonomial), intent(in) :: this
  integer                         :: output
  
  output = sum(this%modes_%total_power())
end function

! Evaluate the contribution to the energy from
!    a univariate, monomial or polynomial at a given displacement.
impure elemental function energy_RealUnivariate(this,displacement) &
   & result(output)
  implicit none
  
  class(RealUnivariate),         intent(in) :: this
  class(RealSingleDisplacement), intent(in) :: displacement
  real(dp)                                  :: output
  
  if (displacement%id==this%id) then
    output = displacement%magnitude**this%power
  elseif (displacement%id==this%paired_id) then
    output = displacement%magnitude**this%paired_power
  else
    call print_line(CODE_ERROR//': Trying to evaluate a univariate at an &
       &incompatible displacement.')
    call err()
  endif
end function

impure elemental function energy_RealMonomial(this,displacement) &
   & result(output)
  implicit none
  
  class(RealMonomial),         intent(in) :: this
  class(RealModeDisplacement), intent(in) :: displacement
  real(dp)                                :: output
  
  integer, allocatable :: ids(:)
  
  integer :: i,j,k
  
  output = this%coefficient
  
  do i=1,size(this)
    ! Find the mode in the displacement which matches that in the monomial.
    if (this%modes_(i)%id==this%modes_(i)%paired_id) then
      if (this%modes_(i)%power/=0) then
        ids = [this%modes_(i)%id]
      endif
    else
      if (this%modes_(i)%power/=0) then
        ids = [this%modes_(i)%id]
      endif
      if (this%modes_(i)%paired_power/=0) then
        ids = [this%modes_(i)%paired_id]
      endif
    endif
    
    do j=1,size(ids)
      k = first(displacement%vectors%id==ids(j), default=0) 
      
      ! If the mode is not present in the displacement,
      !    then the displacement along that mode is zero.
      ! As such, the monomial is zero. (0^n=0 if n>0).
      if (k==0) then
        output = 0.0_dp
        return
      endif
      
      ! If the mode is present in both,
      !    evaluate the univariate at the displacement.
      output = output * this%modes_(i)%energy(displacement%vectors(k))
    enddo
  enddo
end function

impure elemental function energy_RealPolynomial(this,displacement) &
   & result(output)
  implicit none
  
  class(RealPolynomial),       intent(in) :: this
  class(RealModeDisplacement), intent(in) :: displacement
  real(dp)                                :: output
  
  output = sum(this%terms%energy(displacement))
end function

! Evaluate the contribution to the force from
!    a univariate, monomial or polynomial at a given displacement.
! -d/d{u_i} ({u_i}^n) evaluated at u_i=U is -n*U^{n-1}
impure elemental function force_RealUnivariate(this,displacement) &
   & result(output)
  implicit none
  
  class(RealUnivariate),         intent(in) :: this
  class(RealSingleDisplacement), intent(in) :: displacement
  type(RealSingleForce)                     :: output
  
  real(dp) :: force
  
  if (displacement%id==this%id) then
    if (this%power==0) then
      force = 0.0_dp
    elseif (this%power==1) then
      force = -1.0_dp
    else
      force = -this%power * displacement%magnitude**(this%power-1)
    endif
    output = RealSingleForce(id=this%id, magnitude=force)
  elseif (displacement%id==this%paired_id) then
    if (this%paired_power==0) then
      force = 0.0_dp
    elseif (this%paired_power==1) then
      force = -1.0_dp
    else
      force = -this%paired_power &
          & * displacement%magnitude**(this%paired_power-1)
    endif
    output = RealSingleForce(id=this%paired_id, magnitude=force)
  else
    call print_line(CODE_ERROR//': Trying to take the derivative of a &
       & univariate at an incompatible displacement.')
    call err()
  endif
end function

! -d/d{u_i} (c*prod_j[ {u_j}^{n_j} ]) evaluated at {u_i=U_i} is
!    -c*prod_{j/=i}[ {U_j}^{n_j} ] * n_i * {U_i}^{n_i-1}
impure elemental function force_RealMonomial(this,displacement) result(output)
  implicit none
  
  class(RealMonomial),         intent(in) :: this
  class(RealModeDisplacement), intent(in) :: displacement
  type(RealModeForce)                     :: output
  
  integer, allocatable :: id(:)
  integer, allocatable :: power(:)
  
  integer               :: displacement_id
  real(dp)              :: energy
  type(RealSingleForce) :: force
  
  integer,               allocatable :: powers(:)
  integer,               allocatable :: displacement_ids(:)
  real(dp),              allocatable :: energies(:)
  type(RealSingleForce), allocatable :: forces(:)
  
  type(RealSingleForce), allocatable :: components(:)
  
  integer :: i,j,ialloc
  
  ! Match the components of the displacement along each mode with the
  !    univariates making up the monomial.
  ! Evaluate and take the derivative of each univariate at the one-mode
  !    component of the vector.
  powers = [integer::]
  displacement_ids = [integer::]
  energies = [real(dp)::]
  forces = [RealSingleForce::]
  do i=1,size(this)
    id = [integer::]
    power = [integer::]
    if (this%modes_(i)%id==this%modes_(i)%paired_id) then
      if (this%modes_(i)%power>0) then
        id = [id, this%modes_(i)%id]
        power = [power, this%modes_(i)%power]
      endif
    else
      if (this%modes_(i)%power>0) then
        id = [id, this%modes_(i)%id]
        power = [power, this%modes_(i)%power]
      endif
      
      if (this%modes_(i)%paired_power>0) then
        id = [id, this%modes_(i)%paired_id]
        power = [power, this%modes_(i)%paired_power]
      endif
    endif
    
    do j=1,size(id)
      ! Identify the displacement corresponding to each mode.
      ! If displacement_ids(i)=0 then U_i=0.
      displacement_id = first( displacement%vectors%id==id(j), &
                             & default=0)
      
      ! Calculate {U_i}^{n_i}
      if (displacement_id==0) then
        energy = 0.0_dp
      else
        energy = this%modes_(i)%energy(            &
           & displacement%vectors(displacement_id) )
      endif
      
      ! Calculate -d/d{u_i} ({u_i}^{n_i}) evaluated at U_i.
      if (displacement_id==0) then
        if (power(j)==1) then
          force = RealSingleForce(id=id(j), magnitude=-1.0_dp )
        else
          force = RealSingleForce(id=id(j), magnitude=0.0_dp )
        endif
      else
        force = this%modes_(i)%force(              &
           & displacement%vectors(displacement_id) )
      endif
      
      powers = [powers, power(j)]
      displacement_ids = [displacement_ids, displacement_id]
      energies = [energies, energy]
      forces = [forces, force]
    enddo
  enddo
  
  ! Use the Univariate terms to calculate forces along each mode.
  if (count(displacement_ids==0)>1) then
    ! If U_i is zero for more than one i, then
    !    prod_{j/=i}[ {U_j}^{n_j} ] is always zero,
    !    so all derivatives are zero.
    components = [RealSingleForce::]
  elseif (count(displacement_ids==0)==1) then
    i = first(displacement_ids==0, default=0)
    if (powers(i)>1) then
      ! If n_i>1, then the derivative along u_i is also zero.
      components = [RealSingleForce::]
    else
      ! If n_i=1, then the derivative is simply c*prod_{j/=i}[ {U_j}^{n_j}
      components = [ this%coefficient                       &
                 & * product( energies,                     &
                 &            dim=1,                        &
                 &            mask=[(j/=i,j=1,size(this))]) &
                 & * forces(i)                              &
                 & ]
    endif
  else
    ! If no U_i are zero, then the monomial has non-zero derivatives along
    !    each of its modes.
    allocate(components(size(this)), stat=ialloc); call err(ialloc)
    do i=1,size(this)
      components(i) = this%coefficient                       &
                  & * product( energies,                     &
                  &            dim=1,                        &
                  &            mask=[(j/=i,j=1,size(this))]) &
                  & * forces(i)
    enddo
  endif
  
  ! Construct output from components along each mode.
  output = RealModeForce(components)
end function

impure elemental function force_RealPolynomial(this,displacement) &
   & result(output)
  implicit none
  
  class(RealPolynomial),       intent(in) :: this
  class(RealModeDisplacement), intent(in) :: displacement
  type(RealModeForce)                     :: output
  
  if (size(this)==0) then
    output = RealModeForce([RealSingleForce::])
  else
    output = sum(this%terms%force(displacement))
  endif
end function

! Multiplication and division by scalars.
impure elemental function multiply_RealMonomial_real(this,that) &
   & result(output)
  implicit none
  
  type(RealMonomial), intent(in) :: this
  real(dp),           intent(in) :: that
  type(RealMonomial)             :: output
  
  output = RealMonomial(this%coefficient*that, this%modes_)
end function

impure elemental function multiply_real_RealMonomial(this,that) &
   & result(output)
  implicit none
  
  real(dp),           intent(in) :: this
  type(RealMonomial), intent(in) :: that
  type(RealMonomial)             :: output
  
  output = RealMonomial(this*that%coefficient, that%modes_)
end function

impure elemental function multiply_RealPolynomial_real(this,that) &
   & result(output)
  implicit none
  
  type(RealPolynomial), intent(in) :: this
  real(dp),             intent(in) :: that
  type(RealPolynomial)             :: output
  
  output = RealPolynomial(this%terms * that)
end function

impure elemental function multiply_real_RealPolynomial(this,that) &
   & result(output)
  implicit none
  
  real(dp),             intent(in) :: this
  type(RealPolynomial), intent(in) :: that
  type(RealPolynomial)             :: output
  
  output = RealPolynomial(this * that%terms)
end function

impure elemental function divide_RealMonomial_real(this,that) result(output)
  implicit none
  
  type(RealMonomial), intent(in) :: this
  real(dp),           intent(in) :: that
  type(RealMonomial)             :: output
  
  output = RealMonomial(this%coefficient/that, this%modes_)
end function

impure elemental function divide_RealPolynomial_real(this,that) result(output)
  implicit none
  
  type(RealPolynomial), intent(in) :: this
  real(dp),             intent(in) :: that
  type(RealPolynomial)             :: output
  
  output = RealPolynomial(this%terms / that)
end function

! Multiplication between Monomial-like types.
! Uses a merge to keep ids in ascending order.
impure elemental function multiply_RealMonomialable_RealMonomialable(this, &
   & that) result(output)
  implicit none
  
  class(RealMonomialable), intent(in) :: this
  class(RealMonomialable), intent(in) :: that
  type(RealMonomial)                  :: output
  
  real(dp)                          :: coefficient
  type(RealUnivariate), allocatable :: modes(:)
  
  type(RealMonomial) :: this_monomial
  type(RealMonomial) :: that_monomial
  
  integer :: i_this,i_that,i_out,ialloc
  
  this_monomial = this%to_RealMonomial()
  that_monomial = that%to_RealMonomial()
  
  coefficient = this_monomial%coefficient * that_monomial%coefficient
  
  if (size(this_monomial)==0) then
    modes = that_monomial%modes_
  elseif (size(that_monomial)==0) then
    modes = this_monomial%modes_
  else
    i_this = 1
    i_that = 1
    i_out = 0
    allocate( modes(size(this_monomial)+size(that_monomial)), &
            & stat=ialloc); call err(ialloc)
    do while(i_this<=size(this_monomial) .or. i_that<=size(that_monomial))
      i_out = i_out + 1
      if (i_this>size(this_monomial)) then
        modes(i_out) = that_monomial%modes_(i_that)
        i_that = i_that + 1
      elseif (i_that>size(that_monomial)) then
        modes(i_out) = this_monomial%modes_(i_this)
        i_this = i_this + 1
      elseif ( this_monomial%modes_(i_this)%id == &
             & that_monomial%modes_(i_that)%id    ) then
        modes(i_out) = RealUnivariate(                                 &
           & id           = this_monomial%modes_(i_this)%id,           &
           & paired_id    = this_monomial%modes_(i_this)%paired_id,    &
           & power        = this_monomial%modes_(i_this)%power         &
           &              + that_monomial%modes_(i_that)%power,        &
           & paired_power = this_monomial%modes_(i_this)%paired_power  &
           &              + that_monomial%modes_(i_that)%paired_power  )
        i_this = i_this + 1
        i_that = i_that + 1
      elseif ( this_monomial%modes_(i_this)%id < &
             & that_monomial%modes_(i_that)%id   ) then
        modes(i_out) = this_monomial%modes_(i_this)
        i_this = i_this + 1
      elseif ( this_monomial%modes_(i_this)%id > &
             & that_monomial%modes_(i_that)%id   ) then
        modes(i_out) = that_monomial%modes_(i_that)
        i_that = i_that + 1
      else
        call err()
      endif
    enddo
    modes = modes(:i_out)
  endif
  
  output = RealMonomial(coefficient, modes)
end function

! Multiplication between polynomials and monomial-like types.
impure elemental function multiply_RealPolynomial_RealMonomialable( &
   & this,that) result(output)
  implicit none
  
  type(RealPolynomial),    intent(in) :: this
  class(RealMonomialable), intent(in) :: that
  type(RealPolynomial)                :: output
  
  output = RealPolynomial(this%terms * that)
end function

impure elemental function multiply_RealMonomialable_RealPolynomial( &
   & this,that) result(output)
  implicit none
  
  class(RealMonomialable), intent(in) :: this
  type(RealPolynomial),    intent(in) :: that
  type(RealPolynomial)                :: output
  
  output = RealPolynomial(this * that%terms)
end function

! Addition between polynomials and polynomial-like types.
impure elemental function add_RealPolynomialable_RealPolynomialable(this, &
   & that) result(output)
  implicit none
  
  class(RealPolynomialable), intent(in) :: this
  class(RealPolynomialable), intent(in) :: that
  type(RealPolynomial)                  :: output
  
  type(RealPolynomial) :: this_polynomial
  type(RealPolynomial) :: that_polynomial
  
  integer :: no_terms
  
  integer :: i,j,ialloc
  
  this_polynomial = this%to_RealPolynomial()
  that_polynomial = that%to_RealPolynomial()
  
  allocate( output%terms(size(this_polynomial)+size(that_polynomial)), &
          & stat=ialloc); call err(ialloc)
  output%terms(:size(this_polynomial)) = this_polynomial%terms
  no_terms = size(this_polynomial)
  do i=1,size(that_polynomial)
    j = first_equivalent( this_polynomial%terms,    &
                        & that_polynomial%terms(i), &
                        & compare_real_monomials,   &
                        & default=0)
    if (j==0) then
      no_terms = no_terms + 1
      output%terms(no_terms) = that_polynomial%terms(i)
    else
      output%terms(j)%coefficient = output%terms(j)%coefficient &
                                & + that_polynomial%terms(i)%coefficient
    endif
  enddo
  output%terms = output%terms(:no_terms)
end function

! The negative of a polynomial or polynomial-like type.
impure elemental function negative_RealPolynomialable(this) result(output)
  implicit none
  
  class(RealPolynomialable), intent(in) :: this
  type(RealPolynomial)                  :: output
  
  output = this%to_RealPolynomial()
  output%terms%coefficient = - output%terms%coefficient
end function

! Subtraction between polynomials and polynomial-like types.
impure elemental function subtract_RealPolynomialable_RealPolynomialable( &
   & this,that) result(output)
  implicit none
  
  class(RealPolynomialable), intent(in) :: this
  class(RealPolynomialable), intent(in) :: that
  type(RealPolynomial)                  :: output
  
  output = this + (-that)
end function

! Sum polynomial-like types.
function sum_RealPolynomialables(input) result(output)
  implicit none
  
  class(RealPolynomialable), intent(in) :: input(:)
  type(RealPolynomial)                  :: output
  
  integer :: i
  
  if (size(input)==0) then
    output = RealPolynomial([RealMonomial::])
  else
    output = input(1)%to_RealPolynomial()
    do i=2,size(input)
      output = output + input(i)
    enddo
  endif
end function

! ----------------------------------------------------------------------
! Select modes, displacements or forces corresponding to
!    a given univariate or univariates.
! ----------------------------------------------------------------------
function select_mode_RealUnivariate(univariate,modes) result(output)
  implicit none
  
  type(RealUnivariate), intent(in) :: univariate
  type(RealMode),       intent(in) :: modes(:)
  type(RealMode)                   :: output
  
  output = modes(first(modes%id==univariate%id))
end function

function select_modes_RealUnivariates(univariates,modes) result(output)
  implicit none
  
  type(RealUnivariate), intent(in) :: univariates(:)
  type(RealMode),       intent(in) :: modes(:)
  type(RealMode), allocatable      :: output(:)
  
  integer :: i,ialloc
  
  allocate(output(size(univariates)), stat=ialloc); call err(ialloc)
  do i=1,size(univariates)
    output(i) = select_mode(univariates(i), modes)
  enddo
end function

function select_displacement_RealUnivariate(univariate,displacements) &
   & result(output)
  implicit none
  
  type(RealUnivariate),         intent(in) :: univariate
  type(RealSingleDisplacement), intent(in) :: displacements(:)
  type(RealSingleDisplacement)             :: output
  
  output = displacements(first(displacements%id==univariate%id))
end function

function select_displacements_RealUnivariates(univariates,displacements) &
   & result(output)
  implicit none
  
  type(RealUnivariate),         intent(in)  :: univariates(:)
  type(RealSingleDisplacement), intent(in)  :: displacements(:)
  type(RealSingleDisplacement), allocatable :: output(:)
  
  integer :: i,ialloc
  
  allocate(output(size(univariates)), stat=ialloc); call err(ialloc)
  do i=1,size(univariates)
    output(i) = select_displacement(univariates(i), displacements)
  enddo
end function

function select_force_RealUnivariate(univariate,forces) &
   & result(output)
  implicit none
  
  type(RealUnivariate),  intent(in) :: univariate
  type(RealSingleForce), intent(in) :: forces(:)
  type(RealSingleForce)             :: output
  
  output = forces(first(forces%id==univariate%id))
end function

function select_forces_RealUnivariates(univariates,forces) &
   & result(output)
  implicit none
  
  type(RealUnivariate),  intent(in)  :: univariates(:)
  type(RealSingleForce), intent(in)  :: forces(:)
  type(RealSingleForce), allocatable :: output(:)
  
  integer :: i,ialloc
  
  allocate(output(size(univariates)), stat=ialloc); call err(ialloc)
  do i=1,size(univariates)
    output(i) = select_force(univariates(i), forces)
  enddo
end function

! ----------------------------------------------------------------------
! Compares two monomials for equality up to coefficient.
! ----------------------------------------------------------------------
function compare_real_monomials(this,that) result(output)
  implicit none
  
  class(*), intent(in) :: this
  class(*), intent(in) :: that
  logical              :: output
  
  select type(this); type is(RealMonomial)
    select type(that); type is(RealMonomial)
      if (size(this)/=size(that)) then
        output = .false.
      else
        output = all( this%modes_%id           == that%modes_%id    .and.  &
                    & this%modes_%power        == that%modes_%power .and.  &
                    & this%modes_%paired_power == that%modes_%paired_power )
      endif
    end select
  end select
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_RealUnivariate(this,input)
  implicit none
  
  class(RealUnivariate), intent(out) :: this
  type(String),          intent(in)  :: input
  
  type(String), allocatable :: line(:)
  integer                   :: id
  integer                   :: paired_id
  integer                   :: power
  integer                   :: paired_power
  
  select type(this); type is(RealUnivariate)
    ! If id=paired_id=5 and power=paired_power=3 then:
    !    input = '(u5^3)'.
    ! If id=5, paired_id=7, power=3 and paired_power=4 then:
    !    input = '(u5^3*u7^4)'.
    
    ! Strip off brackets, and split into mode and paired mode.
    line = split_line( slice(input,2,len(input)-1), &
                     & delimiter='*')
    
    if (size(line)==1) then
      ! line = [ 'u5^3' ]
      ! ID = paired ID.
      
      ! Split into ID and power.
      line = split_line(line(1), delimiter='^')
      id = int(slice(line(1),2,len(line(1))))
      power = int(line(2))
      paired_id = id
      paired_power = power
    elseif (size(line)==2) then
      ! line = [ 'u5^3', 'u7^4' ]
      ! ID /= paired ID.
      
      ! Split into ID, power, paired ID, paired power.
      line = [ split_line(line(1), delimiter='^'), &
             & split_line(line(2), delimiter='^')  ]
      id = int(slice(line(1),2,len(line(1))))
      power = int(line(2))
      paired_id = int(slice(line(3),2,len(line(3))))
      paired_power = int(line(4))
    else
      call print_line(ERROR//': unable to parse RealUnivariate.')
      call err()
    endif
    
    this = RealUnivariate( id           = id,          &
                         & paired_id    = paired_id,   &
                         & power        = power,       &
                         & paired_power = paired_power )
  end select
end subroutine

function write_RealUnivariate(this) result(output)
  implicit none
  
  class(RealUnivariate), intent(in) :: this
  type(String)                      :: output
  
  select type(this); type is(RealUnivariate)
    if (this%id==this%paired_id) then
      output = '(u'//this%id//'^'//this%power//')'
    else
      output = '(u'//this%id//'^'//this%power// &
             & '*u'//this%paired_id//'^'//this%paired_power//')'
    endif
  end select
end function

impure elemental function new_RealUnivariate_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(RealUnivariate)     :: this
  
  call this%read(input)
end function

subroutine read_RealMonomial(this,input)
  implicit none
  
  class(RealMonomial), intent(out) :: this
  type(String),        intent(in)  :: input
  
  type(String),         allocatable :: line(:)
  real(dp)                          :: coefficient
  type(RealUnivariate), allocatable :: modes(:)
  
  integer :: i,ialloc
  
  select type(this); type is(RealMonomial)
    ! Splitting the input by '*' separates the coefficient and the modes,
    !    but also splits some modes in two.
    line = split_line(input,delimiter='*')
    
    coefficient = dble(line(1))
    
    modes = [RealUnivariate::]
    i = 2
    do while (i<=size(line))
      ! Check if line(i) ends in a bracket.
      if (slice(line(i),len(line(i)),len(line(i)))==')') then
        ! line(i) is a mode on its own.
        modes = [modes, RealUnivariate(line(i))]
        i = i+1
      else
        ! line(i) and line(i+1) together make a mode.
        modes = [modes, RealUnivariate(line(i)//'*'//line(i+1))]
        i = i+2
      endif
    enddo
    
    this = RealMonomial(coefficient, modes)
  end select
end subroutine

function write_RealMonomial(this) result(output)
  implicit none
  
  class(RealMonomial), intent(in) :: this
  type(String)                    :: output
  
  select type(this); type is(RealMonomial)
    if (size(this%modes_)>0) then
      output = this%coefficient//'*'//join(this%modes_, delimiter='*')
    else
      output = str(this%coefficient)
    endif
  end select
end function

impure elemental function new_RealMonomial_String(input) result(this)
  implicit none
  
  type(String),   intent(in) :: input
  type(RealMonomial)         :: this
  
  call this%read(input)
end function

subroutine read_RealPolynomial(this,input)
  implicit none
  
  class(RealPolynomial), intent(out) :: this
  type(String),          intent(in)  :: input
  
  type(String),       allocatable :: terms(:)
  type(String)                    :: plus
  type(RealMonomial), allocatable :: monomials(:)
  
  plus = '+'
  select type(this); type is(RealPolynomial)
    terms = split_line(input)
    terms = terms(filter(terms/=plus))
    
    monomials = RealMonomial(terms)
    
    this = RealPolynomial(monomials)
  end select
end subroutine

function write_RealPolynomial(this) result(output)
  implicit none
  
  class(RealPolynomial), intent(in) :: this
  type(String)                      :: output
  
  select type(this); type is(RealPolynomial)
    output = join(this%terms, delimiter=' + ')
  end select
end function

impure elemental function new_RealPolynomial_String(input) result(this)
  implicit none
  
  type(String),   intent(in) :: input
  type(RealPolynomial)       :: this
  
  call this%read(input)
end function
end module
