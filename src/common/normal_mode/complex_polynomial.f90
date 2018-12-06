! ======================================================================
! The building blocks of basis functions in complex co-ordinates.
! ======================================================================
module complex_polynomial_module
  use utils_module
  use structure_module
  
  use complex_mode_module
  use complex_single_mode_displacement_module
  use complex_single_mode_force_module
  use complex_mode_displacement_module
  use complex_mode_force_module
  implicit none
  
  private
  
  public :: ComplexUnivariate
  public :: ComplexMonomial
  public :: ComplexPolynomial
  public :: size
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: operator(-)
  public :: sum
  public :: conjg
  public :: select_modes
  public :: select_displacements
  public :: select_forces
  public :: compare_complex_monomials
  
  ! --------------------------------------------------
  ! Types, and conversions between types.
  ! --------------------------------------------------
  
  ! It is desirable to be able to convert:
  !    univariates -> monomials -> polynomials
  ! This is acheived by extending from the classes ComplexPolynomialable and
  !    ComplexMonomialable, representing types which can be converted to
  !    ComplexPolynomial and ComplexMonomial respectively.
  
  type, abstract, extends(Stringable) :: ComplexPolynomialable
  contains
    procedure(to_ComplexPolynomial_ComplexPolynomialable), deferred, public &
       & :: to_ComplexPolynomial
    
    procedure(wavevector_ComplexPolynomialable), deferred, public &
       & :: wavevector
  end type
  
  type, abstract, extends(ComplexPolynomialable) :: ComplexMonomialable
  contains
    procedure(to_ComplexMonomial_ComplexMonomialable), deferred, public :: &
       & to_ComplexMonomial
  end type
  
  type, extends(ComplexMonomialable) :: ComplexUnivariate
    integer :: id
    integer :: paired_id
    integer :: power
    integer :: paired_power
  contains
    procedure, public :: to_ComplexMonomial   => &
       & to_ComplexMonomial_ComplexUnivariate
    procedure, public :: to_ComplexPolynomial => &
       & to_ComplexPolynomial_ComplexUnivariate
    
    procedure, public :: total_power => total_power_ComplexUnivariate
    
    procedure, public :: wavevector => &
                       & wavevector_ComplexUnivariate
    
    procedure, public :: energy => energy_ComplexUnivariate
    procedure, public :: force  => force_ComplexUnivariate
    
    procedure, public :: read  => read_ComplexUnivariate
    procedure, public :: write => write_ComplexUnivariate
  end type
  
  interface ComplexUnivariate
    module procedure new_ComplexUnivariate
    module procedure new_ComplexUnivariate_ComplexMode
    module procedure new_ComplexUnivariate_String
  end interface
  
  type, extends(ComplexMonomialable) :: ComplexMonomial
    complex(dp)                                   :: coefficient
    ! Modes is private so that it can be guaranteed to be sorted.
    type(ComplexUnivariate), allocatable, private :: modes_(:)
  contains
    procedure, public :: to_ComplexMonomial   => &
       & to_ComplexMonomial_ComplexMonomial
    procedure, public :: to_ComplexPolynomial => &
       & to_ComplexPolynomial_ComplexMonomial
    
    procedure, public :: mode => mode_ComplexMonomial
    procedure, public :: modes => modes_ComplexMonomial
    procedure, public :: id  => id_ComplexMonomial
    procedure, public :: ids => ids_ComplexMonomial
    procedure, public :: paired_id  => paired_id_ComplexMonomial
    procedure, public :: paired_ids => paired_ids_ComplexMonomial
    procedure, public :: power  => power_ComplexMonomial
    procedure, public :: powers => powers_ComplexMonomial
    procedure, public :: paired_power  => paired_power_ComplexMonomial
    procedure, public :: paired_powers => paired_powers_ComplexMonomial
    
    procedure, public :: simplify => simplify_ComplexMonomial
    
    procedure, public :: total_power => total_power_ComplexMonomial
    
    procedure, public :: wavevector => wavevector_ComplexMonomial
    
    procedure, public :: energy => energy_ComplexMonomial
    procedure, public :: force  => force_ComplexMonomial
    
    procedure, public :: read  => read_ComplexMonomial
    procedure, public :: write => write_ComplexMonomial
  end type
  
  interface ComplexMonomial
    module procedure new_ComplexMonomial
    module procedure new_ComplexMonomial_ComplexMonomialable
    module procedure new_ComplexMonomial_String
  end interface
  
  type, extends(ComplexPolynomialable) :: ComplexPolynomial
    type(ComplexMonomial), allocatable :: terms(:)
  contains
    procedure, public :: to_ComplexPolynomial => &
       & to_ComplexPolynomial_ComplexPolynomial
    
    procedure, public :: simplify => simplify_ComplexPolynomial
    
    procedure, public :: wavevector => &
                       & wavevector_ComplexPolynomial
    
    procedure, public :: energy => energy_ComplexPolynomial
    procedure, public :: force  => force_ComplexPolynomial
    
    procedure, public :: read  => read_ComplexPolynomial
    procedure, public :: write => write_ComplexPolynomial
  end type
  
  interface ComplexPolynomial
    module procedure new_ComplexPolynomial
    module procedure new_ComplexPolynomial_ComplexPolynomialable
    module procedure new_ComplexPolynomial_String
  end interface
  
  abstract interface
    function to_ComplexPolynomial_ComplexPolynomialable(this) result(output)
      import ComplexPolynomial
      import ComplexPolynomialable
      implicit none
      
      class(ComplexPolynomialable), intent(in) :: this
      type(ComplexPolynomial)                  :: output
    end function
    
    function wavevector_ComplexPolynomialable(this,modes,qpoints) &
       & result(output)
      import ComplexPolynomialable
      import ComplexMode
      import QpointData
      implicit none
      
      class(ComplexPolynomialable), intent(in) :: this
      type(ComplexMode),            intent(in) :: modes(:)
      type(QpointData),             intent(in) :: qpoints(:)
      type(QpointData)                         :: output
    end function
    
    function to_ComplexMonomial_ComplexMonomialable(this) result(output)
      import ComplexMonomial
      import ComplexMonomialable
      implicit none
      
      class(ComplexMonomialable), intent(in) :: this
      type(ComplexMonomial)                  :: output
    end function
  end interface
  
  ! --------------------------------------------------
  ! Operations involving types.
  ! --------------------------------------------------
  
  interface size
    module procedure size_ComplexMonomial
    module procedure size_ComplexPolynomial
  end interface
  
  interface conjg
    module procedure conjg_ComplexUnivariate
    module procedure conjg_ComplexMonomial
    module procedure conjg_ComplexPolynomial
  end interface
  
  interface operator(*)
    module procedure multiply_ComplexMonomial_real
    module procedure multiply_real_ComplexMonomial
    module procedure multiply_ComplexMonomial_complex
    module procedure multiply_complex_ComplexMonomial
    
    module procedure multiply_ComplexPolynomial_real
    module procedure multiply_real_ComplexPolynomial
    module procedure multiply_ComplexPolynomial_complex
    module procedure multiply_complex_ComplexPolynomial
    
    module procedure multiply_ComplexMonomialable_ComplexMonomialable
    
    module procedure multiply_ComplexPolynomial_ComplexMonomialable
    module procedure multiply_ComplexMonomialable_ComplexPolynomial
  end interface
  
  interface operator(/)
    module procedure divide_ComplexMonomial_real
    module procedure divide_ComplexMonomial_complex
    
    module procedure divide_ComplexPolynomial_real
    module procedure divide_ComplexPolynomial_complex
  end interface
  
  interface operator(+)
    module procedure add_ComplexPolynomialable_ComplexPolynomialable
  end interface
  
  interface operator(-)
    module procedure negative_ComplexPolynomialable
    
    module procedure subtract_ComplexPolynomialable_ComplexPolynomialable
  end interface
  
  interface sum
    module procedure sum_ComplexPolynomialables
  end interface
  
  interface select_modes
    module procedure select_modes_ComplexUnivariate
    module procedure select_modes_ComplexUnivariates
  end interface
  
  interface select_displacements
    module procedure select_displacements_ComplexUnivariate
    module procedure select_displacements_ComplexUnivariates
  end interface
  
  interface select_forces
    module procedure select_forces_ComplexUnivariate
    module procedure select_forces_ComplexUnivariates
  end interface
contains

! ----------------------------------------------------------------------
! Constructors.
! ----------------------------------------------------------------------
impure elemental function new_ComplexUnivariate(id,paired_id,power, &
   & paired_power) result(this)
  implicit none
  
  integer, intent(in)     :: id
  integer, intent(in)     :: paired_id
  integer, intent(in)     :: power
  integer, intent(in)     :: paired_power
  type(ComplexUnivariate) :: this
  
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

function new_ComplexUnivariate_ComplexMode(mode,power,paired_power) &
   & result(this)
  implicit none
  
  type(ComplexMode), intent(in)           :: mode
  integer,           intent(in)           :: power
  integer,           intent(in), optional :: paired_power
  type(ComplexUnivariate)                 :: this
  
  if (present(paired_power)) then
    if (mode%id==mode%paired_id .and. power/=paired_power) then
      call print_line(ERROR//': Mode is its own pair, but power does not &
         & match paired_power.')
      call err()
    endif
    this = ComplexUnivariate( id           = mode%id,        &
                            & paired_id    = mode%paired_id, &
                            & power        = power,          &
                            & paired_power = paired_power    )
  else
    if (mode%id==mode%paired_id) then
      this = ComplexUnivariate( id           = mode%id, &
                              & paired_id    = mode%id, &
                              & power        = power,   &
                              & paired_power = power    )
    else
      this = ComplexUnivariate( id           = mode%id,        &
                              & paired_id    = mode%paired_id, &
                              & power        = power,          &
                              & paired_power = 0               )
    endif
  endif
end function

function new_ComplexMonomial(coefficient,modes) result(this)
  implicit none
  
  complex(dp),             intent(in) :: coefficient
  type(ComplexUnivariate), intent(in) :: modes(:)
  type(ComplexMonomial)               :: this
  
  this%coefficient = coefficient
  this%modes_      = modes(sort(modes%id))
end function

function new_ComplexMonomial_ComplexMonomialable(input) result(this)
  implicit none
  
  class(ComplexMonomialable), intent(in) :: input
  type(ComplexMonomial)                  :: this
  
  this = input%to_ComplexMonomial()
end function

function new_ComplexPolynomial(terms) result(this)
  implicit none
  
  type(ComplexMonomial), intent(in) :: terms(:)
  type(ComplexPolynomial)           :: this
  
  this%terms = terms
end function

function new_ComplexPolynomial_ComplexPolynomialable(input) result(this)
  implicit none
  
  class(ComplexPolynomialable), intent(in) :: input
  type(ComplexPolynomial)                  :: this
  
  this = input%to_ComplexPolynomial()
end function

! ----------------------------------------------------------------------
! Conversions between types.
! ----------------------------------------------------------------------
function to_ComplexMonomial_ComplexUnivariate(this) result(output)
  implicit none
  
  class(ComplexUnivariate), intent(in) :: this
  type(ComplexMonomial)                :: output
  
  select type(this); type is(ComplexUnivariate)
    output = ComplexMonomial( coefficient = cmplx(1.0_dp,0.0_dp,dp), &
                            & modes       = [this])
  class default
    call err()
  end select
end function

function to_ComplexPolynomial_ComplexUnivariate(this) result(output)
  implicit none
  
  class(ComplexUnivariate), intent(in) :: this
  type(ComplexPolynomial)              :: output
  
  select type(this); type is(ComplexUnivariate)
    output = ComplexPolynomial([this%to_ComplexMonomial()])
  class default
    call err()
  end select
end function

function to_ComplexMonomial_ComplexMonomial(this) result(output)
  implicit none
  
  class(ComplexMonomial), intent(in) :: this
  type(ComplexMonomial)              :: output
  
  select type(this); type is(ComplexMonomial)
    output = this
  class default
    call err()
  end select
end function

function to_ComplexPolynomial_ComplexMonomial(this) result(output)
  implicit none
  
  class(ComplexMonomial), intent(in) :: this
  type(ComplexPolynomial)            :: output
  
  select type(this); type is(ComplexMonomial)
    output = ComplexPolynomial([this])
  class default
    call err()
  end select
end function

function to_ComplexPolynomial_ComplexPolynomial(this) result(output)
  implicit none
  
  class(ComplexPolynomial), intent(in) :: this
  type(ComplexPolynomial)              :: output
  
  select type(this); type is(ComplexPolynomial)
    output = this
  class default
    call err()
  end select
end function

! ----------------------------------------------------------------------
! Operations involving types.
! ----------------------------------------------------------------------

! The number of modes in a monomial, or terms in a polynomial.
function size_ComplexMonomial(this) result(output)
  implicit none
  
  type(ComplexMonomial), intent(in) :: this
  integer                           :: output
  
  output = size(this%modes_)
end function

function size_ComplexPolynomial(this) result(output)
  implicit none
  
  class(ComplexPolynomial), intent(in) :: this
  integer                              :: output
  
  output = size(this%terms)
end function

! Getters for monomials.
impure elemental function mode_ComplexMonomial(this,index) result(output)
  implicit none
  
  class(ComplexMonomial), intent(in) :: this
  integer,                intent(in) :: index
  type(ComplexUnivariate)            :: output
  
  output = this%modes_(index)
end function

function modes_ComplexMonomial(this,indices) result(output)
  implicit none
  
  class(ComplexMonomial), intent(in)           :: this
  integer,                intent(in), optional :: indices(:)
  type(ComplexUnivariate), allocatable         :: output(:)
  
  if (present(indices)) then
    output = this%modes_(indices)
  else
    output = this%modes_
  endif
end function

impure elemental function id_ComplexMonomial(this,index) result(output)
  implicit none
  
  class(ComplexMonomial), intent(in) :: this
  integer,                intent(in) :: index
  integer                            :: output
  
  output = this%modes_(index)%id
end function

function ids_ComplexMonomial(this,indices) result(output)
  implicit none
  
  class(ComplexMonomial), intent(in)           :: this
  integer,                intent(in), optional :: indices(:)
  integer, allocatable                         :: output(:)
  
  if (present(indices)) then
    output = this%modes_(indices)%id
  else
    output = this%modes_%id
  endif
end function

impure elemental function paired_id_ComplexMonomial(this,index) result(output)
  implicit none
  
  class(ComplexMonomial), intent(in) :: this
  integer,                intent(in) :: index
  integer                            :: output
  
  output = this%modes_(index)%paired_id
end function

function paired_ids_ComplexMonomial(this,indices) result(output)
  implicit none
  
  class(ComplexMonomial), intent(in)           :: this
  integer,                intent(in), optional :: indices(:)
  integer, allocatable                         :: output(:)
  
  if (present(indices)) then
    output = this%modes_(indices)%paired_id
  else
    output = this%modes_%paired_id
  endif
end function

impure elemental function power_ComplexMonomial(this,index) result(output)
  implicit none
  
  class(ComplexMonomial), intent(in) :: this
  integer,                intent(in) :: index
  integer                            :: output
  
  output = this%modes_(index)%power
end function

function powers_ComplexMonomial(this,indices) result(output)
  implicit none
  
  class(ComplexMonomial), intent(in)           :: this
  integer,                intent(in), optional :: indices(:)
  integer, allocatable                         :: output(:)
  
  if (present(indices)) then
    output = this%modes_(indices)%power
  else
    output = this%modes_%power
  endif
end function

impure elemental function paired_power_ComplexMonomial(this,index) &
   & result(output)
  implicit none
  
  class(ComplexMonomial), intent(in) :: this
  integer,                intent(in) :: index
  integer                            :: output
  
  output = this%modes_(index)%paired_power
end function

function paired_powers_ComplexMonomial(this,indices) result(output)
  implicit none
  
  class(ComplexMonomial), intent(in)           :: this
  integer,                intent(in), optional :: indices(:)
  integer, allocatable                         :: output(:)
  
  if (present(indices)) then
    output = this%modes_(indices)%paired_power
  else
    output = this%modes_%paired_power
  endif
end function

! Simplify a monomial or polynomial.
impure elemental subroutine simplify_ComplexMonomial(this)
  implicit none
  
  class(ComplexMonomial), intent(inout) :: this
  
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

impure elemental subroutine simplify_ComplexPolynomial(this)
  implicit none
  
  class(ComplexPolynomial), intent(inout) :: this
  
  type(ComplexMonomial), allocatable :: equivalent_monomials(:)
  type(ComplexMonomial), allocatable :: monomials(:)
  
  real(dp) :: max_coefficient
  
  integer :: i,ialloc
  
  call this%terms%simplify()
  
  ! Add together monomials with the same powers.
  monomials = this%terms(set(this%terms, compare_complex_monomials))
  do i=1,size(monomials)
    equivalent_monomials = this%terms(                          &
       & filter(this%terms,compare_complex_monomials,monomials(i)) )
    monomials(i)%coefficient = sum(equivalent_monomials%coefficient)
  enddo
  
  ! Remove any monomials with coefficient less than 1e10 times smaller than
  !    the largest coefficient.
  max_coefficient = maxval(abs(monomials%coefficient))
  monomials = monomials(                                          &
     & filter(abs(monomials%coefficient)>=max_coefficient/1e10_dp) )
  
  this%terms = monomials
end subroutine

! Find the conjugate of a univariate or monomial.
impure elemental function conjg_ComplexUnivariate(this) result(output)
  implicit none
  
  type(ComplexUnivariate), intent(in) :: this
  type(ComplexUnivariate)             :: output
  
  output = ComplexUnivariate( id           = this%paired_id,   &
                            & paired_id    = this%id,          &
                            & power        = this%power,       &
                            & paired_power = this%paired_power )
end function

impure elemental function conjg_ComplexMonomial(this) result(output)
  implicit none
  
  type(ComplexMonomial), intent(in) :: this
  type(ComplexMonomial)             :: output
  
  output = ComplexMonomial( coefficient = conjg(this%coefficient), &
                          & modes       = conjg(this%modes_)       )
end function

impure elemental function conjg_ComplexPolynomial(this) result(output)
  implicit none
  
  type(ComplexPolynomial), intent(in) :: this
  type(ComplexPolynomial)             :: output
  
  output = ComplexPolynomial(conjg(this%terms))
end function

! The total power of a univariate or monomial.
impure elemental function total_power_ComplexUnivariate(this) result(output)
  implicit none
  
  class(ComplexUnivariate), intent(in) :: this
  integer                              :: output
  
  if (this%id==this%paired_id) then
    output = this%power
  else
    output = this%power + this%paired_power
  endif
end function

impure elemental function total_power_ComplexMonomial(this) result(output)
  implicit none
  
  class(ComplexMonomial), intent(in) :: this
  integer                            :: output
  
  output = sum(this%modes_%total_power())
end function

! Returns the Bloch wavevector of a univariate, monomial or polynomial.
function wavevector_ComplexUnivariate(this,modes,qpoints) result(output)
  implicit none
  
  class(ComplexUnivariate), intent(in) :: this
  type(ComplexMode),        intent(in) :: modes(:)
  type(QpointData),         intent(in) :: qpoints(:)
  type(QpointData)                     :: output
  
  type(ComplexMode) :: mode
  type(QpointData)  :: qpoint
  
  mode = modes(first(modes%id==this%id))
  qpoint = qpoints(first(qpoints%id==mode%qpoint_id))
  if (this%id==this%paired_id) then
    qpoint%qpoint = qpoint%qpoint * this%power
  else
    qpoint%qpoint = qpoint%qpoint * (this%power-this%paired_power)
  endif
  
  output = qpoints(first(qpoints==qpoint))
end function

function wavevector_ComplexMonomial(this,modes,qpoints) result(output)
  implicit none
  
  class(ComplexMonomial), intent(in) :: this
  type(ComplexMode),      intent(in) :: modes(:)
  type(QpointData),       intent(in) :: qpoints(:)
  type(QpointData)                   :: output
  
  type(QpointData) :: wavevector
  type(QpointData) :: mode_wavevector
  
  integer :: i
  
  wavevector%qpoint = fracvec(zeroes(3))
  do i=1,size(this%modes_)
    mode_wavevector = this%modes_(i)%wavevector(modes,qpoints)
    wavevector%qpoint = wavevector%qpoint + mode_wavevector%qpoint
  enddo
  
  output = qpoints(first(qpoints==wavevector))
end function

function wavevector_ComplexPolynomial(this,modes,qpoints) result(output)
  implicit none
  
  class(ComplexPolynomial), intent(in) :: this
  type(ComplexMode),        intent(in) :: modes(:)
  type(QpointData),         intent(in) :: qpoints(:)
  type(QpointData)                     :: output
  
  integer :: i
  
  if (size(this)==0) then
    output%qpoint = fracvec(zeroes(3))
    output = qpoints(first(qpoints==output))
  else
    output = this%terms(1)%wavevector(modes,qpoints)
    do i=2,size(this)
      if (this%terms(i)%wavevector(modes,qpoints)/=output) then
        call print_line(ERROR//': Complex polynomial has inconsistent &
           & wavevector.')
        call err()
      endif
    enddo
  endif
end function

! Evaluate the contribution to the energy from
!    a univariate, monomial or polynomial at a given displacement.
impure elemental function energy_ComplexUnivariate(this,displacement) &
   & result(output)
  implicit none
  
  class(ComplexUnivariate),         intent(in) :: this
  class(ComplexSingleDisplacement), intent(in) :: displacement
  complex(dp)                                  :: output
  
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

impure elemental function energy_ComplexMonomial(this,displacement) &
   & result(output)
  implicit none
  
  class(ComplexMonomial),         intent(in) :: this
  class(ComplexModeDisplacement), intent(in) :: displacement
  complex(dp)                                :: output
  
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

impure elemental function energy_ComplexPolynomial(this,displacement) &
   & result(output)
  implicit none
  
  class(ComplexPolynomial),       intent(in) :: this
  class(ComplexModeDisplacement), intent(in) :: displacement
  complex(dp)                                :: output
  
  output = sum(this%terms%energy(displacement))
end function

! Evaluate the contribution to the force from
!    a univariate, monomial or polynomial at a given displacement.
! -d/d{u_i} ({u_i}^n) evaluated at u_i=U is -n*U^{n-1}
impure elemental function force_ComplexUnivariate(this,displacement) &
   & result(output)
  implicit none
  
  class(ComplexUnivariate),         intent(in) :: this
  class(ComplexSingleDisplacement), intent(in) :: displacement
  type(ComplexSingleForce)                     :: output
  
  complex(dp) :: force
  
  if (displacement%id==this%id) then
    if (this%power==0) then
      force = 0.0_dp
    elseif (this%power==1) then
      force = -1.0_dp
    else
      force = -this%power * displacement%magnitude**(this%power-1)
    endif
    output = ComplexSingleForce(id=this%id, magnitude=force)
  elseif (displacement%id==this%paired_id) then
    if (this%paired_power==0) then
      force = 0.0_dp
    elseif (this%paired_power==1) then
      force = -1.0_dp
    else
      force = -this%paired_power &
          & * displacement%magnitude**(this%paired_power-1)
    endif
    output = ComplexSingleForce(id=this%paired_id, magnitude=force)
  else
    call print_line(CODE_ERROR//': Trying to take the derivative of a &
       & univariate at an incompatible displacement.')
    call err()
  endif
end function

! -d/d{u_i} (c*prod_j[ {u_j}^{n_j} ]) evaluated at {u_i=U_i} is
!    -c*prod_{j/=i}[ {U_j}^{n_j} ] * n_i * {U_i}^{n_i-1}
impure elemental function force_ComplexMonomial(this,displacement) &
   & result(output)
  implicit none
  
  class(ComplexMonomial),         intent(in) :: this
  class(ComplexModeDisplacement), intent(in) :: displacement
  type(ComplexModeForce)                     :: output
  
  integer, allocatable :: id(:)
  integer, allocatable :: power(:)
  
  integer                  :: displacement_id
  complex(dp)              :: energy
  type(ComplexSingleForce) :: force
  
  integer,                  allocatable :: powers(:)
  integer,                  allocatable :: displacement_ids(:)
  complex(dp),              allocatable :: energies(:)
  type(ComplexSingleForce), allocatable :: forces(:)
  
  type(ComplexSingleForce), allocatable :: components(:)
  
  integer :: i,j,ialloc
  
  ! Match the components of the displacement along each mode with the
  !    univariates making up the monomial.
  ! Evaluate and take the derivative of each univariate at the one-mode
  !    component of the vector.
  powers = [integer::]
  displacement_ids = [integer::]
  energies = [complex(dp)::]
  forces = [ComplexSingleForce::]
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
      ! Identify the vector corresponding to each mode.
      ! If vector_ids(i)=0 then U_i=0.
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
          force = ComplexSingleForce(             &
             & id=id(j),                          &
             & magnitude=cmplx(-1.0_dp,0.0_dp,dp) )
        else
          force = ComplexSingleForce(            &
             & id=id(j),                         &
             & magnitude=cmplx(0.0_dp,0.0_dp,dp) )
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
    components = [ComplexSingleForce::]
  elseif (count(displacement_ids==0)==1) then
    i = first(displacement_ids==0, default=0)
    if (powers(i)>1) then
      ! If n_i>1, then the derivative along u_i is also zero.
      components = [ComplexSingleForce::]
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
    allocate(components(size(forces)), stat=ialloc); call err(ialloc)
    do i=1,size(forces)
      components(i) = this%coefficient                       &
                  & * product( energies,                     &
                  &            dim=1,                        &
                  &            mask=[(j/=i,j=1,size(this))]) &
                  & * forces(i)
    enddo
  endif
  
  ! Construct output from components along each mode.
  output = ComplexModeForce(components)
end function

impure elemental function force_ComplexPolynomial(this,displacement) &
   & result(output)
  implicit none
  
  class(ComplexPolynomial),       intent(in) :: this
  class(ComplexModeDisplacement), intent(in) :: displacement
  type(ComplexModeForce)                     :: output
  
  if (size(this)==0) then
    output = ComplexModeForce([ComplexSingleForce::])
  else
    output = sum(this%terms%force(displacement))
  endif
end function

! Multiplication and division by scalars.
impure elemental function multiply_ComplexMonomial_real(this,that) &
   & result(output)
  implicit none
  
  type(ComplexMonomial), intent(in) :: this
  real(dp),              intent(in) :: that
  type(ComplexMonomial)             :: output
  
  output = ComplexMonomial(this%coefficient*that, this%modes_)
end function

impure elemental function multiply_real_ComplexMonomial(this,that) &
   & result(output)
  implicit none
  
  real(dp),              intent(in) :: this
  type(ComplexMonomial), intent(in) :: that
  type(ComplexMonomial)             :: output
  
  output = ComplexMonomial(this*that%coefficient, that%modes_)
end function

impure elemental function multiply_ComplexMonomial_complex(this,that) &
   & result(output)
  implicit none
  
  type(ComplexMonomial), intent(in) :: this
  complex(dp),           intent(in) :: that
  type(ComplexMonomial)             :: output
  
  output = ComplexMonomial(this%coefficient*that, this%modes_)
end function

impure elemental function multiply_complex_ComplexMonomial(this,that) &
   & result(output)
  implicit none
  
  complex(dp),           intent(in) :: this
  type(ComplexMonomial), intent(in) :: that
  type(ComplexMonomial)             :: output
  
  output = ComplexMonomial(this*that%coefficient, that%modes_)
end function

impure elemental function multiply_ComplexPolynomial_real(this,that) &
   & result(output)
  implicit none
  
  type(ComplexPolynomial), intent(in) :: this
  real(dp),                intent(in) :: that
  type(ComplexPolynomial)             :: output
  
  output = ComplexPolynomial(this%terms * that)
end function

impure elemental function multiply_real_ComplexPolynomial(this,that) &
   & result(output)
  implicit none
  
  real(dp),                intent(in) :: this
  type(ComplexPolynomial), intent(in) :: that
  type(ComplexPolynomial)             :: output
  
  output = ComplexPolynomial(this * that%terms)
end function

impure elemental function multiply_ComplexPolynomial_complex(this,that) &
   & result(output)
  implicit none
  
  type(ComplexPolynomial), intent(in) :: this
  complex(dp),             intent(in) :: that
  type(ComplexPolynomial)             :: output
  
  output = ComplexPolynomial(this%terms * that)
end function

impure elemental function multiply_complex_ComplexPolynomial(this,that) &
   & result(output)
  implicit none
  
  complex(dp),             intent(in) :: this
  type(ComplexPolynomial), intent(in) :: that
  type(ComplexPolynomial)             :: output
  
  output = ComplexPolynomial(this * that%terms)
end function

impure elemental function divide_ComplexMonomial_real(this,that) &
   & result(output)
  implicit none
  
  type(ComplexMonomial), intent(in) :: this
  real(dp),              intent(in) :: that
  type(ComplexMonomial)             :: output
  
  output = ComplexMonomial(this%coefficient/that, this%modes_)
end function

impure elemental function divide_ComplexMonomial_complex(this,that) &
   & result(output)
  implicit none
  
  type(ComplexMonomial), intent(in) :: this
  complex(dp),           intent(in) :: that
  type(ComplexMonomial)             :: output
  
  output = ComplexMonomial(this%coefficient/that, this%modes_)
end function

impure elemental function divide_ComplexPolynomial_real(this,that) &
   & result(output)
  implicit none
  
  type(ComplexPolynomial), intent(in) :: this
  real(dp),                intent(in) :: that
  type(ComplexPolynomial)             :: output
  
  output = ComplexPolynomial(this%terms / that)
end function

impure elemental function divide_ComplexPolynomial_complex(this,that) &
   & result(output)
  implicit none
  
  type(ComplexPolynomial), intent(in) :: this
  complex(dp),             intent(in) :: that
  type(ComplexPolynomial)             :: output
  
  output = ComplexPolynomial(this%terms / that)
end function

! Multiplication between Monomial-like types.
! Uses a merge to keep ids in ascending order.
impure elemental function multiply_ComplexMonomialable_ComplexMonomialable( &
   & this,that) result(output)
  implicit none
  
  class(ComplexMonomialable), intent(in) :: this
  class(ComplexMonomialable), intent(in) :: that
  type(ComplexMonomial)                  :: output
  
  complex(dp)                          :: coefficient
  type(ComplexUnivariate), allocatable :: modes(:)
  
  type(ComplexMonomial) :: this_monomial
  type(ComplexMonomial) :: that_monomial
  
  integer :: i_this,i_that,i_out,ialloc
  
  this_monomial = this%to_ComplexMonomial()
  that_monomial = that%to_ComplexMonomial()
  
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
        modes(i_out) = ComplexUnivariate(                             &
           & id           = this_monomial%modes_(i_this)%id,          &
           & paired_id    = this_monomial%modes_(i_this)%paired_id,   &
           & power        = this_monomial%modes_(i_this)%power        &
           &              + that_monomial%modes_(i_that)%power,       &
           & paired_power = this_monomial%modes_(i_this)%paired_power &
           &              + that_monomial%modes_(i_that)%paired_power )
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
  
  output = ComplexMonomial(coefficient, modes)
end function

! Multiplication between polynomials and monomial-like types.
impure elemental function multiply_ComplexPolynomial_ComplexMonomialable( &
   & this,that) result(output)
  implicit none
  
  type(ComplexPolynomial),    intent(in) :: this
  class(ComplexMonomialable), intent(in) :: that
  type(ComplexPolynomial)                :: output
  
  output = ComplexPolynomial(this%terms * that)
end function

impure elemental function multiply_ComplexMonomialable_ComplexPolynomial( &
   & this,that) result(output)
  implicit none
  
  class(ComplexMonomialable), intent(in) :: this
  type(ComplexPolynomial),    intent(in) :: that
  type(ComplexPolynomial)                :: output
  
  output = ComplexPolynomial(this * that%terms)
end function

! Addition between polynomials and polynomial-like types.
impure elemental function add_ComplexPolynomialable_ComplexPolynomialable( &
   & this,that) result(output)
  implicit none
  
  class(ComplexPolynomialable), intent(in) :: this
  class(ComplexPolynomialable), intent(in) :: that
  type(ComplexPolynomial)                  :: output
  
  type(ComplexPolynomial) :: this_polynomial
  type(ComplexPolynomial) :: that_polynomial
  
  integer :: no_terms
  
  integer :: i,j,ialloc
  
  this_polynomial = this%to_ComplexPolynomial()
  that_polynomial = that%to_ComplexPolynomial()
  
  allocate( output%terms(size(this_polynomial)+size(that_polynomial)), &
          & stat=ialloc); call err(ialloc)
  output%terms(:size(this_polynomial)) = this_polynomial%terms
  no_terms = size(this_polynomial)
  do i=1,size(that_polynomial)
    j = first_equivalent( this_polynomial%terms,     &
                        & that_polynomial%terms(i),  &
                        & compare_complex_monomials, &
                        & default=0                  )
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
impure elemental function negative_ComplexPolynomialable(this) result(output)
  implicit none
  
  class(ComplexPolynomialable), intent(in) :: this
  type(ComplexPolynomial)                  :: output
  
  output = this%to_ComplexPolynomial()
  output%terms%coefficient = - output%terms%coefficient
end function

! Subtraction between polynomials and polynomial-like types.
impure elemental function                                            &
   & subtract_ComplexPolynomialable_ComplexPolynomialable(this,that) &
   & result(output)
  implicit none
  
  class(ComplexPolynomialable), intent(in) :: this
  class(ComplexPolynomialable), intent(in) :: that
  type(ComplexPolynomial)                  :: output
  
  output = this + (-that)
end function

! Sum polynomial-like types.
function sum_ComplexPolynomialables(input) result(output)
  implicit none
  
  class(ComplexPolynomialable), intent(in) :: input(:)
  type(ComplexPolynomial)                  :: output
  
  integer :: i
  
  if (size(input)==0) then
    output = ComplexPolynomial([ComplexMonomial::])
  else
    output = input(1)%to_ComplexPolynomial()
    do i=2,size(input)
      !output = output + input(i)
    enddo
  endif
end function

! ----------------------------------------------------------------------
! Select modes, displacements or forces corresponding to
!    a given univariate or univariates.
! ----------------------------------------------------------------------
function select_modes_ComplexUnivariate(input,modes) result(output)
  implicit none
  
  type(ComplexUnivariate), intent(in) :: input
  type(ComplexMode),       intent(in) :: modes(:)
  type(ComplexMode), allocatable      :: output(:)
  
  if (input%id==input%paired_id) then
    output = [modes(first(modes%id==input%id))]
  else
    output = [ modes(first(modes%id==input%id)),       &
             & modes(first(modes%id==input%paired_id)) ]
  endif
end function

function select_modes_ComplexUnivariates(input,modes) &
   & result(output)
  implicit none
  
  type(ComplexUnivariate), intent(in) :: input(:)
  type(ComplexMode),       intent(in) :: modes(:)
  type(ComplexMode), allocatable      :: output(:)
  
  integer :: i
  
  output = [ComplexMode::]
  do i=1,size(input)
    output = [output, select_modes(input(i), modes)]
  enddo
end function

function select_displacements_ComplexUnivariate(input,displacements) &
   & result(output)
  implicit none
  
  type(ComplexUnivariate),          intent(in) :: input
  type(ComplexSingleDisplacement),  intent(in) :: displacements(:)
  type(ComplexSingleDisplacement), allocatable :: output(:)
  
  if (input%id==input%paired_id) then
    output = [displacements(first(displacements%id==input%id))]
  else
    output = [ displacements(first(displacements%id==input%id)),       &
             & displacements(first(displacements%id==input%paired_id)) ]
  endif
end function

function select_displacements_ComplexUnivariates(input,displacements) &
   & result(output)
  implicit none
  
  type(ComplexUnivariate),         intent(in)  :: input(:)
  type(ComplexSingleDisplacement), intent(in)  :: displacements(:)
  type(ComplexSingleDisplacement), allocatable :: output(:)
  
  integer :: i
  
  output = [ComplexSingleDisplacement::]
  do i=1,size(input)
    output = [output, select_displacements(input(i), displacements)]
  enddo
end function

function select_forces_ComplexUnivariate(input,forces) result(output)
  implicit none
  
  type(ComplexUnivariate),  intent(in)  :: input
  type(ComplexSingleForce), intent(in)  :: forces(:)
  type(ComplexSingleForce), allocatable :: output(:)
  
  if (input%id==input%paired_id) then
    output = [forces(first(forces%id==input%id))]
  else
    output = [ forces(first(forces%id==input%id)),       &
             & forces(first(forces%id==input%paired_id)) ]
  endif
end function

function select_forces_ComplexUnivariates(input,forces) &
   & result(output)
  implicit none
  
  type(ComplexUnivariate),  intent(in)  :: input(:)
  type(ComplexSingleForce), intent(in)  :: forces(:)
  type(ComplexSingleForce), allocatable :: output(:)
  
  integer :: i
  
  output = [ComplexSingleForce::]
  do i=1,size(input)
    output = [output, select_forces(input(i), forces)]
  enddo
end function

! ----------------------------------------------------------------------
! Compares two monomials for equality up to coefficient.
! ----------------------------------------------------------------------
function compare_complex_monomials(this,that) result(output)
  implicit none
  
  class(*), intent(in) :: this
  class(*), intent(in) :: that
  logical              :: output
  
  select type(this); type is(ComplexMonomial)
    select type(that); type is(ComplexMonomial)
      if (size(this)/=size(that)) then
        output = .false.
      else
        output = all( this%modes_%id==that%modes_%id .and.               &
                    & this%modes_%power==that%modes_%power .and.         &
                    & this%modes_%paired_power==that%modes_%paired_power )
      endif
    end select
  end select
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_ComplexUnivariate(this,input)
  implicit none
  
  class(ComplexUnivariate), intent(out) :: this
  type(String),             intent(in) :: input
  
  type(String), allocatable :: line(:)
  integer                   :: id
  integer                   :: paired_id
  integer                   :: power
  integer                   :: paired_power
  
  select type(this); type is(ComplexUnivariate)
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
      call print_line(ERROR//': unable to parse ComplexUnivariate.')
      call err()
    endif
    
    this = ComplexUnivariate( id           = id,          &
                            & paired_id    = paired_id,   &
                            & power        = power,       &
                            & paired_power = paired_power )
  end select
end subroutine

function write_ComplexUnivariate(this) result(output)
  implicit none
  
  class(ComplexUnivariate), intent(in) :: this
  type(String)                         :: output
  
  select type(this); type is(ComplexUnivariate)
    if (this%id==this%paired_id) then
      output = '(u'//this%id//'^'//this%power//')'
    else
      output = '(u'//this%id//'^'//this%power// &
             & '*u'//this%paired_id//'^'//this%paired_power//')'
    endif
  end select
end function

impure elemental function new_ComplexUnivariate_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(ComplexUnivariate)  :: this
  
  call this%read(input)
end function

subroutine read_ComplexMonomial(this,input)
  implicit none
  
  class(ComplexMonomial), intent(out) :: this
  type(String),           intent(in)  :: input
  
  type(String),            allocatable :: line(:)
  complex(dp)                          :: coefficient
  type(ComplexUnivariate), allocatable :: modes(:)
  
  integer :: i,ialloc
  
  select type(this); type is(ComplexMonomial)
    ! Splitting the input by '*' separates the coefficient and the modes,
    !    but also splits some modes in two.
    line = split_line(input,delimiter='*')
    
    coefficient = cmplx(line(1))
    
    modes = [ComplexUnivariate::]
    i = 2
    do while (i<=size(line))
      ! Check if line(i) ends in a bracket.
      if (slice(line(i),len(line(i)),len(line(i)))==')') then
        ! line(i) is a mode on its own.
        modes = [modes, ComplexUnivariate(line(i))]
        i = i+1
      else
        ! line(i) and line(i+1) together make a mode.
        modes = [modes, ComplexUnivariate(line(i)//'*'//line(i+1))]
        i = i+2
      endif
    enddo
    
    this = ComplexMonomial(coefficient, modes)
  end select
end subroutine

function write_ComplexMonomial(this) result(output)
  implicit none
  
  class(ComplexMonomial), intent(in) :: this
  type(String)                       :: output
  
  select type(this); type is(ComplexMonomial)
    if (size(this%modes_)>0) then
      output = this%coefficient//'*'//join(this%modes_, delimiter='*')
    else
      output = str(this%coefficient)
    endif
  end select
end function

impure elemental function new_ComplexMonomial_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(ComplexMonomial)    :: this
  
  call this%read(input)
end function

subroutine read_ComplexPolynomial(this,input)
  implicit none
  
  class(ComplexPolynomial), intent(out) :: this
  type(String),             intent(in)  :: input
  
  type(String),          allocatable :: terms(:)
  type(String)                       :: plus
  type(ComplexMonomial), allocatable :: monomials(:)
  
  plus = '+'
  select type(this); type is(ComplexPolynomial)
    terms = split_line(input)
    terms = terms(filter(terms/=plus))
    
    monomials = ComplexMonomial(terms)
    
    this = ComplexPolynomial(monomials)
  end select
end subroutine

function write_ComplexPolynomial(this) result(output)
  implicit none
  
  class(ComplexPolynomial), intent(in) :: this
  type(String)                         :: output
  
  select type(this); type is(ComplexPolynomial)
    output = join(this%terms, delimiter=' + ')
  end select
end function

impure elemental function new_ComplexPolynomial_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(ComplexPolynomial)  :: this
  
  call this%read(input)
end function
end module
