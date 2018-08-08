! ======================================================================
! The building blocks of basis functions in complex co-ordinates.
! ======================================================================
module complex_polynomial_submodule
  use utils_module
  use structure_module
  
  use complex_mode_submodule
  use complex_single_mode_displacement_submodule
  use complex_single_mode_force_submodule
  use complex_mode_displacement_submodule
  use complex_mode_force_submodule
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
  public :: select_mode
  public :: select_modes
  public :: select_displacement
  public :: select_displacements
  public :: select_force
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
  contains
    procedure, private :: set_paired_id
    
    procedure, public :: to_ComplexMonomial   => &
       & to_ComplexMonomial_ComplexUnivariate
    procedure, public :: to_ComplexPolynomial => &
       & to_ComplexPolynomial_ComplexUnivariate
    
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
    complex(dp)                          :: coefficient
    type(ComplexUnivariate), allocatable :: modes(:)
  contains
    procedure, private :: set_paired_ids => set_paired_ids_ComplexMonomial
    
    procedure, public :: to_ComplexMonomial   => &
       & to_ComplexMonomial_ComplexMonomial
    procedure, public :: to_ComplexPolynomial => &
       & to_ComplexPolynomial_ComplexMonomial
    
    procedure, public :: simplify => simplify_ComplexMonomial
    
    procedure, public :: wavevector => wavevector_ComplexMonomial
    
    procedure, public :: energy => energy_ComplexMonomial
    procedure, public :: force  => force_ComplexMonomial
    
    procedure, public :: read  => read_ComplexMonomial
    procedure, public :: write => write_ComplexMonomial
  end type
  
  interface ComplexMonomial
    module procedure new_ComplexMonomial
    module procedure new_ComplexMonomial_String
  end interface
  
  type, extends(ComplexPolynomialable) :: ComplexPolynomial
    type(ComplexMonomial), allocatable :: terms(:)
  contains
    procedure, private :: set_paired_ids => set_paired_ids_ComplexPolynomial
    
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
  
  interface select_mode
    module procedure select_mode_ComplexUnivariate
  end interface
  
  interface select_modes
    module procedure select_modes_ComplexUnivariates
  end interface
  
  interface select_displacement
    module procedure select_displacement_ComplexUnivariate
  end interface
  
  interface select_displacements
    module procedure select_displacements_ComplexUnivariates
  end interface
  
  interface select_force
    module procedure select_force_ComplexUnivariate
  end interface
  
  interface select_forces
    module procedure select_forces_ComplexUnivariates
  end interface
contains

! ----------------------------------------------------------------------
! Constructors.
! ----------------------------------------------------------------------
function new_ComplexUnivariate(id,paired_id,power) result(this)
  implicit none
  
  integer, intent(in)     :: id
  integer, intent(in)     :: paired_id
  integer, intent(in)     :: power
  type(ComplexUnivariate) :: this
  
  this%id        = id
  this%paired_id = paired_id
  this%power     = power
end function

function new_ComplexUnivariate_ComplexMode(mode,power) result(this)
  implicit none
  
  type(ComplexMode), intent(in) :: mode
  integer,           intent(in) :: power
  type(ComplexUnivariate)       :: this
  
  this = ComplexUnivariate( id        = mode%id,        &
                          & paired_id = mode%paired_id, &
                          & power     = power           )
end function

function new_ComplexMonomial(coefficient,modes) result(this)
  implicit none
  
  complex(dp),             intent(in) :: coefficient
  type(ComplexUnivariate), intent(in) :: modes(:)
  type(ComplexMonomial)               :: this
  
  this%coefficient = coefficient
  this%modes       = modes
end function

function new_ComplexPolynomial(terms) result(this)
  implicit none
  
  type(ComplexMonomial), intent(in) :: terms(:)
  type(ComplexPolynomial)           :: this
  
  this%terms = terms
end function

! ----------------------------------------------------------------------
! Setter for paired_id, using an array of modes.
! ----------------------------------------------------------------------
subroutine set_paired_id(this,modes)
  implicit none
  
  class(ComplexUnivariate), intent(inout) :: this
  type(ComplexMode),        intent(in)    :: modes(:)
  
  type(ComplexMode) :: mode
  
  mode = modes(first(modes%id==this%id))
  this%paired_id = mode%paired_id
end subroutine

subroutine set_paired_ids_ComplexMonomial(this,modes)
  implicit none
  
  class(ComplexMonomial), intent(inout) :: this
  type(ComplexMode),      intent(in)    :: modes(:)
  
  integer :: i
  
  do i=1,size(this%modes)
    call this%modes(i)%set_paired_id(modes)
  enddo
end subroutine

subroutine set_paired_ids_ComplexPolynomial(this,modes)
  implicit none
  
  class(ComplexPolynomial), intent(inout) :: this
  type(ComplexMode),        intent(in)    :: modes(:)
  
  integer :: i
  
  do i=1,size(this%terms)
    call this%terms(i)%set_paired_ids(modes)
  enddo
end subroutine

! ----------------------------------------------------------------------
! Conversions between types.
! ----------------------------------------------------------------------
function to_ComplexMonomial_ComplexUnivariate(this) result(output)
  implicit none
  
  class(ComplexUnivariate), intent(in) :: this
  type(ComplexMonomial)                :: output
  
  output = ComplexMonomial( coefficient = cmplx(1.0_dp,0.0_dp,dp), &
                          & modes       = [this])
end function

function to_ComplexPolynomial_ComplexUnivariate(this) result(output)
  implicit none
  
  class(ComplexUnivariate), intent(in) :: this
  type(ComplexPolynomial)              :: output
  
  output = ComplexPolynomial([this%to_ComplexMonomial()])
end function

function to_ComplexMonomial_ComplexMonomial(this) result(output)
  implicit none
  
  class(ComplexMonomial), intent(in) :: this
  type(ComplexMonomial)              :: output
  
  output = this
end function

function to_ComplexPolynomial_ComplexMonomial(this) result(output)
  implicit none
  
  class(ComplexMonomial), intent(in) :: this
  type(ComplexPolynomial)            :: output
  
  output = ComplexPolynomial([this])
end function

function to_ComplexPolynomial_ComplexPolynomial(this) result(output)
  implicit none
  
  class(ComplexPolynomial), intent(in) :: this
  type(ComplexPolynomial)              :: output
  
  output = this
end function

! ----------------------------------------------------------------------
! Operations involving types.
! ----------------------------------------------------------------------

! The number of modes in a monomial, or terms in a polynomial.
function size_ComplexMonomial(this) result(output)
  implicit none
  
  type(ComplexMonomial), intent(in) :: this
  integer                           :: output
  
  output = size(this%modes)
end function

function size_ComplexPolynomial(this) result(output)
  implicit none
  
  class(ComplexPolynomial), intent(in) :: this
  integer                              :: output
  
  output = size(this%terms)
end function

! Simplify a monomial or polynomial.
impure elemental subroutine simplify_ComplexMonomial(this)
  implicit none
  
  class(ComplexMonomial), intent(inout) :: this
  
  integer :: i
  
  ! Sort modes in ascending order of ID.
  this%modes = this%modes(sort(this%modes%id))
  
  ! Combine modes with the same ID, remove modes with power=0.
  i = 1
  do while(i<=size(this))
    if (this%modes(i)%power<0) then
      call err()
    elseif (this%modes(i)%power==0) then
      this%modes = [this%modes(:i-1), this%modes(i+1:)]
      cycle
    endif
    
    if (i>1) then
      if (this%modes(i)%id==this%modes(i-1)%id) then
        this%modes(i-1)%power = this%modes(i-1)%power + this%modes(i)%power
        this%modes = [this%modes(:i-1), this%modes(i+1:)]
        cycle
      endif
    endif
    
    i = i+1
  enddo
end subroutine

impure elemental subroutine simplify_ComplexPolynomial(this)
  implicit none
  
  class(ComplexPolynomial), intent(inout) :: this
  
  integer,               allocatable :: equivalent_monomial_locs(:)
  type(ComplexMonomial), allocatable :: equivalent_monomials(:)
  type(ComplexMonomial), allocatable :: monomials(:)
  
  real(dp) :: max_coefficient
  
  integer :: i,ialloc
  
  call this%terms%simplify()
  
  ! Add together any equivalent monomials.
  equivalent_monomial_locs = first_equivalent( this%terms,               &
                                             & compare_complex_monomials )
  allocate( monomials(maxval(equivalent_monomial_locs)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(monomials)
    equivalent_monomials = this%terms(filter(equivalent_monomial_locs==i))
    monomials(i) = equivalent_monomials(1)
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
  
  output = ComplexUnivariate( id        = this%paired_id, &
                            & paired_id = this%id,        &
                            & power     = this%power)
end function

impure elemental function conjg_ComplexMonomial(this) result(output)
  implicit none
  
  type(ComplexMonomial), intent(in) :: this
  type(ComplexMonomial)             :: output
  
  output = ComplexMonomial( coefficient = conjg(this%coefficient), &
                          & modes       = conjg(this%modes))
  output%modes = output%modes(sort(output%modes%id))
end function

impure elemental function conjg_ComplexPolynomial(this) result(output)
  implicit none
  
  type(ComplexPolynomial), intent(in) :: this
  type(ComplexPolynomial)             :: output
  
  output = ComplexPolynomial(conjg(this%terms))
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
  qpoint%qpoint = qpoint%qpoint * this%power
  
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
  do i=1,size(this%modes)
    mode_wavevector = this%modes(i)%wavevector(modes,qpoints)
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
  
  if (this%id/=displacement%id) then
    call print_line(CODE_ERROR//': Trying to evaluate a univariate at an &
       &incompatible displacement.')
    call err()
  endif
  
  output = displacement%magnitude**this%power
end function

impure elemental function energy_ComplexMonomial(this,displacement) &
   & result(output)
  implicit none
  
  class(ComplexMonomial),         intent(in) :: this
  class(ComplexModeDisplacement), intent(in) :: displacement
  complex(dp)                                :: output
  
  integer :: i,j
  
  output = this%coefficient
  
  do i=1,size(this)
    ! Find the mode in the displacement which matches that in the monomial.
    j = first(displacement%vectors%id==this%modes(i)%id,default=0)
    
    ! If the mode is not present in the displacement,
    !    then the displacement along that mode is zero.
    ! As such, the monomial is zero. (0^n=0 if n>0).
    if (j==0) then
      output = 0.0_dp
      return
    endif
    
    ! If the mode is present in both,
    !    evaluate the univariate at the displacement.
    output = output * this%modes(i)%energy(displacement%vectors(j))
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
  
  if (this%id/=displacement%id) then
    call print_line(CODE_ERROR//': Trying to take the derivative of a &
       & univariate at an incompatible displacement.')
    call err()
  endif
  
  if (this%power<1) then
    call err()
  elseif (this%power==1) then
    force = -1.0_dp
  else
    force = -this%power * displacement%magnitude**(this%power-1)
  endif
  
  output = ComplexSingleForce(id=this%id, magnitude=force)
end function

! -d/d{u_i} (c*prod_j[ {u_j}^{n_j} ]) evaluated at {u_i=U_i} is
!    -c*prod_{j/=i}[ {U_j}^{n_j} ] * n_i * {U_i}^{n_i-1}
impure elemental function force_ComplexMonomial(this,displacement) &
   & result(output)
  implicit none
  
  class(ComplexMonomial),         intent(in) :: this
  class(ComplexModeDisplacement), intent(in) :: displacement
  type(ComplexModeForce)                     :: output
  
  integer,                  allocatable :: displacement_ids(:)
  complex(dp),              allocatable :: evaluations(:)
  type(ComplexSingleForce), allocatable :: forces(:)
  type(ComplexSingleForce), allocatable :: components(:)
  
  integer :: i,j,ialloc
  
  ! Match the components of the displacement along each mode with the
  !    univariates making up the monomial.
  ! Evaluate and take the derivative of each univariate at the one-mode
  !    component of the vector.
  allocate( displacement_ids(size(this)),  &
          & evaluations(size(this)),       &
          & forces(size(this)),            &
          & stat=ialloc); call err(ialloc)
  do i=1,size(this)
    ! Identify the vector corresponding to each mode.
    ! If vector_ids(i)=0 then U_i=0.
    displacement_ids(i) = first( displacement%vectors%id==this%modes(i)%id, &
                               & default=0)
    
    ! Calculate {U_i}^{n_i}
    if (displacement_ids(i)==0) then
      evaluations(i) = 0.0_dp
    else
      evaluations(i) = this%modes(i)%energy(         &
         & displacement%vectors(displacement_ids(i)) )
    endif
    
    ! Calculate -d/d{u_i} ({u_i}^{n_i}) evaluated at U_i.
    if (displacement_ids(i)==0) then
      if (this%modes(i)%power==1) then
        forces(i) = ComplexSingleForce(         &
           & id=this%modes(i)%id,               &
           & magnitude=cmplx(-1.0_dp,0.0_dp,dp) )
      else
        forces(i) = ComplexSingleForce(        &
           & id=this%modes(i)%id,              &
           & magnitude=cmplx(0.0_dp,0.0_dp,dp) )
      endif
    else
      forces(i) = this%modes(i)%force(               &
         & displacement%vectors(displacement_ids(i)) )
    endif
  enddo
  
  ! Use the Univariate terms to calculate forces along each mode.
  if (count(displacement_ids==0)>1) then
    ! If U_i is zero for more than one i, then
    !    prod_{j/=i}[ {U_j}^{n_j} ] is always zero,
    !    so all derivatives are zero.
    components = [ComplexSingleForce::]
  elseif (count(displacement_ids==0)==1) then
    i = first(displacement_ids==0, default=0)
    if (this%modes(i)%power>1) then
      ! If n_i>1, then the derivative along u_i is also zero.
      components = [ComplexSingleForce::]
    else
      ! If n_i=1, then the derivative is simply c*prod_{j/=i}[ {U_j}^{n_j}
      components = [ this%coefficient                       &
                 & * product( evaluations,                  &
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
                  & * product( evaluations,                  &
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
  
  output = ComplexMonomial(this%coefficient*that, this%modes)
end function

impure elemental function multiply_real_ComplexMonomial(this,that) &
   & result(output)
  implicit none
  
  real(dp),              intent(in) :: this
  type(ComplexMonomial), intent(in) :: that
  type(ComplexMonomial)             :: output
  
  output = ComplexMonomial(this*that%coefficient, that%modes)
end function

impure elemental function multiply_ComplexMonomial_complex(this,that) &
   & result(output)
  implicit none
  
  type(ComplexMonomial), intent(in) :: this
  complex(dp),           intent(in) :: that
  type(ComplexMonomial)             :: output
  
  output = ComplexMonomial(this%coefficient*that, this%modes)
end function

impure elemental function multiply_complex_ComplexMonomial(this,that) &
   & result(output)
  implicit none
  
  complex(dp),           intent(in) :: this
  type(ComplexMonomial), intent(in) :: that
  type(ComplexMonomial)             :: output
  
  output = ComplexMonomial(this*that%coefficient, that%modes)
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
  
  output = ComplexMonomial(this%coefficient/that, this%modes)
end function

impure elemental function divide_ComplexMonomial_complex(this,that) &
   & result(output)
  implicit none
  
  type(ComplexMonomial), intent(in) :: this
  complex(dp),           intent(in) :: that
  type(ComplexMonomial)             :: output
  
  output = ComplexMonomial(this%coefficient/that, this%modes)
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

! Multiplication between Monomials and Monomial-like types.
function multiply_ComplexMonomialable_ComplexMonomialable(this,that) &
   & result(output)
  implicit none
  
  class(ComplexMonomialable), intent(in) :: this
  class(ComplexMonomialable), intent(in) :: that
  type(ComplexMonomial)                  :: output
  
  type(ComplexMonomial) :: this_monomial
  type(ComplexMonomial) :: that_monomial
  
  integer :: i_this,i_that,i_out,ialloc
  
  this_monomial = this%to_ComplexMonomial()
  that_monomial = that%to_ComplexMonomial()
  
  output%coefficient = this_monomial%coefficient * that_monomial%coefficient
  
  if (size(this_monomial)==0) then
    output%modes = that_monomial%modes
  elseif (size(that_monomial)==0) then
    output%modes = this_monomial%modes
  else
    i_this = 1
    i_that = 1
    i_out = 0
    allocate( output%modes(size(this_monomial)+size(that_monomial)), &
            & stat=ialloc); call err(ialloc)
    do while(i_this<=size(this_monomial) .and. i_that<=size(that_monomial))
      i_out = i_out + 1
      if (i_this>size(this_monomial)) then
        output%modes(i_out) = that_monomial%modes(i_that)
        i_that = i_that + 1
      elseif (i_that>size(that_monomial)) then
        output%modes(i_out) = this_monomial%modes(i_this)
        i_this = i_this + 1
      elseif ( this_monomial%modes(i_this)%id == &
             & that_monomial%modes(i_that)%id) then
        output%modes(i_out)%id = this_monomial%modes(i_this)%id
        output%modes(i_out)%power = this_monomial%modes(i_this)%power &
                                & + that_monomial%modes(i_that)%power
        i_this = i_this + 1
        i_that = i_that + 1
      elseif ( this_monomial%modes(i_this)%id < &
             & that_monomial%modes(i_that)%id) then
        output%modes(i_out) = this_monomial%modes(i_this)
        i_this = i_this + 1
      elseif ( this_monomial%modes(i_this)%id > &
             & that_monomial%modes(i_that)%id) then
        output%modes(i_out) = that_monomial%modes(i_that)
        i_that = i_that + 1
      else
        call err()
      endif
    enddo
    output%modes = output%modes(:i_out)
  endif
end function

! Addition between polynomials and polynomial-like types.
function add_ComplexPolynomialable_ComplexPolynomialable(this,that) &
   & result(output)
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
    j = first( this_polynomial%terms,     &
             & compare_complex_monomials, &
             & that_polynomial%terms(i),  &
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
function negative_ComplexPolynomialable(this) result(output)
  implicit none
  
  class(ComplexPolynomialable), intent(in) :: this
  type(ComplexPolynomial)                  :: output
  
  output = this%to_ComplexPolynomial()
  output%terms%coefficient = - output%terms%coefficient
end function

! Subtraction between polynomials and polynomial-like types.
function subtract_ComplexPolynomialable_ComplexPolynomialable(this,that) &
   & result(output)
  implicit none
  
  class(ComplexPolynomialable), intent(in) :: this
  class(ComplexPolynomialable), intent(in) :: that
  type(ComplexPolynomial)                  :: output
  
  output = this + (-that)
end function

! ----------------------------------------------------------------------
! Select modes, displacements or forces corresponding to
!    a given univariate or univariates.
! ----------------------------------------------------------------------
function select_mode_ComplexUnivariate(univariate,modes) result(output)
  implicit none
  
  type(ComplexUnivariate), intent(in) :: univariate
  type(ComplexMode),       intent(in) :: modes(:)
  type(ComplexMode)                   :: output
  
  output = modes(first(modes%id==univariate%id))
end function

function select_modes_ComplexUnivariates(univariates,modes) &
   & result(output)
  implicit none
  
  type(ComplexUnivariate), intent(in) :: univariates(:)
  type(ComplexMode),       intent(in) :: modes(:)
  type(ComplexMode), allocatable      :: output(:)
  
  integer :: i,ialloc
  
  allocate(output(size(univariates)), stat=ialloc); call err(ialloc)
  do i=1,size(univariates)
    output(i) = select_mode(univariates(i), modes)
  enddo
end function

function select_displacement_ComplexUnivariate(univariate,displacements) &
   & result(output)
  implicit none
  
  type(ComplexUnivariate),         intent(in) :: univariate
  type(ComplexSingleDisplacement), intent(in) :: displacements(:)
  type(ComplexSingleDisplacement)             :: output
  
  output = displacements(first(displacements%id==univariate%id))
end function

function select_displacements_ComplexUnivariates(univariates,displacements) &
   & result(output)
  implicit none
  
  type(ComplexUnivariate),         intent(in)  :: univariates(:)
  type(ComplexSingleDisplacement), intent(in)  :: displacements(:)
  type(ComplexSingleDisplacement), allocatable :: output(:)
  
  integer :: i,ialloc
  
  allocate(output(size(univariates)), stat=ialloc); call err(ialloc)
  do i=1,size(univariates)
    output(i) = select_displacement(univariates(i), displacements)
  enddo
end function

function select_force_ComplexUnivariate(univariate,forces) &
   & result(output)
  implicit none
  
  type(ComplexUnivariate),  intent(in) :: univariate
  type(ComplexSingleForce), intent(in) :: forces(:)
  type(ComplexSingleForce)             :: output
  
  output = forces(first(forces%id==univariate%id))
end function

function select_forces_ComplexUnivariates(univariates,forces) &
   & result(output)
  implicit none
  
  type(ComplexUnivariate),  intent(in)  :: univariates(:)
  type(ComplexSingleForce), intent(in)  :: forces(:)
  type(ComplexSingleForce), allocatable :: output(:)
  
  integer :: i,ialloc
  
  allocate(output(size(univariates)), stat=ialloc); call err(ialloc)
  do i=1,size(univariates)
    output(i) = select_force(univariates(i), forces)
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
      if (size(this%modes)/=size(that%modes)) then
        output = .false.
      else
        output = all( this%modes%id==that%modes%id .and. &
                    & this%modes%power==that%modes%power)
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
  
  select type(this); type is(ComplexUnivariate)
    ! If id=5, paired_id=7 and power=3 then:
    ! input = '(u5=u7*)^7'
    
    ! Split off power.
    line = split_line(input,delimiter='^') ! line = ['(u5=u7*)', '7']
    if (size(line)/=2) then
      call print_line(ERROR//': Unable to convert string to univariate.')
      call err()
    endif
    power = int(line(2))
    
    ! Split id and paired_id
    line = split_line(line(1), delimiter='u') ! line = ['(', '5=', '7*)']
    if (size(line)/=3) then
      call print_line(ERROR//': Unable to convert string to univariate.')
      call err()
    endif
    id = int(slice(line(2),1,len(line(2))-1))
    paired_id = int(slice(line(3),1,len(line(3))-2))
    
    this = ComplexUnivariate(id=id, paired_id=paired_id, power=power)
  end select
end subroutine

function write_ComplexUnivariate(this) result(output)
  implicit none
  
  class(ComplexUnivariate), intent(in) :: this
  type(String)                         :: output
  
  select type(this); type is(ComplexUnivariate)
    output = '(u'//this%id//'=u'//this%paired_id//'*)^'//this%power
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
    ! Splitting the input by '*' separates univariates, but also splits
    !    each univariate in two.
    line = split_line(input,delimiter='*')
    
    coefficient = cmplx(line(1))
    
    ! Each univariate must be re-assembled from its two parts.
    allocate(modes((size(line)-1)/2), stat=ialloc); call err(ialloc)
    do i=1,size(modes)
      modes(i) = ComplexUnivariate(line(2*i)//'*'//line(2*i+1))
    enddo
    
    this = ComplexMonomial(coefficient, modes)
  end select
end subroutine

function write_ComplexMonomial(this) result(output)
  implicit none
  
  class(ComplexMonomial), intent(in) :: this
  type(String)                       :: output
  
  select type(this); type is(ComplexMonomial)
    if (size(this%modes)>0) then
      output = this%coefficient//'*'//join(this%modes, delimiter='*')
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
