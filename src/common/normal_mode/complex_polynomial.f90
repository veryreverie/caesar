! ======================================================================
! The building blocks of basis functions in complex co-ordinates.
! ======================================================================
module complex_polynomial_submodule
  use utils_module
  
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
  public :: conjg
  public :: select_mode
  public :: select_modes
  public :: select_displacement
  public :: select_displacements
  public :: select_force
  public :: select_forces
  
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
  end type
  
  type, abstract, extends(ComplexPolynomialable) :: ComplexMonomialable
  contains
    procedure(to_ComplexMonomial_ComplexMonomialable), deferred, public :: &
       & to_ComplexMonomial
  end type
  
  type, extends(ComplexMonomialable) :: ComplexUnivariate
    integer          :: id
    integer, private :: paired_id_
    integer          :: power
  contains
    procedure, public :: paired_id
    
    procedure, public :: to_ComplexMonomial   => &
       & to_ComplexMonomial_ComplexUnivariate
    procedure, public :: to_ComplexPolynomial => &
       & to_ComplexPolynomial_ComplexUnivariate
    
    procedure, public :: energy => energy_ComplexUnivariate
    procedure, public :: force  => force_ComplexUnivariate
    
    procedure, public :: read  => read_ComplexUnivariate
    procedure, public :: write => write_ComplexUnivariate
  end type
  
  interface ComplexUnivariate
    module procedure new_ComplexUnivariate
    module procedure new_ComplexUnivariate_String
  end interface
  
  type, extends(ComplexMonomialable) :: ComplexMonomial
    complex(dp)                          :: coefficient
    type(ComplexUnivariate), allocatable :: modes(:)
  contains
    procedure, public :: to_ComplexMonomial   => &
       & to_ComplexMonomial_ComplexMonomial
    procedure, public :: to_ComplexPolynomial => &
       & to_ComplexPolynomial_ComplexMonomial
    
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
    procedure, public :: to_ComplexPolynomial => &
       & to_ComplexPolynomial_ComplexPolynomial
    
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
    
    module procedure multiply_ComplexMonomialable_ComplexMonomialable
  end interface
  
  interface operator(/)
    module procedure divide_ComplexMonomial_real
    module procedure divide_ComplexMonomial_complex
  end interface
  
  interface operator(+)
    module procedure add_ComplexPolynomialable_ComplexPolynomialable
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
  
  integer, intent(in)           :: id
  integer, intent(in), optional :: paired_id
  integer, intent(in)           :: power
  type(ComplexUnivariate)       :: this
  
  this%id = id
  if (present(paired_id)) then
    this%paired_id_ = paired_id
  else
    this%paired_id_ = 0
  endif
  this%power = power
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
! Getter for paired_id, with error checking.
! ----------------------------------------------------------------------
impure elemental function paired_id(this) result(output)
  implicit none
  
  class(ComplexUnivariate), intent(in) :: this
  integer                              :: output
  
  if (this%paired_id_==0) then
    call print_line(CODE_ERROR//': Trying to use paired_id before it has &
       & been set.')
    call err()
  endif
  
  output = this%paired_id_
end function

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

! Find the conjugate of a univariate or monomial.
impure elemental function conjg_ComplexUnivariate(this) result(output)
  implicit none
  
  type(ComplexUnivariate), intent(in) :: this
  type(ComplexUnivariate)             :: output
  
  output = ComplexUnivariate( id        = this%paired_id(), &
                            & paired_id = this%id,          &
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
  
  output = this
  output%coefficient = output%coefficient * that
end function

impure elemental function multiply_real_ComplexMonomial(this,that) &
   & result(output)
  implicit none
  
  real(dp),              intent(in) :: this
  type(ComplexMonomial), intent(in) :: that
  type(ComplexMonomial)             :: output
  
  output = that
  output%coefficient = output%coefficient * this
end function

impure elemental function multiply_ComplexMonomial_complex(this,that) &
   & result(output)
  implicit none
  
  type(ComplexMonomial), intent(in) :: this
  complex(dp),           intent(in) :: that
  type(ComplexMonomial)             :: output
  
  output = this
  output%coefficient = output%coefficient * that
end function

impure elemental function multiply_complex_ComplexMonomial(this,that) &
   & result(output)
  implicit none
  
  complex(dp),           intent(in) :: this
  type(ComplexMonomial), intent(in) :: that
  type(ComplexMonomial)             :: output
  
  output = that
  output%coefficient = output%coefficient * this
end function

impure elemental function divide_ComplexMonomial_real(this,that) result(output)
  implicit none
  
  type(ComplexMonomial), intent(in) :: this
  real(dp),              intent(in) :: that
  type(ComplexMonomial)             :: output
  
  output = this
  output%coefficient = output%coefficient / that
end function

impure elemental function divide_ComplexMonomial_complex(this,that) &
   & result(output)
  implicit none
  
  type(ComplexMonomial), intent(in) :: this
  complex(dp),           intent(in) :: that
  type(ComplexMonomial)             :: output
  
  output = this
  output%coefficient = output%coefficient / that
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
    j = first( this_polynomial%terms,    &
             & compare_monomial_modes,   &
             & that_polynomial%terms(j), &
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
contains
  ! Lambda for comparing monomials.
  function compare_monomial_modes(this,that) result(output)
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
! I/O.
! ----------------------------------------------------------------------
subroutine read_ComplexUnivariate(this,input)
  implicit none
  
  class(ComplexUnivariate), intent(out) :: this
  type(String),             intent(in) :: input
  
  type(String), allocatable :: line(:)
  integer                   :: id
  integer                   :: power
  
  select type(this); type is(ComplexUnivariate)
    line = split_line(input,delimiter='^')
    if (size(line)/=2) then
      call print_line(ERROR//': Unable to convert string to univariate.')
      call err()
    endif
    
    id = int(slice(line(1),2,len(line(1))))
    power = int(line(2))
    
    this = ComplexUnivariate(id=id,power=power)
  end select
end subroutine

function write_ComplexUnivariate(this) result(output)
  implicit none
  
  class(ComplexUnivariate), intent(in) :: this
  type(String)                         :: output
  
  select type(this); type is(ComplexUnivariate)
    output = 'u'//this%id//'^'//this%power
  end select
end function

impure elemental function new_ComplexUnivariate_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(ComplexUnivariate)  :: this
  
  this = input
end function

subroutine read_ComplexMonomial(this,input)
  implicit none
  
  class(ComplexMonomial), intent(out) :: this
  type(String),           intent(in)  :: input
  
  type(String),            allocatable :: line(:)
  complex(dp)                          :: coefficient
  type(ComplexUnivariate), allocatable :: modes(:)
  
  select type(this); type is(ComplexMonomial)
    line = split_line(input,delimiter='*')
    
    coefficient = cmplx(line(1))
    modes = ComplexUnivariate(line(2:))
    
    this = ComplexMonomial(coefficient, modes)
  end select
end subroutine

function write_ComplexMonomial(this) result(output)
  implicit none
  
  class(ComplexMonomial), intent(in) :: this
  type(String)                       :: output
  
  select type(this); type is(ComplexMonomial)
    output = this%coefficient//'*'//join(this%modes, delimiter='*')
  end select
end function

impure elemental function new_ComplexMonomial_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(ComplexMonomial)    :: this
  
  this = input
end function

subroutine read_ComplexPolynomial(this,input)
  implicit none
  
  class(ComplexPolynomial), intent(out) :: this
  type(String),             intent(in)  :: input
  
  type(String), allocatable :: terms(:)
  type(String)              :: plus
  
  plus = '+'
  select type(this); type is(ComplexPolynomial)
    terms = split_line(input)
    terms = terms(filter(terms/=plus))
    this = ComplexPolynomial(ComplexMonomial(terms))
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
  
  this = input
end function
end module
