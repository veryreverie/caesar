! ======================================================================
! The building blocks of basis functions in complex co-ordinates.
! ======================================================================
module complex_polynomial_module
  use utils_module
  use structure_module
  
  use real_mode_module
  use real_single_mode_displacement_module
  use real_single_mode_force_module
  use real_mode_displacement_module
  use real_mode_force_module
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
    
    procedure, private :: energy_real => &
                        & energy_RealSingleDisplacement_ComplexUnivariate
    procedure, private :: energy_complex => &
                        & energy_ComplexSingleDisplacement_ComplexUnivariate
    
    procedure, private :: force_real => &
                        & force_RealSingleDisplacement_ComplexUnivariate
    procedure, private :: force_complex => &
                        & force_ComplexSingleDisplacement_ComplexUnivariate
    
    procedure, public :: harmonic_expectation => &
                       & harmonic_expectation_ComplexUnivariate
    
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
    
    procedure, public :: set_modes => set_modes_ComplexMonomial
    
    procedure, public :: simplify => simplify_ComplexMonomial
    
    procedure, public :: total_power => total_power_ComplexMonomial
    
    procedure, public :: wavevector => wavevector_ComplexMonomial
    
    generic,   public  :: energy =>                                    &
                        & energy_RealModeDisplacement_ComplexMonomial, &
                        & energy_ComplexModeDisplacement_ComplexMonomial
    procedure, private :: energy_RealModeDisplacement_ComplexMonomial
    procedure, private :: energy_ComplexModeDisplacement_ComplexMonomial
    
    generic,   public  :: force =>                                    &
                        & force_RealModeDisplacement_ComplexMonomial, &
                        & force_ComplexModeDisplacement_ComplexMonomial
    procedure, private :: force_RealModeDisplacement_ComplexMonomial
    procedure, private :: force_ComplexModeDisplacement_ComplexMonomial
    
    procedure, public :: harmonic_expectation => &
                       & harmonic_expectation_ComplexMonomial
    
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
    
    generic,   public  :: energy =>                                      &
                        & energy_RealModeDisplacement_ComplexPolynomial, &
                        & energy_ComplexModeDisplacement_ComplexPolynomial
    procedure, private :: energy_RealModeDisplacement_ComplexPolynomial
    procedure, private :: energy_ComplexModeDisplacement_ComplexPolynomial
    
    generic,   public  :: force =>                                      &
                        & force_RealModeDisplacement_ComplexPolynomial, &
                        & force_ComplexModeDisplacement_ComplexPolynomial
    procedure, private :: force_RealModeDisplacement_ComplexPolynomial
    procedure, private :: force_ComplexModeDisplacement_ComplexPolynomial
    
    procedure, public :: harmonic_expectation => &
                       & harmonic_expectation_ComplexPolynomial
    
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
      import FractionVector
      implicit none
      
      class(ComplexPolynomialable), intent(in) :: this
      type(ComplexMode),            intent(in) :: modes(:)
      type(QpointData),             intent(in) :: qpoints(:)
      type(FractionVector)                     :: output
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
    module procedure negative_ComplexMonomial
    module procedure negative_ComplexPolynomial
    
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

impure elemental function new_ComplexUnivariate_ComplexMode(mode,power, &
   & paired_power) result(this)
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

! Getter for modes. Has several run types:
!    - If no arguments specified, returns all modes.
!    - If indices specified, returns modes at specified indices.
!    - If ids and paired_ids specified, returns modes with specified ids.
!    - If exclude_ids specicified, returns all modes but those excluded.
function modes_ComplexMonomial(this,indices,ids,paired_ids,exclude_ids) &
   & result(output)
  implicit none
  
  class(ComplexMonomial), intent(in)           :: this
  integer,                intent(in), optional :: indices(:)
  integer,                intent(in), optional :: ids(:)
  integer,                intent(in), optional :: paired_ids(:)
  integer,                intent(in), optional :: exclude_ids(:)
  type(ComplexUnivariate), allocatable         :: output(:)
  
  integer, allocatable :: sort_key(:)
  
  integer :: i,j,k,ialloc
  
  if (count([present(indices),present(ids),present(exclude_ids)])>1) then
    call print_line(CODE_ERROR//': too many optional arguments present.')
    call err()
  elseif (present(ids) .neqv. present(paired_ids)) then
    call print_line(CODE_ERROR// &
       & ': Only one of ids and paired_ids is present.')
    call err()
  endif
  
  if (present(indices)) then
    ! Return the modes with the specified indices.
    output = this%modes_(indices)
  elseif (present(ids)) then
    ! Return only modes with given ids, in the order of the ids specified.
    ! If an id is not present in this monomial, instead include u^0.
    sort_key = sort(ids)
    allocate(output(size(ids)), stat=ialloc); call err(ialloc)
    j = 1
    ! Loop over the ids in ascending order (i.e. ids(sort_key)).
    do i=1,size(ids)
      ! Cycle through this%modes_ until the mode with id=ids(sort_key(i))
      !    is found.
      do
        if (j>size(this%modes_)) then
          exit
        elseif (this%modes_(j)%id>=ids(sort_key(i))) then
          exit
        else
          j = j+1
        endif
      enddo
      
      ! If such a mode exists, add it to the output.
      if (j<=size(this%modes_)) then
        if (this%modes_(j)%id == ids(sort_key(i))) then
          output(sort_key(i)) = this%modes_(j)
          cycle
        endif
      endif
      
      ! If no such mode exists, add u^0 to the output.
      output(sort_key(i)) = ComplexUnivariate(     &
         & id           = ids(sort_key(i)),        &
         & paired_id    = paired_ids(sort_key(i)), &
         & power        = 0,                       &
         & paired_power = 0                        )
    enddo
  elseif (present(exclude_ids)) then
    ! Return all modes apart from those in exclude_ids.
    sort_key = sort(exclude_ids)
    
    allocate(output(size(this%modes_)), stat=ialloc); call err(ialloc)
    j = 1
    k = 0
    ! Loop over this%modes_.
    do i=1,size(this%modes_)
      ! Cycle through exclude_ids in ascending order
      !    (i.e. exclude_ids(sort_key)) until
      !    exculde_ids(j)=this%modes_(i)%id.
      do
        if (j>size(exclude_ids)) then
          exit
        elseif (exclude_ids(sort_key(j))>=this%modes_(i)%id) then
          exit
        else
          j = j+1
        endif
      enddo
      
      ! If the mode is in exclude_ids, exclude it.
      if (j<=size(exclude_ids)) then
        if (exclude_ids(sort_key(j))==this%modes_(i)%id) then
          cycle
        endif
      endif
      
      ! If the mode is not in exclude_ids, add it to the output.
      k = k+1
      output(k) = this%modes_(i)
    enddo
    output = output(:k)
  else
    ! Return all modes in the monomial.
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

! Set the modes of a monomial.
subroutine set_modes_ComplexMonomial(this,modes)
  implicit none
  
  class(ComplexMonomial),  intent(inout) :: this
  type(ComplexUnivariate), intent(in)    :: modes(:)
  
  this%modes_ = modes
end subroutine

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
  
  integer,     allocatable :: unique_terms(:)
  complex(dp), allocatable :: coefficients(:)
  
  integer :: i,j
  
  call this%terms%simplify()
  
  ! Add together monomials with the same powers.
  unique_terms = set(this%terms, compare_complex_monomials)
  coefficients = [(cmplx(0.0_dp,0.0_dp,dp), i=1, size(unique_terms))]
  do i=1,size(this%terms)
    do j=1,size(unique_terms)
      if (compare_complex_monomials( this%terms(i),              &
                                   & this%terms(unique_terms(j)) )) then
        coefficients(j) = coefficients(j) + this%terms(i)%coefficient
      endif
    enddo
  enddo
  
  this%terms = this%terms(unique_terms)
  this%terms%coefficient = coefficients
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
  type(FractionVector)                 :: output
  
  type(ComplexMode) :: mode
  type(QpointData)  :: qpoint
  
  mode = modes(first(modes%id==this%id))
  qpoint = qpoints(first(qpoints%id==mode%qpoint_id))
  if (this%id==this%paired_id) then
    output = qpoint%qpoint * this%power
  else
    output = qpoint%qpoint * (this%power-this%paired_power)
  endif
  
  output = vec(modulo(frac(output),1))
end function

function wavevector_ComplexMonomial(this,modes,qpoints) result(output)
  implicit none
  
  class(ComplexMonomial), intent(in) :: this
  type(ComplexMode),      intent(in) :: modes(:)
  type(QpointData),       intent(in) :: qpoints(:)
  type(FractionVector)               :: output
  
  integer :: i
  
  output = fracvec(zeroes(3))
  do i=1,size(this%modes_)
    output = output + this%modes_(i)%wavevector(modes,qpoints)
  enddo
  
  output = vec(modulo(frac(output),1))
end function

function wavevector_ComplexPolynomial(this,modes,qpoints) result(output)
  implicit none
  
  class(ComplexPolynomial), intent(in) :: this
  type(ComplexMode),        intent(in) :: modes(:)
  type(QpointData),         intent(in) :: qpoints(:)
  type(FractionVector)                 :: output
  
  integer :: i
  
  if (size(this)==0) then
    output = fracvec(zeroes(3))
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
impure elemental function energy_RealSingleDisplacement_ComplexUnivariate( &
   & this,displacement,paired_displacement) result(output)
  implicit none
  
  class(ComplexUnivariate),     intent(in)           :: this
  type(RealSingleDisplacement), intent(in), optional :: displacement
  type(RealSingleDisplacement), intent(in), optional :: paired_displacement
  complex(dp)                                        :: output
  
  complex(dp) :: magnitude
  complex(dp) :: paired_magnitude
  
  if (.not. (present(displacement) .or. present(paired_displacement))) then
    call print_line(CODE_ERROR//': Neither displacement passed.')
    call err()
  elseif (present(paired_displacement) .and. this%id==this%paired_id) then
    call print_line(CODE_ERROR//': paired displacement passed to &
       &univariate with id=paired_id.')
    call err()
  endif
  
  if (this%id==this%paired_id) then
    magnitude = displacement%magnitude
  else
    ! Convert from real to complex co-ordinates.
    ! x_+ = (x_c + ix_s) / sqrt(2)
    ! x_- = (x_c - ix_s) / sqrt(2)
    magnitude = 0
    paired_magnitude = 0
    if (present(displacement)) then
      ! x_+ : x_c / sqrt(2)
      magnitude = magnitude &
              & + displacement%magnitude &
              & / sqrt(2.0_dp)
      ! x_- : x_c / sqrt(2)
      paired_magnitude = paired_magnitude &
                     & + displacement%magnitude &
                     & / sqrt(2.0_dp)
    endif
    if (present(paired_displacement)) then
      ! x_+ : ix_s / sqrt(2)
      magnitude = magnitude &
              & + paired_displacement%magnitude &
              & * cmplx(0,1,dp) / sqrt(2.0_dp)
      ! x_- : - ix_s / sqrt(2)
      paired_magnitude = paired_magnitude &
                     & - paired_displacement%magnitude &
                     & * cmplx(0,1,dp) / sqrt(2.0_dp)
    endif
  endif
  
  output = magnitude**this%power
  if (this%id/=this%paired_id) then
    output = output * paired_magnitude**this%paired_power
  endif
end function

impure elemental function energy_ComplexSingleDisplacement_ComplexUnivariate( &
   & this,displacement,paired_displacement) result(output)
  implicit none
  
  class(ComplexUnivariate),         intent(in)           :: this
  class(ComplexSingleDisplacement), intent(in), optional :: displacement
  class(ComplexSingleDisplacement), intent(in), optional :: paired_displacement
  complex(dp)                                            :: output
  
  if (.not. (present(displacement) .or. present(paired_displacement))) then
    call print_line(CODE_ERROR//': Neither displacement passed.')
    call err()
  elseif (present(paired_displacement) .and. this%id==this%paired_id) then
    call print_line(CODE_ERROR//': paired displacement passed to &
       &univariate with id=paired_id.')
    call err()
  endif
  
  output = 1.0_dp
  
  if (present(displacement)) then
    output = output * displacement%magnitude**this%power
  endif
  
  if (present(paired_displacement)) then
    output = output * paired_displacement%magnitude**this%paired_power
  endif
end function

impure elemental function energy_RealModeDisplacement_ComplexMonomial(this, &
   & displacement) result(output)
  implicit none
  
  class(ComplexMonomial),     intent(in) :: this
  type(RealModeDisplacement), intent(in) :: displacement
  complex(dp)                            :: output
  
  integer :: i,j,k
  
  output = this%coefficient
  
  ! Loop over the modes in this monomial,
  !    find the mode (and paired mode if relevant) in the displacement
  !    corresponding to each mode,
  !    and calculate the energy term by term.
  do i=1,size(this)
    associate (mode => this%modes_(i))
      if (mode%id==mode%paired_id) then
        if (mode%power/=0) then
          j = first(displacement%vectors%id==mode%id, default=0)
          
          if (j==0) then
            ! If the mode is not present in the displacement,
            !    then the displacement along that mode is zero.
            ! As such, the monomial is zero. (0^n=0 if n>0).
            output = 0
            return
          else
            output = output &
                 & * mode%energy_real(displacement=displacement%vectors(j))
          endif
        endif
      else
        if (mode%power/=0 .or. mode%paired_power/=0) then
          j = first(displacement%vectors%id==mode%id, default=0)
          k = first(displacement%vectors%id==mode%paired_id, default=0)
          
          if (j==0 .and. k==0) then
            ! If both modes are not present in the displacement,
            !    then the displacement along that mode is zero.
            ! As such, the monomial is zero. (0^n=0 if n>0).
            output = 0
            return
          elseif (j==0) then
            output = output * mode%energy_real(              &
               & paired_displacement=displacement%vectors(k) )
          elseif (k==0) then
            output = output * mode%energy_real(       &
               & displacement=displacement%vectors(j) )
          else
            output = output * mode%energy_real(              &
               & displacement=displacement%vectors(j),       &
               & paired_displacement=displacement%vectors(k) )
          endif
        endif
      endif
    end associate
  enddo
end function

impure elemental function energy_ComplexModeDisplacement_ComplexMonomial( &
   & this,displacement) result(output)
  implicit none
  
  class(ComplexMonomial),        intent(in) :: this
  type(ComplexModeDisplacement), intent(in) :: displacement
  complex(dp)                               :: output
  
  integer :: i,j,k
  
  output = this%coefficient
  
  ! Loop over the modes in this monomial,
  !    find the mode (and paired mode if relevant) in the displacement
  !    corresponding to each mode,
  !    and calculate the energy term by term.
  do i=1,size(this)
    associate (mode => this%modes_(i))
      if (mode%id==mode%paired_id) then
        if (mode%power/=0) then
          j = first(displacement%vectors%id==mode%id, default=0)
          
          if (j==0) then
            ! If the mode is not present in the displacement,
            !    then the displacement along that mode is zero.
            ! As such, the monomial is zero. (0^n=0 if n>0).
            output = 0
            return
          else
            output = output * mode%energy_complex(    &
               & displacement=displacement%vectors(j) )
          endif
        endif
      else
        if (mode%power/=0 .or. mode%paired_power/=0) then
          j = first(displacement%vectors%id==mode%id, default=0)
          k = first(displacement%vectors%id==mode%paired_id, default=0)
          
          if (      (j==0 .and. mode%power/=0)        &
             & .or. (k==0 .and. mode%paired_power/=0) ) then
            ! If either mode is not present in the displacement,
            !    then the displacement along that mode is zero.
            ! As such, the monomial is zero. (0^n=0 if n>0).
            output = 0
            return
          elseif (j==0) then
            output = output * mode%energy_complex(           &
               & paired_displacement=displacement%vectors(k) )
          elseif (k==0) then
            output = output * mode%energy_complex(    &
               & displacement=displacement%vectors(j) )
          else
            output = output * mode%energy_complex(           &
               & displacement=displacement%vectors(j),       &
               & paired_displacement=displacement%vectors(k) )
          endif
        endif
      endif
    end associate
  enddo
end function

impure elemental function energy_RealModeDisplacement_ComplexPolynomial( &
   & this,displacement) result(output)
  implicit none
  
  class(ComplexPolynomial),    intent(in) :: this
  class(RealModeDisplacement), intent(in) :: displacement
  complex(dp)                             :: output
  
  output = sum(this%terms%energy(displacement))
end function

impure elemental function energy_ComplexModeDisplacement_ComplexPolynomial( &
   & this,displacement) result(output)
  implicit none
  
  class(ComplexPolynomial),       intent(in) :: this
  class(ComplexModeDisplacement), intent(in) :: displacement
  complex(dp)                                :: output
  
  output = sum(this%terms%energy(displacement))
end function

! Evaluate the contribution to the force from
!    a univariate, monomial or polynomial at a given displacement.
! -d/d{u_i} ({u_i}^n) evaluated at u_i=U is -n*U^{n-1}
function force_RealSingleDisplacement_ComplexUnivariate(this,displacement, &
   & paired_displacement) result(output)
  implicit none
  
  class(ComplexUnivariate),     intent(in)           :: this
  type(RealSingleDisplacement), intent(in), optional :: displacement
  type(RealSingleDisplacement), intent(in), optional :: paired_displacement
  complex(dp), allocatable                           :: output(:)
  
  complex(dp) :: magnitude
  complex(dp) :: paired_magnitude
  
  complex(dp) :: values(2)
  complex(dp) :: derivatives(2)
  
  if (.not. (present(displacement) .or. present(paired_displacement))) then
    call print_line(CODE_ERROR//': Neither displacement passed.')
    call err()
  elseif (present(paired_displacement) .and. this%id==this%paired_id) then
    call print_line(CODE_ERROR//': paired displacement passed to &
       &univariate with id=paired_id.')
    call err()
  endif
  
  if (this%id==this%paired_id) then
    magnitude = displacement%magnitude
  else
    ! Convert from real to complex co-ordinates.
    ! x_+ = (x_c + ix_s) / sqrt(2)
    ! x_- = (x_c - ix_s) / sqrt(2)
    magnitude = 0
    paired_magnitude = 0
    if (present(displacement)) then
      ! x_+ : x_c / sqrt(2)
      magnitude = magnitude &
              & + displacement%magnitude &
              & / sqrt(2.0_dp)
      ! x_- : x_c / sqrt(2)
      paired_magnitude = paired_magnitude &
                     & + displacement%magnitude &
                     & / sqrt(2.0_dp)
    endif
    if (present(paired_displacement)) then
      ! x_+ : ix_s / sqrt(2)
      magnitude = magnitude &
              & + paired_displacement%magnitude &
              & * cmplx(0,1,dp) / sqrt(2.0_dp)
      ! x_- : - ix_s / sqrt(2)
      paired_magnitude = paired_magnitude &
                     & - paired_displacement%magnitude &
                     & * cmplx(0,1,dp) / sqrt(2.0_dp)
    endif
  endif
  
  if (this%power>1) then
    values(1) = magnitude**this%power
    derivatives(1) = this%power * magnitude**(this%power-1)
  elseif (this%power==1) then
    values(1) = magnitude
    derivatives(1) = 1
  else
    values(1) = 1
    derivatives(1) = 0
  endif
  
  if (this%id==this%paired_id) then
    output = [-derivatives(1)]
    return
  endif
  
  if (this%paired_power>1) then
    values(2) = paired_magnitude**this%paired_power
    derivatives(2) = this%paired_power &
                 & * paired_magnitude**(this%paired_power-1)
  elseif (this%paired_power==1) then
    values(2) = paired_magnitude
    derivatives(2) = 1
  else
    values(2) = 1
    derivatives(2) = 0
  endif
  
  ! Construct the output in complex co-ordinates...
  output = [-derivatives(1)*values(2), -derivatives(2)*values(1)]
  ! ... then convert back to real co-ordinates.
  ! f_c =   (f_+ + f_-) / sqrt(2)
  ! f_s = - (f_+ - f_-) / sqrt(2)i
  output = [ ( output(1)+output(2))/sqrt(2.0_dp),            &
           & (-output(1)+output(2))/cmplx(0,sqrt(2.0_dp),dp) ]
end function

function force_ComplexSingleDisplacement_ComplexUnivariate(this, &
   & displacement,paired_displacement) result(output)
  implicit none
  
  class(ComplexUnivariate),        intent(in)           :: this
  type(ComplexSingleDisplacement), intent(in), optional :: displacement
  type(ComplexSingleDisplacement), intent(in), optional :: paired_displacement
  complex(dp), allocatable                              :: output(:)
  
  complex(dp) :: magnitude
  complex(dp) :: paired_magnitude
  
  complex(dp) :: values(2)
  complex(dp) :: derivatives(2)
  
  if (.not. (present(displacement) .or. present(paired_displacement))) then
    call print_line(CODE_ERROR//': Neither displacement passed.')
    call err()
  elseif (present(paired_displacement) .and. this%id==this%paired_id) then
    call print_line(CODE_ERROR//': paired displacement passed to &
       &univariate with id=paired_id.')
    call err()
  endif
  
  if (present(displacement)) then
    magnitude = displacement%magnitude
  else
    magnitude = 0
  endif
  
  if (this%power>1) then
    values(1) = magnitude**this%power
    derivatives(1) = this%power * magnitude**(this%power-1)
  elseif (this%power==1) then
    values(1) = magnitude
    derivatives(1) = 1
  else
    values(1) = 1
    derivatives(1) = 0
  endif
  
  if (this%id==this%paired_id) then
    output = [-derivatives(1)]
    return
  endif
  
  if (present(paired_displacement)) then
    paired_magnitude = paired_displacement%magnitude
  else
    paired_magnitude = 0
  endif
  
  if (this%paired_power>1) then
    values(2) = paired_magnitude**this%paired_power
    derivatives(2) = this%paired_power &
                 & * paired_magnitude**(this%paired_power-1)
  elseif (this%paired_power==1) then
    values(2) = paired_magnitude
    derivatives(2) = 1
  else
    values(2) = 1
    derivatives(2) = 0
  endif
  
  output = [-derivatives(1)*values(2), -derivatives(2)*values(1)]
end function

! -d/d{u_i} (c*prod_j[ {u_j}^{n_j} ]) evaluated at {u_i=U_i} is
!    -c*prod_{j/=i}[ {U_j}^{n_j} ] * n_i * {U_i}^{n_i-1}
impure elemental function force_RealModeDisplacement_ComplexMonomial(this, &
   & displacement) result(output)
  implicit none
  
  class(ComplexMonomial),     intent(in) :: this
  type(RealModeDisplacement), intent(in) :: displacement
  type(RealModeForce)                    :: output
  
  logical, allocatable :: value_is_zero(:)
  logical, allocatable :: force_is_zero(:)
  
  type(RealSingleDisplacement), allocatable :: displacements(:)
  type(RealSingleDisplacement), allocatable :: paired_displacements(:)
  
  complex(dp)              :: energy
  integer,     allocatable :: ids(:)
  complex(dp), allocatable :: forces(:)
  
  integer :: i,j,k,ialloc
  
  integer :: no_modes
  
  ! Match each mode in the monomial with the matching mode(s)
  !    in the displacement.
  allocate( displacements(size(this)),        &
          & paired_displacements(size(this)), &
          & value_is_zero(size(this)),        &
          & force_is_zero(size(this)),        &
          & stat=ialloc); call err(ialloc)
  no_modes = 0
  do i=1,size(this)
    associate (mode => this%modes_(i))
      j = first(displacement%vectors%id==mode%id, default=0)
      if (j/=0) then
        displacements(i) = displacement%vectors(j)
      else
        displacements(i) = RealSingleDisplacement(mode%id, 0.0_dp)
      endif
      if (mode%id==mode%paired_id) then
        value_is_zero(i) = j==0 .and. mode%power/=0
        force_is_zero(i) = j==0 .and. mode%power/=1
        if (.not. force_is_zero(i)) then
          no_modes = no_modes+1
        endif
      else
        k = first(displacement%vectors%id==mode%paired_id, default=0)
        if (k/=0) then
          paired_displacements(i) = displacement%vectors(k)
        else
          paired_displacements(i) = RealSingleDisplacement( mode%paired_id, &
                                                          & 0.0_dp          )
        endif
        value_is_zero(i) = (j==0 .and. k==0) &
                   & .and. (mode%power/=0 .or. mode%paired_power/=0)
        force_is_zero(i) = (j==0 .and. mode%power/=1) &
                   & .and. (k==0 .and. mode%paired_power/=1)
        if (.not. force_is_zero(i)) then
          no_modes = no_modes+2
        endif
      endif
    end associate
  enddo
  
  if (count(value_is_zero)>1) then
    ! If U_i is zero for more than one i, then
    !    prod_{j/=i}[ {U_j}^{n_j} ] is always zero,
    !    so all derivatives are zero.
    allocate( ids(0),    &
            & forces(0), &
            & stat=ialloc); call err(ialloc)
  else
    i = first(value_is_zero, default=0)
    if (i/=0) then
      ! Only one U_i is zero.
      if (force_is_zero(i)) then
        ! The derivative is also zero along mode i.
        allocate( ids(0),    &
                & forces(0), &
                & stat=ialloc); call err(ialloc)
      else
        ! The derivative is not zero along mode i.
        associate (mode => this%modes_(i))
        if (mode%id==mode%paired_id) then
          ids = [mode%id]
          forces = mode%force_real(displacements(i))
        else
          ids = [mode%id, mode%paired_id]
          forces = mode%force_real(displacements(i), paired_displacements(i))
        endif
        end associate
        do j=1,size(this)
          if (j/=i) then
            associate (mode => this%modes_(j))
              if (mode%id==mode%paired_id) then
                forces = forces &
                     & * mode%energy_real(displacements(j))
              else
                forces = forces                                    &
                     & * mode%energy_real( displacements(j),       &
                     &                     paired_displacements(j) )
              endif
            end associate
          endif
        enddo
      endif
    else
      ! The force is non-zero along multiple modes.
      allocate( ids(no_modes),    &
              & forces(no_modes), &
              & stat=ialloc); call err(ialloc)
      forces = 1.0_dp
      j = 0
      do i=1,size(this)
        associate (mode => this%modes_(i))
          if (mode%id==mode%paired_id) then
            energy = mode%energy_real(displacements(i))
          else
            energy = mode%energy_real( displacements(i),       &
                                     & paired_displacements(i) )
          endif
          
          if (force_is_zero(i)) then
            forces = forces * energy
          else
            if (mode%id==mode%paired_id) then
              j = j+1
              ids(j) = mode%id
              forces(j:j) = forces(j:j) * mode%force_real(displacements(i))
              forces(:j-1) = forces(:j-1) * energy
              forces(j+1:) = forces(j+1:) * energy
            else
              j = j+2
              ids(j-1) = mode%id
              ids(j) = mode%paired_id
              forces(j-1:j) = mode%force_real( displacements(i),       &
                                             & paired_displacements(i) )
              forces(:j-2) = forces(:j-2) * energy
              forces(j+1:) = forces(j+1:) * energy
            endif
          endif
        end associate
      enddo
    endif
  endif
  
  ! Construct output from components along each mode.
  output = RealModeForce(RealSingleForce(ids, real(forces*this%coefficient)))
end function

impure elemental function force_ComplexModeDisplacement_ComplexMonomial( &
   & this,displacement) result(output)
  implicit none
  
  class(ComplexMonomial),        intent(in) :: this
  type(ComplexModeDisplacement), intent(in) :: displacement
  type(ComplexModeForce)                    :: output
  
  logical, allocatable :: value_is_zero(:)
  logical, allocatable :: force_is_zero(:)
  
  type(ComplexSingleDisplacement), allocatable :: displacements(:)
  type(ComplexSingleDisplacement), allocatable :: paired_displacements(:)
  
  complex(dp)              :: energy
  integer,     allocatable :: ids(:)
  complex(dp), allocatable :: forces(:)
  
  integer :: i,j,k,ialloc
  
  integer :: no_modes
  
  ! Match each mode in the monomial with the matching mode(s)
  !    in the displacement.
  allocate( displacements(size(this)),        &
          & paired_displacements(size(this)), &
          & value_is_zero(size(this)),        &
          & force_is_zero(size(this)),        &
          & stat=ialloc); call err(ialloc)
  no_modes = 0
  do i=1,size(this)
    associate (mode => this%modes_(i))
      j = first(displacement%vectors%id==mode%id, default=0)
      if (j/=0) then
        displacements(i) = displacement%vectors(j)
      else
        displacements(i) = ComplexSingleDisplacement( &
                            & mode%id,                &
                            & cmplx(0.0_dp,0.0_dp,dp) )
      endif
      if (mode%id==mode%paired_id) then
        value_is_zero(i) = j==0 .and. mode%power/=0
        force_is_zero(i) = j==0 .and. mode%power/=1
        if (.not. force_is_zero(i)) then
          no_modes = no_modes+1
        endif
      else
        k = first(displacement%vectors%id==mode%paired_id, default=0)
        if (k/=0) then
          paired_displacements(i) = displacement%vectors(k)
        else
          paired_displacements(i) = ComplexSingleDisplacement( &
                                     & mode%paired_id,         &
                                     & cmplx(0.0_dp,0.0_dp,dp) )
        endif
        value_is_zero(i) = (j==0 .and. mode%power/=0) &
                    & .or. (k==0 .and. mode%paired_power/=0)
        force_is_zero(i) = (j==0 .and. mode%power/=1) &
                   & .and. (k==0 .and. mode%paired_power/=1)
        if (.not. force_is_zero(i)) then
          no_modes = no_modes+2
        endif
      endif
    end associate
  enddo
  
  if (count(value_is_zero)>1) then
    ! If U_i is zero for more than one i, then
    !    prod_{j/=i}[ {U_j}^{n_j} ] is always zero,
    !    so all derivatives are zero.
    allocate( ids(0),    &
            & forces(0), &
            & stat=ialloc); call err(ialloc)
  else
    i = first(value_is_zero, default=0)
    if (i/=0) then
      ! Only one U_i is zero.
      if (force_is_zero(i)) then
        ! The derivative is also zero along mode i.
        allocate( ids(0),    &
                & forces(0), &
                & stat=ialloc); call err(ialloc)
      else
        ! The derivative is not zero along mode i.
        associate (mode => this%modes_(i))
        if (mode%id==mode%paired_id) then
          ids = [mode%id]
          forces = mode%force_complex(displacements(i))
        else
          ids = [mode%id, mode%paired_id]
          forces = mode%force_complex( displacements(i),       &
                                     & paired_displacements(i) )
        endif
        end associate
        do j=1,size(this)
          if (j/=i) then
            associate (mode => this%modes_(j))
              if (mode%id==mode%paired_id) then
                forces = forces &
                     & * mode%energy_complex(displacements(j))
              else
                forces = forces &
                     & * mode%energy_complex( displacements(j),       &
                     &                        paired_displacements(j) )
              endif
            end associate
          endif
        enddo
      endif
    else
      ! The force is non-zero along multiple modes.
      allocate( ids(no_modes),    &
              & forces(no_modes), &
              & stat=ialloc); call err(ialloc)
      forces = 1.0_dp
      j = 0
      do i=1,size(this)
        associate (mode => this%modes_(i))
          if (mode%id==mode%paired_id) then
            energy = mode%energy_complex(displacements(i))
          else
            energy = mode%energy_complex( displacements(i),       &
                                        & paired_displacements(i) )
          endif
          
          if (force_is_zero(i)) then
            forces = forces * energy
          else
            if (mode%id==mode%paired_id) then
              j = j+1
              ids(j) = mode%id
              forces(j:j) = forces(j:j) * mode%force_complex(displacements(i))
              forces(:j-1) = forces(:j-1) * energy
              forces(j+1:) = forces(j+1:) * energy
            else
              j = j+2
              ids(j-1) = mode%id
              ids(j) = mode%paired_id
              forces(j-1:j) = mode%force_complex( displacements(i),       &
                                                & paired_displacements(i) )
              forces(:j-2) = forces(:j-2) * energy
              forces(j+1:) = forces(j+1:) * energy
            endif
          endif
        end associate
      enddo
    endif
  endif
  
  ! Construct output from components along each mode.
  output = ComplexModeForce(ComplexSingleForce(ids, forces*this%coefficient))
end function

impure elemental function force_RealModeDisplacement_ComplexPolynomial(this, &
   & displacement) result(output)
  implicit none
  
  class(ComplexPolynomial),   intent(in) :: this
  type(RealModeDisplacement), intent(in) :: displacement
  type(RealModeForce)                    :: output
  
  type(RealSingleForce) :: zero_force(0)
  
  if (size(this)==0) then
    output = RealModeForce(zero_force)
  else
    output = sum(this%terms%force(displacement))
  endif
end function

impure elemental function force_ComplexModeDisplacement_ComplexPolynomial( &
   & this,displacement) result(output)
  implicit none
  
  class(ComplexPolynomial),      intent(in) :: this
  type(ComplexModeDisplacement), intent(in) :: displacement
  type(ComplexModeForce)                    :: output
  
  type(ComplexSingleForce) :: zero_force(0)
  
  if (size(this)==0) then
    output = ComplexModeForce(zero_force)
  else
    output = sum(this%terms%force(displacement))
  endif
end function

! The harmonic expectation of a monomial or polynomial.
impure elemental function harmonic_expectation_ComplexPolynomial(this, &
   & frequency,thermal_energy,supercell_size) result(output)
  implicit none
  
  class(ComplexPolynomial), intent(in) :: this
  real(dp),                 intent(in) :: frequency
  real(dp),                 intent(in) :: thermal_energy
  integer,                  intent(in) :: supercell_size
  real(dp)                             :: output
  
  output = sum(this%terms%harmonic_expectation( frequency,      &
                                              & thermal_energy, &
                                              & supercell_size  ))
end function

impure elemental function harmonic_expectation_ComplexMonomial(this, &
   & frequency,thermal_energy,supercell_size) result(output)
  implicit none
  
  class(ComplexMonomial), intent(in) :: this
  real(dp),               intent(in) :: frequency
  real(dp),               intent(in) :: thermal_energy
  integer,                intent(in) :: supercell_size
  real(dp)                           :: output
  
  output = real( this%coefficient                                             &
       &       * product(this%modes_%harmonic_expectation( frequency,         &
       &                                                   thermal_energy,    &
       &                                                   supercell_size  )) )
end function

impure elemental function harmonic_expectation_ComplexUnivariate(this, &
   & frequency,thermal_energy,supercell_size) result(output)
  implicit none
  
  class(ComplexUnivariate), intent(in) :: this
  real(dp),                 intent(in) :: frequency
  real(dp),                 intent(in) :: thermal_energy
  integer,                  intent(in) :: supercell_size
  real(dp)                             :: output
  
  integer :: n
  
  real(dp) :: factor
  
  ! If the mode is its own conjugate, then
  !            f(n) (1+e^(-w/T))^n
  ! <(u)^2n> = -------------------
  !            (2Nw(1-e^(-w/T)))^n
  ! where f(n) is the odd factorial.
  
  ! If the mode is not its own conjugate, then
  !               n! (1+e^(-w/T))^n
  ! <(u+u-)^n> = -------------------
  !              (2Nw(1-e^(-w/T)))^n
  
  ! Calculate n, and return if the answer is 0 or 1.
  if (this%id==this%paired_id) then
    if (modulo(this%power,2)==1) then
      output = 0
      return
    endif
    n = this%power/2
  else
    if (this%power/=this%paired_power) then
      output = 0
      return
    endif
    n = this%power
  endif
  
  ! <X^0> = <1> = 1.
  if (n==0) then
    output = 1
    return
  endif
  
  ! Calculate (1+e^(-w/T))/(1-e^(-w/T)).
  if (frequency < 1e-10_dp*thermal_energy) then
    ! High temperature regime (N.B. this is unbounded for w=0).
    ! If w/T<1e-10, (w/T)^2 < 1e-20, smaller than floating point error.
    !
    ! e^(-w/T) = 1 - w/T + 0.5(w/T)^2 + ...
    !
    ! (1+e^(-w/T))/(1-e^(-w/T)) = (2 - w/T + ...) / (w/T - 0.5(w/T)^2 + ...)
    !                           = (2T/w - 1 + ...) / (1 - 0.5w/T + ...)
    !                           = (2T/w - 1 + ...) * (1 + 0.5w/T + ...)
    !                           = 2T/w + 0 + ...
    factor = 2*thermal_energy/frequency
  elseif (frequency < 23*thermal_energy) then
    ! Normal temperature regime.
    factor = (1+exp(-frequency/thermal_energy)) &
         & / (1-exp(-frequency/thermal_energy))
  elseif (frequency < 690*thermal_energy) then
    ! Low temperature regime.
    ! If w/T>23, exp(-w/T)<1e-9, so exp(-2w/T) < 1e-20.
    ! 
    ! (1+e^(-w/T))/(1-e^(-w/T)) = (1 + e^(-w/T)) * (1 + e^(-w/T) + ...)
    !                           = 1 + 2e^(-w/T) + ...
    factor = 1+2*exp(-frequency/thermal_energy)
  else
    ! Very low temperature regime.
    ! If w/T>690, exp(-w/T)<1e-300, smaller than floating point can store.
    factor = 1
  endif
  
  ! Calculate the output.
  if (this%id==this%paired_id) then
    output = odd_factorial(n) &
         & * (factor/(2.0_dp*supercell_size*frequency))**n
  else
    output = factorial(n) &
         & * (factor/(2.0_dp*supercell_size*frequency))**n
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

! The negative of a monomial or polynomial.
impure elemental function negative_ComplexMonomial(this) result(output)
  implicit none
  
  class(ComplexMonomial), intent(in) :: this
  type(ComplexMonomial)              :: output
  
  output = ComplexMonomial( modes       = this%modes_,      &
                          & coefficient = -this%coefficient )
end function

impure elemental function negative_ComplexPolynomial(this) result(output)
  implicit none
  
  class(ComplexPolynomial), intent(in) :: this
  type(ComplexPolynomial)              :: output
  
  output = ComplexPolynomial(-this%terms)
end function

! Subtraction between polynomials and polynomial-like types.
impure elemental function                                            &
   & subtract_ComplexPolynomialable_ComplexPolynomialable(this,that) &
   & result(output)
  implicit none
  
  class(ComplexPolynomialable), intent(in) :: this
  class(ComplexPolynomialable), intent(in) :: that
  type(ComplexPolynomial)                  :: output
  
  output = this + (-that%to_ComplexPolynomial())
end function

! Sum polynomial-like types.
function sum_ComplexPolynomialables(input) result(output)
  implicit none
  
  class(ComplexPolynomialable), intent(in) :: input(:)
  type(ComplexPolynomial)                  :: output
  
  type(ComplexMonomial) :: zero_monomial(0)
  
  integer :: i
  
  if (size(input)==0) then
    output = ComplexPolynomial(zero_monomial)
  else
    output = input(1)%to_ComplexPolynomial()
    do i=2,size(input)
      output = output + input(i)
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
  
  integer :: i,ialloc
  
  allocate(output(0), stat=ialloc); call err(ialloc)
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
  
  integer :: i,ialloc
  
  allocate(output(0), stat=ialloc); call err(ialloc)
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
  
  integer :: i,ialloc
  
  allocate(output(0), stat=ialloc); call err(ialloc)
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
    
    allocate(modes(0), stat=ialloc); call err(ialloc)
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
