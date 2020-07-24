! ======================================================================
! Polynomials in the complex representation, but with monomials stored as
!    (m+m*) and i(m-m*) pairs with real coefficients.
! ======================================================================
module paired_polynomial_module
  use utils_module
  use structure_module
  
  use real_mode_module
  use real_single_mode_displacement_module
  use real_single_mode_force_module
  use real_mode_displacement_module
  use real_mode_force_module
  use complex_polynomial_module
  implicit none
  
  private
  
  public :: PairedMonomial
  public :: PairedPolynomial
  public :: size
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: operator(-)
  public :: sum
  public :: conjg
  public :: compare_paired_monomials
  public :: matching_pair
  
  ! --------------------------------------------------
  ! Types, and conversions between types.
  ! --------------------------------------------------
  ! It is desirable to be able to convert monomials to polynomials
  ! This is acheived by extending PairedMonomial from PairedPolynomialable.
  
  type, abstract, extends(Stringable) :: PairedPolynomialable
  contains
    procedure(to_PairedPolynomial_PairedPolynomialable), deferred, public &
       & :: to_PairedPolynomial
  end type
  
  type, extends(PairedPolynomialable) :: PairedMonomial
    real(dp) :: coefficient
    ! Modes is private so that it can be guaranteed to be sorted.
    type(ComplexUnivariate), allocatable, private :: modes_(:)
    ! pair_ = : 0 if m
    !           1 if m+m*
    !          -1 if i(u-u*).
    integer, private :: pair_
  contains
    procedure, public :: to_PairedPolynomial => &
                       & to_PairedPolynomial_PairedMonomial
    
    procedure, public :: mode => mode_PairedMonomial
    procedure, public :: modes => modes_PairedMonomial
    procedure, public :: id => id_PairedMonomial
    procedure, public :: ids => ids_PairedMonomial
    procedure, public :: paired_id  => paired_id_PairedMonomial
    procedure, public :: paired_ids => paired_ids_PairedMonomial
    procedure, public :: power  => power_PairedMonomial
    procedure, public :: powers => powers_PairedMonomial
    procedure, public :: paired_power  => paired_power_PairedMonomial
    procedure, public :: paired_powers => paired_powers_PairedMonomial
    
    procedure, public :: set_modes => set_modes_PairedMonomial
    
    procedure, public :: simplify => simplify_PairedMonomial
    
    procedure, public :: total_power => total_power_PairedMonomial
    
    generic,   public  :: energy => &
                        & energy_RealModeDisplacement_PairedMonomial
    procedure, private :: energy_RealModeDisplacement_PairedMonomial
    
    generic,   public  :: force => &
                        & force_RealModeDisplacement_PairedMonomial
    procedure, private :: force_RealModeDisplacement_PairedMonomial
    
    procedure, public :: harmonic_expectation => &
                       & harmonic_expectation_PairedMonomial
    
    procedure, public :: read  => read_PairedMonomial
    procedure, public :: write => write_PairedMonomial
  end type
  
  interface PairedMonomial
    module procedure new_PairedMonomial
    module procedure new_PairedMonomial_String
  end interface
  
  type, extends(PairedPolynomialable) :: PairedPolynomial
    type(PairedMonomial), allocatable :: terms(:)
  contains
    procedure, public :: to_PairedPolynomial => &
                       & to_PairedPolynomial_PairedPolynomial
    
    procedure, public :: simplify => simplify_PairedPolynomial
    
    generic,   public  :: energy => &
                        & energy_RealModeDisplacement_PairedPolynomial
    procedure, private :: energy_RealModeDisplacement_PairedPolynomial
    
    generic,   public  :: force => &
                        & force_RealModeDisplacement_PairedPolynomial
    procedure, private :: force_RealModeDisplacement_PairedPolynomial
    
    procedure, public :: harmonic_expectation => &
                       & harmonic_expectation_PairedPolynomial
    
    procedure, public :: read  => read_PairedPolynomial
    procedure, public :: write => write_PairedPolynomial
  end type
  
  interface PairedPolynomial
    module procedure new_PairedPolynomial
    module procedure new_PairedPolynomial_PairedPolynomialable
    module procedure new_PairedPolynomial_String
  end interface
  
  abstract interface
    function to_PairedPolynomial_PairedPolynomialable(this) result(output)
      import PairedPolynomial
      import PairedPolynomialable
      implicit none
      
      class(PairedPolynomialable), intent(in) :: this
      type(PairedPolynomial)                  :: output
    end function
  end interface
  
  ! --------------------------------------------------
  ! Operations involving types.
  ! --------------------------------------------------
  interface size
    module procedure size_PairedMonomial
    module procedure size_PairedPolynomial
  end interface
  
  interface conjg
    module procedure conjg_PairedMonomial
    module procedure conjg_PairedPolynomial
  end interface
  
  interface operator(*)
    module procedure multiply_PairedMonomial_real
    module procedure multiply_real_PairedMonomial
    
    module procedure multiply_PairedPolynomial_real
    module procedure multiply_real_PairedPolynomial
  end interface
  
  interface operator(/)
    module procedure divide_PairedMonomial_real
    
    module procedure divide_PairedPolynomial_real
  end interface
  
  interface operator(+)
    module procedure add_PairedPolynomialable_PairedPolynomialable
  end interface
  
  interface operator(-)
    module procedure negative_PairedMonomial
    module procedure negative_PairedPolynomial
    
    module procedure subtract_PairedPolynomialable_PairedPolynomialable
  end interface
  
  interface sum
    module procedure sum_PairedPolynomialables
  end interface
contains

! ----------------------------------------------------------------------
! Constructors.
! ----------------------------------------------------------------------
function new_PairedMonomial(coefficient,modes,pair) result(this)
  implicit none
  
  real(dp),                intent(in) :: coefficient
  type(ComplexUnivariate), intent(in) :: modes(:)
  integer,                 intent(in) :: pair
  type(PairedMonomial)                :: this
  
  this%coefficient = coefficient
  this%modes_      = modes(sort(modes%id))
  this%pair_       = pair
end function

function new_PairedPolynomial(terms) result(this)
  implicit none
  
  type(PairedMonomial), intent(in) :: terms(:)
  type(PairedPolynomial)           :: this
  
  this%terms = terms
end function

function new_PairedPolynomial_PairedPolynomialable(input) result(this)
  implicit none
  
  class(PairedPolynomialable), intent(in) :: input
  type(PairedPolynomial)                  :: this
  
  this = input%to_PairedPolynomial()
end function

! ----------------------------------------------------------------------
! Conversions between types.
! ----------------------------------------------------------------------
function to_PairedPolynomial_PairedMonomial(this) result(output)
  implicit none
  
  class(PairedMonomial), intent(in) :: this
  type(PairedPolynomial)            :: output
  
  select type(this); type is(PairedMonomial)
    output = PairedPolynomial([this])
  class default
    call err()
  end select
end function

function to_PairedPolynomial_PairedPolynomial(this) result(output)
  implicit none
  
  class(PairedPolynomial), intent(in) :: this
  type(PairedPolynomial)              :: output
  
  select type(this); type is(PairedPolynomial)
    output = this
  class default
    call err()
  end select
end function

! ----------------------------------------------------------------------
! Operations involving types.
! ----------------------------------------------------------------------

! The number of modes in a monomial, or terms in a polynomial.
function size_PairedMonomial(this) result(output)
  implicit none
  
  type(PairedMonomial), intent(in) :: this
  integer                          :: output
  
  output = size(this%modes_)
end function

function size_PairedPolynomial(this) result(output)
  implicit none
  
  class(PairedPolynomial), intent(in) :: this
  integer                             :: output
  
  output = size(this%terms)
end function

! Getters for monomials.
impure elemental function mode_PairedMonomial(this,index) result(output)
  implicit none
  
  class(PairedMonomial), intent(in) :: this
  integer,               intent(in) :: index
  type(ComplexUnivariate)           :: output
  
  output = this%modes_(index)
end function

! Getter for modes. Has several run types:
!    - If no arguments specified, returns all modes.
!    - If indices specified, returns modes at specified indices.
!    - If ids and paired_ids specified, returns modes with specified ids.
!    - If exclude_ids specicified, returns all modes but those excluded.
function modes_PairedMonomial(this,indices,ids,paired_ids,exclude_ids) &
   & result(output)
  implicit none
  
  class(PairedMonomial), intent(in)           :: this
  integer,               intent(in), optional :: indices(:)
  integer,               intent(in), optional :: ids(:)
  integer,               intent(in), optional :: paired_ids(:)
  integer,               intent(in), optional :: exclude_ids(:)
  type(ComplexUnivariate), allocatable        :: output(:)
  
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

impure elemental function id_PairedMonomial(this,index) result(output)
  implicit none
  
  class(PairedMonomial), intent(in) :: this
  integer,               intent(in) :: index
  integer                           :: output
  
  output = this%modes_(index)%id
end function

function ids_PairedMonomial(this,indices) result(output)
  implicit none
  
  class(PairedMonomial), intent(in)           :: this
  integer,               intent(in), optional :: indices(:)
  integer, allocatable                        :: output(:)
  
  if (present(indices)) then
    output = this%modes_(indices)%id
  else
    output = this%modes_%id
  endif
end function

impure elemental function paired_id_PairedMonomial(this,index) result(output)
  implicit none
  
  class(PairedMonomial), intent(in) :: this
  integer,               intent(in) :: index
  integer                           :: output
  
  output = this%modes_(index)%paired_id
end function

function paired_ids_PairedMonomial(this,indices) result(output)
  implicit none
  
  class(PairedMonomial), intent(in)           :: this
  integer,               intent(in), optional :: indices(:)
  integer, allocatable                        :: output(:)
  
  if (present(indices)) then
    output = this%modes_(indices)%paired_id
  else
    output = this%modes_%paired_id
  endif
end function

impure elemental function power_PairedMonomial(this,index) result(output)
  implicit none
  
  class(PairedMonomial), intent(in) :: this
  integer,               intent(in) :: index
  integer                           :: output
  
  output = this%modes_(index)%power
end function

function powers_PairedMonomial(this,indices) result(output)
  implicit none
  
  class(PairedMonomial), intent(in)           :: this
  integer,               intent(in), optional :: indices(:)
  integer, allocatable                        :: output(:)
  
  if (present(indices)) then
    output = this%modes_(indices)%power
  else
    output = this%modes_%power
  endif
end function

impure elemental function paired_power_PairedMonomial(this,index) &
   & result(output)
  implicit none
  
  class(PairedMonomial), intent(in) :: this
  integer,               intent(in) :: index
  integer                           :: output
  
  output = this%modes_(index)%paired_power
end function

function paired_powers_PairedMonomial(this,indices) result(output)
  implicit none
  
  class(PairedMonomial), intent(in)           :: this
  integer,               intent(in), optional :: indices(:)
  integer, allocatable                        :: output(:)
  
  if (present(indices)) then
    output = this%modes_(indices)%paired_power
  else
    output = this%modes_%paired_power
  endif
end function

! Set the modes of a monomial.
subroutine set_modes_PairedMonomial(this,modes,pair)
  implicit none
  
  class(PairedMonomial),   intent(inout)        :: this
  type(ComplexUnivariate), intent(in)           :: modes(:)
  integer,                 intent(in), optional :: pair
  
  this%modes_ = modes
  if (present(pair)) then
    this%pair_ = pair
  endif
end subroutine

! Simplify a monomial or polynomial.
impure elemental subroutine simplify_PairedMonomial(this)
  implicit none
  
  class(PairedMonomial), intent(inout) :: this
  
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

impure elemental subroutine simplify_PairedPolynomial(this)
  implicit none
  
  class(PairedPolynomial), intent(inout) :: this
  
  integer,  allocatable :: unique_terms(:)
  real(dp), allocatable :: coefficients(:)
  
  integer :: i,j
  
  call this%terms%simplify()
  
  ! Add together monomials with the same powers.
  unique_terms = set(this%terms, compare_paired_monomials)
  coefficients = [(0.0_dp, i=1, size(unique_terms))]
  do i=1,size(this%terms)
    do j=1,size(unique_terms)
      if (compare_paired_monomials( this%terms(i),              &
                                  & this%terms(unique_terms(j)) )) then
        coefficients(j) = coefficients(j) + this%terms(i)%coefficient
      endif
    enddo
  enddo
  
  this%terms = this%terms(unique_terms)
  this%terms%coefficient = coefficients
end subroutine

! Find the conjugate of a monomial or polynomial.
impure elemental function conjg_PairedMonomial(this) result(output)
  implicit none
  
  type(PairedMonomial), intent(in) :: this
  type(PairedMonomial)             :: output
  
  output = this
  
  if (this%pair_==-1) then
    output%coefficient = -output%coefficient
  endif
end function

impure elemental function conjg_PairedPolynomial(this) result(output)
  implicit none
  
  type(PairedPolynomial), intent(in) :: this
  type(PairedPolynomial)             :: output
  
  output = PairedPolynomial(conjg(this%terms))
end function

! The total power of a monomial.
impure elemental function total_power_PairedMonomial(this) result(output)
  implicit none
  
  class(PairedMonomial), intent(in) :: this
  integer                           :: output
  
  output = sum(this%modes_%total_power())
end function

! Evaluate the contribution to the energy from a monomial or polynomial
!    at a given displacement.
impure elemental function energy_RealModeDisplacement_PairedMonomial(this, &
   & displacement) result(output)
  implicit none
  
  class(PairedMonomial),      intent(in) :: this
  type(RealModeDisplacement), intent(in) :: displacement
  real(dp)                               :: output
  
  complex(dp) :: complex_energy
  
  integer :: i,j,k
  
  complex_energy = this%coefficient
  
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
            complex_energy = complex_energy                          &
                         & * mode%energy_real(                       &
                         &      displacement=displacement%vectors(j) )
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
            complex_energy = complex_energy * mode%energy_real( &
                  & paired_displacement=displacement%vectors(k) )
          elseif (k==0) then
            complex_energy = complex_energy * mode%energy_real( &
                         & displacement=displacement%vectors(j) )
          else
            complex_energy = complex_energy * mode%energy_real( &
                  & displacement=displacement%vectors(j),       &
                  & paired_displacement=displacement%vectors(k) )
          endif
        endif
      endif
    end associate
  enddo
  
  if (this%pair_==0) then
    output = real(complex_energy)
  elseif (this%pair_==1) then
    output = 2*real(complex_energy)
  else
    output = -2*aimag(complex_energy)
  endif
end function

impure elemental function energy_RealModeDisplacement_PairedPolynomial( &
   & this,displacement) result(output)
  implicit none
  
  class(PairedPolynomial),     intent(in) :: this
  class(RealModeDisplacement), intent(in) :: displacement
  real(dp)                                :: output
  
  output = sum(this%terms%energy(displacement))
end function

! Evaluate the contribution to the force from
!    a monomial or polynomial at a given displacement.

! -d/d{u_i} (c*prod_j[ {u_j}^{n_j} ]) evaluated at {u_i=U_i} is
!    -c*prod_{j/=i}[ {U_j}^{n_j} ] * n_i * {U_i}^{n_i-1}
impure elemental function force_RealModeDisplacement_PairedMonomial(this, &
   & displacement) result(output)
  implicit none
  
  class(PairedMonomial),      intent(in) :: this
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
  if (this%pair_==0) then
    output = RealModeForce(RealSingleForce( ids, &
                                          & real(forces*this%coefficient)))
  elseif (this%pair_==1) then
    output = RealModeForce(RealSingleForce( ids, &
                                          & 2*real(forces*this%coefficient)))
  else
    output = RealModeForce(RealSingleForce( ids, &
                                          & -2*aimag(forces*this%coefficient)))
  endif
end function

impure elemental function force_RealModeDisplacement_PairedPolynomial(this, &
   & displacement) result(output)
  implicit none
  
  class(PairedPolynomial),    intent(in) :: this
  type(RealModeDisplacement), intent(in) :: displacement
  type(RealModeForce)                    :: output
  
  type(RealSingleForce) :: zero_force(0)
  
  if (size(this)==0) then
    output = RealModeForce(zero_force)
  else
    output = sum(this%terms%force(displacement))
  endif
end function

! The harmonic expectation of a monomial or polynomial.
impure elemental function harmonic_expectation_PairedMonomial(this, &
   & frequency,thermal_energy,supercell_size) result(output)
  implicit none
  
  class(PairedMonomial), intent(in) :: this
  real(dp),              intent(in) :: frequency
  real(dp),              intent(in) :: thermal_energy
  integer,               intent(in) :: supercell_size
  real(dp)                          :: output
  
  if (this%pair_==-1) then
    output = 0
    return
  endif
  
  output = this%coefficient                                          &
       & * product(this%modes_%harmonic_expectation( frequency,      &
       &                                             thermal_energy, &
       &                                             supercell_size  ))
  
  if (this%pair_==1) then
    output = output*2
  endif
end function

impure elemental function harmonic_expectation_PairedPolynomial(this, &
   & frequency,thermal_energy,supercell_size) result(output)
  implicit none
  
  class(PairedPolynomial), intent(in) :: this
  real(dp),                intent(in) :: frequency
  real(dp),                intent(in) :: thermal_energy
  integer,                 intent(in) :: supercell_size
  real(dp)                            :: output
  
  output = sum(this%terms%harmonic_expectation( frequency,      &
                                              & thermal_energy, &
                                              & supercell_size  ))
end function

! Multiplication and division by scalars.
impure elemental function multiply_PairedMonomial_real(this,that) &
   & result(output)
  implicit none
  
  type(PairedMonomial), intent(in) :: this
  real(dp),             intent(in) :: that
  type(PairedMonomial)             :: output
  
  output = PairedMonomial(this%coefficient*that, this%modes_, this%pair_)
end function

impure elemental function multiply_real_PairedMonomial(this,that) &
   & result(output)
  implicit none
  
  real(dp),             intent(in) :: this
  type(PairedMonomial), intent(in) :: that
  type(PairedMonomial)             :: output
  
  output = PairedMonomial(this*that%coefficient, that%modes_, that%pair_)
end function

impure elemental function multiply_PairedPolynomial_real(this,that) &
   & result(output)
  implicit none
  
  type(PairedPolynomial), intent(in) :: this
  real(dp),               intent(in) :: that
  type(PairedPolynomial)             :: output
  
  output = PairedPolynomial(this%terms * that)
end function

impure elemental function multiply_real_PairedPolynomial(this,that) &
   & result(output)
  implicit none
  
  real(dp),               intent(in) :: this
  type(PairedPolynomial), intent(in) :: that
  type(PairedPolynomial)             :: output
  
  output = PairedPolynomial(this * that%terms)
end function

impure elemental function divide_PairedMonomial_real(this,that) &
   & result(output)
  implicit none
  
  type(PairedMonomial), intent(in) :: this
  real(dp),             intent(in) :: that
  type(PairedMonomial)             :: output
  
  output = PairedMonomial(this%coefficient/that, this%modes_, this%pair_)
end function

impure elemental function divide_PairedPolynomial_real(this,that) &
   & result(output)
  implicit none
  
  type(PairedPolynomial), intent(in) :: this
  real(dp),               intent(in) :: that
  type(PairedPolynomial)             :: output
  
  output = PairedPolynomial(this%terms / that)
end function

! Addition between polynomials and polynomial-like types.
impure elemental function add_PairedPolynomialable_PairedPolynomialable( &
   & this,that) result(output)
  implicit none
  
  class(PairedPolynomialable), intent(in) :: this
  class(PairedPolynomialable), intent(in) :: that
  type(PairedPolynomial)                  :: output
  
  type(PairedPolynomial) :: this_polynomial
  type(PairedPolynomial) :: that_polynomial
  
  integer :: no_terms
  
  integer :: i,j,ialloc
  
  this_polynomial = this%to_PairedPolynomial()
  that_polynomial = that%to_PairedPolynomial()
  
  allocate( output%terms(size(this_polynomial)+size(that_polynomial)), &
          & stat=ialloc); call err(ialloc)
  output%terms(:size(this_polynomial)) = this_polynomial%terms
  no_terms = size(this_polynomial)
  do i=1,size(that_polynomial)
    j = first_equivalent( this_polynomial%terms,    &
                        & that_polynomial%terms(i), &
                        & compare_paired_monomials, &
                        & default=0                 )
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
impure elemental function negative_PairedMonomial(this) result(output)
  implicit none
  
  class(PairedMonomial), intent(in) :: this
  type(PairedMonomial)              :: output
  
  output = PairedMonomial( modes       = this%modes_,       &
                         & coefficient = -this%coefficient, &
                         & pair        = this%pair_         )
end function

impure elemental function negative_PairedPolynomial(this) result(output)
  implicit none
  
  class(PairedPolynomial), intent(in) :: this
  type(PairedPolynomial)              :: output
  
  output = PairedPolynomial(-this%terms)
end function

! Subtraction between polynomials and polynomial-like types.
impure elemental function                                            &
   & subtract_PairedPolynomialable_PairedPolynomialable(this,that) &
   & result(output)
  implicit none
  
  class(PairedPolynomialable), intent(in) :: this
  class(PairedPolynomialable), intent(in) :: that
  type(PairedPolynomial)                  :: output
  
  output = this + (-that%to_PairedPolynomial())
end function

! Sum polynomial-like types.
function sum_PairedPolynomialables(input) result(output)
  implicit none
  
  class(PairedPolynomialable), intent(in) :: input(:)
  type(PairedPolynomial)                  :: output
  
  type(PairedMonomial) :: zero_monomial(0)
  
  integer :: i
  
  if (size(input)==0) then
    output = PairedPolynomial(zero_monomial)
  else
    output = input(1)%to_PairedPolynomial()
    do i=2,size(input)
      output = output + input(i)
    enddo
  endif
end function

! ----------------------------------------------------------------------
! Compares two monomials for equality up to coefficient.
! ----------------------------------------------------------------------
function compare_paired_monomials(this,that) result(output)
  implicit none
  
  class(*), intent(in) :: this
  class(*), intent(in) :: that
  logical              :: output
  
  select type(this); type is(PairedMonomial)
    select type(that); type is(PairedMonomial)
      if (size(this)/=size(that)) then
        output = .false.
      elseif (this%pair_/=that%pair_) then
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
! Returns .true. if the arguments have matching modes, and:
!    this%pair_ ==  0 and that%pair_ == 0
! or this%pair_ ==  1 and that%pair_== -1
! or this%pair_ == -1 and that%pair_==  1
! ----------------------------------------------------------------------
function matching_pair(this,that) result(output)
  implicit none
  
  class(*), intent(in) :: this
  class(*), intent(in) :: that
  logical              :: output
  
  select type(this); type is(PairedMonomial)
    select type(that); type is(PairedMonomial)
      if (size(this)/=size(that)) then
        output = .false.
      elseif (this%pair_+that%pair_/=0) then
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
subroutine read_PairedMonomial(this,input)
  implicit none
  
  class(PairedMonomial), intent(out) :: this
  type(String),           intent(in)  :: input
  
  type(String),            allocatable :: line(:)
  real(dp)                             :: coefficient
  type(ComplexUnivariate), allocatable :: modes(:)
  integer                              :: pair
  
  integer :: i,ialloc
  
  select type(this); type is(PairedMonomial)
    ! Splitting the input by '*' separates the coefficient and the modes,
    !    but also splits some modes in two.
    line = split_line(input,delimiter='*')
    
    ! Check for '+cc', and remove it from parsing.
    i = size(line)
    if (slice(line(i),len(line(i))-2,len(line(i)))=='+cc') then
      pair = 1
      line(i) = slice(line(i),1,len(line(i))-3)
    endif
    
    ! Check for an 'i' after the coefficient, and remove it from parsing.
    ! N.B. setting pair=-1 intentionally overwrites the above pair=1.
    if (slice(line(1),len(line(1)),len(line(1)))=='i') then
      pair = -1
      line(1) = slice(line(1),1,len(line(1))-1)
    endif
    
    coefficient = dble(line(1))
    
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
    
    this = PairedMonomial(coefficient, modes, pair)
  end select
end subroutine

function write_PairedMonomial(this) result(output)
  implicit none
  
  class(PairedMonomial), intent(in) :: this
  type(String)                       :: output
  
  select type(this); type is(PairedMonomial)
    if (size(this%modes_)==0) then
      output = ''
    else
      output = '*'//join(this%modes_, delimiter='*')
    endif
    
    if (this%pair_==0) then
      output = this%coefficient//output
    elseif (this%pair_==1) then
      output = this%coefficient//output//'+cc'
    else
      output = this%coefficient//'i'//output//'+cc'
    endif
  end select
end function

impure elemental function new_PairedMonomial_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(PairedMonomial)     :: this
  
  call this%read(input)
end function

subroutine read_PairedPolynomial(this,input)
  implicit none
  
  class(PairedPolynomial), intent(out) :: this
  type(String),            intent(in)  :: input
  
  type(String),         allocatable :: terms(:)
  type(PairedMonomial), allocatable :: monomials(:)
  
  type(String) :: plus
  
  plus = '+'
  select type(this); type is(PairedPolynomial)
    terms = split_line(input)
    terms = terms(filter(terms/=plus))
    
    monomials = PairedMonomial(terms)
    
    this = PairedPolynomial(monomials)
  end select
end subroutine

function write_PairedPolynomial(this) result(output)
  implicit none
  
  class(PairedPolynomial), intent(in) :: this
  type(String)                         :: output
  
  select type(this); type is(PairedPolynomial)
    output = join(this%terms, delimiter=' + ')
  end select
end function

impure elemental function new_PairedPolynomial_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(PairedPolynomial)   :: this
  
  call this%read(input)
end function
end module
