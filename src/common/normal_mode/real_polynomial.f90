! ======================================================================
! The building blocks of basis functions in real co-ordinates.
! ======================================================================
module caesar_real_polynomial_module
  use caesar_utils_module
  
  use caesar_real_mode_module
  use caesar_real_single_mode_displacement_module
  use caesar_real_single_mode_force_module
  use caesar_real_mode_displacement_module
  use caesar_real_mode_force_module
  use caesar_complex_mode_module
  use caesar_complex_single_mode_displacement_module
  use caesar_complex_single_mode_force_module
  use caesar_complex_mode_displacement_module
  use caesar_complex_mode_force_module
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
    
    procedure, private :: energy_real => &
                        & energy_RealSingleDisplacement_RealUnivariate
    procedure, private :: energy_complex => &
                        & energy_ComplexSingleDisplacement_RealUnivariate
    
    procedure, private :: force_real => &
                        & force_RealSingleDisplacement_RealUnivariate
    procedure, private :: force_complex => &
                        & force_ComplexSingleDisplacement_RealUnivariate
    
    procedure, public :: read  => read_RealUnivariate
    procedure, public :: write => write_RealUnivariate
  end type
  
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
    
    generic,   public  :: energy =>                                 &
                        & energy_RealModeDisplacement_RealMonomial, &
                        & energy_ComplexModeDisplacement_RealMonomial
    procedure, private :: energy_RealModeDisplacement_RealMonomial
    procedure, private :: energy_ComplexModeDisplacement_RealMonomial
    
    generic,   public  :: force =>                                 &
                        & force_RealModeDisplacement_RealMonomial, &
                        & force_ComplexModeDisplacement_RealMonomial
    procedure, private :: force_RealModeDisplacement_RealMonomial
    procedure, private :: force_ComplexModeDisplacement_RealMonomial
    
    procedure, public :: read  => read_RealMonomial
    procedure, public :: write => write_RealMonomial
  end type
  
  type, extends(RealPolynomialable) :: RealPolynomial
    type(RealMonomial), allocatable :: terms(:)
  contains
    procedure, public :: to_RealPolynomial => to_RealPolynomial_RealPolynomial
    
    procedure, public :: simplify => simplify_RealPolynomial
    
    generic,   public  :: energy =>                                   &
                        & energy_RealModeDisplacement_RealPolynomial, &
                        & energy_ComplexModeDisplacement_RealPolynomial
    procedure, private :: energy_RealModeDisplacement_RealPolynomial
    procedure, private :: energy_ComplexModeDisplacement_RealPolynomial
    
    generic,   public  :: force =>                                   &
                        & force_RealModeDisplacement_RealPolynomial, &
                        & force_ComplexModeDisplacement_RealPolynomial
    procedure, private :: force_RealModeDisplacement_RealPolynomial
    procedure, private :: force_ComplexModeDisplacement_RealPolynomial
    
    procedure, public :: read  => read_RealPolynomial
    procedure, public :: write => write_RealPolynomial
  end type
  
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
  
  interface RealUnivariate
    ! ----------------------------------------------------------------------
    ! Constructors.
    ! ----------------------------------------------------------------------
    impure elemental module function new_RealUnivariate(id,paired_id,power, &
       & paired_power) result(this) 
      integer, intent(in)  :: id
      integer, intent(in)  :: paired_id
      integer, intent(in)  :: power
      integer, intent(in)  :: paired_power
      type(RealUnivariate) :: this
    end function
  
    module function new_RealUnivariate_RealMode(mode,power,paired_power) &
       & result(this) 
      type(RealMode), intent(in)           :: mode
      integer,        intent(in)           :: power
      integer,        intent(in), optional :: paired_power
      type(RealUnivariate)                 :: this
    end function
  end interface
  
  interface RealMonomial
    module function new_RealMonomial(coefficient,modes) result(this) 
      real(dp),             intent(in) :: coefficient
      type(RealUnivariate), intent(in) :: modes(:)
      type(RealMonomial)               :: this
    end function
  
    module function new_RealMonomial_RealMonomialable(input) result(this) 
      class(RealMonomialable), intent(in) :: input
      type(RealMonomial)                  :: this
    end function
  end interface
  
  interface RealPolynomial
    module function new_RealPolynomial(terms) result(this) 
      type(RealMonomial), intent(in) :: terms(:)
      type(RealPolynomial)           :: this
    end function
  
    module function new_RealPolynomial_RealPolynomialable(input) result(this) 
      class(RealPolynomialable), intent(in) :: input
      type(RealPolynomial)                  :: this
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Conversions between types.
    ! ----------------------------------------------------------------------
    module function to_RealMonomial_RealUnivariate(this) result(output) 
      class(RealUnivariate), intent(in) :: this
      type(RealMonomial)                :: output
    end function
  end interface
  
  interface
    module function to_RealPolynomial_RealUnivariate(this) result(output) 
      class(RealUnivariate), intent(in) :: this
      type(RealPolynomial)              :: output
    end function
  end interface
  
  interface
    module function to_RealMonomial_RealMonomial(this) result(output) 
      class(RealMonomial), intent(in) :: this
      type(RealMonomial)              :: output
    end function
  end interface
  
  interface
    module function to_RealPolynomial_RealMonomial(this) result(output) 
      class(RealMonomial), intent(in) :: this
      type(RealPolynomial)            :: output
    end function
  end interface
  
  interface
    module function to_RealPolynomial_RealPolynomial(this) result(output) 
      class(RealPolynomial), intent(in) :: this
      type(RealPolynomial)              :: output
    end function
  end interface
  
  interface size
    ! ----------------------------------------------------------------------
    ! Operations involving types.
    ! ----------------------------------------------------------------------
    
    ! The number of modes in a monomial, or terms in a polynomial.
    module function size_RealMonomial(this) result(output) 
      type(RealMonomial), intent(in) :: this
      integer                        :: output
    end function
  
    module function size_RealPolynomial(this) result(output) 
      class(RealPolynomial), intent(in) :: this
      integer                           :: output
    end function
  end interface
  
  interface
    ! Getters for monomials.
    impure elemental module function mode_RealMonomial(this,index) &
       & result(output) 
      class(RealMonomial), intent(in) :: this
      integer,             intent(in) :: index
      type(RealUnivariate)            :: output
    end function
  end interface
  
  interface
    module function modes_RealMonomial(this,indices) result(output) 
      class(RealMonomial), intent(in)           :: this
      integer,             intent(in), optional :: indices(:)
      type(RealUnivariate), allocatable         :: output(:)
    end function
  end interface
  
  interface
    impure elemental module function id_RealMonomial(this,index) &
       & result(output) 
      class(RealMonomial), intent(in) :: this
      integer,             intent(in) :: index
      integer                         :: output
    end function
  end interface
  
  interface
    module function ids_RealMonomial(this,indices) result(output) 
      class(RealMonomial), intent(in)           :: this
      integer,             intent(in), optional :: indices(:)
      integer, allocatable                      :: output(:)
    end function
  end interface
  
  interface
    impure elemental module function paired_id_RealMonomial(this,index) &
       & result(output) 
      class(RealMonomial), intent(in) :: this
      integer,             intent(in) :: index
      integer                         :: output
    end function
  end interface
  
  interface
    module function paired_ids_RealMonomial(this,indices) result(output) 
      class(RealMonomial), intent(in)           :: this
      integer,             intent(in), optional :: indices(:)
      integer, allocatable                      :: output(:)
    end function
  end interface
  
  interface
    impure elemental module function power_RealMonomial(this,index) &
       & result(output) 
      class(RealMonomial), intent(in) :: this
      integer,             intent(in) :: index
      integer                         :: output
    end function
  end interface
  
  interface
    module function powers_RealMonomial(this,indices) result(output) 
      class(RealMonomial), intent(in)           :: this
      integer,             intent(in), optional :: indices(:)
      integer, allocatable                      :: output(:)
    end function
  end interface
  
  interface
    impure elemental module function paired_power_RealMonomial(this,index) &
       & result(output) 
      class(RealMonomial), intent(in) :: this
      integer,             intent(in) :: index
      integer                         :: output
    end function
  end interface
  
  interface
    module function paired_powers_RealMonomial(this,indices) result(output) 
      class(RealMonomial), intent(in)           :: this
      integer,             intent(in), optional :: indices(:)
      integer, allocatable                      :: output(:)
    end function
  end interface
  
  interface
    ! Simplify a monomial or polynomial.
    impure elemental module subroutine simplify_RealMonomial(this) 
      class(RealMonomial), intent(inout) :: this
    end subroutine
  end interface
  
  interface
    impure elemental module subroutine simplify_RealPolynomial(this) 
      class(RealPolynomial), intent(inout) :: this
    end subroutine
  end interface
  
  interface
    ! The total power of a univariate or monomial.
    impure elemental module function total_power_RealUnivariate(this) &
       & result(output) 
      class(RealUnivariate), intent(in) :: this
      integer                           :: output
    end function
  end interface
  
  interface
    impure elemental module function total_power_RealMonomial(this) &
       & result(output) 
      class(RealMonomial), intent(in) :: this
      integer                         :: output
    end function
  end interface
  
  interface
    ! Evaluate the contribution to the energy from
    !    a univariate, monomial or polynomial at a given displacement.
    impure elemental module function energy_RealSingleDisplacement_RealUnivariate(this,displacement,paired_displacement) result(output) 
      class(RealUnivariate),        intent(in)           :: this
      type(RealSingleDisplacement), intent(in), optional :: displacement
      type(RealSingleDisplacement), intent(in), optional :: paired_displacement
      real(dp)                                           :: output
    end function
  end interface
  
  interface
    impure elemental module function energy_ComplexSingleDisplacement_RealUnivariate(   this,displacement,paired_displacement) result(output) 
      class(RealUnivariate),           intent(in)           :: this
      type(ComplexSingleDisplacement), intent(in), optional :: displacement
      type(ComplexSingleDisplacement), intent(in), optional :: paired_displacement
      complex(dp)                                           :: output
    end function
  end interface
  
  interface
    impure elemental module function energy_RealModeDisplacement_RealMonomial(this,displacement) result(output) 
      class(RealMonomial),        intent(in) :: this
      type(RealModeDisplacement), intent(in) :: displacement
      real(dp)                               :: output
    end function
  end interface
  
  interface
    impure elemental module function energy_ComplexModeDisplacement_RealMonomial(this,displacement) result(output) 
      class(RealMonomial),           intent(in) :: this
      type(ComplexModeDisplacement), intent(in) :: displacement
      complex(dp)                               :: output
    end function
  end interface
  
  interface
    impure elemental module function energy_RealModeDisplacement_RealPolynomial(this,displacement) result(output) 
      class(RealPolynomial),      intent(in) :: this
      type(RealModeDisplacement), intent(in) :: displacement
      real(dp)                               :: output
    end function
  end interface
  
  interface
    impure elemental module function energy_ComplexModeDisplacement_RealPolynomial(   this,displacement) result(output) 
      class(RealPolynomial),         intent(in) :: this
      type(ComplexModeDisplacement), intent(in) :: displacement
      complex(dp)                               :: output
    end function
  end interface
  
  interface
    ! Evaluate the contribution to the force from
    !    a univariate, monomial or polynomial at a given displacement.
    ! -d/d{u_i} ({u_i}^n) evaluated at u_i=U is -n*U^{n-1}
    module function force_RealSingleDisplacement_RealUnivariate(this, &
       & displacement,paired_displacement) result(output) 
      class(RealUnivariate),        intent(in)           :: this
      type(RealSingleDisplacement), intent(in), optional :: displacement
      type(RealSingleDisplacement), intent(in), optional :: paired_displacement
      real(dp), allocatable                              :: output(:)
    end function
  end interface
  
  interface
    module function force_ComplexSingleDisplacement_RealUnivariate(this, &
       & displacement,paired_displacement) result(output) 
      class(RealUnivariate),           intent(in)           :: this
      type(ComplexSingleDisplacement), intent(in), optional :: displacement
      type(ComplexSingleDisplacement), intent(in), optional :: paired_displacement
      complex(dp), allocatable                              :: output(:)
    end function
  end interface
  
  interface
    ! -d/d{u_i} (c*prod_j[ {u_j}^{n_j} ]) evaluated at {u_i=U_i} is
    !    -c*prod_{j/=i}[ {U_j}^{n_j} ] * n_i * {U_i}^{n_i-1}
    impure elemental module function force_RealModeDisplacement_RealMonomial(this,displacement) result(output) 
      class(RealMonomial),        intent(in) :: this
      type(RealModeDisplacement), intent(in) :: displacement
      type(RealModeForce)                    :: output
    end function
  end interface
  
  interface
    impure elemental module function force_ComplexModeDisplacement_RealMonomial(this,displacement) result(output) 
      class(RealMonomial),           intent(in) :: this
      type(ComplexModeDisplacement), intent(in) :: displacement
      type(ComplexModeForce)                    :: output
    end function
  end interface
  
  interface
    impure elemental module function force_RealModeDisplacement_RealPolynomial(this,displacement) result(output) 
      class(RealPolynomial),      intent(in) :: this
      type(RealModeDisplacement), intent(in) :: displacement
      type(RealModeForce)                    :: output
    end function
  end interface
  
  interface
    impure elemental module function force_ComplexModeDisplacement_RealPolynomial(this,displacement) result(output) 
      class(RealPolynomial),         intent(in) :: this
      type(ComplexModeDisplacement), intent(in) :: displacement
      type(ComplexModeForce)                    :: output
    end function
  end interface
  
  interface operator(*)
    ! Multiplication and division by scalars.
    impure elemental module function multiply_RealMonomial_real(this,that) &
       & result(output) 
      type(RealMonomial), intent(in) :: this
      real(dp),           intent(in) :: that
      type(RealMonomial)             :: output
    end function
  
    impure elemental module function multiply_real_RealMonomial(this,that) &
       & result(output) 
      real(dp),           intent(in) :: this
      type(RealMonomial), intent(in) :: that
      type(RealMonomial)             :: output
    end function
  
    impure elemental module function multiply_RealPolynomial_real(this,that) &
       & result(output) 
      type(RealPolynomial), intent(in) :: this
      real(dp),             intent(in) :: that
      type(RealPolynomial)             :: output
    end function
  
    impure elemental module function multiply_real_RealPolynomial(this,that) &
       & result(output) 
      real(dp),             intent(in) :: this
      type(RealPolynomial), intent(in) :: that
      type(RealPolynomial)             :: output
    end function
  end interface
  
  interface operator(/)
    impure elemental module function divide_RealMonomial_real(this,that) &
       & result(output) 
      type(RealMonomial), intent(in) :: this
      real(dp),           intent(in) :: that
      type(RealMonomial)             :: output
    end function
  
    impure elemental module function divide_RealPolynomial_real(this,that) &
       & result(output) 
      type(RealPolynomial), intent(in) :: this
      real(dp),             intent(in) :: that
      type(RealPolynomial)             :: output
    end function
  end interface
  
  interface operator(*)
    ! Multiplication between Monomial-like types.
    ! Uses a merge to keep ids in ascending order.
    impure elemental module function multiply_RealMonomialable_RealMonomialable(this,that) result(output) 
      class(RealMonomialable), intent(in) :: this
      class(RealMonomialable), intent(in) :: that
      type(RealMonomial)                  :: output
    end function
  
    ! Multiplication between polynomials and monomial-like types.
    impure elemental module function multiply_RealPolynomial_RealMonomialable(   this,that) result(output) 
      type(RealPolynomial),    intent(in) :: this
      class(RealMonomialable), intent(in) :: that
      type(RealPolynomial)                :: output
    end function
  
    impure elemental module function multiply_RealMonomialable_RealPolynomial(   this,that) result(output) 
      class(RealMonomialable), intent(in) :: this
      type(RealPolynomial),    intent(in) :: that
      type(RealPolynomial)                :: output
    end function
  end interface
  
  interface operator(+)
    ! Addition between polynomials and polynomial-like types.
    impure elemental module function add_RealPolynomialable_RealPolynomialable(this,that) result(output) 
      class(RealPolynomialable), intent(in) :: this
      class(RealPolynomialable), intent(in) :: that
      type(RealPolynomial)                  :: output
    end function
  end interface
  
  interface operator(-)
    ! The negative of a polynomial or polynomial-like type.
    impure elemental module function negative_RealPolynomialable(this) &
       & result(output) 
      class(RealPolynomialable), intent(in) :: this
      type(RealPolynomial)                  :: output
    end function
  
    ! Subtraction between polynomials and polynomial-like types.
    impure elemental module function subtract_RealPolynomialable_RealPolynomialable(   this,that) result(output) 
      class(RealPolynomialable), intent(in) :: this
      class(RealPolynomialable), intent(in) :: that
      type(RealPolynomial)                  :: output
    end function
  end interface
  
  interface sum
    ! Sum polynomial-like types.
    module function sum_RealPolynomialables(input) result(output) 
      class(RealPolynomialable), intent(in) :: input(:)
      type(RealPolynomial)                  :: output
    end function
  end interface
  
  interface select_mode
    ! ----------------------------------------------------------------------
    ! Select modes, displacements or forces corresponding to
    !    a given univariate or univariates.
    ! ----------------------------------------------------------------------
    module function select_mode_RealUnivariate(univariate,modes) &
       & result(output) 
      type(RealUnivariate), intent(in) :: univariate
      type(RealMode),       intent(in) :: modes(:)
      type(RealMode)                   :: output
    end function
  end interface
  
  interface select_modes
    module function select_modes_RealUnivariates(univariates,modes) &
       & result(output) 
      type(RealUnivariate), intent(in) :: univariates(:)
      type(RealMode),       intent(in) :: modes(:)
      type(RealMode), allocatable      :: output(:)
    end function
  end interface
  
  interface select_displacement
    module function select_displacement_RealUnivariate(univariate, &
       & displacements) result(output) 
      type(RealUnivariate),         intent(in) :: univariate
      type(RealSingleDisplacement), intent(in) :: displacements(:)
      type(RealSingleDisplacement)             :: output
    end function
  end interface
  
  interface select_displacements
    module function select_displacements_RealUnivariates(univariates, &
       & displacements) result(output) 
      type(RealUnivariate),         intent(in)  :: univariates(:)
      type(RealSingleDisplacement), intent(in)  :: displacements(:)
      type(RealSingleDisplacement), allocatable :: output(:)
    end function
  end interface
  
  interface select_force
    module function select_force_RealUnivariate(univariate,forces) &
       & result(output) 
      type(RealUnivariate),  intent(in) :: univariate
      type(RealSingleForce), intent(in) :: forces(:)
      type(RealSingleForce)             :: output
    end function
  end interface
  
  interface select_forces
    module function select_forces_RealUnivariates(univariates,forces) &
       & result(output) 
      type(RealUnivariate),  intent(in)  :: univariates(:)
      type(RealSingleForce), intent(in)  :: forces(:)
      type(RealSingleForce), allocatable :: output(:)
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Compares two monomials for equality up to coefficient.
    ! ----------------------------------------------------------------------
    module function compare_real_monomials(this,that) result(output) 
      class(*), intent(in) :: this
      class(*), intent(in) :: that
      logical              :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_RealUnivariate(this,input) 
      class(RealUnivariate), intent(out) :: this
      type(String),          intent(in)  :: input
    end subroutine
  end interface
  
  interface
    module function write_RealUnivariate(this) result(output) 
      class(RealUnivariate), intent(in) :: this
      type(String)                      :: output
    end function
  end interface
  
  interface RealUnivariate
    impure elemental module function new_RealUnivariate_String(input) &
       & result(this) 
      type(String), intent(in) :: input
      type(RealUnivariate)     :: this
    end function
  end interface
  
  interface
    module subroutine read_RealMonomial(this,input) 
      class(RealMonomial), intent(out) :: this
      type(String),        intent(in)  :: input
    end subroutine
  end interface
  
  interface
    module function write_RealMonomial(this) result(output) 
      class(RealMonomial), intent(in) :: this
      type(String)                    :: output
    end function
  end interface
  
  interface RealMonomial
    impure elemental module function new_RealMonomial_String(input) &
       & result(this) 
      type(String),   intent(in) :: input
      type(RealMonomial)         :: this
    end function
  end interface
  
  interface
    module subroutine read_RealPolynomial(this,input) 
      class(RealPolynomial), intent(out) :: this
      type(String),          intent(in)  :: input
    end subroutine
  end interface
  
  interface
    module function write_RealPolynomial(this) result(output) 
      class(RealPolynomial), intent(in) :: this
      type(String)                      :: output
    end function
  end interface
  
  interface RealPolynomial
    impure elemental module function new_RealPolynomial_String(input) &
       & result(this) 
      type(String),   intent(in) :: input
      type(RealPolynomial)       :: this
    end function
  end interface
end module
