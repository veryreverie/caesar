! ======================================================================
! Polynomials in the complex representation, but with monomials stored as
!    (m+m*) and i(m-m*) pairs with real coefficients.
! ======================================================================
module caesar_paired_polynomial_module
  use caesar_utils_module
  use caesar_structure_module
  
  use caesar_real_mode_module
  use caesar_real_single_mode_displacement_module
  use caesar_real_single_mode_force_module
  use caesar_real_mode_displacement_module
  use caesar_real_mode_force_module
  use caesar_complex_polynomial_module
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
  
  abstract interface
    function to_PairedPolynomial_PairedPolynomialable(this) result(output)
      import PairedPolynomial
      import PairedPolynomialable
      implicit none
      
      class(PairedPolynomialable), intent(in) :: this
      type(PairedPolynomial)                  :: output
    end function
  end interface
  
  interface PairedMonomial
    ! ----------------------------------------------------------------------
    ! Constructors.
    ! ----------------------------------------------------------------------
    module function new_PairedMonomial(coefficient,modes,pair) result(this) 
      real(dp),                intent(in) :: coefficient
      type(ComplexUnivariate), intent(in) :: modes(:)
      integer,                 intent(in) :: pair
      type(PairedMonomial)                :: this
    end function
  end interface
  
  interface PairedPolynomial
    module function new_PairedPolynomial(terms) result(this) 
      type(PairedMonomial), intent(in) :: terms(:)
      type(PairedPolynomial)           :: this
    end function
  
    module function new_PairedPolynomial_PairedPolynomialable(input) &
       & result(this) 
      class(PairedPolynomialable), intent(in) :: input
      type(PairedPolynomial)                  :: this
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Conversions between types.
    ! ----------------------------------------------------------------------
    module function to_PairedPolynomial_PairedMonomial(this) result(output) 
      class(PairedMonomial), intent(in) :: this
      type(PairedPolynomial)            :: output
    end function
  end interface
  
  interface
    module function to_PairedPolynomial_PairedPolynomial(this) result(output) 
      class(PairedPolynomial), intent(in) :: this
      type(PairedPolynomial)              :: output
    end function
  end interface
  
  interface size
    ! ----------------------------------------------------------------------
    ! Operations involving types.
    ! ----------------------------------------------------------------------
    
    ! The number of modes in a monomial, or terms in a polynomial.
    module function size_PairedMonomial(this) result(output) 
      type(PairedMonomial), intent(in) :: this
      integer                          :: output
    end function
  
    module function size_PairedPolynomial(this) result(output) 
      class(PairedPolynomial), intent(in) :: this
      integer                             :: output
    end function
  end interface
  
  interface
    ! Getters for monomials.
    impure elemental module function mode_PairedMonomial(this,index) &
       & result(output) 
      class(PairedMonomial), intent(in) :: this
      integer,               intent(in) :: index
      type(ComplexUnivariate)           :: output
    end function
  end interface
  
  interface
    ! Getter for modes. Has several run types:
    !    - If no arguments specified, returns all modes.
    !    - If indices specified, returns modes at specified indices.
    !    - If ids and paired_ids specified, returns modes with specified ids.
    !    - If exclude_ids specicified, returns all modes but those excluded.
    module function modes_PairedMonomial(this,indices,ids,paired_ids, &
       & exclude_ids) result(output) 
      class(PairedMonomial), intent(in)           :: this
      integer,               intent(in), optional :: indices(:)
      integer,               intent(in), optional :: ids(:)
      integer,               intent(in), optional :: paired_ids(:)
      integer,               intent(in), optional :: exclude_ids(:)
      type(ComplexUnivariate), allocatable        :: output(:)
    end function
  end interface
  
  interface
    impure elemental module function id_PairedMonomial(this,index) &
       & result(output) 
      class(PairedMonomial), intent(in) :: this
      integer,               intent(in) :: index
      integer                           :: output
    end function
  end interface
  
  interface
    module function ids_PairedMonomial(this,indices) result(output) 
      class(PairedMonomial), intent(in)           :: this
      integer,               intent(in), optional :: indices(:)
      integer, allocatable                        :: output(:)
    end function
  end interface
  
  interface
    impure elemental module function paired_id_PairedMonomial(this,index) &
       & result(output) 
      class(PairedMonomial), intent(in) :: this
      integer,               intent(in) :: index
      integer                           :: output
    end function
  end interface
  
  interface
    module function paired_ids_PairedMonomial(this,indices) result(output) 
      class(PairedMonomial), intent(in)           :: this
      integer,               intent(in), optional :: indices(:)
      integer, allocatable                        :: output(:)
    end function
  end interface
  
  interface
    impure elemental module function power_PairedMonomial(this,index) &
       & result(output) 
      class(PairedMonomial), intent(in) :: this
      integer,               intent(in) :: index
      integer                           :: output
    end function
  end interface
  
  interface
    module function powers_PairedMonomial(this,indices) result(output) 
      class(PairedMonomial), intent(in)           :: this
      integer,               intent(in), optional :: indices(:)
      integer, allocatable                        :: output(:)
    end function
  end interface
  
  interface
    impure elemental module function paired_power_PairedMonomial(this,index) &
       & result(output) 
      class(PairedMonomial), intent(in) :: this
      integer,               intent(in) :: index
      integer                           :: output
    end function
  end interface
  
  interface
    module function paired_powers_PairedMonomial(this,indices) result(output) 
      class(PairedMonomial), intent(in)           :: this
      integer,               intent(in), optional :: indices(:)
      integer, allocatable                        :: output(:)
    end function
  end interface
  
  interface
    ! Set the modes of a monomial.
    module subroutine set_modes_PairedMonomial(this,modes,pair) 
      class(PairedMonomial),   intent(inout)        :: this
      type(ComplexUnivariate), intent(in)           :: modes(:)
      integer,                 intent(in), optional :: pair
    end subroutine
  end interface
  
  interface
    ! Simplify a monomial or polynomial.
    impure elemental module subroutine simplify_PairedMonomial(this) 
      class(PairedMonomial), intent(inout) :: this
    end subroutine
  end interface
  
  interface
    impure elemental module subroutine simplify_PairedPolynomial(this) 
      class(PairedPolynomial), intent(inout) :: this
    end subroutine
  end interface
  
  interface conjg
    ! Find the conjugate of a monomial or polynomial.
    impure elemental module function conjg_PairedMonomial(this) result(output) 
      type(PairedMonomial), intent(in) :: this
      type(PairedMonomial)             :: output
    end function
  
    impure elemental module function conjg_PairedPolynomial(this) &
       & result(output) 
      type(PairedPolynomial), intent(in) :: this
      type(PairedPolynomial)             :: output
    end function
  end interface
  
  interface
    ! The total power of a monomial.
    impure elemental module function total_power_PairedMonomial(this) &
       & result(output) 
      class(PairedMonomial), intent(in) :: this
      integer                           :: output
    end function
  end interface
  
  interface
    ! Evaluate the contribution to the energy from a monomial or polynomial
    !    at a given displacement.
    impure elemental module function energy_RealModeDisplacement_PairedMonomial(this,displacement) result(output) 
      class(PairedMonomial),      intent(in) :: this
      type(RealModeDisplacement), intent(in) :: displacement
      real(dp)                               :: output
    end function
  end interface
  
  interface
    impure elemental module function energy_RealModeDisplacement_PairedPolynomial(   this,displacement) result(output) 
      class(PairedPolynomial),     intent(in) :: this
      class(RealModeDisplacement), intent(in) :: displacement
      real(dp)                                :: output
    end function
  end interface
  
  interface
    ! Evaluate the contribution to the force from
    !    a monomial or polynomial at a given displacement.
    
    ! -d/d{u_i} (c*prod_j[ {u_j}^{n_j} ]) evaluated at {u_i=U_i} is
    !    -c*prod_{j/=i}[ {U_j}^{n_j} ] * n_i * {U_i}^{n_i-1}
    impure elemental module function force_RealModeDisplacement_PairedMonomial(this,displacement) result(output) 
      class(PairedMonomial),      intent(in) :: this
      type(RealModeDisplacement), intent(in) :: displacement
      type(RealModeForce)                    :: output
    end function
  end interface
  
  interface
    impure elemental module function force_RealModeDisplacement_PairedPolynomial(this,displacement) result(output) 
      class(PairedPolynomial),    intent(in) :: this
      type(RealModeDisplacement), intent(in) :: displacement
      type(RealModeForce)                    :: output
    end function
  end interface
  
  interface
    ! The harmonic expectation of a monomial or polynomial.
    impure elemental module function harmonic_expectation_PairedMonomial(this,frequency,thermal_energy,supercell_size) result(output) 
      class(PairedMonomial), intent(in) :: this
      real(dp),              intent(in) :: frequency
      real(dp),              intent(in) :: thermal_energy
      integer,               intent(in) :: supercell_size
      real(dp)                          :: output
    end function
  end interface
  
  interface
    impure elemental module function harmonic_expectation_PairedPolynomial(this,frequency,thermal_energy,supercell_size) result(output) 
      class(PairedPolynomial), intent(in) :: this
      real(dp),                intent(in) :: frequency
      real(dp),                intent(in) :: thermal_energy
      integer,                 intent(in) :: supercell_size
      real(dp)                            :: output
    end function
  end interface
  
  interface operator(*)
    ! Multiplication and division by scalars.
    impure elemental module function multiply_PairedMonomial_real(this,that) &
       & result(output) 
      type(PairedMonomial), intent(in) :: this
      real(dp),             intent(in) :: that
      type(PairedMonomial)             :: output
    end function
  
    impure elemental module function multiply_real_PairedMonomial(this,that) &
       & result(output) 
      real(dp),             intent(in) :: this
      type(PairedMonomial), intent(in) :: that
      type(PairedMonomial)             :: output
    end function
  
    impure elemental module function multiply_PairedPolynomial_real(this, &
       & that) result(output) 
      type(PairedPolynomial), intent(in) :: this
      real(dp),               intent(in) :: that
      type(PairedPolynomial)             :: output
    end function
  
    impure elemental module function multiply_real_PairedPolynomial(this, &
       & that) result(output) 
      real(dp),               intent(in) :: this
      type(PairedPolynomial), intent(in) :: that
      type(PairedPolynomial)             :: output
    end function
  end interface
  
  interface operator(/)
    impure elemental module function divide_PairedMonomial_real(this,that) &
       & result(output) 
      type(PairedMonomial), intent(in) :: this
      real(dp),             intent(in) :: that
      type(PairedMonomial)             :: output
    end function
  
    impure elemental module function divide_PairedPolynomial_real(this,that) &
       & result(output) 
      type(PairedPolynomial), intent(in) :: this
      real(dp),               intent(in) :: that
      type(PairedPolynomial)             :: output
    end function
  end interface
  
  interface operator(+)
    ! Addition between polynomials and polynomial-like types.
    impure elemental module function add_PairedPolynomialable_PairedPolynomialable(   this,that) result(output) 
      class(PairedPolynomialable), intent(in) :: this
      class(PairedPolynomialable), intent(in) :: that
      type(PairedPolynomial)                  :: output
    end function
  end interface
  
  interface operator(-)
    ! The negative of a monomial or polynomial.
    impure elemental module function negative_PairedMonomial(this) &
       & result(output) 
      class(PairedMonomial), intent(in) :: this
      type(PairedMonomial)              :: output
    end function
  
    impure elemental module function negative_PairedPolynomial(this) &
       & result(output) 
      class(PairedPolynomial), intent(in) :: this
      type(PairedPolynomial)              :: output
    end function
  
    ! Subtraction between polynomials and polynomial-like types.
    impure elemental module function                                              subtract_PairedPolynomialable_PairedPolynomialable(this,that) result(output) 
      class(PairedPolynomialable), intent(in) :: this
      class(PairedPolynomialable), intent(in) :: that
      type(PairedPolynomial)                  :: output
    end function
  end interface
  
  interface sum
    ! Sum polynomial-like types.
    module function sum_PairedPolynomialables(input) result(output) 
      class(PairedPolynomialable), intent(in) :: input(:)
      type(PairedPolynomial)                  :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Compares two monomials for equality up to coefficient.
    ! ----------------------------------------------------------------------
    module function compare_paired_monomials(this,that) result(output) 
      class(*), intent(in) :: this
      class(*), intent(in) :: that
      logical              :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Returns .true. if the arguments have matching modes, and:
    !    this%pair_ ==  0 and that%pair_ == 0
    ! or this%pair_ ==  1 and that%pair_== -1
    ! or this%pair_ == -1 and that%pair_==  1
    ! ----------------------------------------------------------------------
    module function matching_pair(this,that) result(output) 
      class(*), intent(in) :: this
      class(*), intent(in) :: that
      logical              :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_PairedMonomial(this,input) 
      class(PairedMonomial), intent(out) :: this
      type(String),           intent(in)  :: input
    end subroutine
  end interface
  
  interface
    module function write_PairedMonomial(this) result(output) 
      class(PairedMonomial), intent(in) :: this
      type(String)                       :: output
    end function
  end interface
  
  interface PairedMonomial
    impure elemental module function new_PairedMonomial_String(input) &
       & result(this) 
      type(String), intent(in) :: input
      type(PairedMonomial)     :: this
    end function
  end interface
  
  interface
    module subroutine read_PairedPolynomial(this,input) 
      class(PairedPolynomial), intent(out) :: this
      type(String),            intent(in)  :: input
    end subroutine
  end interface
  
  interface
    module function write_PairedPolynomial(this) result(output) 
      class(PairedPolynomial), intent(in) :: this
      type(String)                         :: output
    end function
  end interface
  
  interface PairedPolynomial
    impure elemental module function new_PairedPolynomial_String(input) &
       & result(this) 
      type(String), intent(in) :: input
      type(PairedPolynomial)   :: this
    end function
  end interface
end module
