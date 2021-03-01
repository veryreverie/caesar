! ======================================================================
! The building blocks of basis functions in complex co-ordinates.
! ======================================================================
module caesar_complex_polynomial_module
  use caesar_utils_module
  use caesar_structure_module
  
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
  end type
  
  type, abstract, extends(ComplexPolynomialable) :: ComplexMonomialable
  contains
    procedure(to_ComplexMonomial_ComplexMonomialable), deferred, public :: &
       & to_ComplexMonomial
    
    procedure(wavevector_ComplexMonomialable), deferred, public &
       & :: wavevector
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
    
    procedure, public :: energy_real => &
                       & energy_RealSingleDisplacement_ComplexUnivariate
    procedure, public :: energy_complex => &
                       & energy_ComplexSingleDisplacement_ComplexUnivariate
    
    procedure, public :: force_real => &
                       & force_RealSingleDisplacement_ComplexUnivariate
    procedure, public :: force_complex => &
                       & force_ComplexSingleDisplacement_ComplexUnivariate
    
    procedure, public :: harmonic_expectation => &
                       & harmonic_expectation_ComplexUnivariate
    
    procedure, public :: read  => read_ComplexUnivariate
    procedure, public :: write => write_ComplexUnivariate
  end type
  
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
  
  type, extends(ComplexPolynomialable) :: ComplexPolynomial
    type(ComplexMonomial), allocatable :: terms(:)
  contains
    procedure, public :: to_ComplexPolynomial => &
       & to_ComplexPolynomial_ComplexPolynomial
    
    procedure, public :: simplify => simplify_ComplexPolynomial
    
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
    
    function wavevector_ComplexMonomialable(this,modes,qpoints) &
       & result(output)
      import ComplexMonomialable
      import ComplexMode
      import QpointData
      import FractionVector
      implicit none
      
      class(ComplexMonomialable), intent(in) :: this
      type(ComplexMode),            intent(in) :: modes(:)
      type(QpointData),             intent(in) :: qpoints(:)
      type(FractionVector)                     :: output
    end function
  end interface
  
  interface ComplexUnivariate
    ! ----------------------------------------------------------------------
    ! Constructors.
    ! ----------------------------------------------------------------------
    impure elemental module function new_ComplexUnivariate(id,paired_id, &
       & power,paired_power) result(this) 
      integer, intent(in)     :: id
      integer, intent(in)     :: paired_id
      integer, intent(in)     :: power
      integer, intent(in)     :: paired_power
      type(ComplexUnivariate) :: this
    end function
  
    impure elemental module function new_ComplexUnivariate_ComplexMode(mode, &
       & power,paired_power) result(this) 
      type(ComplexMode), intent(in)           :: mode
      integer,           intent(in)           :: power
      integer,           intent(in), optional :: paired_power
      type(ComplexUnivariate)                 :: this
    end function
  end interface
  
  interface ComplexMonomial
    module function new_ComplexMonomial(coefficient,modes) result(this) 
      complex(dp),             intent(in) :: coefficient
      type(ComplexUnivariate), intent(in) :: modes(:)
      type(ComplexMonomial)               :: this
    end function
  
    module function new_ComplexMonomial_ComplexMonomialable(input) &
       & result(this) 
      class(ComplexMonomialable), intent(in) :: input
      type(ComplexMonomial)                  :: this
    end function
  end interface
  
  interface ComplexPolynomial
    module function new_ComplexPolynomial(terms) result(this) 
      type(ComplexMonomial), intent(in) :: terms(:)
      type(ComplexPolynomial)           :: this
    end function
  
    module function new_ComplexPolynomial_ComplexPolynomialable(input) &
       & result(this) 
      class(ComplexPolynomialable), intent(in) :: input
      type(ComplexPolynomial)                  :: this
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Conversions between types.
    ! ----------------------------------------------------------------------
    module function to_ComplexMonomial_ComplexUnivariate(this) result(output) 
      class(ComplexUnivariate), intent(in) :: this
      type(ComplexMonomial)                :: output
    end function
  end interface
  
  interface
    module function to_ComplexPolynomial_ComplexUnivariate(this) &
       & result(output) 
      class(ComplexUnivariate), intent(in) :: this
      type(ComplexPolynomial)              :: output
    end function
  end interface
  
  interface
    module function to_ComplexMonomial_ComplexMonomial(this) result(output) 
      class(ComplexMonomial), intent(in) :: this
      type(ComplexMonomial)              :: output
    end function
  end interface
  
  interface
    module function to_ComplexPolynomial_ComplexMonomial(this) result(output) 
      class(ComplexMonomial), intent(in) :: this
      type(ComplexPolynomial)            :: output
    end function
  end interface
  
  interface
    module function to_ComplexPolynomial_ComplexPolynomial(this) &
       & result(output) 
      class(ComplexPolynomial), intent(in) :: this
      type(ComplexPolynomial)              :: output
    end function
  end interface
  
  interface size
    ! ----------------------------------------------------------------------
    ! Operations involving types.
    ! ----------------------------------------------------------------------
    
    ! The number of modes in a monomial, or terms in a polynomial.
    module function size_ComplexMonomial(this) result(output) 
      type(ComplexMonomial), intent(in) :: this
      integer                           :: output
    end function
  
    module function size_ComplexPolynomial(this) result(output) 
      class(ComplexPolynomial), intent(in) :: this
      integer                              :: output
    end function
  end interface
  
  interface
    ! Getters for monomials.
    impure elemental module function mode_ComplexMonomial(this,index) &
       & result(output) 
      class(ComplexMonomial), intent(in) :: this
      integer,                intent(in) :: index
      type(ComplexUnivariate)            :: output
    end function
  end interface
  
  interface
    ! Getter for modes. Has several run types:
    !    - If no arguments specified, returns all modes.
    !    - If indices specified, returns modes at specified indices.
    !    - If ids and paired_ids specified, returns modes with specified ids.
    !    - If exclude_ids specicified, returns all modes but those excluded.
    module function modes_ComplexMonomial(this,indices,ids,paired_ids, &
       & exclude_ids) result(output) 
      class(ComplexMonomial), intent(in)           :: this
      integer,                intent(in), optional :: indices(:)
      integer,                intent(in), optional :: ids(:)
      integer,                intent(in), optional :: paired_ids(:)
      integer,                intent(in), optional :: exclude_ids(:)
      type(ComplexUnivariate), allocatable         :: output(:)
    end function
  end interface
  
  interface
    impure elemental module function id_ComplexMonomial(this,index) &
       & result(output) 
      class(ComplexMonomial), intent(in) :: this
      integer,                intent(in) :: index
      integer                            :: output
    end function
  end interface
  
  interface
    module function ids_ComplexMonomial(this,indices) result(output) 
      class(ComplexMonomial), intent(in)           :: this
      integer,                intent(in), optional :: indices(:)
      integer, allocatable                         :: output(:)
    end function
  end interface
  
  interface
    impure elemental module function paired_id_ComplexMonomial(this,index) &
       & result(output) 
      class(ComplexMonomial), intent(in) :: this
      integer,                intent(in) :: index
      integer                            :: output
    end function
  end interface
  
  interface
    module function paired_ids_ComplexMonomial(this,indices) result(output) 
      class(ComplexMonomial), intent(in)           :: this
      integer,                intent(in), optional :: indices(:)
      integer, allocatable                         :: output(:)
    end function
  end interface
  
  interface
    impure elemental module function power_ComplexMonomial(this,index) &
       & result(output) 
      class(ComplexMonomial), intent(in) :: this
      integer,                intent(in) :: index
      integer                            :: output
    end function
  end interface
  
  interface
    module function powers_ComplexMonomial(this,indices) result(output) 
      class(ComplexMonomial), intent(in)           :: this
      integer,                intent(in), optional :: indices(:)
      integer, allocatable                         :: output(:)
    end function
  end interface
  
  interface
    impure elemental module function paired_power_ComplexMonomial(this,index) &
       & result(output) 
      class(ComplexMonomial), intent(in) :: this
      integer,                intent(in) :: index
      integer                            :: output
    end function
  end interface
  
  interface
    module function paired_powers_ComplexMonomial(this,indices) result(output) 
      class(ComplexMonomial), intent(in)           :: this
      integer,                intent(in), optional :: indices(:)
      integer, allocatable                         :: output(:)
    end function
  end interface
  
  interface
    ! Set the modes of a monomial.
    module subroutine set_modes_ComplexMonomial(this,modes) 
      class(ComplexMonomial),  intent(inout) :: this
      type(ComplexUnivariate), intent(in)    :: modes(:)
    end subroutine
  end interface
  
  interface
    ! Simplify a monomial or polynomial.
    impure elemental module subroutine simplify_ComplexMonomial(this) 
      class(ComplexMonomial), intent(inout) :: this
    end subroutine
  end interface
  
  interface
    impure elemental module subroutine simplify_ComplexPolynomial(this) 
      class(ComplexPolynomial), intent(inout) :: this
    end subroutine
  end interface
  
  interface conjg
    ! Find the conjugate of a univariate, monomial or polynomial.
    impure elemental module function conjg_ComplexUnivariate(this) &
       & result(output) 
      type(ComplexUnivariate), intent(in) :: this
      type(ComplexUnivariate)             :: output
    end function
  
    impure elemental module function conjg_ComplexMonomial(this) &
       & result(output) 
      type(ComplexMonomial), intent(in) :: this
      type(ComplexMonomial)             :: output
    end function
  
    impure elemental module function conjg_ComplexPolynomial(this) &
       & result(output) 
      type(ComplexPolynomial), intent(in) :: this
      type(ComplexPolynomial)             :: output
    end function
  end interface
  
  interface
    ! The total power of a univariate or monomial.
    impure elemental module function total_power_ComplexUnivariate(this) &
       & result(output) 
      class(ComplexUnivariate), intent(in) :: this
      integer                              :: output
    end function
  end interface
  
  interface
    impure elemental module function total_power_ComplexMonomial(this) &
       & result(output) 
      class(ComplexMonomial), intent(in) :: this
      integer                            :: output
    end function
  end interface
  
  interface
    ! Returns the Bloch wavevector of a univariate or monomial.
    module function wavevector_ComplexUnivariate(this,modes,qpoints) &
       & result(output) 
      class(ComplexUnivariate), intent(in) :: this
      type(ComplexMode),        intent(in) :: modes(:)
      type(QpointData),         intent(in) :: qpoints(:)
      type(FractionVector)                 :: output
    end function
  end interface
  
  interface
    module function wavevector_ComplexMonomial(this,modes,qpoints) &
       & result(output) 
      class(ComplexMonomial), intent(in) :: this
      type(ComplexMode),      intent(in) :: modes(:)
      type(QpointData),       intent(in) :: qpoints(:)
      type(FractionVector)               :: output
    end function
  end interface
  
  interface
    ! Evaluate the contribution to the energy from
    !    a univariate, monomial or polynomial at a given displacement.
    impure elemental module function energy_RealSingleDisplacement_ComplexUnivariate(   this,displacement,paired_displacement) result(output) 
      class(ComplexUnivariate),     intent(in)           :: this
      type(RealSingleDisplacement), intent(in), optional :: displacement
      type(RealSingleDisplacement), intent(in), optional :: paired_displacement
      complex(dp)                                        :: output
    end function
  end interface
  
  interface
    impure elemental module function energy_ComplexSingleDisplacement_ComplexUnivariate(   this,displacement,paired_displacement) result(output) 
      class(ComplexUnivariate),         intent(in)           :: this
      class(ComplexSingleDisplacement), intent(in), optional :: displacement
      class(ComplexSingleDisplacement), intent(in), optional :: paired_displacement
      complex(dp)                                            :: output
    end function
  end interface
  
  interface
    impure elemental module function energy_RealModeDisplacement_ComplexMonomial(this,displacement) result(output) 
      class(ComplexMonomial),     intent(in) :: this
      type(RealModeDisplacement), intent(in) :: displacement
      complex(dp)                            :: output
    end function
  end interface
  
  interface
    impure elemental module function energy_ComplexModeDisplacement_ComplexMonomial(   this,displacement) result(output) 
      class(ComplexMonomial),        intent(in) :: this
      type(ComplexModeDisplacement), intent(in) :: displacement
      complex(dp)                               :: output
    end function
  end interface
  
  interface
    impure elemental module function energy_RealModeDisplacement_ComplexPolynomial(   this,displacement) result(output) 
      class(ComplexPolynomial),    intent(in) :: this
      class(RealModeDisplacement), intent(in) :: displacement
      complex(dp)                             :: output
    end function
  end interface
  
  interface
    impure elemental module function energy_ComplexModeDisplacement_ComplexPolynomial(   this,displacement) result(output) 
      class(ComplexPolynomial),       intent(in) :: this
      class(ComplexModeDisplacement), intent(in) :: displacement
      complex(dp)                                :: output
    end function
  end interface
  
  interface
    ! Evaluate the contribution to the force from
    !    a univariate, monomial or polynomial at a given displacement.
    ! -d/d{u_i} ({u_i}^n) evaluated at u_i=U is -n*U^{n-1}
    module function force_RealSingleDisplacement_ComplexUnivariate(this, &
       & displacement,paired_displacement) result(output) 
      class(ComplexUnivariate),     intent(in)           :: this
      type(RealSingleDisplacement), intent(in), optional :: displacement
      type(RealSingleDisplacement), intent(in), optional :: paired_displacement
      complex(dp), allocatable                           :: output(:)
    end function
  end interface
  
  interface
    module function force_ComplexSingleDisplacement_ComplexUnivariate(this, &
       & displacement,paired_displacement) result(output) 
      class(ComplexUnivariate),        intent(in)           :: this
      type(ComplexSingleDisplacement), intent(in), optional :: displacement
      type(ComplexSingleDisplacement), intent(in), optional :: paired_displacement
      complex(dp), allocatable                              :: output(:)
    end function
  end interface
  
  interface
    ! -d/d{u_i} (c*prod_j[ {u_j}^{n_j} ]) evaluated at {u_i=U_i} is
    !    -c*prod_{j/=i}[ {U_j}^{n_j} ] * n_i * {U_i}^{n_i-1}
    impure elemental module function force_RealModeDisplacement_ComplexMonomial(this,displacement) result(output) 
      class(ComplexMonomial),     intent(in) :: this
      type(RealModeDisplacement), intent(in) :: displacement
      type(RealModeForce)                    :: output
    end function
  end interface
  
  interface
    impure elemental module function force_ComplexModeDisplacement_ComplexMonomial(   this,displacement) result(output) 
      class(ComplexMonomial),        intent(in) :: this
      type(ComplexModeDisplacement), intent(in) :: displacement
      type(ComplexModeForce)                    :: output
    end function
  end interface
  
  interface
    impure elemental module function force_RealModeDisplacement_ComplexPolynomial(this,displacement) result(output) 
      class(ComplexPolynomial),   intent(in) :: this
      type(RealModeDisplacement), intent(in) :: displacement
      type(RealModeForce)                    :: output
    end function
  end interface
  
  interface
    impure elemental module function force_ComplexModeDisplacement_ComplexPolynomial(   this,displacement) result(output) 
      class(ComplexPolynomial),      intent(in) :: this
      type(ComplexModeDisplacement), intent(in) :: displacement
      type(ComplexModeForce)                    :: output
    end function
  end interface
  
  interface
    ! The harmonic expectation of a monomial or polynomial.
    impure elemental module function harmonic_expectation_ComplexPolynomial(this,frequency,thermal_energy,supercell_size) result(output) 
      class(ComplexPolynomial), intent(in) :: this
      real(dp),                 intent(in) :: frequency
      real(dp),                 intent(in) :: thermal_energy
      integer,                  intent(in) :: supercell_size
      real(dp)                             :: output
    end function
  end interface
  
  interface
    impure elemental module function harmonic_expectation_ComplexMonomial(this,frequency,thermal_energy,supercell_size) result(output) 
      class(ComplexMonomial), intent(in) :: this
      real(dp),               intent(in) :: frequency
      real(dp),               intent(in) :: thermal_energy
      integer,                intent(in) :: supercell_size
      real(dp)                           :: output
    end function
  end interface
  
  interface
    impure elemental module function harmonic_expectation_ComplexUnivariate(this,frequency,thermal_energy,supercell_size) result(output) 
      class(ComplexUnivariate), intent(in) :: this
      real(dp),                 intent(in) :: frequency
      real(dp),                 intent(in) :: thermal_energy
      integer,                  intent(in) :: supercell_size
      real(dp)                             :: output
    end function
  end interface
  
  interface operator(*)
    ! Multiplication and division by scalars.
    impure elemental module function multiply_ComplexMonomial_real(this,that) &
       & result(output) 
      type(ComplexMonomial), intent(in) :: this
      real(dp),              intent(in) :: that
      type(ComplexMonomial)             :: output
    end function
  
    impure elemental module function multiply_real_ComplexMonomial(this,that) &
       & result(output) 
      real(dp),              intent(in) :: this
      type(ComplexMonomial), intent(in) :: that
      type(ComplexMonomial)             :: output
    end function
  
    impure elemental module function multiply_ComplexMonomial_complex(this, &
       & that) result(output) 
      type(ComplexMonomial), intent(in) :: this
      complex(dp),           intent(in) :: that
      type(ComplexMonomial)             :: output
    end function
  
    impure elemental module function multiply_complex_ComplexMonomial(this, &
       & that) result(output) 
      complex(dp),           intent(in) :: this
      type(ComplexMonomial), intent(in) :: that
      type(ComplexMonomial)             :: output
    end function
  
    impure elemental module function multiply_ComplexPolynomial_real(this, &
       & that) result(output) 
      type(ComplexPolynomial), intent(in) :: this
      real(dp),                intent(in) :: that
      type(ComplexPolynomial)             :: output
    end function
  
    impure elemental module function multiply_real_ComplexPolynomial(this, &
       & that) result(output) 
      real(dp),                intent(in) :: this
      type(ComplexPolynomial), intent(in) :: that
      type(ComplexPolynomial)             :: output
    end function
  
    impure elemental module function multiply_ComplexPolynomial_complex(this,that) result(output) 
      type(ComplexPolynomial), intent(in) :: this
      complex(dp),             intent(in) :: that
      type(ComplexPolynomial)             :: output
    end function
  
    impure elemental module function multiply_complex_ComplexPolynomial(this,that) result(output) 
      complex(dp),             intent(in) :: this
      type(ComplexPolynomial), intent(in) :: that
      type(ComplexPolynomial)             :: output
    end function
  end interface
  
  interface operator(/)
    impure elemental module function divide_ComplexMonomial_real(this,that) &
       & result(output) 
      type(ComplexMonomial), intent(in) :: this
      real(dp),              intent(in) :: that
      type(ComplexMonomial)             :: output
    end function
  
    impure elemental module function divide_ComplexMonomial_complex(this, &
       & that) result(output) 
      type(ComplexMonomial), intent(in) :: this
      complex(dp),           intent(in) :: that
      type(ComplexMonomial)             :: output
    end function
  
    impure elemental module function divide_ComplexPolynomial_real(this,that) &
       & result(output) 
      type(ComplexPolynomial), intent(in) :: this
      real(dp),                intent(in) :: that
      type(ComplexPolynomial)             :: output
    end function
  
    impure elemental module function divide_ComplexPolynomial_complex(this, &
       & that) result(output) 
      type(ComplexPolynomial), intent(in) :: this
      complex(dp),             intent(in) :: that
      type(ComplexPolynomial)             :: output
    end function
  end interface
  
  interface operator(*)
    ! Multiplication between Monomial-like types.
    ! Uses a merge to keep ids in ascending order.
    impure elemental module function multiply_ComplexMonomialable_ComplexMonomialable(   this,that) result(output) 
      class(ComplexMonomialable), intent(in) :: this
      class(ComplexMonomialable), intent(in) :: that
      type(ComplexMonomial)                  :: output
    end function
  
    ! Multiplication between polynomials and monomial-like types.
    impure elemental module function multiply_ComplexPolynomial_ComplexMonomialable(   this,that) result(output) 
      type(ComplexPolynomial),    intent(in) :: this
      class(ComplexMonomialable), intent(in) :: that
      type(ComplexPolynomial)                :: output
    end function
  
    impure elemental module function multiply_ComplexMonomialable_ComplexPolynomial(   this,that) result(output) 
      class(ComplexMonomialable), intent(in) :: this
      type(ComplexPolynomial),    intent(in) :: that
      type(ComplexPolynomial)                :: output
    end function
  end interface
  
  interface operator(+)
    ! Addition between polynomials and polynomial-like types.
    impure elemental module function add_ComplexPolynomialable_ComplexPolynomialable(   this,that) result(output) 
      class(ComplexPolynomialable), intent(in) :: this
      class(ComplexPolynomialable), intent(in) :: that
      type(ComplexPolynomial)                  :: output
    end function
  end interface
  
  interface operator(-)
    ! The negative of a monomial or polynomial.
    impure elemental module function negative_ComplexMonomial(this) &
       & result(output) 
      class(ComplexMonomial), intent(in) :: this
      type(ComplexMonomial)              :: output
    end function
  
    impure elemental module function negative_ComplexPolynomial(this) &
       & result(output) 
      class(ComplexPolynomial), intent(in) :: this
      type(ComplexPolynomial)              :: output
    end function
  
    ! Subtraction between polynomials and polynomial-like types.
    impure elemental module function                                              subtract_ComplexPolynomialable_ComplexPolynomialable(this,that) result(output) 
      class(ComplexPolynomialable), intent(in) :: this
      class(ComplexPolynomialable), intent(in) :: that
      type(ComplexPolynomial)                  :: output
    end function
  end interface
  
  interface sum
    ! Sum polynomial-like types.
    module function sum_ComplexPolynomialables(input) result(output) 
      class(ComplexPolynomialable), intent(in) :: input(:)
      type(ComplexPolynomial)                  :: output
    end function
  end interface
  
  interface select_modes
    ! ----------------------------------------------------------------------
    ! Select modes, displacements or forces corresponding to
    !    a given univariate or univariates.
    ! ----------------------------------------------------------------------
    module function select_modes_ComplexUnivariate(input,modes) result(output) 
      type(ComplexUnivariate), intent(in) :: input
      type(ComplexMode),       intent(in) :: modes(:)
      type(ComplexMode), allocatable      :: output(:)
    end function
  
    module function select_modes_ComplexUnivariates(input,modes) &
       & result(output) 
      type(ComplexUnivariate), intent(in) :: input(:)
      type(ComplexMode),       intent(in) :: modes(:)
      type(ComplexMode), allocatable      :: output(:)
    end function
  end interface
  
  interface select_displacements
    module function select_displacements_ComplexUnivariate(input, &
       & displacements) result(output) 
      type(ComplexUnivariate),          intent(in) :: input
      type(ComplexSingleDisplacement),  intent(in) :: displacements(:)
      type(ComplexSingleDisplacement), allocatable :: output(:)
    end function
  
    module function select_displacements_ComplexUnivariates(input, &
       & displacements) result(output) 
      type(ComplexUnivariate),         intent(in)  :: input(:)
      type(ComplexSingleDisplacement), intent(in)  :: displacements(:)
      type(ComplexSingleDisplacement), allocatable :: output(:)
    end function
  end interface
  
  interface select_forces
    module function select_forces_ComplexUnivariate(input,forces) &
       & result(output) 
      type(ComplexUnivariate),  intent(in)  :: input
      type(ComplexSingleForce), intent(in)  :: forces(:)
      type(ComplexSingleForce), allocatable :: output(:)
    end function
  
    module function select_forces_ComplexUnivariates(input,forces) &
       & result(output) 
      type(ComplexUnivariate),  intent(in)  :: input(:)
      type(ComplexSingleForce), intent(in)  :: forces(:)
      type(ComplexSingleForce), allocatable :: output(:)
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Compares two monomials for equality up to coefficient.
    ! ----------------------------------------------------------------------
    module function compare_complex_monomials(this,that) result(output) 
      class(*), intent(in) :: this
      class(*), intent(in) :: that
      logical              :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_ComplexUnivariate(this,input) 
      class(ComplexUnivariate), intent(out) :: this
      type(String),             intent(in) :: input
    end subroutine
  end interface
  
  interface
    module function write_ComplexUnivariate(this) result(output) 
      class(ComplexUnivariate), intent(in) :: this
      type(String)                         :: output
    end function
  end interface
  
  interface ComplexUnivariate
    impure elemental module function new_ComplexUnivariate_String(input) &
       & result(this) 
      type(String), intent(in) :: input
      type(ComplexUnivariate)  :: this
    end function
  end interface
  
  interface
    module subroutine read_ComplexMonomial(this,input) 
      class(ComplexMonomial), intent(out) :: this
      type(String),           intent(in)  :: input
    end subroutine
  end interface
  
  interface
    module function write_ComplexMonomial(this) result(output) 
      class(ComplexMonomial), intent(in) :: this
      type(String)                       :: output
    end function
  end interface
  
  interface ComplexMonomial
    impure elemental module function new_ComplexMonomial_String(input) &
       & result(this) 
      type(String), intent(in) :: input
      type(ComplexMonomial)    :: this
    end function
  end interface
  
  interface
    module subroutine read_ComplexPolynomial(this,input) 
      class(ComplexPolynomial), intent(out) :: this
      type(String),             intent(in)  :: input
    end subroutine
  end interface
  
  interface
    module function write_ComplexPolynomial(this) result(output) 
      class(ComplexPolynomial), intent(in) :: this
      type(String)                         :: output
    end function
  end interface
  
  interface ComplexPolynomial
    impure elemental module function new_ComplexPolynomial_String(input) &
       & result(this) 
      type(String), intent(in) :: input
      type(ComplexPolynomial)  :: this
    end function
  end interface
end module
