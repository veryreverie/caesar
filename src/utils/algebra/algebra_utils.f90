!> Algebra-related utilities.
module caesar_algebra_utils_module
  use caesar_foundations_module
  use caesar_io_module
  
  use caesar_mathematical_constants_module
  use caesar_linear_algebra_module
  implicit none
  
  private
  
  public :: gcd
  public :: lcm
  public :: triple_product
  public :: l2_norm
  public :: sum_squares
  public :: exp_2pii
  public :: cos_2pi
  public :: sin_2pi
  public :: factorial
  public :: real_factorial
  public :: log_factorial
  public :: odd_factorial
  public :: real_odd_factorial
  public :: log_odd_factorial
  public :: binomial
  public :: real_binomial
  public :: log_binomial
  public :: multinomial
  public :: real_multinomial
  public :: log_multinomial
  public :: int_sqrt
  public :: bound
  
  interface l2_norm
    !> Returns the L2 norm of the vector `input`.
    module function l2_norm_reals(input) result(output) 
      real(dp), intent(in) :: input(:)
      real(dp)             :: output
    end function

    !> Returns the L2 norm of the vector `input`.
    module function l2_norm_complexes(input) result(output) 
      complex(dp), intent(in) :: input(:)
      real(dp)                :: output
    end function
  end interface
  
  interface sum_squares
    !> Returns the sum of the norm squared elements of `input`.
    impure elemental module function sum_squares_RealVector(input) &
       & result(output) 
      type(RealVector), intent(in) :: input
      real(dp)                     :: output
    end function

    !> Returns the sum of the norm squared elements of `input`.
    impure elemental module function sum_squares_RealMatrix(input) &
       & result(output) 
      type(RealMatrix), intent(in) :: input
      real(dp)                     :: output
    end function

    !> Returns the sum of the norm squared elements of `input`.
    impure elemental module function sum_squares_ComplexVector(input) &
       & result(output) 
      type(ComplexVector), intent(in) :: input
      real(dp)                        :: output
    end function

    !> Returns the sum of the norm squared elements of `input`.
    impure elemental module function sum_squares_ComplexMatrix(input) &
       & result(output) 
      type(ComplexMatrix), intent(in) :: input
      real(dp)                        :: output
    end function
  end interface
  
  interface triple_product
    !> Returns the triple product of three vectors, `a^b.c`.
    !> The input vectors must be three-dimensional.
    module function triple_product_IntVector(a,b,c) result(output) 
      type(IntVector) :: a
      type(IntVector) :: b
      type(IntVector) :: c
      integer         :: output
    end function

    !> Returns the triple product of three vectors, `a^b.c`.
    !> The input vectors must be three-dimensional.
    module function triple_product_RealVector(a,b,c) result(output) 
      type(RealVector) :: a
      type(RealVector) :: b
      type(RealVector) :: c
      real(dp)         :: output
    end function
  end interface
  
  interface
    !> Returns the factorial of `input`.
    impure elemental module function factorial(input) result(output) 
      integer, intent(in) :: input
      integer             :: output
    end function
    
    !> Returns the factorial of `input`.
    impure elemental module function real_factorial(input) result(output) 
      integer, intent(in) :: input
      real(dp)            :: output
    end function
  
    !> Returns the natural log of the factorial of `input`, `ln(n!)`.
    impure elemental module function log_factorial(input) result(output) 
      integer, intent(in) :: input
      real(dp)            :: output
    end function
  
    !> Returns the product of the first `input` odd numbers.
    !> `odd_factorial(n)` = prod_{k=1}^n[ 2*k-1 ] = (2n-1)!/(2^n * (n-1)!).
    impure elemental module function odd_factorial(input) result(output) 
      integer, intent(in) :: input
      integer             :: output
    end function
  
    !> Returns the product of the first `input` odd numbers.
    !> `real_odd_factorial(n)` = prod_{k=1}^n[ 2*k-1 ]
    !>                         = (2n-1)!/(2^n * (n-1)!).
    impure elemental module function real_odd_factorial(input) result(output) 
      integer, intent(in) :: input
      real(dp)            :: output
    end function

    !> Returns the natural log of [[odd_factorial]] of `input`.
    impure elemental module function log_odd_factorial(input) result(output) 
      integer, intent(in) :: input
      real(dp)            :: output
    end function
  
    !> Returns the binomial coefficient of `(top bottom)`.
    !> `binomial(a,b)` = a!/(b!*(a-b)!)
    impure elemental module function binomial(top,bottom) result(output) 
      integer, intent(in) :: top
      integer, intent(in) :: bottom
      integer             :: output
    end function

    !> Returns the binomial coefficient of `(top bottom)`.
    !> `real_binomial(a,b)` = a!/(b!*(a-b)!)
    impure elemental module function real_binomial(top,bottom) result(output) 
      integer, intent(in) :: top
      integer, intent(in) :: bottom
      real(dp)            :: output
    end function

    !> Returns the natural log of the [[binomial]] of `(top,bottom)`.
    impure elemental module function log_binomial(top,bottom) result(output) 
      integer, intent(in) :: top
      integer, intent(in) :: bottom
      real(dp)            :: output
    end function

    !> Returns the multinomial coefficient of a numerator and a set of
    !>    denominators.
    !> `multinomial(a, [b,c,d,...])` = a!/(b!c!d!...).
    module function multinomial(top,bottom) result(output) 
      integer, intent(in) :: top
      integer, intent(in) :: bottom(:)
      integer             :: output
    end function

    !> Returns the multinomial coefficient of a numerator and a set of
    !>    denominators.
    !> `real_multinomial(a, [b,c,d,...])` = a!/(b!c!d!...).
    module function real_multinomial(top,bottom) result(output) 
      integer, intent(in) :: top
      integer, intent(in) :: bottom(:)
      real(dp)            :: output
    end function

    !> Returns the natural log of the [[monomial]] of `(top,bottom)`.
    module function log_multinomial(top,bottom) result(output) 
      integer, intent(in) :: top
      integer, intent(in) :: bottom(:)
      real(dp)            :: output
    end function

    !> Returns the Square-root of an integer `input`.
    !> Throws an error if `input` is not square.
    module function int_sqrt(input) result(output) 
      integer, intent(in) :: input
      integer             :: output
    end function
  end interface
  
  ! Greatest common denominator.
  interface gcd
    !> Returns the greatest common divisor of the moduli of two integers.
    !> N.B. `gcd(a,0)=gcd(0,a)=a`.
    module function gcd_2(int_1,int_2) result(output) 
      integer, intent(in) :: int_1
      integer, intent(in) :: int_2
      integer             :: output
    end function

    !> Returns the greatest common divisor of
    !>    the moduli of three integers.
    module function gcd_3(int_1,int_2,int_3) result(output) 
      integer, intent(in) :: int_1
      integer, intent(in) :: int_2
      integer, intent(in) :: int_3
      integer             :: output
    end function

    !> Returns the greatest common divisor of the moduli of four integers.
    module function gcd_4(int_1,int_2,int_3,int_4) result(output) 
      integer, intent(in) :: int_1
      integer, intent(in) :: int_2
      integer, intent(in) :: int_3
      integer, intent(in) :: int_4
      integer             :: output
    end function

    !> Returns the greatest common divisor of
    !>    the moduli of an array of integers.
    !> N.B. `gcd([])=0` and `gcd([a])=a` to preserve associativity.
    recursive module function gcd_integers(input) result(output) 
      integer, intent(in) :: input(:)
      integer             :: output
    end function
  end interface
  
  interface lcm
    !> Returns the lowest common multiple of two integers.
    !> N.B. `lcm(a,0)=lcm(0,a)=0`.
    module function lcm_2(int_1,int_2) result(output) 
      integer, intent(in) :: int_1
      integer, intent(in) :: int_2
      integer             :: output
    end function

    !> Returns the lowest common multiple of three integers.
    module function lcm_3(int_1,int_2,int_3) result(output) 
      integer, intent(in) :: int_1
      integer, intent(in) :: int_2
      integer, intent(in) :: int_3
      integer             :: output
    end function

    !> Returns the lowest common multiple of four integers.
    module function lcm_4(int_1,int_2,int_3,int_4) result(output) 
      integer, intent(in) :: int_1
      integer, intent(in) :: int_2
      integer, intent(in) :: int_3
      integer, intent(in) :: int_4
      integer             :: output
    end function

    !> Returns the lowest common multiple of an array of integers.
    !> N.B. `lcm([])=1` and `lcm([a])=a` to preserve associativity.
    recursive module function lcm_integers(input) result(output) 
      integer, intent(in) :: input(:)
      integer             :: output
    end function
  end interface
  
  interface exp_2pii
    !> Returns `exp(2*pi*i*input)`.
    impure elemental module function exp_2pii_real(input) result(output) 
      real(dp), intent(in) :: input
      complex(dp)          :: output
    end function
  end interface
  
  interface cos_2pi
    !> Returns `cos(2*pi*input)`.
    impure elemental module function cos_2pi_real(input) result(output) 
      real(dp), intent(in) :: input
      real(dp)             :: output
    end function
  end interface
  
  interface sin_2pi
    !> Returns `sin(2*pi*input)`.
    impure elemental module function sin_2pi_real(input) result(output) 
      real(dp), intent(in) :: input
      real(dp)             :: output
    end function
  end interface
  
  interface
    !> Returns `input` if `input` is between `lower_bound` and `upper_bound`,
    !>    or the relevant bound if not.
    !> `bound(a,b,c) = max(a,min(b,c))`.
    impure elemental module function bound(input,lower_bound,upper_bound) &
       & result(output) 
      real(dp), intent(in) :: input
      real(dp), intent(in) :: lower_bound
      real(dp), intent(in) :: upper_bound
      real(dp)             :: output
    end function
  end interface
end module
