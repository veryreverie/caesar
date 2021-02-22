! Algebra-related utilities.
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
  
  ! The L2 norm of an array of reals.
  interface l2_norm
    ! ----------------------------------------------------------------------
    ! Vector L2 norm.
    ! Replicates norm2() from f2008 standard.
    ! ----------------------------------------------------------------------
    module function l2_norm_reals(input) result(output) 
      real(dp), intent(in) :: input(:)
      real(dp)             :: output
    end function

    module function l2_norm_complexes(input) result(output) 
      complex(dp), intent(in) :: input(:)
      real(dp)                :: output
    end function
  end interface
  
  ! The sum of the elements squared of a matrix.
  interface sum_squares
    ! ----------------------------------------------------------------------
    ! The sum of the squares of the elements of a matrix.
    ! ----------------------------------------------------------------------
    impure elemental module function sum_squares_RealVector(input) &
       & result(output) 
      type(RealVector), intent(in) :: input
      real(dp)                     :: output
    end function

    impure elemental module function sum_squares_RealMatrix(input) &
       & result(output) 
      type(RealMatrix), intent(in) :: input
      real(dp)                     :: output
    end function

    impure elemental module function sum_squares_ComplexVector(input) &
       & result(output) 
      type(ComplexVector), intent(in) :: input
      real(dp)                        :: output
    end function

    impure elemental module function sum_squares_ComplexMatrix(input) &
       & result(output) 
      type(ComplexMatrix), intent(in) :: input
      real(dp)                        :: output
    end function
  end interface
  
  ! The triple product of three vectors.
  interface triple_product
    ! ----------------------------------------------------------------------
    ! Vector triple product.
    ! ----------------------------------------------------------------------
    module function triple_product_IntVector(a,b,c) result(output) 
      type(IntVector) :: a
      type(IntVector) :: b
      type(IntVector) :: c
      integer         :: output
    end function

    module function triple_product_RealVector(a,b,c) result(output) 
      type(RealVector) :: a
      type(RealVector) :: b
      type(RealVector) :: c
      real(dp)         :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Factorial. factorial(n)=n! = prod_{k=1}^n[ k ].
    ! Stores calculated results to avoid excess computation.
    ! ----------------------------------------------------------------------
    impure elemental module function factorial(input) result(output) 
      integer, intent(in) :: input
      integer             :: output
    end function

    module subroutine calculate_factorials(factorials) 
      integer, intent(inout), allocatable :: factorials(:)
    end subroutine
  
    ! ----------------------------------------------------------------------
    ! Factorial, but giving the answer as a real.
    ! Stores calculated results to avoid excess computation.
    ! ----------------------------------------------------------------------
    impure elemental module function real_factorial(input) result(output) 
      integer, intent(in) :: input
      real(dp)            :: output
    end function

    module subroutine calculate_real_factorials(real_factorials,input) 
      real(dp), intent(inout), allocatable :: real_factorials(:)
      integer,  intent(in)                 :: input
    end subroutine
  
    ! ----------------------------------------------------------------------
    ! Calculates ln(n!), the log of the factorial of the input.
    ! ----------------------------------------------------------------------
    impure elemental module function log_factorial(input) result(output) 
      integer, intent(in) :: input
      real(dp)            :: output
    end function

    module subroutine calculate_log_factorials(log_factorials, &
       & log_factorials_calculated,input) 
      real(dp), intent(inout), allocatable :: log_factorials(:)
      integer,  intent(inout)              :: log_factorials_calculated
      integer,  intent(in)                 :: input
    end subroutine
  
    ! ----------------------------------------------------------------------
    ! The product of the first n odd numbers, rather than the first n numbers.
    ! odd_factorial(n) = prod_{k=1}^n[ 2*k-1 ] = (2n-1)!/(2^n * (n-1)!).
    ! Stores calculated results to avoid excess computation.
    ! ----------------------------------------------------------------------
    impure elemental module function odd_factorial(input) result(output) 
      integer, intent(in) :: input
      integer             :: output
    end function

    module subroutine calculate_odd_factorials(odd_factorials) 
      integer, intent(inout), allocatable :: odd_factorials(:)
    end subroutine
  
    ! ----------------------------------------------------------------------
    ! Odd factorial, but giving the answer as a real.
    ! Stores calculated results to avoid excess computation.
    ! ----------------------------------------------------------------------
    impure elemental module function real_odd_factorial(input) result(output) 
      integer, intent(in) :: input
      real(dp)            :: output
    end function

    module subroutine calculate_real_odd_factorials(real_odd_factorials,input) 
      real(dp), intent(inout), allocatable :: real_odd_factorials(:)
      integer,  intent(in)                 :: input
    end subroutine

    ! ----------------------------------------------------------------------
    ! Calculates ln(real_odd_factorial(n)),
    !    the log of the odd factorial of the input.
    ! ----------------------------------------------------------------------
    impure elemental module function log_odd_factorial(input) result(output) 
      integer, intent(in) :: input
      real(dp)            :: output
    end function

    module subroutine calculate_log_odd_factorials(log_odd_factorials, &
       & log_odd_factorials_calculated,input) 
      real(dp), intent(inout), allocatable :: log_odd_factorials(:)
      integer,  intent(inout)              :: log_odd_factorials_calculated
      integer,  intent(in)                 :: input
    end subroutine
  
    ! ----------------------------------------------------------------------
    ! Calculates the binomial coefficient of two integers.
    ! binomial(a,b) = a!/(b!*(a-b)!)
    ! ----------------------------------------------------------------------
    impure elemental module function binomial(top,bottom) result(output) 
      integer, intent(in) :: top
      integer, intent(in) :: bottom
      integer             :: output
    end function

    module subroutine calculate_binomials(binomials,binomials_calculated, &
       & top,bottom) 
      integer, intent(inout), allocatable :: binomials(:,:)
      integer, intent(inout), allocatable :: binomials_calculated(:)
      integer, intent(in)                 :: top
      integer, intent(in)                 :: bottom
    end subroutine

    impure elemental module function real_binomial(top,bottom) result(output) 
      integer, intent(in) :: top
      integer, intent(in) :: bottom
      real(dp)            :: output
    end function

    module subroutine calculate_real_binomials(real_binomials, &
       & real_binomials_calculated,top,bottom) 
      real(dp), intent(inout), allocatable :: real_binomials(:,:)
      integer,  intent(inout), allocatable :: real_binomials_calculated(:)
      integer,  intent(in)                 :: top
      integer,  intent(in)                 :: bottom
    end subroutine

    ! ----------------------------------------------------------------------
    ! Calculates ln(bin(top,bottom)) = ln(top!/(bottom!(top-bottom)!)),
    !    the log of the binomial of two integers.
    ! ----------------------------------------------------------------------
    impure elemental module function log_binomial(top,bottom) result(output) 
      integer, intent(in) :: top
      integer, intent(in) :: bottom
      real(dp)            :: output
    end function

    ! ----------------------------------------------------------------------
    ! Calculates the multinomial coefficient of a numerator and a set of
    !    denominators.
    ! multinomial(a,[b,c,d,...]) = a!/(b!c!d!...)
    ! ----------------------------------------------------------------------
    module function multinomial(top,bottom) result(output) 
      integer, intent(in) :: top
      integer, intent(in) :: bottom(:)
      integer             :: output
    end function

    module function real_multinomial(top,bottom) result(output) 
      integer, intent(in) :: top
      integer, intent(in) :: bottom(:)
      real(dp)            :: output
    end function

    ! ----------------------------------------------------------------------
    ! Calculates ln(multinomial(top,bottom)) = ln(top!/product(bottom!)),
    !    the log of the multinomial of a numerator and set of denominators.
    ! ----------------------------------------------------------------------
    module function log_multinomial(top,bottom) result(output) 
      integer, intent(in) :: top
      integer, intent(in) :: bottom(:)
      real(dp)            :: output
    end function

    ! ----------------------------------------------------------------------
    ! Square-root of an integer.
    ! Returns an error if the integer is not square.
    ! ----------------------------------------------------------------------
    module function int_sqrt(input) result(output) 
      integer, intent(in) :: input
      integer             :: output
    end function
  end interface
  
  ! Greatest common denominator.
  interface gcd
    ! ----------------------------------------------------------------------
    ! Calculate the greatest common divisor of two non-negative integers using
    ! Euclid's algorithm.
    ! ----------------------------------------------------------------------
    ! For convenience, it is defined that
    !    - gcd(0,b)=b, and gcd(a,0)=a.
    !    - gcd(-a,b)=gcd(a,-b)=gcd(-a,-b)=gcd(a,b).
    ! a=Ac, b=Bc, A<=B. c is the gcd of a and b.
    ! gcd(a,b) = gcd(Ac,Bc) = gcd(Ac,(B-nA)c).
    module function gcd_2(int_1,int_2) result(output) 
      integer, intent(in) :: int_1
      integer, intent(in) :: int_2
      integer             :: output
    end function

    ! Calculate the gcd of more than two integers, by recursively finding
    !    pairwise gcds.
    module function gcd_3(int_1,int_2,int_3) result(output) 
      integer, intent(in) :: int_1
      integer, intent(in) :: int_2
      integer, intent(in) :: int_3
      integer             :: output
    end function

    module function gcd_4(int_1,int_2,int_3,int_4) result(output) 
      integer, intent(in) :: int_1
      integer, intent(in) :: int_2
      integer, intent(in) :: int_3
      integer, intent(in) :: int_4
      integer             :: output
    end function

    ! Calculate the gcd of an array of integers, by recusively finding
    !    pairwise gcds.
    ! For convenience, gcd([a])=a if a/=0, and gcd([0])=1.
    recursive module function gcd_integers(input) result(output) 
      integer, intent(in) :: input(:)
      integer             :: output
    end function
  end interface
  
  ! Lowest common multiple.
  interface lcm
    ! ----------------------------------------------------------------------
    ! Calculate the lowest common multiple of two integers.
    ! lcm(a,b) = |a*b|/gcd(a,b)
    ! ----------------------------------------------------------------------
    module function lcm_2(int_1,int_2) result(output) 
      integer, intent(in) :: int_1
      integer, intent(in) :: int_2
      integer             :: output
    end function

    ! Calculate the lcm of more than two integers, by recursively finding
    !    pairwise lcms.
    module function lcm_3(int_1,int_2,int_3) result(output) 
      integer, intent(in) :: int_1
      integer, intent(in) :: int_2
      integer, intent(in) :: int_3
      integer             :: output
    end function

    module function lcm_4(int_1,int_2,int_3,int_4) result(output) 
      integer, intent(in) :: int_1
      integer, intent(in) :: int_2
      integer, intent(in) :: int_3
      integer, intent(in) :: int_4
      integer             :: output
    end function

    ! Calculate the lcm of an array of integers, by recusively finding
    !    pairwise lcms.
    ! For convenience, lcm([a])=a if |a|.
    recursive module function lcm_integers(input) result(output) 
      integer, intent(in) :: input(:)
      integer             :: output
    end function
  end interface
  
  ! exp(2*pi*i*input), cos(2*pi*input) and sin(2*pi*input).
  interface exp_2pii
    ! ----------------------------------------------------------------------
    ! exp(2*pi*i*input), cos(2*pi*input) and sin(2*pi*input).
    ! ----------------------------------------------------------------------
    impure elemental module function exp_2pii_real(input) result(output) 
      real(dp), intent(in) :: input
      complex(dp)          :: output
    end function
  end interface
  
  interface cos_2pi
    impure elemental module function cos_2pi_real(input) result(output) 
      real(dp), intent(in) :: input
      real(dp)             :: output
    end function
  end interface
  
  interface sin_2pi
    impure elemental module function sin_2pi_real(input) result(output) 
      real(dp), intent(in) :: input
      real(dp)             :: output
    end function
  end interface
  
  interface
    ! bound(input,lower,upper) = input if lower < input < upper
    !                            lower if input < lower
    !                            upper if upper < input
    impure elemental module function bound(input,lower_bound,upper_bound) &
       & result(output) 
      real(dp), intent(in) :: input
      real(dp), intent(in) :: lower_bound
      real(dp), intent(in) :: upper_bound
      real(dp)             :: output
    end function
  end interface
end module
