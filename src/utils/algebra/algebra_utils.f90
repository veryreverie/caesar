! ======================================================================
! various utilities
! ======================================================================
module algebra_utils_submodule
  use precision_module
  use io_module
  
  use mathematical_constants_submodule
  use linear_algebra_submodule
  implicit none
  
  private
  
  public :: gcd
  public :: lcm
  public :: triple_product
  public :: l2_norm
  public :: sum_squares
  public :: exp_2pii
  public :: factorial
  public :: binomial
  public :: int_sqrt
  
  ! Greatest common denominator.
  interface gcd
    module procedure gcd_2
    module procedure gcd_3
    module procedure gcd_4
    module procedure gcd_integers
  end interface
  
  ! Lowest common multiple.
  interface lcm
    module procedure lcm_2
    module procedure lcm_3
    module procedure lcm_4
    module procedure lcm_integers
  end interface
  
  ! The triple product of three vectors.
  interface triple_product
    module procedure triple_product_IntVector
    module procedure triple_product_RealVector
  end interface
  
  ! The L2 norm of an array of reals.
  interface l2_norm
    module procedure l2_norm_reals
  end interface
  
  ! The sum of the elements squared of a matrix.
  interface sum_squares
    module procedure sum_squares_RealVector
    module procedure sum_squares_RealMatrix
    module procedure sum_squares_ComplexVector
    module procedure sum_squares_ComplexMatrix
  end interface
  
  ! exp(2*pi*i*input).
  interface exp_2pii
    module procedure exp_2pii_real
    module procedure exp_2pii_integer
  end interface
contains

! ----------------------------------------------------------------------
! Vector L2 norm.
! Replicates norm2() from f2008 standard.
! ----------------------------------------------------------------------
function l2_norm_reals(input) result(output)
  implicit none
  
  real(dp), intent(in) :: input(:)
  real(dp)             :: output
  
  output = sqrt(dot_product(input,input))
end function

! ----------------------------------------------------------------------
! The sum of the squares of the elements of a matrix.
! ----------------------------------------------------------------------
impure elemental function sum_squares_RealVector(input) result(output)
  implicit none
  
  type(RealVector), intent(in) :: input
  real(dp)                     :: output
  
  output = sum(dble(input)*dble(input))
end function

impure elemental function sum_squares_RealMatrix(input) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: input
  real(dp)                     :: output
  
  output = sum(dble(input)*dble(input))
end function

impure elemental function sum_squares_ComplexVector(input) result(output)
  implicit none
  
  type(ComplexVector), intent(in) :: input
  real(dp)                        :: output
  
  output = sum(real(cmplx(input)*conjg(cmplx(input))))
end function

impure elemental function sum_squares_ComplexMatrix(input) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: input
  real(dp)                        :: output
  
  output = sum(real(cmplx(input)*conjg(cmplx(input))))
end function

! ----------------------------------------------------------------------
! Vector triple product.
! ----------------------------------------------------------------------
function triple_product_IntVector(a,b,c) result(output)
  implicit none
  
  type(IntVector) :: a
  type(IntVector) :: b
  type(IntVector) :: c
  integer         :: output
  
  integer :: matrix(3,3)
  
  matrix(:,1) = int(a)
  matrix(:,2) = int(b)
  matrix(:,3) = int(c)
  
  output = determinant(matrix)
end function

function triple_product_RealVector(a,b,c) result(output)
  implicit none
  
  type(RealVector) :: a
  type(RealVector) :: b
  type(RealVector) :: c
  real(dp)         :: output
  
  real(dp) :: matrix(3,3)
  
  matrix(:,1) = dble(a)
  matrix(:,2) = dble(b)
  matrix(:,3) = dble(c)
  
  output = determinant(matrix)
end function

! ----------------------------------------------------------------------
! Factorial.
! Stores calculated results to avoid excess computation.
! ----------------------------------------------------------------------
function factorial(input) result(output)
  implicit none
  
  integer, intent(in) :: input
  integer             :: output
  
  ! Lookup is saved between calls of this function.
  integer, allocatable, save :: factorials(:)
  integer, allocatable       :: temp(:)
  
  integer :: i,j,ialloc
  
  ! The first time this function is called,
  !    calculate all storable factorials, {i!} s.t. i!<=huge(0).
  if (.not. allocated(factorials)) then
    ! Find the largest integer i s.t. i! is too large to be
    !    stored as an integer.
    i = 0
    j = 1
    do
      i = i+1
      if (j>huge(0)/i) then
        allocate(factorials(i), stat=ialloc); call err(ialloc)
        exit
      endif
      j = j*i ! j = i!.
    enddo
    
    ! Calculate and store all storable factorials.
    factorials(1) = 1
    do i=1,size(factorials)-1
      factorials(i+1) = factorials(i)*i ! factorials(i+1) = i!.
    enddo
  endif
  
  ! Check that the input integer, n, obeys n>=0,
  !    and small enough that n! is storable as an integer.
  if (input<0) then
    call print_line( ERROR//': Trying to calculate the factorial of an &
       &integer less than 0.')
    call err()
  elseif (input>=size(factorials)) then
    call print_line( ERROR//': Trying to calculate the factorial of an &
       &integer greater than '//size(factorials)-1//'.')
    call err()
  endif
  
  ! factorials(1) = 0!, so factorials(n+1) = n!.
  output = factorials(input+1)
end function

! ----------------------------------------------------------------------
! Calculates the binomial coefficient of two integers.
! binomial(a,b) = a!/(b!*(a-b)!)
! ----------------------------------------------------------------------
function binomial(top,bottom) result(output)
  implicit none
  
  integer, intent(in) :: top
  integer, intent(in) :: bottom
  integer             :: output
  
  output = factorial(top) / (factorial(bottom)*factorial(top-bottom))
end function

! ----------------------------------------------------------------------
! Square-root of an integer.
! Returns an error if the integer is not square.
! ----------------------------------------------------------------------
function int_sqrt(input) result(output)
  implicit none
  
  integer, intent(in) :: input
  integer             :: output
  
  output = nint(sqrt(real(input,dp)))
  if (output*output/=input) then
    call err()
  endif
end function

! ----------------------------------------------------------------------
! Calculate the greatest common divisor of two non-negative integers using
! Euclid's algorithm.
! ----------------------------------------------------------------------
! For convenience, it is defined that
!    - gcd(0,b)=b, and gcd(a,0)=a.
!    - gcd(-a,b)=gcd(a,-b)=gcd(-a,-b)=gcd(a,b).
! a=Ac, b=Bc, A<=B. c is the gcd of a and b.
! gcd(a,b) = gcd(Ac,Bc) = gcd(Ac,(B-nA)c).
function gcd_2(int_1,int_2) result(output)
  implicit none
  
  integer, intent(in) :: int_1
  integer, intent(in) :: int_2
  integer             :: output
  
  integer :: a_b(2) ! a_b(1)=a, a_b(2)=b.
  
  a_b = [ min(abs(int_1), abs(int_2)), &
        & max(abs(int_1), abs(int_2))  ]
  
  ! Loop until the g.c.d. is found, obeying the following:
  !    gcd(0,b) = b
  !    gcd(a,b) = gcd(modulo(b,a), a)
  do while (a_b(1)/=0)
    a_b = [modulo(a_b(2),a_b(1)), a_b(1)]
  enddo
  
  output = a_b(2)
end function

! Calculate the gcd of more than two integers, by recursively finding
!    pairwise gcds.
function gcd_3(int_1,int_2,int_3) result(output)
  implicit none
  
  integer, intent(in) :: int_1
  integer, intent(in) :: int_2
  integer, intent(in) :: int_3
  integer             :: output
  
  output = gcd(gcd(int_1,int_2),int_3)
end function

function gcd_4(int_1,int_2,int_3,int_4) result(output)
  implicit none
  
  integer, intent(in) :: int_1
  integer, intent(in) :: int_2
  integer, intent(in) :: int_3
  integer, intent(in) :: int_4
  integer             :: output
  
  output = gcd(gcd(int_1,int_2,int_3),int_4)
end function

! Calculate the gcd of an array of integers, by recusively finding
!    pairwise gcds.
! For convenience, gcd([a])=a if a/=0, and gcd([0])=1.
recursive function gcd_integers(input) result(output)
  implicit none
  
  integer, intent(in) :: input(:)
  integer             :: output
  
  if (size(input)==0) then
    call print_line(CODE_ERROR//': called gcd on an empty array.')
    call err()
  elseif (size(input)==1) then
    if (input(1)==0) then
      output = 1
    else
      output = abs(input(1))
    endif
  else
    output = gcd([gcd(input(1),input(2)), input(3:)])
  endif
end function

! ----------------------------------------------------------------------
! Calculate the lowest common multiple of two integers.
! lcm(a,b) = |a*b|/gcd(a,b)
! ----------------------------------------------------------------------
function lcm_2(int_1,int_2) result(output)
  implicit none
  
  integer, intent(in) :: int_1
  integer, intent(in) :: int_2
  integer             :: output
  
  output = abs(int_1*int_2)/gcd(int_1,int_2)
end function

! Calculate the lcm of more than two integers, by recursively finding
!    pairwise lcms.
function lcm_3(int_1,int_2,int_3) result(output)
  implicit none
  
  integer, intent(in) :: int_1
  integer, intent(in) :: int_2
  integer, intent(in) :: int_3
  integer             :: output
  
  output = lcm(lcm(int_1,int_2),int_3)
end function

function lcm_4(int_1,int_2,int_3,int_4) result(output)
  implicit none
  
  integer, intent(in) :: int_1
  integer, intent(in) :: int_2
  integer, intent(in) :: int_3
  integer, intent(in) :: int_4
  integer             :: output
  
  output = lcm(lcm(int_1,int_2,int_3),int_4)
end function

! Calculate the lcm of an array of integers, by recusively finding
!    pairwise lcms.
! For convenience, lcm([a])=a if |a|.
recursive function lcm_integers(input) result(output)
  implicit none
  
  integer, intent(in) :: input(:)
  integer             :: output
  
  if (size(input)==0) then
    call print_line(CODE_ERROR//': called lcm on an empty array.')
    call err()
  elseif (size(input)==1) then
    output = abs(input(1))
  else
    output = lcm([lcm(input(1),input(2)), input(3:)])
  endif
end function

! ----------------------------------------------------------------------
! Returns exp(2*pi*i*input).
! ----------------------------------------------------------------------
impure elemental function exp_2pii_real(input) result(output)
  implicit none
  
  real(dp), intent(in) :: input
  complex(dp)          :: output
  
  real(dp) :: exponent
  
  exponent = 2*PI*input
  output = cmplx(cos(exponent),sin(exponent),dp)
end function

impure elemental function exp_2pii_integer(input) result(output)
  implicit none
  
  integer, intent(in) :: input
  complex(dp)         :: output
  
  output = exp_2pii(real(input,dp))
end function
end module
