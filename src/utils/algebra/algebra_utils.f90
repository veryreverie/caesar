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
  public :: int_sqrt
  
  ! Greatest common denominator.
  interface gcd
    module procedure gcd_2
    module procedure gcd_3
    module procedure gcd_4
  end interface
  
  ! Lowest common multiple.
  interface lcm
    module procedure lcm_2
    module procedure lcm_3
    module procedure lcm_4
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
  integer, allocatable, save :: lookup(:)
  integer, allocatable       :: temp(:)
  
  integer :: i,ialloc
  
  ! Initialise lookup the first time factorial is called.
  if (.not. allocated(lookup)) then
    lookup = [1]
  endif
  
  ! If the requested factorial has not previously been calculated,
  !    extends lookup to include it.
  if (input>=size(lookup)) then
    temp = lookup
    deallocate(lookup, stat=ialloc); call err(ialloc)
    allocate(lookup(input+1), stat=ialloc); call err(ialloc)
    lookup(:size(temp)) = temp
    do i=size(temp)+1,size(lookup)
      lookup(i) = (i-1)*lookup(i-1)
    enddo
  endif
  
  ! lookup(1) = 0!, so lookup(n+1) = n!.
  output = lookup(input+1)
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

! ----------------------------------------------------------------------
! Calculate the lowest common multiple of two non-negative integers.
! lcm(a,b) = a*b/gcd(a,b)
! ----------------------------------------------------------------------
function lcm_2(int_1,int_2) result(output)
  implicit none
  
  integer, intent(in) :: int_1
  integer, intent(in) :: int_2
  integer             :: output
  
  output = int_1*int_2/gcd(int_1,int_2)
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
