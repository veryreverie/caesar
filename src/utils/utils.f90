! ======================================================================
! various utilities
! ======================================================================
module utils_module
  use constants_module, only : dp
  use string_module
  use io_module
  
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
  
  ! The sum of the elements squared of a matrix.
  interface sum_squares
    module procedure sum_squares_RealVector
    module procedure sum_squares_RealMatrix
    module procedure sum_squares_ComplexVector
    module procedure sum_squares_ComplexMatrix
  end interface
contains

! ----------------------------------------------------------------------
! Returns an array containing the command line arguments.
! ----------------------------------------------------------------------
function command_line_args() result(args)
  implicit none

  type(String), allocatable :: args(:)
  
  ! Temporary variables.
  integer         :: i,ialloc
  integer         :: arg_count
  character(1000) :: temp_char

  arg_count = command_argument_count()
  allocate (args(arg_count+1), stat=ialloc); call err(ialloc)
  do i=0,arg_count
    call get_command_argument(i, temp_char)
    args(i+1) = trim(temp_char)
  enddo
end function

! ----------------------------------------------------------------------
! Make a directory, if it doesn't already exist.
! ----------------------------------------------------------------------
subroutine mkdir(dirname)
  implicit none
  
  type(String), intent(in) :: dirname
  
  integer :: result_code
  
  result_code = system_call( &
     & 'if [ ! -e '//dirname//' ]; then mkdir '//dirname//'; fi')
  if (result_code/=0) then
    call print_line('Error: failed to make directory: '//dirname)
    call err()
  endif
end subroutine

! ----------------------------------------------------------------------
! Vector L2 norm.
! Replicates norm2() from f2008 standard.
! ----------------------------------------------------------------------
pure function l2_norm(input) result(output)
  implicit none
  
  real(dp), intent(in) :: input(:)
  real(dp)             :: output
  
  output = sqrt(dot_product(input,input))
end function

! ----------------------------------------------------------------------
! The sum of the squares of the elements of a matrix.
! ----------------------------------------------------------------------
function sum_squares_RealVector(input) result(output)
  use linear_algebra_module
  implicit none
  
  type(RealVector), intent(in) :: input
  real(dp)                     :: output
  
  output = sum(dble(input)*dble(input))
end function

function sum_squares_RealMatrix(input) result(output)
  use linear_algebra_module
  implicit none
  
  type(RealMatrix), intent(in) :: input
  real(dp)                     :: output
  
  output = sum(dble(input)*dble(input))
end function

function sum_squares_ComplexVector(input) result(output)
  use linear_algebra_module
  implicit none
  
  type(ComplexVector), intent(in) :: input
  real(dp)                        :: output
  
  output = sum(real(cmplx(input)*conjg(cmplx(input))))
end function

function sum_squares_ComplexMatrix(input) result(output)
  use linear_algebra_module
  implicit none
  
  type(ComplexMatrix), intent(in) :: input
  real(dp)                        :: output
  
  output = sum(real(cmplx(input)*conjg(cmplx(input))))
end function

! ----------------------------------------------------------------------
! Vector outer product.
! ----------------------------------------------------------------------
function outer_product(input1,input2) result(output)
  implicit none
  
  real(dp), intent(in)  :: input1(:)
  real(dp), intent(in)  :: input2(:)
  real(dp), allocatable :: output(:,:)
  
  ! Temporary variables.
  integer :: i,ialloc
  
  allocate(output(size(input1), size(input2)), stat=ialloc); call err(ialloc)
  do i=1,size(input2)
    output(:,i) = input2(i) * input1
  enddo
end function

! ----------------------------------------------------------------------
! Vector triple product.
! ----------------------------------------------------------------------
function triple_product_IntVector(a,b,c) result(output)
  use linear_algebra_module
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
  use linear_algebra_module
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
! a=Ac, b=Bc, A>=B. c is the gcd of a and b.
! gcd(a,b) = gcd(Ac,Bc) = gcd((A-nB)c,Bc).
recursive function gcd_helper(a,b) result(output)
  implicit none
  
  integer, intent(in) :: a
  integer, intent(in) :: b
  integer             :: output
  
  if (b==0) then
    output = a
  else
    output = gcd_helper(b,modulo(a,b))
  endif
end function

function gcd_2(int_1,int_2) result(output)
  implicit none
  
  integer, intent(in) :: int_1
  integer, intent(in) :: int_2
  integer             :: output
  
  output = gcd_helper( max(abs(int_1),abs(int_2)), &
                     & min(abs(int_1),abs(int_2))  )
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
end module
