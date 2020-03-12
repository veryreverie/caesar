! ======================================================================
! various utilities
! ======================================================================
module algebra_utils_module
  use precision_module
  use io_module
  
  use mathematical_constants_module
  use linear_algebra_module
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
    module procedure l2_norm_complexes
  end interface
  
  ! The sum of the elements squared of a matrix.
  interface sum_squares
    module procedure sum_squares_RealVector
    module procedure sum_squares_RealMatrix
    module procedure sum_squares_ComplexVector
    module procedure sum_squares_ComplexMatrix
  end interface
  
  ! exp(2*pi*i*input), cos(2*pi*input) and sin(2*pi*input).
  interface exp_2pii
    module procedure exp_2pii_real
  end interface
  
  interface cos_2pi
    module procedure cos_2pi_real
  end interface
  
  interface sin_2pi
    module procedure sin_2pi_real
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

function l2_norm_complexes(input) result(output)
  implicit none
  
  complex(dp), intent(in) :: input(:)
  real(dp)                :: output
  
  output = sqrt(abs(dot_product(input,input)))
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
! Factorial. factorial(n)=n! = prod_{k=1}^n[ k ].
! Stores calculated results to avoid excess computation.
! ----------------------------------------------------------------------
impure elemental function factorial(input) result(output)
  implicit none
  
  integer, intent(in) :: input
  integer             :: output
  
  ! An array of factorials, stored to avoid re-calculating.
  ! Calculated the first time a related function is called.
  integer, allocatable, save :: factorials(:)
  
  ! Check that n>=0.
  if (input<0) then
    call print_line( ERROR//': Trying to calculate the factorial of an &
       &integer less than 0.')
    call err()
  endif
  
  call calculate_factorials(factorials)
  
  ! Check that n is small enough that n! is storable as an integer.
  if (input>=size(factorials)) then
    call print_line( ERROR//': Trying to calculate the factorial of an &
       &integer greater than '//(size(factorials)-1)//'. This is too large &
       &to fit in the standard integer type.')
    call err()
  endif
  
  ! factorials(1) = 0!, so factorials(n+1) = n!.
  output = factorials(input+1)
end function

subroutine calculate_factorials(factorials)
  implicit none
  
  integer, intent(inout), allocatable :: factorials(:)
  
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
end subroutine

! ----------------------------------------------------------------------
! Factorial, but giving the answer as a real.
! Stores calculated results to avoid excess computation.
! ----------------------------------------------------------------------
impure elemental function real_factorial(input) result(output)
  implicit none
  
  integer, intent(in) :: input
  real(dp)            :: output
  
  ! An array of factorials, stored to avoid re-calculating.
  ! Calculated the first time a related function is called.
  real(dp), allocatable, save :: real_factorials(:)
  
  ! Check that n>=0.
  if (input<0) then
    call print_line( ERROR//': Trying to calculate the factorial of an &
       &integer less than 0.')
    call err()
  endif
  
  call calculate_real_factorials(real_factorials,input)
  
  ! Check that n is small enough that n! is storable as an dp float.
  if (input>=size(real_factorials)) then
    call print_line( ERROR//': Trying to calculate the factorial of an &
       &integer greater than '//(size(real_factorials)-1)//'. This is too &
       &large to fit in the double-precision real type.')
    call err()
  endif
  
  ! real_factorials(1) = 0!, so real_factorials(n+1) = n!.
  output = real_factorials(input+1)
end function

subroutine calculate_real_factorials(real_factorials,input)
  implicit none
  
  real(dp), intent(inout), allocatable :: real_factorials(:)
  integer,  intent(in)                 :: input
  
  integer  :: i,ialloc
  real(dp) :: j
  
  ! The first time this function is called,
  !    calculate all storable factorials, {i!} s.t. i!<=huge(0).
  if (.not. allocated(real_factorials)) then
    ! Find the largest integer i s.t. i! is too large to be
    !    stored as a real.
    i = 0
    j = 1
    do
      i = i+1
      if (j>huge(0.0_dp)/i) then
        allocate(real_factorials(i), stat=ialloc); call err(ialloc)
        exit
      endif
      j = j*i ! j = i!.
    enddo
    
    ! Calculate and store all storable factorials.
    real_factorials(1) = 1.0_dp
    do i=1,size(real_factorials)-1
      real_factorials(i+1) = real_factorials(i)*i ! real_factorials(i+1) = i!.
    enddo
  endif
end subroutine

! ----------------------------------------------------------------------
! Calculates ln(n!), the log of the factorial of the input.
! ----------------------------------------------------------------------
impure elemental function log_factorial(input) result(output)
  implicit none
  
  integer, intent(in) :: input
  real(dp)            :: output
  
  ! An array of ln(n!), stored to avoid re-calculating.
  real(dp), allocatable, save :: log_factorials(:)
  integer,               save :: log_factorials_calculated = -1
  
  if (input<0) then
    call print_line( ERROR//': Trying to calculate the factorial of an &
       &integer less than 0.')
    call err()
  endif
  
  if (input > log_factorials_calculated) then
    call calculate_log_factorials( log_factorials,            &
                                 & log_factorials_calculated, &
                                 & input                      )
  endif
  
  ! log_factorials(1) = ln(0!), so ln(n!) = log_factorials(n+1).
  output = log_factorials(input+1)
end function

subroutine calculate_log_factorials(log_factorials,log_factorials_calculated, &
   & input)
  implicit none
  
  real(dp), intent(inout), allocatable :: log_factorials(:)
  integer,  intent(inout)              :: log_factorials_calculated
  integer,  intent(in)                 :: input
  
  real(dp), allocatable :: old_log_factorials(:)
  
  integer :: n,ialloc
  
  ! Check if log_factorials needs expanding.
  ! Uses a doubling strategy to minimise the number of re-allocations.
  if (.not. allocated(log_factorials)) then
    allocate(log_factorials(input+1), stat=ialloc); call err(ialloc)
    log_factorials(1) = 0
    log_factorials_calculated = 0
  elseif (size(log_factorials)<=input) then
    old_log_factorials = log_factorials
    deallocate(log_factorials, stat=ialloc); call err(ialloc)
    allocate( log_factorials(max(input+1,2*size(old_log_factorials))), &
            & stat=ialloc); call err(ialloc)
    log_factorials(:size(old_log_factorials)) = old_log_factorials
  endif
  
  ! Calculate additional ln(n!) as needed.
  ! Uses ln(n!) = ln((n-1)!) + ln(n).
  ! N.B. log_factorials(1) = ln(0!), so ln(n!) = log_factorials(n+1).
  do n=log_factorials_calculated+1,input
    log_factorials(n+1) = log_factorials(n) + log(real(n,dp))
  enddo
  log_factorials_calculated = max(log_factorials_calculated,input)
end subroutine

! ----------------------------------------------------------------------
! The product of the first n odd numbers, rather than the first n numbers.
! odd_factorial(n) = prod_{k=1}^n[ 2*k-1 ] = (2n-1)!/(2^n * (n-1)!).
! Stores calculated results to avoid excess computation.
! ----------------------------------------------------------------------
impure elemental function odd_factorial(input) result(output)
  implicit none
  
  integer, intent(in) :: input
  integer             :: output
  
  ! An array of odd_factorials, stored to avoid re-calculating.
  ! Calculated the first time a related function is called.
  integer, allocatable, save :: odd_factorials(:)
  
  ! Check that n>=0.
  if (input<0) then
    call print_line( ERROR//': Trying to calculate the factorial of an &
       &integer less than 0.')
    call err()
  endif
  
  call calculate_odd_factorials(odd_factorials)
  
  ! Check that n is small enough that n! is storable as an integer.
  if (input>=size(odd_factorials)) then
    call print_line( ERROR//': Trying to calculate the odd factorial of an &
       &integer greater than '//(size(odd_factorials)-1)//'. This is too &
       &large to fit in the standard integer type.')
    call err()
  endif
  
  ! odd_factorials(1) = odd_factorial(0),
  !    so odd_factorials(n+1) = odd_factorial(n).
  output = odd_factorials(input+1)
end function

subroutine calculate_odd_factorials(odd_factorials)
  implicit none
  
  integer, intent(inout), allocatable :: odd_factorials(:)
  
  integer :: i,j,ialloc
  
  ! The first time this function is called,
  !    calculate all storable factorials, {f(i)} s.t. f(i)<=huge(0).
  if (.not. allocated(odd_factorials)) then
    ! Find the largest integer i s.t. f(i) is too large to be
    !    stored as an integer.
    i = 0
    j = 1
    do
      i = i+1
      if (j>huge(0)/(2*i-1)) then
        allocate(odd_factorials(i), stat=ialloc); call err(ialloc)
        exit
      endif
      j = j*(2*i-1) ! j = f(i).
    enddo
    
    ! Calculate and store all storable odd factorials.
    odd_factorials(1) = 1
    do i=1,size(odd_factorials)-1
      ! oddfactorials(i+1) = f(i).
      odd_factorials(i+1) = odd_factorials(i)*(2*i-1)
    enddo
  endif
end subroutine

! ----------------------------------------------------------------------
! Odd factorial, but giving the answer as a real.
! Stores calculated results to avoid excess computation.
! ----------------------------------------------------------------------
impure elemental function real_odd_factorial(input) result(output)
  implicit none
  
  integer, intent(in) :: input
  real(dp)            :: output
  
  ! An array of odd factorials, stored to avoid re-calculating.
  ! Calculated the first time a related function is called.
  real(dp), allocatable, save :: real_odd_factorials(:)
  
  ! Check that n>=0.
  if (input<0) then
    call print_line( ERROR//': Trying to calculate the odd factorial of an &
       &integer less than 0.')
    call err()
  endif
  
  call calculate_real_odd_factorials(real_odd_factorials, input)
  
  ! Check that n is small enough that n! is storable as an dp float.
  if (input>=size(real_odd_factorials)) then
    call print_line( ERROR//': Trying to calculate the odd factorial of an &
       &integer greater than '//(size(real_odd_factorials)-1)//'. This is too &
       &large to fit in the double-precision real type.')
    call err()
  endif
  
  ! real_odd_factorials(1) = real_odd_factorial(0),
  !    so real_odd_factorials(n+1) = real_odd_factorial(n).
  output = real_odd_factorials(input+1)
end function

subroutine calculate_real_odd_factorials(real_odd_factorials,input)
  implicit none
  
  real(dp), intent(inout), allocatable :: real_odd_factorials(:)
  integer,  intent(in)                 :: input
  
  integer  :: i,ialloc
  real(dp) :: j
  
  ! The first time this function is called,
  !    calculate all storable odd factorials, {f(i)} s.t. f(i)<=huge(0).
  if (.not. allocated(real_odd_factorials)) then
    ! Find the largest integer i s.t. f(i) is too large to be
    !    stored as a real.
    i = 0
    j = 1
    do
      i = i+1
      if (j>huge(0.0_dp)/(2*i-1)) then
        allocate(real_odd_factorials(i), stat=ialloc); call err(ialloc)
        exit
      endif
      j = j*(2*i-1) ! j = f(i).
    enddo
    
    ! Calculate and store all storable odd factorials.
    real_odd_factorials(1) = 1.0_dp
    do i=1,size(real_odd_factorials)-1
      ! real_odd_factorials(i+1) = f(i).
      real_odd_factorials(i+1) = real_odd_factorials(i)*(2*i-1)
    enddo
  endif
end subroutine

! ----------------------------------------------------------------------
! Calculates ln(real_odd_factorial(n)),
!    the log of the odd factorial of the input.
! ----------------------------------------------------------------------
impure elemental function log_odd_factorial(input) result(output)
  implicit none
  
  integer, intent(in) :: input
  real(dp)            :: output
  
  ! An array of ln(f(n)), stored to avoid re-calculating.
  real(dp), allocatable, save :: log_odd_factorials(:)
  integer,               save :: log_odd_factorials_calculated
  
  if (input<0) then
    call print_line( ERROR//': Trying to calculate the odd factorial of an &
       &integer less than 0.')
    call err()
  endif
  
  call calculate_log_odd_factorials( log_odd_factorials,            &
                                   & log_odd_factorials_calculated, &
                                   & input                          )
  
  ! log_odd_factorials(1) = ln(f(0)), so ln(f(n)) = log_odd_factorials(n+1).
  output = log_odd_factorials(input+1)
end function

subroutine calculate_log_odd_factorials(log_odd_factorials, &
   & log_odd_factorials_calculated,input)
  implicit none
  
  real(dp), intent(inout), allocatable :: log_odd_factorials(:)
  integer,  intent(inout)              :: log_odd_factorials_calculated
  integer,  intent(in)                 :: input
  
  real(dp), allocatable :: old_log_odd_factorials(:)
  
  integer :: n,ialloc
  
  ! Check if log_odd_factorials needs expanding.
  ! Uses a doubling strategy to minimise the number of re-allocations.
  if (.not. allocated(log_odd_factorials)) then
    allocate(log_odd_factorials(input+1), stat=ialloc); call err(ialloc)
    log_odd_factorials(1) = 0
    log_odd_factorials_calculated = 0
  elseif (size(log_odd_factorials)<=input) then
    old_log_odd_factorials = log_odd_factorials
    deallocate(log_odd_factorials, stat=ialloc); call err(ialloc)
    allocate( log_odd_factorials(max( input+1,                           &
            &                         2*size(old_log_odd_factorials) )), &
            & stat=ialloc); call err(ialloc)
    log_odd_factorials(:size(old_log_odd_factorials)) = old_log_odd_factorials
  endif
  
  ! Calculate additional ln(f(n)) as needed.
  ! Uses ln(f(n)) = ln(f(n-1)) + ln(2*n-1).
  ! N.B. log_odd_factorials(1) = ln(f(0)),
  !    so ln(f(n)) = log_odd_factorials(n+1).
  do n=log_odd_factorials_calculated+1,input
    log_odd_factorials(n+1) = log_odd_factorials(n) + log(real((2*n-1),dp))
  enddo
  log_odd_factorials_calculated = max(log_odd_factorials_calculated,input)
end subroutine

! ----------------------------------------------------------------------
! Calculates the binomial coefficient of two integers.
! binomial(a,b) = a!/(b!*(a-b)!)
! ----------------------------------------------------------------------
impure elemental function binomial(top,bottom) result(output)
  implicit none
  
  integer, intent(in) :: top
  integer, intent(in) :: bottom
  integer             :: output
  
  integer :: smaller_bottom
  
  ! An array of binomials, stored to avoid re-calculating.
  ! Calculated the first time a related function is called.
  integer, allocatable, save :: binomials(:,:)
  integer, allocatable, save :: binomials_calculated(:)
  
  if (top<0) then
    call err()
  elseif (bottom<0) then
    call err()
  elseif (bottom>top) then
    call err()
  endif
  
  ! binomial(n,k)=binomial(n,n-k). Pick whichever of k and n-k is smallest.
  smaller_bottom = min(bottom,top-bottom)
  
  if (smaller_bottom==0) then
    output = 1
  else
    call calculate_binomials(binomials,binomials_calculated,top,smaller_bottom)
    
    ! N.B. unlike factorials, binomial(n,k)=binomials(n,k),
    !    since binomial(n,0)=1, and does not need to be stored.
    output = binomials(top,smaller_bottom)
  endif
end function

subroutine calculate_binomials(binomials,binomials_calculated,top,bottom)
  implicit none
  
  integer, intent(inout), allocatable :: binomials(:,:)
  integer, intent(inout), allocatable :: binomials_calculated(:)
  integer, intent(in)                 :: top
  integer, intent(in)                 :: bottom
  
  integer, allocatable :: old_binomials(:,:)
  
  integer :: old_n,old_k
  integer :: max_n,max_k
  
  integer :: n,k,ialloc
  
  ! If the binomials array is too small, reallocate it.
  ! Uses a doubling strategy to reduce the number of reallocations.
  if (.not. allocated(binomials)) then
    ! Allocate binomials and binomials_calculated, and initalise:
    ! binomials = (1         )  binomials_calculated=(1, 2, ..., k)
    !             (2, 1      )
    !             (3, ?, 1   )
    !             (.  .  . . )
    !             (k ...    1)
    !             (.  .  .  .)
    !             (n, ?, ?, ?)
    allocate( binomials(top,bottom),        &
            & binomials_calculated(bottom), &
            & stat=ialloc); call err(ialloc)
    binomials(:,1) = [(n,n=1,top)]
    binomials_calculated(1) = top
    do k=2,bottom
      binomials(k,k) = 1
      binomials_calculated(k) = k
    enddo
  elseif (top>size(binomials,1) .or. bottom>size(binomials,2)) then
    ! Reallocate binomials, and copy over the old contents.
    old_n = size(binomials,1)
    old_k = size(binomials,2)
    max_n = max(top,2*old_n)
    max_k = max(bottom,2*old_k)
    old_binomials = binomials
    deallocate(binomials, stat=ialloc); call err(ialloc)
    allocate(binomials(max_n,max_k), stat=ialloc); call err(ialloc)
    binomials(:old_n,:old_k) = old_binomials
    
    ! Set the (n,1) column to n.
    binomials(old_n+1:,1) = [(n,n=old_n+1,max_n)]
    
    ! Set the (k,k) diagonal to 1.
    do k=old_k+1,max_k
      binomials(k,k) = 1
    enddo
    
    ! Update binomials_calculated.
    binomials_calculated(1) = max_n
    binomials_calculated = [binomials_calculated, [(k,k=old_k+1,max_k)]]
  endif
  
  ! If the requested binomial has not been calculated, calculate it.
  ! Calculates binomials using the recurrence relation
  !    b(n,k) = b(n-1,k) + b(n-1,k-1)
  !    b(n,1) = n
  !    b(n,n) = 1
  if (binomials_calculated(bottom)<top) then
    ! Calculate b(n,0) up to b(top-bottom,0) and b(n,n) up to b(bottom,bottom).
    do k=2,bottom
      do n=binomials_calculated(k)+1,top-bottom+k
        if (huge(0)-binomials(n-1,k)<binomials(n-1,k-1)) then
          call print_line(ERROR//': Binomial calculation overflowing. &
             &Binomial('//top//','//bottom//') is too large.')
          call err()
        endif
        binomials(n,k) = binomials(n-1,k)+binomials(n-1,k-1)
      enddo
      binomials_calculated(k) = top-bottom+k
    enddo
  endif
end subroutine

impure elemental function real_binomial(top,bottom) result(output)
  implicit none
  
  integer, intent(in) :: top
  integer, intent(in) :: bottom
  real(dp)            :: output
  
  integer :: smaller_bottom
  
  ! An array of binomials, stored to avoid re-calculating.
  ! Calculated the first time a related function is called.
  real(dp), allocatable, save :: real_binomials(:,:)
  integer,  allocatable, save :: real_binomials_calculated(:)
  
  if (top<0) then
    call err()
  elseif (bottom<0) then
    call err()
  elseif (bottom>top) then
    call err()
  endif
  
  ! binomial(n,k)=binomial(n,n-k). Pick whichever of k and n-k is smallest.
  smaller_bottom = min(bottom,top-bottom)
  
  if (smaller_bottom==0) then
    output = 1
  else
    call calculate_real_binomials( real_binomials,            &
                                 & real_binomials_calculated, &
                                 & top,                       &
                                 & smaller_bottom             )
    
    ! N.B. unlike factorials, binomial(n,k)=BINOMIALS(n,k),
    !    since binomial(n,0)=1, and does not need to be stored.
    output = real_binomials(top,smaller_bottom)
  endif
end function

subroutine calculate_real_binomials(real_binomials,real_binomials_calculated, &
   & top,bottom)
  implicit none
  
  real(dp), intent(inout), allocatable :: real_binomials(:,:)
  integer,  intent(inout), allocatable :: real_binomials_calculated(:)
  integer,  intent(in)                 :: top
  integer,  intent(in)                 :: bottom
  
  real(dp), allocatable :: old_binomials(:,:)
  
  integer :: old_n,old_k
  integer :: max_n,max_k
  
  integer :: n,k,ialloc
  
  ! If the real_binomials array is too small, reallocate it.
  ! Uses a doubling strategy to reduce the number of reallocations.
  if (.not. allocated(real_binomials)) then
    ! Allocate real_binomials and real_binomials_calculated, and initalise:
    ! real_binomials = (1         )  real_binomials_calculated=(1, 2, ..., k)
    !                  (2, 1      )
    !                  (3, ?, 1   )
    !                  (.  .  . . )
    !                  (k ...    1)
    !                  (.  .  .  .)
    !                  (n, ?, ?, ?)
    allocate( real_binomials(top,bottom),        &
            & real_binomials_calculated(bottom), &
            & stat=ialloc); call err(ialloc)
    real_binomials(:,1) = [(n,n=1,top)]
    real_binomials_calculated(1) = top
    do k=2,bottom
      real_binomials(k,k) = 1
      real_binomials_calculated(k) = k
    enddo
  elseif (top>size(real_binomials,1) .or. bottom>size(real_binomials,2)) then
    ! Reallocate real_binomials, and copy over the old contents.
    old_n = size(real_binomials,1)
    old_k = size(real_binomials,2)
    max_n = max(top,2*old_n)
    max_k = max(bottom,2*old_k)
    old_binomials = real_binomials
    deallocate(real_binomials, stat=ialloc); call err(ialloc)
    allocate(real_binomials(max_n,max_k), stat=ialloc); call err(ialloc)
    real_binomials(:old_n,:old_k) = old_binomials
    
    ! Set the (n,1) column to n.
    real_binomials(old_n+1:,1) = [(n,n=old_n+1,max_n)]
    
    ! Set the (k,k) diagonal to 1.
    do k=old_k+1,max_k
      real_binomials(k,k) = 1
    enddo
    
    ! Update real_binomials_calculated.
    real_binomials_calculated(1) = max_n
    real_binomials_calculated = [ real_binomials_calculated, &
                                & [(k,k=old_k+1,max_k)]      ]
  endif
  
  ! If the requested binomial has not been calculated, calculate it.
  ! Calculates binomials using the recurrence relation
  !    b(n,k) = b(n-1,k) + b(n-1,k-1)
  !    b(n,1) = n
  !    b(n,n) = 1
  if (real_binomials_calculated(bottom)<top) then
    ! Calculate b(n,0) up to b(top-bottom,0) and b(n,n) up to b(bottom,bottom).
    do k=2,bottom
      do n=real_binomials_calculated(k)+1,top-bottom+k
        if (huge(0)-real_binomials(n-1,k)<real_binomials(n-1,k-1)) then
          call print_line(ERROR//': Binomial calculation overflowing. &
             &Binomial('//top//','//bottom//') is too large.')
          call err()
        endif
        real_binomials(n,k) = real_binomials(n-1,k)+real_binomials(n-1,k-1)
      enddo
      real_binomials_calculated(k) = top-bottom+k
    enddo
  endif
end subroutine

! ----------------------------------------------------------------------
! Calculates ln(bin(top,bottom)) = ln(top!/(bottom!(top-bottom)!)),
!    the log of the binomial of two integers.
! ----------------------------------------------------------------------
impure elemental function log_binomial(top,bottom) result(output)
  implicit none
  
  integer, intent(in) :: top
  integer, intent(in) :: bottom
  real(dp)            :: output
  
  output = log_factorial(top)    &
       & - log_factorial(bottom) &
       & - log_factorial(top-bottom)
end function

! ----------------------------------------------------------------------
! Calculates the multinomial coefficient of a numerator and a set of
!    denominators.
! multinomial(a,[b,c,d,...]) = a!/(b!c!d!...)
! ----------------------------------------------------------------------
function multinomial(top,bottom) result(output)
  implicit none
  
  integer, intent(in) :: top
  integer, intent(in) :: bottom(:)
  integer             :: output
  
  if (top<0) then
    call err()
  elseif (any(bottom<0)) then
    call err()
  elseif (sum(bottom)>top) then
    call err()
  endif
  
  output = factorial(top) / product(factorial(bottom))
end function

function real_multinomial(top,bottom) result(output)
  implicit none
  
  integer, intent(in) :: top
  integer, intent(in) :: bottom(:)
  real(dp)            :: output
  
  if (top<0) then
    call err()
  elseif (any(bottom<0)) then
    call err()
  elseif (sum(bottom)>top) then
    call err()
  endif
  
  output = exp(log_factorial(top)-sum(log_factorial(bottom)))
end function

! ----------------------------------------------------------------------
! Calculates ln(multinomial(top,bottom)) = ln(top!/product(bottom!)),
!    the log of the multinomial of a numerator and set of denominators.
! ----------------------------------------------------------------------
function log_multinomial(top,bottom) result(output)
  implicit none
  
  integer, intent(in) :: top
  integer, intent(in) :: bottom(:)
  real(dp)            :: output
  
  output = log_factorial(top) - sum(log_factorial(bottom))
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
  
  integer :: a
  integer :: b
  integer :: a_old
  
  a = min(abs(int_1), abs(int_2))
  b = max(abs(int_1), abs(int_2))
  
  ! Loop until the g.c.d. is found, obeying the following:
  !    gcd(0,b) = b
  !    gcd(a,b) = gcd(modulo(b,a), a)
  do while (a/=0)
    a_old = a
    a = modulo(b,a)
    b = a_old
  enddo
  
  output = b
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
! exp(2*pi*i*input), cos(2*pi*input) and sin(2*pi*input).
! ----------------------------------------------------------------------
impure elemental function exp_2pii_real(input) result(output)
  implicit none
  
  real(dp), intent(in) :: input
  complex(dp)          :: output
  
  real(dp) :: exponent
  
  exponent = 2*PI*input
  output = cmplx(cos(exponent),sin(exponent),dp)
end function

impure elemental function cos_2pi_real(input) result(output)
  implicit none
  
  real(dp), intent(in) :: input
  real(dp)             :: output
  
  output = cos(2*PI*input)
end function

impure elemental function sin_2pi_real(input) result(output)
  implicit none
  
  real(dp), intent(in) :: input
  real(dp)             :: output
  
  output = sin(2*PI*input)
end function
end module
