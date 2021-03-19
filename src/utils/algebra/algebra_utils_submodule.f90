submodule (caesar_algebra_utils_module) caesar_algebra_utils_submodule
  use caesar_algebra_module
contains

module procedure l2_norm_reals
  output = sqrt(dot_product(input,input))
end procedure

module procedure l2_norm_complexes
  output = sqrt(abs(dot_product(input,input)))
end procedure

module procedure sum_squares_RealVector
  output = sum(dble(input)*dble(input))
end procedure

module procedure sum_squares_RealMatrix
  output = sum(dble(input)*dble(input))
end procedure

module procedure sum_squares_ComplexVector
  output = sum(real(cmplx(input)*conjg(cmplx(input))))
end procedure

module procedure sum_squares_ComplexMatrix
  output = sum(real(cmplx(input)*conjg(cmplx(input))))
end procedure

module procedure triple_product_IntVector
  integer :: matrix(3,3)
  
  matrix(:,1) = int(a)
  matrix(:,2) = int(b)
  matrix(:,3) = int(c)
  
  output = determinant(matrix)
end procedure

module procedure triple_product_RealVector
  real(dp) :: matrix(3,3)
  
  matrix(:,1) = dble(a)
  matrix(:,2) = dble(b)
  matrix(:,3) = dble(c)
  
  output = determinant(matrix)
end procedure

module procedure factorial
  ! An array of factorials, stored to avoid re-calculating.
  ! Calculated the first time a related module procedure is called.
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
end procedure

subroutine calculate_factorials(factorials)
  integer, intent(inout), allocatable :: factorials(:)
  
  integer :: i,j,ialloc
  
  ! The first time this module procedure is called,
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

module procedure real_factorial
  ! An array of factorials, stored to avoid re-calculating.
  ! Calculated the first time a related module procedure is called.
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
end procedure

subroutine calculate_real_factorials(real_factorials,input)
  real(dp), intent(inout), allocatable :: real_factorials(:)
  integer,  intent(in)                 :: input
  
  integer  :: i,ialloc
  real(dp) :: j
  
  ! The first time this module procedure is called,
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

module procedure log_factorial
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
end procedure

subroutine calculate_log_factorials(log_factorials,log_factorials_calculated, &
   & input)
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

module procedure odd_factorial
  ! An array of odd_factorials, stored to avoid re-calculating.
  ! Calculated the first time a related module procedure is called.
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
end procedure

subroutine calculate_odd_factorials(odd_factorials)
  integer, intent(inout), allocatable :: odd_factorials(:)
  
  integer :: i,j,ialloc
  
  ! The first time this module procedure is called,
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

module procedure real_odd_factorial
  ! An array of odd factorials, stored to avoid re-calculating.
  ! Calculated the first time a related module procedure is called.
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
end procedure

subroutine calculate_real_odd_factorials(real_odd_factorials,input)
  real(dp), intent(inout), allocatable :: real_odd_factorials(:)
  integer,  intent(in)                 :: input
  
  integer  :: i,ialloc
  real(dp) :: j
  
  ! The first time this module procedure is called,
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

module procedure log_odd_factorial
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
end procedure

subroutine calculate_log_odd_factorials(log_odd_factorials, &
   & log_odd_factorials_calculated,input)
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

module procedure binomial
  integer :: smaller_bottom
  
  ! An array of binomials, stored to avoid re-calculating.
  ! Calculated the first time a related module procedure is called.
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
end procedure

subroutine calculate_binomials(binomials,binomials_calculated,top,bottom)
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

module procedure real_binomial
  integer :: smaller_bottom
  
  ! An array of binomials, stored to avoid re-calculating.
  ! Calculated the first time a related module procedure is called.
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
end procedure

module subroutine calculate_real_binomials(real_binomials, &
   & real_binomials_calculated,top,bottom)
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

module procedure log_binomial
  output = log_factorial(top)    &
       & - log_factorial(bottom) &
       & - log_factorial(top-bottom)
end procedure

module procedure multinomial
  if (top<0) then
    call err()
  elseif (any(bottom<0)) then
    call err()
  elseif (sum(bottom)>top) then
    call err()
  endif
  
  output = factorial(top) / product(factorial(bottom))
end procedure

module procedure real_multinomial
  if (top<0) then
    call err()
  elseif (any(bottom<0)) then
    call err()
  elseif (sum(bottom)>top) then
    call err()
  endif
  
  output = exp(log_factorial(top)-sum(log_factorial(bottom)))
end procedure

module procedure log_multinomial
  output = log_factorial(top) - sum(log_factorial(bottom))
end procedure

module procedure int_sqrt
  output = nint(sqrt(real(input,dp)))
  if (output*output/=input) then
    call err()
  endif
end procedure

module procedure gcd_2
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
end procedure

module procedure gcd_3
  output = gcd(gcd(int_1,int_2),int_3)
end procedure

module procedure gcd_4
  output = gcd(gcd(int_1,int_2,int_3),int_4)
end procedure

module procedure gcd_integers
  if (size(input)==0) then
    output = 0
  elseif (size(input)==1) then
    output = abs(input(1))
  else
    output = gcd([gcd(input(1),input(2)), input(3:)])
  endif
end procedure

module procedure lcm_2
  integer :: product
  
  product = abs(int_1*int_2)
  
  if (product>0) then
    output = product/gcd(int_1,int_2)
  else
    output = 0
  endif
end procedure

module procedure lcm_3
  output = lcm(lcm(int_1,int_2),int_3)
end procedure

module procedure lcm_4
  output = lcm(lcm(int_1,int_2,int_3),int_4)
end procedure

module procedure lcm_integers
  if (size(input)==0) then
    output = 1
  elseif (size(input)==1) then
    output = abs(input(1))
  else
    output = lcm([lcm(input(1),input(2)), input(3:)])
  endif
end procedure

module procedure exp_2pii_real
  real(dp) :: exponent
  
  exponent = 2*PI*input
  output = cmplx(cos(exponent),sin(exponent),dp)
end procedure

module procedure cos_2pi_real
  output = cos(2*PI*input)
end procedure

module procedure sin_2pi_real
  output = sin(2*PI*input)
end procedure

module procedure bound
  if (lower_bound>upper_bound) then
    call print_line(ERROR//': Upper bound < lower bound.')
    call err()
  endif
  
  output = max(lower_bound,min(input,upper_bound))
end procedure
end submodule
