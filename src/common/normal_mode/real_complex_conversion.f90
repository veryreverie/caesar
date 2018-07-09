! ======================================================================
! Converts objects between complex and real co-ordinates.
! ======================================================================

! Complex modes with ids i and j, with i<j, are paired under the transformation
!    q -> -q if they transform as u_i -> u_j and u_j -> u_i.
! N.B. u_i = u_j*.
! These modes can then be transformed to and from the real modes i and j
!    as follows.
! The complex u_i and u_j will be written u+ and u-,
!    and the real u_i and u_j will be written c and s respectively.
!
! N.B. Taking q to be the q-point of u+, u+ looks like e^{+2pi*i*q.r},
!                                        u- looks like e^{-2pi*i*q.r},
!                                        c  looks like cos(2pi*q.r),
!                                        s  looks like sin(2pi*q.r).
! 
! N.B. the choice of r_i as the cos mode and r_j as the sin mode is arbitrary,
!    but will be used consistently throughout.
!
! u+ = (c + i*s) /  sqrt(2)
! u- = (c - i*s) /  sqrt(2)
! c  = (u+ + u-) /  sqrt(2)    = real(u+)*sqrt(2) =  real(u-)*sqrt(2)
! s  = (u+ - u-) / (sqrt(2)*i) = imag(u+)*sqrt(2) = -imag(u-)*sqrt(2)
!
! u+ = (a+ib)e^{2pi i q.R}
! u- = (a-ib)e^{2pi i q.R}
! c  = sqrt(2)(a*cos(2pi q.R) - b*sin(2pi q.R))
! s  = sqrt(2)(a*sin(2pi q.R) + b*cos(2pi q.R))
! 
! Under q -> -q : u+ ->   u-
!                 u- ->   u+
!                 c  ->   c
!                 s  -> - s
!
! N.B. if i=j, then the mode is real (2q=G, so e^{+iq.r}=e^{-iq.r}=cos(q.r)).
! In this case, u+ = c, and u- = s = 0.

module real_complex_conversion_submodule
  use utils_module
  
  use complex_mode_submodule
  use real_mode_submodule
  use complex_polynomial_submodule
  use real_polynomial_submodule
  implicit none
  
  private
  
  public :: complex_to_real
  public :: real_to_complex
  public :: conversion_matrix
  
  interface complex_to_real
    module procedure complex_to_real_Modes
  end interface
  
  interface real_to_complex
    module procedure real_to_complex_Modes
  end interface
  
  interface conversion_matrix
    module procedure conversion_matrix_ComplexMonomials_RealMonomials
    module procedure conversion_matrix_RealMonomials_ComplexMonomials
  end interface
  
  interface element
    module procedure element_ComplexMonomial_RealMonomial
    module procedure element_RealMonomial_ComplexMonomial
    
    module procedure element_ComplexUnivariates_RealUnivariates
    module procedure element_RealUnivariates_ComplexUnivariates
  end interface
contains

! ----------------------------------------------------------------------
! Functions for converting between real and complex co-ordinates.
! ----------------------------------------------------------------------
! Construct a set of real modes from a set of complex modes, and vice versa.
function complex_to_real_Modes(input) result(output)
  implicit none
  
  type(ComplexMode), intent(in) :: input(:)
  type(RealMode), allocatable   :: output(:)
  
  type(RealVector), allocatable :: a(:)
  type(RealVector), allocatable :: b(:)
  
  type(RealVector), allocatable :: cos_vector(:)
  type(RealVector), allocatable :: sin_vector(:)
  
  integer :: qpoint_id_plus
  integer :: qpoint_id_minus
  
  integer :: mode,ialloc
  
  allocate(output(size(input)), stat=ialloc); call err(ialloc)
  do mode=1,size(input)
    if (input(mode)%id==input(mode)%paired_id) then
      ! output(mode) is its own pair. It is already real.
      ! c = u+, s = u- = 0.
      qpoint_id_plus = input(mode)%qpoint_id
      qpoint_id_minus = qpoint_id_plus
      a = real(input(mode)%unit_vector)
      b = a * 0.0_dp
      cos_vector = a
      sin_vector = b
    elseif (input(mode)%id<input(mode)%paired_id) then
      ! output(mode) is c.
      ! input(mode) is u+.
      ! u+ = (a+ib)e^{2pi i q.R}
      ! c  = sqrt(2)(a*cos(2pi q.R) - b*sin(2pi q.R))
      qpoint_id_plus = input(mode)%qpoint_id
      qpoint_id_minus = input(mode)%paired_qpoint_id
      a = real( input(mode)%unit_vector)
      b = aimag(input(mode)%unit_vector)
      cos_vector =  sqrt(2.0_dp) * a
      sin_vector = -sqrt(2.0_dp) * b
    elseif (input(mode)%id>input(mode)%paired_id) then
      ! output(mode) is s.
      ! input(mode) is u-.
      ! u- = (a-ib)e^{2pi i q.R}
      ! s  = sqrt(2)(a*sin(2pi q.R) + b*cos(2pi q.R))
      a =  real( input(mode)%unit_vector)
      b = -aimag(input(mode)%unit_vector)
      cos_vector = sqrt(2.0_dp) * b
      sin_vector = sqrt(2.0_dp) * a
    else
      call err()
    endif
    
    output(mode) = RealMode(                                    &
       & id                   = input(mode)%id,                 &
       & paired_id            = input(mode)%paired_id,          &
       & frequency            = input(mode)%frequency,          &
       & spring_constant      = input(mode)%spring_constant,    &
       & soft_mode            = input(mode)%soft_mode,          &
       & translational_mode   = input(mode)%translational_mode, &
       & cos_vector           = cos_vector,                     &
       & sin_vector           = sin_vector,                     &
       & qpoint_id_plus       = qpoint_id_plus,                 &
       & qpoint_id_minus      = qpoint_id_minus,                &
       & subspace_id          = input(mode)%subspace_id)
  enddo
end function

function real_to_complex_Modes(input) result(output)
  implicit none
  
  type(RealMode), intent(in)     :: input(:)
  type(ComplexMode), allocatable :: output(:)
  
  type(RealVector), allocatable :: a(:)
  type(RealVector), allocatable :: b(:)
  
  type(ComplexVector), allocatable :: unit_vector(:)
  
  integer :: qpoint_id
  integer :: paired_qpoint_id
  
  integer :: i,ialloc
  
  allocate(output(size(input)), stat=ialloc); call err(ialloc)
  do i=1,size(input)
    ! Construct vectors and q-point ids.
    if (input(i)%id==input(i)%paired_id) then
      ! output(i) is its own pair. It is already real.
      ! c = u+, s = u- = 0.
      qpoint_id = input(i)%qpoint_id_plus
      paired_qpoint_id = qpoint_id
      a = input(i)%cos_vector
      unit_vector = cmplxvec(a)
    elseif (input(i)%id<input(i)%paired_id) then
      ! output(i) is u+.
      ! input(i) is c.
      ! c  = sqrt(2)(a*cos(2pi q.R) - b*sin(2pi q.R))
      ! u+ = (a+ib)e^{2pi i q.R}
      qpoint_id = input(i)%qpoint_id_plus
      paired_qpoint_id = input(i)%qpoint_id_minus
      a =  input(i)%cos_vector / sqrt(2.0_dp)
      b = -input(i)%sin_vector / sqrt(2.0_dp)
      unit_vector = cmplxvec(a,b)
    elseif (input(i)%paired_id<input(i)%id) then
      ! output(i) is u-.
      ! input(i) is s.
      ! s  = sqrt(2)(a*sin(2pi q.R) + b*cos(2pi q.R))
      ! u- = (a-ib)e^{2pi i q.R}
      qpoint_id = input(i)%qpoint_id_minus
      paired_qpoint_id = input(i)%qpoint_id_plus
      a = input(i)%sin_vector / sqrt(2.0_dp)
      b = input(i)%cos_vector / sqrt(2.0_dp)
      unit_vector = cmplxvec(a,-b)
    else
      call err()
    endif
    
    output(i) = ComplexMode(                                 &
       & id                   = input(i)%id,                 &
       & paired_id            = input(i)%paired_id,          &
       & frequency            = input(i)%frequency,          &
       & spring_constant      = input(i)%spring_constant,    &
       & soft_mode            = input(i)%soft_mode,          &
       & translational_mode   = input(i)%translational_mode, &
       & unit_vector          = unit_vector,                 &
       & qpoint_id            = qpoint_id,                   &
       & paired_qpoint_id     = paired_qpoint_id,            &
       & subspace_id          = input(i)%subspace_id)
  enddo
end function

! Construct the matrix M which converts from one co-ordinate set to the other,
!    s.t. this = M * that.
function conversion_matrix_ComplexMonomials_RealMonomials(this,that) &
   & result(output)
  implicit none
  
  type(ComplexMonomial), intent(in) :: this(:)
  type(RealMonomial),    intent(in) :: that(:)
  type(ComplexMatrix)               :: output
  
  complex(dp), allocatable :: matrix(:,:)
  
  integer :: i,j,ialloc
  
  allocate(matrix(size(this),size(that)), stat=ialloc); call err(ialloc)
  do i=1,size(that)
    do j=1,size(this)
      matrix(j,i) = element(this(j),that(i))
    enddo
  enddo
  
  output = matrix
end function

function conversion_matrix_RealMonomials_ComplexMonomials(this,that) &
   & result(output)
  implicit none
  
  type(RealMonomial),    intent(in) :: this(:)
  type(ComplexMonomial), intent(in) :: that(:)
  type(ComplexMatrix)               :: output
  
  complex(dp), allocatable :: matrix(:,:)
  
  integer :: i,j,ialloc
  
  allocate(matrix(size(this),size(that)), stat=ialloc); call err(ialloc)
  do i=1,size(that)
    do j=1,size(this)
      matrix(j,i) = element(this(j),that(i))
    enddo
  enddo
  
  output = matrix
end function

! ----------------------------------------------------------------------
! Private helper functions.
! ----------------------------------------------------------------------

! Calculate an element for conversion_matrix.
!
! Calculates the coefficient of 'that' in the expansion of 'this'.
! i.e. 'this' = sum[ element(this,that) * 'that' ].
impure elemental function element_ComplexMonomial_RealMonomial(this,that) &
   & result(output)
  implicit none
  
  type(ComplexMonomial), intent(in) :: this
  type(RealMonomial),    intent(in) :: that
  complex(dp)                       :: output
  
  integer :: i,j,k,l
  
  if (sum(this%modes%power)/=sum(that%modes%power)) then
    output = 0
    return
  endif
  
  output = this%coefficient / that%coefficient
  ! Loop over modes, multiplying the output by the overlap of each mode pair
  !    in 'this' with the eqivalent mode pair in 'that'.
  do i=1,size(this)
    ! Find the location of the mode paired to this one in 'this'.
    ! Also find the location of this mode and its pair in 'that'.
    j = first(this%modes%paired_id() == this%modes(i)%id, default=0)
    k = first(that%modes%id          == this%modes(i)%id, default=0)
    l = first(that%modes%paired_id() == this%modes(i)%id, default=0)
    
    if (j/=0 .and. j<i) then
      ! This pairing has already been accounted for, when the paired mode
      !    came up.
      cycle
    endif
    
    if (k==0 .and. l==0) then
      ! Neither this mode nor its pair is present in 'that'.
      output = 0
      return
    elseif (j==0 .and. k==0) then
      output = output * element( this=this%modes(i), &
                               & that=that%modes(l))
    elseif (j==0 .and. l==0) then
      output = output * element( this=this%modes(i), &
                               & that=that%modes(k))
    elseif (j==0) then
      output = output * element( this     =this%modes(i), &
                               & that     =that%modes(k), &
                               & that_pair=that%modes(l))
    elseif (k==0) then
      output = output * element( this     =this%modes(i), &
                               & this_pair=this%modes(j), &
                               & that     =that%modes(l))
    elseif (l==0) then
      output = output * element( this     =this%modes(i), &
                               & this_pair=this%modes(j), &
                               & that     =that%modes(k))
    else
      output = output * element( this     =this%modes(i), &
                               & this_pair=this%modes(j), &
                               & that     =that%modes(k), &
                               & that_pair=that%modes(l))
    endif
  enddo
end function

! Calculate an element for conversion_matrix.
!
! Calculates the coefficient of 'that' in the expansion of 'this'.
! i.e. 'this' = sum[ element(this,that) * 'that' ].
impure elemental function element_RealMonomial_ComplexMonomial(this,that) &
   & result(output)
  implicit none
  
  type(RealMonomial),    intent(in) :: this
  type(ComplexMonomial), intent(in) :: that
  complex(dp)                       :: output
  
  integer :: i,j,k,l
  
  if (sum(this%modes%power)/=sum(that%modes%power)) then
    output = 0
    return
  endif
  
  output = this%coefficient / that%coefficient
  ! Loop over modes, multiplying the output by the overlap of each mode pair
  !    in 'this' with the eqivalent mode pair in 'that'.
  do i=1,size(this)
    ! Find the location of the mode paired to this one in 'this'.
    ! Also find the location of this mode and its pair in 'that'.
    j = first(this%modes%paired_id() == this%modes(i)%id, default=0)
    k = first(that%modes%id          == this%modes(i)%id, default=0)
    l = first(that%modes%paired_id() == this%modes(i)%id, default=0)
    
    if (j/=0 .and. j<i) then
      ! This pairing has already been accounted for, when the paired mode
      !    came up.
      cycle
    endif
    
    if (k==0 .and. l==0) then
      ! Neither this mode nor its pair is present in 'that'.
      output = 0
      return
    elseif (j==0 .and. k==0) then
      output = output * element( this=this%modes(i), &
                               & that=that%modes(l))
    elseif (j==0 .and. l==0) then
      output = output * element( this=this%modes(i), &
                               & that=that%modes(k))
    elseif (j==0) then
      output = output * element( this     =this%modes(i), &
                               & that     =that%modes(k), &
                               & that_pair=that%modes(l))
    elseif (k==0) then
      output = output * element( this     =this%modes(i), &
                               & this_pair=this%modes(j), &
                               & that     =that%modes(l))
    elseif (l==0) then
      output = output * element( this     =this%modes(i), &
                               & this_pair=this%modes(j), &
                               & that     =that%modes(k))
    else
      output = output * element( this     =this%modes(i), &
                               & this_pair=this%modes(j), &
                               & that     =that%modes(k), &
                               & that_pair=that%modes(l))
    endif
  enddo
end function

impure elemental function element_ComplexUnivariates_RealUnivariates &
   & (this,this_pair,that,that_pair) result(output)
  implicit none
  
  type(ComplexUnivariate), intent(in)           :: this
  type(ComplexUnivariate), intent(in), optional :: this_pair
  type(RealUnivariate),    intent(in)           :: that
  type(RealUnivariate),    intent(in), optional :: that_pair
  complex(dp)                                   :: output
  
  type(ComplexUnivariate) :: up ! u+ = (c  + i*s)/sqrt(2)
  type(ComplexUnivariate) :: um ! u- = (c  - i*s)/sqrt(2)
  type(RealUnivariate)    :: c  ! c  = (u+ + u- )/sqrt(2)
  type(RealUnivariate)    :: s  ! s  = (u+ - u- )/sqrt(2)i
  
  integer :: j,k,l,n
  
  ! Check for optional arguments, and identify which mode is i and which is j.
  if (present(this_pair)) then
    if (this%id<this%paired_id()) then
      up = this
      um = this_pair
    else
      up = this_pair
      um = this
    endif
  else
    if (this%id<this%paired_id()) then
      up = this
      um = ComplexUnivariate(this%paired_id(),this%id,power=0)
    else
      up = ComplexUnivariate(this%paired_id(),this%id,power=0)
      um = this
    endif
  endif
  
  if (present(that_pair)) then
    if (that%id<that%paired_id()) then
      c = that
      s = that_pair
    else
      c = that_pair
      s = that
    endif
  else
    if (that%id<that%paired_id()) then
      c = that
      s = RealUnivariate(that%paired_id(),that%id,power=0)
    else
      c = RealUnivariate(that%paired_id(),that%id,power=0)
      s = that
    endif
  endif
  
  ! Check that inputs are paired up.
  if (up%paired_id()/=um%id .or. up%id/=um%paired_id()) then
    call print_line(CODE_ERROR//': Unpaired modes passed to function.')
    call err()
  elseif (c%paired_id()/=s%id .or. c%id/=s%paired_id()) then
    call print_line(CODE_ERROR//': Unpaired modes passed to function.')
    call err()
  endif
  
  ! If the modes don't match, or the powers don't match, then the element
  !    is zero.
  if (up%id/=c%id .or. um%id/=s%id) then
    output = 0
    return
  elseif (up%power+um%power /= c%power+s%power) then
    output = 0
    return
  endif
  
  ! Check if the mode is real. If so, the overlap is 1.
  if (up%id==um%id) then
    output = 1
    return
  endif
  
  ! (u+)^j * (u-)^{n-j} = sum_{k=0}^n M(j,k,n) c^k * s^{n-k}.
  !
  ! (u+)^j * (u-)^{n-j} = (c+is)^j * (c-is)^{n-j} / sqrt(2)^n
  !
  !    = sum_{l=0}^j     [ bin(j,l)   * c^l * ( is)^{j-l}   ]
  !    * sum_{m=0}^{n-j} [ bin(n-j,m) * c^m * (-is)^{n-j-m} ]
  !    / sqrt(2)^n
  !
  !    = sum_{l,m} [ bin(j,l) bin(n-j,m) c^{l+m} (is)^{n-l-m} (-1)^{n-j-m} ]
  !    / sqrt(2)^n
  !
  ! k = l+m
  !
  ! 0 <= l <= j
  !
  !    0      <= m   <= n-j
  ! => 0      <= k-l <= n-j
  ! => -k     <= -l  <= n-j-k
  ! => -n+j+k <= l   <= k
  !
  ! => max(0,-n+j+k) <= l <= min(j,k)
  !
  ! => sum_{l,m} = sum_{k=0}^n sum_{l=max(0,-n+j+k)}^{min(j,k)}
  !
  ! => (u+)^j*(u-)^{n-j} =
  !       sum_{k,l} [ bin(j,l)*bin(n-j,k-l) c^k (is)^{n-k} (-1)^{n-j-k+l} ]
  !       / sqrt(2)^n
  !
  ! => M(j,k,n) = i^{n-k} * (-1)^{n-j-k} / sqrt(2)^n
  !             * sum_{l=max(0,j+k-n)}^{min(j,k)} bin(j,l)*bin(n-j,k-l)*(-1)^l
  j = up%power
  k = c%power
  n = up%power + um%power
  output = 0
  do l=max(0,j+k-n),min(j,k)
    output = output + binomial(j,l)*binomial(n-j,k-l)*(-1)**l
  enddo
  output = output                         &
       & * cmplx(0.0_dp,1.0_dp,dp)**(n-k) &
       & * (-1.0_dp)**(n-j-k)             &
       & / sqrt(2.0_dp)**n
end function

impure elemental function element_RealUnivariates_ComplexUnivariates &
   & (this,this_pair,that,that_pair) result(output)
  implicit none
  
  type(RealUnivariate),    intent(in)           :: this
  type(RealUnivariate),    intent(in), optional :: this_pair
  type(ComplexUnivariate), intent(in)           :: that
  type(ComplexUnivariate), intent(in), optional :: that_pair
  complex(dp)                                   :: output
  
  type(ComplexUnivariate) :: up ! u+ = (c  + i*s)/sqrt(2)
  type(ComplexUnivariate) :: um ! u- = (c  - i*s)/sqrt(2)
  type(RealUnivariate)    :: c  ! c  = (u+ + u- )/sqrt(2)
  type(RealUnivariate)    :: s  ! s  = (u+ - u- )/sqrt(2)i
  
  integer :: j,k,l,n
  
  ! Check for optional arguments, and identify which mode is i and which is j.
  if (present(this_pair)) then
    if (this%id<this%paired_id()) then
      c = this
      s = this_pair
    else
      c = this_pair
      s = this
    endif
  else
    if (this%id<this%paired_id()) then
      c = this
      s = RealUnivariate(this%paired_id(),this%id,power=0)
    else
      c = RealUnivariate(this%paired_id(),this%id,power=0)
      s = this
    endif
  endif
  
  if (present(that_pair)) then
    if (that%id<that%paired_id()) then
      up = that
      um = that_pair
    else
      up = that_pair
      um = that
    endif
  else
    if (that%id<that%paired_id()) then
      up = that
      um = ComplexUnivariate(that%paired_id(),that%id,power=0)
    else
      up = ComplexUnivariate(that%paired_id(),that%id,power=0)
      um = that
    endif
  endif
  
  ! Check that inputs are paired up.
  if (up%paired_id()/=um%id .or. up%id/=um%paired_id()) then
    call print_line(CODE_ERROR//': Unpaired modes passed to function.')
    call err()
  elseif (c%paired_id()/=s%id .or. c%id/=s%paired_id()) then
    call print_line(CODE_ERROR//': Unpaired modes passed to function.')
    call err()
  endif
  
  ! If the modes don't match, or the powers don't match, then the element
  !    is zero.
  if (up%id/=c%id .or. um%id/=s%id) then
    output = 0
    return
  elseif (up%power+um%power /= c%power+s%power) then
    output = 0
    return
  endif
  
  ! Check if the mode is real. If so, the overlap is 1.
  if (up%id==um%id) then
    output = 1
    return
  endif
  
  ! c^j * s*{n-j} = sum_{k=0}^n M(j,k,n) (u+)^k * (u-)^{n-k}.
  !
  ! c^j * s^{n-j} = (u+ + u-)^j * ((u+ - u-)/i)^{n-j} / sqrt(2)^n
  !
  !    = sum_{l=0}^j     [ bin(j,l)   * (u+)^l * ( u-)^{j-l}   ]
  !    * sum_{m=0}^{n-j} [ bin(n-j,m) * (u+)^m * (-u-)^{n-j-m} ]
  !    / (sqrt(2)^n * i^{n-j})
  !
  !    = sum_{l,m} [ bin(j,l) bin(n-j,m) (u+)^{l+m} (u-)^{n-l-m} (-1)^{n-j-m} ]
  !    / (sqrt(2)^n * i^{n-j})
  !
  ! k = l+m
  !
  ! 0 <= l <= j
  !
  !    0      <= m   <= n-j
  ! => 0      <= k-l <= n-j
  ! => -k     <= -l  <= n-j-k
  ! => -n+j+k <= l   <= k
  !
  ! => max(0,-n+j+k) <= l <= min(j,k)
  !
  ! => sum_{l,m} = sum_{k=0}^n sum_{l=max(0,-n+j+k)}^{min(j,k)}
  !
  ! => c^j*s^{n-j} =
  !      sum_{k,l} [ bin(j,l) bin(n-j,k-l) (u+)^k (u-)^{n-k} (-1)^{n-j-k+l} ]
  !      / (sqrt(2)^n * i^{n-j})
  !
  ! => M(j,k,n) = i^{j-n} * (-1)^{n-j-k} / sqrt(2)^n
  !             * sum_{l=max(0,j+k-n)}^{min(j,k)} bin(j,l)*bin(n-j,k-l)*(-1)^l
  j = c%power
  k = up%power
  n = c%power + s%power
  output = 0
  do l=max(0,j+k-n),min(j,k)
    output = output + binomial(j,l)*binomial(n-j,k-l)*(-1)**l
  enddo
  output = output                         &
       & * cmplx(0.0_dp,1.0_dp,dp)**(j-n) &
       & * (-1.0_dp)**(n-j-k)             &
       & / sqrt(2.0_dp)**n
end function
end module
