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
! The choice of which mode is u+/c and which is u-/s is arbitrary,
!    but the choice of the mode with the smaller id is u+/c is used throughout.
! 
! u+ = (c  + i*s) /  sqrt(2)
! u- = (c  - i*s) /  sqrt(2)
! c  = (u+ +  u-) /  sqrt(2)
! s  = (u+ -  u-) / (sqrt(2)*i)
!
! u+ = (a+ib)e^{ 2pi i q.R}
! u- = (a-ib)e^{-2pi i q.R}
! c  = sqrt(2)(a*cos(2pi q.R) - b*sin(2pi q.R))
! s  = sqrt(2)(a*sin(2pi q.R) + b*cos(2pi q.R))
! 
! Under q -> -q : u+ ->   u-
!                 u- ->   u+
!                 c  ->   c
!                 s  -> - s
!
! N.B. if i=j, then the mode is real (2q=G, so e^{+iq.r}=e^{-iq.r}=cos(q.r)).
! In this case, u+ = c, and the mode has no pair (u- = s = 0).

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
    module procedure complex_to_real_mode
  end interface
  
  interface real_to_complex
    module procedure real_to_complex_mode
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
impure elemental function complex_to_real_mode(input) result(output)
  implicit none
  
  type(ComplexMode), intent(in) :: input
  type(RealMode)                :: output
  
  type(RealVector), allocatable :: a(:)
  type(RealVector), allocatable :: b(:)
  
  type(RealVector), allocatable :: cos_vector(:)
  type(RealVector), allocatable :: sin_vector(:)
  
  integer :: qpoint_id_plus
  integer :: qpoint_id_minus
  
  if (input%id==input%paired_id) then
    ! output = input is its own pair. It is already real.
    ! c = u+, s = u- = 0.
    qpoint_id_plus = input%qpoint_id
    qpoint_id_minus = qpoint_id_plus
    a = real(input%unit_vector)
    b = a * 0.0_dp
    cos_vector = a
    sin_vector = b
  elseif (input%id<input%paired_id) then
    ! output is c.
    ! input is u+.
    ! u+ = (a+ib)e^{ 2pi i q.R}
    ! c  = sqrt(2)(a*cos(2pi q.R) - b*sin(2pi q.R))
    qpoint_id_plus = input%qpoint_id
    qpoint_id_minus = input%paired_qpoint_id
    a = real( input%unit_vector)
    b = aimag(input%unit_vector)
    cos_vector =  sqrt(2.0_dp) * a
    sin_vector = -sqrt(2.0_dp) * b
  elseif (input%id>input%paired_id) then
    ! output is s.
    ! input is u-.
    ! u- = (a-ib)e^{-2pi i q.R}
    ! s  = sqrt(2)(a*sin(2pi q.R) + b*cos(2pi q.R))
    qpoint_id_plus = input%paired_qpoint_id
    qpoint_id_minus = input%qpoint_id
    a =  real( input%unit_vector)
    b = -aimag(input%unit_vector)
    cos_vector = sqrt(2.0_dp) * b
    sin_vector = sqrt(2.0_dp) * a
  else
    call err()
  endif
  
  output = RealMode(                                  &
     & id                 = input%id,                 &
     & paired_id          = input%paired_id,          &
     & frequency          = input%frequency,          &
     & spring_constant    = input%spring_constant,    &
     & soft_mode          = input%soft_mode,          &
     & translational_mode = input%translational_mode, &
     & cos_vector         = cos_vector,               &
     & sin_vector         = sin_vector,               &
     & qpoint_id_plus     = qpoint_id_plus,           &
     & qpoint_id_minus    = qpoint_id_minus,          &
     & subspace_id        = input%subspace_id         )
end function

impure elemental function real_to_complex_mode(input) result(output)
  implicit none
  
  type(RealMode), intent(in) :: input
  type(ComplexMode)          :: output
  
  type(RealVector), allocatable :: a(:)
  type(RealVector), allocatable :: b(:)
  
  type(ComplexVector), allocatable :: unit_vector(:)
  
  integer :: qpoint_id
  integer :: paired_qpoint_id
  
  ! Construct vectors and q-point ids.
  if (input%id==input%paired_id) then
    ! output = input is its own pair. It is already real.
    ! c = u+, s = u- = 0.
    qpoint_id = input%qpoint_id_plus
    paired_qpoint_id = qpoint_id
    a = input%cos_vector
    unit_vector = cmplxvec(a)
  elseif (input%id<input%paired_id) then
    ! output is u+.
    ! input is c.
    ! c  = sqrt(2)(a*cos(2pi q.R) - b*sin(2pi q.R))
    ! u+ = (a+ib)e^{2pi i q.R}
    qpoint_id = input%qpoint_id_plus
    paired_qpoint_id = input%qpoint_id_minus
    a =  input%cos_vector / sqrt(2.0_dp)
    b = -input%sin_vector / sqrt(2.0_dp)
    unit_vector = cmplxvec(a,b)
  elseif (input%paired_id<input%id) then
    ! output is u-.
    ! input is s.
    ! s  = sqrt(2)(a*sin(2pi q.R) + b*cos(2pi q.R))
    ! u- = (a-ib)e^{-2pi i q.R}
    qpoint_id = input%qpoint_id_minus
    paired_qpoint_id = input%qpoint_id_plus
    a = input%sin_vector / sqrt(2.0_dp)
    b = input%cos_vector / sqrt(2.0_dp)
    unit_vector = cmplxvec(a,-b)
  else
    call err()
  endif
  
  output = ComplexMode(                               &
     & id                 = input%id,                 &
     & paired_id          = input%paired_id,          &
     & frequency          = input%frequency,          &
     & spring_constant    = input%spring_constant,    &
     & soft_mode          = input%soft_mode,          &
     & translational_mode = input%translational_mode, &
     & unit_vector        = unit_vector,              &
     & qpoint_id          = qpoint_id,                &
     & paired_qpoint_id   = paired_qpoint_id,         &
     & subspace_id        = input%subspace_id         )
end function

! Construct the matrix M which converts from one basis set to the other,
!    s.t. this = M * that.
! If include_coefficients is .true. then the coefficients of the monomials are
!    included in the definition of the basis. If .false., then the coefficients
!    are not included.
! e.g. if the real basis is A*u+ and B*u-, and the complex basis is C*c and
!    D*s, where c=(u+)+i(u-) and s=(u+)-i(u-), then:
!    - if .false. then the conversion matrix is (c) = ( 1  i) . (u+)
!                                               (s)   ( 1 -i)   (u-)
!    - if .true. then the conversion matrix is (Cc) = ( C/A  Ci/A) . (Au+)
!                                              (Ds)   ( D/B -Di/B)   (Bu-)
function conversion_matrix_ComplexMonomials_RealMonomials(this,that, &
   & include_coefficients) result(output)
  implicit none
  
  type(ComplexMonomial), intent(in) :: this(:)
  type(RealMonomial),    intent(in) :: that(:)
  logical,               intent(in) :: include_coefficients
  type(ComplexMatrix)               :: output
  
  complex(dp), allocatable :: matrix(:,:)
  
  integer :: i,j,ialloc
  
  allocate(matrix(size(this),size(that)), stat=ialloc); call err(ialloc)
  do i=1,size(that)
    do j=1,size(this)
      matrix(j,i) = element(this(j),that(i),include_coefficients)
    enddo
  enddo
  
  output = matrix
end function

function conversion_matrix_RealMonomials_ComplexMonomials(this,that, &
   & include_coefficients) result(output)
  implicit none
  
  type(RealMonomial),    intent(in) :: this(:)
  type(ComplexMonomial), intent(in) :: that(:)
  logical,               intent(in) :: include_coefficients
  type(ComplexMatrix)               :: output
  
  complex(dp), allocatable :: matrix(:,:)
  
  integer :: i,j,ialloc
  
  allocate(matrix(size(this),size(that)), stat=ialloc); call err(ialloc)
  do i=1,size(that)
    do j=1,size(this)
      matrix(j,i) = element(this(j),that(i),include_coefficients)
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
impure elemental function element_ComplexMonomial_RealMonomial(this,that, &
   & include_coefficients) result(output)
  implicit none
  
  type(ComplexMonomial), intent(in) :: this
  type(RealMonomial),    intent(in) :: that
  logical,               intent(in) :: include_coefficients
  complex(dp)                       :: output
  
  logical, allocatable :: that_mode_accounted(:)
  
  integer :: i,j,k,l
  
  if ( sum(this%modes%power+this%modes%paired_power) /= &
     & sum(that%modes%power+that%modes%paired_power)    ) then
    output = 0
    return
  endif
  
  output = 1.0_dp
  
  that_mode_accounted = [(.false., i=1, size(that))]
  
  ! Loop over modes, multiplying the output by the overlap of each mode pair
  !    in 'this' with the eqivalent mode pair in 'that'.
  do i=1,size(this)
    ! Find the location of this%modes(i) in 'that'.
    j = first(that%modes%id==this%modes(i)%id, default=0)
    
    if (j==0) then
      ! If the modes doesn't exist in that, and the mode in 'this' has
      !    non-zero power, then the overlap between 'this' and 'that' is zero.
      ! If the mode in 'this' has zero power, then the mode has zero power
      !    in both, and can be neglected.
      if (this%modes(i)%power/=0 .or. this%modes(i)%paired_power/=0) then
        output = 0
        return
      endif
    else
      ! If the mode exists in both, account for it in the overlap.
      output = output * element(this%modes(i), that%modes(j))
      that_mode_accounted(j) = .true.
    endif
  enddo
  
  ! Check remaining modes in 'that'. If any have non-zero power then the
  !    overlap is zero.
  do i=1,size(that)
    if (.not. that_mode_accounted(j)) then
      if (that%modes(j)%power/=0 .or. that%modes(j)%paired_power/=0) then
        output = 0
        return
      endif
    endif
  enddo
  
  if (include_coefficients) then
    output = output * this%coefficient / that%coefficient
  endif
end function

! Calculate an element for conversion_matrix.
!
! Calculates the coefficient of 'that' in the expansion of 'this'.
! i.e. 'this' = sum[ element(this,that) * 'that' ].
impure elemental function element_RealMonomial_ComplexMonomial(this,that, &
   & include_coefficients) result(output)
  implicit none
  
  type(RealMonomial),    intent(in) :: this
  type(ComplexMonomial), intent(in) :: that
  logical,               intent(in) :: include_coefficients
  complex(dp)                       :: output
  
  logical, allocatable :: that_mode_accounted(:)
  
  integer :: i,j,k,l
  
  if ( sum(this%modes%power+this%modes%paired_power) /= &
     & sum(that%modes%power+that%modes%paired_power)    ) then
    output = 0
    return
  endif
  
  output = 1.0_dp
  
  that_mode_accounted = [(.false., i=1, size(that))]
  
  ! Loop over modes, multiplying the output by the overlap of each mode pair
  !    in 'this' with the eqivalent mode pair in 'that'.
  do i=1,size(this)
    ! Find the location of this%modes(i) in 'that'.
    j = first(that%modes%id==this%modes(i)%id, default=0)
    
    if (j==0) then
      ! If the modes doesn't exist in that, and the mode in 'this' has
      !    non-zero power, then the overlap between 'this' and 'that' is zero.
      ! If the mode in 'this' has zero power, then the mode has zero power
      !    in both, and can be neglected.
      if (this%modes(i)%power/=0 .or. this%modes(i)%paired_power/=0) then
        output = 0
        return
      endif
    else
      ! If the mode exists in both, account for it in the overlap.
      output = output * element(this%modes(i), that%modes(j))
      that_mode_accounted(j) = .true.
    endif
  enddo
  
  ! Check remaining modes in 'that'. If any have non-zero power then the
  !    overlap is zero.
  do i=1,size(that)
    if (.not. that_mode_accounted(j)) then
      if (that%modes(j)%power/=0 .or. that%modes(j)%paired_power/=0) then
        output = 0
        return
      endif
    endif
  enddo
  
  if (include_coefficients) then
    output = output * this%coefficient / that%coefficient
  endif
end function

impure elemental function element_ComplexUnivariates_RealUnivariates(this, &
   & that) result(output)
  implicit none
  
  type(ComplexUnivariate), intent(in) :: this
  type(RealUnivariate),    intent(in) :: that
  complex(dp)                         :: output
  
  integer :: j,k,l,n
  
  ! Check that inputs are paired up.
  if (this%id/=that%id) then
    call print_line(CODE_ERROR//': Unpaired modes passed to function.')
    call err()
  endif
  
  ! If the powers don't match, then the overlap is zero.
  if (this%power+this%paired_power /= that%power+that%paired_power) then
    output = 0
    return
  endif
  
  ! Check if the mode is real. If so, the overlap is 1.
  if (this%id==this%paired_id) then
    output = 1
    return
  endif
  
  ! this = (u+)^j * (u-)^{n-j}.
  ! that = (c )^k * (s )^{n-k}.
  !
  ! (u+)^j * (u-)^{n-j} = sum_{k=0}^n M(j,k,n) c^k * s*{n-k}.
  !
  ! (u+)^j * (u-)^{n-j} = (c+is)^j * (c-is)^(n-j) / sqrt(2)^n
  !
  ! => this = sum_{l=0}^j sum_{m=0}^{n-j}
  !           [ bin(j,l) bin(n-j,m) i^(-n+2j-l+m) c^{l+m} s^{n-l-m} ]
  !           / sqrt(2)^n
  !
  ! k = l+m
  !
  ! => sum_{l=0}^n sum_{m=0}^{n-j} = sum_{k=0}^n sum_{l=max(0,j+k-n)}^{min(j,k)}
  !
  ! => that = sum_{k,l} [ bin(j,l)*bin(n-j,k-l) i^(-n+2j+k-2l) c^k s^{n-k} ]
  !         / sqrt(2)^n
  !
  ! => M(j,k,n) = i^{-n+2j+k} / sqrt(2)^n
  !             * sum_{l=max(0,j+k-n)}^{min(j,k)} bin(j,l)*bin(n-j,k-l)*(-1)^l
  j = this%power
  k = that%power
  n = this%power + this%paired_power
  output = 0
  do l=max(0,j+k-n),min(j,k)
    output = output + binomial(j,l)*binomial(n-j,k-l)*(-1)**l
  enddo
  output = output                              &
       & * cmplx(0.0_dp,1.0_dp,dp)**(-n+2*j+k) &
       & / sqrt(2.0_dp)**n
end function

impure elemental function element_RealUnivariates_ComplexUnivariates(this, &
   & that) result(output)
  implicit none
  
  type(RealUnivariate),    intent(in) :: this
  type(ComplexUnivariate), intent(in) :: that
  complex(dp)                         :: output
  
  integer :: j,k,l,n
  
  ! Check that inputs are paired up.
  if (this%id/=that%id) then
    call print_line(CODE_ERROR//': Unpaired modes passed to function.')
    call err()
  endif
  
  ! If the powers don't match, then the overlap is zero.
  if (this%power+this%paired_power /= that%power+that%paired_power) then
    output = 0
    return
  endif
  
  ! Check if the mode is real. If so, the overlap is 1.
  if (this%id==this%paired_id) then
    output = 1
    return
  endif
  
  ! this = (c )^j * (s )^{n-j}.
  ! that = (u+)^k * (u-)^{n-k}.
  !
  ! c^j * s^{n-j} = sum_{k=0}^n M(j,k,n) (u+)^k * (u-)^{n-k}.
  !
  ! c^j * s^{n-j} = (u+ + u-)^j * ((u+ - u-)/i)^{n-j} / sqrt(2)^n
  !
  !    = sum_{l=0}^j sum_{m=0}^{n-j}
  !      [ bin(j,l) bin(n-j,m) (u+)^{l+m} (u-)^{n-l-m} (-1)^{n-j-m} ]
  !    / (sqrt(2)^n * i^{n-j})
  !
  ! k = l+m
  !
  ! => sum_{l=0}^n sum_{m=0}^{n-j} = sum_{k=0}^n sum_{l=max(0,j+k-n)}^{min(j,k)}
  !
  ! => that = sum_{k,l}[ bin(j,l) bin(n-j,k-l) u+^k u-^{n-k} (-1)^{n-j-k+l} ]
  !      / (sqrt(2)^n * i^{n-j})
  !
  ! => M(j,k,n) = i^{n-j-2k} / sqrt(2)^n
  !             * sum_{l=max(0,j+k-n)}^{min(j,k)} bin(j,l)*bin(n-j,k-l)*(-1)^l
  j = this%power
  k = that%power
  n = this%power + this%paired_power
  output = 0
  do l=max(0,j+k-n),min(j,k)
    output = output + binomial(j,l)*binomial(n-j,k-l)*(-1)**l
  enddo
  output = output                             &
       & * cmplx(0.0_dp,1.0_dp,dp)**(n-j-2*k) &
       & / sqrt(2.0_dp)**n
end function
end module
