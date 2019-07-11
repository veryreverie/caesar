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

module real_complex_conversion_module
  use utils_module
  
  use complex_mode_module
  use real_mode_module
  use complex_polynomial_module
  use real_polynomial_module
  implicit none
  
  private
  
  public :: complex_to_real
  public :: real_to_complex
  public :: basis_conversion_matrix
  public :: coefficient_conversion_matrix
  
  interface complex_to_real
    module procedure complex_to_real_mode
  end interface
  
  interface real_to_complex
    module procedure real_to_complex_mode
  end interface
  
  interface basis_conversion_matrix
    module procedure basis_conversion_matrix_ComplexMonomials_RealMonomials
    module procedure basis_conversion_matrix_RealMonomials_ComplexMonomials
  end interface
  
  interface coefficient_conversion_matrix
    module procedure coefficient_conversion_matrix_ComplexMonomials_RealMonomials
    module procedure coefficient_conversion_matrix_RealMonomials_ComplexMonomials
  end interface
  
  interface element
    module procedure element_ComplexMonomial_RealMonomial
    module procedure element_RealMonomial_ComplexMonomial
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
function basis_conversion_matrix_ComplexMonomials_RealMonomials(this,that, &
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
  
  output = mat(matrix)
end function

function basis_conversion_matrix_RealMonomials_ComplexMonomials(this,that, &
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
  
  output = mat(matrix)
end function

! ----------------------------------------------------------------------
! Conversion matrices for coefficients rather than basis functions.
! ----------------------------------------------------------------------
! Coefficients transform in the opposite manner to basis functions.
! If f(x) = sum_i a_ig_i(x) = sum_j b_jh_j(x) then if
!    f_i(x) = sum_j C_ij h_j(x)
!    h_j(x) = sum_i D_ji g_j(x)
! where sum_j C_ij D_jk is the identity, then
!    a_i    = sum_j D_ji b_j
!    b_j    = sum_j C_ij a_i
! So the conversion matrices are D^T and C^T rather than C and D.
function coefficient_conversion_matrix_ComplexMonomials_RealMonomials(this, &
   & that,include_coefficients) result(output)
  implicit none
  
  type(ComplexMonomial), intent(in) :: this(:)
  type(RealMonomial),    intent(in) :: that(:)
  logical,               intent(in) :: include_coefficients
  type(ComplexMatrix)               :: output
  
  output = transpose(basis_conversion_matrix(that,this,include_coefficients))
end function

function coefficient_conversion_matrix_RealMonomials_ComplexMonomials(this, &
   & that,include_coefficients) result(output)
  implicit none
  
  type(RealMonomial),    intent(in) :: this(:)
  type(ComplexMonomial), intent(in) :: that(:)
  logical,               intent(in) :: include_coefficients
  type(ComplexMatrix)               :: output
  
  output = transpose(basis_conversion_matrix(that,this,include_coefficients))
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
  
  type(ComplexUnivariate) :: this_mode
  type(RealUnivariate)    :: that_mode
  
  logical, allocatable :: that_mode_accounted(:)
  
  complex(dp) :: overlap
  
  integer :: i,j,p,q,l,n
  
  if (this%total_power()/=that%total_power()) then
    output = 0
    return
  endif
  
  output = 1.0_dp
  
  that_mode_accounted = [(.false., i=1, size(that))]
  
  ! Loop over modes, multiplying the output by the overlap of each mode pair
  !    in 'this' with the eqivalent mode pair in 'that'.
  do i=1,size(this)
    this_mode = this%mode(i)
    
    ! Find the location of this%modes(i) in 'that'.
    j = first_equivalent(that%ids(), this%id(i), default=0, sorted=.true.)
    
    if (j==0) then
      ! If the modes doesn't exist in that, and the mode in 'this' has
      !    non-zero power, then the overlap between 'this' and 'that' is zero.
      ! If the mode in 'this' has zero power, then the mode has zero power
      !    in both, and can be neglected.
      if (this_mode%power/=0 .or. this_mode%paired_power/=0) then
        output = 0
        return
      endif
    else
      that_mode = that%mode(j)
      ! If the mode exists in both, account for it in the overlap.
      if (this_mode%total_power()/=that_mode%total_power()) then
        ! There is no overlap between this and that.
        output = 0
        return
      elseif (this_mode%id==this_mode%paired_id) then
        ! The mode is real. The overlap is 1.
        continue
      else
        ! this_mode = (u+)^p * (u-)^{n-p}.
        ! that_mode = (c )^q * (s )^{n-q}.
        !
        ! (u+)^p * (u-)^{n-p} = sum_{q=0}^n M(p,q,n) c^q * s*{n-q}.
        !
        ! (u+)^p * (u-)^{n-p} = (c+is)^p * (c-is)^(n-p) / sqrt(2)^n
        !
        ! => this = sum_{l=0}^p sum_{m=0}^{n-p}
        !           [ bin(p,l) bin(n-p,m) i^(-n+2j-l+m) c^{l+m} s^{n-l-m} ]
        !           / sqrt(2)^n
        !
        ! q = l+m
        !
        ! => sum_{l=0}^n sum_{m=0}^{n-p} = sum_{q=0}^n sum_{l=max(0,p+q-n)}^{min(p,q)}
        !
        ! => that = sum_{q,l} [ bin(p,l)*bin(n-p,q-l) i^(-n+2j+q-2l) c^q s^{n-q} ]
        !         / sqrt(2)^n
        !
        ! => M(p,q,n) = i^{-n+2j+q} / sqrt(2)^n
        !             * sum_{l=max(0,p+q-n)}^{min(p,q)} bin(p,l)*bin(n-p,q-l)*(-1)^l
        p = this_mode%power
        q = that_mode%power
        n = this_mode%total_power()
        overlap = 0
        do l=max(0,p+q-n),min(p,q)
          overlap = overlap                                      &
                & + exp(log_binomial(p,l)+log_binomial(n-p,q-l)) &
                & * (-1)**l
        enddo
        overlap = overlap                            &
             & * cmplx(0.0_dp,1.0_dp,dp)**(-n+2*p+q) &
             & / sqrt(2.0_dp)**n
        output = output * overlap
      endif
      
      that_mode_accounted(j) = .true.
    endif
  enddo
  
  ! Check remaining modes in 'that'. If any have non-zero power then the
  !    overlap is zero.
  do i=1,size(that)
    if (.not. that_mode_accounted(j)) then
      if (that%power(j)/=0 .or. that%paired_power(j)/=0) then
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
  
  type(RealUnivariate)    :: this_mode
  type(ComplexUnivariate) :: that_mode
  
  logical, allocatable :: that_mode_accounted(:)
  
  complex(dp) :: overlap
  
  integer :: i,j,p,q,l,n
  
  if (this%total_power()/=that%total_power()) then
    output = 0
    return
  endif
  
  output = 1.0_dp
  
  that_mode_accounted = [(.false., i=1, size(that))]
  
  ! Loop over modes, multiplying the output by the overlap of each mode pair
  !    in 'this' with the eqivalent mode pair in 'that'.
  do i=1,size(this)
    this_mode = this%mode(i)
    
    ! Find the location of this%modes(i) in 'that'.
    j = first_equivalent(that%ids(), this%id(i), default=0, sorted=.true.)
    
    if (j==0) then
      ! If the modes doesn't exist in that, and the mode in 'this' has
      !    non-zero power, then the overlap between 'this' and 'that' is zero.
      ! If the mode in 'this' has zero power, then the mode has zero power
      !    in both, and can be neglected.
      if (this%power(i)/=0 .or. this%paired_power(i)/=0) then
        output = 0
        return
      endif
    else
      that_mode = that%mode(j)
      
      ! If the mode exists in both, account for it in the overlap.
      if (this_mode%total_power()/=that_mode%total_power()) then
        ! There is no overlap between this and that.
        output = 0
        return
      elseif (this_mode%id==this_mode%paired_id) then
        ! The mode is real. The overlap is 1.
        continue
      else
        ! this_mode = (c )^p * (s )^{n-p}.
        ! that_mode = (u+)^q * (u-)^{n-q}.
        !
        ! c^p * s^{n-p} = sum_{q=0}^n M(p,q,n) (u+)^q * (u-)^{n-q}.
        !
        ! c^p * s^{n-p} = (u+ + u-)^p * ((u+ - u-)/i)^{n-p} / sqrt(2)^n
        !
        !    = sum_{l=0}^p sum_{m=0}^{n-p}
        !      [ bin(p,l) bin(n-p,m) (u+)^{l+m} (u-)^{n-l-m} (-1)^{n-p-m} ]
        !    / (sqrt(2)^n * i^{n-p})
        !
        ! q = l+m
        !
        ! => sum_{l=0}^n sum_{m=0}^{n-p} = sum_{q=0}^n sum_{l=max(0,p+q-n)}^{min(p,q)}
        !
        ! => that = sum_{q,l}[ bin(p,l) bin(n-p,q-l) u+^q u-^{n-q} (-1)^{n-p-q+l} ]
        !      / (sqrt(2)^n * i^{n-p})
        !
        ! => M(p,q,n) = i^{n-p-2k} / sqrt(2)^n
        !             * sum_{l=max(0,p+q-n)}^{min(p,q)} bin(p,l)*bin(n-p,q-l)*(-1)^l
        p = this_mode%power
        q = that_mode%power
        n = this_mode%total_power()
        overlap = 0
        do l=max(0,p+q-n),min(p,q)
          overlap = overlap                                      &
                & + exp(log_binomial(p,l)+log_binomial(n-p,q-l)) &
                & * (-1)**l
        enddo
        overlap = overlap                             &
             & * cmplx(0.0_dp,1.0_dp,dp)**(n-p-2*q) &
             & / sqrt(2.0_dp)**n
        output = output * overlap
      endif
      
      that_mode_accounted(j) = .true.
    endif
  enddo
  
  ! Check remaining modes in 'that'. If any have non-zero power then the
  !    overlap is zero.
  do i=1,size(that)
    if (.not. that_mode_accounted(j)) then
      if (that%power(j)/=0 .or. that%paired_power(j)/=0) then
        output = 0
        return
      endif
    endif
  enddo
  
  if (include_coefficients) then
    output = output * this%coefficient / that%coefficient
  endif
end function
end module
