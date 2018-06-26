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
! Under q -> -q : u+ ->   u-
!                 u- ->   u+
!                 c  ->   c
!                 s  -> - s
!
! N.B. if i=j, then the mode is real (2q=G, so e^{+iq.r}=e^{-iq.r}=cos(q.r)).

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
  
  type(RealVector), allocatable :: mass_weighted_vector(:)
  type(RealVector), allocatable :: cartesian_vector(:)
  
  integer :: qpoint_id_plus
  integer :: qpoint_id_minus
  
  integer :: mode,ialloc
  
  allocate(output(size(input)), stat=ialloc); call err(ialloc)
  do mode=1,size(input)
    if (input(mode)%id==input(mode)%paired_id) then
      ! output(mode) is its own pair. It is already real.
      ! c = u+.
      qpoint_id_plus = input(mode)%qpoint_id
      qpoint_id_minus = qpoint_id_plus
      mass_weighted_vector = real(input(mode)%mass_weighted_vector)
      cartesian_vector = real(input(mode)%cartesian_vector)
    elseif (input(mode)%id<input(mode)%paired_id) then
      ! output(mode) is c.
      ! input(mode) is u+.
      ! c = real(u+)*sqrt(2).
      qpoint_id_plus = input(mode)%qpoint_id
      qpoint_id_minus = input(mode)%paired_qpoint_id
      mass_weighted_vector = real(input(mode)%mass_weighted_vector) &
                         & * sqrt(2.0_dp)
      cartesian_vector = real(input(mode)%cartesian_vector) * sqrt(2.0_dp)
    elseif (input(mode)%id>input(mode)%paired_id) then
      ! output(mode) is s.
      ! input(mode) is u-.
      ! s = -imag(u-)*sqrt(2).
      mass_weighted_vector = -aimag(input(mode)%mass_weighted_vector) &
                         & * sqrt(2.0_dp)
      cartesian_vector = -aimag(input(mode)%cartesian_vector) &
                     & * sqrt(2.0_dp)
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
       & mass_weighted_vector = mass_weighted_vector,           &
       & cartesian_vector     = cartesian_vector,               &
       & qpoint_id_plus       = qpoint_id_plus,                 &
       & qpoint_id_minus      = qpoint_id_minus,                &
       & subspace_id          = input(mode)%subspace_id)
  enddo
end function

function real_to_complex_Modes(input) result(output)
  implicit none
  
  type(RealMode), intent(in)     :: input(:)
  type(ComplexMode), allocatable :: output(:)
  
  type(ComplexVector), allocatable :: mass_weighted_vector(:)
  type(ComplexVector), allocatable :: cartesian_vector(:)
  
  integer :: qpoint_id
  integer :: paired_qpoint_id
  
  integer :: mode,pair,ialloc
  
  allocate(output(size(input)), stat=ialloc); call err(ialloc)
  do mode=1,size(input)
    ! Get the location of the paired mode.
    pair = first(input%id==input(mode)%paired_id)
    
    ! Construct vectors and q-point ids.
    if (input(mode)%id==input(pair)%id) then
      ! output(mode) is its own pair. It is real.
      ! u+ = c.
      qpoint_id = input(mode)%qpoint_id_plus
      paired_qpoint_id = qpoint_id
      mass_weighted_vector = cmplxvec(input(mode)%mass_weighted_vector)
      cartesian_vector = cmplxvec(input(mode)%cartesian_vector)
    elseif (input(mode)%id<input(pair)%id) then
      ! output(mode) is u+.
      ! input(mode) is c, and input(pair) is s.
      ! u+ = (c + i*s)/sqrt(2).
      qpoint_id = input(mode)%qpoint_id_plus
      paired_qpoint_id = input(mode)%qpoint_id_minus
      mass_weighted_vector =                           &
         & cmplxvec( input(mode)%mass_weighted_vector, &
         &           input(pair)%mass_weighted_vector) &
         & / sqrt(2.0_dp)
      cartesian_vector =                           &
         & cmplxvec( input(mode)%cartesian_vector, &
         &           input(pair)%cartesian_vector) &
         & / sqrt(2.0_dp)
    elseif (input(mode)%paired_id<input(mode)%id) then
      ! output(mode) is u-.
      ! input(mode) is s, and input(pair) is c.
      ! u- = (c - i*s)/sqrt(2).
      qpoint_id = input(mode)%qpoint_id_minus
      paired_qpoint_id = input(mode)%qpoint_id_plus
      mass_weighted_vector =                           &
         & cmplxvec( input(pair)%mass_weighted_vector, &
         &         - input(mode)%mass_weighted_vector) &
         & / sqrt(2.0_dp)
      cartesian_vector =                           &
         & cmplxvec( input(pair)%cartesian_vector, &
         &         - input(mode)%cartesian_vector) &
         & / sqrt(2.0_dp)
    else
      call err()
    endif
    
    output(mode) = ComplexMode(                                 &
       & id                   = input(mode)%id,                 &
       & paired_id            = input(mode)%paired_id,          &
       & frequency            = input(mode)%frequency,          &
       & spring_constant      = input(mode)%spring_constant,    &
       & soft_mode            = input(mode)%soft_mode,          &
       & translational_mode   = input(mode)%translational_mode, &
       & mass_weighted_vector = mass_weighted_vector,           &
       & cartesian_vector     = cartesian_vector,               &
       & qpoint_id            = qpoint_id,                      &
       & paired_qpoint_id     = paired_qpoint_id,               &
       & subspace_id          = input(mode)%subspace_id)
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

impure elemental function element_RealMonomial_ComplexMonomial(this,that) &
   & result(output)
  implicit none
  
  type(RealMonomial),    intent(in) :: this
  type(ComplexMonomial), intent(in) :: that
  complex(dp)                       :: output
  
  output = conjg(element(that,this))
end function

! Calculates a single-mode or two-mode factor.
!
! The single-mode factor is the Univariant equivalent of the Monomial function
!    above, i.e. 'this' = sum[ element(this,that) * 'that' ].
!
! The two-mode factor takes the Univariates representing the mode and its
!    pair under q -> -q on both sides, so
! 'this'*'this_pair' = sum[ element(this,pair,that,pair) * 'that'*'that_pair' ]

!impure elemental function element_ComplexUnivariate_RealUnivariate(this,that) &
!   & result(output)
!  implicit none
!  
!  type(ComplexUnivariate), intent(in) :: this
!  type(RealUnivariate),    intent(in) :: that
!  complex(dp)                         :: output
!  
!  ! Check for consistency between modes.
!  ! The array [this%id,this%paired_id,that%id,that%paired_id] must either be
!  !    [x,x,x,x] for 'this' and 'that' are the same real mode,
!  !    [x,y,x,y] or [x,y,y,x] for 'this' and 'that' are the same complex mode,
!  ! or [x,x,y,y] or [x,y,z,a] for 'this' and 'that' are different modes.
!  if (this%id==that%id .or. this%id==that%paired_id) then
!    if (this%id==this%paired_id .neqv. that%id==that%paired_id) then
!      call print_line(CODE_ERROR//': Modes paired inconsistently.')
!      call err()
!    elseif (this%id==that%id .neqv. this%paired_id==that%paired_id) then
!      call print_line(CODE_ERROR//': Modes paired inconsistently.')
!      call err()
!    elseif (this%id==that%paired_id .neqv. this%paired_id==that%id) then
!      call print_line(CODE_ERROR//': Modes paired inconsistently.')
!      call err()
!    endif
!  endif
!  
!  if (all(this%id/=[that%id,that%paired_id]) .or. this%power/=that%power) then
!    ! The modes don't match, or their powers don't match.
!    ! Either way, the overlap is 0.
!    output = 0
!  elseif (this%id==this%paired_id) then
!    ! The mode is real. The overlap is 1.
!    output = 1
!  elseif (that%id<that%paired_id) then
!    ! 'that' is a cos mode.
!    ! c_i = (r_i + i.r_j)/sqrt(2).
!    ! c_j = (r_i - i.r_j)/sqrt(2).
!    ! In both cases the r_i coefficient is 1/sqrt(2).
!    ! Thus the (r_i)^power coefficient is (1/sqrt(2))^power.
!    output = 1/sqrt(2.0_dp)**this%power
!  elseif (that%id>that%paired_id .and. this%id<that%paired_id) then
!    ! 'that' is a sin mode and 'this' is an e^{+iq.r} mode.
!    ! c_i = (r_i + i.r_j)/sqrt(2).
!    ! The r_j coefficient is i/sqrt(2).
!    ! Thus the (r_j)^power coefficient is (i/sqrt(2))^power.
!    output = cmplx(0.0_dp,1/sqrt(2.0_dp),dp)**this%power
!  elseif (that%id>that%paired_id .and. this%id>that%paired_id) then
!    ! 'that' is a sin mode and 'this' is an e^{-iq.r} mode.
!    ! c_j = (r_i - i.r_j)/sqrt(2).
!    ! The r_j coefficient is -i/sqrt(2).
!    ! Thus the (r_j)^power coefficient is (-i/sqrt(2))^power.
!    output = cmplx(0.0_dp,-1/sqrt(2.0_dp),dp)**this%power
!  else
!    ! The above options should be exhaustive. Something has gone wrong.
!    call err()
!  endif
!end function

impure elemental function element_ComplexUnivariates_RealUnivariates &
   & (this,this_pair,that,that_pair) result(output)
  implicit none
  
  type(ComplexUnivariate), intent(in)           :: this
  type(ComplexUnivariate), intent(in), optional :: this_pair
  type(RealUnivariate),    intent(in)           :: that
  type(RealUnivariate),    intent(in), optional :: that_pair
  complex(dp)                                   :: output
  
  type(ComplexUnivariate) :: c_i
  type(ComplexUnivariate) :: c_j
  type(RealUnivariate)    :: r_i
  type(RealUnivariate)    :: r_j
  
  integer :: a,b,c,n
  
  ! Check for optional arguments, and identify which mode is i and which is j.
  if (present(this_pair)) then
    if (this%id<this%paired_id()) then
      c_i = this
      c_j = this_pair
    else
      c_i = this_pair
      c_j = this
    endif
  else
    if (this%id<this%paired_id()) then
      c_i = this
      c_j = ComplexUnivariate(this%paired_id(),this%id,power=0)
    else
      c_i = ComplexUnivariate(this%paired_id(),this%id,power=0)
      c_j = this
    endif
  endif
  
  if (present(that_pair)) then
    if (that%id<that%paired_id()) then
      r_i = that
      r_j = that_pair
    else
      r_i = that_pair
      r_j = that
    endif
  else
    if (that%id<that%paired_id()) then
      r_i = that
      r_j = RealUnivariate(that%paired_id(),that%id,power=0)
    else
      r_i = RealUnivariate(that%paired_id(),that%id,power=0)
      r_j = that
    endif
  endif
  
  ! Check that inputs are paired up.
  if (c_i%paired_id()/=c_j%id .or. c_i%id/=c_j%paired_id()) then
    call print_line(CODE_ERROR//': Unpaired modes passed to function.')
    call err()
  elseif (r_i%paired_id()/=r_j%id .or. r_i%id/=r_j%paired_id()) then
    call print_line(CODE_ERROR//': Unpaired modes passed to function.')
    call err()
  endif
  
  ! If the modes don't match, or the powers don't match, then the element
  !    is zero.
  if (c_i%id/=r_i%id .or. c_j%id/=r_j%id) then
    output = 0
    return
  elseif (c_i%power+c_j%power /= r_i%power+r_j%power) then
    output = 0
    return
  endif
  
  ! Check if the mode is real. If so, the overlap is 1.
  if (c_i%id==c_j%id) then
    output = 1
    return
  endif
  
  ! (c_i)^a(c_j)^{n-a} = sum_{b=0}^n M(a,b,n) (r_i)^b(r_j)^{n-b}.
  !
  ! M(a,b,n) = i^{-(n-b)}
  !          * (1/sqrt(2))^n
  !          * sum_{c-max(0,a+b-n)}^{min(a,b)} bin(a,c)*bin(n-a,b-c)*(-1)^{a-c}
  a = c_i%power
  b = r_i%power
  n = c_i%power + c_j%power
  output = 0
  do c=max(0,a+b-n),min(a,b)
    output = output + binomial(a,c)*binomial(n-a,b-c)*(-1)**(a-c)
  enddo
  output = output                            &
       & * cmplx(0.0_dp,1.0_dp,dp)**(-(n-b)) &
       & / sqrt(2.0_dp)**n
end function

impure elemental function element_RealUnivariates_ComplexUnivariates(this, &
   & this_pair,that,that_pair) result(output)
  implicit none
  
  type(RealUnivariate),    intent(in)           :: this
  type(RealUnivariate),    intent(in), optional :: this_pair
  type(ComplexUnivariate), intent(in)           :: that
  type(ComplexUnivariate), intent(in), optional :: that_pair
  complex(dp)                                   :: output
  
  output = conjg(element(that,that_pair,this,this_pair))
end function
end module
