! ======================================================================
! Converts objects between complex and real co-ordinates.
! ======================================================================

! Complex modes with ids i and j, with i<j, are paired under the transformation
!    q -> -q if they transform as u_i -> u_j and u_j -> u_i.
! N.B. u_i = u_j*.
! These modes can then be transformed to and from the real modes i and j
!    as follows.
! The complex u_i will be written c_i, and the real u_i will be written r_i.
! 
! c_i = (r_i + i.r_j) /  sqrt(2)
! c_j = (r_i - i.r_j) /  sqrt(2)
! r_i = (c_i +   c_j) /  sqrt(2)    = real(c_i)*sqrt(2) =  real(c_j)*sqrt(2)
! r_j = (c_i -   c_j) / (sqrt(2)*i) = imag(c_i)*sqrt(2) = -imag(c_j)*sqrt(2)
! 
! Under q -> -q : (c_i) -> (0  1) . (c_i)
!                 (c_j)    (1  0)   (c_j)
!
!                 (r_i) -> (1  0) . (r_i)
!                 (r_j)    (0 -1)   (r_j)
!
! N.B. Taking q to be the q-point of c_i, c_i looks like e^{+iq.r},
!                                         c_j looks like e^{-iq.r},
!                                         r_i looks like cos(q.r),
!                                         r_j looks like sin(q.r).
!
! N.B. the choice of r_i as the cos mode and r_j as the sin mode is arbitrary,
!    but will be used consistently throughout.
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
  
  integer :: mode,ialloc
  
  ! Allocate output, and copy over everything which is unchanged by the
  !    transformation to real co-ordinates.
  allocate(output(size(input)), stat=ialloc); call err(ialloc)
  output%id                 = input%id
  output%paired_id          = input%paired_id
  output%frequency          = input%frequency
  output%soft_mode          = input%soft_mode
  output%translational_mode = input%translational_mode
  output%qpoint_id          = input%qpoint_id
  output%subspace_id        = input%subspace_id
  
  ! Convert displacements to real co-ordinates.
  do mode=1,size(input)
    if (input(mode)%id==input(mode)%paired_id) then
      ! This mode is its own pair. It is already real.
      ! mode = i = j.
      ! r_i = c_i.
      output(mode)%primitive_displacements = &
         & real(input(mode)%primitive_displacements)
    elseif (input(mode)%id<input(mode)%paired_id) then
      ! Construct the cos(q.r) mode, r_i.
      ! mode = i.
      ! r_i = real(c_i)*sqrt(2).
      output(mode)%primitive_displacements = &
         & real(input(mode)%primitive_displacements)*sqrt(2.0_dp)
    elseif (input(mode)%id>input(mode)%paired_id) then
      ! Construct the sin(q.r) mode, r_j.
      ! mode = j.
      ! r_j = -imag(c_j)*sqrt(2).
      output(mode)%primitive_displacements = &
         & -aimag(input(mode)%primitive_displacements)*sqrt(2.0_dp)
    else
      call err()
    endif
  enddo
end function

function real_to_complex_Modes(input) result(output)
  implicit none
  
  type(RealMode), intent(in)     :: input(:)
  type(ComplexMode), allocatable :: output(:)
  
  integer :: i,j,mode,pair,ialloc
  
  ! Allocate output, and copy over everything which is unchanged by the
  !    transformation to real co-ordinates.
  allocate(output(size(input)), stat=ialloc); call err(ialloc)
  output%id                 = input%id
  output%paired_id          = input%paired_id
  output%frequency          = input%frequency
  output%soft_mode          = input%soft_mode
  output%translational_mode = input%translational_mode
  output%qpoint_id          = input%qpoint_id
  output%subspace_id        = input%subspace_id
  
  ! Convert displacements to complex co-ordinates.
  do mode=1,size(input)
    ! Get the location of the paired mode with id j.
    pair = first(input%id==input(mode)%paired_id)
    if (input(mode)%id==input(pair)%id) then
      ! This mode is its own pair. It is real.
      ! mode = pair = i = j.
      ! c_i = r_i.
      output(mode)%primitive_displacements = &
         & cmplxvec(input(mode)%primitive_displacements)
    elseif (input(mode)%id<input(pair)%id) then
      ! Construct the e^{+iq.r} mode, c_i.
      ! mode = i, pair = j.
      ! c_i = (r_i + i.r_j)/sqrt(2).
      output(mode)%primitive_displacements =              &
         & cmplxvec( input(mode)%primitive_displacements, &
         &           input(pair)%primitive_displacements) &
         & / sqrt(2.0_dp)
    elseif (input(mode)%paired_id<input(mode)%id) then
      ! Construct the e^{-iq.r} mode, c_j.
      ! mode = j, pair = i.
      ! c_j = (r_i - i.r_j)/sqrt(2).
      output(mode)%primitive_displacements =              &
         & cmplxvec( input(pair)%primitive_displacements, &
         &          -input(mode)%primitive_displacements) &
         & / sqrt(2.0_dp)
    else
      call err()
    endif
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
  
  integer :: c_i_power
  integer :: c_j_power
  integer :: r_i_power
  integer :: r_j_power
  
  integer :: i,j,k,l,ialloc
  
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
    j = first(this%modes%paired_id == this%modes(i)%id, default=0)
    k = first(that%modes%id        == this%modes(i)%id, default=0)
    l = first(that%modes%paired_id == this%modes(i)%id, default=0)
    
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
    if (this%id<this%paired_id) then
      c_i = this
      c_j = this_pair
    else
      c_i = this_pair
      c_j = this
    endif
  else
    if (this%id<this%paired_id) then
      c_i = this
      c_j = ComplexUnivariate(this%paired_id,this%id,0)
    else
      c_i = ComplexUnivariate(this%paired_id,this%id,0)
      c_j = this
    endif
  endif
  
  if (present(that_pair)) then
    if (that%id<that%paired_id) then
      r_i = that
      r_j = that_pair
    else
      r_i = that_pair
      r_j = that
    endif
  else
    if (that%id<that%paired_id) then
      r_i = that
      r_j = RealUnivariate(that%paired_id,that%id,0)
    else
      r_i = RealUnivariate(that%paired_id,that%id,0)
      r_j = that
    endif
  endif
  
  ! Check that inputs are paired up.
  if (c_i%paired_id/=c_j%id .or. c_i%id/=c_j%paired_id) then
    call print_line(CODE_ERROR//': Unpaired modes passed to function.')
    call err()
  elseif (r_i%paired_id/=r_j%id .or. r_i%id/=r_j%paired_id) then
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
