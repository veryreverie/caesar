! ======================================================================
! Calculates and stores 1-D mode overlaps for interpolating polynomials
!    from a coarse q-point grid to a fine q-point grid.
! ======================================================================
! The absolute overlap for atom j and modes u1 and u2 is just <u1|j><j|u2>,
!    i.e. the inner product of the projections of u1 and u2 onto atom j.
! The relative overlap for atom j and modes u1 and u2 is
!    sum_k <u1|k>e^{2 pi i (q1-q2).R_{jk}}<k|u2>, where R_{jk} is the minimum
!    image R-vector between atoms j and k.
! The overlap between two products of modes {u1} and {u2} is
!    sum_atom sum_j [ absolute(u1j,u2j,atom) 
!                   * product_{k/=j} [ relative(u1k,u2k,atom) ] ]
! N.B. the product overlap must be averaged over permutations of {u1} against
!    {u2}.
module polynomial_interpolator_module
  use common_module
  use anharmonic_common_module
  use permutation_module
  implicit none
  
  private
  
  public :: PolynomialInterpolator
  public :: construct_fine_monomials
  
  type, extends(NoDefaultConstructor) :: PolynomialInterpolator
    complex(dp), allocatable, private :: absolute_(:,:,:)
    complex(dp), allocatable, private :: relative_(:,:,:)
    integer,     allocatable, private :: fine_map_(:)
    integer,     allocatable, private :: coarse_map_(:)
  contains
    procedure, private :: absolute_overlap
    procedure, private :: relative_overlap
    procedure, private :: product_overlap
    
    generic,   public  :: overlap =>                               &
                        & overlap_ComplexMonomial_ComplexMonomial, &
                        & overlap_ComplexMonomial_ComplexPolynomial
    procedure, private :: overlap_ComplexMonomial_ComplexMonomial
    procedure, private :: overlap_ComplexMonomial_ComplexPolynomial
  end type
  
  public :: new_PolynomialInterpolator
  interface PolynomialInterpolator
    module procedure new_PolynomialInterpolator
  end interface
contains

! Constructor.
function new_PolynomialInterpolator(min_images,anharmonic_data, &
   & interpolated_anharmonic_data) result(this) 
  implicit none
  
  type(MinImages),      intent(in) :: min_images(:,:)
  type(AnharmonicData), intent(in) :: anharmonic_data
  type(AnharmonicData), intent(in) :: interpolated_anharmonic_data
  type(PolynomialInterpolator)     :: this
  
  type(ComplexMode), allocatable :: fine_modes(:)
  type(RealVector),  allocatable :: fine_qpoints(:)
  type(ComplexMode), allocatable :: coarse_modes(:)
  type(Realvector),  allocatable :: coarse_qpoints(:)
  
  type(AtomData) :: atom_i
  type(AtomData) :: atom_l
  
  type(RealVector) :: q
  
  integer :: i,j,k,l,ialloc
  
  type(IntVector), allocatable :: rvectors(:)
  
  coarse_modes = anharmonic_data%complex_modes
  fine_modes = interpolated_anharmonic_data%complex_modes
  fine_qpoints = dblevec(interpolated_anharmonic_data%qpoints%qpoint)
  
  coarse_qpoints = [(                                                       &
     & dblevec(anharmonic_data%qpoints(first(                               &
     &    anharmonic_data%qpoints%id==coarse_modes(i)%qpoint_id ))%qpoint), &
     & i=1,                                                                 &
     & size(coarse_modes)                                                   )]
  fine_qpoints = [(                                                         &
     & dblevec(interpolated_anharmonic_data%qpoints(first(                  &
     &    interpolated_anharmonic_data%qpoints%id==fine_modes(i)%qpoint_id  &
     &                                                     ))%qpoint),      &
     & i=1,                                                                 &
     & size(fine_modes)                                                     )]
  
  this%fine_map_ = [(0,i=1,maxval(fine_modes%id))]
  do i=1,size(fine_modes)
    this%fine_map_(fine_modes(i)%id) = i
  enddo
  
  this%coarse_map_ = [(0,i=1,maxval(coarse_modes%id))]
  do i=1,size(coarse_modes)
    this%coarse_map_(coarse_modes(i)%id) = i
  enddo
  
  allocate( this%absolute_( size(fine_modes),                       &
          &                 size(coarse_modes),                     &
          &                 size(anharmonic_data%structure%atoms)), &
          & this%relative_( size(fine_modes),                       &
          &                 size(coarse_modes),                     &
          &                 size(anharmonic_data%structure%atoms)), &
          & stat=ialloc); call err(ialloc)
  this%relative_ = 0
  do i=1,size(anharmonic_data%structure%atoms)
    atom_i = anharmonic_data%structure%atoms(i)
    do j=1,size(coarse_modes)
      do k=1,size(fine_modes)
        this%absolute_(k, j, i) = coarse_modes(j)%unit_vector(atom_i%id()) &
                              & * conjg(fine_modes(k)%unit_vector(atom_i%id()))
        
        q = coarse_qpoints(j)-fine_qpoints(k)
        do l=1,size(anharmonic_data%anharmonic_supercell%atoms)
          atom_l = anharmonic_data%anharmonic_supercell%atoms(l)
          rvectors = min_images(i,l)%image_rvectors
          this%relative_(k, j, i) =                                 &
             &   this%relative_(k, j, i)                            &
             & + coarse_modes(j)%unit_vector(atom_l%prim_id())      &
             & * conjg(fine_modes(k)%unit_vector(atom_l%prim_id())) &
             & * sum(exp_2pii(q*rvectors))                          &
             & / size(rvectors)
        enddo
      enddo
    enddo
  enddo
  
  this%absolute_ =                                                   &
     & ( this%absolute_                                              &
     & * interpolated_anharmonic_data%anharmonic_supercell%sc_size ) &
     & / anharmonic_data%anharmonic_supercell%sc_size
  
  this%relative_ = this%relative_ &
               & / anharmonic_data%anharmonic_supercell%sc_size
end function

! Calculates <u1|j><j|u2>.
impure elemental function absolute_overlap(this,fine_mode,coarse_mode,atom) &
   & result(output)
  implicit none
  
  class(PolynomialInterpolator), intent(in) :: this
  integer,                       intent(in) :: fine_mode
  integer,                       intent(in) :: coarse_mode
  integer,                       intent(in) :: atom
  complex(dp)                               :: output
  
  output = this%absolute_( this%fine_map_(fine_mode),     &
                         & this%coarse_map_(coarse_mode), &
                         & atom                           )
end function

! Calculates sum_k <u1|k>e^{2 pi i (q1-q2).R_{jk}}<k|u2>.
impure elemental function relative_overlap(this,fine_mode,coarse_mode,atom) &
   & result(output)
  implicit none
  
  class(PolynomialInterpolator), intent(in) :: this
  integer,                       intent(in) :: fine_mode
  integer,                       intent(in) :: coarse_mode
  integer,                       intent(in) :: atom
  complex(dp)                               :: output
  
  output = this%relative_( this%fine_map_(fine_mode),     &
                         & this%coarse_map_(coarse_mode), &
                         & atom                           )
end function

! Calculates <{u1}|{u2}> for a single permutation of {u1} and {u2}.
function product_overlap(this,fine_modes,coarse_modes) result(output)
  implicit none
  
  class(PolynomialInterpolator), intent(in) :: this
  integer,                       intent(in) :: fine_modes(:)
  integer,                       intent(in) :: coarse_modes(:)
  complex(dp)                               :: output
  
  complex(dp), allocatable :: absolute_overlaps(:)
  complex(dp), allocatable :: relative_overlaps(:)
  
  integer :: i,j,k
  
  output = 0
  ! Sum across atoms.
  do i=1,size(this%relative_,3)
    absolute_overlaps = this%absolute_overlap(fine_modes, coarse_modes, i)
    relative_overlaps = this%relative_overlap(fine_modes, coarse_modes, i)
    
    ! The term adds
    !    sum_j[absolute(j)*product_{k/=j}[relative(k)]]
    !    / sum_j[1].
    ! If all relative(k) /= 0 for all k, this is equivalent to
    !    product_j[relative(j)]*sum_j[absolute(j)/relative(j)]
    !    / sum_j[1].
    ! If relative(k)==0 for k=j only, this is
    !    absolute(j)*product_{k/=j}[relative(j)] / sum_j[1]
    ! If relative(k)==0 for more than one value of k, this is zero.
    
    j = first(abs(relative_overlaps)<1e-200_dp, default=0)
    if (j==0) then
      ! All relative(k) /= 0.
      output = output                     &
           & + product(relative_overlaps) &
           & * sum(absolute_overlaps/relative_overlaps)
    else
      k = first(abs(relative_overlaps(j+1:))<1e-200_dp, default=0)
      if (k==0) then
        ! relative(k)==0 only for k=j.
        output = output                           &
             & + absolute_overlaps(j)             &
             & * product(relative_overlaps(:j-1)) &
             & * product(relative_overlaps(j+1:))
      else
        ! relative(k)==0 for more than one value of k.
        continue
      endif
    endif
  enddo
  
  output = output / size(fine_modes)
end function

! Calculates the overlap between a monomial on the fine grid and a monomial
!    on the coarse grid.
impure elemental function overlap_ComplexMonomial_ComplexMonomial(this,fine, &
   & coarse) result(output)
  implicit none
  
  class(PolynomialInterpolator), intent(in) :: this
  type(ComplexMonomial),         intent(in) :: fine
  type(ComplexMonomial),         intent(in) :: coarse
  complex(dp)                               :: output
  
  integer, allocatable :: fine_ids(:)
  integer, allocatable :: coarse_ids(:)
  
  type(PermutationData) :: permutation
  
  integer :: no_permutations
  
  if (fine%total_power()/=coarse%total_power()) then
    output = 0
    return
  endif
  
  fine_ids = to_ids(fine)
  coarse_ids = to_ids(coarse)
  
  permutation = PermutationData(coarse_ids, fine_ids)
  
  no_permutations = 0
  output = 0
  do
    fine_ids = permutation%b()
    coarse_ids = permutation%a()
    
    output = output                                    &
         & + this%product_overlap(fine_ids,coarse_ids) &
         & * exp(permutation%log_no_permutations())
    
    call permutation%next_permutation()
    
    no_permutations = no_permutations + 1
    
    if (permutation%all_permutations_done()) then
      exit
    endif
  enddo
  
  output = output             &
       & * coarse%coefficient &
       & * exp(log_no_permutations(fine) - log_no_permutations(coarse))
end function

function to_ids(monomial) result(output)
  implicit none
  
  type(ComplexMonomial), intent(in) :: monomial
  integer, allocatable              :: output(:)
  
  integer :: i,j,ialloc
  
  allocate(output(monomial%total_power()), stat=ialloc); call err(ialloc)
  j = 0
  do i=1,size(monomial)
    output(j+1:j+monomial%power(i)) = monomial%id(i)
    j = j+monomial%power(i)
    if (monomial%id(i)/=monomial%paired_id(i)) then
      output(j+1:j+monomial%paired_power(i)) = monomial%paired_id(i)
      j = j+monomial%paired_power(i)
    endif
  enddo
end function

! Calculates the number of ways a given monomial can be permuted.
function log_no_permutations(monomial) result(output)
  implicit none
  
  type(ComplexMonomial), intent(in) :: monomial
  real(dp)                          :: output
  
  integer :: power
  integer :: total
  
  integer :: i
  
  output = 0
  
  total = 0
  do i=1,size(monomial)
    power = monomial%power(i)
    total = total + power
    output = output - log_factorial(power)
    if (monomial%id(i)/=monomial%paired_id(i)) then
      power = monomial%paired_power(i)
      total = total + power
      output = output - log_factorial(power)
    endif
  enddo
  
  output = output + log_factorial(total)
end function

! Calculates the overlap between an monomial on the fine grid and a polynomial
!    on the coarse grid.
impure elemental function overlap_ComplexMonomial_ComplexPolynomial( &
   & this,fine,coarse) result(output)
  implicit none
  
  class(PolynomialInterpolator), intent(in) :: this
  type(ComplexMonomial),         intent(in) :: fine
  type(ComplexPolynomial),       intent(in) :: coarse
  complex(dp)                               :: output
  
  output = sum(this%overlap(fine,coarse%terms))
end function

! Construct the monomials at a q-point on the fine grid.
function construct_fine_monomials(power,modes) result(output)
  implicit none
  
  integer,           intent(in)      :: power
  type(ComplexMode), intent(in)      :: modes(:)
  type(ComplexMonomial), allocatable :: output(:)
  
  type(ComplexUnivariate)            :: no_modes(0)
  type(ComplexMonomial), allocatable :: old(:)
  
  integer :: i,j,k,l
  
  ! Construct the monomial containing no modes.
  old = [ComplexMonomial( coefficient = cmplx(1.0_dp,0.0_dp,dp), &
                        & modes       = no_modes                 )]
  
  ! Loop over each mode in modes.
  ! For each mode, loop over all monomials containing the previous modes,
  !    and for each monomial append the new mode with all valid powers and
  !    paired_powers such that sum(monomial%powers()) and
  !    sum(monomial%paired_powers()) are both between 0 and power inclusive.
  do i=1,size(modes)-1
    output = [(                                                       &
       & (                                                            &
       &   ( ComplexMonomial(                                         &
       &        coefficient = cmplx(1.0_dp,0.0_dp,dp),                &
       &        modes       = [ old(l)%modes(),                       &
       &                        ComplexUnivariate(                    &
       &                           id           = modes(i)%id,        &
       &                           paired_id    = modes(i)%paired_id, &
       &                           power        = j,                  &
       &                           paired_power = k            )] ),  &
       &     j=0,                                                     &
       &     power-sum(old(l)%powers()) ),                            &
       &   k=0,                                                       &
       &   power-sum(old(l)%paired_powers()) ),                       &
       & l=1,                                                         &
       & size(old)                                                    )]
    old = output
  enddo
  
  ! Add powers and paired_powers of the final mode such that
  !    sum(monomial%powers()) = sum(monomial%paired_powers()) = power.
  output = [(                                                        &
     & ComplexMonomial(                                                 &
     &    coefficient = cmplx(1.0_dp,0.0_dp,dp),                        &
     &    modes       = [                                               &
     &       old(i)%modes(),                                            &
     &       ComplexUnivariate(                                         &
     &          id           = modes(size(modes))%id,                   &
     &          paired_id    = modes(size(modes))%paired_id,            &
     &          power        = power-sum(old(i)%powers()),              &
     &          paired_power = power-sum(old(i)%paired_powers()) ) ] ), &
     & i=1,                                                             &
     & size(old)                                                        )]
  
  call output%simplify()
end function
end module
