! ======================================================================
! Sampling points for sampling a given set of basis functions.
! ======================================================================
module sampling_points_module
  use common_module
  
  use anharmonic_common_module
  
  use basis_function_module
  implicit none
  
  private
  
  public :: SamplingPoints
  public :: generate_sampling_points
  public :: size
  public :: construct_sample_matrix
  public :: make_sample_vector
  
  type, extends(Stringsable) :: SamplingPoints
    type(RealModeDisplacement), allocatable :: points(:)
  contains
    procedure, public :: read  => read_SamplingPoints
    procedure, public :: write => write_SamplingPoints
  end type
  
  interface SamplingPoints
    module procedure new_SamplingPoints
    module procedure new_SamplingPoints_Strings
    module procedure new_SamplingPoints_StringArray
  end interface
  
  interface size
    module procedure size_SamplingPoints
  end interface
contains

! Constructor and size() function.
function new_SamplingPoints(points) result(this)
  implicit none
  
  type(RealModeDisplacement), intent(in) :: points(:)
  type(SamplingPoints)                   :: this
  
  this%points = points
end function

function size_SamplingPoints(this) result(output)
  implicit none
  
  type(SamplingPoints), intent(in) :: this
  integer                          :: output
  
  output = size(this%points)
end function

! Generates a set of sampling points for sampling a set of basis functions.
function generate_sampling_points1(basis_functions,potential_expansion_order, &
   & maximum_weighted_displacement,frequency_of_max_displacement,real_modes, &
   & energy_to_force_ratio) result(output)
  implicit none
  
  type(BasisFunctionsAndUniqueTerms), intent(in) :: basis_functions(:)
  integer,                            intent(in) :: potential_expansion_order
  real(dp),                           intent(in) :: maximum_weighted_displacement
  real(dp),                           intent(in) :: frequency_of_max_displacement
  type(RealMode),                     intent(in) :: real_modes(:)
  real(dp),                           intent(in) :: energy_to_force_ratio
  type(SamplingPoints)                           :: output
  
  type(RealMonomial), allocatable :: unique_terms(:)
  
  type(RealMonomial), allocatable :: unique_bases(:)
  type(RealMonomial), allocatable :: matching_bases(:)
  
  type(RealModeDisplacement), allocatable :: points(:)
  
  integer :: i,ialloc
  
  unique_terms = [( basis_functions(i)%unique_terms, &
                  & i=1,                             &
                  & size(basis_functions)            )]
  
  ! Sampling points are generated in groups, where each group contains all
  !    basis functions with the same modes with non-zero power.
  ! e.g. u1^2, u1^3 and u1^4 are all in one group,
  !    and u1^1*u3^2, u1^4&u3^1 and u1^5*u3^7 are all in one group.
  
  ! Identify one basis function in each group.
  unique_bases = unique_terms(set(unique_terms, compare_modes))
  
  ! Loop over the groups.
  allocate(points(0), stat=ialloc); call err(ialloc)
  do i=1,size(unique_bases)
    ! Identify the basis functions in the group.
    ! These are the basis functions equivalent to the representative function.
    matching_bases = unique_terms(filter( unique_terms,   &
                                        & compare_modes,  &
                                        & unique_bases(i) ))
    
    ! Generate the sampling points for the group of basis functions.
    points = [                                                             &
       & points,                                                           &
       & generate_sampling_points_helper( matching_bases,                  &
       &                                  potential_expansion_order,       &
       &                                  maximum_weighted_displacement,   &
       &                                  frequency_of_max_displacement,   &
       &                                  real_modes                     ) ]
  enddo
  
  output = SamplingPoints(points)
contains
  ! Lambda for comparing if two RealMonomials contain the same set of modes
  !    with non-zero power.
  function compare_modes(this,that) result(output)
    implicit none
    
    class(*), intent(in) :: this
    class(*), intent(in) :: that
    logical              :: output
    
    integer, allocatable :: this_ids(:)
    integer, allocatable :: that_ids(:)
    
    integer :: i,ialloc
    
    select type(this); type is(RealMonomial)
      select type(that); type is(RealMonomial)
        allocate(this_ids(0), stat=ialloc); call err(ialloc)
        do i=1,size(this)
          if (this%power(i)>0) then
            this_ids = [this_ids, this%id(i)]
          endif
          
          if (this%paired_id(i)/=this%id(i)) then
            if (this%paired_power(i)>0) then
              this_ids = [this_ids, this%paired_id(i)]
            endif
          endif
        enddo
        
        allocate(that_ids(0), stat=ialloc); call err(ialloc)
        do i=1,size(that)
          if (that%power(i)>0) then
            that_ids = [that_ids, that%id(i)]
          endif
          
          if (that%paired_id(i)/=that%id(i)) then
            if (that%paired_power(i)>0) then
              that_ids = [that_ids, that%paired_id(i)]
            endif
          endif
        enddo
        
        this_ids = this_ids(set(this_ids))
        that_ids = that_ids(set(that_ids))
        
        if (size(this_ids)/=size(that_ids)) then
          output = .false.
        else
          this_ids = this_ids(sort(this_ids))
          that_ids = that_ids(sort(that_ids))
          output = all(this_ids==that_ids)
        endif
      end select
    end select
  end function
end function

! Helper function for generate_sampling_points.
function generate_sampling_points_helper(monomials,potential_expansion_order, &
   & maximum_weighted_displacement,frequency_of_max_displacement,real_modes)  &
   & result(output)
  implicit none
  
  type(RealMonomial), intent(in)          :: monomials(:)
  integer,            intent(in)          :: potential_expansion_order
  real(dp),           intent(in)          :: maximum_weighted_displacement
  real(dp),           intent(in)          :: frequency_of_max_displacement
  type(RealMode),     intent(in)          :: real_modes(:)
  type(RealModeDisplacement), allocatable :: output(:)
  
  integer, allocatable :: ids(:)
  integer, allocatable :: powers(:)
  
  integer                                   :: sum_powers
  type(RealSingleDisplacement), allocatable :: vectors(:)
  type(RealMode)                            :: mode
  real(dp)                                  :: mode_power
  real(dp)                                  :: mode_frequency
  real(dp)                                  :: magnitude
  
  integer :: i,j,ialloc
  
  allocate(output(0), stat=ialloc); call err(ialloc)
  
  do i=1,size(monomials)
    allocate( ids(0),    &
            & powers(0), &
            & stat=ialloc); call err(ialloc)
    do j=1,size(monomials(i))
      if (monomials(i)%power(j)>0) then
        ids = [ids, monomials(i)%id(j)]
        powers = [powers, monomials(i)%power(j)]
      endif
      
      if ( monomials(i)%paired_id(j)/=monomials(i)%id(j) .and. &
         & monomials(i)%paired_power(j)>0                      ) then
        ids = [ids, monomials(i)%paired_id(j)]
        powers = [powers, monomials(i)%paired_power(j)]
      endif
    enddo
    
    sum_powers = sum(powers)
    
    allocate(vectors(size(ids)), stat=ialloc); call err(ialloc)
    do j=1,size(ids)
      ! Calculate the displacement along mode j in monomial i.
      ! This is equal to d_j = d_max *     ( p_i   / p_max          )
      !                              * sqrt( w_min / max(w_j,w_min) )
      !                              * sqrt( p_j   / p_i)
      !    - d_max is the maximum mass-weighted displacement.
      !    - p_j   is the power of mode j in monomial i.
      !    - p_i   is the sum of the powers of the modes in monomial i.
      !    - p_max is the maximum sum of powers in any monomial.
      !    - w_j   is the frequency of mode j in monomial i.
      !    - w_min is the frequency cutoff for displacement scaling.
      !
      ! As such the displacement corresponding to a given monomial
      !    (the L2 sum across d_j for that monomial), is at most d_max.
      mode = real_modes(first(real_modes%id==ids(j)))
      mode_power = powers(j)
      mode_frequency = max(mode%frequency, frequency_of_max_displacement)
      magnitude = maximum_weighted_displacement       &
              & * sqrt( mode_power                    &
              &       * sum_powers                    &
              &       * frequency_of_max_displacement &
              &       / mode_frequency )              &
              & / potential_expansion_order
      
      vectors(j) = RealSingleDisplacement(mode, magnitude)
    enddo
    
    output = [output, RealModeDisplacement(vectors)]
    output = [output, -RealModeDisplacement(vectors)]
    deallocate(ids, powers, vectors, stat=ialloc); call err(ialloc)
  enddo
end function

! Construct the matrix necessary for L2 fitting of basis function coefficients.
! The elements are the energies and forces at a given sampling point due to
!    a given basis function.
function construct_sample_matrix(basis_functions,sampling_points, &
   & modes,energy_force_ratio) result(output)
  implicit none
  
  type(BasisFunction),        intent(in) :: basis_functions(:)
  type(RealModeDisplacement), intent(in) :: sampling_points(:)
  type(RealMode),             intent(in) :: modes(:)
  real(dp),                   intent(in) :: energy_force_ratio
  real(dp), allocatable                  :: output(:,:)
  
  real(dp)            :: energy
  type(RealModeForce) :: forces
  
  integer :: dims
  
  integer :: i,j,ialloc
  
  dims = 1+size(modes)
  
  allocate( output( size(sampling_points)*dims,    &
          &         size(basis_functions)       ), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(basis_functions)
    do j=1,size(sampling_points)
      energy = basis_functions(i)%energy(sampling_points(j))
      forces = basis_functions(i)%force(sampling_points(j))
      
      output((j-1)*dims+1:j*dims, i) = make_sample_vector( energy,            &
                                                         & forces,            &
                                                         & modes,             &
                                                         & energy_force_ratio )
    enddo
  enddo
end function

! Convert an energy and force into a single vector which can be inserted
!    into a matrix for passing to LAPACK.
function make_sample_vector(energy,force,modes,energy_force_ratio) &
   & result(output)
  implicit none
  
  real(dp),            intent(in) :: energy
  type(RealModeForce), intent(in) :: force
  type(RealMode),      intent(in) :: modes(:)
  real(dp),            intent(in) :: energy_force_ratio
  real(dp), allocatable           :: output(:)
  
  integer :: i,ialloc
  
  allocate(output(1+size(modes)), stat=ialloc); call err(ialloc)
  output(1) = energy / energy_force_ratio
  do i=1,size(modes)
    output(1+i) = force%force(modes(i))
  enddo
end function

! Generates a set of sampling points for sampling a set of basis functions.
function generate_sampling_points(basis_functions_,potential_expansion_order, &
   & maximum_weighted_displacement,frequency_of_max_displacement,real_modes, &
   & energy_to_force_ratio) result(output)
  implicit none
  
  type(BasisFunctionsAndUniqueTerms), intent(in) :: basis_functions_(:)
  integer,                            intent(in) :: potential_expansion_order
  real(dp),                           intent(in) :: maximum_weighted_displacement
  real(dp),                           intent(in) :: frequency_of_max_displacement
  type(RealMode),                     intent(in) :: real_modes(:)
  real(dp),                           intent(in) :: energy_to_force_ratio
  type(SamplingPoints)                           :: output
  
  type(BasisFunction), allocatable :: basis_functions(:)
  
  type(RealModeDisplacement), allocatable :: basis_points(:)
  logical,                    allocatable :: basis_points_used(:)
  real(dp),                   allocatable :: determinants(:)
  
  type(RealModeDisplacement), allocatable :: sampling_points(:)
  
  integer :: points_per_basis_function
  
  integer :: i,j,k,l,m,ialloc
  
  points_per_basis_function = 2
  
  basis_functions = [( basis_functions_(i)%basis_functions, &
                     & i=1,                                 &
                     & size(basis_functions_)               )]
  
  allocate( sampling_points(points_per_basis_function*size(basis_functions)), &
          & stat=ialloc); call err(ialloc)
  l = 0
  do i=1,size(basis_functions)
    ! Generate the possible sampling points for this basis function.
    basis_points = generate_basis_points( basis_functions(i),            &
                                        & potential_expansion_order,     &
                                        & maximum_weighted_displacement, &
                                        & frequency_of_max_displacement, &
                                        & real_modes                     )
    basis_points_used = [(.false., j=1, size(basis_points))]
    determinants = [(0.0_dp, j=1, size(basis_points))]
    
    ! Find the sampling point which most resolves the first i basis functions.
    ! Repeat until points_per_basis_function sampling points have been found.
    do j=1,points_per_basis_function
      l = l+1
      do k=1,size(basis_points)
        if (.not. basis_points_used(k)) then
          sampling_points(l) = basis_points(k)
          determinants(k) = calculate_average_eigenvalue( &
                                  & basis_functions(:i),  &
                                  & sampling_points(:l),  &
                                  & real_modes,           &
                                  & energy_to_force_ratio )
        endif
      enddo
      m = maxloc(abs(determinants), 1, mask=.not.basis_points_used)
      sampling_points(l) = basis_points(m)
      basis_points_used(m) = .true.
      
      if (i==size(basis_functions) .and. j==points_per_basis_function) then
        call print_line( 'Average |eigenvalue| of fitting matrix: '// &
                       & determinants(m)*energy_to_force_ratio//' (Ha).')
      endif
    enddo
  enddo
  
  output = SamplingPoints(sampling_points)
end function

function generate_basis_points(basis_function,potential_expansion_order, &
   & maximum_weighted_displacement,frequency_of_max_displacement,modes)  &
   & result(output)
  implicit none
  
  type(BasisFunction), intent(in)         :: basis_function
  integer,             intent(in)         :: potential_expansion_order
  real(dp),            intent(in)         :: maximum_weighted_displacement
  real(dp),            intent(in)         :: frequency_of_max_displacement
  type(RealMode),      intent(in)         :: modes(:)
  type(RealModeDisplacement), allocatable :: output(:)
  
  type(ComplexMonomial), allocatable :: monomials(:)
  
  type(RealSingleDisplacement), allocatable :: displacement(:)
  
  type(RealMode) :: mode
  real(dp)       :: frequency
  real(dp)       :: magnitude
  
  integer :: i,j,k,ialloc
  
  monomials = basis_function%terms()
  allocate(output(2*size(monomials)), stat=ialloc); call err(ialloc)
  do i=1,size(monomials)
    allocate(displacement(0), stat=ialloc); call err(ialloc)
    do j=1,size(monomials(i))
      if (monomials(i)%power(j)>0) then
        mode = modes(first(modes%id==monomials(i)%id(j)))
        frequency = max(mode%frequency, frequency_of_max_displacement)
        magnitude = maximum_weighted_displacement         &
                & * sqrt( monomials(i)%power(j)           &
                &       * frequency_of_max_displacement   &
                &       / frequency                     ) &
                & / potential_expansion_order
        displacement = [displacement, RealSingleDisplacement(mode, magnitude)]
      endif
      if (monomials(i)%id(j)/=monomials(i)%paired_id(j)) then
        if (monomials(i)%paired_power(j)>0) then
          mode = modes(first(modes%id==monomials(i)%paired_id(j)))
          frequency = max(mode%frequency, frequency_of_max_displacement)
          magnitude = maximum_weighted_displacement         &
                  & * sqrt( monomials(i)%paired_power(j)    &
                  &       * frequency_of_max_displacement   &
                  &       / frequency                     ) &
                  & / potential_expansion_order
          displacement = [ displacement,                           &
                         & RealSingleDisplacement(mode, magnitude) ]
        endif
      endif
    enddo
    output(2*i-1:2*i) = [  RealModeDisplacement(displacement), &
                        & -RealModeDisplacement(displacement)  ]
    deallocate(displacement, stat=ialloc); call err(ialloc)
  enddo
end function

! Calculates the geometric average of the absolute eigenvalues of A^T . A,
!    the matrix which must be inverted to calculate
!    basis function coefficients.
function calculate_average_eigenvalue(basis_functions,sampling_points,modes, &
   & energy_force_ratio) result(output)
  implicit none
  
  type(BasisFunction),        intent(in) :: basis_functions(:)
  type(RealModeDisplacement), intent(in) :: sampling_points(:)
  type(RealMode),             intent(in) :: modes(:)
  real(dp),                   intent(in) :: energy_force_ratio
  real(dp)                               :: output
  
  type(RealMatrix)          :: matrix
  type(RealQRDecomposition) :: qr
  real(dp), allocatable     :: evals(:)
  
  integer :: i
  
  matrix = mat( construct_sample_matrix( basis_functions,     &
              &                          sampling_points,     &
              &                          modes,               &
              &                          energy_force_ratio ) )
  
  qr = qr_decomposition(transpose(matrix)*matrix)
  
  ! The eigenvalues are the diagonal entries of the R matrix
  !    in a QR decomposition.
  evals = [(qr%r(i,i), i=1, size(qr%r,1))]
  
  evals = abs(evals)
  
  if (any(evals<1e-200_dp)) then
    output = 0.0_dp
    return
  endif
  
  output = product(evals**(1.0_dp/size(qr%r,1)))
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_SamplingPoints(this,input)
  implicit none
  
  class(SamplingPoints), intent(out) :: this
  type(String),          intent(in)  :: input(:)
  
  select type(this); type is(SamplingPoints)
    this = SamplingPoints(RealModeDisplacement(split_into_sections(input)))
  class default
    call err()
  end select
end subroutine

function write_SamplingPoints(this) result(output)
  implicit none
  
  class(SamplingPoints), intent(in) :: this
  type(String), allocatable         :: output(:)
  
  select type(this); type is(SamplingPoints)
    output = str(this%points,separating_line='')
  class default
    call err()
  end select
end function

function new_SamplingPoints_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(SamplingPoints)     :: this
  
  call this%read(input)
end function

impure elemental function new_SamplingPoints_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(SamplingPoints)          :: this
  
  this = SamplingPoints(str(input))
end function
end module
