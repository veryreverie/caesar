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
  
  integer :: i,j,ialloc
  
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
