! ======================================================================
! Sampling points for sampling a given set of basis functions.
! ======================================================================
module sampling_points_module
  use common_module
  
  use anharmonic_common_module
  
  use basis_function_module
  use basis_functions_module
  implicit none
  
  private
  
  public :: SamplingPoints
  public :: generate_sampling_points
  public :: size
  
  type, extends(Stringsable) :: SamplingPoints
    type(RealModeDisplacement), allocatable :: points(:)
  contains
    procedure, public :: read  => read_SamplingPoints
    procedure, public :: write => write_SamplingPoints
  end type
  
  interface SamplingPoints
    module procedure new_SamplingPoints
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
function generate_sampling_points(basis_functions,potential_expansion_order, &
   & maximum_weighted_displacement,frequency_of_max_displacement,real_modes) &
   & result(output)
  implicit none
  
  type(BasisFunction), intent(in) :: basis_functions(:)
  integer,             intent(in) :: potential_expansion_order
  real(dp),            intent(in) :: maximum_weighted_displacement
  real(dp),            intent(in) :: frequency_of_max_displacement
  type(RealMode),      intent(in) :: real_modes(:)
  type(SamplingPoints)            :: output
  
  type(ModeCoupling), allocatable :: couplings(:)
  type(ModeCoupling), allocatable :: unique_couplings(:)
  
  type(RealMonomial),  allocatable :: unique_terms(:)
  
  type(RealModeDisplacement), allocatable :: points(:)
  
  integer, allocatable :: matching_couplings(:)
  
  integer :: i
  
  ! Construct the mode coupling corresponding to the unique term in each
  !    basis function.
  ! e.g. (u7)^4*(u9)^3*(u11)^1 => [7,9,11].
  couplings = ModeCoupling(basis_functions%unique_term)
  
  ! De-duplicate the couplings.
  unique_couplings = couplings(set(couplings,compare_ModeCoupling))
  
  points = [RealModeDisplacement::]
  do i=1,size(unique_couplings)
    ! Gather together all unique terms with couplings which are the same as
    !    unique_couplings(i).
    matching_couplings = filter(couplings==unique_couplings(i))
    unique_terms = basis_functions(matching_couplings)%unique_term
    
    ! Construct an array of sampling points for the unique terms,
    !    and append this to the output array.
    points = [                                                           &
       & points,                                                         &
       & generate_sampling_points_helper( unique_terms,                  &
       &                                  potential_expansion_order,     &
       &                                  maximum_weighted_displacement, &
       &                                  frequency_of_max_displacement, &
       &                                  real_modes)                    &
       & ]
  enddo
  
  output = SamplingPoints(points)
contains
  ! Lambda for comparing ModeCoupling.
  function compare_ModeCoupling(this,that) result(output)
    implicit none
    
    class(*), intent(in) :: this
    class(*), intent(in) :: that
    logical              :: output
    
    select type(this); type is(ModeCoupling)
      select type(that); type is(ModeCoupling)
        output = this==that
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
  
  integer                                   :: sum_powers
  type(RealSingleDisplacement), allocatable :: vectors(:)
  type(RealMode)                            :: mode
  real(dp)                                  :: mode_power
  real(dp)                                  :: mode_frequency
  real(dp)                                  :: magnitude
  
  integer :: i,j,ialloc
  
  output = [RealModeDisplacement::]
  
  do i=1,size(output)
    sum_powers = sum(monomials(i)%modes%power)
    allocate(vectors(size(monomials(i))), stat=ialloc); call err(ialloc)
    do j=1,size(monomials(i))
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
      mode = real_modes(first(real_modes%id==monomials(i)%modes(j)%id))
      mode_power = monomials(i)%modes(j)%power
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
    deallocate(vectors, stat=ialloc); call err(ialloc)
  enddo
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
  end select
end subroutine

function write_SamplingPoints(this) result(output)
  implicit none
  
  class(SamplingPoints), intent(in) :: this
  type(String), allocatable         :: output(:)
  
  select type(this); type is(SamplingPoints)
    output = str(this%points,separating_line='')
  end select
end function

impure elemental function new_SamplingPoints_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(SamplingPoints)          :: this
  
  this = input
end function
end module
