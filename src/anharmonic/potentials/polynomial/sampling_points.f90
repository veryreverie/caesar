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
function generate_sampling_points(basis_functions,potential_expansion_order, &
   & maximum_weighted_displacement,frequency_of_max_displacement,real_modes) &
   & result(output)
  implicit none
  
  type(BasisFunctions), intent(in) :: basis_functions
  integer,              intent(in) :: potential_expansion_order
  real(dp),             intent(in) :: maximum_weighted_displacement
  real(dp),             intent(in) :: frequency_of_max_displacement
  type(RealMode),       intent(in) :: real_modes(:)
  type(SamplingPoints)             :: output
  
  type(RealMonomial), allocatable :: unique_bases(:)
  type(RealMonomial), allocatable :: matching_bases(:)
  
  type(RealModeDisplacement), allocatable :: points(:)
  
  integer :: i
  
  ! Sampling points are generated in groups, where each group contains all
  !    basis functions with the same modes with non-zero power.
  ! e.g. u1^2, u1^3 and u1^4 are all in one group,
  !    and u1^1*u3^2, u1^4&u3^1 and u1^5*u3^7 are all in one group.
  
  ! Identify one basis function in each group.
  unique_bases = basis_functions%unique_terms(         &
     & set(basis_functions%unique_terms,compare_modes) )
  
  ! Loop over the groups.
  points = [RealModeDisplacement::]
  do i=1,size(unique_bases)
    ! Identify the basis functions in the group.
    ! These are the basis functions equivalent to the representative function.
    matching_bases = basis_functions%unique_terms(                          &
       & filter(basis_functions%unique_terms,compare_modes,unique_bases(i)) )
    
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
    
    integer :: i
    
    select type(this); type is(RealMonomial)
      select type(that); type is(RealMonomial)
        this_ids = [integer::]
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
        
        that_ids = [integer::]
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
  
  output = [RealModeDisplacement::]
  
  do i=1,size(monomials)
    ids = [integer::]
    powers = [integer::]
    do j=1,size(monomials(i))
      if (monomials(i)%power(j)>0) then
        ids = [ids, monomials%id(j)]
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
