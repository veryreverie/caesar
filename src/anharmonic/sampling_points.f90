! ======================================================================
! Sampling points for sampling a given set of basis functions.
! ======================================================================
module sampling_points_module
  use common_module
  
  use coupled_modes_module
  use mode_monomial_module
  use basis_function_module
  use basis_functions_module
  implicit none
  
  private
  
  public :: SamplingPoints
  public :: generate_sampling_points
  
  type :: SamplingPoints
    type(RealModeDisplacement), allocatable :: points(:)
  contains
    procedure, public :: write_file => write_file_SamplingPoints
    procedure, public :: read_file  => read_file_SamplingPoints
  end type
  
  interface size
    module procedure size_SamplingPoints
  end interface
contains

function size_SamplingPoints(this) result(output)
  implicit none
  
  type(SamplingPoints), intent(in) :: this
  integer                          :: output
  
  output = size(this%points)
end function

function generate_sampling_points(basis_functions,potential_expansion_order, &
   & maximum_weighted_displacement) result(output)
  implicit none
  
  type(BasisFunction), intent(in) :: basis_functions(:)
  integer,             intent(in) :: potential_expansion_order
  real(dp),            intent(in) :: maximum_weighted_displacement
  type(SamplingPoints)            :: output
  
  type(CoupledModes), allocatable :: couplings(:)
  type(CoupledModes), allocatable :: unique_couplings(:)
  
  type(RealMonomial),  allocatable :: unique_terms(:)
  
  type(RealModeDisplacement), allocatable :: points(:)
  
  integer :: i,ialloc
  
  ! Construct the mode coupling corresponding to the unique term in each
  !    basis function.
  ! e.g. (u7)^4*(u9)^3*(u11)^1 => [7,9,11].
  couplings = construct_coupled_modes(basis_functions%unique_term)
  
  ! De-duplicate the couplings.
  unique_couplings = couplings(set(couplings,compare_CoupledModes))
  
  points = [RealModeDisplacement::]
  do i=1,size(unique_couplings)
    ! Gather together all unique terms with couplings which are the same as
    !    unique_couplings(i).
    unique_terms = &
       & basis_functions(filter(couplings==unique_couplings(i)))%unique_term
    
    ! Construct an array of sampling points for the unique terms,
    !    and append this to the output array.
    points = [                                                           &
       & points,                                                         &
       & generate_sampling_points_helper( unique_terms,                  &
       &                                  potential_expansion_order,     &
       &                                  maximum_weighted_displacement) &
       & ]
  enddo
  
  output = SamplingPoints(points)
contains
  ! Lambda for comparing CoupledModes.
  function compare_CoupledModes(this,that) result(output)
    implicit none
    
    class(*), intent(in) :: this
    class(*), intent(in) :: that
    logical              :: output
    
    select type(this); type is(CoupledModes)
      select type(that); type is(CoupledModes)
        output = this==that
      end select
    end select
  end function
end function

! Helper function for generate_sampling_points.
function generate_sampling_points_helper(monomials,potential_expansion_order, &
   & maximum_weighted_displacement) result(output)
  implicit none
  
  type(RealMonomial), intent(in)          :: monomials(:)
  integer,            intent(in)          :: potential_expansion_order
  real(dp),           intent(in)          :: maximum_weighted_displacement
  type(RealModeDisplacement), allocatable :: output(:)
  
  integer               :: power
  real(dp), allocatable :: fractions(:)
  
  integer :: i,j,ialloc
  
  allocate(output(size(monomials)), stat=ialloc); call err(ialloc)
  do i=1,size(output)
    power = sum(monomials(i)%modes%power)
    allocate( output(i)%displacements(size(monomials(i))), &
            & stat=ialloc); call err(ialloc)
    output(i)%displacements%id = monomials(i)%modes%id
    fractions = monomials(i)%modes%power / real(power,dp)
    output(i)%displacements%displacement = maximum_weighted_displacement &
                                       & * sqrt(fractions)
  enddo
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine write_file_SamplingPoints(this,filename)
  implicit none
  
  class(SamplingPoints), intent(in) :: this
  type(String),          intent(in) :: filename
  
  type(OFile) :: file
  
  integer :: i
  
  file = OFile(filename)
  do i=1,size(this)
    call file%print_line('Sampling point '//i)
    call file%print_lines(this%points(i))
    call file%print_line('')
  enddo
end subroutine

subroutine read_file_SamplingPoints(this,filename)
  implicit none
  
  class(SamplingPoints), intent(out) :: this
  type(String),          intent(in)  :: filename
  
  type(IFile)                    :: file
  type(StringArray), allocatable :: sections(:)
  
  integer :: i,ialloc
  
  file = IFile(filename)
  sections = split(file%lines())
  allocate(this%points(size(sections)), stat=ialloc); call err(ialloc)
  do i=1,size(sections)
    this%points(i) = RealModeDisplacement(sections(i)%strings(2:))
  enddo
end subroutine
end module
