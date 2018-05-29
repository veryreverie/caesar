! ======================================================================
! A wrapped polymorphic pointer to a potential.
! ======================================================================
! Wraps all of PotentialData's methods,
!    calling them on the pointed-to potential.
! See example module below for how to use this type.
module potential_pointer_module
  use common_module
  
  use coupled_subspaces_module
  use potential_module
  implicit none
  
  private
  
  public :: PotentialPointer
  public :: assignment(=)
  
  type, extends(PotentialData) :: PotentialPointer
    class(PotentialData), allocatable :: potential
  contains
    procedure, public :: generate_sampling_points => &
       & generate_sampling_points_PotentialPointer
  end type
  
  interface assignment(=)
    module procedure assign_PotentialPointer_PotentialData
  end interface
contains

! Assign a PotentialData from any type which extends PotentialData.
subroutine assign_PotentialPointer_PotentialData(output,input)
  implicit none
  
  type(PotentialPointer), intent(out) :: output
  class(PotentialData),   intent(in)  :: input
  
  integer :: ialloc
  
  select type(input); class is(PotentialPointer)
    allocate( output%potential, source=input%potential, &
            & stat=ialloc); call err(ialloc)
  class default
    allocate( output%potential, source=input, &
            & stat=ialloc); call err(ialloc)
  end select
end subroutine

! Wrappers for all of PotentialData's methods.
subroutine generate_sampling_points_PotentialPointer(this,coupled_subspaces)
  implicit none
  
  class(PotentialPointer), intent(inout) :: this
  type(CoupledSubspaces),  intent(in)    :: coupled_subspaces(:)
  
  if (.not. allocated(this%potential)) then
    call print_line(CODE_ERROR//': Trying to use a PotentialPointer before &
       &it has been allocated.')
    call err()
  endif
  
  call this%potential%generate_sampling_points( coupled_subspaces)
end subroutine
end module

! ======================================================================
! An example module, showing how to use PotentialData and PotentialPointer.
! ======================================================================
module potential_example_module
  use common_module
  
  use coupled_subspaces_module
  use potential_module
  use potential_pointer_module
  implicit none
  
  private
  
  public :: potential_example_subroutine
  
  type, extends(PotentialData) :: PotentialDataExample
    type(String) :: example_contents
  contains
    procedure, public :: generate_sampling_points => &
       & generate_sampling_points_PotentialDataExample
  end type
  
  interface PotentialDataExample
    module procedure new_PotentialDataExample
  end interface
contains

! Constructor for example class.
! This is where any PotentialDataExample-specific data is input.
function new_PotentialDataExample(example_contents) result(this)
  implicit none
  
  type(String), intent(in)   :: example_contents
  type(PotentialDataExample) :: this
  
  this%example_contents = example_contents
end function

! Overloads of PotentialData's methods.
subroutine generate_sampling_points_PotentialDataExample(this, &
   & coupled_subspaces)
  implicit none
  
  class(PotentialDataExample), intent(inout) :: this
  type(CoupledSubspaces),      intent(in)    :: coupled_subspaces(:)
  
  ! Code to generate sampling points goes here.
  call print_line('PotentialDataExample: generating sampling points.')
end subroutine

! The class in use.
subroutine potential_example_subroutine(working_directory,coupled_subspaces)
  implicit none
  
  type(String),           intent(in) :: working_directory
  type(CoupledSubspaces), intent(in) :: coupled_subspaces(:)
  
  type(String) :: example_contents
  
  type(PotentialPointer) :: potential
  
  ! Set the pointer to point to a PotentialDataExample type.
  ! This is any PotentialDataExample-specific data is input,
  !    in this case the variable example_contents.
  example_contents = 'example'
  potential = PotentialDataExample(example_contents)
  
  ! Now PotentialData's methods can be called.
  ! They will all be forwareded to the PotentialDataExample instance.
  call potential%generate_sampling_points(coupled_subspaces)
end subroutine
end module
