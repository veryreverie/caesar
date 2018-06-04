! ======================================================================
! A wrapped polymorphic pointer to a potential.
! ======================================================================
! Wraps all of PotentialData's methods,
!    calling them on the pointed-to potential.
! See example module below for how to use this type.
module potential_pointer_module
  use common_module
  
  use anharmonic_data_module
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
    procedure, public :: generate_potential => &
       & generate_potential_PotentialPointer
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
subroutine generate_sampling_points_PotentialPointer(this,inputs, &
   & sampling_points_dir,logfile,write_lambda)
  implicit none
  
  class(PotentialPointer), intent(inout) :: this
  type(AnharmonicData),    intent(in)    :: inputs
  type(String),            intent(in)    :: sampling_points_dir
  type(OFile),             intent(inout) :: logfile
  procedure(WriteLambda)                 :: write_lambda
  
  if (.not. allocated(this%potential)) then
    call print_line(CODE_ERROR//': Trying to use a PotentialPointer before &
       &it has been allocated.')
    call err()
  endif
  
  call this%potential%generate_sampling_points( inputs,              &
                                              & sampling_points_dir, &
                                              & logfile,             &
                                              & write_lambda)
end subroutine

subroutine generate_potential_PotentialPointer(this,inputs, &
   & sampling_points_dir,logfile,read_lambda)
  implicit none
  
  class(PotentialPointer), intent(inout) :: this
  type(AnharmonicData),    intent(in)    :: inputs
  type(String),            intent(in)    :: sampling_points_dir
  type(OFile),             intent(inout) :: logfile
  procedure(ReadLambda)                  :: read_lambda
  
  if (.not. allocated(this%potential)) then
    call print_line(CODE_ERROR//': Trying to use a PotentialPointer before &
       &it has been allocated.')
    call err()
  endif
  
  call this%potential%generate_potential( inputs,              &
                                        & sampling_points_dir, &
                                        & logfile,             &
                                        & read_lambda)
end subroutine
end module

! ======================================================================
! An example module, showing how to use PotentialData and PotentialPointer.
! ======================================================================
module potential_example_module
  use common_module
  
  use anharmonic_data_module
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
    procedure, public :: generate_potential => &
       & generate_potential_PotentialDataExample
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
subroutine generate_sampling_points_PotentialDataExample(this,inputs, &
   & sampling_points_dir,logfile,write_lambda)
  implicit none
  
  class(PotentialDataExample), intent(inout) :: this
  type(AnharmonicData),        intent(in)    :: inputs
  type(String),                intent(in)    :: sampling_points_dir
  type(OFile),                 intent(inout) :: logfile
  procedure(WriteLambda)                     :: write_lambda
  
  call print_line('PotentialDataExample: generating sampling points.')
  call print_line('Example contents = '//this%example_contents)
  
  ! Code to generate sampling points goes here.
end subroutine

subroutine generate_potential_PotentialDataExample(this,inputs, &
   & sampling_points_dir,logfile,read_lambda)
  implicit none
  
  class(PotentialDataExample), intent(inout) :: this
  type(AnharmonicData),        intent(in)    :: inputs
  type(String),                intent(in)    :: sampling_points_dir
  type(OFile),                 intent(inout) :: logfile
  procedure(ReadLambda)                      :: read_lambda
  
  call print_line('PotentialDataExample: generating potential.')
  call print_line('Example contents = '//this%example_contents)
  
  ! Code to generate sampling points goes here.
end subroutine

! The class in use.
subroutine potential_example_subroutine()
  implicit none
  
  ! A polymorphic pointer, which can store an object of any type which
  !    extends PotentialData.
  type(PotentialPointer) :: potential
  
  ! An example variable with which to initialise a PotentialDataExample.
  type(String) :: example_contents
  
  ! Variables for generate_sampling points.
  type(AnharmonicData) :: anharmonic_data
  type(String)         :: sampling_points_dir
  type(OFile)          :: logfile
  
  ! Set the pointer to point to a PotentialDataExample type.
  ! This is where any PotentialDataExample-specific data is input,
  !    in this case the variable example_contents.
  example_contents = 'example'
  potential = PotentialDataExample(example_contents)
  
  ! Now PotentialData's methods can be called.
  ! They will all be forwareded to the PotentialDataExample instance.
  
  ! Generates sampling points, in a manner specific to the representation.
  call potential%generate_sampling_points( anharmonic_data,     &
                                         & sampling_points_dir, &
                                         & logfile,             &
                                         & write_structure_file_lambda)
  
  ! Code to run electronic structure goes here.
  
  ! Generates the potential, in a manner specific to the representation.
  call potential%generate_potential( anharmonic_data,     &
                                   & sampling_points_dir, &
                                   & logfile,             &
                                   & read_electronic_structure_lambda)
contains
  ! Lambda of type WriteLambda.
  subroutine write_structure_file_lambda(structure,directory)
    implicit none
    
    type(StructureData), intent(in) :: structure
    type(String),        intent(in) :: directory
    
    call print_line('Writing structure in directory '//directory)
    
    ! Code to write structure goes here.
  end subroutine
  
  ! Lambda of type ReadLambda.
  function read_electronic_structure_lambda(directory) result(output)
    implicit none
    
    type(String), intent(in)  :: directory
    type(ElectronicStructure) :: output
    
    call print_line('Reading electronic structure from directory '//directory)
    
    ! Code to read electronic structure goes here.
  end function
end subroutine
end module
