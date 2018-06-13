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
    
    procedure, public :: energy => energy_PotentialPointer
    procedure, public :: force  => force_PotentialPointer
    
    procedure, private :: check => check_PotentialPointer
    
    procedure, public :: read  => read_PotentialPointer
    procedure, public :: write => write_PotentialPointer
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

! Checks that the pointer has been allocated before it is used.
subroutine check_PotentialPointer(this)
  implicit none
  
  class(PotentialPointer), intent(in) :: this
  
  if (.not. allocated(this%potential)) then
    call print_line(CODE_ERROR//': Trying to use a PotentialPointer before &
       &it has been allocated.')
    call err()
  endif
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
  
  call this%check()
  
  call this%potential%generate_sampling_points( inputs,              &
                                              & sampling_points_dir, &
                                              & logfile,             &
                                              & write_lambda)
end subroutine

subroutine generate_potential_PotentialPointer(this,inputs, &
   & weighted_energy_force_ratio,sampling_points_dir,logfile,read_lambda)
  implicit none
  
  class(PotentialPointer), intent(inout) :: this
  type(AnharmonicData),    intent(in)    :: inputs
  real(dp),                intent(in)    :: weighted_energy_force_ratio
  type(String),            intent(in)    :: sampling_points_dir
  type(OFile),             intent(inout) :: logfile
  procedure(ReadLambda)                  :: read_lambda
  
  call this%check()
  
  call this%potential%generate_potential( inputs,                      &
                                        & weighted_energy_force_ratio, &
                                        & sampling_points_dir,         &
                                        & logfile,                     &
                                        & read_lambda)
end subroutine

impure elemental function energy_PotentialPointer(this,displacement) &
   & result(output)
  implicit none
  
  class(PotentialPointer),    intent(in) :: this
  type(RealModeDisplacement), intent(in) :: displacement
  real(dp)                               :: output
  
  call this%check()
  
  output = this%potential%energy(displacement)
end function

impure elemental function force_PotentialPointer(this,displacement) &
   & result(output)
  implicit none
  
  class(PotentialPointer),    intent(in) :: this
  type(RealModeDisplacement), intent(in) :: displacement
  type(RealModeForce)                    :: output
  
  call this%check()
  
  output = this%potential%force(displacement)
end function

subroutine read_PotentialPointer(this,input)
  implicit none
  
  class(PotentialPointer), intent(out) :: this
  type(String),            intent(in)  :: input(:)
  
  ! Since this%potential is polymorphic, it is not possible to determine its
  !    type in an intent(out) subroutine.
  ! As such, a PotentialPointer(StringArray) constructor is not provided.
  call print_line(CODE_ERROR//' Unable to read a PotentialPointer.')
  call err()
end subroutine

function write_PotentialPointer(this) result(output)
  implicit none
  
  class(PotentialPointer), intent(in) :: this
  type(String), allocatable           :: output(:)
  
  select type(this); type is(PotentialPointer)
    output = str(this%potential)
  end select
end function
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
    
    procedure, public :: energy => energy_PotentialDataExample
    procedure, public :: force  => force_PotentialDataExample
    
    procedure, public :: read  => read_PotentialDataExample
    procedure, public :: write => write_PotentialDataExample
  end type
  
  interface PotentialDataExample
    module procedure new_PotentialDataExample
    module procedure new_PotentialDataExample_StringArray
  end interface
contains

! --------------------------------------------------
! Constructor for example class.
! --------------------------------------------------
! This is where any PotentialDataExample-specific data is input.
function new_PotentialDataExample(example_contents) result(this)
  implicit none
  
  type(String), intent(in)   :: example_contents
  type(PotentialDataExample) :: this
  
  this%example_contents = example_contents
end function

! --------------------------------------------------
! Overloads of PotentialData's methods.
! --------------------------------------------------
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
   & weighted_energy_force_ratio,sampling_points_dir,logfile,read_lambda)
  implicit none
  
  class(PotentialDataExample), intent(inout) :: this
  type(AnharmonicData),        intent(in)    :: inputs
  real(dp),                    intent(in)    :: weighted_energy_force_ratio
  type(String),                intent(in)    :: sampling_points_dir
  type(OFile),                 intent(inout) :: logfile
  procedure(ReadLambda)                      :: read_lambda
  
  call print_line('PotentialDataExample: generating potential.')
  call print_line('Example contents = '//this%example_contents)
  
  ! Code to generate sampling points goes here.
end subroutine

impure elemental function energy_PotentialDataExample(this,displacement) &
   & result(output)
  implicit none
  
  class(potentialDataExample), intent(in) :: this
  type(RealModeDisplacement),  intent(in) :: displacement
  real(dp)                                :: output
  
  call print_line('PotentialDataExample: calculating energy.')
  call print_line('Example contents = '//this%example_contents)
  
  ! Code to calculate energies goes here.
end function

impure elemental function force_PotentialDataExample(this,displacement) &
   & result(output)
  implicit none
  
  class(potentialDataExample), intent(in) :: this
  type(RealModeDisplacement),  intent(in) :: displacement
  type(RealModeForce)                     :: output
  
  call print_line('PotentialDataExample: calculating force.')
  call print_line('Example contents = '//this%example_contents)
  
  ! Code to calculate forces goes here.
end function

! --------------------------------------------------
! I/O.
! --------------------------------------------------
subroutine read_PotentialDataExample(this,input)
  implicit none
  
  class(PotentialDataExample), intent(out) :: this
  type(String),                intent(in)  :: input(:)
  
  ! Code to read potential from strings goes here.
end subroutine

function write_PotentialDataExample(this) result(output)
  implicit none
  
  class(PotentialDataExample), intent(in) :: this
  type(String), allocatable               :: output(:)
  
  ! Code to write potential to strings goes here.
end function

impure elemental function new_PotentialDataExample_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(PotentialDataExample)    :: this
  
  this = input
end function

! The class in use.
subroutine potential_example_subroutine(wd)
  implicit none
  
  type(String), intent(in) :: wd
  
  ! A polymorphic pointer, which can store an object of any type which
  !    extends PotentialData.
  type(PotentialPointer) :: potential
  
  ! An example variable with which to initialise a PotentialDataExample.
  type(String) :: example_contents
  
  ! Variables for generate_sampling points.
  type(AnharmonicData) :: anharmonic_data
  type(String)         :: sampling_points_dir
  type(OFile)          :: logfile
  
  ! Variables for generate_potential.
  real(dp) :: weighted_energy_force_ratio
  
  ! Variables for energy and force.
  type(RealModeDisplacement) :: displacement
  real(dp)                   :: energy
  type(RealModeForce)        :: force
  
  ! Files.
  type(OFile) :: output_file
  type(IFile) :: input_file
  
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
  call potential%generate_potential( anharmonic_data,             &
                                   & weighted_energy_force_ratio, &
                                   & sampling_points_dir,         &
                                   & logfile,                     &
                                   & read_electronic_structure_lambda)
  
  ! Once the potential has been generated, it can be used to calculate
  !    energies and forces.
  energy = potential%energy(displacement)
  force  = potential%force(displacement)
  
  ! The potential can be written to file directly from the potential pointer's
  !    write method, but must be read using the specific potential's
  !    constructor.
  output_file = OFile(wd//'example_potential.file')
  call output_file%print_lines(potential)
  
  input_file = IFile(wd//'example_potential.file')
  potential = PotentialDataExample(StringArray(input_file%lines()))
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
