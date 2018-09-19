! ======================================================================
! A class which reads the results of electronic structure calculations.
! ======================================================================
module calculation_reader_submodule
  use utils_module
  
  use structure_module
  
  use structure_file_submodule
  use electronic_structure_data_submodule
  use electronic_structure_file_submodule
  implicit none
  
  private
  
  public :: CalculationReader
  
  type, extends(NoDefaultConstructor) :: CalculationReader
    type(String), private              :: working_directory_
    type(String), private, allocatable :: directories_(:)
  contains
    procedure, public :: directories_read
    procedure, public :: read_calculation
    procedure, public :: read_calculations
  end type
  
  interface CalculationReader
    module procedure new_CalculationReader
  end interface
contains

! Constructor.
function new_CalculationReader(working_directory) result(this)
  implicit none
  
  type(String), intent(in) :: working_directory
  type(CalculationReader)  :: this
  
  this%working_directory_ = working_directory
  this%directories_       = [String::]
end function

! Return a list of the directories from which calculations have been read.
function directories_read(this) result(output)
  implicit none
  
  class(CalculationReader), intent(in) :: this
  type(String), allocatable            :: output(:)
  
  output = this%directories_
end function

! Read the results of an electronic structure calculation from the given
!    directory, and record the directory.
function read_calculation(this,directory) result(output)
  implicit none
  
  class(CalculationReader), intent(inout) :: this
  type(String),             intent(in)    :: directory
  type(ElectronicStructure)               :: output
  
  type(String) :: absolute_directory
  type(IFile)  :: electronic_structure_file
  
  absolute_directory = this%working_directory_//'/'//directory
  electronic_structure_file = IFile(                   &
     & absolute_directory//'/electronic_structure.dat' )
  output = ElectronicStructure(electronic_structure_file%lines())
  
  ! Record the directory.
  this%directories_ = [this%directories_, directory]
end function

! Read the results of a set of electronic structure calculations from the given
!    directories, and record the directories.
subroutine read_calculations(this,directories)
  implicit none
  
  class(CalculationReader), intent(inout) :: this
  type(String),             intent(in)    :: directories(:)
  type(ElectronicStructure), allocatable  :: output(:)
  
  integer :: i,ialloc
  
  ! Read the calculations.
  allocate(output(size(directories)), stat=ialloc); call err(ialloc)
  do i=1,size(directories)
    output(i) = this%read_calculation(directories(i))
  enddo
  
  ! Record the directories.
  this%directories_ = [this%directories_, directories]
end subroutine
end module
