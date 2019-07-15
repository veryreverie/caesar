! ======================================================================
! A class which writes and keeps track of electronic structure calculation
!    input directories.
! ======================================================================
module calculation_writer_module
  use utils_module
  
  use structure_module
  
  use electronic_structure_common_module
  
  use electronic_structure_file_module
  implicit none
  
  private
  
  public :: CalculationWriter
  
  type, extends(NoDefaultConstructor) :: CalculationWriter
    type(String), private              :: file_type_
    type(String), private              :: seedname_
    type(String), private              :: input_filename_
    type(String), private, allocatable :: directories_(:)
  contains
    procedure, public :: directories_written
    procedure, public :: write_calculation
  end type
  
  interface CalculationWriter
    module procedure new_CalculationWriter
  end interface
contains

! Constructor.
function new_CalculationWriter(file_type,seedname) result(this)
  implicit none
  
  type(String), intent(in) :: file_type
  type(String), intent(in) :: seedname
  type(CalculationWriter)  :: this
  
  integer :: ialloc
  
  this%file_type_      = file_type
  this%seedname_       = seedname
  this%input_filename_ = make_input_filename(file_type,seedname)
  allocate(this%directories_(0), stat=ialloc); call err(ialloc)
  
  ! Check that the input file exists.
  if (.not. file_exists(this%input_filename_)) then
    call print_line(ERROR//': The input file '//this%input_filename_// &
                   &' does not exist in the working directory.'        )
    call quit()
  endif
end function

! Return a list of the directories written.
function directories_written(this) result(output)
  implicit none
  
  class(CalculationWriter), intent(in) :: this
  type(String), allocatable            :: output(:)
  
  output = this%directories_
end function

! Write a calculation directory, and record the directory written.
subroutine write_calculation(this,structure,directory)
  implicit none
  
  class(CalculationWriter), intent(inout) :: this
  type(StructureData),      intent(in)    :: structure
  type(String),             intent(in)    :: directory
  
  type(OFile) :: structure_file
  
  ! Check that the directory has not already been written to by this class.
  if (any(this%directories_==directory)) then
    call print_line(CODE_ERROR//': Trying to write an electronic structure &
       &calculation to a directory in which an electronic structure &
       &calculation has already been written by this execution of Caesar.')
    call print_line('Directory: '//directory)
    call err()
  endif
  
  ! Make the directory, and add a structure.dat file and
  !    an electronic structure input file.
  call mkdir(directory)
  structure_file = OFile(directory//'/structure.dat')
  call structure_file%print_lines(structure)
  call StructureData_to_input_file(              &
     & this%file_type_,                          &
     & structure,                                &
     & this%input_filename_,                     &
     & directory//'/'//this%input_filename_      )
  
  ! Record the directory.
  this%directories_ = [this%directories_, directory]
end subroutine
end module
