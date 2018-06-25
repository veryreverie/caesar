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
    type(String), private              :: file_type_
    type(String), private              :: seedname_
    type(String), private              :: calculation_type_
    real(dp),     private              :: symmetry_precision_
    type(String), private              :: filename_
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
function new_CalculationReader(working_directory,file_type,seedname, &
   & calculation_type,symmetry_precision) result(this)
  implicit none
  
  type(String), intent(in) :: working_directory
  type(String), intent(in) :: file_type
  type(String), intent(in) :: seedname
  type(String), intent(in) :: calculation_type
  real(dp),     intent(in) :: symmetry_precision
  type(CalculationReader)  :: this
    
  this%working_directory_  = working_directory
  this%file_type_          = file_type
  this%seedname_           = seedname
  this%calculation_type_   = calculation_type
  this%symmetry_precision_ = symmetry_precision
  
  ! If the calculation type is 'script' then the electronic structure
  !    calculation has already been run, so the output file should be read
  ! If the calculation type is 'quip' then the calculation is still to be run,
  !    so the input file should be read.
  if (calculation_type=='script') then
    this%filename_ = make_output_filename(file_type,seedname)
  elseif (calculation_type=='quip') then
    this%filename_ = make_input_filename(file_type,seedname)
  else
    call print_line(ERROR//': calculation_type must be either "script" or &
       & "quip.')
    call err()
  endif
  
  this%directories_ = [String::]
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
  
  type(StructureData) :: structure
  
  type(OFile) :: electronic_structure_file
  
  ! Read in structure.
  structure = read_structure_file( directory//'/structure.dat', &
                                 & this%symmetry_precision_,    &
                                 & calculate_symmetry=.false.   )
  
  ! Read the calculation.
  output = read_output_file( this%file_type_,                &
                           & directory//'/'//this%filename_, &
                           & structure,                      &
                           & this%working_directory_,        &
                           & this%seedname_,                 &
                           & this%symmetry_precision_,       &
                           & this%calculation_type_          )
  
  ! Write an electronic_structure.dat file.
  electronic_structure_file = OFile(directory//'/electronic_structure.dat')
  call electronic_structure_file%print_lines(output)
  
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
