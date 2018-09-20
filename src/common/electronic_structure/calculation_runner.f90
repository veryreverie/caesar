! ======================================================================
! A class which runs electronic structure calculations.
! ======================================================================
module calculation_runner_submodule
  use utils_module
  
  use structure_module
  
  use structure_file_submodule
  use electronic_structure_data_submodule
  use electronic_structure_file_submodule
  implicit none
  
  private
  
  public :: CalculationRunner
  
  type, extends(NoDefaultConstructor) :: CalculationRunner
    type(String), private              :: file_type_
    type(String), private              :: seedname_
    type(String), private              :: run_script_
    integer,      private              :: no_cores_
    type(String), private              :: calculation_type_
    type(String), private              :: filename_
    type(String), private, allocatable :: directories_(:)
  contains
    procedure, public :: directories_run
    procedure, public :: run_calculation
    procedure, public :: run_calculations
  end type
  
  interface CalculationRunner
    module procedure new_CalculationRunner
  end interface
contains

! Constructor.
function new_CalculationRunner(file_type,seedname,run_script,no_cores, &
   & calculation_type) result(this)
  implicit none
  
  type(String), intent(in) :: file_type
  type(String), intent(in) :: seedname
  type(String), intent(in) :: run_script
  integer,      intent(in) :: no_cores
  type(String), intent(in) :: calculation_type
  type(CalculationRunner)  :: this
    
  this%file_type_        = file_type
  this%seedname_         = seedname
  this%run_script_       = run_script
  this%no_cores_         = no_cores
  this%calculation_type_ = calculation_type
  
  ! If the calculation type is 'none' then the electronic structure
  !    calculation has already been run, so the output file should be read
  ! If the calculation type is 'quip' then the calculation is still to be run,
  !    so the input file should be read.
  if (calculation_type=='none') then
    this%filename_ = make_output_filename(file_type,seedname)
  elseif (calculation_type=='quip') then
    this%filename_ = make_input_filename(file_type,seedname)
  else
    call print_line(ERROR//': calculation_type must be either "none" or &
       & "quip.')
    call err()
  endif
  
  this%directories_ = [String::]
end function

! Return a list of the directories in which calculations have been run.
function directories_run(this) result(output)
  implicit none
  
  class(CalculationRunner), intent(in) :: this
  type(String), allocatable            :: output(:)
  
  output = this%directories_
end function

! Run a calculation in the given directory, and record the directory.
subroutine run_calculation(this,directory)
  implicit none
  
  class(CalculationRunner), intent(inout) :: this
  type(String),             intent(in)    :: directory
  
  type(String) :: working_directory
  
  integer                   :: result_code
  type(IFile)               :: structure_file
  type(StructureData)       :: structure
  type(ElectronicStructure) :: electronic_structure
  type(OFile)               :: electronic_structure_file
  
  working_directory = format_path('.')
  
  ! Run the calculation.
  call print_line('Running calculation in directory '//directory)
  result_code = system_call( 'cd '//working_directory//';' //' '// &
                           & this%run_script_              //' '// &
                           & this%file_type_               //' '// &
                           & directory                     //' '// &
                           & this%no_cores_                //' '// &
                           & this%seedname_                        )
  call print_line('Result code: '//result_code)
  
  ! Convert the electronic structure result into an ElectronicStructure.
  structure_file = IFile(directory//'/structure.dat')
  structure = StructureData(structure_file%lines())
  electronic_structure = read_output_file( &
         & this%file_type_,                &
         & directory//'/'//this%filename_, &
         & structure,                      &
         & directory,                      &
         & this%seedname_,                 &
         & this%calculation_type_          )
  electronic_structure_file = OFile(directory//'/electronic_structure.dat' )
  call electronic_structure_file%print_lines(electronic_structure)
  
  ! Record the directory.
  this%directories_ = [this%directories_, directory]
end subroutine

! Run a calculation in the given directories, and record the directories.
subroutine run_calculations(this,directories)
  implicit none
  
  class(CalculationRunner), intent(inout) :: this
  type(String),             intent(in)    :: directories(:)
  
  integer :: i
  
  do i=1,size(directories)
    call print_line(i//' of '//size(directories)//':')
    call this%run_calculation(directories(i))
  enddo
end subroutine
end module
