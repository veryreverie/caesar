! ======================================================================
! A class which runs electronic structure calculations.
! ======================================================================
module calculation_runner_module
  use utils_module
  
  use structure_module
  
  use structure_file_module
  use electronic_structure_data_module
  use electronic_structure_file_module
  implicit none
  
  private
  
  public :: CalculationRunner
  
  type, extends(NoDefaultConstructor) :: CalculationRunner
    type(String), private              :: file_type_
    type(String), private              :: seedname_
    type(String), private              :: run_script_
    integer,      private              :: no_cores_
    integer,      private              :: no_nodes_
    type(String), private              :: run_script_data_
    type(String), private              :: calculation_type_
    logical,      private              :: exit_on_error_
    logical,      private              :: repeat_calculations_
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
   & no_nodes,run_script_data,calculation_type,exit_on_error,          &
   & repeat_calculations) result(this)
  implicit none
  
  type(String), intent(in) :: file_type
  type(String), intent(in) :: seedname
  type(String), intent(in) :: run_script
  integer,      intent(in) :: no_cores
  integer,      intent(in) :: no_nodes
  type(String), intent(in) :: run_script_data
  type(String), intent(in) :: calculation_type
  logical,      intent(in) :: exit_on_error
  logical,      intent(in) :: repeat_calculations
  type(CalculationRunner)  :: this
    
  this%file_type_           = file_type
  this%seedname_            = seedname
  this%run_script_          = run_script
  this%no_cores_            = no_cores
  this%no_nodes_            = no_nodes
  this%run_script_data_     = run_script_data
  this%calculation_type_    = calculation_type
  this%exit_on_error_       = exit_on_error
  this%repeat_calculations_ = repeat_calculations
  
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
  
  ! Check that this directory has not already been run.
  if (any(this%directories_==directory)) then
    call print_line(CODE_ERROR//': Trying to run an electronic structure &
       &calculation in a directory in which an electronic structure &
       &calculation has already been run by this execution of Caesar.')
    call print_line('Directory: '//directory)
    call err()
  endif
  
  ! Record the directory.
  this%directories_ = [this%directories_, directory]
  
  ! Check if the calculation has been run succesfully before.
  if (.not. this%repeat_calculations_) then
    if (file_exists(directory//'/electronic_structure.dat')) then
      call print_line('Skipping successful calculation in directory '// &
         & directory)
      return
    endif
  endif
  
  ! Run the calculation.
  call print_line('Running calculation in directory '//directory)
  result_code = system_call( 'cd '//working_directory//';' //' '// &
                           & this%run_script_              //' '// &
                           & this%file_type_               //' '// &
                           & directory                     //' '// &
                           & this%no_cores_                //' '// &
                           & this%no_nodes_                //' '// &
                           & this%seedname_                //' '// &
                           & this%run_script_data_                 )
  call print_line('Result code: '//result_code)
  
  ! Check the result code.
  if (result_code/=0) then
    if (this%exit_on_error_) then
      stop
    else
      return
    endif
  endif
  
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
