! ======================================================================
! A class which runs electronic structure calculations.
! ======================================================================
module caesar_calculation_runner_module
  use caesar_utils_module
  
  use caesar_structure_module
  
  use caesar_electronic_structure_data_module
  use caesar_electronic_structure_common_module
  
  use caesar_electronic_structure_file_module
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
    logical,      private              :: use_forces_
    logical,      private              :: use_hessians_
    logical,      private              :: calculate_stress_
    logical,      private              :: exit_on_error_
    logical,      private              :: repeat_calculations_
    type(String), private, allocatable :: directories_(:)
  contains
    procedure, public :: directories_run
    procedure, public :: run_calculation
    procedure, public :: run_calculations
  end type
  
  interface CalculationRunner
    ! Constructor.
    module function new_CalculationRunner(file_type,seedname,run_script, &
       & no_cores,no_nodes,run_script_data,calculation_type,             &
       & calculate_stress,use_forces,use_hessians,exit_on_error,         &
       & repeat_calculations) result(this) 
      type(String), intent(in) :: file_type
      type(String), intent(in) :: seedname
      type(String), intent(in) :: run_script
      integer,      intent(in) :: no_cores
      integer,      intent(in) :: no_nodes
      type(String), intent(in) :: run_script_data
      type(String), intent(in) :: calculation_type
      logical,      intent(in) :: use_forces
      logical,      intent(in) :: use_hessians
      logical,      intent(in) :: calculate_stress
      logical,      intent(in) :: exit_on_error
      logical,      intent(in) :: repeat_calculations
      type(CalculationRunner)  :: this
    end function
  end interface
  
  interface
    ! Return a list of the directories in which calculations have been run.
    module function directories_run(this) result(output) 
      class(CalculationRunner), intent(in) :: this
      type(String), allocatable            :: output(:)
    end function
  end interface
  
  interface
    ! Run a calculation in the given directory, and record the directory.
    module subroutine run_calculation(this,directory) 
      class(CalculationRunner), intent(inout) :: this
      type(String),             intent(in)    :: directory
    end subroutine
  end interface
  
  interface
    ! Run a calculation in the given directories, and record the directories.
    module subroutine run_calculations(this,directories) 
      class(CalculationRunner), intent(inout) :: this
      type(String),             intent(in)    :: directories(:)
    end subroutine
  end interface
end module
