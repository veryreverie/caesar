! ======================================================================
! A class which reads the results of electronic structure calculations.
! ======================================================================
module caesar_calculation_reader_module
  use caesar_utils_module
  
  use caesar_structure_module
  use caesar_normal_mode_module
  
  use caesar_electronic_structure_data_module
  use caesar_electronic_structure_common_module
  
  use caesar_electronic_structure_file_module
  implicit none
  
  private
  
  public :: CalculationReader
  
  type, extends(NoDefaultConstructor) :: CalculationReader
    type(FractionVector), private, allocatable :: loto_direction_
    type(String),         private, allocatable :: directories_(:)
  contains
    procedure, public :: directories_read
    procedure, public :: read_calculation
    procedure, public :: read_calculations
  end type
  
  interface CalculationReader
    ! Constructor.
    module function new_CalculationReader(loto_direction) result(this) 
      type(FractionVector), intent(in), optional :: loto_direction
      type(CalculationReader)                    :: this
    end function
  end interface
  
  interface
    ! Return a list of the directories from which calculations have been read.
    module function directories_read(this) result(output) 
      class(CalculationReader), intent(in) :: this
      type(String), allocatable            :: output(:)
    end function
  end interface
  
  interface
    ! Read the results of an electronic structure calculation from the given
    !    directory, and record the directory.
    module function read_calculation(this,directory,displacement) &
       & result(output) 
      class(CalculationReader),    intent(inout)        :: this
      type(String),                intent(in)           :: directory
      type(CartesianDisplacement), intent(in), optional :: displacement
      type(ElectronicStructure)                         :: output
    end function
  end interface
  
  interface
    ! Read the results of a set of electronic structure calculations from the given
    !    directories, and record the directories.
    module subroutine read_calculations(this,directories) 
      class(CalculationReader), intent(inout) :: this
      type(String),             intent(in)    :: directories(:)
    end subroutine
  end interface
end module
