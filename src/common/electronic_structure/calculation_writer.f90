! ======================================================================
! A class which writes and keeps track of electronic structure calculation
!    input directories.
! ======================================================================
module caesar_calculation_writer_module
  use caesar_utils_module
  
  use caesar_structure_module
  
  use caesar_electronic_structure_common_module
  
  use caesar_electronic_structure_file_module
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
    ! Constructor.
    module function new_CalculationWriter(file_type,seedname) result(this) 
      type(String), intent(in) :: file_type
      type(String), intent(in) :: seedname
      type(CalculationWriter)  :: this
    end function
  end interface
  
  interface
    ! Return a list of the directories written.
    module function directories_written(this) result(output) 
      class(CalculationWriter), intent(in) :: this
      type(String), allocatable            :: output(:)
    end function
  end interface
  
  interface
    ! Write a calculation directory, and record the directory written.
    module subroutine write_calculation(this,structure,directory) 
      class(CalculationWriter), intent(inout) :: this
      type(StructureData),      intent(in)    :: structure
      type(String),             intent(in)    :: directory
    end subroutine
  end interface
end module
