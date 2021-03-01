! ======================================================================
! Converts input files to and from StructureData.
! ======================================================================
module caesar_electronic_structure_file_module
  use caesar_utils_module
  
  use caesar_structure_module
  
  use caesar_electronic_structure_data_module
  use caesar_electronic_structure_common_module
  use caesar_quip_module
  use caesar_castep_module
  use caesar_qe_module
  use caesar_vasp_module
  use caesar_xyz_module
  implicit none
  
  private
  
  public :: make_input_filename
  public :: make_output_filename
  
  public :: input_file_to_StructureData
  public :: StructureData_to_input_file
  
  public :: read_output_file
  
  interface
    ! Converts a file seedname into the appropriate dft input filename.
    module function make_input_filename(file_type,seedname) result(output) 
      type(String), intent(in) :: file_type
      type(String), intent(in) :: seedname
      type(String)             :: output
    end function
  end interface
  
  interface
    module function make_input_filename_caesar() result(output) 
      type(String) :: output
    end function
  end interface
  
  interface
    module function make_output_filename(file_type,seedname) result(output) 
      type(String), intent(in) :: file_type
      type(String), intent(in) :: seedname
      type(String)             :: output
    end function
  end interface
  
  interface
    module function make_output_filename_caesar() result(output) 
      type(String) :: output
    end function
  end interface
  
  interface
    ! Reads an input file, and constructs a StructureData.
    module function input_file_to_StructureData(file_type,filename) &
       & result(output) 
      type(String), intent(in) :: file_type
      type(String), intent(in) :: filename
      type(StructureData)      :: output
    end function
  end interface
  
  interface
    ! Writes an input file from a StructureData.
    module subroutine StructureData_to_input_file(file_type,structure, &
       & input_filename,output_filename) 
      type(String),        intent(in)           :: file_type
      type(StructureData), intent(in)           :: structure
      type(String),        intent(in), optional :: input_filename
      type(String),        intent(in)           :: output_filename
    end subroutine
  end interface
  
  interface
    module function read_output_file(file_type,structure,directory,seedname, &
       & calculation_type,use_forces,use_hessians,calculate_stress)          &
       & result(output) 
      type(String),        intent(in) :: file_type
      type(StructureData), intent(in) :: structure
      type(String),        intent(in) :: directory
      type(String),        intent(in) :: seedname
      type(String),        intent(in) :: calculation_type
      logical,             intent(in) :: use_forces
      logical,             intent(in) :: use_hessians
      logical,             intent(in) :: calculate_stress
      type(ElectronicStructure)       :: output
    end function
  end interface
end module
