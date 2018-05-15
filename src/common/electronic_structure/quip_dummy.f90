! ======================================================================
! A dummy module for when QUIP is not being linked against.
! ======================================================================
module quip_wrapper_submodule
  use utils_module
  
  use structure_module
  
  use electronic_structure_data_submodule
  implicit none
  
  private
  
  public :: make_input_filename_xyz
  public :: read_input_file_xyz
  public :: write_input_file_xyz
  public :: run_quip_on_structure
  public :: QUIP_LINKED
  
  logical, parameter :: QUIP_LINKED = .false.
contains

function make_input_filename_xyz(seedname) result(output)
  implicit none
  
  type(String), intent(in) :: seedname
  type(String)             :: output
  
  call print_line(ERROR//': Cannot read .xyz files because Caesar has not &
     &been linked against QUIP. Please use -DPATH_TO_QUIP when running cmake.')
  stop
end function

function read_input_file_xyz(filename) result(output)
  implicit none
  
  type(String), intent(in)  :: filename
  type(BasicStructure)      :: output
  
  call print_line(ERROR//': Cannot read .xyz files because Caesar has not &
     &been linked against QUIP. Please use -DPATH_TO_QUIP when running cmake.')
  stop
end function

subroutine write_input_file_xyz(structure,input_filename,output_filename)
  implicit none
  
  type(StructureData), intent(in)           :: structure
  type(String),        intent(in), optional :: input_filename
  type(String),        intent(in)           :: output_filename
  
  call print_line(ERROR//': Cannot write .xyz files because Caesar has not &
     &been linked against QUIP. Please use -DPATH_TO_QUIP when running cmake.')
  stop
end subroutine

function run_quip_on_structure(structure,dir,seedname) result(output)
  implicit none
  
  type(StructureData), intent(in) :: structure
  type(String),        intent(in) :: dir
  type(String),        intent(in) :: seedname
  type(ElectronicStructure)       :: output
  
  call print_line(ERROR//': Cannot run calculations using Quip because Caesar &
     &has not been linked against QUIP. Please use -DPATH_TO_QUIP when &
     &running cmake.')
  stop
end function
end module
