! ======================================================================
! Converts input files to and from StructureData.
! ======================================================================
module electronic_structure_file_submodule
  use utils_module
  
  use structure_module
  
  use electronic_structure_data_submodule
  use structure_file_submodule
  use castep_wrapper_submodule
  use qe_wrapper_submodule
  use vasp_wrapper_submodule
  use quip_wrapper_submodule
  implicit none
  
  private
  
  public :: make_input_filename
  public :: make_output_filename
  
  public :: input_file_to_StructureData
  public :: StructureData_to_input_file
  
  public :: read_output_file
contains

! Converts a file seedname into the appropriate dft input filename.
function make_input_filename(file_type,seedname) result(output)
  implicit none
  
  type(String), intent(in) :: file_type
  type(String), intent(in) :: seedname
  type(String)             :: output
  
  if (file_type == 'caesar') then
    output = make_input_filename_caesar()
  elseif (file_type == 'castep' .or. file_type == 'quip') then
    output = make_input_filename_castep(seedname)
  elseif (file_type == 'vasp') then
    output = make_input_filename_vasp(seedname)
  elseif (file_type == 'qe') then
    output = make_input_filename_qe(seedname)
  elseif (file_type == 'xyz') then
    output = make_input_filename_xyz(seedname)
  else
    call print_line('Unrecognised input file type: '//file_type)
    call err()
  endif
end function

function make_output_filename(file_type,seedname) result(output)
  implicit none
  
  type(String), intent(in) :: file_type
  type(String), intent(in) :: seedname
  type(String)             :: output
  
  if (file_type == 'castep') then
    output = make_output_filename_castep(seedname)
  elseif (file_type == 'vasp') then
    output = make_output_filename_vasp(seedname)
  elseif (file_type == 'qe') then
    output = make_output_filename_qe(seedname)
  else
    call print_line('Unrecognised output file type: '//file_type)
    call err()
  endif
end function

! Reads an input file, and constructs a StructureData.
function input_file_to_StructureData(file_type,filename, &
   & symmetry_precision,calculate_symmetry) result(output)
  implicit none
  
  type(String), intent(in)           :: file_type
  type(String), intent(in)           :: filename
  real(dp),     intent(in)           :: symmetry_precision
  logical,      intent(in), optional :: calculate_symmetry
  type(StructureData)                :: output
  
  type(BasicStructure) :: basic_structure
  
  if (file_type == 'caesar') then
    output = read_structure_file( filename,           &
                                & symmetry_precision, &
                                & calculate_symmetry)
  else
    if (file_type == 'castep') then
      basic_structure = read_input_file_castep(filename)
    elseif (file_type == 'xyz') then
      basic_structure = read_input_file_xyz(filename)
    else
      call print_line('Reading '//file_type//' input files not yet supported.')
      call err()
    endif
    output = StructureData( basic_structure,    &
                          & symmetry_precision, &
                          & calculate_symmetry=calculate_symmetry)
  endif
end function

! Writes an input file from a StructureData.
subroutine StructureData_to_input_file(file_type,structure,input_filename, &
   & output_filename)
  implicit none
  
  type(String),        intent(in)           :: file_type
  type(StructureData), intent(in)           :: structure
  type(String),        intent(in), optional :: input_filename
  type(String),        intent(in)           :: output_filename
  
  if (file_type=='caesar') then
    call write_structure_file(structure,output_filename)
  elseif (file_type=='castep') then
    call write_input_file_castep(structure,input_filename,output_filename)
  elseif (file_type=='xyz') then
    call write_input_file_xyz(structure,input_filename,output_filename)
  else
    call print_line('Writing '//file_type//' input files not yet supported.')
    call err()
  endif
end subroutine

function read_output_file(file_type,filename,structure,dir,seedname, &
   & symmetry_precision,calculation_type) result(output)
  implicit none
  
  type(String),        intent(in) :: file_type
  type(String),        intent(in) :: filename
  type(StructureData), intent(in) :: structure
  type(String),        intent(in) :: dir
  type(String),        intent(in) :: seedname
  real(dp),            intent(in) :: symmetry_precision
  type(String),        intent(in) :: calculation_type
  type(ElectronicStructure)       :: output
  
  type(StructureData) :: displaced_structure
  
  if (calculation_type=='script') then
    if (file_type=='castep') then
      output = read_output_file_castep(filename,structure)
    elseif (file_type=='qe') then
      output = read_output_file_qe(filename,structure)
    elseif (file_type=='vasp') then
      output = read_output_file_vasp(filename,structure)
    else
      call print_line(ERROR//': Unrecognised output file type: '//file_type)
      call err()
    endif
  elseif (calculation_type=='quip') then
    displaced_structure = input_file_to_StructureData( &
                                 & file_type,          &
                                 & filename,           &
                                 & symmetry_precision, &
                                 & calculate_symmetry=.false.)
    output = run_quip_on_structure( displaced_structure, &
                                  & dir,                 &
                                  & seedname)
  else
    call print_line(ERROR//': calculation_type must be either "script" or &
       & "quip".')
    call err()
  endif
end function
end module
