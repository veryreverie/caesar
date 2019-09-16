! ======================================================================
! Converts input files to and from StructureData.
! ======================================================================
module electronic_structure_file_module
  use utils_module
  
  use structure_module
  
  use electronic_structure_data_module
  use electronic_structure_common_module
  use quip_module
  use castep_module
  use qe_module
  use vasp_module
  use xyz_module
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
  elseif (file_type == 'quantum_espresso') then
    output = make_input_filename_qe(seedname)
  elseif (file_type == 'xyz') then
    output = make_input_filename_xyz(seedname)
  else
    call print_line('Unrecognised input file type: '//file_type)
    call err()
  endif
end function

function make_input_filename_caesar() result(output)
  implicit none
  
  type(String) :: output
  
  output = 'structure.dat'
end function

function make_output_filename(file_type,seedname) result(output)
  implicit none
  
  type(String), intent(in) :: file_type
  type(String), intent(in) :: seedname
  type(String)             :: output
  
  if (file_type == 'caesar') then
    output = make_output_filename_caesar()
  elseif (file_type == 'castep') then
    output = make_output_filename_castep(seedname)
  elseif (file_type == 'vasp') then
    output = make_output_filename_vasp(seedname)
  elseif (file_type == 'quantum_espresso') then
    output = make_output_filename_qe(seedname)
  elseif (file_type == 'xyz') then
    output = make_output_filename_caesar()
  else
    call print_line('Unrecognised output file type: '//file_type)
    call err()
  endif
end function

function make_output_filename_caesar() result(output)
  implicit none
  
  type(String) :: output
  
  output = 'electronic_structure.dat'
end function

! Reads an input file, and constructs a StructureData.
function input_file_to_StructureData(file_type,filename) result(output)
  implicit none
  
  type(String), intent(in) :: file_type
  type(String), intent(in) :: filename
  type(StructureData)      :: output
  
  type(IFile)          :: structure_file
  type(BasicStructure) :: basic_structure
  
  if (file_type == 'caesar') then
    structure_file = IFile(filename)
    output = StructureData(structure_file%lines())
  else
    if (file_type == 'castep') then
      basic_structure = read_input_file_castep(filename)
    elseif (file_type == 'quantum_espresso') then
      basic_structure = read_input_file_qe(filename)
    elseif (file_type == 'xyz') then
      basic_structure = read_input_file_xyz(filename)
    else
      call print_line('Reading '//file_type//' input files not yet supported.')
      call err()
    endif
    output = StructureData(basic_structure)
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
  
  type(OFile) :: structure_file
  
  if (file_type=='caesar') then
    structure_file = OFile(output_filename)
    call structure_file%print_lines(structure)
  elseif (file_type=='castep') then
    call write_input_file_castep(structure,input_filename,output_filename)
  elseif (file_type=='quantum_espresso') then
    call write_input_file_qe(structure,input_filename,output_filename)
  elseif (file_type=='xyz') then
    call write_input_file_xyz( BasicStructure(structure),     &
                             & input_filename,output_filename )
  else
    call print_line('Writing '//file_type//' input files not yet supported.')
    call err()
  endif
end subroutine

function read_output_file(file_type,structure,directory,seedname, &
   & calculation_type,use_forces,use_hessians,calculate_stress) result(output)
  implicit none
  
  type(String),        intent(in) :: file_type
  type(StructureData), intent(in) :: structure
  type(String),        intent(in) :: directory
  type(String),        intent(in) :: seedname
  type(String),        intent(in) :: calculation_type
  logical,             intent(in) :: use_forces
  logical,             intent(in) :: use_hessians
  logical,             intent(in) :: calculate_stress
  type(ElectronicStructure)       :: output
  
  type(IFile)         :: electronic_structure_file
  type(IFile)         :: displaced_structure_file
  type(StructureData) :: displaced_structure
  
  if (calculation_type=='none') then
    if (file_type=='caesar' .or. file_type=='xyz') then
      electronic_structure_file = IFile(directory//'/electronic_structure.dat')
      output = ElectronicStructure(electronic_structure_file%lines())
    elseif (file_type=='castep') then
      output = read_output_file_castep( directory,       &
                                      & seedname,        &
                                      & structure,       &
                                      & use_forces,      &
                                      & use_hessians,    &
                                      & calculate_stress )
    elseif (file_type=='quantum_espresso') then
      output = read_output_file_qe( directory,       &
                                  & seedname,        &
                                  & structure,       &
                                  & use_forces,      &
                                  & use_hessians,    &
                                  & calculate_stress )
    elseif (file_type=='vasp') then
      output = read_output_file_vasp( directory,       &
                                    & seedname,        &
                                    & structure,       &
                                    & use_forces,      &
                                    & use_hessians,    &
                                    & calculate_stress )
    else
      call print_line(ERROR//': Unrecognised output file type: '//file_type)
      call err()
    endif
  elseif (calculation_type=='quip') then
    displaced_structure_file = IFile(directory//'/structure.dat')
    displaced_structure = StructureData(displaced_structure_file%lines())
    output = run_quip_on_structure( BasicStructure(displaced_structure), &
                                  & seedname,                            &
                                  & use_forces,                          &
                                  & use_hessians,                        &
                                  & calculate_stress                     )
  else
    call print_line(ERROR//': calculation_type must be either "script" or &
       & "quip".')
    call err()
  endif
end function
end module
