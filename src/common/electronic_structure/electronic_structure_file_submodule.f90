submodule (caesar_electronic_structure_file_module) caesar_electronic_structure_file_submodule
  use caesar_electronic_structure_module
contains

module procedure make_input_filename
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
end procedure

module procedure make_input_filename_caesar
  output = 'structure.dat'
end procedure

module procedure make_output_filename
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
end procedure

module procedure make_output_filename_caesar
  output = 'electronic_structure.dat'
end procedure

module procedure input_file_to_StructureData
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
end procedure

module procedure StructureData_to_input_file
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
end procedure

module procedure read_output_file
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
end procedure
end submodule
