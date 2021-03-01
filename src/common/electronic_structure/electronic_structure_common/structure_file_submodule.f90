submodule (caesar_structure_file_module) caesar_structure_file_submodule
  use caesar_electronic_structure_common_module
contains

module procedure read_structure_file
  type(IFile) :: structure_file
  
  structure_file = IFile(filename)
  this = StructureData(structure_file%lines())
end procedure

module procedure write_structure_file
  type(OFile) :: structure_file
  
  structure_file = OFile(filename)
  call structure_file%print_lines(this)
end procedure
end submodule
