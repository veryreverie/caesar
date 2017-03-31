! ----------------------------------------------------------------------
! Calculates the symmetries of the structure.
! ----------------------------------------------------------------------
module calculate_symmetry_module
contains
subroutine calculate_symmetry(this,temp_cell_filename,temp_symm_filename)
  use string_module
  use file_module
  use structure_module
  use structure_to_dft_module
  implicit none
  
  type(StructureData), intent(inout) :: this
  type(String),        intent(in)    :: temp_cell_filename
  type(String),        intent(in)    :: temp_symm_filename
  
  call structure_to_dft( dft_code        = str('castep'), &
                       & structure_sc    = this,          &
                       & output_filename = temp_cell_filename )
  
  call system_call( &
     & 'cellsym --symmetry '//temp_cell_filename//' > '//temp_symm_filename)
  
  call read_symmetry_file(this, temp_symm_filename)
  
  call system_call('rm '//temp_cell_filename)
  call system_call('rm '//temp_symm_filename)
end subroutine
end module
