module caesar_atom_mapping_module
  use caesar_common_module
  implicit none
  
  interface
    ! ----------------------------------------------------------------------
    ! Finds a mapping between the atoms in two versions of the same structure.
    ! ----------------------------------------------------------------------
    module function atom_mapping(structure_a,structure_b) result(output) 
      type(StructureData), intent(in) :: structure_a
      type(StructureData), intent(in) :: structure_b
      type(Group)                     :: output
    end function
  end interface
end module
