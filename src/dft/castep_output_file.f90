module caesar_castep_output_file_module
  use caesar_common_module
  implicit none
  
  type, extends(NoDefaultConstructor) :: CastepOutputFile
    real(dp)             :: energy
    type(CartesianForce) :: forces
    integer              :: no_kpoints
    integer              :: kpoints_mp_grid(3)
  end type
  
  interface
    module function read_castep_output_file(structure,directory,seedname) &
       & result(output) 
      type(StructureData), intent(in) :: structure
      type(String),        intent(in) :: directory
      type(String),        intent(in) :: seedname
      type(CastepOutputFile)          :: output
    end function
  end interface
end module
