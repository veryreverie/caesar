module caesar_castep_output_file_module
  use caesar_common_module
  implicit none
  
  type, extends(NoDefaultConstructor) :: CastepOutputFile
    real(dp)             :: energy
    type(CartesianForce) :: forces
    integer              :: no_kpoints
    integer              :: kpoints_mp_grid(3)
  end type
contains

function read_castep_output_file(structure,directory,seedname) result(output)
  implicit none
  
  type(StructureData), intent(in) :: structure
  type(String),        intent(in) :: directory
  type(String),        intent(in) :: seedname
  type(CastepOutputFile)          :: output
  
  type(String)              :: filename
  type(ElectronicStructure) :: output_file
  type(IFile)               :: castep_output_file
  
  type(String), allocatable :: line(:)
  integer                   :: i
  
  output_file = read_output_file( file_type        = str('castep'), &
                                & structure        = structure,     &
                                & directory        = directory,     &
                                & seedname         = seedname,      &
                                & calculation_type = str('script'), &
                                & use_forces       = .true.,        &
                                & use_hessians     = .false.,       &
                                & calculate_stress = .false.        )
  
  output%energy = output_file%energy()
  output%forces = output_file%forces()
  
  filename = directory//'/'//make_output_filename( file_type = str('castep'), &
                                                 & seedname  = seedname       )
  castep_output_file = IFile(filename)
  do i=1,size(castep_output_file)
    line = split_line(lower_case(castep_output_file%line(i)))
    if (size(line)>=6) then
      if (join(line(:3))=='number of kpoints') then
        output%no_kpoints = int(line(6))
      elseif (join(line(:3))=='mp grid size') then
        output%kpoints_mp_grid = int(line(size(line)-2:))
      endif
    endif
  enddo
end function
end module
