module castep_output_file_module
  use common_module
  implicit none
  
  type, extends(ElectronicStructure) :: CastepOutputFile
    integer :: no_kpoints
    integer :: kpoints_mp_grid(3)
  end type
contains

function read_castep_output_file(filename,structure,dir,seedname, &
   & symmetry_precision) result(output)
  implicit none
  
  type(String),        intent(in) :: filename
  type(StructureData), intent(in) :: structure
  type(String),        intent(in) :: dir
  type(String),        intent(in) :: seedname
  real(dp),            intent(in) :: symmetry_precision
  type(CastepOutputFile)          :: output
  
  type(ElectronicStructure) :: output_file
  type(IFile)               :: castep_output_file
  
  type(String), allocatable :: line(:)
  integer                   :: i
  
  output_file = read_output_file( str('castep'),      &
                                & filename,           &
                                & structure,          &
                                & dir,                &
                                & seedname,           &
                                & symmetry_precision, &
                                & str('script'))
  
  output%energy = output_file%energy
  output%forces = output_file%forces
  
  castep_output_file = filename
  do i=1,size(castep_output_file)
    line = split(lower_case(castep_output_file%line(i)))
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
