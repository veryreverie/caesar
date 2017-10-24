module castep_output_file_module
  use constants_module, only : dp
  use string_module
  use io_module
  
  use dft_output_file_module
  implicit none
  
  type, extends(DftOutputFile) :: CastepOutputFile
    integer :: no_kpoints
    integer :: kpoints_mp_grid(3)
  end type
contains

function read_castep_output_file(filename) result(output)
  use ifile_module
  implicit none
  
  type(String), intent(in) :: filename
  type(CastepOutputFile)   :: output
  
  type(DftOutputFile) :: dft_output_file
  type(IFile)         :: castep_output_file
  
  type(String), allocatable :: line(:)
  integer                   :: i
  
  dft_output_file = read_dft_output_file(str('castep'),filename)
  
  output%no_atoms = dft_output_file%no_atoms
  output%species = dft_output_file%species
  output%energy = dft_output_file%energy
  output%forces = dft_output_file%forces
  
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
