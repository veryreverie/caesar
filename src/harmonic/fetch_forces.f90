! Program to fetch forces from qe output
module fetch_forces_module
contains

subroutine fetch_forces(dft_code,dft_dir,seedname,atom,disp,forces_filename)
  use string_module
  use dft_output_file_module
  use file_module
  implicit none
  
  type(String), intent(in) :: dft_code
  type(String), intent(in) :: dft_dir
  type(String), intent(in) :: seedname
  integer,      intent(in) :: atom
  integer,      intent(in) :: disp
  type(String), intent(in) :: forces_filename
  
  ! file contents
  type(DftOutputFile) :: dft_output
  
  ! file units
  integer :: forces_file
  
  ! temporary variables
  integer  :: i,j
  
  dft_output = read_dft_output_file(dft_code,dft_dir,seedname)
  
  forces_file = open_write_file(forces_filename)
  do i=1,dft_output%no_atoms
    do j=1,3
      write(forces_file,*) char(str(atom)//' '//disp//' '//i//' '//j//' '//&
        & dft_output%forces(j,i))
    enddo
  enddo
  close(forces_file)
end subroutine
end module
