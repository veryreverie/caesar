! Program to fetch forces from qe output
module fetch_forces_qe_module
contains

subroutine fetch_forces_qe(qe_output_filename,atom,disp,forces_filename)
  use string_module
  use qe_output_file_module
  use constants, only : Ry, bohr
  use file_io,   only : open_write_file
  implicit none
  
  type(String), intent(in) :: qe_output_filename
  integer,      intent(in) :: atom
  integer,      intent(in) :: disp
  type(String), intent(in) :: forces_filename
  
  ! file contents
  type(QeOutputFile) :: qe_output
  
  ! file units
  integer :: forces_file
  
  ! temporary variables
  integer :: i,j
  
  qe_output = read_qe_output_file(qe_output_filename)
  
  forces_file = open_write_file(forces_filename)
  do i=1,qe_output%no_atoms
    do j=1,3
      write(forces_file,*) char(str(atom)//' '//disp//' '//i//' '//j//' '//&
        & qe_output%forces%forces(j,i)*Ry/bohr)
    enddo
  enddo
  close(forces_file)
end subroutine
end module
