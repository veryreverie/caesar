! Program to fetch forces from qe output
module fetch_forces_module
contains

subroutine fetch_forces(codename,dft_output_filename,atom,disp,forces_filename)
  use string_module
  use dft_output_file_module
  use constants, only : dp, Ry, bohr
  use file_module,   only : open_write_file
  implicit none
  
  type(String), intent(in) :: codename
  type(String), intent(in) :: dft_output_filename
  integer,      intent(in) :: atom
  integer,      intent(in) :: disp
  type(String), intent(in) :: forces_filename
  
  ! file contents
  type(DftOutputFile) :: dft_output
  
  ! file units
  integer :: forces_file
  
  ! temporary variables
  integer  :: i,j
  real(dp) :: conversion
  
  if (codename=="castep") then
    conversion = 1.d0
    dft_output = read_castep_output_file(dft_output_filename)
  elseif (codename=="qe") then
    conversion = Ry/bohr
    dft_output = read_qe_output_file(dft_output_filename)
  endif
  
  forces_file = open_write_file(forces_filename)
  do i=1,dft_output%no_atoms
    do j=1,3
      write(forces_file,*) char(str(atom)//' '//disp//' '//i//' '//j//' '//&
        & dft_output%forces(j,i)*conversion)
    enddo
  enddo
  close(forces_file)
end subroutine
end module
