module generate_sc_path_module
  implicit none
contains

subroutine generate_sc_path(filenames)
  use constants, only : dp
  use utils,     only : lower_case
  use file_module
  use string_module
  implicit none
  
  type(String), intent(in) :: filenames(:)
  
  integer               :: i
  integer               :: supercell(3,3)
  integer               :: no_points
  real(dp), allocatable :: path(:,:)
  
  ! line numbers
  integer :: path_file_length
  integer :: path_start_line
  integer :: path_end_line
  
  ! file units
  integer :: supercell_file
  integer :: path_file
  integer :: sc_bs_path_file
  
  ! file names
  type(String) :: supercell_filename
  type(String) :: path_filename
  type(String) :: sc_bs_path_filename
  
  character(100) :: line
  
  supercell_filename = filenames(1)
  path_filename = filenames(2)
  sc_bs_path_filename = filenames(3)

  ! Read in supercell matrix
  supercell_file = open_read_file(supercell_filename)
  do i=1,3
    read(supercell_file,*) supercell(i,:)
  enddo
  close(supercell_file)
  
  ! Parse path file
  path_file_length = count_lines(path_filename)
  path_file = open_read_file(path_filename)
  do i=1,path_file_length
    read(path_file,"(a)") line
    line = lower_case(line)
    if (line(1:5)=="block") then
      path_start_line = i
    elseif (line(1:8)=="endblock") then
      path_end_line = i
    endif
  enddo
  
  no_points = path_end_line-path_start_line-1
  allocate(path(3,no_points))
  rewind(path_file)
  
  do i=1,path_file_length
    read(path_file,"(a)") line
    if(i>path_start_line .and. i<path_end_line) then
      read(line,*) path(:,i-path_start_line)
    endif
  enddo
  close(path_file)

  ! Write output
  sc_bs_path_file = open_write_file(filenames(3))
  do i=1,no_points
    write(sc_bs_path_file,*) matmul(supercell,path(:,i))
  enddo
  close(sc_bs_path_file)
end subroutine
end module
