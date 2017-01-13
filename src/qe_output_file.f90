! A class for holding the information in a qe .out file
module qe_output_file_module
  use constants, only : dp
  implicit none
  
  private
  
  public :: QeOutputFile
  public :: read_qe_output_file
  public :: new
  public :: drop
  
  type QeOutputForces
    integer,  allocatable :: atom_ids(:)
    integer,  allocatable :: atom_types(:)
    real(dp), allocatable :: forces(:,:)
  end type
  
  type QeOutputFile
    integer              :: no_atoms
    type(QeOutputForces) :: forces
  end type
  
  interface read_qe_output_file
    module procedure read_qe_output_file_character
    module procedure read_qe_output_file_string
  end interface
  
  interface new
    module procedure new_QeOutputForces
  end interface
  
  interface drop
    module procedure drop_QeOutputForces
    module procedure drop_QeOutputFile
  end interface
    
contains

function read_qe_output_file_character(filename) result(output)
  use file_io
  use utils, only : lower_case
  implicit none
  
  character(*), intent(in) :: filename
  type(QeOutputFile)       :: output
  
  integer        :: qe_file
  integer        :: file_length
  integer        :: i
  character(100) :: line
  integer        :: forces_start_line
  integer        :: forces_end_line
  character(100) :: temp_char
  
  qe_file = open_read_file(filename)
  file_length = count_lines(qe_file)
  
  do i=1,file_length
    read(qe_file,"(a)") line
    line = lower_case(trim(line))
    if (line(1:13)=="forces acting") then
      forces_start_line = i
    elseif (line(1:11)=="total force") then
      forces_end_line = i
    endif
  enddo
  
  output%no_atoms = forces_end_line-forces_start_line-3
  call new(output%forces,output%no_atoms)
  
  rewind(qe_file)
  do i=1,forces_start_line+1
    read(qe_file,*)
  enddo
  do i=1,output%no_atoms
    read(qe_file,*) temp_char,                   &
                  & output%forces%atom_ids(i),   &
                  & temp_char,                   &
                  & output%forces%atom_types(i), &
                  & temp_char,                   &
                  & temp_char,                   &
                  & output%forces%forces(:,i)
  enddo
  
  close(qe_file)
end function

function read_qe_output_file_string(filename) result(output)
  use string_module
  implicit none
  
  type(String), intent(in) :: filename
  type(QeOutputFile)       :: output
  
  output = read_qe_output_file(char(filename))
end function

subroutine new_QeOutputForces(this, no_atoms)
  implicit none
  
  type(QeOutputForces), intent(out) :: this
  integer,              intent(in)  :: no_atoms
  
  allocate(this%atom_ids(no_atoms))
  allocate(this%atom_types(no_atoms))
  allocate(this%forces(3,no_atoms))
end subroutine

subroutine drop_QeOutputForces(this)
  implicit none
  
  type(QeOutputForces), intent(inout) :: this
  
  deallocate(this%atom_ids)
  deallocate(this%atom_types)
  deallocate(this%forces)
end subroutine

subroutine drop_QeOutputFile(this)
  implicit none
  
  type(QeOutputFile), intent(inout) :: this
  
  call drop(this%forces)
end subroutine

end module
