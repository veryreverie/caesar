! A class for holding the information in a castep .castep or qe .out file
module dft_output_file_module
  use constants, only : dp
  implicit none
  
  private
  
  public :: DftOutputFile
  public :: read_castep_output_file
  public :: read_qe_output_file
  public :: new
  public :: drop
  
  type DftOutputFile
    integer                   :: no_atoms
    character(2), allocatable :: species(:)
    real(dp),     allocatable :: forces(:,:)
  end type
  
  interface read_castep_output_file
    module procedure read_castep_output_file_character
    module procedure read_castep_output_file_string
  end interface
  
  interface read_qe_output_file
    module procedure read_qe_output_file_character
    module procedure read_qe_output_file_string
  end interface
  
  interface new
    module procedure new_DftOutputFile
  end interface
  
  interface drop
    module procedure drop_DftOutputFile
  end interface
    
contains

function read_castep_output_file_character(filename) result(output)
  use file_io
  use utils, only : lower_case
  implicit none
  
  character(*), intent(in) :: filename
  type(DftOutputFile)      :: output
  
  ! file unit
  integer        :: castep_file
  
  ! line numbers
  integer        :: file_length
  integer        :: forces_start_line
  integer        :: forces_end_line
  
  ! temporary variables
  integer        :: i
  character(100) :: line
  character(100) :: temp_char
  
  castep_file = open_read_file(filename)
  file_length = count_lines(castep_file)
  
  forces_start_line = 0
  do i=1,file_length
    read(castep_file,"(a)") line
    line = lower_case(trim(line))
    if (line=="*********************** Forces ***********************") then
      forces_start_line = i
    elseif (forces_start_line/=0 .and. &
       & line=="******************************************************") then
      forces_end_line = i
    endif
  enddo
  
  call new(output,forces_end_line-forces_start_line-7)
  
  rewind(castep_file)
  
  do i=1,forces_start_line+5
    read(castep_file,*)
  enddo
  do i=1,output%no_atoms
    read(castep_file,*) temp_char,         &
                      & output%species(i), &
                      & temp_char,         &
                      & output%forces(:,i)
  enddo
  
  close(castep_file)
end function

function read_castep_output_file_string(filename) result(output)
  use string_module
  implicit none
  
  type(String), intent(in) :: filename
  type(DftOutputFile)      :: output
  
  output = read_castep_output_file(char(filename))
end function

function read_qe_output_file_character(filename) result(output)
  use file_io
  use utils, only : lower_case
  implicit none
  
  character(*), intent(in) :: filename
  type(DftOutputFile)      :: output
  
  ! file unit
  integer        :: qe_file
  
  ! line numbers
  integer        :: file_length
  integer        :: species_start_line
  integer        :: species_end_line
  integer        :: forces_start_line
  integer        :: forces_end_line
  
  ! qe "type" to species conversion
  integer                   :: no_species
  character(2), allocatable :: species(:)
  integer                   :: species_type
  
  ! temporary variables
  integer        :: i
  character(100) :: line
  character(100) :: temp_char
  
  qe_file = open_read_file(filename)
  file_length = count_lines(qe_file)
  
  species_start_line = 0
  do i=1,file_length
    read(qe_file,"(a)") line
    line = lower_case(trim(line))
    if (line(1:14)=="atomic species") then
      species_start_line = i
    elseif (species_start_line/=0 .and. line=="") then
      species_end_line = i
    elseif (line(1:13)=="forces acting") then
      forces_start_line = i
    elseif (line(1:11)=="total force") then
      forces_end_line = i
    endif
  enddo
  
  no_species = species_end_line-species_start_line-1
  allocate(species(no_species))
  
  call new(output,forces_end_line-forces_start_line-3)
  
  rewind(qe_file)
  
  do i=1,species_start_line
    read(qe_file,*)
  enddo
  do i=1,no_species
    read(qe_file,*) species(i)
  enddo
  do i=species_end_line,forces_start_line+1
    read(qe_file,*)
  enddo
  do i=1,output%no_atoms
    read(qe_file,*) temp_char,    &
                  & temp_char,    &
                  & temp_char,    &
                  & species_type, &
                  & temp_char,    &
                  & temp_char,    &
                  & output%forces(:,i)
    output%species(i) = species(species_type)
  enddo
  
  close(qe_file)
end function

function read_qe_output_file_string(filename) result(output)
  use string_module
  implicit none
  
  type(String), intent(in) :: filename
  type(DftOutputFile)      :: output
  
  output = read_qe_output_file(char(filename))
end function

subroutine new_DftOutputFile(this, no_atoms)
  implicit none
  
  type(DftOutputFile), intent(out) :: this
  integer,             intent(in)  :: no_atoms
  
  allocate(this%species(no_atoms))
  allocate(this%forces(3,no_atoms))
end subroutine

subroutine drop_DftOutputFile(this)
  implicit none
  
  type(DftOutputFile), intent(inout) :: this
  
  deallocate(this%species)
  deallocate(this%forces)
end subroutine
end module
