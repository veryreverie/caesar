! A class for holding the information in a castep .castep or qe .out file
module dft_output_file_module
  use constants, only : dp
  implicit none
  
  private
  
  public :: DftOutputFile
  public :: read_dft_output_file
  public :: new
  public :: drop
  
  type DftOutputFile
    integer                   :: no_atoms
    character(2), allocatable :: species(:)
    real(dp)                  :: energy
    real(dp),     allocatable :: forces(:,:)
  end type
  
  interface new
    module procedure new_DftOutputFile
  end interface
  
  interface drop
    module procedure drop_DftOutputFile
  end interface
    
contains

function read_castep_output_file(filename) result(output)
  use utils, only : lower_case
  use string_module
  use file_module
  implicit none
  
  type(String), intent(in) :: filename
  type(DftOutputFile)      :: output
  
  ! file unit
  integer        :: castep_file
  
  ! line numbers
  integer        :: file_length
  integer        :: energy_line
  integer        :: forces_start_line
  integer        :: forces_end_line
  
  ! temporary variables
  integer        :: i,j
  character(100) :: line
  character(100) :: dump
  
  file_length = count_lines(filename)
  castep_file = open_read_file(filename)
  
  ! Work out line numbers
  forces_start_line = 0
  do i=1,file_length
    read(castep_file,"(a)") line
    line = lower_case(trim(line))
    ! energy
    if (line(1:12)=="final energy") then
      energy_line = i
    ! forces
    elseif (line=="*********************** Forces ***********************") then
      forces_start_line = i
    elseif (forces_start_line/=0 .and. &
          & line=="******************************************************") then
      forces_end_line = i
    endif
  enddo
  
  ! Allocate output
  call new(output,forces_end_line-forces_start_line-7)
  
  rewind(castep_file)
  
  ! Read data
  do i=1,file_length
    read(castep_file,"(a)") line
    ! energy
    if (i==energy_line) then
      read(line,*) dump,dump,dump,dump,output%energy
    ! forces
    elseif (i>forces_start_line+5 .and. i<forces_end_line) then
      j = i-(forces_start_line+5)
      read(line,*) dump, output%species(j), dump, output%forces(:,j)
    endif
  enddo
  
  close(castep_file)
end function

function read_qe_output_file(filename) result(output)
  use constants, only : Ry,bohr
  use utils,     only : lower_case
  use string_module
  use file_module
  implicit none
  
  type(String), intent(in) :: filename
  type(DftOutputFile)      :: output
  
  ! file unit
  integer        :: qe_file
  
  ! line numbers
  integer        :: file_length
  integer        :: species_start_line
  integer        :: species_end_line
  integer        :: energy_line
  integer        :: forces_start_line
  integer        :: forces_end_line
  
  ! qe "type" to species conversion
  integer                   :: no_species
  character(2), allocatable :: species(:)
  integer                   :: species_type
  
  ! temporary variables
  integer        :: i,j
  character(100) :: line
  character(100) :: dump
  
  file_length = count_lines(filename)
  qe_file = open_read_file(filename)
  
  ! Work out line numbers
  species_start_line = 0
  do i=1,file_length
    read(qe_file,"(a)") line
    line = lower_case(trim(line))
    ! species
    if (line(1:14)=="atomic species") then
      species_start_line = i
    elseif (species_start_line/=0 .and. line=="") then
      species_end_line = i
    ! energy
    elseif (line(1:1)=="!") then
      energy_line=i
    ! forces
    elseif (line(1:13)=="forces acting") then
      forces_start_line = i
    elseif (line(1:11)=="total force") then
      forces_end_line = i
    endif
  enddo
  
  no_species = species_end_line-species_start_line-1
  allocate(species(no_species))
  
  ! Allocate output
  allocate(output%species(forces_end_line-forces_start_line-3))
  allocate(output%forces(3,forces_end_line-forces_start_line-3))
  
  rewind(qe_file)
  
  ! Read data
  do i=1,file_length
    read(qe_file,"(a)") line
    ! species
    if (i>species_start_line .and. i<species_end_line) then
      j = i-species_start_line
      read(line,*) species(j)
    ! energy
    elseif (i==energy_line) then
      read(line,*) dump,dump,dump,dump, output%energy
    ! forces
    elseif (i>forces_start_line+1 .and. i<forces_end_line-1) then
      j = i-(forces_start_line+1)
      read(line,*) dump,dump,dump, species_type, dump,dump, output%forces(:,j)
      output%species(j) = species(species_type)
    endif
  enddo
  
  output%forces = output%forces*Ry/bohr
  
  close(qe_file)
end function

function read_dft_output_file(dft_code,dft_dir,seedname) result(output)
  use string_module
  implicit none
  
  type(String), intent(in) :: dft_code
  type(String), intent(in) :: dft_dir
  type(String), intent(in) :: seedname
  type(DftOutputFile)      :: output
  
  if (dft_code=="castep") then
    output = read_castep_output_file(dft_dir//'/'//seedname//'.castep')
  elseif (dft_code=="qe") then
    output = read_qe_output_file(dft_dir//'/'//seedname//'.out')
  endif
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
