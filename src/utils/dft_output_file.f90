! A class for holding the information in a castep .castep or qe .out file
module dft_output_file_module
  use constants, only : dp
  use string_module
  implicit none
  
  private
  
  public :: DftOutputFile
  public :: read_dft_output_file
  public :: new
  public :: drop
  
  type DftOutputFile
    integer                   :: no_atoms
    type(String), allocatable :: species(:)
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

subroutine new_DftOutputFile(this, no_atoms)
  implicit none
  
  type(DftOutputFile), intent(out) :: this
  integer,             intent(in)  :: no_atoms
  
  this%no_atoms = no_atoms
  allocate(this%species(no_atoms))
  allocate(this%forces(3,no_atoms))
end subroutine

subroutine drop_DftOutputFile(this)
  implicit none
  
  type(DftOutputFile), intent(inout) :: this
  
  deallocate(this%species)
  deallocate(this%forces)
end subroutine

function read_castep_output_file(filename) result(output)
  use string_module
  use file_module
  implicit none
  
  type(String), intent(in) :: filename
  type(DftOutputFile)      :: output
  
  ! file contents
  type(String), allocatable :: castep_file(:)
  type(String), allocatable :: line(:)
  
  ! line numbers
  integer        :: energy_line
  integer        :: forces_start_line
  integer        :: forces_end_line
  
  ! temporary variables
  integer        :: i
  
  castep_file = read_lines(filename)
  
  ! Work out line numbers
  energy_line = 0
  forces_start_line = 0
  forces_end_line = 0
  do i=1,size(castep_file)
    line = split(lower_case(castep_file(i)))
    ! energy
    if (size(line)>=2) then
      if (line(1)=="final" .and. line(2)=="energy,") then
        energy_line = i
      endif
    endif
    ! forces
    if (size(line)>=2) then
      if (line(1)=="***********************" .and. line(2)=="forces") then
        forces_start_line = i
      endif
    endif
    if (size(line)>=1) then
      if(forces_start_line/=0 .and. &
       &line(1)=="******************************************************") then
        forces_end_line = i
      endif
    endif
  enddo
  
  if (energy_line==0) then
    call print_line("Error: Energy not found in "//char(filename))
  endif
  if (forces_start_line==0) then
    call print_line("Error: Start of forces not found in "//char(filename))
  endif
  if (forces_end_line==0) then
    call print_line("Error: End of forces not found in "//char(filename))
  endif
  
  ! Allocate output
  call new(output,forces_end_line-forces_start_line-7)
  
  ! Read data
  line = split(castep_file(energy_line))
  output%energy = dble(line(5))
  
  do i=1,output%no_atoms
    line = split(castep_file(forces_start_line+5+i))
    output%species(i) = line(2)
    output%forces(:,i) = dble(line(4:6))
  enddo
end function

function read_qe_output_file(filename) result(output)
  use constants, only : ev_per_rydberg,angstrom_per_bohr
  use string_module
  use file_module
  implicit none
  
  type(String), intent(in) :: filename
  type(DftOutputFile)      :: output
  
  ! file contents
  type(String), allocatable :: qe_file(:)
  type(String), allocatable :: line(:)
  
  ! line numbers
  integer        :: species_start_line
  integer        :: species_end_line
  integer        :: energy_line
  integer        :: forces_start_line
  integer        :: forces_end_line
  
  ! qe "type" to species conversion
  integer                   :: no_species
  type(String), allocatable :: species(:)
  
  ! temporary variables
  integer        :: i
  
  qe_file = read_lines(filename)
  
  ! Work out line numbers
  species_start_line = 0
  do i=1,size(qe_file)
    line = split(lower_case(qe_file(i)))
    ! species
    if (line(1)=="atomic" .and. line(2)=="species") then
      species_start_line = i
    elseif ( species_start_line/=0 .and. &
           & species_end_line==0   .and. &
           & size(line)==0) then
      species_end_line = i
    ! energy
    elseif (line(1)=="!") then
      energy_line=i
    ! forces
    elseif (line(1)=="forces" .and. line(2)=="acting") then
      forces_start_line = i
    elseif (line(1)=="total" .and. line(2)=="force") then
      forces_end_line = i
    endif
  enddo
  
  no_species = species_end_line-species_start_line-1
  allocate(species(no_species))
  
  ! Allocate output
  allocate(output%species(forces_end_line-forces_start_line-3))
  allocate(output%forces(3,forces_end_line-forces_start_line-3))
  
  ! Read data
  do i=1,species_end_line-species_start_line-1
    line = split(qe_file(species_start_line+1))
    species(i) = line(1)
  enddo
  
  line = split(qe_file(energy_line))
  output%energy = dble(line(5))
  
  do i=1,forces_end_line-forces_start_line-3
    line = split(qe_file(forces_start_line+1+i))
    output%species(i) = species(int(line(4)))
    output%forces(:,i) = dble(line(7:9))
  enddo
  
  output%forces = output%forces*ev_per_rydberg/angstrom_per_bohr
end function

function read_dft_output_file(dft_code,filename) result(output)
  use string_module
  implicit none
  
  type(String), intent(in) :: dft_code
  type(String), intent(in) :: filename
  type(DftOutputFile)      :: output
  
  if (dft_code=="castep") then
    output = read_castep_output_file(filename)
  elseif (dft_code=="qe") then
    output = read_qe_output_file(filename)
  endif
end function
end module
