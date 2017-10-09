! ======================================================================
! A class for holding the information in a castep .castep or qe .out file.
! ======================================================================
module dft_output_file_module
  use constants_module, only : dp
  use string_module
  use io_module
  use linear_algebra_module
  implicit none
  
  private
  
  public :: make_dft_output_filename
  public :: DftOutputFile
  public :: read_dft_output_file
  public :: new
  
  type DftOutputFile
    integer                       :: no_atoms
    type(String), allocatable     :: species(:)
    real(dp)                      :: energy
    type(RealVector), allocatable :: forces(:)
  end type
  
  interface new
    module procedure new_DftOutputFile
  end interface
  
contains

! ----------------------------------------------------------------------
! Converts a file seedname into the appropriate dft input or output filename.
! ----------------------------------------------------------------------
function make_dft_output_filename(dft_code,seedname) result(output)
  implicit none
  
  type(String), intent(in) :: dft_code
  type(String), intent(in) :: seedname
  type(String)             :: output
  
  if (dft_code == 'castep') then
    output = seedname//'.castep'
  elseif (dft_code == 'vasp') then
    output = 'OUTCAR'
  elseif (dft_code == 'qe') then
    output = seedname//'.out'
  else
    call print_line('Unrecognised dft code: '//dft_code)
    call err()
  endif
end function

subroutine new_DftOutputFile(this, no_atoms)
  implicit none
  
  type(DftOutputFile), intent(out) :: this
  integer,             intent(in)  :: no_atoms
  
  this%no_atoms = no_atoms
  allocate(this%species(no_atoms))
  allocate(this%forces(no_atoms))
end subroutine

function read_castep_output_file(filename) result(output)
  use constants_module, only : angstrom_per_bohr, ev_per_hartree
  use ifile_module
  implicit none
  
  type(String), intent(in) :: filename
  type(DftOutputFile)      :: output
  
  ! file contents
  type(IFile)               :: castep_file
  type(String), allocatable :: line(:)
  
  ! line numbers
  integer        :: energy_line
  integer        :: forces_start_line
  integer        :: forces_end_line
  
  ! temporary variables
  integer        :: i
  
  castep_file = filename
  
  ! Work out line numbers
  energy_line = 0
  forces_start_line = 0
  forces_end_line = 0
  do i=1,size(castep_file)
    line = split(lower_case(castep_file%line(i)))
    ! energy
    if (size(line)>=2) then
      if (line(1)=='final' .and. line(2)=='energy,') then
        energy_line = i
      endif
    endif
    ! forces
    if (size(line)>=2) then
      if ( (line(1)=='***********************' .and. line(2)=='forces') .or. &
         & (line(1)=='*****************' .and. line(2)=='symmetrised')  .or. &
         & (line(1)=='*******************' .and. line(2)=='unconstrained')   &
         & ) then
        forces_start_line = i
      endif
    endif
    if (size(line)==1) then
      if ( forces_start_line/=0 .and. &
         & forces_end_line==0   .and. &
         & len(line(1))>=5) then
        if (slice(line(1),1,5)=='*****') then
          forces_end_line = i
        endif
      endif
    endif
  enddo
  
  if (energy_line==0) then
    call print_line('Error: Energy not found in '//char(filename))
    stop
  endif
  if (forces_start_line==0) then
    call print_line('Error: Start of forces not found in '//char(filename))
    stop
  endif
  if (forces_end_line==0) then
    call print_line('Error: End of forces not found in '//char(filename))
    stop
  endif
  
  ! Allocate output
  call new(output,forces_end_line-forces_start_line-7)
  
  ! Read data
  line = split(castep_file%line(energy_line))
  output%energy = dble(line(5)) / ev_per_hartree
  
  do i=1,output%no_atoms
    line = split(castep_file%line(forces_start_line+5+i))
    output%species(i) = line(2)
    output%forces(i) = dble(line(4:6)) * angstrom_per_bohr / ev_per_hartree
  enddo
end function

function read_qe_output_file(filename) result(output)
  use constants_module, only : ev_per_rydberg, ev_per_hartree
  use ifile_module
  implicit none
  
  type(String), intent(in) :: filename
  type(DftOutputFile)      :: output
  
  ! file contents
  type(IFile)               :: qe_file
  type(String), allocatable :: line(:)
  
  ! line numbers
  integer :: species_start_line
  integer :: species_end_line
  integer :: energy_line
  integer :: forces_start_line
  integer :: forces_end_line
  
  ! qe 'type' to species conversion
  integer                   :: no_species
  type(String), allocatable :: species(:)
  
  ! temporary variables
  integer        :: i
  
  qe_file = filename
  
  ! Work out line numbers
  species_start_line = 0
  do i=1,size(qe_file)
    line = split(lower_case(qe_file%line(i)))
    ! species
    if (line(1)=='atomic' .and. line(2)=='species') then
      species_start_line = i
    elseif ( species_start_line/=0 .and. &
           & species_end_line==0   .and. &
           & size(line)==0) then
      species_end_line = i
    ! energy
    elseif (line(1)=='!') then
      energy_line=i
    ! forces
    elseif (line(1)=='forces' .and. line(2)=='acting') then
      forces_start_line = i
    elseif (line(1)=='total' .and. line(2)=='force') then
      forces_end_line = i
    endif
  enddo
  
  no_species = species_end_line-species_start_line-1
  allocate(species(no_species))
  
  ! Allocate output
  allocate(output%species(forces_end_line-forces_start_line-3))
  allocate(output%forces(forces_end_line-forces_start_line-3))
  
  ! Read data
  do i=1,species_end_line-species_start_line-1
    line = split(qe_file%line(species_start_line+1))
    species(i) = line(1)
  enddo
  
  line = split(qe_file%line(energy_line))
  output%energy = dble(line(5)) * ev_per_rydberg / ev_per_hartree
  
  do i=1,forces_end_line-forces_start_line-3
    line = split(qe_file%line(forces_start_line+1+i))
    output%species(i) = species(int(line(4)))
    output%forces(i) = dble(line(7:9)) * ev_per_rydberg / ev_per_hartree
  enddo
end function

function read_dft_output_file(dft_code,filename) result(output)
  implicit none
  
  type(String), intent(in) :: dft_code
  type(String), intent(in) :: filename
  type(DftOutputFile)      :: output
  
  if (dft_code=='castep') then
    output = read_castep_output_file(filename)
  elseif (dft_code=='qe') then
    output = read_qe_output_file(filename)
  endif
end function
end module
