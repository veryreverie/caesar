! ======================================================================
! A class representing the contents of structure.dat
! ======================================================================
module structure_module
  use constants, only : dp
  use string_module
  implicit none
  
  ! the structure class
  type StructureData
    real(dp)                  :: lattice(3,3)             ! lattice.dat
    real(dp)                  :: recip_lattice(3,3)
    integer                   :: no_atoms                 ! equilibrium.dat
    character(2), allocatable :: species(:)               ! equilibrium.dat
    real(dp),     allocatable :: mass(:)                  ! equilibrium.dat
    real(dp),     allocatable :: atoms(:,:)               ! equilibrium.dat
    integer                   :: no_modes
    integer                   :: no_symmetries            ! symmetry.dat
    real(dp),     allocatable :: symmetries(:,:)          ! symmetry.dat
    real(dp),     allocatable :: rotation_matrices(:,:,:)
    real(dp),     allocatable :: offsets(:,:)
  end type
  
  interface read_structure_file
    module procedure read_structure_file_character
    module procedure read_structure_file_string
  end interface

  interface drop
    module procedure drop_structure
  end interface
  
contains

! reads structure.dat
function read_structure_file_character(filename) result(output)
  use utils,          only : lower_case
  use file_io,        only : open_read_file, count_lines
  use linear_algebra, only : inv_33
  implicit none
  
  character(*), intent(in) :: filename
  type(StructureData)      :: output
  
  integer        :: file_length
  integer        :: file_unit
  character(100) :: line
  integer        :: i
  
  ! line numbers
  integer :: lattice_line  ! The line "Lattice"
  integer :: atoms_line    ! The line "Atoms"
  integer :: symmetry_line ! The line "Symmetry"
  integer :: end_line      ! The line "End"
  
  ! initialise line numbers
  lattice_line = 0
  atoms_line = 0
  symmetry_line = 0
  end_line = 0
  
  ! work out layout of file
  file_unit = open_read_file(filename)
  file_length = count_lines(file_unit)
  do i=1,file_length
    read(file_unit,"(a)") line
    line = lower_case(trim(line))
    if (line=="lattice") then
      lattice_line = i
    elseif (line=="atoms") then
      atoms_line = i
    elseif (line=="symmetry") then
      symmetry_line = i
    elseif (line=="end") then
      end_line = i
    endif
  enddo
  
  ! check layout is consistent with input file
  if (lattice_line==0) then
    write(*,*) "structure.dat does not contain the line 'Lattice'"
    stop
  elseif (atoms_line==0) then
    write(*,*) "Structure.dat does not contain the line 'Atoms'"
    stop
  elseif (symmetry_line==0) then
    write(*,*) "Structure.dat does not contain the line 'Symmetry'"
    stop
  elseif (end_line==0) then
    write(*,*) "Structure.dat does not contain the line 'End'"
    stop
  elseif (lattice_line/=1) then
    write(*,*) "Line 1 of structure.dat is not 'Lattice'"
    stop
  elseif (atoms_line/=4) then
    write(*,*) "Line 5 of structure.dat is not 'Atoms'"
    stop
  endif
  
  ! set counts, and allocate arrays
  output%no_atoms = symmetry_line-atoms_line-1
  output%no_symmetries = (end_line-symmetry_line-1)/4
  allocate(output%species(output%no_atoms))
  allocate(output%mass(output%no_atoms))
  allocate(output%atoms(3,output%no_atoms))
  allocate(output%symmetries(3,output%no_symmetries*4))
  allocate(output%rotation_matrices(3,3,output%no_symmetries))
  allocate(output%offsets(3,output%no_symmetries))
  
  ! read file into arrays
  rewind(file_unit)
  read(file_unit,*) ! lattice_line
  do i=1,3
    read(file_unit,*) output%lattice(i,:)
  enddo
  read(file_unit,*) ! atoms_line
  do i=1,output%no_atoms
    read(file_unit,*) output%species(i),output%mass(i),output%atoms(:,i)
  enddo
  read(file_unit,*) ! symmetry_line
  do i=1,4*output%no_symmetries
    read(file_unit,*) output%symmetries(:,i)
  enddo
  
  close(file_unit)
  
  ! calculate derived quantities
  output%recip_lattice = transpose(inv_33(output%lattice))
  
  output%no_modes = output%no_atoms*3
  
  do i=1,output%no_symmetries
    output%rotation_matrices(1,:,i) = output%symmetries(:,4*i-3)
    output%rotation_matrices(2,:,i) = output%symmetries(:,4*i-2)
    output%rotation_matrices(3,:,i) = output%symmetries(:,4*i-1)
    output%offsets(:,i)             = output%symmetries(:,4*i)
  enddo
end function

function read_structure_file_string(filename) result(output)
  implicit none
  
  type(String), intent(in) :: filename
  type(StructureData)      :: output
  
  output = read_structure_file(char(filename))
end function

! Deallocates a Structure
subroutine drop_structure(this)
  implicit none
  
  type(StructureData), intent(inout) :: this
  
  deallocate(this%species)
  deallocate(this%atoms)
  deallocate(this%symmetries)
  deallocate(this%rotation_matrices)
end subroutine
end module
