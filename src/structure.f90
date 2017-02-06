! ======================================================================
! A class representing the contents of structure.dat
! ======================================================================
module structure_module
  use constants, only : dp
  use string_module
  implicit none
  
  ! the structure class
  type StructureData
    ! Lattice data
    real(dp)                  :: lattice(3,3)
    real(dp)                  :: recip_lattice(3,3)
    real(dp)                  :: volume
    ! Atom data
    integer                   :: no_atoms
    integer                   :: no_modes
    character(2), allocatable :: species(:)
    real(dp),     allocatable :: mass(:)
    real(dp),     allocatable :: atoms(:,:)
    real(dp),     allocatable :: frac_atoms(:,:)
    real(dp),     allocatable :: cart_atoms(:,:)
    ! Symmetry data
    integer                   :: no_symmetries
    real(dp),     allocatable :: symmetries(:,:)
    real(dp),     allocatable :: rotation_matrices(:,:,:)
    real(dp),     allocatable :: offsets(:,:)
    real(dp),     allocatable :: offsets_cart(:,:)
  end type
  
  interface new
    module procedure new_StructureData
  end interface
  
  interface drop
    module procedure drop_StructureData
  end interface
  
  interface read_structure_file
    module procedure read_structure_file_character
    module procedure read_structure_file_string
  end interface
  
  interface write_structure_file
    module procedure write_structure_file_character
    module procedure write_structure_file_string
  end interface
  
  interface read_symmetry_file
    module procedure read_symmetry_file_character
    module procedure read_symmetry_file_string
  end interface

contains

subroutine new_StructureData(this,no_atoms,no_symmetries)
  implicit none
  
  type(StructureData) :: this
  integer             :: no_atoms
  integer             :: no_symmetries
  
  this%no_atoms = no_atoms
  this%no_modes = no_atoms*3
  allocate(this%species(no_atoms))
  allocate(this%mass(no_atoms))
  allocate(this%atoms(3,no_atoms))
  allocate(this%frac_atoms(3,no_atoms))
  allocate(this%cart_atoms(3,no_atoms))
  
  this%no_symmetries = no_symmetries
  if (no_symmetries /= 0) then
    allocate(this%symmetries(3,no_symmetries*4))
    allocate(this%rotation_matrices(3,3,no_symmetries))
    allocate(this%offsets(3,no_symmetries))
    allocate(this%offsets_cart(3,no_symmetries))
  endif
end subroutine

! Deallocates a Structure
subroutine drop_StructureData(this)
  implicit none
  
  type(StructureData), intent(inout) :: this
  
  deallocate(this%species)
  deallocate(this%mass)
  deallocate(this%atoms)
  deallocate(this%frac_atoms)
  deallocate(this%cart_atoms)
  
  if (this%no_symmetries /= 0) then
    deallocate(this%symmetries)
    deallocate(this%rotation_matrices)
    deallocate(this%offsets)
    deallocate(this%offsets_cart)
  endif
end subroutine

! reads structure.dat
function read_structure_file_character(filename) result(this)
  use utils,          only : lower_case
  use file_module
  use linear_algebra, only : inv_33, determinant33
  implicit none
  
  character(*), intent(in) :: filename
  type(StructureData)      :: this
  
  integer        :: file_length
  integer        :: file_unit
  character(100) :: line
  integer        :: i, j
  integer        :: no_atoms
  integer        :: no_symmetries
  
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
  file_length = count_lines(filename)
  file_unit = open_read_file(filename)
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
  if (lattice_line/=1) then
    write(*,*) "Line 1 of structure.dat is not 'Lattice'"
    stop
  elseif (atoms_line/=4) then
    write(*,*) "Line 5 of structure.dat is not 'Atoms'"
    stop
  elseif (end_line==0) then
    write(*,*) "Structure.dat does not contain the line 'End'"
    stop
  endif
  
  ! set counts, and allocate structure
  if (symmetry_line==0) then
    ! structure.dat does not contain symmetries
    no_atoms = end_line-atoms_line-1
    no_symmetries = 0
  else
    no_atoms = symmetry_line-atoms_line-1
    no_symmetries = (end_line-symmetry_line-1)/4
  endif
  call new(this,no_atoms,no_symmetries)
  
  ! read file into arrays
  rewind(file_unit)
  
  do i=1,file_length
    read(file_unit,"(a)") line
    if (i>lattice_line .and. i<=lattice_line+3) then
      ! Read in lattice
      j=i-lattice_line
      read(line,*) this%lattice(j,:)
    elseif (i>atoms_line .and. i<=atoms_line+this%no_atoms) then
      ! Read in atoms
      j=i-atoms_line
      read(line,*) this%species(j), this%mass(j), this%atoms(:,j)
    elseif (i>symmetry_line .and. i<= symmetry_line+this%no_symmetries) then
      ! Read in symmetries
      j=i-symmetry_line
      read(line,*) this%symmetries(:,i)
    endif
  enddo
  
  close(file_unit)
  
  ! calculate derived quantities
  this%recip_lattice = transpose(inv_33(this%lattice))
  this%volume        = abs(determinant33(this%lattice))
  this%frac_atoms    = matmul(this%recip_lattice,this%atoms)
  this%cart_atoms    = matmul(this%lattice,this%frac_atoms)
  
  do i=1,no_symmetries
    this%rotation_matrices(1,:,i) = this%symmetries(:,4*i-3)
    this%rotation_matrices(2,:,i) = this%symmetries(:,4*i-2)
    this%rotation_matrices(3,:,i) = this%symmetries(:,4*i-1)
    this%offsets(:,i)             = this%symmetries(:,4*i)
  enddo
  
  this%offsets_cart = matmul(this%lattice,this%offsets)
end function

function read_structure_file_string(filename) result(this)
  implicit none
  
  type(String), intent(in) :: filename
  type(StructureData)      :: this
  
  this = read_structure_file(char(filename))
end function

subroutine write_structure_file_character(this,filename)
  use string_module
  use file_module
  implicit none
  
  type(StructureData), intent(in) :: this
  character(*),        intent(in) :: filename
  
  integer :: structure_file
  integer :: i
  
  structure_file = open_write_file(filename)
  write(structure_file,"(a)") "Lattice"
  do i=1,3
    write(structure_file,*) this%lattice(i,:)
  enddo
  write(structure_file,"(a)") "Atoms"
  do i=1,this%no_atoms
    write(structure_file,*) this%species(i), &
                          & this%mass(i),    &
                          & this%atoms(:,i)
  enddo
  if (this%no_symmetries/=0) then
    write(structure_file,"(a)") "Symmetry"
  endif
  do i=1,4*this%no_symmetries
    write(structure_file,*) this%symmetries(:,i)
  enddo
  write(structure_file,"(a)") "End"
  close(structure_file)
end subroutine

subroutine write_structure_file_string(this,filename)
  implicit none
  
  type(StructureData), intent(in) :: this
  type(String),        intent(in) :: filename
  
  call write_structure_file(this,char(filename))
end subroutine

! ----------------------------------------------------------------------
! Reads the symmetry operations from the output of cellsym
! ----------------------------------------------------------------------
subroutine read_symmetry_file_character(this,filename)
  use utils, only : lower_case
  use file_module
  implicit none
  
  type(StructureData), intent(inout) :: this
  character(*),        intent(in)    :: filename
  
  ! line numbers
  integer :: file_length
  integer :: start_line
  integer :: end_line
  
  ! file unit
  integer :: symmetry_file
  
  ! working variables
  integer        :: i, j
  character(100) :: line
  
  file_length = count_lines(filename)
  symmetry_file = open_read_file(filename)
  
  ! Work out line numbers
  do i=1,file_length
    read(symmetry_file,"(a)") line
    line = lower_case(trim(line))
    if (line(1:19)=="%block symmetry_ops") then
      start_line = i
    elseif (line(1:22)=="%endblock symmetry_ops") then
      end_line = i
    endif
  enddo
  
  this%no_symmetries = (end_line-start_line)/5
  allocate(this%symmetries(3,4*this%no_symmetries))
  
  rewind(symmetry_file)
  
  ! Read in symmetries
  do i=1,file_length
    read(symmetry_file,"(a)") line
    if (i>start_line .and. i<end_line) then
      j = i-start_line
      if (modulo(j,5)/=0) then
        read(line,*) this%symmetries(:,(j/5)*4+modulo(j,5))
      endif
    endif
  enddo
  
  close(symmetry_file)
  
  do i=1,this%no_symmetries
    this%rotation_matrices(1,:,i) = this%symmetries(:,4*i-3)
    this%rotation_matrices(2,:,i) = this%symmetries(:,4*i-2)
    this%rotation_matrices(3,:,i) = this%symmetries(:,4*i-1)
    this%offsets(:,i)             = this%symmetries(:,4*i)
  enddo
end subroutine

subroutine read_symmetry_file_String(this,filename)
  implicit none
  
  type(StructureData), intent(inout) :: this
  type(String),        intent(in)    :: filename
  
  call read_symmetry_file(this,char(filename))
end subroutine
end module
