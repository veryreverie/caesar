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
    deallocate(this%rotation_matrices)
    deallocate(this%offsets)
    deallocate(this%offsets_cart)
  endif
end subroutine

! reads structure.dat
function read_structure_file_character(filename) result(this)
  use file_module
  use linear_algebra, only : inv_33, determinant33
  use string_module
  implicit none
  
  character(*), intent(in) :: filename
  type(StructureData)      :: this
  
  type(String), allocatable :: contents(:)
  type(String)              :: line
  character(100)            :: line_char
  integer                   :: i,j
  integer                   :: no_atoms
  integer                   :: no_symmetries
  
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
  contents = read_lines(filename)
  do i=1,size(contents)
    line = lower_case(contents(i))
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
  elseif (atoms_line/=5) then
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
  do i=1,3
    line_char = char(contents(lattice_line+i))
    read(line_char,*) this%lattice(i,:)
  enddo
  
  do i=1,this%no_atoms
    line_char = char(contents(atoms_line+i))
    read(line_char,*) this%species(i), this%mass(i), this%atoms(:,i)
  enddo
  
  do i=1,this%no_symmetries*4
    do j=1,3
      line_char = char(contents(symmetry_line+(i-1)*4+j))
      read(line_char,*) this%rotation_matrices(j,:,i)
    enddo
    line_char = char(contents(symmetry_line+(i-1)*4+4))
    read(line_char,*) this%offsets(:,i)
  enddo
  
  ! calculate derived quantities
  call calculate_derived_atom_quantities(this)
  
  if (this%no_symmetries>0) then
    call calculate_derived_symmetry_quantities(this)
  endif
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
  integer :: i,j
  
  structure_file = open_write_file(filename)
  write(structure_file,"(a)") "Lattice"
  do i=1,3
    write(structure_file,*) this%lattice(i,:)
  enddo
  write(structure_file,"(a)") "Atoms"
  do i=1,this%no_atoms
    write(structure_file,*) this%species(i), this%mass(i), this%atoms(:,i)
  enddo
  if (this%no_symmetries/=0) then
    write(structure_file,"(a)") "Symmetry"
  endif
  do i=1,4*this%no_symmetries
    do j=1,3
      write(structure_file,*) this%rotation_matrices(j,:,i)
    enddo
    write(structure_file,*) this%offsets(:,i)
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
  
  ! File contents
  type(String), allocatable :: contents(:)
  
  ! line numbers
  integer :: file_length
  integer :: start_line
  integer :: end_line
  
  ! working variables
  integer        :: i, j
  character(100) :: line_char
  
  contents = read_lines(filename)
  file_length = count_lines(filename)
  
  ! Work out line numbers
  do i=1,size(contents)
    line_char = char(lower_case(contents(i)))
    if (line_char(1:19)=="%block symmetry_ops") then
      start_line = i
    elseif (line_char(1:22)=="%endblock symmetry_ops") then
      end_line = i
    endif
  enddo
  
  this%no_symmetries = (end_line-start_line)/5
  allocate(this%rotation_matrices(3,3,this%no_symmetries))
  allocate(this%offsets(3,this%no_symmetries))
  
  ! Read in symmetries
  do i=1,this%no_symmetries
    do j=1,3
      line_char = char(contents(start_line+(i-1)*5+j))
      read(line_char,*) this%rotation_matrices(j,:,i)
    enddo
    line_char = char(contents(start_line+(i-1)*5+4))
    read(line_char,*) this%offsets(:,i)
  enddo
  
  call calculate_derived_symmetry_quantities(this)
end subroutine

subroutine read_symmetry_file_String(this,filename)
  implicit none
  
  type(StructureData), intent(inout) :: this
  type(String),        intent(in)    :: filename
  
  call read_symmetry_file(this,char(filename))
end subroutine

! calculate derived quantities relating to lattice and atoms
subroutine calculate_derived_atom_quantities(this)
  use linear_algebra, only : inv_33, determinant33
  implicit none
  
  type(StructureData), intent(inout) :: this
  
  this%recip_lattice = transpose(inv_33(this%lattice))
  this%volume        = abs(determinant33(this%lattice))
  this%frac_atoms    = matmul(this%recip_lattice,this%atoms)
  this%cart_atoms    = matmul(this%lattice,this%frac_atoms)
end subroutine

! calculate derived quantities relating to symmetries
subroutine calculate_derived_symmetry_quantities(this)
  implicit none
  
  type(StructureData), intent(inout) :: this
  
  this%offsets_cart = matmul(this%lattice,this%offsets)
end subroutine
end module
