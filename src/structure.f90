! ======================================================================
! A class representing the contents of structure.dat
! ======================================================================
module structure_module
  use constants, only : dp
  use supercell_module
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
    ! Symmetry data
    integer                   :: no_symmetries
    real(dp),     allocatable :: rotation_matrices(:,:,:)
    real(dp),     allocatable :: offsets(:,:)
    real(dp),     allocatable :: offsets_cart(:,:)
    ! Superell data
    type(SupercellData)       :: supercell
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

subroutine new_StructureData(this,no_atoms,no_symmetries,sc_size)
  implicit none
  
  type(StructureData), intent(out) :: this
  integer,             intent(in)  :: no_atoms
  integer,             intent(in)  :: no_symmetries
  integer,             intent(in)  :: sc_size
  
  this%no_atoms = no_atoms
  this%no_modes = no_atoms*3
  allocate(this%species(no_atoms))
  allocate(this%mass(no_atoms))
  allocate(this%atoms(3,no_atoms))
  
  this%no_symmetries = no_symmetries
  if (no_symmetries /= 0) then
    allocate(this%rotation_matrices(3,3,no_symmetries))
    allocate(this%offsets(3,no_symmetries))
    allocate(this%offsets_cart(3,no_symmetries))
  endif
  
  call new(this%supercell,sc_size)
end subroutine

! Deallocates a Structure
subroutine drop_StructureData(this)
  implicit none
  
  type(StructureData), intent(inout) :: this
  
  deallocate(this%species)
  deallocate(this%mass)
  deallocate(this%atoms)
  
  if (this%no_symmetries /= 0) then
    deallocate(this%rotation_matrices)
    deallocate(this%offsets)
    deallocate(this%offsets_cart)
  endif
  
  call drop(this%supercell)
end subroutine

! reads structure.dat
function read_structure_file_character(filename,supercell) result(this)
  use file_module
  use linear_algebra, only : invert, invert_int
  use string_module
  implicit none
  
  character(*),        intent(in) :: filename
  type(SupercellData), intent(in) :: supercell
  type(StructureData)             :: this
  
  type(String), allocatable :: structure_file(:)
  type(String), allocatable :: line(:)
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
  structure_file = read_lines(filename)
  do i=1,size(structure_file)
    line = split(lower_case(structure_file(i)))
    if (line(1)=="lattice") then
      lattice_line = i
    elseif (line(1)=="atoms") then
      atoms_line = i
    elseif (line(1)=="symmetry") then
      symmetry_line = i
    elseif (line(1)=="end") then
      end_line = i
    endif
  enddo
  
  ! check layout is consistent with input file
  if (lattice_line/=1) then
    call print_line("Error: line 1 of "//filename//" is not 'Lattice'")
    stop
  elseif (atoms_line/=5) then
    call print_line("Error: line 5 of "//filename//" is not 'Atoms'")
    stop
  elseif (end_line/=size(structure_file)) then
    call print_line("Error: the last line of "//filename//" is not 'End'")
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
  
  call new(this,no_atoms,no_symmetries,supercell%sc_size)
  
  ! Store supercell data
  ! TODO: change input file format to include supercell data
  this%supercell = supercell
  
  ! read file into arrays
  do i=1,3
    this%lattice(i,:) = dble(split(structure_file(lattice_line+i)))
  enddo
  
  do i=1,this%no_atoms
    line = split(structure_file(atoms_line+i))
    this%species(i) = char(line(1))
    this%mass(i) = dble(line(2))
    this%atoms(:,i) = dble(line(3:5))
  enddo
  
  do i=1,this%no_symmetries
    do j=1,3
      line = split(structure_file(symmetry_line+(i-1)*4+j))
      this%rotation_matrices(j,:,i) = dble(line)
    enddo
    line = split(structure_file(symmetry_line+(i-1)*4+4))
    this%offsets(:,i) = dble(line)
  enddo
  
  ! calculate derived quantities
  call calculate_derived_atom_quantities(this)
  
  if (this%no_symmetries>0) then
    call calculate_derived_symmetry_quantities(this)
  endif
end function

function read_structure_file_string(filename,supercell) result(this)
  use string_module
  implicit none
  
  type(String),        intent(in) :: filename
  type(SupercellData), intent(in) :: supercell
  type(StructureData)             :: this
  
  this = read_structure_file(char(filename),supercell)
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
  
  call print_line(structure_file,'Lattice')
  do i=1,3
    call print_line(structure_file, this%lattice(i,:))
  enddo
  
  call print_line(structure_file,'Atoms')
  do i=1,this%no_atoms
    call print_line(structure_file, this%species(i)//' '// &
                                  & this%mass(i)//' '//    &
                                  & this%atoms(:,i))
  enddo
  
  if (this%no_symmetries/=0) then
    call print_line(structure_file,'Symmetry')
    do i=1,this%no_symmetries
      do j=1,3
        call print_line(structure_file, this%rotation_matrices(j,:,i))
      enddo
      call print_line(structure_file, this%offsets(:,i))
    enddo
  endif
  
  call print_line(structure_file,'End')
  
  close(structure_file)
end subroutine

subroutine write_structure_file_string(this,filename)
  use string_module
  implicit none
  
  type(StructureData), intent(in) :: this
  type(String),        intent(in) :: filename
  
  call write_structure_file(this,char(filename))
end subroutine

! ----------------------------------------------------------------------
! Reads the symmetry operations from the output of cellsym
! ----------------------------------------------------------------------
subroutine read_symmetry_file_character(this,filename)
  use string_module
  use file_module
  implicit none
  
  type(StructureData), intent(inout) :: this
  character(*),        intent(in)    :: filename
  
  ! File contents
  type(String), allocatable :: structure_file(:)
  type(String), allocatable :: line(:)
  
  ! line numbers
  integer :: start_line
  integer :: end_line
  
  ! working variables
  integer        :: i, j
  
  structure_file = read_lines(filename)
  
  ! Work out line numbers
  do i=1,size(structure_file)
    line = split(lower_case(structure_file(i)))
    if (size(line)==2) then
      if (line(1)=="%block" .and. line(2)=="symmetry_ops") then
        start_line = i
      elseif (line(1)=="%endblock" .and. line(2)=="symmetry_ops") then
        end_line = i
      endif
    endif
  enddo
  
  this%no_symmetries = (end_line-start_line)/5
  allocate(this%rotation_matrices(3,3,this%no_symmetries))
  allocate(this%offsets(3,this%no_symmetries))
  allocate(this%offsets_cart(3,this%no_symmetries))
  
  ! Read in symmetries
  do i=1,this%no_symmetries
    do j=1,3
      line = split(structure_file(start_line+(i-1)*5+j))
      this%rotation_matrices(j,:,i) = dble(line)
    enddo
    line = split(structure_file(start_line+(i-1)*5+4))
    this%offsets(:,i) = dble(line)
  enddo
  
  call calculate_derived_symmetry_quantities(this)
end subroutine

subroutine read_symmetry_file_String(this,filename)
  use string_module
  implicit none
  
  type(StructureData), intent(inout) :: this
  type(String),        intent(in)    :: filename
  
  call read_symmetry_file(this,char(filename))
end subroutine

! calculate derived quantities relating to lattice and atoms
subroutine calculate_derived_atom_quantities(this)
  use linear_algebra, only : invert, determinant
  implicit none
  
  type(StructureData), intent(inout) :: this
  
  this%recip_lattice = transpose(invert(this%lattice))
  this%volume        = abs(determinant(this%lattice))
end subroutine

! calculate derived quantities relating to symmetries
subroutine calculate_derived_symmetry_quantities(this)
  implicit none
  
  type(StructureData), intent(inout) :: this
  
  this%offsets_cart = matmul(this%lattice,this%offsets)
end subroutine
end module
