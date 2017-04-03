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
    type(String), allocatable :: species(:)
    real(dp),     allocatable :: mass(:)
    real(dp),     allocatable :: atoms(:,:)
    ! Conversions between atom representations.
    !    [1...no_atoms] vs [1...no_atoms_in_prim]*[sc_size].
    integer,      allocatable :: atom_to_prim(:)
    integer,      allocatable :: atom_to_gvec(:)
    integer,      allocatable :: gvec_and_prim_to_atom(:,:)
    ! Symmetry data
    integer                   :: no_symmetries
    real(dp),     allocatable :: rotation_matrices(:,:,:)
    real(dp),     allocatable :: offsets(:,:)
    real(dp),     allocatable :: offsets_cart(:,:)
    ! Superell data
    integer                   :: sc_size
    integer                   :: supercell(3,3)
    integer                   :: recip_supercell(3,3)
    integer,      allocatable :: gvectors(:,:)
    integer,      allocatable :: paired_gvec(:)
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
  
  integer :: prim,gvec,atom
  
  this%no_atoms = no_atoms
  this%no_modes = no_atoms*3
  allocate(this%species(no_atoms))
  allocate(this%mass(no_atoms))
  allocate(this%atoms(3,no_atoms))
  
  allocate(this%atom_to_prim(no_atoms))
  allocate(this%atom_to_gvec(no_atoms))
  allocate(this%gvec_and_prim_to_atom(no_atoms/sc_size,sc_size))
  do gvec=1,sc_size
    do prim=1,no_atoms/sc_size
      atom = (prim-1)*sc_size + gvec
      
      this%atom_to_prim(atom) = prim
      this%atom_to_gvec(atom) = gvec
      this%gvec_and_prim_to_atom(prim,gvec) = atom
    enddo
  enddo
  
  this%no_symmetries = no_symmetries
  if (no_symmetries /= 0) then
    allocate(this%rotation_matrices(3,3,no_symmetries))
    allocate(this%offsets(3,no_symmetries))
    allocate(this%offsets_cart(3,no_symmetries))
  endif
  
  this%sc_size = sc_size
  allocate(this%gvectors(3,sc_size))
  allocate(this%paired_gvec(sc_size))
end subroutine

! Deallocates a Structure
subroutine drop_StructureData(this)
  implicit none
  
  type(StructureData), intent(inout) :: this
  
  deallocate(this%species)
  deallocate(this%mass)
  deallocate(this%atoms)
  
  deallocate(this%atom_to_prim)
  deallocate(this%atom_to_gvec)
  deallocate(this%gvec_and_prim_to_atom)
  
  if (this%no_symmetries /= 0) then
    deallocate(this%rotation_matrices)
    deallocate(this%offsets)
    deallocate(this%offsets_cart)
  endif
  
  deallocate(this%gvectors)
  deallocate(this%paired_gvec)
end subroutine

! reads structure.dat
function read_structure_file_character(filename) result(this)
  use constants,      only : identity
  use linear_algebra, only : invert, invert_int
  use string_module
  use file_module
  implicit none
  
  character(*),        intent(in) :: filename
  type(StructureData)             :: this
  
  type(String), allocatable :: structure_file(:)
  type(String), allocatable :: line(:)
  integer                   :: i,j
  integer                   :: no_atoms
  integer                   :: no_symmetries
  integer                   :: sc_size
  
  ! line numbers
  integer :: lattice_line   ! The line "Lattice"
  integer :: atoms_line     ! The line "Atoms"
  integer :: symmetry_line  ! The line "Symmetry"
  integer :: supercell_line ! The line "Supercell"
  integer :: gvectors_line  ! The line "G-vectors"
  integer :: end_line       ! The line "End"
  
  ! ------------------------------
  ! Initialise line numbers.
  ! ------------------------------
  lattice_line = 0
  atoms_line = 0
  symmetry_line = 0
  supercell_line = 0
  end_line = 0
  
  ! ------------------------------
  ! Work out layout of file.
  ! ------------------------------
  structure_file = read_lines(filename)
  do i=1,size(structure_file)
    line = split(lower_case(structure_file(i)))
    
    if (size(line)==0) then
      cycle
    endif
    
    if (line(1)=="lattice") then
      lattice_line = i
    elseif (line(1)=="atoms") then
      atoms_line = i
    elseif (line(1)=="symmetry") then
      symmetry_line = i
    elseif (line(1)=="supercell") then
      supercell_line = i
    elseif (line(1)=="g-vectors") then
      gvectors_line = i
    elseif (line(1)=="end") then
      end_line = i
    endif
  enddo
  
  ! ------------------------------
  ! Check layout is as expected.
  ! ------------------------------
  if (lattice_line/=1) then
    call print_line("Error: line 1 of "//filename//" is not 'Lattice'")
    call err()
  elseif (atoms_line/=5) then
    call print_line("Error: line 5 of "//filename//" is not 'Atoms'")
    call err()
  elseif (end_line/=size(structure_file)) then
    call print_line("Error: the last line of "//filename//" is not 'End'")
    call err()
  endif
  
  if (supercell_line/=0 .and. gvectors_line/=0) then
    if (gvectors_line-supercell_line/=4) then
      call print_line('Error: the lines "Supercell" and "G-vectors" in '// &
         & filename//' are not four lines apart.')
    endif
  endif
  
  ! ------------------------------
  ! Set counts.
  ! ------------------------------
  if (symmetry_line==0 .and. supercell_line==0) then
    ! structure.dat does not contain symmetries or supercell data
    no_atoms = end_line-atoms_line-1
    no_symmetries = 0
    sc_size = 1
  elseif (symmetry_line==0) then
    ! structure.dat does not contain symmetries
    no_atoms = supercell_line-atoms_line-1
    no_symmetries = 0
    sc_size = end_line-gvectors_line-1
  elseif (supercell_line==0) then
    ! structure.dat does not contain supercell data
    no_atoms = symmetry_line-atoms_line-1
    no_symmetries = (end_line-symmetry_line-1)/5
    sc_size = 1
  else
    no_atoms = symmetry_line-atoms_line-1
    no_symmetries = (supercell_line-symmetry_line-1)/5
    sc_size = (end_line-gvectors_line-1)
  endif
  
  ! ------------------------------
  ! Allocate structure.
  ! ------------------------------
  call new(this,no_atoms,no_symmetries,sc_size)
  
  ! ------------------------------
  ! Read file into arrays.
  ! ------------------------------
  do i=1,3
    this%lattice(i,:) = dble(split(structure_file(lattice_line+i)))
  enddo
  
  do i=1,this%no_atoms
    line = split(structure_file(atoms_line+i))
    this%species(i) = line(1)
    this%mass(i) = dble(line(2))
    this%atoms(:,i) = dble(line(3:5))
  enddo
  
  do i=1,this%no_symmetries
    do j=1,3
      line = split(structure_file(symmetry_line+(i-1)*5+j))
      this%rotation_matrices(j,:,i) = dble(line)
    enddo
    line = split(structure_file(symmetry_line+(i-1)*5+4))
    this%offsets(:,i) = dble(line)
  enddo
  
  if (supercell_line==0) then
    this%supercell = identity
    this%gvectors(:,1) = 0
  else
    do i=1,3
      this%supercell(i,:) = int(split(structure_file(supercell_line+i)))
    enddo
    
    do i=1,sc_size
      this%gvectors(:,i) = int(split(structure_file(gvectors_line+i)))
    enddo
  endif
  
  ! calculate derived quantities
  call calculate_derived_atom_quantities(this)
  
  if (this%no_symmetries>0) then
    call calculate_derived_symmetry_quantities(this)
  endif
  
  call calculate_derived_supercell_quantities(this)
end function

function read_structure_file_string(filename) result(this)
  use string_module
  implicit none
  
  type(String),        intent(in) :: filename
  type(StructureData)             :: this
  
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
      call print_line(structure_file, '')
    enddo
  endif
  
  call print_line(structure_file,'Supercell')
  do i=1,3
    call print_line(structure_file, this%supercell(i,:))
  enddo
  call print_line(structure_file,'G-vectors')
  do i=1,this%sc_size
    call print_line(structure_file, this%gvectors(:,i))
  enddo
  
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

! ----------------------------------------------------------------------
! calculate derived quantities relating to lattice and atoms
! ----------------------------------------------------------------------
subroutine calculate_derived_atom_quantities(this)
  use linear_algebra, only : invert, determinant
  implicit none
  
  type(StructureData), intent(inout) :: this
  
  this%recip_lattice = transpose(invert(this%lattice))
  this%volume        = abs(determinant(this%lattice))
end subroutine

! ----------------------------------------------------------------------
! calculate derived quantities relating to symmetries
! ----------------------------------------------------------------------
subroutine calculate_derived_symmetry_quantities(this)
  implicit none
  
  type(StructureData), intent(inout) :: this
  
  this%offsets_cart = matmul(this%lattice,this%offsets)
end subroutine

! ----------------------------------------------------------------------
! Calculate the derived quantities relating to the supercell and G-vectors.
! ----------------------------------------------------------------------
subroutine calculate_derived_supercell_quantities(this)
  use linear_algebra, only : invert_int
  implicit none
  
  type(StructureData), intent(inout) :: this
  
  integer :: i,j
  
  this%recip_supercell = invert_int(this%supercell)
  
  do i=1,this%sc_size
    do j=1,i
      if (all(modulo( this%gvectors(:,i)+this%gvectors(:,j), &
                    & this%sc_size)==0)) then
        this%paired_gvec(i) = j
        this%paired_gvec(j) = i
      endif
    enddo
  enddo
end subroutine

! ----------------------------------------------------------------------
! Calculate the relationships between G-vectors, modulo the supercell.
!    so if gvec(:,i)+gvec(:,j)=gvec(:,k) then operate(output(i),j)=k
! ----------------------------------------------------------------------
function calculate_gvector_group(this) result(output)
  use group_module
  implicit none
  
  type(StructureData), intent(in) :: this
  type(Group), allocatable        :: output(:)
  
  integer, allocatable :: operation(:)
  
  integer :: gvector_k(3)
  
  integer :: i,j,k
  
  allocate(operation(this%sc_size))
  allocate(output(this%sc_size))
  do i=1,this%sc_size
    do j=1,this%sc_size
      gvector_k = this%gvectors(:,i)+this%gvectors(:,j)
      do k=1,this%sc_size
        if (all(modulo(gvector_k-this%gvectors(:,k),this%sc_size)==0)) then
          operation(j) = k
        endif
      enddo
    enddo
    output(i) = operation
  enddo
end function
end module
