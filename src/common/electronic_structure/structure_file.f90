! ======================================================================
! Reads and writes structure.dat files from StructureData.
! ======================================================================
module structure_file_submodule
  use utils_module
  
  use structure_module
  implicit none
  
  private
  
  public :: make_input_filename_caesar
  public :: read_structure_file
  public :: write_structure_file
contains

function make_input_filename_caesar() result(output)
  implicit none
  
  type(String) :: output
  
  output = 'structure.dat'
end function

function read_structure_file(filename) result(this)
  implicit none
  
  type(String), intent(in) :: filename
  type(StructureData)      :: this
  
  type(IFile)               :: structure_file
  type(String), allocatable :: lines(:)
  
  ! line numbers
  integer :: lattice_line      ! The line "Lattice"
  integer :: atoms_line        ! The line "Atoms"
  integer :: symmetry_line     ! The line "Symmetry"
  integer :: symmetry_end_line ! The last line of the Symmetry block.
  integer :: supercell_line    ! The line "Supercell"
  integer :: rvectors_line     ! The line "R-vectors"
  integer :: gvectors_line     ! The line "G-vectors"
  integer :: end_line          ! The line "End"
  
  ! Counts.
  integer     :: no_atoms
  integer     :: sc_size
  integer     :: no_atoms_prim
  
  ! Lattice.
  type(RealMatrix) :: lattice_matrix
  
  ! Supercell, R-vectors and G-vectors.
  type(IntMatrix)              :: supercell_matrix
  type(IntVector), allocatable :: rvectors(:)
  type(IntVector), allocatable :: gvectors(:)
  
  ! Atom species, masses and cartesian positions.
  type(String),     allocatable :: species(:)
  type(String),     allocatable :: species2(:)
  real(dp),         allocatable :: masses(:)
  real(dp),         allocatable :: masses2(:)
  type(RealVector), allocatable :: positions(:,:)
  type(RealVector), allocatable :: positions2(:)
  integer,          allocatable :: atom_rvector_ids(:)
  integer,          allocatable :: atom_prim_ids(:)
  
  ! Symmetry data.
  real(dp)                         :: symmetry_precision
  type(BasicSymmetry), allocatable :: symmetries(:)
  
  ! Temporary variables.
  integer                        :: i,ialloc
  type(StringArray), allocatable :: sections(:)
  type(String),      allocatable :: line(:)
  integer                        :: atom,rvec,prim
  
  ! ------------------------------
  ! Read in file.
  ! ------------------------------
  structure_file = IFile(filename)
  lines = structure_file%lines()
  this = StructureData(lines)
  
  !! ------------------------------
  !! Initialise line numbers.
  !! ------------------------------
  !lattice_line = 0
  !atoms_line = 0
  !symmetry_line = 0
  !symmetry_end_line = 0
  !supercell_line = 0
  !rvectors_line = 0
  !gvectors_line = 0
  !end_line = 0
  !
  !! ------------------------------
  !! Work out layout of file.
  !! ------------------------------
  !do i=1,size(lines)
  !  line = split_line(lower_case(lines(i)))
  !  
  !  if (size(line)==0) then
  !    cycle
  !  endif
  !  
  !  if (line(1)=="lattice") then
  !    lattice_line = i
  !  elseif (line(1)=="atoms") then
  !    atoms_line = i
  !  elseif (line(1)=="symmetry") then
  !    symmetry_line = i
  !  elseif (line(1)=="supercell") then
  !    supercell_line = i
  !  elseif (line(1)=="r-vectors") then
  !    rvectors_line = i
  !  elseif (line(1)=="g-vectors") then
  !    gvectors_line = i
  !  elseif (line(1)=="end") then
  !    end_line = i
  !  endif
  !enddo
  !
  !! ------------------------------
  !! Check layout is as expected.
  !! ------------------------------
  !if (lattice_line/=1) then
  !  call print_line('Error: line 1 of '//filename//' is not "Lattice"')
  !  call err()
  !elseif (atoms_line/=5) then
  !  call print_line('Error: line 5 of '//filename//' is not "Atoms"')
  !  call err()
  !elseif (end_line/=size(lines)) then
  !  call print_line('Error: the last line of '//filename//' is not "End"')
  !  call err()
  !elseif ( any([supercell_line,rvectors_line,gvectors_line]==0) .and. &
  !       & any([supercell_line,rvectors_line,gvectors_line]/=0)) then
  !  call print_line('Error: some but not all of Supercell, R-vectors and &
  !     &G-vectors are present in '//filename//'.')
  !  call err()
  !endif
  !
  !if (supercell_line/=0) then
  !  if (rvectors_line-supercell_line/=4) then
  !    call print_line('Error: the lines "Supercell" and "R-vectors" in '// &
  !       & filename//' are not four lines apart.')
  !  endif
  !endif
  !
  !! ------------------------------
  !! Set counts.
  !! ------------------------------
  !if (symmetry_line==0 .and. supercell_line==0) then
  !  ! structure.dat does not contain symmetries or supercell data
  !  no_atoms = end_line-atoms_line-1
  !  sc_size = 1
  !elseif (symmetry_line==0) then
  !  ! structure.dat does not contain symmetries
  !  no_atoms = supercell_line-atoms_line-1
  !  sc_size = end_line-gvectors_line-1
  !elseif (supercell_line==0) then
  !  ! structure.dat does not contain supercell data
  !  no_atoms = symmetry_line-atoms_line-1
  !  symmetry_end_line = end_line-1
  !  sc_size = 1
  !else
  !  no_atoms = symmetry_line-atoms_line-1
  !  symmetry_end_line = supercell_line-1
  !  sc_size = end_line-gvectors_line-1
  !endif
  !
  !if (modulo(no_atoms,sc_size)/=0) then
  !  call print_line(ERROR//' The number of atoms is not a multiple of the &
  !     &number of R-vectors')
  !  call err()
  !endif
  !
  !no_atoms_prim = no_atoms/sc_size
  !
  !! ------------------------------
  !! Read in lattice.
  !! ------------------------------
  !lattice_matrix = &
  !   & RealMatrix(lines(lattice_line+1:lattice_line+3))
  !
  !! ------------------------------
  !! Read in supercell, R-vectors and G-vectors.
  !! ------------------------------
  !if (supercell_line==0) then
  !  supercell_matrix = make_identity_matrix(3)
  !  rvectors = [vec([0,0,0])]
  !  gvectors = [vec([0,0,0])]
  !else
  !  supercell_matrix = &
  !     & IntMatrix(IntArray(lines(supercell_line+1:supercell_line+3)))
  !  
  !  allocate( rvectors(sc_size), &
  !          & gvectors(sc_size), &
  !          & stat=ialloc); call err(ialloc)
  !  do i=1,sc_size
  !    rvectors(i) = int(split_line(lines(rvectors_line+i)))
  !    gvectors(i) = int(split_line(lines(gvectors_line+i)))
  !  enddo
  !endif
  !
  !! ------------------------------
  !! Read in atomic information.
  !! ------------------------------
  !allocate( species(no_atoms_prim),            &
  !        & species2(no_atoms),                &
  !        & masses(no_atoms_prim),             &
  !        & masses2(no_atoms),                 &
  !        & positions(no_atoms_prim, sc_size), &
  !        & positions2(no_atoms),              &
  !        & atom_rvector_ids(no_atoms),        &
  !        & atom_prim_ids(no_atoms),           &
  !        & stat=ialloc); call err(ialloc)
  !do atom=1,no_atoms
  !  line = split_line(lines(atoms_line+atom))
  !  
  !  rvec = (atom-1)/no_atoms_prim + 1
  !  prim = modulo(atom-1,no_atoms_prim) + 1
  !  
  !  if (rvec==1) then
  !    species(prim) = line(1)
  !    masses(prim) = dble(line(2))
  !  else
  !    if (species(prim) /= line(1)) then
  !      call print_line(ERROR//': The species of the same atom at different &
  !         &R-vectors does not match.')
  !      call err()
  !    endif
  !  endif
  !  positions(prim,rvec)   = dble(line(3:5))
  !  species2(atom) = line(1)
  !  masses2(atom) = dble(line(2))
  !  positions2(atom)       = dble(line(3:5))
  !  atom_rvector_ids(atom) = rvec
  !  atom_prim_ids(atom) = prim
  !enddo
  !
  !! ------------------------------
  !! Construct output without symmetry.
  !! ------------------------------
  !this = StructureData(                                      &
  !   & basic_structure = BasicStructure( lattice_matrix,     &
  !   &                                   species2,           &
  !   &                                   masses2,            &
  !   &                                   positions2      ),  &
  !   & basic_supercell = BasicSupercell( supercell_matrix,   &
  !   &                                   rvectors,           &
  !   &                                   gvectors,           &
  !   &                                   atom_rvector_ids,   &
  !   &                                   atom_prim_ids     ) )
  !
  !! ------------------------------
  !! Read in symmetries if present, and add symmetry to output.
  !! ------------------------------
  !if (symmetry_line/=0) then
  !  line = split_line(lines(symmetry_line+1))
  !  symmetry_precision = dble(line(2))
  !  
  !  sections = split_into_sections(lines(symmetry_line+2:symmetry_end_line))
  !  symmetries = BasicSymmetry(sections)
  !  
  !  call this%calculate_symmetry(symmetry_precision, symmetries)
  !endif
end function

subroutine write_structure_file(this,filename)
  implicit none
  
  type(StructureData), intent(in) :: this
  type(String),        intent(in) :: filename
  
  type(OFile) :: structure_file
  
  type(BasicSymmetry), allocatable :: symmetries(:)
  
  integer :: i
  
  structure_file = OFile(filename)
  
  call structure_file%print_lines(this)
  
  !call structure_file%print_line('Lattice')
  !call structure_file%print_lines(this%lattice)
  !call structure_file%print_line('Atoms')
  !do i=1,this%no_atoms
  !  call structure_file%print_line( this%atoms(i)%species() //' '// &
  !                                & this%atoms(i)%mass()    //' '// &
  !                                & this%atoms(i)%cartesian_position())
  !enddo
  !
  !if (size(this%symmetries)/=0) then
  !  call structure_file%print_line('Symmetry')
  !  call structure_file%print_line('Precision: '//this%symmetry_precision)
  !  symmetries = BasicSymmetry(this%symmetries)
  !  call structure_file%print_lines( symmetries,        &
  !                                 & separating_line='' )
  !endif
  !
  !call structure_file%print_line('Supercell')
  !call structure_file%print_lines(this%supercell)
  !call structure_file%print_line('R-vectors')
  !do i=1,this%sc_size
  !  call structure_file%print_line(this%rvectors(i))
  !enddo
  !call structure_file%print_line('G-vectors')
  !do i=1,this%sc_size
  !  call structure_file%print_line(this%gvectors(i))
  !enddo
  !
  !call structure_file%print_line('End')
end subroutine
end module
