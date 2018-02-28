! ======================================================================
! All data relating to a given atomic configuration.
! ======================================================================
module structure_module
  use constants_module, only : dp
  use string_module
  use io_module
  
  use linear_algebra_module
  use fraction_algebra_module
  use group_module
  use atom_module
  use basic_symmetry_module
  use symmetry_module
  implicit none
  
  private
  
  public :: StructureData
  public :: read_structure_file
  public :: write_structure_file
  public :: construct_supercell
  
  ! ----------------------------------------------------------------------
  ! The structure type.
  ! ----------------------------------------------------------------------
  type :: StructureData
    ! ------------------------------
    ! Lattice data
    ! ------------------------------
    type(RealMatrix) :: lattice
    type(RealMatrix) :: recip_lattice
    real(dp)         :: volume
    
    ! ------------------------------
    ! Atom data
    ! ------------------------------
    integer :: no_atoms
    integer :: no_atoms_prim
    integer :: no_modes
    integer :: no_modes_prim
    
    type(AtomData), allocatable :: atoms(:)
    
    ! ------------------------------
    ! Symmetry data (in fractional co-ordinates).
    ! ------------------------------
    type(SymmetryOperator), allocatable :: symmetries(:)
    
    ! ------------------------------
    ! Superell data
    ! ------------------------------
    ! The number of primitive cells in the supercell.
    integer :: sc_size
    
    ! The lattice vectors of the supercell,
    !    in fractional primitive cell co-ordinates.
    type(IntMatrix) :: supercell
    
    ! invert(transpose(supercell))
    type(FractionMatrix) :: recip_supercell
    
    ! The R-vectors of the primitive cell which are not related by supercell
    !    lattice vectors.
    type(IntVector), allocatable :: rvectors(:)
    ! The G-vectors of the reciprocal supercell which are not related by
    !    primitive reciprocal cell vectors.
    type(IntVector), allocatable :: gvectors(:)
    
    ! --------------------------------------------------
    ! Groups describing rotation by symmetry or addition of vectors.
    ! --------------------------------------------------
    ! The groups describing symmetry operations.
    ! e.g. if S(i)*S(j)=S(k) then symmetry_group(i)*j = k.
    type(Group), allocatable :: symmetry_group(:)
    ! The groups describing R-vector and G-vector operations.
    ! e.g. if Rvec(i)+Rvec(j)=Rvec(k) then rvec_group(i)*j=k.
    type(Group), allocatable :: rvector_group(:)
    type(Group), allocatable :: gvector_group(:)
    
    ! --------------------------------------------------
    ! The IDs of inverse symmetries and paired R-vectors and G-vectors.
    ! Used to provide inverse and pair functions.
    ! --------------------------------------------------
    ! S(i)*S(inverse_symmetry(i)) = identity.
    integer, private, allocatable :: symmetry_inverse_ids_(:)
    ! The ID of the R-vector, j, s.t. rvectors(i) + rvectors(j) = 0.
    integer, allocatable :: rvector_paired_ids_(:)
    ! The ID of the G-vector, j, s.t. gvectors(i) + gvectors(j) = 0.
    integer, allocatable :: gvector_paired_ids_(:)
  contains
    ! Return inverse symmetries or paired R-vectors / G-vectors.
    procedure, public :: inverse_symmetries
    procedure, public :: paired_rvectors
    procedure, public :: paired_gvectors
    
    ! Return groups corresponding to inverse symmetries or paired vectors.
    procedure, public :: inverse_symmetry_group
    procedure, public :: paired_rvector_group
    procedure, public :: paired_gvector_group
    
    ! Procedures involved in constructing a StructureData.
    procedure, private :: calculate_symmetries
    procedure, private :: calculate_symmetry_group
    procedure, private :: calculate_symmetry_inverses
    procedure, private :: calculate_rvector_group
    procedure, private :: calculate_gvector_group
  end type
  
  interface StructureData
    module procedure new_StructureData
  end interface
contains

! ----------------------------------------------------------------------
! Functions to provide inverse symmetries and paired R-vectors and G-vectors,
!    and the groups corresponding to those objects.
! ----------------------------------------------------------------------
! Returns the inverse of a symmetry.
function inverse_symmetries(this,symmetry_id) result(output)
  implicit none
  
  class(StructureData), intent(in) :: this
  integer,              intent(in) :: symmetry_id
  type(SymmetryOperator)           :: output
  
  output = this%symmetries(this%symmetry_inverse_ids_(symmetry_id))
end function

! Returns the pair of an R-vector.
! i.e. paired_rvectors(i) + rvectors(i) = 0, modulo supercell R-vectors.
function paired_rvectors(this,rvector_id) result(output)
  implicit none
  
  class(StructureData), intent(in) :: this
  integer,              intent(in) :: rvector_id
  type(IntVector)                  :: output
  
  output = this%rvectors(this%rvector_paired_ids_(rvector_id))
end function

! Returns the pair of an G-vector.
! i.e. paired_gvectors(i) + gvectors(i) = 0, modulo primitive G-vectors.
function paired_gvectors(this,gvector_id) result(output)
  implicit none
  
  class(StructureData), intent(in) :: this
  integer,              intent(in) :: gvector_id
  type(IntVector)                  :: output
  
  output = this%gvectors(this%gvector_paired_ids_(gvector_id))
end function

! Return the symmetry group corresponding to the inverse of symmetry i.
function inverse_symmetry_group(this,symmetry_id) result(output)
  implicit none
  
  class(StructureData), intent(in) :: this
  integer,              intent(in) :: symmetry_id
  type(Group)                      :: output
  
  output = this%symmetry_group(this%symmetry_inverse_ids_(symmetry_id))
end function

! Return the group describing the subtraction of R-vector i.
function paired_rvector_group(this,rvector_id) result(output)
  implicit none
  
  class(StructureData), intent(in) :: this
  integer,              intent(in) :: rvector_id
  type(Group)                      :: output
  
  output = this%rvector_group(this%rvector_paired_ids_(rvector_id))
end function

! Return the group describing the subtraction of G-vector i.
function paired_gvector_group(this,gvector_id) result(output)
  implicit none
  
  class(StructureData), intent(in) :: this
  integer,              intent(in) :: gvector_id
  type(Group)                      :: output
  
  output = this%gvector_group(this%gvector_paired_ids_(gvector_id))
end function

! ----------------------------------------------------------------------
! Allocates all arrays, and sets no_ variables.
! ----------------------------------------------------------------------
function new_StructureData(lattice_matrix,supercell_matrix,rvectors,gvectors, &
   & species_primitive_cell,masses_primitive_cell,cartesian_positions,        &
   & basic_symmetries_in,symmetry_precision) result(this)
  use linear_algebra_module
  use basic_symmetry_module
  implicit none
  
  ! --------------------------------------------------
  ! Function arguments and Output.
  ! --------------------------------------------------
  ! Lattice and Supercell inputs.
  type(RealMatrix), intent(in) :: lattice_matrix
  type(IntMatrix),  intent(in) :: supercell_matrix
  type(IntVector),  intent(in) :: rvectors(:)
  type(IntVector),  intent(in) :: gvectors(:)
  
  ! Atom inputs.
  type(String),     intent(in) :: species_primitive_cell(:)
  real(dp),         intent(in) :: masses_primitive_cell(:)
  type(RealVector), intent(in) :: cartesian_positions(:,:)
  
  ! Symmetry inputs.
  ! Either basic_symmetries_in should be supplied (as read from a file),
  !    or symmetry_precision should be supplied (to calculate fresh symmetries).
  type(BasicSymmetry), intent(in), optional :: basic_symmetries_in(:)
  real(dp),            intent(in), optional :: symmetry_precision
  
  ! Output.
  type(StructureData) :: this
  
  ! --------------------------------------------------
  ! Working variables.
  ! --------------------------------------------------
  ! Symmetries.
  type(BasicSymmetry), allocatable :: basic_symmetries(:)
  
  ! Temporary variables.
  integer :: id,prim_id,rvec_id
  integer :: i,j,ialloc
  
  ! Fill out |S| and numbers of atoms and modes.
  this%sc_size = abs(determinant(supercell_matrix))
  this%no_atoms_prim = size(species_primitive_cell)
  this%no_atoms = this%no_atoms_prim*this%sc_size
  this%no_modes_prim = this%no_atoms_prim*3
  this%no_modes = this%no_atoms*3
  
  ! Check that inputs are consistent with |S|, no_atoms and one another.
  if (size(rvectors)/=this%sc_size) then
    call print_line(ERROR//': The size of the supercell, |S|, does not match &
       &the number of R-vectors provided.')
    call err()
  elseif (size(gvectors)/=this%sc_size) then
    call print_line(ERROR//': The size of the supercell, |S|, does not match &
       &the number of R-vectors provided.')
    call err()
  elseif (size(masses_primitive_cell)/=this%no_atoms_prim) then
    call print_line(ERROR//': The number of species provided does not match &
       &the number of masses provided.')
    call err()
  elseif (size(cartesian_positions,1)/=this%no_atoms_prim) then
    call print_line(ERROR//': The number of positions provided does not match &
       &the number of masses provided.')
    call err()
  elseif (size(cartesian_positions,2)/=this%sc_size) then
    call print_line(ERROR//': The number of positions provided does not match &
       &the number of R-vectors provided.')
  endif
  
  ! Fill out lattice, L, and supercell, S, matrix information.
  this%lattice         = lattice_matrix
  this%recip_lattice   = transpose(invert(lattice_matrix))
  this%supercell       = supercell_matrix
  this%recip_supercell = transpose(invert(this%supercell))
  this%volume          = abs(determinant(this%lattice))
  
  ! Fill out R-vectors (multiples of the primitive lattice which lie within
  !    the supercell lattice), and G-vectors (multiples of the reciprocal
  !    supercell lattice which lie within the reciprocal primitive lattice).
  this%rvectors = rvectors
  this%gvectors = gvectors
  
  allocate( this%rvector_paired_ids_(this%sc_size), &
          & this%gvector_paired_ids_(this%sc_size), &
          & stat=ialloc); call err(ialloc)
  
  this%rvector_paired_ids_ = 0
  do i=1,this%sc_size
    do j=1,this%sc_size
      if (is_int( this%recip_supercell &
              & * (this%rvectors(i)+this%rvectors(j)))) then
        if ( this%rvector_paired_ids_(i)==0 .and. &
           & this%rvector_paired_ids_(j)==0) then
          this%rvector_paired_ids_(i) = j
          this%rvector_paired_ids_(j) = i
        elseif ( this%rvector_paired_ids_(i)/=j .or. &
               & this%rvector_paired_ids_(j)/=i) then
          call print_line(ERROR//': Pairing of R-vectors inconsistent.')
          call err()
        endif
      endif
    enddo
  enddo
  if (any(this%rvector_paired_ids_==0)) then
    call print_line(ERROR//': not all paired R-vectors found.')
    call err()
  endif
  
  this%gvector_paired_ids_ = 0
  do i=1,this%sc_size
    do j=1,this%sc_size
      if (is_int( transpose(this%recip_supercell) &
              & * (this%gvectors(i)+this%gvectors(j)))) then
        if ( this%gvector_paired_ids_(i)==0 .and. &
           & this%gvector_paired_ids_(j)==0) then
          this%gvector_paired_ids_(i) = j
          this%gvector_paired_ids_(j) = i
        elseif ( this%gvector_paired_ids_(i)/=j .or. &
               & this%gvector_paired_ids_(j)/=i) then
          call print_line(ERROR//': Pairing of F-vectors inconsistent.')
          call err()
        endif
      endif
    enddo
  enddo
  if (any(this%gvector_paired_ids_==0)) then
    call print_line(ERROR//': not all paired G-vectors found.')
    call err()
  endif
  
  ! Fill out atom data.
  allocate( this%atoms(this%no_atoms), &
          & stat=ialloc); call err(ialloc)
  do id=1,this%no_atoms
    prim_id = modulo(id-1,this%no_atoms_prim) + 1
    rvec_id = (id-1)/this%no_atoms_prim + 1
    this%atoms(id) = AtomData( species_primitive_cell(prim_id),      &
                             & masses_primitive_cell(prim_id),       &
                             & cartesian_positions(prim_id,rvec_id), &
                             & this%lattice,                         &
                             & this%recip_lattice,                   &
                             & id,                                   &
                             & prim_id,                              &
                             & rvec_id)
  enddo
  
  ! Fill out symmetry information.
  if (present(symmetry_precision)) then
    if (present(basic_symmetries_in)) then
      call print_line(CODE_ERROR//': Both symmetry_precision and symmetries &
         &provided.')
      call err()
    endif
  else
    if (.not. present(basic_symmetries_in)) then
      call print_line(CODE_ERROR//': Neither symmetry_precision nor &
         &symmetries provided.')
      call err()
    endif
  endif
  
  if (present(symmetry_precision)) then
    basic_symmetries = calculate_basic_symmetries( this%lattice, &
                                                 & this%atoms, &
                                                 & symmetry_precision)
  else
    basic_symmetries = basic_symmetries_in
  endif
  
  if (size(basic_symmetries)/=0) then
    this%symmetries = this%calculate_symmetries(basic_symmetries)
    this%symmetry_group = this%calculate_symmetry_group()
    this%symmetry_inverse_ids_ =  this%calculate_symmetry_inverses()
  else
    this%symmetries = [SymmetryOperator::]
    this%symmetry_group = [Group::]
    this%symmetry_inverse_ids_ = [integer::]
  endif
  
  call this%calculate_rvector_group()
  call this%calculate_gvector_group()
end function

! ----------------------------------------------------------------------
! Reads structure.dat
! ----------------------------------------------------------------------
function read_structure_file(filename) result(this)
  use ifile_module
  implicit none
  
  type(String), intent(in) :: filename
  type(StructureData)      :: this
  
  type(IFile) :: structure_file
  
  ! line numbers
  integer :: lattice_line   ! The line "Lattice"
  integer :: atoms_line     ! The line "Atoms"
  integer :: symmetry_line  ! The line "Symmetry"
  integer :: supercell_line ! The line "Supercell"
  integer :: rvectors_line  ! The line "R-vectors"
  integer :: gvectors_line  ! The line "G-vectors"
  integer :: end_line       ! The line "End"
  
  ! Counts.
  integer     :: no_atoms
  integer     :: no_symmetries
  integer     :: sc_size
  integer     :: no_atoms_prim
  
  ! Lattice.
  real(dp) :: lattice_matrix(3,3)
  
  ! Supercell, R-vectors and G-vectors.
  integer                      :: supercell_matrix(3,3)
  type(IntVector), allocatable :: rvectors(:)
  type(IntVector), allocatable :: gvectors(:)
  
  ! Atom species, masses and cartesian positions.
  type(String),     allocatable :: species(:)
  real(dp),         allocatable :: masses(:)
  type(RealVector), allocatable :: positions(:,:)
  
  ! Symmetry data.
  integer                          :: rotation_matrix(3,3)
  type(BasicSymmetry), allocatable :: symmetries(:)
  
  ! Temporary variables.
  integer                   :: i,j,ialloc
  type(String), allocatable :: line(:)
  integer                   :: atom,rvec,prim
  
  ! ------------------------------
  ! Initialise line numbers.
  ! ------------------------------
  lattice_line = 0
  atoms_line = 0
  symmetry_line = 0
  supercell_line = 0
  rvectors_line = 0
  gvectors_line = 0
  end_line = 0
  
  ! ------------------------------
  ! Work out layout of file.
  ! ------------------------------
  structure_file = filename
  do i=1,size(structure_file)
    line = split(lower_case(structure_file%line(i)))
    
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
    elseif (line(1)=="r-vectors") then
      rvectors_line = i
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
    call print_line('Error: line 1 of '//filename//' is not "Lattice"')
    call err()
  elseif (atoms_line/=5) then
    call print_line('Error: line 5 of '//filename//' is not "Atoms"')
    call err()
  elseif (end_line/=size(structure_file)) then
    call print_line('Error: the last line of '//filename//' is not "End"')
    call err()
  elseif ( any([supercell_line,rvectors_line,gvectors_line]==0) .and. &
         & any([supercell_line,rvectors_line,gvectors_line]/=0)) then
    call print_line('Error: some but not all of Supercell, R-vectors and &
       &G-vectors are present in '//filename//'.')
    call err()
  endif
  
  if (supercell_line/=0) then
    if (rvectors_line-supercell_line/=4) then
      call print_line('Error: the lines "Supercell" and "R-vectors" in '// &
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
    no_symmetries = (end_line-symmetry_line-1)/7
    sc_size = 1
  else
    no_atoms = symmetry_line-atoms_line-1
    no_symmetries = (supercell_line-symmetry_line-1)/7
    sc_size = end_line-gvectors_line-1
  endif
  
  if (modulo(no_atoms,sc_size)/=0) then
    call print_line(ERROR//' The number of atoms is not a multiple of the &
       &number of R-vectors')
    call err()
  endif
  
  no_atoms_prim = no_atoms/sc_size
  
  ! ------------------------------
  ! Read in lattice.
  ! ------------------------------
  do i=1,3
    lattice_matrix(i,:) = dble(split(structure_file%line(lattice_line+i)))
  enddo
  
  ! ------------------------------
  ! Read in supercell, R-vectors and G-vectors.
  ! ------------------------------
  if (supercell_line==0) then
    supercell_matrix = int(make_identity_matrix(3))
    rvectors = [vec([0,0,0])]
    gvectors = [vec([0,0,0])]
  else
    do i=1,3
      supercell_matrix(i,:) = int(split(structure_file%line(supercell_line+i)))
    enddo
    
    allocate( rvectors(sc_size), &
            & gvectors(sc_size), &
            & stat=ialloc); call err(ialloc)
    do i=1,sc_size
      rvectors(i) = int(split(structure_file%line(rvectors_line+i)))
      gvectors(i) = int(split(structure_file%line(gvectors_line+i)))
    enddo
  endif
  
  ! ------------------------------
  ! Read in atomic information.
  ! ------------------------------
  allocate( species(no_atoms_prim),            &
          & masses(no_atoms_prim),             &
          & positions(no_atoms_prim, sc_size), &
          & stat=ialloc); call err(ialloc)
  do atom=1,no_atoms
    line = split(structure_file%line(atoms_line+atom))
    
    rvec = (atom-1)/no_atoms_prim + 1
    prim = modulo(atom-1,no_atoms_prim) + 1
    
    if (rvec==1) then
      species(prim) = line(1)
      masses(prim) = dble(line(2))
    else
      if (species(prim) /= line(1)) then
        call print_line(ERROR//': The species of the same atom at different &
           &R-vectors does not match.')
        call err()
      endif
    endif
    positions(prim,rvec) = dble(line(3:5))
  enddo
  
  ! ------------------------------
  ! Read in symmetries.
  ! ------------------------------
  allocate(symmetries(no_symmetries), stat=ialloc); call err(ialloc)
  do i=1,no_symmetries
    do j=1,3
      rotation_matrix(j,:) = int(split( &
         & structure_file%line(symmetry_line+(i-1)*7+j+1)))
    enddo
    symmetries(i)%rotation = rotation_matrix
    symmetries(i)%translation = &
       & dble(split(structure_file%line(symmetry_line+(i-1)*7+6)))
  enddo
  
  ! ------------------------------
  ! Allocate structure.
  ! ------------------------------
  this = StructureData( mat(lattice_matrix),   &
                      & mat(supercell_matrix), &
                      & rvectors,              &
                      & gvectors,              &
                      & species,               &
                      & masses,                &
                      & positions,             &
                      & basic_symmetries_in=symmetries)
end function

subroutine write_structure_file(this,filename)
  use ofile_module
  implicit none
  
  type(StructureData), intent(in) :: this
  type(String),        intent(in) :: filename
  
  type(OFile) :: structure_file
  
  integer :: i
  
  structure_file = filename
  
  call structure_file%print_line('Lattice')
  call structure_file%print_line(this%lattice)
  call structure_file%print_line('Atoms')
  do i=1,this%no_atoms
    call structure_file%print_line( this%atoms(i)%species() //' '// &
                                  & this%atoms(i)%mass()    //' '// &
                                  & this%atoms(i)%cartesian_position())
  enddo
  
  if (size(this%symmetries)/=0) then
    call structure_file%print_line('Symmetry')
    do i=1,size(this%symmetries)
      call structure_file%print_line('Rotation:')
      call structure_file%print_line(this%symmetries(i)%rotation)
      call structure_file%print_line('Translation:')
      call structure_file%print_line(this%symmetries(i)%translation)
      call structure_file%print_line('')
    enddo
  endif
  
  call structure_file%print_line('Supercell')
  call structure_file%print_line(this%supercell)
  call structure_file%print_line('R-vectors')
  do i=1,this%sc_size
    call structure_file%print_line(this%rvectors(i))
  enddo
  call structure_file%print_line('G-vectors')
  do i=1,this%sc_size
    call structure_file%print_line(this%gvectors(i))
  enddo
  
  call structure_file%print_line('End')
end subroutine

! ----------------------------------------------------------------------
! Calculate the relationships between R-vectors, modulo the supercell.
!    so if rvec(:,i)+rvec(:,j)=rvec(:,k) then output(i)*j=k
! ----------------------------------------------------------------------
subroutine calculate_rvector_group(this)
  use group_module
  implicit none
  
  class(StructureData), intent(inout) :: this
  
  integer, allocatable :: operation(:)
  
  type(IntVector) :: rvector_k
  
  integer :: i,j,k,ialloc
  
  allocate( operation(this%sc_size),          &
          & this%rvector_group(this%sc_size), &
          & stat=ialloc); call err(ialloc)
  do i=1,this%sc_size
    operation = 0
    do j=1,this%sc_size
      rvector_k = this%rvectors(i)+this%rvectors(j)
      do k=1,this%sc_size
        if (is_int(this%recip_supercell*(rvector_k-this%rvectors(k)))) then
          operation(j) = k
        endif
      enddo
    enddo
    if (any(operation==0)) then
      call print_line('Error: R-vector group incomplete.')
      call err()
    endif
    this%rvector_group(i) = Group(operation)
  enddo
end subroutine

! ----------------------------------------------------------------------
! Calculate the relationships between G-vectors, modulo the reciprocal
!    primitive lattice.
! so if gvec(:,i)+gvec(:,j)=gvec(:,k) then output(i)*j=k
! ----------------------------------------------------------------------
subroutine calculate_gvector_group(this)
  use group_module
  implicit none
  
  class(StructureData), intent(inout) :: this
  
  integer, allocatable :: operation(:)
  
  type(IntVector) :: gvector_k
  
  integer :: i,j,k,ialloc
  
  allocate( operation(this%sc_size), &
          & this%gvector_group(this%sc_size),    &
          & stat=ialloc); call err(ialloc)
  do i=1,this%sc_size
    operation = 0
    do j=1,this%sc_size
      gvector_k = this%gvectors(i)+this%gvectors(j)
      do k=1,this%sc_size
        if (is_int( transpose(this%recip_supercell) &
                & * (gvector_k-this%gvectors(k)))) then
          operation(j) = k
        endif
      enddo
    enddo
    if (any(operation==0)) then
      call print_line('Error: G-vector group incomplete.')
      call err()
    endif
    this%gvector_group(i) = Group(operation)
  enddo
end subroutine

! ----------------------------------------------------------------------
! Calculates the set of vectors which are not related to one another by
!    lattice vectors.
! ----------------------------------------------------------------------
! This is a helper function for construct_supercell.
! This is used to either calculate the R-vectors of the primitive cell which
!    are unique in the supercell, or the G-vectors of the reciprocal supercell
!    which are unique in the reciprocal primitive cell.
function calculate_unique_vectors(lattice,centre_on_origin) result(output)
  use linear_algebra_module
  implicit none
  
  ! Lattice vectors are the rows of lattice.
  type(IntMatrix), intent(in)  :: lattice
  logical,         intent(in)  :: centre_on_origin
  type(IntVector), allocatable :: output(:)
  
  integer         :: lattice_size
  integer         :: no_vectors
  type(IntVector) :: frac_vec
  type(IntVector) :: prim_vec
  integer         :: i,j,k,ialloc
  
  if (size(lattice,1)/=3 .or. size(lattice,2)/=3) then
    call print_line('Error: lattice is not 3x3.')
    call err()
  endif
  
  lattice_size = abs(determinant(lattice))
  
  allocate(output(lattice_size), stat=ialloc); call err(ialloc)
  
  no_vectors = 0
  do i=0,lattice_size-1
    do j=0,lattice_size-1
      do k=0,lattice_size-1
        
        ! Construct vectors in scaled fractional primitive co-ordinates.
        ! (scaled by lattice_size, to preserve integer representation).
        if (centre_on_origin) then
          frac_vec = [i-lattice_size/2,j-lattice_size/2,k-lattice_size/2]
        else
          frac_vec = [i,j,k]
        endif
        
        ! Transform to scaled fractional supercell co-ordinates.
        prim_vec = transpose(lattice)*frac_vec
        
        ! Check if the scaled co-ordinate is scaling*(integer co-ordinate).
        if (all(modulo(int(prim_vec),lattice_size)==0)) then
          no_vectors = no_vectors+1
          output(no_vectors) = prim_vec / lattice_size
        endif
      enddo
    enddo
  enddo
end function

! ----------------------------------------------------------------------
! Makes a supercell from a given primitive cell and supercell matrix.
! ----------------------------------------------------------------------
function construct_supercell(this,supercell_matrix,calculate_symmetries, &
   & symmetry_precision) result(supercell)
  use linear_algebra_module
  implicit none
  
  type(StructureData), intent(in)           :: this
  type(IntMatrix),     intent(in)           :: supercell_matrix
  logical,             intent(in), optional :: calculate_symmetries
  real(dp),            intent(in), optional :: symmetry_precision
  type(StructureData)                       :: supercell
  
  ! R-vector and G-vector information.
  type(IntVector), allocatable :: rvectors(:)
  type(IntVector), allocatable :: gvectors(:)
  
  ! Atomic information.
  type(String),     allocatable :: species(:)
  real(dp),         allocatable :: masses(:)
  type(RealVector), allocatable :: positions(:,:)
  
  ! Temporary variables
  logical             :: calculate_symmetries_flag
  integer             :: sc_size
  type(BasicSymmetry) :: empty_symmetries(0)
  
  integer :: prim,rvec,ialloc
  
  if (present(calculate_symmetries)) then
    calculate_symmetries_flag = calculate_symmetries
  else
    calculate_symmetries_flag = .true.
  endif
  
  ! Generate supercell and lattice data.
  sc_size = abs(determinant(supercell_matrix))
  
  ! Construct R-vectors and G-vectors.
  rvectors = calculate_unique_vectors(supercell_matrix, .false.)
  gvectors = calculate_unique_vectors(transpose(supercell_matrix), .true.)
  
  ! Check that rvectors(1) is the Gamma-point.
  if (rvectors(1)/=vec([0,0,0])) then
    call print_line(CODE_ERROR//': The first R-vector is not the Gamma point.')
    call err()
  endif
  
  ! Construct atomic information.
  allocate( species(this%no_atoms),           &
          & masses(this%no_atoms),            &
          & positions(this%no_atoms,sc_size), &
          & stat=ialloc); call err(ialloc)
  do prim=1,this%no_atoms
    species(prim) = this%atoms(prim)%species()
    masses(prim) = this%atoms(prim)%mass()
    do rvec=1,sc_size
      positions(prim,rvec) = transpose(this%lattice)                  &
                         & * ( this%atoms(prim)%fractional_position() &
                         &   + rvectors(rvec))
    enddo
  enddo
  
  ! Check symmetry precision has been supplied if required.
  if (calculate_symmetries_flag .and. .not. present(symmetry_precision)) then
    call print_line('Symmetry requested but no precision specified.')
    call err()
  endif
  
  ! Construct output.
  if (calculate_symmetries_flag) then
    supercell = StructureData(                                          &
       & supercell_matrix * this%lattice,                               &
       & supercell_matrix,                                              &
       & calculate_unique_vectors(supercell_matrix, .false.),           &
       & calculate_unique_vectors(transpose(supercell_matrix), .true.), &
       & species,                                                       &
       & masses,                                                        &
       & positions,                                                     &
       & symmetry_precision=symmetry_precision)
  else
    ! A compiler bug in gfortran 5.4 prevents [SymmetryOperator::] being passed
    !    directly to the function.
    empty_symmetries = [BasicSymmetry::]
    supercell = StructureData(                                          &
       & supercell_matrix * this%lattice,                               &
       & supercell_matrix,                                              &
       & calculate_unique_vectors(supercell_matrix, .false.),           &
       & calculate_unique_vectors(transpose(supercell_matrix), .true.), &
       & species,                                                       &
       & masses,                                                        &
       & positions,                                                     &
       & basic_symmetries_in=empty_symmetries)
  endif
end function

! ----------------------------------------------------------------------
! Calculates symmetry operations, in various representations.
! ----------------------------------------------------------------------
function calculate_symmetries(this,symmetries) result(output)
  use atom_module
  use basic_symmetry_module
  implicit none
  
  class(StructureData), intent(in)    :: this
  type(BasicSymmetry),  intent(in)    :: symmetries(:)
  type(SymmetryOperator), allocatable :: output(:)
  
  ! Atom group variables.
  type(RealVector)             :: transformed_position
  type(RealVector)             :: distance
  real(dp),        allocatable :: offsets(:)
  integer,         allocatable :: operations(:,:)
  type(IntVector), allocatable :: rvectors(:,:)
  integer,         allocatable :: prim_operations(:,:)
  type(IntVector), allocatable :: prim_rvectors(:,:)
  
  type(AtomData) :: atom_j
  type(AtomData) :: atom_k
  
  integer :: no_symmetries
  
  ! Temporary variables.
  integer :: i,j,k,ialloc
  
  no_symmetries = size(symmetries)
  
  allocate(output(no_symmetries), stat=ialloc); call err(ialloc)
  
  ! Returns if the symmetry is only a dummy placeholder.
  if (no_symmetries==0) then
    return
  endif

  ! --------------------------------------------------
  ! Calculates how symmetries map atoms onto other atoms.
  ! --------------------------------------------------
  ! If symmetry i maps atoms(j)%frac_pos to atoms(k)%frac_pos + R then
  !    - output(i)%atom_group * j = k.
  !    - output(i)%rvector(j)     = R.
  
  ! Work out which atoms map to which atoms under each symmetry operation.
  allocate( offsets(this%no_atoms),                             &
          & operations(this%no_atoms, no_symmetries),           &
          & rvectors(this%no_atoms, no_symmetries),             &
          & prim_operations(this%no_atoms_prim, no_symmetries), &
          & prim_rvectors(this%no_atoms_prim, no_symmetries),   &
          & stat=ialloc); call err(ialloc)
  prim_operations = 0
  do i=1,no_symmetries
    do j=1,this%no_atoms
      atom_j = this%atoms(j)
      
      ! Calculate the position of the transformed atom.
      transformed_position = symmetries(i)%rotation       &
                         & * atom_j%fractional_position() &
                         & + symmetries(i)%translation
      
      ! Identify which atom is closest to the transformed position,
      !    modulo supercell lattice vectors.
      do k=1,this%no_atoms
        distance = transformed_position - this%atoms(k)%fractional_position()
        offsets(k) = l2_norm(distance - vec(nint(distance)))
      enddo
      atom_k = this%atoms(minloc(offsets,1))
      operations(j,i) = atom_k%id()
      
      if (prim_operations(atom_j%prim_id(),i)==0) then
        prim_operations(atom_j%prim_id(),i) = atom_k%prim_id()
      else
        if (prim_operations(atom_j%prim_id(),i)/=atom_k%prim_id()) then
          call print_line( CODE_ERROR// &
                         & ': Transformation of atoms inconsistent.')
          call err()
        endif
      endif
      
      ! R * position(i) + t = position(j) + an R-vector.
      ! Identify this R-vector, and record it.
      rvectors(j,i) = nint( transformed_position &
                        & - atom_k%fractional_position())
      
      ! If atom j is in the primitive cell,
      !    record the primitive R-vector change.
      if (atom_j%prim_id()==atom_j%id()) then
        prim_rvectors(atom_j%id(),i) = this%rvectors(atom_k%rvec_id()) &
                                   & + transpose(this%supercell)       &
                                   & * rvectors(j,i)
      endif
      
      ! Check that the transformed atom is acceptably close to its image.
      if (offsets(operations(j,i))>1.0e-10_dp) then
        call print_line(CODE_ERROR//': Error mapping atoms under symmetry.')
        call err()
      elseif (l2_norm( transformed_position           &
                   & - atom_k%fractional_position() &
                   & - rvectors(j,i) )>1e-10_dp) then
        call print_line(CODE_ERROR//': Error mapping atoms under symmetry.')
        call err()
      endif
    enddo
  enddo
  
  ! Check that each symmetry is one-to-one, and that mapped atoms are of the
  !    same species.
  do i=1,no_symmetries
    do j=1,this%no_atoms
      if (count(operations(:,i)==j)/=1) then
        call print_line('Error: symmetry operation not one-to-one.')
        call err()
      endif
      
      if (this%atoms(operations(j,i))%species()/=this%atoms(j)%species()) then
        call print_line('Error: symmetry operation between different species.')
        call err()
      endif
    enddo
  enddo
  
  ! --------------------------------------------------
  ! Record basic symmetry properties.
  ! --------------------------------------------------
  do i=1,no_symmetries
    output(i) = SymmetryOperator( symmetries(i),               &
                                & this%lattice,                &
                                & this%recip_lattice,          &
                                & Group(operations(:,i)),      &
                                & rvectors(:,i),               &
                                & Group(prim_operations(:,i)), &
                                & prim_rvectors(:,i))
  enddo
end function

! --------------------------------------------------
! Calculates how symmetries map symmetries onto other symmetries.
! --------------------------------------------------
! symmetry_group(i) * j = k if symmetry i * symmetry j = symmetry k,
!    modulo translations by lattice vectors.
function calculate_symmetry_group(this) result(output)
  implicit none
  
  class(StructureData), intent(in) :: this
  type(Group), allocatable         :: output(:)
  
  integer, allocatable :: symmetry_group(:)
  type(IntMatrix)      :: rotation_ij
  type(Group)          :: atom_group_ij
  integer              :: identity
  
  integer :: i,j,k,ialloc
  
  allocate( symmetry_group(size(this%symmetries)), &
          & output(size(this%symmetries)),         &
          & stat=ialloc); call err(ialloc)
  do i=1,size(this%symmetries)
    symmetry_group = 0
    do j=1,size(this%symmetries)
      rotation_ij = this%symmetries(i)%rotation &
               & * this%symmetries(j)%rotation
      atom_group_ij = this%symmetries(i)%atom_group &
                & * this%symmetries(j)%atom_group
      do k=1,size(this%symmetries)
        if ( rotation_ij==this%symmetries(k)%rotation .and. &
           & atom_group_ij==this%symmetries(k)%atom_group) then
          if (symmetry_group(j)==0) then
            symmetry_group(j) = k
          elseif (symmetry_group(j)/=k) then
            call print_line(ERROR//': Symmetry group inconsistent.')
            call err()
          endif
        endif
      enddo
      
      if (symmetry_group(j)==0) then
        call print_line(ERROR//': symmetry '//i//' times symmetry '//j//' is &
           &not itself a symmetry.')
        call err()
      endif
    enddo
    
    output(i) = Group(symmetry_group)
  enddo
end function

! --------------------------------------------------
! Calculates the inverse of each symmetry.
! --------------------------------------------------
! If symmetry i * symmetry j = I then this%symmetry_inverse_ids_(i)=j.
function calculate_symmetry_inverses(this) result(output)
  implicit none
  
  class(StructureData), intent(in) :: this
  integer, allocatable             :: output(:)
  
  integer, allocatable :: symmetry_group(:)
  type(IntMatrix)      :: rotation_ij
  type(Group)          :: atom_group_ij
  integer              :: identity
  
  integer :: i,j,k,ialloc
  
  ! Locate the identity operator. S*S=S iff S=I.
  identity = 0
  do i=1,size(this%symmetries)
    if (this%symmetry_group(i)*i == i) then
      if (identity==0) then
        identity = i
      else
        call print_line(ERROR//': The identity symmetry has been found &
           &twice.')
        call err()
      endif
    endif
  enddo
  
  if (identity==0) then
    call print_line(ERROR//': The identity symmetry has not been found.')
    call err()
  endif
  
  ! Locate the inverse of each operator.
  allocate(output(size(this%symmetries)), stat=ialloc); call err(ialloc)
  output = 0
  do i=1,size(this%symmetries)
    do j=1,size(this%symmetries)
      if (this%symmetry_group(i)*j==identity) then
        if (output(i)==0 .and. output(j)==0) then
          output(i) = j
          output(j) = i
        elseif (output(i)/=j .or. output(j)/=i) then
          call print_line(ERROR//': Symmetry inverses are not consistent.')
          call err()
        endif
      endif
    enddo
  enddo
  
  if (any(output==0)) then
    call print_line(ERROR//': Unable to find all symmetry inverses.')
    call err()
  endif
end function
end module
