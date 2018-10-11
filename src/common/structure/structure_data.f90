! ======================================================================
! All data relating to a given atomic configuration.
! ======================================================================
module structure_data_module
  use utils_module
  
  use atom_module
  use spglib_module
  
  use basic_symmetry_module
  use basic_structure_module
  use symmetry_module
  implicit none
  
  private
  
  public :: StructureData
  public :: BasicStructure
  
  ! ----------------------------------------------------------------------
  ! The structure type.
  ! ----------------------------------------------------------------------
  type, extends(Stringsable) :: StructureData
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
    real(dp)                            :: symmetry_precision
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
    ! Groups describing addition of vectors.
    ! --------------------------------------------------
    ! The groups describing R-vector and G-vector operations.
    ! e.g. if Rvec(i)+Rvec(j)=Rvec(k) then rvec_group(i)*j=k.
    type(Group), allocatable :: rvector_group(:)
    type(Group), allocatable :: gvector_group(:)
    
    ! --------------------------------------------------
    ! The IDs of paired R-vectors and G-vectors.
    ! Used to provide inverse and pair functions.
    ! --------------------------------------------------
    ! The ID of the R-vector, j, s.t. rvectors(i) + rvectors(j) = 0.
    integer, allocatable :: rvector_paired_ids_(:)
    ! The ID of the G-vector, j, s.t. gvectors(i) + gvectors(j) = 0.
    integer, allocatable :: gvector_paired_ids_(:)
  contains
    ! Calculate the symmetry operators of the structure.
    procedure, public :: calculate_symmetry
    
    ! Return inverse symmetries or paired R-vectors / G-vectors.
    procedure, public :: paired_rvectors
    procedure, public :: paired_gvectors
    
    ! Return groups corresponding to paired vectors.
    procedure, public :: paired_rvector_group
    procedure, public :: paired_gvector_group
    
    ! Procedures involved in constructing a StructureData.
    procedure, private :: calculate_rvector_group
    procedure, private :: calculate_gvector_group
    
    ! I/O.
    procedure, public :: read  => read_StructureData
    procedure, public :: write => write_StructureData
  end type
  
  interface StructureData
    module procedure new_StructureData
    module procedure new_StructureData_Strings
    module procedure new_StructureData_StringArray
  end interface
  
  interface BasicStructure
    module procedure new_BasicStructure_StructureData
  end interface
contains

! ----------------------------------------------------------------------
! Functions to provide paired R-vectors and G-vectors,
!    and the groups corresponding to those objects.
! ----------------------------------------------------------------------
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
function new_StructureData(basic_structure,basic_supercell, &
   & recip_lattice_matrix,recip_supercell_matrix) result(this)
  implicit none
  
  type(BasicStructure), intent(in)           :: basic_structure
  type(BasicSupercell), intent(in), optional :: basic_supercell
  type(RealMatrix),     intent(in), optional :: recip_lattice_matrix
  type(FractionMatrix), intent(in), optional :: recip_supercell_matrix
  type(StructureData)                        :: this
  
  ! Working variables.
  integer,             allocatable :: atom_rvector_ids(:)
  integer,             allocatable :: atom_prim_ids(:)
  
  ! Temporary variables.
  integer :: id
  integer :: i,j,ialloc
  
  ! Copy structure information: lattice vectors and numbers of atoms.
  this%lattice = basic_structure%lattice_matrix
  if (present(recip_lattice_matrix)) then
    this%recip_lattice = recip_lattice_matrix
  else
    this%recip_lattice = transpose(invert(this%lattice))
  endif
  this%volume = abs(determinant(this%lattice))
  
  this%no_atoms = size(basic_structure%atoms)
  this%no_modes = this%no_atoms*3
  
  ! Copy supercell information: the supercell matrix, R-vectors and G-vectors.
  ! Assume this is the primitive cell if this information is not given.
  if (present(basic_supercell)) then
    this%sc_size = size(basic_supercell%rvectors)
    
    if (size(basic_supercell%atom_rvector_ids)/=this%no_atoms) then
      call print_line(ERROR//': incompatible structure and supercell &
         & information.')
      call err()
    elseif (modulo(this%no_atoms,this%sc_size)/=0) then
      call print_line(ERROR//': incompatible structure and supercell &
         & information.')
    endif
    
    this%supercell = basic_supercell%supercell_matrix
    if (present(recip_supercell_matrix)) then
      this%recip_supercell = recip_supercell_matrix
    else
      this%recip_supercell = transpose(invert(this%supercell))
    endif
    this%rvectors    = basic_supercell%rvectors
    this%gvectors    = basic_supercell%gvectors
    atom_rvector_ids = basic_supercell%atom_rvector_ids
    atom_prim_ids    = basic_supercell%atom_prim_ids
  else
    this%sc_size = 1
    
    this%supercell       = make_identity_matrix(3)
    this%recip_supercell = fracmat(this%supercell)
    this%rvectors        = [zeroes(3)]
    this%gvectors        = [zeroes(3)]
    atom_rvector_ids     = [(1,i=1,size(basic_structure%atoms))]
    atom_prim_ids        = [(i,i=1,size(basic_structure%atoms))]
  endif
  
  this%no_atoms_prim   = this%no_atoms/this%sc_size
  this%no_modes_prim   = this%no_atoms_prim*3
  
  ! Pair up R-vectors and G-vectors.
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
  
  call this%calculate_rvector_group()
  call this%calculate_gvector_group()
  
  ! Fill out atom data: atomic species, masses and positions.
  allocate( this%atoms(this%no_atoms), &
          & stat=ialloc); call err(ialloc)
  do id=1,this%no_atoms
    this%atoms(id) = AtomData( basic_structure%atoms(id), &
                             & this%lattice,              &
                             & this%recip_lattice,        &
                             & id,                        &
                             & atom_prim_ids(id),         &
                             & atom_rvector_ids(id))
  enddo
  
  ! Fill out symmetry information with dummy data.
  this%symmetry_precision = 0.0_dp
  this%symmetries = [SymmetryOperator::]
end function

! ----------------------------------------------------------------------
! Convert to a BasicStructure.
! ----------------------------------------------------------------------
impure elemental function new_BasicStructure_StructureData(this) result(output)
  implicit none
  
  type(StructureData), intent(in) :: this
  type(BasicStructure)            :: output
  
  output = BasicStructure(this%lattice, BasicAtom(this%atoms))
end function

! ----------------------------------------------------------------------
! Calculate symmetries.
! ----------------------------------------------------------------------
subroutine calculate_symmetry(this,symmetry_precision,symmetries, &
   & loto_direction)
  implicit none
  
  class(StructureData), intent(inout)        :: this
  real(dp),             intent(in)           :: symmetry_precision
  type(BasicSymmetry),  intent(in), optional :: symmetries(:)
  type(FractionVector), intent(in), optional :: loto_direction
  
  type(SpglibSymmetries)           :: spglib_symmetries
  type(BasicSymmetry), allocatable :: basic_symmetries(:)
  
  type(Group), allocatable :: symmetry_groups(:)
  
  integer :: no_symmetries
  
  ! Variables for calculating symmetry groups.
  integer, allocatable :: symmetry_group(:)
  type(IntMatrix)      :: tensor_ij
  type(Group)          :: atom_group_ij
  
  ! Variables for calculating symmetry inverses.
  integer              :: identity_id
  integer, allocatable :: inverse_symmetry_ids(:)
  
  integer :: i,j,ialloc
  
  this%symmetry_precision = symmetry_precision
  
  if (present(symmetries)) then
    if (size(symmetries)==0) then
      call print_line(ERROR//': No symmetries given. There must be at least &
         &the identity symmetry.')
      call err()
    endif
    basic_symmetries = symmetries
    if (present(loto_direction)) then
      if (any(loto_direction*basic_symmetries%tensor/=loto_direction)) then
        call print_line(ERROR//': A symmetry does not leave the LO/TO &
           &direction invariant.')
        call err()
      endif
    endif
  else
    spglib_symmetries = SpglibSymmetries( this%lattice,      &
                                        & this%atoms,        &
                                        & symmetry_precision )
    basic_symmetries = BasicSymmetry( spglib_symmetries, &
                                    & this%atoms,        &
                                    & symmetry_precision )
    if (present(loto_direction)) then
      basic_symmetries = basic_symmetries(filter(                 &
         & loto_direction*basic_symmetries%tensor==loto_direction ))
    endif
  endif
  
  no_symmetries = size(basic_symmetries)
  
  ! Returns if the symmetry is only a dummy placeholder.
  if (no_symmetries==0) then
    return
  endif

  ! --------------------------------------------------
  ! Calculates how symmetries map symmetries onto other symmetries.
  ! --------------------------------------------------
  
  ! Calculate symmetry groups.
  allocate( symmetry_group(no_symmetries),  &
          & symmetry_groups(no_symmetries), &
          & stat=ialloc); call err(ialloc)
  do i=1,no_symmetries
    symmetry_group = 0
    do j=1,no_symmetries
      tensor_ij = basic_symmetries(i)%tensor &
              & * basic_symmetries(j)%tensor
      atom_group_ij = basic_symmetries(i)%atom_group &
                  & * basic_symmetries(j)%atom_group
      
      if (count( basic_symmetries%tensor==tensor_ij .and.   &
               & basic_symmetries%atom_group==atom_group_ij )/=1) then
        call print_line(ERROR//': symmetry group inconsistent.')
        call err()
      endif
      
      symmetry_group(j) = first( basic_symmetries%tensor==tensor_ij .and.   &
                               & basic_symmetries%atom_group==atom_group_ij )
    enddo
    
    symmetry_groups(i) = Group(symmetry_group)
  enddo
  
  ! Locate the identity operator. S*S=S iff S=I.
  identity_id = 0
  do i=1,no_symmetries
    if (symmetry_groups(i)*i == i) then
      if (identity_id==0) then
        identity_id = i
      else
        call print_line(ERROR//': The identity symmetry has been found &
           &twice.')
        call err()
      endif
    endif
  enddo
  
  if (identity_id==0) then
    call print_line(ERROR//': The identity symmetry has not been found.')
    call err()
  endif
  
  ! Locate the inverse of each operator.
  allocate(inverse_symmetry_ids(no_symmetries), stat=ialloc); call err(ialloc)
  inverse_symmetry_ids = 0
  do i=1,no_symmetries
    do j=1,no_symmetries
      if (symmetry_groups(i)*j==identity_id) then
        if ( inverse_symmetry_ids(i)==0 .and. &
           & inverse_symmetry_ids(j)==0       ) then
          inverse_symmetry_ids(i) = j
          inverse_symmetry_ids(j) = i
        elseif ( inverse_symmetry_ids(i)/=j .or. &
               & inverse_symmetry_ids(j)/=i      ) then
          call print_line(ERROR//': Symmetry inverses are not consistent.')
          call err()
        endif
      endif
    enddo
  enddo
  
  if (any(inverse_symmetry_ids==0)) then
    call print_line(ERROR//': Unable to find all symmetry inverses.')
    call err()
  endif
  
  ! --------------------------------------------------
  ! Record basic symmetry properties.
  ! --------------------------------------------------
  deallocate(this%symmetries, stat=ialloc); call err(ialloc)
  allocate(this%symmetries(no_symmetries), stat=ialloc); call err(ialloc)
  do i=1,no_symmetries
    this%symmetries(i) = SymmetryOperator( basic_symmetries(i),    &
                                         & this%lattice,           &
                                         & this%recip_lattice,     &
                                         & symmetry_groups(i),     &
                                         & inverse_symmetry_ids(i) )
  enddo
end subroutine

! ----------------------------------------------------------------------
! Calculate the relationships between R-vectors, modulo the supercell.
!    so if rvec(:,i)+rvec(:,j)=rvec(:,k) then output(i)*j=k
! ----------------------------------------------------------------------
subroutine calculate_rvector_group(this)
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
! I/O.
! ----------------------------------------------------------------------
subroutine read_StructureData(this,input)
  implicit none
  
  class(StructureData), intent(out) :: this
  type(String),         intent(in)  :: input(:)
  
  ! line numbers
  integer :: lattice_line         ! The line "Lattice"
  integer :: recip_lattice_line   ! The line "Reciprocal Lattice"
  integer :: atoms_line           ! The line "Atoms"
  integer :: symmetry_line        ! The line "Symmetry"
  integer :: symmetry_end_line    ! The last line of the Symmetry block.
  integer :: supercell_line       ! The line "Supercell"
  integer :: recip_supercell_line ! The line "Reciprocal Supercell"
  integer :: rvectors_line        ! The line "R-vectors"
  integer :: gvectors_line        ! The line "G-vectors"
  integer :: end_line             ! The line "End"
  
  ! Counts.
  integer     :: no_atoms
  integer     :: sc_size
  integer     :: no_atoms_prim
  
  ! Lattice.
  type(RealMatrix) :: lattice_matrix
  type(RealMatrix) :: recip_lattice_matrix
  
  ! Supercell, R-vectors and G-vectors.
  type(IntMatrix)              :: supercell_matrix
  type(FractionMatrix)         :: recip_supercell_matrix
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
  
  ! Output data.
  type(BasicStructure) :: basic_structure
  type(BasicSupercell) :: basic_supercell
  
  ! Temporary variables.
  integer                        :: i,ialloc
  type(StringArray), allocatable :: sections(:)
  type(String),      allocatable :: line(:)
  integer                        :: atom,rvec,prim
  
  select type(this); type is(StructureData)
    ! ------------------------------
    ! Initialise line numbers.
    ! ------------------------------
    lattice_line = 0
    recip_lattice_line = 0
    atoms_line = 0
    symmetry_line = 0
    symmetry_end_line = 0
    supercell_line = 0
    recip_supercell_line = 0
    rvectors_line = 0
    gvectors_line = 0
    end_line = 0
    
    ! ------------------------------
    ! Work out layout of file.
    ! ------------------------------
    do i=1,size(input)
      line = split_line(lower_case(input(i)))
      
      if (size(line)==0) then
        cycle
      endif
      
      if (line(1)=='lattice') then
        lattice_line = i
      elseif (line(1)=='atoms') then
        atoms_line = i
      elseif (line(1)=='symmetry') then
        symmetry_line = i
      elseif (line(1)=='supercell') then
        supercell_line = i
      elseif (line(1)=='r-vectors') then
        rvectors_line = i
      elseif (line(1)=='g-vectors') then
        gvectors_line = i
      elseif (line(1)=='end') then
        end_line = i
      elseif (line(1)=='reciprocal') then
        if (size(line)>1) then
          if (line(2)=='lattice') then
            recip_lattice_line = i
          elseif (line(2)=='supercell') then
            recip_supercell_line = i
          endif
        endif
      endif
    enddo
    
    ! ------------------------------
    ! Check layout is as expected.
    ! ------------------------------
    if (lattice_line/=1) then
      call print_line(ERROR//': line 1 of structure.dat is not "Lattice"')
      call err()
    elseif (recip_lattice_line/=0 .and. recip_lattice_line/=5) then
      call print_line(ERROR//': the line "Reciprocal Lattice" is not line 5 &
         &of structure.dat')
      call err()
    elseif (recip_lattice_line==0 .and. atoms_line/=5) then
      call print_line(ERROR//': line 5 of structure.dat is not four lines &
         &after "Lattice".')
      call err()
    elseif (recip_lattice_line/=0 .and. atoms_line/=9) then
      call print_line(ERROR//': line 5 of structure.dat is not four lines &
         &after "Reciprocal Lattice".')
      call err()
    elseif (end_line/=size(input)) then
      call print_line(ERROR//': the last line of structure.dat is not "End"')
      call err()
    elseif ( any([supercell_line,rvectors_line,gvectors_line]==0) .and. &
           & any([supercell_line,rvectors_line,gvectors_line]/=0)) then
      call print_line(ERROR//': some but not all of Supercell, R-vectors and &
         &G-vectors are present in structure.dat.')
      call err()
    elseif (supercell_line==0 .and. recip_supercell_line/=0) then
      call print_line(ERROR//': Reciprocal Supercell is present in &
         &structure.dat, but Supercell is not.')
      call err()
    elseif ( recip_supercell_line/=0 .and.          &
           & recip_supercell_line-supercell_line/=4 ) then
      call print_line(ERROR//': The line "Reciprocal Supercell" is not four &
         &lines after the line "Supercell" in structure.dat.')
      call err()
    endif
    
    if (supercell_line/=0) then
      if (recip_supercell_line==0) then
        if (rvectors_line-supercell_line/=4) then
          call print_line(ERROR//': the lines "Supercell" and "R-vectors" in &
             &structure.dat are not four lines apart.')
        endif
      else
        if (rvectors_line-recip_supercell_line/=4) then
          call print_line(ERROR//': the lines "Reciprocal Supercell" and &
             &"R-vectors" in structure.dat are not four lines apart.')
        endif
      endif
    endif
    
    ! ------------------------------
    ! Set counts.
    ! ------------------------------
    if (symmetry_line==0 .and. supercell_line==0) then
      ! structure.dat does not contain symmetries or supercell data
      no_atoms = end_line-atoms_line-1
      sc_size = 1
    elseif (symmetry_line==0) then
      ! structure.dat does not contain symmetries
      no_atoms = supercell_line-atoms_line-1
      sc_size = end_line-gvectors_line-1
    elseif (supercell_line==0) then
      ! structure.dat does not contain supercell data
      no_atoms = symmetry_line-atoms_line-1
      symmetry_end_line = end_line-1
      sc_size = 1
    else
      no_atoms = symmetry_line-atoms_line-1
      symmetry_end_line = supercell_line-1
      sc_size = end_line-gvectors_line-1
    endif
    
    if (modulo(no_atoms,sc_size)/=0) then
      call print_line(ERROR//' The number of atoms is not a multiple of the &
         &number of R-vectors')
      call err()
    endif
    
    no_atoms_prim = no_atoms/sc_size
    
    ! ------------------------------
    ! Read in lattice (and reciprocal lattice if present).
    ! ------------------------------
    lattice_matrix = RealMatrix(input(lattice_line+1:lattice_line+3))
    if (recip_lattice_line/=0) then
      recip_lattice_matrix = RealMatrix(                    &
         & input(recip_lattice_line+1:recip_lattice_line+3) )
    endif
    
    ! ------------------------------
    ! Read in supercell, reciprocal supercell, R-vectors and G-vectors.
    ! ------------------------------
    if (supercell_line==0) then
      supercell_matrix = make_identity_matrix(3)
      rvectors = [vec([0,0,0])]
      gvectors = [vec([0,0,0])]
    else
      supercell_matrix = IntMatrix(input(supercell_line+1:supercell_line+3))
      
      if (recip_supercell_line/=0) then
        recip_supercell_matrix = FractionMatrix(                  &
           & input(recip_supercell_line+1:recip_supercell_line+3) )
      endif
      
      allocate( rvectors(sc_size), &
              & gvectors(sc_size), &
              & stat=ialloc); call err(ialloc)
      do i=1,sc_size
        rvectors(i) = int(split_line(input(rvectors_line+i)))
        gvectors(i) = int(split_line(input(gvectors_line+i)))
      enddo
    endif
    
    ! ------------------------------
    ! Read in atomic information.
    ! ------------------------------
    allocate( species(no_atoms_prim),            &
            & species2(no_atoms),                &
            & masses(no_atoms_prim),             &
            & masses2(no_atoms),                 &
            & positions(no_atoms_prim, sc_size), &
            & positions2(no_atoms),              &
            & atom_rvector_ids(no_atoms),        &
            & atom_prim_ids(no_atoms),           &
            & stat=ialloc); call err(ialloc)
    do atom=1,no_atoms
      line = split_line(input(atoms_line+atom))
      
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
      positions(prim,rvec)   = dble(line(3:5))
      species2(atom) = line(1)
      masses2(atom) = dble(line(2))
      positions2(atom)       = dble(line(3:5))
      atom_rvector_ids(atom) = rvec
      atom_prim_ids(atom) = prim
    enddo
    
    ! ------------------------------
    ! Construct output without symmetry.
    ! ------------------------------
    basic_structure = BasicStructure( lattice_matrix, &
                                      species2,       &
                                      masses2,        &
                                      positions2      )
    basic_supercell = BasicSupercell( supercell_matrix, &
                                      rvectors,         &
                                      gvectors,         &
                                      atom_rvector_ids, &
                                      atom_prim_ids     )
    if (recip_lattice_line==0 .and. recip_supercell_line==0) then
      this = StructureData(basic_structure, basic_supercell)
    elseif (recip_lattice_line==0) then
      this = StructureData( basic_structure,                                &
                          & basic_supercell,                                &
                          & recip_supercell_matrix = recip_supercell_matrix )
    elseif (recip_supercell_line==0) then
      this = StructureData( basic_structure,                            &
                          & basic_supercell,                            &
                          & recip_lattice_matrix = recip_lattice_matrix )
    else
      this = StructureData( basic_structure,                                &
                          & basic_supercell,                                &
                          & recip_lattice_matrix = recip_lattice_matrix,    &
                          & recip_supercell_matrix = recip_supercell_matrix )
    endif
    
    ! ------------------------------
    ! Read in symmetries if present, and add symmetry to output.
    ! ------------------------------
    if (symmetry_line/=0) then
      line = split_line(input(symmetry_line+1))
      symmetry_precision = dble(line(2))
      
      sections = split_into_sections(input(symmetry_line+2:symmetry_end_line))
      symmetries = BasicSymmetry(sections)
      
      call this%calculate_symmetry(symmetry_precision, symmetries)
    endif
  class default
    call err()
  end select
end subroutine

function write_StructureData(this) result(output)
  implicit none
  
  class(StructureData), intent(in) :: this
  type(String), allocatable        :: output(:)
  
  type(BasicSymmetry), allocatable :: symmetries(:)
  
  integer :: i
  
  select type(this); type is(StructureData)
    output = [                            &
             & str('Lattice'),            &
             & str(this%lattice),         &
             & str('Reciprocal Lattice'), &
             & str(this%recip_lattice),   &
             & str('Atoms')               ]
    do i=1,this%no_atoms
      output = [ output, &
               & this%atoms(i)%species() //' '//    &
               & this%atoms(i)%mass()    //' '//    &
               & this%atoms(i)%cartesian_position() ]
    enddo
    
    if (size(this%symmetries)/=0) then
      symmetries = BasicSymmetry(this%symmetries)
      output = [ output,                                      &
               & str('Symmetry'),                             &
               & str('Precision: '//this%symmetry_precision), &
               & str(symmetries, separating_line='')          ]
    endif
    
    output = [ output,                      &
             & str('Supercell'),            &
             & str(this%supercell),         &
             & str('Reciprocal Supercell'), &
             & str(this%recip_supercell),   &
             & str('R-vectors'),            &
             & str(this%rvectors),          &
             & str('G-vectors'),            &
             & str(this%gvectors),          &
             & str('End')                   ]
  class default
    call err()
  end select
end function

function new_StructureData_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(StructureData)      :: this
  
  call this%read(input)
end function

impure elemental function new_StructureData_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(StructureData)           :: this
  
  this = StructureData(str(input))
end function
end module
