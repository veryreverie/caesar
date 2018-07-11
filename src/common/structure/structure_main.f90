! ======================================================================
! All data relating to a given atomic configuration.
! ======================================================================
module structure_submodule
  use utils_module
  
  use basic_symmetry_submodule
  use basic_structure_submodule
  use atom_submodule
  use calculate_symmetry_submodule
  use symmetry_submodule
  implicit none
  
  private
  
  public :: StructureData
  
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
    ! Calculate the symmetry operators of the structure.
    procedure, public :: calculate_symmetry
    
    ! Return inverse symmetries or paired R-vectors / G-vectors.
    procedure, public :: inverse_symmetries
    procedure, public :: paired_rvectors
    procedure, public :: paired_gvectors
    
    ! Return whether to treat a symmetry as sin or cos when converting to
    !    real co-ordinates.
    procedure, public :: symmetry_is_sin
    
    ! Return groups corresponding to inverse symmetries or paired vectors.
    procedure, public :: inverse_symmetry_group
    procedure, public :: paired_rvector_group
    procedure, public :: paired_gvector_group
    
    ! Procedures involved in constructing a StructureData.
    procedure, private :: calculate_symmetry_operators
    procedure, private :: calculate_symmetry_group
    procedure, private :: calculate_symmetry_inverses
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

! Return whether to treat the symmetry as sin or cos when converting to
!    real co-ordinates.
! The choice is arbitrary, as long as one of each pair is sent each way,
!    and symmetries which are their own pair are sent to cos.
impure elemental function symmetry_is_sin(this,symmetry_id) result(output)
  implicit none
  
  class(StructureData), intent(in) :: this
  integer,              intent(in) :: symmetry_id
  logical                          :: output
  
  output = this%symmetry_inverse_ids_(symmetry_id)<symmetry_id
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
function new_StructureData(basic_structure,basic_supercell) result(this)
  implicit none
  
  type(BasicStructure), intent(in)           :: basic_structure
  type(BasicSupercell), intent(in), optional :: basic_supercell
  type(StructureData)                        :: this
  
  ! Working variables.
  integer,             allocatable :: atom_rvector_ids(:)
  integer,             allocatable :: atom_prim_ids(:)
  type(BasicSymmetry), allocatable :: symmetries(:)
  
  ! Temporary variables.
  integer :: id
  integer :: i,j,ialloc
  
  ! Copy structure information: lattice vectors and numbers of atoms.
  this%lattice         = basic_structure%lattice_matrix
  this%recip_lattice   = transpose(invert(this%lattice))
  this%volume          = abs(determinant(this%lattice))
  
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
    
    this%supercell       = basic_supercell%supercell_matrix
    this%recip_supercell = transpose(invert(this%supercell))
    this%rvectors        = basic_supercell%rvectors
    this%gvectors        = basic_supercell%gvectors
    atom_rvector_ids     = basic_supercell%atom_rvector_ids
    atom_prim_ids        = basic_supercell%atom_prim_ids
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
  this%symmetry_group = [Group::]
  this%symmetry_inverse_ids_ = [integer::]
end function

! ----------------------------------------------------------------------
! Calculate symmetries.
! ----------------------------------------------------------------------
subroutine calculate_symmetry(this,symmetry_precision,symmetries)
  implicit none
  
  class(StructureData), intent(inout)        :: this
  real(dp),             intent(in)           :: symmetry_precision
  type(BasicSymmetry),  intent(in), optional :: symmetries(:)
  
  type(BasicSymmetry), allocatable :: basic_symmetries(:)
  
  this%symmetry_precision = symmetry_precision
  
  if (present(symmetries)) then
    if (size(symmetries)==0) then
      call print_line(ERROR//': No symmetries given. There must be at least &
         &the identity symmetry.')
      call err()
    endif
    basic_symmetries = symmetries
  else
    basic_symmetries = calculate_basic_symmetries( this%lattice,      &
                                                 & this%atoms,        &
                                                 & symmetry_precision )
  endif
  
  this%symmetries = this%calculate_symmetry_operators( basic_symmetries,  &
                                                     & symmetry_precision )
  this%symmetry_group = this%calculate_symmetry_group()
  this%symmetry_inverse_ids_ = this%calculate_symmetry_inverses()
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
! Calculates symmetry operations, in various representations.
! ----------------------------------------------------------------------
function calculate_symmetry_operators(this,symmetries,symmetry_precision) &
   & result(output)
  implicit none
  
  class(StructureData), intent(in)    :: this
  type(BasicSymmetry),  intent(in)    :: symmetries(:)
  real(dp),             intent(in)    :: symmetry_precision
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
      if (offsets(operations(j,i))>symmetry_precision) then
        call print_line(CODE_ERROR//': An acepted symmetry maps atom '//i// &
           &' to further than symmetry_precision from atom '//operations(j,i))
        call err()
      elseif (l2_norm( transformed_position           &
                   & - atom_k%fractional_position() &
                   & - rvectors(j,i) )>symmetry_precision) then
        call print_line(CODE_ERROR//': An acepted symmetry maps atom '//i// &
           &' to further than symmetry_precision from atom '//operations(j,i))
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
  
  integer :: identity
  
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

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_StructureData(this,input)
  implicit none
  
  class(StructureData), intent(out) :: this
  type(String),         intent(in)  :: input(:)
  
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
  
  select type(this); type is(StructureData)
    ! ------------------------------
    ! Initialise line numbers.
    ! ------------------------------
    lattice_line = 0
    atoms_line = 0
    symmetry_line = 0
    symmetry_end_line = 0
    supercell_line = 0
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
      call print_line('Error: line 1 of structure.dat is not "Lattice"')
      call err()
    elseif (atoms_line/=5) then
      call print_line('Error: line 5 of structure.dat is not "Atoms"')
      call err()
    elseif (end_line/=size(input)) then
      call print_line('Error: the last line of structure.dat is not "End"')
      call err()
    elseif ( any([supercell_line,rvectors_line,gvectors_line]==0) .and. &
           & any([supercell_line,rvectors_line,gvectors_line]/=0)) then
      call print_line('Error: some but not all of Supercell, R-vectors and &
         &G-vectors are present in structure.dat.')
      call err()
    endif
    
    if (supercell_line/=0) then
      if (rvectors_line-supercell_line/=4) then
        call print_line('Error: the lines "Supercell" and "R-vectors" in &
           &structure.dat are not four lines apart.')
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
    ! Read in lattice.
    ! ------------------------------
    lattice_matrix = RealMatrix(input(lattice_line+1:lattice_line+3))
    
    ! ------------------------------
    ! Read in supercell, R-vectors and G-vectors.
    ! ------------------------------
    if (supercell_line==0) then
      supercell_matrix = make_identity_matrix(3)
      rvectors = [vec([0,0,0])]
      gvectors = [vec([0,0,0])]
    else
      supercell_matrix = IntMatrix(input(supercell_line+1:supercell_line+3))
      
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
    this = StructureData(                                      &
       & basic_structure = BasicStructure( lattice_matrix,     &
       &                                   species2,           &
       &                                   masses2,            &
       &                                   positions2      ),  &
       & basic_supercell = BasicSupercell( supercell_matrix,   &
       &                                   rvectors,           &
       &                                   gvectors,           &
       &                                   atom_rvector_ids,   &
       &                                   atom_prim_ids     ) )
    
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
    output = [                    &
             & str('Lattice'),    &
             & str(this%lattice), &
             & str('Atoms')       ]
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
    
    output = [ output,              &
             & str('Supercell'),    &
             & str(this%supercell), &
             & str('R-vectors'),    &
             & str(this%rvectors),  &
             & str('G-vectors'),    &
             & str(this%gvectors),  &
             & str('End')           ]
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
