! ======================================================================
! All data relating to a given atomic configuration.
! ======================================================================
module structure_module
  use constants_module, only : dp
  use string_module
  use io_module
  
  use linear_algebra_module
  use group_module
  implicit none
  
  private
  
  public :: SymmetryOperator
  public :: StructureData
  public :: new
  public :: read_structure_file
  public :: write_structure_file
  
  ! ----------------------------------------------------------------------
  ! A symmetry operation.
  ! ----------------------------------------------------------------------
  type :: SymmetryOperator
    ! The rotation and translation in fractional co-ordinates.
    ! R and T.
    type(IntMatrix)  :: rotation
    type(RealVector) :: translation
    
    ! The mapping from atoms to other atoms.
    ! If R * x_i + T = x_j then atom_group*i = j,
    !    where x_i is the equilibrium position of atom i.
    type(Group) :: atom_group
    
    ! The mapping from symmetry operators to other symmetry operators.
    ! If S * S_i = S_j then symmetry_group*i = j,
    !    where S is this operator, and S_i and S_j are other operators.
    type(Group) :: operator_group
    
    ! The id of the operator S_j, s.g. S*S_j = S_j*S = I.
    integer :: inverse
  end type
  
  ! ----------------------------------------------------------------------
  ! The structure type.
  ! ----------------------------------------------------------------------
  type StructureData
    ! ------------------------------
    ! Lattice data
    ! ------------------------------
    type(RealMatrix) :: lattice
    type(RealMatrix) :: recip_lattice
    real(dp)         :: volume
    
    ! ------------------------------
    ! Atom data
    ! ------------------------------
    integer                       :: no_atoms
    integer                       :: no_atoms_prim
    integer                       :: no_modes
    integer                       :: no_modes_prim
    type(String),     allocatable :: species(:)
    real(dp),         allocatable :: mass(:)
    type(RealVector), allocatable :: atoms(:)
    
    ! ------------------------------
    ! Conversions between atom representations.
    !    [1...no_atoms] vs [1...no_atoms_in_prim]*[sc_size].
    ! ------------------------------
    ! Mapping from an atom in the supercell to
    !    the equivalent atom in the primitive cell.
    integer,      allocatable :: atom_to_prim(:)
    ! Mapping from an atom in the supercell to the R-vector of the primitive
    !    cell containing that atom.
    integer,      allocatable :: atom_to_rvec(:)
    ! The inverse of the above mappings.
    integer,      allocatable :: rvec_and_prim_to_atom(:,:)
    
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
    
    ! invert_int(transpose(supercell))
    type(IntMatrix) :: recip_supercell
    
    ! The R-vectors of the primitive cell which are not related by supercell
    !    lattice vectors.
    type(IntVector), allocatable :: rvectors(:)
    ! The ID of the R-vector, j, s.t. rvectors(:,i) + rvectors(:,j) = 0.
    integer, allocatable :: paired_rvec(:)
    ! The G-vectors of the reciprocal supercell which are not related by
    !    primitive reciprocal cell vectors.
    type(IntVector), allocatable :: gvectors(:)
    ! The ID of the G-vector, j, s.t. gvectors(:,i) + gvectors(:,j) = 0.
    integer, allocatable :: paired_gvec(:)
  contains
    procedure, public  :: calculate_derived_quantities
    procedure, public  :: calculate_symmetry
    procedure, public  :: calculate_rvector_group
    procedure, public  :: calculate_gvector_group
    procedure, public  :: calculate_cartesian_rotations
    procedure, private :: calculate_atom_symmetry_group
    procedure, private :: calculate_operator_symmetry_group
    procedure, private :: calculate_operator_inverses
  end type
  
  interface new
    module procedure new_StructureData
  end interface
  
  ! ----------------------------------------------------------------------
  ! spglib interface.
  ! ----------------------------------------------------------------------
  ! spglib_calculate_symmetries calculates symmetry operations, but has
  !    nowhere to store them. Space should be allocated, and passed into
  !    spglib_retrieve_symmetries, which will copy the data from
  !    C-style malloc memory to Fortran allocatable memory. Finally,
  !    drop_spg_dataset deallocates the malloc memory.
  interface
    subroutine spglib_calculate_symmetries( &
       & lattice,position,types,num_atom,symprec, &
       & success,spg_dataset_pointer,n_operations) &
       & bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      
      ! Inputs.
      real(kind=c_double),  intent(in) :: lattice(3,3)
      real(kind=c_double),  intent(in) :: position(3,*)
      integer(kind=c_int),  intent(in) :: types(*)
      integer(kind=c_int),  intent(in) :: num_atom
      real(kind=c_double),  intent(in) :: symprec
      
      ! Output.
      logical(kind=c_bool), intent(out) :: success
      type(c_ptr),          intent(out) :: spg_dataset_pointer
      integer(kind=c_int),  intent(out) :: n_operations
    end subroutine
  end interface
  
  interface
    subroutine spglib_retrieve_symmetries(spg_dataset_pointer, &
       & rotations,translations) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      
      ! Inputs.
      type(c_ptr), intent(in) :: spg_dataset_pointer
      
      ! Outputs.
      integer(kind=c_int), intent(out) :: rotations(3,3,*)
      real(kind=c_double), intent(out) :: translations(3,*)
    end subroutine
  end interface
  
  interface
    subroutine drop_spg_dataset(spg_dataset_pointer) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      
      type(c_ptr), intent(inout) :: spg_dataset_pointer
    end subroutine
  end interface
contains

! ----------------------------------------------------------------------
! Allocates all arrays, and sets no_ variables.
! ----------------------------------------------------------------------
subroutine new_StructureData(this,no_atoms,no_symmetries,sc_size)
  implicit none
  
  type(StructureData), intent(out) :: this
  integer,             intent(in)  :: no_atoms
  integer,             intent(in)  :: no_symmetries
  integer,             intent(in)  :: sc_size
  
  integer :: prim,gvec,atom
  integer :: ialloc
  
  this%no_atoms = no_atoms
  this%no_atoms_prim = no_atoms/sc_size
  this%no_modes = no_atoms*3
  this%no_modes_prim = no_atoms*3/sc_size
  allocate( this%species(no_atoms),                               &
          & this%mass(no_atoms),                                  &
          & this%atoms(no_atoms),                                 &
          & this%atom_to_prim(no_atoms),                          &
          & this%atom_to_rvec(no_atoms),                          &
          & this%rvec_and_prim_to_atom(no_atoms/sc_size,sc_size), &
          & stat=ialloc); call err(ialloc)
  do gvec=1,sc_size
    do prim=1,no_atoms/sc_size
      atom = (prim-1)*sc_size + gvec
      
      this%atom_to_prim(atom) = prim
      this%atom_to_rvec(atom) = gvec
      this%rvec_and_prim_to_atom(prim,gvec) = atom
    enddo
  enddo
  
  if (no_symmetries /= 0) then
    allocate( this%symmetries(no_symmetries), &
            & stat=ialloc); call err(ialloc)
  endif
  
  this%sc_size = sc_size
  allocate( this%rvectors(sc_size),    &
          & this%paired_rvec(sc_size), &
          & this%gvectors(sc_size),    &
          & this%paired_gvec(sc_size), &
          & stat=ialloc); call err(ialloc)
end subroutine

! ----------------------------------------------------------------------
! Reads structure.dat
! ----------------------------------------------------------------------
function read_structure_file(filename) result(this)
  use ifile_module
  implicit none
  
  type(String), intent(in) :: filename
  type(StructureData)      :: this
  
  type(IFile) :: structure_file
  integer     :: no_atoms
  integer     :: no_symmetries
  integer     :: sc_size
  
  ! line numbers
  integer :: lattice_line   ! The line "Lattice"
  integer :: atoms_line     ! The line "Atoms"
  integer :: symmetry_line  ! The line "Symmetry"
  integer :: supercell_line ! The line "Supercell"
  integer :: rvectors_line  ! The line "R-vectors"
  integer :: gvectors_line  ! The line "G-vectors"
  integer :: end_line       ! The line "End"
  
  ! Temproary variables.
  integer                   :: i,j
  type(String), allocatable :: line(:)
  real(dp)                  :: temp_real(3,3)
  integer                   :: temp_int(3,3)
  
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
    no_symmetries = (end_line-symmetry_line-1)/5
    sc_size = 1
  else
    no_atoms = symmetry_line-atoms_line-1
    no_symmetries = (supercell_line-symmetry_line-1)/11
    sc_size = end_line-gvectors_line-1
  endif
  
  ! ------------------------------
  ! Allocate structure.
  ! ------------------------------
  call new(this,no_atoms,no_symmetries,sc_size)
  
  ! ------------------------------
  ! Read file into arrays.
  ! ------------------------------
  do i=1,3
    temp_real(i,:) = dble(split(structure_file%line(lattice_line+i)))
  enddo
  this%lattice = temp_real
  
  do i=1,this%no_atoms
    line = split(structure_file%line(atoms_line+i))
    this%species(i) = line(1)
    this%mass(i) = dble(line(2))
    this%atoms(i) = dble(line(3:5))
  enddo
  
  do i=1,size(this%symmetries)
    do j=1,3
      temp_int(j,:) = int(split( &
         & structure_file%line(symmetry_line+(i-1)*11+j+1)))
    enddo
    this%symmetries(i)%rotation = temp_int
    this%symmetries(i)%translation = dble(split( &
       & structure_file%line(symmetry_line+(i-1)*11+6)))
    this%symmetries(i)%atom_group = int(split( &
       & structure_file%line(symmetry_line+(i-1)*11+8)))
    this%symmetries(i)%operator_group = int(split( &
       & structure_file%line(symmetry_line+(i-1)*11+10)))
  enddo
  
  if (supercell_line==0) then
    this%supercell = make_identity_matrix(3)
    this%gvectors(1) = [ 0, 0, 0 ]
  else
    do i=1,3
      temp_int(i,:) = int(split(structure_file%line(supercell_line+i)))
    enddo
    this%supercell = temp_int
    
    do i=1,sc_size
      this%rvectors(i) = int(split(structure_file%line(rvectors_line+i)))
    enddo
    
    do i=1,sc_size
      this%gvectors(i) = int(split(structure_file%line(gvectors_line+i)))
    enddo
  endif
  
  ! calculate derived quantities
  call calculate_derived_quantities(this)
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
    call structure_file%print_line( this%species(i)//' '// &
                                  & this%mass(i)//' '//    &
                                  & this%atoms(i))
  enddo
  
  if (size(this%symmetries)/=0) then
    call structure_file%print_line('Symmetry')
    do i=1,size(this%symmetries)
      call structure_file%print_line('Rotation:')
      call structure_file%print_line(this%symmetries(i)%rotation)
      call structure_file%print_line('Translation:')
      call structure_file%print_line(this%symmetries(i)%translation)
      call structure_file%print_line('Atom group:')
      call structure_file%print_line(this%symmetries(i)%atom_group)
      call structure_file%print_line('Operator group:')
      call structure_file%print_line(this%symmetries(i)%operator_group)
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
! Calculate derived quantities.
! ----------------------------------------------------------------------
subroutine calculate_derived_quantities(this)
  use linear_algebra_module, only : invert, determinant
  implicit none
  
  class(StructureData), intent(inout) :: this
  
  integer :: i,j
  
  this%recip_supercell = transpose(invert_int(this%supercell))
  this%recip_lattice   = transpose(invert(this%lattice))
  this%volume          = abs(determinant(this%lattice))
  
  this%paired_rvec = 0
  this%paired_gvec = 0
  
  do i=1,this%sc_size
    do j=1,i
      if (all(modulo( int( this%recip_supercell                  &
                    &    * (this%rvectors(i)+this%rvectors(j))), &
                    & this%sc_size) == 0)) then
        this%paired_rvec(i) = j
        this%paired_rvec(j) = i
      endif
      
      if (all(modulo( int( transpose(this%recip_supercell)       &
                &        * (this%gvectors(i)+this%gvectors(j))), &
                & this%sc_size) == 0)) then
        this%paired_gvec(i) = j
        this%paired_gvec(j) = i
      endif
    enddo
  enddo
  
  if (any(this%paired_rvec==0)) then
    call print_line('Error: not all paired R-vectors found.')
    call err()
  elseif (any(this%paired_gvec==0)) then
    call print_line('Error: not all paired G-vectors found.')
    call err()
  endif
end subroutine

! ----------------------------------------------------------------------
! Calculate the relationships between R-vectors, modulo the supercell.
!    so if rvec(:,i)+rvec(:,j)=rvec(:,k) then output(i)*j=k
! ----------------------------------------------------------------------
function calculate_rvector_group(this) result(output)
  use group_module
  implicit none
  
  class(StructureData), intent(in) :: this
  type(Group), allocatable         :: output(:)
  
  integer, allocatable :: operation(:)
  
  type(IntVector) :: rvector_k
  
  integer :: i,j,k,ialloc
  
  allocate( operation(this%sc_size), &
          & output(this%sc_size),    &
          & stat=ialloc); call err(ialloc)
  do i=1,this%sc_size
    operation = 0
    do j=1,this%sc_size
      rvector_k = this%rvectors(i)+this%rvectors(j)
      do k=1,this%sc_size
        if (all(modulo(int(this%recip_supercell*(rvector_k-this%rvectors(k))),&
                  & this%sc_size) == 0)) then
          operation(j) = k
        endif
      enddo
    enddo
    if (any(operation==0)) then
      call print_line('Error: R-vector group incomplete.')
      call err()
    endif
    output(i) = operation
  enddo
end function

! ----------------------------------------------------------------------
! Calculate the relationships between G-vectors, modulo the reciprocal
!    primitive lattice.
! so if gvec(:,i)+gvec(:,j)=gvec(:,k) then output(i)*j=k
! ----------------------------------------------------------------------
function calculate_gvector_group(this) result(output)
  use group_module
  implicit none
  
  class(StructureData), intent(in) :: this
  type(Group), allocatable         :: output(:)
  
  integer, allocatable :: operation(:)
  
  type(IntVector) :: gvector_k
  
  integer :: i,j,k,ialloc
  
  allocate( operation(this%sc_size), &
          & output(this%sc_size),    &
          & stat=ialloc); call err(ialloc)
  do i=1,this%sc_size
    operation = 0
    do j=1,this%sc_size
      gvector_k = this%gvectors(i)+this%gvectors(j)
      do k=1,this%sc_size
        if (all(modulo( int(   transpose(this%recip_supercell) &
                      &      * (gvector_k-this%gvectors(k))),  &
                      & this%sc_size) == 0)) then
          operation(j) = k
        endif
      enddo
    enddo
    if (any(operation==0)) then
      call print_line('Error: G-vector group incomplete.')
      call err()
    endif
    output(i) = operation
  enddo
end function

! ----------------------------------------------------------------------
! Transforms rotations into cartesian co-ordinates.
! ----------------------------------------------------------------------
function calculate_cartesian_rotations(this) result(output)
  implicit none
  
  class(StructureData), intent(in) :: this
  type(RealMatrix), allocatable    :: output(:)
  
  integer :: i,ialloc
  
  allocate(output(size(this%symmetries)), stat=ialloc); call err(ialloc)
  do i=1,size(this%symmetries)
    output(i) = transpose(this%lattice)     &
            & * this%symmetries(i)%rotation &
            & * this%recip_lattice
  enddo
end function

! ----------------------------------------------------------------------
! Calls spglib, and calculates all symmetry operations.
! ----------------------------------------------------------------------
subroutine calculate_symmetry(this,symmetry_precision)
  use iso_c_binding
  implicit none
  
  class(StructureData), intent(inout) :: this
  real(dp),             intent(in)    :: symmetry_precision
  
  ! spglib inputs.
  integer             :: atom_types(this%no_atoms)
  
  ! spglib data.
  logical(kind=c_bool) :: spglib_success
  type(c_ptr)          :: spg_dataset_pointer
  integer              :: no_symmetries
  
  ! Temporary matrix storage, for passing to and from C.
  real(dp), allocatable :: atoms(:,:)
  integer,  allocatable :: rotations(:,:,:)
  real(dp), allocatable :: translations(:,:)
  
  ! Temporary variables.
  integer :: i,j,ialloc
  integer :: atom_type
  
  atom_type = 0
  do_i : do i=1,this%no_atoms
    do j=1,i-1
      if (this%species(i)==this%species(j)) then
        atom_types(i) = atom_types(j)
        cycle do_i
      endif
    enddo
    atom_types(i) = atom_type
    atom_type = atom_type+1
  enddo do_i
  
  ! Calculate symmetries, without space to store them.
  allocate(atoms(3,this%no_atoms), stat=ialloc); call err(ialloc)
  do i=1,this%no_atoms
    atoms(:,i) = dble(this%recip_lattice*this%atoms(i))
  enddo
  call spglib_calculate_symmetries( dble(this%lattice),  &
                                  & atoms,               &
                                  & atom_types,          &
                                  & this%no_atoms,       &
                                  & symmetry_precision,  &
                                  & spglib_success,      &
                                  & spg_dataset_pointer, &
                                  & no_symmetries)
  
  ! Check for errors.
  if (.not. spglib_success) then
    call print_line('Error: spglib symmetry finder failed. &
       &Try increasing symmetry_precision.')
    call err()
  endif
  
  ! Allocate space for symmetries.
  allocate( rotations(3,3,no_symmetries),     &
          & translations(3,no_symmetries),    &
          & this%symmetries(no_symmetries), &
          & stat=ialloc); call err(ialloc)
  
  ! Retrieve symmetries into allocated space.
  call spglib_retrieve_symmetries( spg_dataset_pointer, &
                                 & rotations,      &
                                 & translations)
  do i=1,size(this%symmetries)
    this%symmetries(i)%rotation = rotations(:,:,i)
    this%symmetries(i)%translation = translations(:,i)
  enddo
  
  ! Deallocate C memory.
  call drop_spg_dataset(spg_dataset_pointer)
  
  ! Calculate further symmetry properties.
  call this%calculate_atom_symmetry_group()
  call this%calculate_operator_symmetry_group()
  call this%calculate_operator_inverses()
end subroutine

! ----------------------------------------------------------------------
! Calculates how symmetries map atoms onto other atoms.
! ----------------------------------------------------------------------
! If symmetry i maps atom j to atom k then output(i)*j=k.
subroutine calculate_atom_symmetry_group(this)
  use group_module
  use linear_algebra_module
  implicit none
  
  class(StructureData), intent(inout) :: this
  
  integer,  allocatable :: operations(:,:)
  
  ! Objects in fractional supercell co-ordinates.
  !   n.b. in these co-ordintates, supercell lattice vectors are unit vectors.
  !   This is a different convention to the scaled co-ordinates used elsewhere.
  type(RealVector), allocatable :: atom_pos_frac(:)
  type(RealVector)              :: transformed_pos_frac
  
  ! Distances between atoms and transformed atoms.
  type(RealVector)      :: delta
  real(dp), allocatable :: distances(:)
  
  ! Temporary variables.
  integer :: i,j,k
  
  allocate(atom_pos_frac(this%no_atoms))
  allocate(operations(this%no_atoms,size(this%symmetries)))
  allocate(distances(this%no_atoms))
  
  ! Transform atom positions into fractional supercell co-ordinates.
  do i=1,this%no_atoms
    atom_pos_frac(i) = this%recip_lattice * this%atoms(i)
  enddo
  
  ! Work out which atoms map to which atoms under each symmetry operation.
  do i=1,size(this%symmetries)
    do j=1,this%no_atoms
      ! Calculate the position of the transformed atom.
      transformed_pos_frac = this%symmetries(i)%rotation &
                         & * atom_pos_frac(j)                 &
                         & + this%symmetries(i)%translation
      
      ! Identify which atom is closest to the transformed position,
      !    modulo supercell lattice vectors.
      do k=1,this%no_atoms
        delta = transformed_pos_frac - atom_pos_frac(k)
        distances(k) = l2_norm(delta-vec(nint(dble(delta))))
      enddo
      operations(j,i) = minloc(distances,1)
      
      ! Check that the transformed atom is acceptably close to its image.
      if (distances(operations(j,i))>1.0e-10_dp) then
        call err()
      endif
    enddo
  enddo
  
  ! Check that each symmetry is one-to-one, and that mapped atoms are of the
  !    same species.
  do i=1,size(this%symmetries)
    do j=1,this%no_atoms
      if (count(operations(:,i)==j)/=1) then
        call print_line('Error: symmetry operation not one-to-one.')
        call err()
      endif
      
      if (this%species(operations(j,i))/=this%species(j)) then
        call print_line('Error: symmetry operation between different species.')
        call err()
      endif
    enddo
  enddo
  
  do i=1,size(this%symmetries)
    this%symmetries(i)%atom_group = operations(:,i)
  enddo
end subroutine

! ----------------------------------------------------------------------
! Calculates how symmetries map symmetries onto other symmetries.
! ----------------------------------------------------------------------
! If symmetry i * symmetry j = symmetry k then output(i)*j=k.
subroutine calculate_operator_symmetry_group(this)
  use group_module
  use linear_algebra_module
  implicit none
  
  class(StructureData), intent(inout) :: this
  
  type(IntMatrix) :: rotation_k
  type(Group)     :: operation_k
  
  ! Temporary variables.
  integer              :: i,j,k,ialloc
  integer, allocatable :: temp_operator_group(:)
  
  allocate( temp_operator_group(size(this%symmetries)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(this%symmetries)
    do_j : do j=1,size(this%symmetries)
      rotation_k = this%symmetries(i)%rotation &
               & * this%symmetries(j)%rotation
      operation_k = this%symmetries(i)%atom_group &
                & * this%symmetries(j)%atom_group
      do k=1,size(this%symmetries)
        if ( rotation_k==this%symmetries(k)%rotation .and. &
           & operation_k==this%symmetries(k)%atom_group) then
          temp_operator_group(j) = k
          cycle do_j
        endif
      enddo
      
      call print_line('Error: symmetry '//i//' times symmetry '//j//' is not &
         &itself a symmetry.')
      call err()
    enddo do_j
    
    this%symmetries(i)%operator_group = temp_operator_group
  enddo
end subroutine

! ----------------------------------------------------------------------
! Calculates the inverse of each symmetry.
! ----------------------------------------------------------------------
! If symmetry i * symmetry j = I then output(i)=j.
subroutine calculate_operator_inverses(this)
  use group_module
  implicit none
  
  class(StructureData), intent(inout) :: this
  
  type(Group) :: identity_group
  
  ! Temporary variables.
  integer :: i,j,ialloc
  integer :: identity
  
  ! Locate the identity operator. S*S=S iff S=I.
  identity = 0
  do i=1,size(this%symmetries)
    if (this%symmetries(i)%operator_group*i == i) then
      identity = i
    endif
  enddo
  
  if (identity==0) then
    call print_line('Error: The identity symmetry has not been found.')
    call err()
  endif
  
  ! Locate the inverse of each operator.
  do_i : do i=1,size(this%symmetries)
    do j=1,size(this%symmetries)
      if (this%symmetries(i)%operator_group*j==identity) then
        this%symmetries(i)%inverse = j
        
        if (j<i) then
          if (this%symmetries(j)%inverse/=i) then
            call print_line('Error: operator inverses do not match.')
            call err()
          endif
        endif
        
        cycle do_i
      endif
    enddo
    
    call print_line('Error: operator inverse not found.')
    call err()
  enddo do_i
end subroutine
end module
