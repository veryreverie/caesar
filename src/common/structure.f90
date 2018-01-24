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
    ! The ID of the R-vector, j, s.t. rvectors(:,i) + rvectors(:,j) = 0.
    integer, allocatable :: paired_rvec(:)
    ! The G-vectors of the reciprocal supercell which are not related by
    !    primitive reciprocal cell vectors.
    type(IntVector), allocatable :: gvectors(:)
    ! The ID of the G-vector, j, s.t. gvectors(:,i) + gvectors(:,j) = 0.
    integer, allocatable :: paired_gvec(:)
    
    ! The groups describing R-vector and G-vector operations.
    ! e.g. if Rvec(i)+Rvec(j)=Rvec(k) then rvec_group(i)*j=k.
    type(Group), allocatable :: rvector_group(:)
    type(Group), allocatable :: gvector_group(:)
  contains
    procedure, private  :: calculate_rvector_group
    procedure, private  :: calculate_gvector_group
  end type
  
  interface StructureData
    module procedure new_StructureData
  end interface
contains

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
  
  allocate( this%paired_rvec(this%sc_size), &
          & this%paired_gvec(this%sc_size), &
          & stat=ialloc); call err(ialloc)
  
  this%paired_rvec = 0
  do i=1,this%sc_size
    do j=1,i
      if (is_int( this%recip_supercell &
              & * (this%rvectors(i)+this%rvectors(j)))) then
        this%paired_rvec(i) = j
        this%paired_rvec(j) = i
      endif
    enddo
  enddo
  if (any(this%paired_rvec==0)) then
    call print_line(ERROR//': not all paired R-vectors found.')
    call err()
  endif
  
  this%paired_gvec = 0
  do i=1,this%sc_size
    do j=1,i
      if (is_int( transpose(this%recip_supercell) &
              & * (this%gvectors(i)+this%gvectors(j)))) then
        this%paired_gvec(i) = j
        this%paired_gvec(j) = i
      endif
    enddo
  enddo
  if (any(this%paired_gvec==0)) then
    call print_line('Error: not all paired G-vectors found.')
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
  
  this%symmetries = calculate_symmetries( basic_symmetries,   &
                                        & this%lattice,       &
                                        & this%recip_lattice, &
                                        & this%atoms)
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
    this%rvector_group(i) = operation
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
    this%gvector_group(i) = operation
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
function construct_supercell(structure,supercell_matrix,calculate_symmetries, &
   & symmetry_precision) result(supercell)
  use linear_algebra_module
  implicit none
  
  type(StructureData), intent(in)           :: structure
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
  logical :: calculate_symmetries_flag
  integer :: sc_size
  
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
  allocate( species(structure%no_atoms),           &
          & masses(structure%no_atoms),            &
          & positions(structure%no_atoms,sc_size), &
          & stat=ialloc); call err(ialloc)
  do prim=1,structure%no_atoms
    species(prim) = structure%atoms(prim)%species()
    masses(prim) = structure%atoms(prim)%mass()
    do rvec=1,sc_size
      positions(prim,rvec) = transpose(structure%lattice)                  &
                         & * ( structure%atoms(prim)%fractional_position() &
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
    supercell = StructureData(                                  &
       & supercell_matrix * structure%lattice,                  &
       & supercell_matrix,                                      &
       & calculate_unique_vectors(supercell_matrix, .false.),   &
       & calculate_unique_vectors( transpose(supercell_matrix), &
       &                     .true.),                           &
       & species,                                               &
       & masses,                                                &
       & positions,                                             &
       & symmetry_precision=symmetry_precision)
  else
    supercell = StructureData(                                  &
       & supercell_matrix * structure%lattice,                  &
       & supercell_matrix,                                      &
       & calculate_unique_vectors(supercell_matrix, .false.),   &
       & calculate_unique_vectors( transpose(supercell_matrix), &
       &                     .true.),                           &
       & species,                                               &
       & masses,                                                &
       & positions,                                             &
       & basic_symmetries_in=[BasicSymmetry::])
  endif
end function
end module
