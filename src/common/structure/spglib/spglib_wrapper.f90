! ======================================================================
! Uses spglib to calculate symmetry operations.
! ======================================================================
module spglib_wrapper_module
  use utils_module
  
  use atom_module
  
  use spglib_symmetries_module
  implicit none
  
  private
  
  public :: SpglibSymmetries
  public :: snap_to_symmetry
  public :: SPGLIB_LINKED
  
  logical, parameter :: SPGLIB_LINKED = .true.
  
  interface SpglibSymmetries
    module procedure new_SpglibSymmetries_calculated
  end interface
  
  interface BasicStructure
    module procedure new_BasicStructure_spglib
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
      
      real(kind=c_double),  intent(in)  :: lattice(3,3)
      real(kind=c_double),  intent(in)  :: position(3,*)
      integer(kind=c_int),  intent(in)  :: types(*)
      integer(kind=c_int),  intent(in)  :: num_atom
      real(kind=c_double),  intent(in)  :: symprec
      logical(kind=c_bool), intent(out) :: success
      type(c_ptr),          intent(out) :: spg_dataset_pointer
      integer(kind=c_int),  intent(out) :: n_operations
    end subroutine
  end interface
  
  interface
    subroutine spglib_retrieve_symmetries(spg_dataset_pointer,               &
       & spacegroup_number,international_symbol,transformation,origin_shift, &
       & n_operations,tensors,translations,n_atoms,pointgroup_symbol) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      
      type(c_ptr),            intent(in)  :: spg_dataset_pointer
      integer(kind=c_int),    intent(out) :: spacegroup_number
      character(kind=c_char), intent(out) :: international_symbol(11)
      real(kind=c_double),    intent(out) :: transformation(3,3)
      real(kind=c_double),    intent(out) :: origin_shift(3)
      integer(kind=c_int),    intent(out) :: n_operations
      integer(kind=c_int),    intent(out) :: tensors(3,3,*)
      real(kind=c_double),    intent(out) :: translations(3,*)
      integer(kind=c_int),    intent(out) :: n_atoms
      character(kind=c_char), intent(out) :: pointgroup_symbol(6)
    end subroutine
  end interface
  
  interface
    subroutine drop_spg_dataset(spg_dataset_pointer) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      
      type(c_ptr), intent(inout) :: spg_dataset_pointer
    end subroutine
  end interface
  
  interface
    function spglib_standardize_cell(lattice,position,types,num_atom, &
       & to_primitive,no_idealize,symprec) result(output) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      
      real(kind=c_double), intent(inout) :: lattice(3,3)
      real(kind=c_double), intent(inout) :: position(3,*)
      integer(kind=c_int), intent(inout) :: types(*)
      integer(kind=c_int), intent(in)    :: num_atom
      integer(kind=c_int), intent(in)    :: to_primitive
      integer(kind=c_int), intent(in)    :: no_idealize
      real(kind=c_double), intent(in)    :: symprec
      integer(kind=c_int)                :: output
    end function
  end interface
contains

function new_SpglibSymmetries_calculated(lattice,atoms,symmetry_precision) &
   & result(this)
  use, intrinsic :: iso_c_binding
  implicit none
  
  type(RealMatrix), intent(in) :: lattice
  type(AtomData),   intent(in) :: atoms(:)
  real(dp),         intent(in) :: symmetry_precision
  type(SpglibSymmetries)       :: this
  
  ! Spglib input variables.
  integer,  allocatable :: atom_types(:)
  real(dp), allocatable :: positions(:,:)
  
  ! Spglib output variables.
  integer               :: spacegroup_number
  character(11)         :: international_symbol
  real(dp)              :: transformation(3,3)
  real(dp)              :: origin_shift(3)
  integer               :: n_operations
  integer,  allocatable :: tensors(:,:,:)
  real(dp), allocatable :: translations(:,:)
  integer               :: n_atoms
  character(6)          :: pointgroup_symbol
  
  ! spglib data.
  logical(kind=c_bool) :: spglib_success
  type(c_ptr)          :: spg_dataset_pointer
  
  ! Output variables.
  type(RealMatrix)              :: output_transformation
  type(Realvector)              :: output_origin_shift
  type(IntMatrix),  allocatable :: output_tensors(:)
  type(RealVector), allocatable :: output_translations(:)
  
  type(String) :: output_international_symbol
  type(String) :: output_pointgroup_symbol
  
  ! Temporary variables.
  integer :: i,ialloc
  
  ! Calculate spglib_atom_types array.
  atom_types = calculate_spglib_atom_types(atoms)
  
  ! Calculate symmetries, without space to store them.
  allocate(positions(3,size(atoms)), stat=ialloc); call err(ialloc)
  do i=1,size(atoms)
    positions(:,i) = dble(atoms(i)%fractional_position())
  enddo
  call spglib_calculate_symmetries( dble(lattice),       &
                                  & positions,           &
                                  & atom_types,          &
                                  & size(atoms),         &
                                  & symmetry_precision,  &
                                  & spglib_success,      &
                                  & spg_dataset_pointer, &
                                  & n_operations         )
  
  ! Check for errors.
  if (.not. spglib_success) then
    call print_line('Error: spglib symmetry finder failed. &
       &Try increasing symmetry_precision.')
    call err()
  endif
  
  ! Allocate space for symmetries.
  allocate( tensors(3,3,n_operations),    &
          & translations(3,n_operations), &
          & stat=ialloc); call err(ialloc)
  
  ! Retrieve symmetries into allocated space.
  call spglib_retrieve_symmetries( spg_dataset_pointer,  &
                                 & spacegroup_number,    &
                                 & international_symbol, &
                                 & transformation,       &
                                 & origin_shift,         &
                                 & n_operations,         &
                                 & tensors,              &
                                 & translations,         &
                                 & n_atoms,              &
                                 & pointgroup_symbol     )
  
  ! Deallocate C memory.
  call drop_spg_dataset(spg_dataset_pointer)
  
  ! Convert from arrays to vectors/matrices.
  output_transformation = mat(transformation)
  output_origin_shift = vec(origin_shift)
  output_tensors = [( mat(tensors(:,:,i)), i=1, n_operations )]
  output_translations = [( vec(translations(:,i)), i=1, n_operations )]
  
  ! Convert from character arrays to strings.
  output_international_symbol = parse_c_string(international_symbol)
  output_pointgroup_symbol = parse_c_string(pointgroup_symbol)
  
  ! Construct output.
  this = SpglibSymmetries( spacegroup_number,           &
                         & output_international_symbol, &
                         & output_transformation,       &
                         & output_origin_shift,         &
                         & n_operations,                &
                         & output_tensors,              &
                         & output_translations,         &
                         & n_atoms,                     &
                         & output_pointgroup_symbol     )
end function

! Snap a structure to symmetry.
! The input lattice is I.
! The output lattice is O.
! Spglib constructs two lattices:
!   unsnapped = transformation . input
!   snapped   = transformation . output . rotation
! Where 'rotation' is a rotation or improper rotation matrix,
!    and 'transformation' is an integer matrix with determinant +/- 1.
! transformation^T = input^{-1} . unsnapped^T
! -> transformation   =  unsnapped . input^{-T}
function snap_to_symmetry(lattice,atoms,symmetry_precision) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: lattice
  type(AtomData),   intent(in) :: atoms(:)
  real(dp),         intent(in) :: symmetry_precision
  type(BasicStructure)         :: output
  
  real(dp)              :: old_lattice(3,3)
  real(dp), allocatable :: old_positions(:,:)
  integer,  allocatable :: old_types(:)
  
  real(dp)              :: spg_lattice(3,3)
  real(dp), allocatable :: spg_positions(:,:)
  integer,  allocatable :: spg_types(:)
  
  type(RealMatrix) :: snapped_lattice
  type(RealMatrix) :: unsnapped_lattice
  type(RealMatrix) :: rotation
  type(RealMatrix) :: real_transformation
  type(IntMatrix)  :: transformation
  
  type(RealMatrix) :: new_lattice
  
  type(RealVector), allocatable :: old_frac_positions(:)
  type(RealVector), allocatable :: check_frac_positions(:)
  type(RealVector), allocatable :: new_frac_positions(:)
  integer,          allocatable :: mapping(:)
  
  real(dp) :: max_lattice_change
  real(dp) :: max_atom_change
  
  integer :: return_code
  
  integer :: i,ialloc
  
  ! Construct Spglib types.
  old_lattice = dble(lattice)
  allocate(old_positions(3,size(atoms)), stat=ialloc); call err(ialloc)
  do i=1,size(atoms)
    old_positions(:,i) = dble(atoms(i)%fractional_position())
  enddo
  old_types = calculate_spglib_atom_types(atoms)
  
  ! Call Spglib such that the atoms are snapped but the lattice vectors are
  !    unchanged (except by linear combinations).
  spg_lattice = old_lattice
  spg_positions = old_positions
  spg_types = old_types
  return_code = spglib_standardize_cell( lattice      = spg_lattice,       &
                                       & position     = spg_positions,     &
                                       & types        = spg_types,         &
                                       & num_atom     = size(atoms),       &
                                       & to_primitive = 1,                 &
                                       & no_idealize  = 1,                 &
                                       & symprec      = symmetry_precision )
  if (return_code==0) then
    call print_line(ERROR//': spglib standardize_cell failed.')
    call err()
  endif
  ! unsnapped = transformation . input.
  unsnapped_lattice = mat(spg_lattice)
  check_frac_positions = [(vec(spg_positions(:,i)), i=1, size(atoms))]
  
  ! Call Spglib such that the lattice vectors are also snapped.
  ! N.B. this will in general rotate the structure by an arbitrary amount.
  ! N.B. this may also invert and re-order the atomic positions.
  spg_lattice = old_lattice
  spg_positions = old_positions
  spg_types = old_types
  return_code = spglib_standardize_cell( lattice      = spg_lattice,       &
                                       & position     = spg_positions,     &
                                       & types        = spg_types,         &
                                       & num_atom     = size(atoms),       &
                                       & to_primitive = 1,                 &
                                       & no_idealize  = 0,                 &
                                       & symprec      = symmetry_precision )
  if (return_code==0) then
    call print_line(ERROR//': spglib standardize_cell failed.')
    call err()
  endif
  ! snapped = transformation . output . rotation.
  snapped_lattice    = mat(spg_lattice)
  new_frac_positions = [(vec(spg_positions(:,i)), i=1, size(atoms))]
  
  ! Check if the positions have been inverted by spglib,
  !    and undo the inversion if they have.
  if (   sum(l2_norm( new_frac_positions                              &
     &              - check_frac_positions                            &
     &              - nint( new_frac_positions                        &
     &                    - check_frac_positions ) ))                 &
     & > sum(l2_norm( new_frac_positions(size(atoms):1:-1)            &
     &              + check_frac_positions                            &
     &              - nint( new_frac_positions(size(atoms):1:-1)      &
     &                    + check_frac_positions                 ) )) ) then
    new_frac_positions = - new_frac_positions
  endif
  
  ! Construct (transformation . output), i.e. the lattice which has
  !    the lengths and angles of the snapped lattice,
  !    but the orientation of the un-snapped lattice.
  new_lattice = standardize_lattice(snapped_lattice) &
            & * standard_transformation(unsnapped_lattice)
  
  ! Construct the old atomic co-ordinates in the fractional co-ordinates of the
  !    new lattice.
  old_frac_positions = [(    invert(transpose(new_lattice)) &
                        &  * transpose(lattice)             &
                        &  * vec(old_positions(:,i)),       &
                        & i=1,                              &
                        & size(atoms)                       )]
  
  ! Construct the mapping from the new positions to the old positions.
  mapping = [(                                                               &
     & minloc( l2_norm( new_frac_positions-old_frac_positions(i)             &
     &                - nint(new_frac_positions-old_frac_positions(i)) ),    &
     &         1                                                          ), &
     & i=1,                                                                  &
     & size(atoms)                                                           )]
  if (size(mapping)/=size(set(mapping))) then
    call print_line(CODE_ERROR//': spglib atom mapping is not one-to-one.')
    call err()
  endif
  
  ! Translate the atomic positions by the correct lattice vectors,
  !    and transform the atomic positions into cartesian co-ordinates.
  new_frac_positions = [(   new_frac_positions(mapping(i))          &
                        & - nint( new_frac_positions(mapping(i))    &
                        &       - old_frac_positions(i)          ), &
                        & i=1,                                      &
                        & size(atoms)                               )]
  
  ! Construct the 'transformation' matrix.
  real_transformation = transpose(lattice * invert(unsnapped_lattice))
  call check_int(real_transformation, 'spglib transformation matrix')
  transformation = nint(real_transformation)
  if (determinant(transformation)/=1) then
    call print_line(ERROR//': Spglib transformation matrix does not have &
       &determinant 1.')
    call err()
  endif
  
  ! Transform the lattice and atoms by this transformation, to return to the
  !    original cell.
  output = BasicStructure(                                              &
     & atoms,                                                           &
     & transpose(transformation)*new_lattice,                           &
     & [(transpose(new_lattice)*new_frac_positions(i),i=1,size(atoms))] )
  
  ! Calculate changes.
  max_lattice_change = maxval([(                             &
     & l2_norm(lattice%row(i)-output%lattice_matrix%row(i)), &
     & i=1,                                                  &
     & 3                                                     )])
  max_atom_change =                                    &
     & maxval(l2_norm( atoms%cartesian_position()      &
     &               - output%atoms%cartesian_position ))
  call print_line('Snapping to symmetry.')
  call print_line('Maximum lattice vector change  : '// &
                 & max_lattice_change//' (bohr)')
  call print_line('Maximum atomic position change : '// &
                 & max_atom_change//' (bohr)')
end function

! Construct the orthonormal transformation which transforms the lattice to
!    the standard lattice, as defined by standardize_lattice.
function standard_transformation(lattice) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: lattice
  type(RealMatrix)             :: output
  
  output = transpose(standardize_lattice(lattice)) &
       & * transpose(invert(lattice))
end function

! Rotate a lattice so that a is along x and b is in the x-y plane.
! N.B. this will always form a right-handed lattice,
!    even if the input is left-handed.
function standardize_lattice(lattice) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: lattice
  type(RealMatrix)             :: output
  
  real(dp) :: dot_products(3,3)
  real(dp) :: ax
  real(dp) :: bx, by
  real(dp) :: cx, cy, cz
  real(dp) :: norms(3)
  real(dp) :: matrix(3,3)
  
  integer :: i
  
  !                   ( a.a a.b a.c )
  ! Construct L.T^T = ( b.a b.b b.c )
  !                   ( c.a c.b c.c )
  dot_products = dble(lattice*transpose(lattice))
  
  ! a = (a.x,0,0) -> a.x = |a| = sqrt(a.a).
  ax = sqrt(dot_products(1,1))
  ! b.x = b.a / |a|.
  bx = dot_products(1,2) / ax
  ! b.y = sqrt(|b|^2 - (b.x)^2)
  by = sqrt(dot_products(2,2) - bx**2)
  ! c.x = c.x / |a|.
  cx = dot_products(1,3) / ax
  ! c.b = c.x*b.x+c.y*b.y -> c.y = (c.b-c.x*b.x)/b.y.
  cy = (dot_products(2,3)-cx*bx) / by
  ! c.z = sqrt(|c|^2 - (c.x)^2 - (c.y)^2)
  cz = sqrt(dot_products(3,3) - cx**2 - cy**2)
  
  output = mat( [ ax, 0.0_dp, 0.0_dp,    &
              &   bx, by    , 0.0_dp,    &
              &   cx, cy    , cz      ], &
              & 3,3                      )
end function

function new_BasicStructure_spglib(atoms,lattice,positions) result(output)
  implicit none
  
  type(AtomData),   intent(in) :: atoms(:)
  type(RealMatrix), intent(in) :: lattice
  type(RealVector), intent(in) :: positions(:)
  type(BasicStructure)         :: output
  
  type(BasicAtom),  allocatable :: output_atoms(:)
  
  output_atoms = BasicAtom(atoms%species(), atoms%mass(), positions)
  output = BasicStructure( lattice_matrix = lattice,     &
                         & atoms          = output_atoms )
end function

! Calculate spglib_atom_types array.
! spglib_atom_types(i) == spglib_atom_types(j) if and only if
!    atoms(i)%species() == atoms(j)%species().
function calculate_spglib_atom_types(atoms) result(output)
  implicit none
  
  type(AtomData), intent(in) :: atoms(:)
  integer, allocatable       :: output(:)
  
  type(String) :: species
  integer      :: no_atom_types
  
  integer :: i,j,ialloc

  allocate(output(size(atoms)), stat=ialloc); call err(ialloc)
  no_atom_types = 0
  do i=1,size(atoms)
    species = atoms(i)%species()
    j = first(atoms%species()==species)
    if (j==i) then
      no_atom_types = no_atom_types + 1
      output(i) = no_atom_types
    else
      output(i) = output(j)
    endif
  enddo
end function
end module
