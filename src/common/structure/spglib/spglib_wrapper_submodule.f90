submodule (caesar_spglib_wrapper_module) caesar_spglib_wrapper_submodule
  use caesar_spglib_module
contains

module procedure new_SpglibSymmetries_calculated
  use, intrinsic :: iso_c_binding
  
  ! Spglib input variables.
  real(dp)              :: real_lattice(3,3)
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
  real_lattice = dble(lattice)
  call spglib_calculate_symmetries( real_lattice,        &
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
end procedure

module procedure snap_to_symmetry
  real(dp)              :: old_lattice(3,3)
  real(dp), allocatable :: old_positions(:,:)
  integer,  allocatable :: old_types(:)
  
  real(dp)              :: spg_lattice(3,3)
  real(dp), allocatable :: spg_positions(:,:)
  integer,  allocatable :: spg_types(:)
  
  type(RealMatrix) :: snapped_lattice
  type(RealMatrix) :: unsnapped_lattice
  type(RealMatrix) :: real_transformation
  type(IntMatrix)  :: transformation
  
  type(RealMatrix) :: new_lattice
  
  type(RealVector), allocatable :: old_frac_positions(:)
  type(RealVector), allocatable :: new_frac_positions(:)
  
  type(RealVector), allocatable :: new_positions(:)
  type(BasicAtom),  allocatable :: output_atoms(:)
  
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
  
  ! Call Spglib such that the lattice vectors are also snapped.
  ! N.B. this will in general rotate the structure by an arbitrary amount.
  ! N.B. this may also in general translate the structure by an arbitrary
  !    amount.
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
  
  ! Construct (transformation . output), i.e. the lattice which has
  !    the lengths and angles of the snapped lattice,
  !    but the orientation of the un-snapped lattice.
  new_lattice = standardize_lattice(snapped_lattice) &
            & * standard_transformation(unsnapped_lattice)
  
  ! Construct the 'transformation' matrix.
  real_transformation = transpose(lattice * invert(unsnapped_lattice))
  call check_int(real_transformation, 'spglib transformation matrix')
  transformation = nint(real_transformation)
  if (abs(determinant(transformation))/=1) then
    call print_line(ERROR//': Spglib transformation matrix does not have &
       &determinant 1.')
    call print_line('This may be because the input structure is not a &
       &primitive cell.')
    call print_line('If the structure is known to be a primitive cell of &
       &the system, try lowering symmetry_precision.')
    call print_line('N.B. running phonon calculations on a system with a unit &
       &cell which is x*y*z primitive cells is equivalent to running the &
       &calculation on that primitive cell but with an x*y*z q-point grid. &
       &Please increase the q-point grid rather than running calculations on &
       &a supercell.')
    call err()
  endif
  
  ! Construct the old atomic co-ordinates in the fractional co-ordinates of the
  !    new lattice.
  old_frac_positions = [(    invert(transpose(new_lattice)) &
                        &  * transpose(lattice)             &
                        &  * vec(old_positions(:,i)),       &
                        & i=1,                              &
                        & size(atoms)                       )]
  
  ! Translate the entire cell by three translations:
  !    - The first ensures that new_frac_positions(i) and old_frac_positions(i)
  !         are approximately separated by a lattice vector.
  !    - The second calculates that lattice vector, and subtracts it from
  !         new_positions(i), such that it approximately equals
  !         old_frac_positions(i).
  !    - The third moves the new frac positions to be as close as possible to
  !         the old frac positions.
  new_frac_positions = [(   new_frac_positions(i)                          &
                        & - (new_frac_positions(1)-old_frac_positions(1)), &
                        & i=1,                                             &
                        & size(atoms)                                      )]
  
  new_frac_positions = [(   new_frac_positions(i)          &
                        & - nint( new_frac_positions(i)    &
                        &       - old_frac_positions(i) ), &
                        & i=1,                             &
                        & size(atoms)                      )]
  
  new_frac_positions = [(   new_frac_positions(i)                      &
                        & - sum(new_frac_positions-old_frac_positions) &
                        & / size(atoms),                               &
                        & i=1,                                         &
                        & size(atoms)                                  )]
  
  ! Transform the lattice and atoms by this transformation, to return to the
  !    original cell.
  new_positions = [( transpose(new_lattice)*new_frac_positions(i), &
                   & i=1,                                          &
                   & size(atoms)                                   )]
  
  output_atoms = BasicAtom(atoms%species(), atoms%mass(), new_positions)
  output = BasicStructure(                                     &
     & lattice_matrix = transpose(transformation)*new_lattice, &
     & atoms          = output_atoms                           )
  
  ! Calculate changes.
  max_lattice_change = 0
  do i=1,3
    max_lattice_change = max(                                 &
       & max_lattice_change,                                  &
       & l2_norm(lattice%row(i)-output%lattice_matrix%row(i)) )
  enddo
  max_atom_change =                                    &
     & maxval(l2_norm( atoms%cartesian_position()      &
     &               - output%atoms%cartesian_position ))
  call print_line('Snapping to symmetry.')
  call print_line('Maximum lattice vector change  : '// &
                 & max_lattice_change//' (bohr)')
  call print_line('Maximum atomic position change : '// &
                 & max_atom_change//' (bohr)')
end procedure

module procedure standard_transformation
  output = transpose(standardize_lattice(lattice)) &
       & * transpose(invert(lattice))
end procedure

module procedure standardize_lattice
  real(dp) :: dot_products(3,3)
  real(dp) :: ax
  real(dp) :: bx, by
  real(dp) :: cx, cy, cz
  
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
              & [3,3]                    )
end procedure

module procedure calculate_spglib_atom_types
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
end procedure
end submodule
