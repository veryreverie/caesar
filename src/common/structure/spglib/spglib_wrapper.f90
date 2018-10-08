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
  public :: SPGLIB_LINKED
  
  logical, parameter :: SPGLIB_LINKED = .true.
  
  interface SpglibSymmetries
    module procedure new_SpglibSymmetries_calculated
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
    subroutine spglib_retrieve_symmetries(spg_dataset_pointer, &
       & spacegroup_number,hall_number,international_symbol,hall_symbol, &
       & choice,transformation,origin_shift,n_operations,tensors, &
       & translations,n_atoms,pointgroup_symbol) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      
      type(c_ptr),            intent(in)  :: spg_dataset_pointer
      integer(kind=c_int),    intent(out) :: spacegroup_number
      integer(kind=c_int),    intent(out) :: hall_number
      character(kind=c_char), intent(out) :: international_symbol(11)
      character(kind=c_char), intent(out) :: hall_symbol(17)
      character(kind=c_char), intent(out) :: choice(6)
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
contains

function new_SpglibSymmetries_calculated(lattice,atoms,symmetry_precision) &
   & result(this)
  use, intrinsic :: iso_c_binding
  implicit none
  
  type(RealMatrix), intent(in) :: lattice
  type(AtomData),   intent(in) :: atoms(:)
  real(dp),         intent(in) :: symmetry_precision
  type(SpglibSymmetries)       :: this
  
  type(String) :: species
  integer      :: no_atom_types
  
  ! Spglib input variables.
  integer,  allocatable :: atom_types(:)
  real(dp), allocatable :: positions(:,:)
  
  ! Spglib output variables.
  integer               :: spacegroup_number
  integer               :: hall_number
  character(11)         :: international_symbol
  character(17)         :: hall_symbol
  character(6)          :: choice
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
  
  ! Temporary variables.
  integer :: i,j,ialloc
  
  ! Calculate spglib_atom_types array.
  ! spglib_atom_types(i) == spglib_atom_types(j) if and only if
  !    atoms(i)%species() == atoms(j)%species().
  allocate(atom_types(size(atoms)), stat=ialloc); call err(ialloc)
  no_atom_types = 0
  do i=1,size(atoms)
    species = atoms(i)%species()
    j = first(atoms%species()==species)
    if (j==i) then
      no_atom_types = no_atom_types + 1
      atom_types(i) = no_atom_types
    else
      atom_types(i) = atom_types(j)
    endif
  enddo
  
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
                                 & hall_number,          &
                                 & international_symbol, &
                                 & hall_symbol,          &
                                 & choice,               &
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
  
  this = SpglibSymmetries( spacegroup_number,         &
                         & hall_number,               &
                         & str(international_symbol), &
                         & str(hall_symbol),          &
                         & str(choice),               &
                         & output_transformation,     &
                         & output_origin_shift,       &
                         & n_operations,              &
                         & output_tensors,            &
                         & output_translations,       &
                         & n_atoms,                   &
                         & str(pointgroup_symbol)     )
end function
end module
