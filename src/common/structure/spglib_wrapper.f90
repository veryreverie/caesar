! ======================================================================
! Uses spglib to calculate symmetry operations.
! ======================================================================
module spglib_wrapper_submodule
  use utils_module
  
  use basic_symmetry_submodule
  use basic_structure_submodule
  use atom_submodule
  implicit none
  
  private
  
  public :: calculate_basic_symmetries
  
  public :: SPGLIB_LINKED
  
  logical, parameter :: SPGLIB_LINKED = .true.
  
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
    subroutine spglib_retrieve_symmetries(spg_dataset_pointer,tensors, &
       & translations) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      
      ! Inputs.
      type(c_ptr), intent(in) :: spg_dataset_pointer
      
      ! Outputs.
      integer(kind=c_int), intent(out) :: tensors(3,3,*)
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
! Calls spglib, and calculates all symmetry operations.
! ----------------------------------------------------------------------
function calculate_basic_symmetries(lattice,atoms,symmetry_precision) &
   & result(output)
  use, intrinsic :: iso_c_binding
  implicit none
  
  type(RealMatrix), intent(in)     :: lattice
  type(AtomData),   intent(in)     :: atoms(:)
  real(dp),         intent(in)     :: symmetry_precision
  type(BasicSymmetry), allocatable :: output(:)
  
  type(String) :: species
  integer      :: no_atom_types
  
  ! spglib variables.
  integer,  allocatable :: spglib_atom_types(:)
  real(dp), allocatable :: spglib_positions(:,:)
  integer,  allocatable :: spglib_tensors(:,:,:)
  real(dp), allocatable :: spglib_translations(:,:)
  
  ! spglib data.
  logical(kind=c_bool) :: spglib_success
  type(c_ptr)          :: spg_dataset_pointer
  integer              :: no_symmetries
  
  ! Output variables.
  type(IntMatrix),  allocatable :: tensors(:)
  type(RealVector), allocatable :: translations(:)
  type(Group),      allocatable :: atom_groups(:)
  type(IntVector),  allocatable :: rvectors(:,:)
  
  ! Variables for calculating atom groups and R-vectors.
  type(RealVector)              :: transformed_position
  type(RealVector), allocatable :: kj_displacements(:)
  type(IntVector),  allocatable :: kj_rvectors(:)
  real(dp),         allocatable :: kj_distances(:)
  integer,          allocatable :: atom_group(:)
  
  ! Temporary variables.
  integer :: i,j,k,ialloc
  
  ! --------------------------------------------------
  ! Calculate tensors and translations.
  ! --------------------------------------------------
  
  ! Calculate spglib_atom_types array.
  ! spglib_atom_types(i) == spglib_atom_types(j) if and only if
  !    atoms(i)%species() == atoms(j)%species().
  allocate(spglib_atom_types(size(atoms)), stat=ialloc); call err(ialloc)
  no_atom_types = 0
  do i=1,size(atoms)
    species = atoms(i)%species()
    j = first(atoms%species()==species)
    if (j==i) then
      no_atom_types = no_atom_types + 1
      spglib_atom_types(i) = no_atom_types
    else
      spglib_atom_types(i) = spglib_atom_types(j)
    endif
  enddo
  
  ! Calculate symmetries, without space to store them.
  allocate(spglib_positions(3,size(atoms)), stat=ialloc); call err(ialloc)
  do i=1,size(atoms)
    spglib_positions(:,i) = dble(atoms(i)%fractional_position())
  enddo
  call spglib_calculate_symmetries( dble(lattice),       &
                                  & spglib_positions,    &
                                  & spglib_atom_types,   &
                                  & size(atoms),         &
                                  & symmetry_precision,  &
                                  & spglib_success,      &
                                  & spg_dataset_pointer, &
                                  & no_symmetries        )
  
  ! Check for errors.
  if (.not. spglib_success) then
    call print_line('Error: spglib symmetry finder failed. &
       &Try increasing symmetry_precision.')
    call err()
  endif
  
  ! Allocate space for symmetries.
  allocate( spglib_tensors(3,3,no_symmetries),    &
          & spglib_translations(3,no_symmetries), &
          & stat=ialloc); call err(ialloc)
  
  ! Retrieve symmetries into allocated space.
  call spglib_retrieve_symmetries( spg_dataset_pointer, &
                                 & spglib_tensors,      &
                                 & spglib_translations  )
  
  ! Deallocate C memory.
  call drop_spg_dataset(spg_dataset_pointer)
  
  ! Construct tensor and translation arrays.
  allocate( tensors(no_symmetries),      &
          & translations(no_symmetries), &
          & stat=ialloc); call err(ialloc)
  do i=1,no_symmetries
    tensors(i) = mat(spglib_tensors(:,:,i))
    translations(i) = vec(spglib_translations(:,i))
  enddo
  
  ! --------------------------------------------------
  ! Calculate atom group and R-vectors.
  ! --------------------------------------------------
  ! The symmetry S transforms atomic equilibrium positions, {r_i} in
  !    fractional co-ordinates as:
  ! r_i -> r_j + R
  !
  ! atom_group * i = j
  ! rvectors(i) = R
  allocate( atom_groups(no_symmetries),           &
          & rvectors(size(atoms), no_symmetries), &
          & kj_displacements(size(atoms)),        &
          & kj_rvectors(size(atoms)),             &
          & kj_distances(size(atoms)),            &
          & atom_group(size(atoms)),              &
          & stat=ialloc); call err(ialloc)
  do i=1,no_symmetries
    do j=1,size(atoms)
      ! Calculate the position of the transformed atom.
      transformed_position = tensors(i) * atoms(j)%fractional_position() &
                         & + translations(i)
      
      do k=1,size(atoms)
        ! Calculate the displacement from atom k to the transfomed atom j.
        kj_displacements(k) = transformed_position &
                          & - atoms(k)%fractional_position()
        ! Calculate the nearest R-vector to the displacement.
        kj_rvectors(k) = nint(kj_displacements(k))
        ! Calculate the distance between the displacement and the nearest
        !    R-vector.
        kj_distances(k) = l2_norm(kj_displacements(k)-kj_rvectors(k))
      enddo
      
      ! Check that exactly one k->j distance is small.
      if (count(kj_distances<symmetry_precision)/=1) then
        call print_line(ERROR//': Symmetries do not map atoms as expected.')
        call err()
      endif
      
      ! Identify the atom which symmetry i maps atom j to.
      k = minloc(kj_distances,1)
      
      if (atoms(k)%species()/=atoms(j)%species()) then
        call print_line(ERROR//': symmetry operation maps between atoms of &
           &different species.')
        call err()
      endif
      
      atom_group(j) = atoms(k)%id()
      rvectors(j,i) = kj_rvectors(k)
    enddo
    atom_groups(i) = Group(atom_group)
  enddo
  
  ! --------------------------------------------------
  ! Construct output.
  ! --------------------------------------------------
  
  allocate(output(no_symmetries), stat=ialloc); call err(ialloc)
  do i=1,no_symmetries
    output(i) = BasicSymmetry( id          = i,               &
                             & tensor      = tensors(i),      &
                             & translation = translations(i), &
                             & atom_group  = atom_groups(i),  &
                             & rvectors    = rvectors(:,i)    )
  enddo
end function
end module
