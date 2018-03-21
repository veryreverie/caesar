! ======================================================================
! Uses spglib to calculate symmetry operations.
! ======================================================================
module calculate_symmetry_submodule
  use utils_module
  
  use basic_structure_submodule
  use atom_submodule
  implicit none
  
  private
  
  public :: calculate_basic_symmetries
  
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
  
  ! spglib inputs.
  integer, allocatable :: atom_types(:)
  
  ! spglib data.
  logical(kind=c_bool) :: spglib_success
  type(c_ptr)          :: spg_dataset_pointer
  integer              :: no_symmetries
  
  ! Temporary matrix storage, for passing to and from C.
  real(dp), allocatable :: positions(:,:)
  integer,  allocatable :: rotations(:,:,:)
  real(dp), allocatable :: translations(:,:)
  
  integer :: atom_type
  
  ! Temporary variables.
  integer :: i,j,ialloc
  
  allocate(atom_types(size(atoms)), stat=ialloc); call err(ialloc)
  atom_type = 0
  do_i : do i=1,size(atoms)
    do j=1,i-1
      if (atoms(i)%species()==atoms(j)%species()) then
        atom_types(i) = atom_types(j)
        cycle do_i
      endif
    enddo
    atom_types(i) = atom_type
    atom_type = atom_type+1
  enddo do_i
  
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
          & output(no_symmetries), &
          & stat=ialloc); call err(ialloc)
  
  ! Retrieve symmetries into allocated space.
  call spglib_retrieve_symmetries( spg_dataset_pointer, &
                                 & rotations,      &
                                 & translations)
  do i=1,no_symmetries
    output(i)%rotation = rotations(:,:,i)
    output(i)%translation = translations(:,i)
  enddo
  
  ! Deallocate C memory.
  call drop_spg_dataset(spg_dataset_pointer)
end function
end module
