! ----------------------------------------------------------------------
! Calculates the symmetries of the structure.
! ----------------------------------------------------------------------
module calculate_symmetry_module
  use iso_c_binding
  implicit none
  
  ! --------------------------------------------------
  ! spglib interface.
  ! --------------------------------------------------
  ! spglib_calculate_symmetries calculates symmetry operations, but has
  !    nowhere to store them. Space should be allocated, and passed into
  !    spglib_retrieve_symmetries, which will copy the data from
  !    C-style malloc memory to Fortran allocatable memory.
  interface
    subroutine spglib_calculate_symmetries( &
       & lattice,position,types,num_atom,symprec, &
       & spg_dataset_pointer,n_operations) &
       & bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      
      ! Inputs.
      real(kind=c_double), intent(in) :: lattice(3,3)
      real(kind=c_double), intent(in) :: position(3,*)
      integer(kind=c_int), intent(in) :: types(*)
      integer(kind=c_int), intent(in) :: num_atom
      real(kind=c_double), intent(in) :: symprec
      
      ! Output.
      type(c_ptr),         intent(out) :: spg_dataset_pointer
      integer(kind=c_int), intent(out) :: n_operations
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
contains

subroutine calculate_symmetry(this)
  use string_module
  use file_module
  use structure_module
  implicit none
  
  type(StructureData), intent(inout) :: this
  
  ! spglib inputs.
  real(dp), parameter :: symmetry_precision = 1.0e-1_dp
  integer             :: atom_types(this%no_atoms)
  
  ! spglib data.
  type(c_ptr) :: spg_dataset_pointer
  
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
  call spglib_calculate_symmetries( transpose(this%lattice), &
                                  & matmul(this%recip_lattice,this%atoms), &
                                  & atom_types,              &
                                  & this%no_atoms,           &
                                  & symmetry_precision,      &
                                  & spg_dataset_pointer,     &
                                  & this%no_symmetries)
  
  ! Allocate space for symmetries.
  allocate( this%rotations(3,3,this%no_symmetries),  &
          & this%offsets(3,this%no_symmetries),      &
          & stat=ialloc); call err(ialloc)
  
  ! Retrieve symmetries into allocated space.
  call spglib_retrieve_symmetries( spg_dataset_pointer, &
                                 & this%rotations,      &
                                 & this%offsets)
end subroutine
end module
