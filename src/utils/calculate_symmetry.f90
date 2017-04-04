! ----------------------------------------------------------------------
! Calculates the symmetries of the structure.
! ----------------------------------------------------------------------
module calculate_symmetry_module
  use iso_c_binding
  implicit none
  
  interface calculate_symmetry
    module procedure calculate_symmetry_old
    module procedure calculate_symmetry_new
  end interface
  
  ! spglib interface.
  interface
    subroutine spglib_get_no_symmetries( &
       & lattice,position,types,num_atom,symprec, &
       & n_operations) &
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
      integer(kind=c_int), intent(out) :: n_operations
    end subroutine
  end interface
  
  interface
    subroutine spglib_get_symmetries( &
       & lattice,position,types,num_atom,symprec, &
       & rotations,translations) &
       & bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      
      ! Inputs.
      real(kind=c_double), intent(in) :: lattice(3,3)
      real(kind=c_double), intent(in) :: position(3,*)
      integer(kind=c_int), intent(in) :: types(*)
      integer(kind=c_int), intent(in) :: num_atom
      real(kind=c_double), intent(in) :: symprec
      
      ! Outputs.
      integer(kind=c_int), intent(out) :: rotations(3,3,*)
      real(kind=c_double), intent(out) :: translations(3,*)
    end subroutine
  end interface
contains

subroutine calculate_symmetry_new(this)
  use string_module
  use file_module
  use structure_module
  implicit none
  
  type(StructureData), intent(inout) :: this
  
  ! spglib inputs
  real(dp), parameter :: symmetry_precision = 1.0e-2_dp
  integer             :: atom_types(this%no_atoms)
  
  ! spglib outputs
  integer               :: no_symmetries
  integer,  allocatable :: rotations(:,:,:)
  real(dp), allocatable :: translations(:,:)
  
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
  
  call spglib_get_no_symmetries( transpose(this%lattice), &
                               & this%atoms,              &
                               & atom_types,              &
                               & this%no_atoms,           &
                               & symmetry_precision,      &
                               & no_symmetries)
  
  allocate( rotations(3,3,no_symmetries),  &
          & translations(3,no_symmetries), &
          & stat=ialloc); call err(ialloc)
  
  call spglib_get_symmetries( transpose(this%lattice), &
                            & this%atoms,              &
                            & atom_types,              &
                            & this%no_atoms,           &
                            & symmetry_precision,      &
                            & rotations,               &
                            & translations)
end subroutine

subroutine calculate_symmetry_old(this,temp_cell_filename,temp_symm_filename)
  use string_module
  use file_module
  use structure_module
  use structure_to_dft_module
  implicit none
  
  type(StructureData), intent(inout) :: this
  type(String),        intent(in)    :: temp_cell_filename
  type(String),        intent(in)    :: temp_symm_filename
  
  call structure_to_dft( dft_code        = str('castep'), &
                       & structure_sc    = this,          &
                       & output_filename = temp_cell_filename )
  
  call system_call( &
     & 'cellsym --symmetry '//temp_cell_filename//' > '//temp_symm_filename)
  
  call read_symmetry_file(this, temp_symm_filename)
  
  call system_call('rm '//temp_cell_filename)
  call system_call('rm '//temp_symm_filename)
end subroutine
end module
