! ======================================================================
! Uses spglib to calculate symmetry operations.
! ======================================================================
module caesar_spglib_wrapper_module
  use caesar_utils_module
  
  use caesar_atom_module
  
  use caesar_spglib_symmetries_module
  implicit none
  
  private
  
  public :: SpglibSymmetries
  public :: snap_to_symmetry
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
  
  interface SpglibSymmetries
    module function new_SpglibSymmetries_calculated(lattice,atoms, &
       & symmetry_precision) result(this) 
      use, intrinsic :: iso_c_binding
      implicit none
      
      type(RealMatrix), intent(in) :: lattice
      type(AtomData),   intent(in) :: atoms(:)
      real(dp),         intent(in) :: symmetry_precision
      type(SpglibSymmetries)       :: this
    end function
  end interface
  
  interface
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
    module function snap_to_symmetry(lattice,atoms,symmetry_precision) &
       & result(output) 
      type(RealMatrix), intent(in) :: lattice
      type(AtomData),   intent(in) :: atoms(:)
      real(dp),         intent(in) :: symmetry_precision
      type(BasicStructure)         :: output
    end function
  end interface
  
  interface
    ! Construct the orthonormal transformation which transforms the lattice to
    !    the standard lattice, as defined by standardize_lattice.
    module function standard_transformation(lattice) result(output) 
      type(RealMatrix), intent(in) :: lattice
      type(RealMatrix)             :: output
    end function
  end interface
  
  interface
    ! Rotate a lattice so that a is along x and b is in the x-y plane.
    ! N.B. this will always form a right-handed lattice,
    !    even if the input is left-handed.
    module function standardize_lattice(lattice) result(output) 
      type(RealMatrix), intent(in) :: lattice
      type(RealMatrix)             :: output
    end function
  end interface
  
  interface
    ! Calculate spglib_atom_types array.
    ! spglib_atom_types(i) == spglib_atom_types(j) if and only if
    !    atoms(i)%species() == atoms(j)%species().
    module function calculate_spglib_atom_types(atoms) result(output) 
      type(AtomData), intent(in) :: atoms(:)
      integer, allocatable       :: output(:)
    end function
  end interface
end module
