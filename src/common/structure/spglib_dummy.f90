! ======================================================================
! A dummy module for when spglib is not being linked.
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
  
  logical, parameter :: SPGLIB_LINKED = .false.
contains

! ----------------------------------------------------------------------
! Replaces calculate_basic_symmetries in spglib_wrapper.f90.
! ----------------------------------------------------------------------
function calculate_basic_symmetries(lattice,atoms,symmetry_precision) &
   & result(output)
  implicit none
  
  type(RealMatrix), intent(in)     :: lattice
  type(AtomData),   intent(in)     :: atoms(:)
  real(dp),         intent(in)     :: symmetry_precision
  type(BasicSymmetry), allocatable :: output(:)
  
  call print_line(ERROR//': Cannot calculate symmetries because Caesar has &
     &not been linked against spglib. Please use the CMake flag &
     &-DLINK_TO_SPGLIB:LOGICAL=true to link against spglib.')
  stop
end function
end module
