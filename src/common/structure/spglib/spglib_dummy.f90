! ======================================================================
! A dummy module for when spglib is not being linked.
! ======================================================================
module spglib_wrapper_module
  use utils_module
  
  use atom_module
  
  use spglib_symmetries_module
  implicit none
  
  private
  
  public :: SpglibSymmetries
  public :: SPGLIB_LINKED
  
  logical, parameter :: SPGLIB_LINKED = .false.
  
  interface SpglibSymmetries
    module procedure new_SpglibSymmetries_calculated
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
  
  call print_line(ERROR//': Cannot calculate symmetries because Caesar has &
     &not been linked against spglib. Please use the CMake flag &
     &-DLINK_TO_SPGLIB:LOGICAL=true to link against spglib.')
  stop
end function
end module
