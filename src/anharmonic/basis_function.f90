! ======================================================================
! Generates basis functions.
! ======================================================================
module basis_function_module
  use constants_module, only : dp
  use string_module
  use io_module
  implicit none
  
  private
  
  public :: generate_basis_functions
contains

function generate_basis_functions(coupling,normal_modes,qpoints, &
   & subspaces,symmetries) result(output)
  use coupling_module
  use normal_mode_module
  use qpoints_module
  use integer_arrays_module
  use polynomial_module
  use logic_module
  use degeneracy_module
  implicit none
  
  type(CoupledSubspaces), intent(in) :: coupling
  type(ComplexMode),      intent(in) :: normal_modes(:)
  type(QpointData),       intent(in) :: qpoints(:)
  type(DegenerateModes),  intent(in) :: subspaces(:)
  type(ComplexMatrix),    intent(in) :: symmetries(:,:)
  type(Polynomial), allocatable      :: output(:)
  
  type(DegenerateModes), allocatable :: coupled_subspaces(:)
  
  integer :: i,j,ialloc
  
  do i=1,size(coupling)
    coupled_subspaces = coupling%coupled_subspaces(subspaces)
  enddo
end function

end module
