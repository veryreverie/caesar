! ======================================================================
! A polynomial representation of a potential.
! ======================================================================
module polynomial_potential_module
  use common_module
  
  use anharmonic_common_module
  
  use basis_function_module
  use basis_functions_module
  use sampling_points_module
  implicit none
  
  private
  
  public :: PolynomialPotential
  
  type, extends(PotentialData) :: PolynomialPotential
    integer                           :: potential_expansion_order
    type(BasisFunctions), allocatable :: basis_functions(:)
    type(SamplingPoints), allocatable :: sampling_points(:)
  contains
    procedure, public :: generate_sampling_points => &
       & generate_sampling_points_PolynomialPotential
  end type
  
  interface PolynomialPotential
    module procedure new_PolynomialPotential
  end interface
contains

! Constructor.
function new_PolynomialPotential(potential_expansion_order) result(this)
  implicit none
  
  integer, intent(in)       :: potential_expansion_order
  type(PolynomialPotential) :: this
  
  this%potential_expansion_order = potential_expansion_order
end function

! Generate sampling points.
subroutine generate_sampling_points_PolynomialPotential(this,coupled_subspaces)
  implicit none
  
  class(PolynomialPotential), intent(inout) :: this
  type(CoupledSubspaces),     intent(in)    :: coupled_subspaces(:)
  
  type(SubspaceMonomial), allocatable :: subspace_monomials(:)
  
  integer :: i
  
  do i=1,size(coupled_subspaces)
    ! Generate the set of subspace monomials corresponding to the subspace
    !    coupling.
    ! e.g. the coupling [1,2] might have monomials [1,2], [1,1,2] and [1,2,2].
    subspace_monomials = generate_subspace_monomials( &
                              & coupled_subspaces(i), &
                              & this%potential_expansion_order)
    
  enddo
end subroutine
end module
