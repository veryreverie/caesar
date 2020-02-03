! ======================================================================
! A density matrix, in the harmonic basis.
! ======================================================================
! N.B. While in general density matrices are not sparse, the operators
!    are sparse in the harmonic basis.
! The only use of the density matrix is to calculate the thermal expectation,
!    trace[DO], where D is the density matrix and O is the operator.
! As such, the density matrices are calculated and stored in the same sparse
!    representation as the Hamiltonian.
module density_matrix_module
  use common_module
  
  use coupled_states_module
  implicit none
  
  private
  
  public :: DensityMatrix
  
  ! Each state i couples with states in harmonic_couplings_(i)%ids(),
  !    where harmonic_couplings_ is a set of arrays stored in WavevectorBasis.
  ! If size(harmonic_couplings_(i)%ids())==n, then these states correspond to
  !    values(keys(i)+1:keys(i)+n).
  type, extends(NoDefaultConstructor) :: DensityMatrix
    type(FractionVector)  :: wavevector
    integer,  allocatable :: keys(:)
    real(dp), allocatable :: values(:)
  end type
end module
