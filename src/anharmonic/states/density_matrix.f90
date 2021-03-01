! ======================================================================
! A density matrix, in the harmonic basis.
! ======================================================================
! N.B. While in general density matrices are not sparse, the operators
!    are sparse in the harmonic basis.
! The only use of the density matrix is to calculate the thermal expectation,
!    trace[DO], where D is the density matrix and O is the operator.
! As such, the density matrices are calculated and stored in the same sparse
!    representation as the Hamiltonian.
module caesar_density_matrix_module
  use caesar_common_module
  
  use caesar_anharmonic_common_module
  
  use caesar_wavevector_state_module
  use caesar_coupled_states_module
  implicit none
  
  private
  
  public :: DensityMatrix
  
  ! Each state i couples with states in harmonic_couplings_(i)%ids(),
  !    where harmonic_couplings_ is a set of arrays stored in WavevectorBasis.
  ! If size(harmonic_couplings_(i)%ids())==n, then these states correspond to
  !    values(keys(i)+1:keys(i)+n).
  type, extends(NoDefaultConstructor) :: DensityMatrix
    type(FractionVector)  :: wavevector
    integer,  allocatable :: state_ids(:)
    integer,  allocatable :: bra_ids(:)
    integer,  allocatable :: ket_ids(:)
    real(dp), allocatable :: values(:)
  end type
  
  interface DensityMatrix
    ! Generate the density matrices associated with a WavevectorStates.
    ! This function has been somewhat optimised,
    !    as it takes up a large proportion of the total runtime.
    module function new_DensityMatrix(wavevector,couplings,states,weights) &
       & result(output) 
      type(FractionVector),  intent(in) :: wavevector
      type(CoupledStates),   intent(in) :: couplings(:)
      type(WavevectorState), intent(in) :: states(:)
      real(dp),              intent(in) :: weights(:)
      type(DensityMatrix)               :: output
    end function
  end interface
end module
