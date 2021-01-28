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
    module procedure new_DensityMatrix
  end interface
contains

! Generate the density matrices associated with a WavevectorStates.
! This function has been somewhat optimised, as it takes up a large proportion
!    of the total runtime.
function new_DensityMatrix(wavevector,couplings,states,weights) result(output)
  implicit none
  
  type(FractionVector),  intent(in) :: wavevector
  type(CoupledStates),   intent(in) :: couplings(:)
  type(WavevectorState), intent(in) :: states(:)
  real(dp),              intent(in) :: weights(:)
  type(DensityMatrix)               :: output
  
  integer               :: no_states
  real(dp), allocatable :: coefficients(:,:)
  
  type(IntArray1d), allocatable :: selected_couplings(:)
  
  integer :: no_elements
  
  integer :: i,j,k,ialloc
  
  output%wavevector = wavevector
  
  if (size(states)==0) then
    return
  endif
  
  no_states = size(states(1)%coefficients)
  
  ! Convert state coefficients to a matrix, weighted by the state weights.
  ! This matrix is constructed as the transpose of the eigenvector matrix,
  !    to allow for fast access when calculating matrix elements.
  allocate(coefficients(size(states),no_states), stat=ialloc); call err(ialloc)
  do i=1,size(states)
    coefficients(i,:) = states(i)%coefficients * sqrt(weights(i))
  enddo
  
  output%state_ids = states(1)%state_ids
  selected_couplings = selected_states_couplings(couplings, output%state_ids)
  no_elements = sum([( size(selected_couplings(i)), &
                     & i=1,                         &
                     & size(selected_couplings)     )])
  allocate( output%bra_ids(no_elements), &
          & output%ket_ids(no_elements), &
          & output%values(no_elements),  &
          & stat=ialloc); call err(ialloc)
  k = 0
  do i=1,size(selected_couplings)
    do j=1,size(selected_couplings(i))
      k = k+1
      output%bra_ids(k) = output%state_ids(i)
      output%ket_ids(k) = output%state_ids(selected_couplings(i)%i(j))
      output%values(k) = dot_product(                  &
         & coefficients(:,selected_couplings(i)%i(j)), &
         & coefficients(:,i)                           )
    enddo
  enddo
end function
end module
