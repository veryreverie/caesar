! ======================================================================
! The SCF cycle of VSCF.
! ======================================================================
module scf_module
  use constants_module, only : dp
  use string_module
  use io_module
  use linear_algebra_module
  
  type :: ScfResult
    real(dp)                      :: energy_change
    type(RealMatrix), allocatable :: eigenstates(:)
  end type

contains

function scf(harmonic_couplings,potential,eigenstates) result(output)
  use potential_basis_module
  implicit none
  
  type(RealMatrix), intent(in) :: harmonic_couplings(:,:)
  type(Monomial),   intent(in) :: potential(:)
  type(RealMatrix), intent(in) :: eigenstates(:)
  type(ScfResult)              :: output
  
  ! Array sizes.
  integer :: no_modes
  integer :: no_states
  
  ! Working variables.
  type(Monomial), allocatable :: mean_field(:)
  type(RealMatrix)            :: hamiltonian
  type(RealEigenstuff)        :: eigenstuff
  
  ! Temporary variables.
  real(dp), allocatable :: zero(:,:)
  integer :: i,j,ialloc
  
  no_modes = size(eigenstates)
  no_states = size(eigenstates(1),1)
  
  allocate(zero(no_states,no_states), stat=ialloc); call err(ialloc)
  zero = 0
  
  do i=1,no_modes
    ! Integrate over all other modes to get mean field potential.
    mean_field = potential
    do j=1,no_modes
      if (j/=i) then
        mean_field = integrate_over_mode_average( mean_field,              &
                                                & i,                       &
                                                & harmonic_couplings(:,i), &
                                                & eigenstates(i))
      endif
    enddo
    
    ! Construct Hamiltonian in terms of harmonic eigenfunctions.
    hamiltonian = zero
    do j=1,size(mean_field)
      hamiltonian = hamiltonian               &
                & + mean_field(j)%coefficient &
                & * harmonic_couplings(mean_field(j)%powers(i),i)
    enddo
    
    ! Diagonalise Hamiltonian to get new eigenstates.
    eigenstuff = calculate_eigenstuff(hamiltonian)
  enddo
end function

end module
