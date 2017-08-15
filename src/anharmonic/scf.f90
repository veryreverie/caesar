! ======================================================================
! The SCF cycle of VSCF.
! ======================================================================
module scf_module
  use constants_module, only : dp
  use string_module
  use io_module
  use linear_algebra_module
contains

function scf(potential,eigenstuff) result(output)
  use potential_module
  implicit none
  
  type(PolynomialPotential), intent(in)  :: potential
  type(RealEigenstuff),      intent(in)  :: eigenstuff(:)
  type(RealEigenstuff), allocatable      :: output(:)
  
  ! Array sizes.
  integer :: no_modes
  integer :: no_states
  
  ! Working variables.
  type(PolynomialPotential) :: mean_field
  type(RealMatrix)          :: hamiltonian
  
  ! Temporary variables.
  integer :: i,j,ialloc
  
  no_modes = size(eigenstuff)
  no_states = size(eigenstuff(1)%evals)
  
  allocate(output(no_modes), stat=ialloc); call err(ialloc)
  
  do i=1,no_modes
    ! Integrate over all other modes to get mean field potential.
    mean_field = potential
    do j=1,no_modes
      if (j/=i) then
        call mean_field%integrate_over_mode_average(i,eigenstuff(i))
      endif
    enddo
    
    ! Construct Hamiltonian in terms of harmonic eigenfunctions.
    hamiltonian = mean_field%construct_hamiltonian(i)
    
    ! Diagonalise Hamiltonian to get new eigenstates.
    output(i) = calculate_eigenstuff(hamiltonian)
  enddo
end function
end module
