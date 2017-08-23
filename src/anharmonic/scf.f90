! ======================================================================
! The SCF cycle of VSCF.
! ======================================================================
module scf_module
  use constants_module, only : dp
  use string_module
  use io_module
  use linear_algebra_module
contains

function scf(potential,harmonic_states,eigenstuff) result(output)
  use potential_module
  use eigenstates_module
  implicit none
  
  type(PolynomialPotential), intent(in) :: potential
  type(SingleModeState),     intent(in) :: harmonic_states(:,:)
  type(RealEigenstuff),      intent(in) :: eigenstuff(:)
  type(RealEigenstuff), allocatable     :: output(:)
  
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
        call mean_field%integrate_over_mode_average( j,                    &
                                                   & harmonic_states(:,j), &
                                                   & eigenstuff(j))
      endif
    enddo
    
    ! Construct Hamiltonian in terms of harmonic eigenfunctions.
    hamiltonian = mean_field%construct_hamiltonian(i, harmonic_states(:,i))
    
    ! Diagonalise Hamiltonian to get new eigenstates.
    output(i) = calculate_eigenstuff(hamiltonian)
  enddo
end function

function vscf(harmonic_states,potential,max_scf_cycles, &
   & scf_convergence_threshold) result(output)
  use eigenstates_module
  use potential_module
  implicit none
  
  type(SingleModeState),     intent(in) :: harmonic_states(:,:)
  type(PolynomialPotential), intent(in) :: potential
  integer,                   intent(in) :: max_scf_cycles
  real(dp),                  intent(in) :: scf_convergence_threshold
  type(SingleModeState), allocatable    :: output(:,:)
  
  integer :: no_harmonic_states
  integer :: no_modes
  
  integer                           :: scf_step
  type(RealEigenstuff), allocatable :: vscf_eigenstuff(:)
  real(dp),             allocatable :: energy_change(:,:)
  
  integer :: i,j,k,ialloc
  
  no_harmonic_states = size(harmonic_states,1)
  no_modes = size(harmonic_states,2)
  
  ! Initialise eigenstuff to harmonic values.
  allocate( vscf_eigenstuff(no_modes), &
          & stat=ialloc); call err(ialloc)
  do i=1,no_modes
    vscf_eigenstuff(i)%evecs = real(int(identity(no_harmonic_states)))
    do j=1,no_harmonic_states
      vscf_eigenstuff(i)%evals(j) = (j-0.5_dp)*harmonic_states(j,i)%frequency
    enddo
  enddo
  
  ! Initialise energy change array.
  allocate( energy_change(no_harmonic_states, no_modes), &
          & stat=ialloc); call err(ialloc)
  
  ! SCF cycles.
  do scf_step=1,max_scf_cycles
    ! Record old eigenvalues.
    do i=1,no_modes
      energy_change(:,i) = vscf_eigenstuff(j)%evals
    enddo
    
    ! Run SCF calculations.
    vscf_eigenstuff = scf(potential,harmonic_states,vscf_eigenstuff)
    
    ! Calculate L2-norm change in eigenvalues.
    do i=1,no_modes
      energy_change(:,i) = l2_norm(vec( vscf_eigenstuff(j)%evals &
                                      & - energy_change(:,i) ))
    enddo
    
    ! Check for convergence
    if (sum(energy_change)<scf_convergence_threshold) then
      exit
    elseif (scf_step==max_scf_cycles) then
      call err()
    endif
  enddo
  
  ! Construct VSCF basis.
  allocate( output(no_harmonic_states,no_modes), &
          & stat=ialloc); call err(ialloc)
  do i=1,no_modes
    do j=1,no_harmonic_states
      output(j,i) = harmonic_states(k,j)
      output(j,i)%coefficients = 0
      do k=1,size(vscf_eigenstuff(1))
        output(j,i)%coefficients = output(k,j)%coefficients &
                               & + vscf_eigenstuff(i)%evecs(j,k) &
                               & * harmonic_states(k,i)%coefficients
      enddo
    enddo
  enddo
end function
end module
