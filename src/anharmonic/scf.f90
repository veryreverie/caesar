! ======================================================================
! The SCF cycle of VSCF.
! ======================================================================
module scf_module
  use constants_module, only : dp
  use string_module
  use io_module
  
  use linear_algebra_module
  use eigenstates_module
  
  private
  
  type, public :: ScfOutput
    type(SingleModeState), allocatable :: states(:,:)
    real(dp),              allocatable :: energies(:,:)
  end type
  
  public :: vscf
contains

function scf(potential,states) result(output)
  use potential_module
  use eigenstates_module
  implicit none
  
  type(PolynomialPotential), intent(in) :: potential
  type(SingleModeState),     intent(in) :: states(:,:)
  type(ScfOutput)                       :: output
  
  ! Array sizes.
  integer :: no_modes
  integer :: no_states
  integer :: no_coefficients
  
  ! Working variables.
  type(PolynomialPotential) :: mean_field
  type(RealMatrix)          :: hamiltonian
  type(RealEigenstuff)      :: eigenstuff
  
  ! Temporary variables.
  integer :: i,j,k,ialloc
  
  no_modes = size(states,2)
  no_states = size(states,1)
  no_coefficients = size(states(1,1)%coefficients)
  
  allocate( output%states(no_states,no_modes),   &
          & output%energies(no_states,no_modes), &
          & stat=ialloc); call err(ialloc)
  
  do i=1,no_modes
    ! Integrate over all other modes to get mean field potential.
    mean_field = potential
    do j=1,no_modes
      if (j/=i) then
        call mean_field%integrate_over_mode_average(j,states(:,j))
      endif
    enddo
    
    ! Construct Hamiltonian in terms of harmonic eigenfunctions.
    hamiltonian = mean_field%construct_hamiltonian(i,states(:,i))
    
    ! Diagonalise Hamiltonian to get new states in terms of old states.
    eigenstuff = calculate_eigenstuff(hamiltonian)
    
    ! Copy energies into output.
    output%energies(:,i) = eigenstuff%evals
    
    ! Calculate new states.
    do j=1,no_states
      output%states(j,i) = states(j,i) * 0.0_dp
      do k=1,no_states
        output%states(j,i) = output%states(j,i) &
                         & + eigenstuff%evecs(j,k) * states(k,i)
      enddo
    enddo
  enddo
end function

function vscf(potential,harmonic_states,max_scf_cycles, &
   & scf_convergence_threshold) result(output)
  use eigenstates_module
  use potential_module
  implicit none
  
  type(PolynomialPotential), intent(in) :: potential
  type(SingleModeState),     intent(in) :: harmonic_states(:,:)
  integer,                   intent(in) :: max_scf_cycles
  real(dp),                  intent(in) :: scf_convergence_threshold
  type(ScfOutput)                       :: output
  
  ! Array sizes.
  integer :: no_modes
  integer :: no_states
  
  ! Working variables.
  integer               :: scf_step
  real(dp), allocatable :: old_energies(:,:)
  
  ! Temporary variables.
  integer :: i,j,ialloc
  
  no_modes = size(harmonic_states,2)
  no_states = size(harmonic_states,1)
  
  ! Initialise output to harmonic values.
  output%states = harmonic_states
  allocate( output%energies(no_states,no_modes), &
          & old_energies(no_states,no_modes),    &
          & stat=ialloc); call err(ialloc)
  do i=1,no_modes
    do j=1,no_states
      output%energies(j,i) = (j-0.5_dp)*harmonic_states(j,i)%frequency
    enddo
  enddo
  
  ! SCF cycles.
  do scf_step=1,max_scf_cycles
    ! Record old energies.
    old_energies = output%energies
    
    ! Run SCF calculations.
    output = scf(potential,output%states)
    
    ! Check for convergence
    if (all(abs(output%energies-old_energies)<scf_convergence_threshold)) then
      exit
    elseif (scf_step==max_scf_cycles) then
      call err()
    endif
  enddo
end function
end module
