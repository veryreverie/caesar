! ======================================================================
! The SCF cycle of VSCF.
! ======================================================================
module scf_module
  use constants_module, only : dp
  use string_module
  use io_module
  
  private
  
  public :: vscf
contains

! ----------------------------------------------------------------------
! Updates the VSCF states along each mode to be consistent with those
!    along all other modes.
! ----------------------------------------------------------------------
function scf(potential,input,harmonic_states) result(output)
  use linear_algebra_module
  use potential_module
  use harmonic_states_module
  use vscf_states_module
  implicit none
  
  type(PolynomialPotential), intent(in) :: potential
  type(VscfStates),          intent(in) :: input(:)
  type(HarmonicStates),      intent(in) :: harmonic_states(:)
  type(VscfStates), allocatable         :: output(:)
  
  ! Working variables.
  type(PolynomialPotential) :: mean_field
  type(RealMatrix)          :: hamiltonian
  
  ! Temporary variables.
  integer :: i,j
  
  output = input
  
  do i=1,size(input)
    
    call print_line('')
    call print_line('V:')
    call print_line(potential)
    
    ! Integrate over all other modes to get mean field potential.
    mean_field = potential
    do j=1,size(input)
      if (j/=i) then
        call mean_field%integrate_over_mode_average(input(j))
        
        call print_line('')
        call print_line('<'//j//'|V|'//j//'>:')
        call print_line(mean_field)
      endif
    enddo
    
    ! Construct Hamiltonian in terms of harmonic eigenfunctions.
    hamiltonian = mean_field%construct_hamiltonian(harmonic_states(i))
    
    call print_line(hamiltonian)
    
    ! Diagonalise Hamiltonian to get new states.
    call output(i)%update_hamiltonian(hamiltonian)
  enddo
end function

! ----------------------------------------------------------------------
! Constructs VSCF single-mode states from the potential and harmonic states.
! ----------------------------------------------------------------------
function vscf(potential,harmonic_states,max_scf_cycles, &
   & scf_convergence_threshold) result(output)
  use potential_module
  use harmonic_states_module
  use vscf_states_module
  implicit none
  
  type(PolynomialPotential), intent(in) :: potential
  type(HarmonicStates),      intent(in) :: harmonic_states(:)
  integer,                   intent(in) :: max_scf_cycles
  real(dp),                  intent(in) :: scf_convergence_threshold
  type(VscfStates), allocatable         :: output(:)
  
  ! Array sizes.
  integer :: no_modes
  
  ! Working variables.
  integer                       :: scf_step
  type(VscfStates), allocatable :: old_states(:)
  
  real(dp) :: max_energy_change
  
  ! Temporary variables.
  integer :: i,j,ialloc
  
  no_modes = size(harmonic_states)
  
  ! Initialise output to harmonic values.
  allocate(output(no_modes), stat=ialloc); call err(ialloc)
  do i=1,no_modes
    output(i) = construct_harmonic_vscf_states(harmonic_states(i))
  enddo
  
  ! SCF cycles.
  do scf_step=1,max_scf_cycles
    ! Record old states.
    old_states = output
    
    ! Run SCF calculations.
    output = scf(potential,output,harmonic_states)
    
    ! Check for convergence
    max_energy_change = 0
    do i=1,no_modes
      do j=0,output(i)%cutoff()
        max_energy_change = max( max_energy_change,           &
                               & abs( old_states(i)%energy(j) &
                               &    - output(i)%energy(j)))
      enddo
    enddo
    call print_line(max_energy_change)
    if (max_energy_change<scf_convergence_threshold) then
      exit
    elseif (scf_step==max_scf_cycles) then
      call err()
    endif
  enddo
end function
end module
