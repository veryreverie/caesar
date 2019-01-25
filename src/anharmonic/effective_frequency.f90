! ======================================================================
! Uses the self-consistent harmonic approximation to find the effective
!    harmonic potential which most closely matches the given single-subspace
!    potential at the given temperature.
! ======================================================================
module effective_frequency_module
  use common_module
  
  use states_module
  use anharmonic_common_module
  use potentials_module
  implicit none
  
  private
  
  public :: calculate_effective_frequency
contains

! An effective harmonic potential Vh(w) with frequency w is defined as
!    V(w) = sum_i 0.5 w^2 (u_i)^2, where u_i is the ith degenerate mode.
! The effective frequency of the VSCF potential Vv at thermal energy KbT is
!    given by the frequency of the effective potential Vh(w) whose states
!    minimise the free energy of the VSCF Hamiltonian.
function calculate_effective_frequency(potential,subspace,anharmonic_data,   &
   & thermal_energy,initial_frequency,no_basis_states,frequency_convergence) &
   & result(output)
  implicit none
  
  class(PotentialData),     intent(in) :: potential
  type(DegenerateSubspace), intent(in) :: subspace
  type(AnharmonicData),     intent(in) :: anharmonic_data
  real(dp),                 intent(in) :: thermal_energy
  real(dp),                 intent(in) :: initial_frequency
  integer,                  intent(in) :: no_basis_states
  real(dp),                 intent(in) :: frequency_convergence
  real(dp)                             :: output
  
  real(dp) :: frequency
  real(dp) :: frequencies(3)
  real(dp) :: free_energies(3)
  
  type(NewtonRaphson) :: solver
  
  integer :: i
  
  frequency = initial_frequency
  
  solver = NewtonRaphson(                                   &
     & starting_value      = initial_frequency,             &
     & finite_displacement = 0.01_dp*frequency_convergence, &
     & convergence_threshold = frequency_convergence,       &
     & lower_bound           = 0.0_dp                       )
  do
    frequencies = solver%get_inputs()
    
    do i=1,3
      free_energies(i) = calculate_free_energy( potential,       &
                                              & frequencies(i),  &
                                              & thermal_energy,  &
                                              & no_basis_states, &
                                              & subspace,        &
                                              & anharmonic_data  )
    enddo
    
    call solver%set_outputs(free_energies)
    
    if (solver%converged()) then
      output = solver%solution()
      exit
    endif
  enddo
end function

! ----------------------------------------------------------------------
! Calculate the free energy of the given potential in the given effective
!    harmonic basis.
! ----------------------------------------------------------------------
! If Fa(b) is the free energy of the Hamiltonian a in the basis of Hamiltonian
!    b, then Fh(h) is the free energy of the effective harmonic system,
!    and Fv(v) is the free energy of the VSCF system.
! The Gibbs-Bogoliubov inequality states:
!    Fv(v) <= Fv(h)
! If Pa_i is the bose factor of the ith state in the basis of Hamiltonian a,
!
!    Fh(h) = sum_i(Ph_i <ih|T+Vh|ih>) + KbT sum_i(Ph_i ln(Ph_i))
!
!    Fv(h) = sum_i(Ph_i <ih|T+Vv|ih>) + KbT sum_i(Ph_i ln(Ph_i))
!
! -> Fv(h) = Fh(h) + sum_i(Ph_i <ih|Vv-Vh|ih>)
function calculate_free_energy(potential,frequency,thermal_energy,no_states, &
   & subspace,anharmonic_data) result(output)
  implicit none
  
  class(PotentialData),     intent(in) :: potential
  real(dp),                 intent(in) :: frequency
  real(dp),                 intent(in) :: thermal_energy
  integer,                  intent(in) :: no_states
  type(DegenerateSubspace), intent(in) :: subspace
  type(AnharmonicData),     intent(in) :: anharmonic_data
  real(dp)                             :: output
  
  type(StructureData)     :: supercell
  type(ThermodynamicData) :: harmonic_thermodynamics
  real(dp)                :: harmonic_free_energy
  real(dp)                :: harmonic_potential_energy_expectation
  real(dp)                :: anharmonic_potential_energy_expectation
  
  supercell = anharmonic_data%anharmonic_supercell
  
  ! Calculate the free energy per primitive cell
  !    of the harmonic system in the harmonic basis.
  harmonic_thermodynamics = ThermodynamicData(thermal_energy, frequency)
  harmonic_free_energy = harmonic_thermodynamics%free_energy &
                     & * size(subspace)                      &
                     & / real(supercell%sc_size,dp)
  
  ! Calculate <V> for the harmonic and anharmonic potentials.
  
  ! U = <V> + <T>. Under the harmonic approximation, <V>=<T>, so <V>=U/2.
  harmonic_potential_energy_expectation = harmonic_thermodynamics%energy &
                                      & * size(subspace)                 &
                                      & / (2.0_dp * real(supercell%sc_size,dp))
  anharmonic_potential_energy_expectation =            &
     & potential%harmonic_expectation( frequency,      &
     &                                 thermal_energy, &
     &                                 no_states,      &
     &                                 subspace,       &
     &                                 anharmonic_data )
  
  output = harmonic_free_energy                  &
       & - harmonic_potential_energy_expectation &
       & + anharmonic_potential_energy_expectation
end function
end module
