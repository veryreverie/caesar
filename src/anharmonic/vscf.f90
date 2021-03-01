! ======================================================================
! Run VSCF on a potential to generate VSCF ground states.
! ======================================================================
! Each step in the routine:
!    1) starts from a set of single-subspace potentials, P.
!    2) uses P to generate a set of single-subspace states, S(P),
!         such that the states diagonalise the potentials.
!    3) uses S(P) and P to generate the free energy, F(P,S(P)).
!    4) checks for convergence of F(P,S(P)), and returns if converged.
!    5) uses S(P) to generate a new set of single-subspace potentials, P(S(P)).
!    6) uses a Pulay scheme to generate the next P, which attempts to minimise
!          F(P,S(P)) subject to P(S(P))=P.
!
! The Pulay scheme sees P as a vector of real coefficients. It is down to the
!    implementation of the potential to define how P is generated.
module caesar_vscf_module
  use caesar_common_module
  
  use caesar_states_module
  
  use caesar_anharmonic_common_module
  use caesar_potentials_module
  
  use caesar_generate_subspace_potentials_module
  implicit none
  
  private
  
  public :: VscfOutput
  
  public :: run_vscf
  
  type, extends(NoDefaultConstructor) :: VscfOutput
    type(PotentialPointer)           :: potential
    type(BasisStatesPointer)         :: states
    type(StressPointer), allocatable :: stress
  end type
  
  interface VscfOutput
    ! Constructor.
    impure elemental module function new_VscfOutput(potential,states,stress) &
       & result(this) 
      class(PotentialData), intent(in)           :: potential
      Class(BasisStates),   intent(in)           :: states
      class(StressData),    intent(in), optional :: stress
      type(VscfOutput)                           :: this
    end function
  end interface
  
  interface
    ! Vscf routine.
    module function run_vscf(potential,stress,subspaces,subspace_bases,   &
       & thermal_energy,state_energy_cutoff,frequencies,convergence_data, &
       & anharmonic_data,random_generator,starting_configuration,         &
       & convergence_file) result(output) 
      class(PotentialData),     intent(in)              :: potential
      class(StressData),        intent(in),    optional :: stress
      type(DegenerateSubspace), intent(in)              :: subspaces(:)
      class(SubspaceBasis),     intent(in)              :: subspace_bases(:)
      real(dp),                 intent(in)              :: thermal_energy
      real(dp),                 intent(in)              :: state_energy_cutoff
      real(dp),                 intent(in)              :: frequencies(:)
      type(ConvergenceData),    intent(in)              :: convergence_data
      type(AnharmonicData),     intent(in)              :: anharmonic_data
      type(RandomReal),         intent(in)              :: random_generator
      type(VscfOutput),         intent(in),    optional :: starting_configuration(:)
      type(OFile),              intent(inout), optional :: convergence_file
      type(VscfOutput), allocatable                     :: output(:)
    end function
  end interface
end module
