! ======================================================================
! Calculates potential energy or potential stress.
! ======================================================================
module braket_module
  use common_module
  
  use anharmonic_data_module
  use stress_prefactors_module
  use abstract_classes_module
  implicit none
  
  private
  
  public :: potential_energy
  
  public :: potential_stress
  
  interface potential_energy
    module procedure potential_energy_state
    module procedure potential_energy_state_state
  end interface
  
  interface potential_stress
    module procedure potential_stress_state
    module procedure potential_stress_state_state
  end interface
contains

! ----------------------------------------------------------------------
! The bra-ket of an arbitrary potential.
! Integrates across all dimensions and returns a real scalar.
! ----------------------------------------------------------------------
! Calculates <state|V|state> as a constant.
recursive function potential_energy_state(state,potential,subspace, &
   & subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceState),     intent(in) :: state
  class(PotentialData),     intent(in) :: potential
  type(DegenerateSubspace), intent(in) :: subspace
  class(SubspaceBasis),     intent(in) :: subspace_basis
  type(AnharmonicData),     intent(in) :: anharmonic_data
  real(dp)                             :: output
  
  type(PotentialPointer), allocatable :: integrated_potential
  
  integrated_potential = PotentialPointer(potential)
  call integrated_potential%braket( state,                            &
                                  & subspace        = subspace,       &
                                  & subspace_basis  = subspace_basis, &
                                  & anharmonic_data = anharmonic_data )
  output = integrated_potential%undisplaced_energy()
end function

! Calculates <bra|V|ket> as a constant.
recursive function potential_energy_state_state(bra,potential,ket,subspace, &
   & subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceState),     intent(in) :: bra
  class(PotentialData),     intent(in) :: potential
  class(SubspaceState),     intent(in) :: ket
  type(DegenerateSubspace), intent(in) :: subspace
  class(SubspaceBasis),     intent(in) :: subspace_basis
  type(AnharmonicData),     intent(in) :: anharmonic_data
  real(dp)                             :: output
  
  type(PotentialPointer), allocatable :: integrated_potential
  
  integrated_potential = PotentialPointer(potential)
  call integrated_potential%braket( bra,            &
                                  & ket,            &
                                  & subspace,       &
                                  & subspace_basis, &
                                  & anharmonic_data )
  output = integrated_potential%undisplaced_energy()
end function

! ----------------------------------------------------------------------
! The bra-ket of an arbitrary stress.
! Integrates across all dimensions and returns a constant tensor.
! ----------------------------------------------------------------------
! Calculates <state|stress|state> as a constant, for the potential stress.
recursive function potential_stress_state(state,stress,subspace, &
   & subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceState),     intent(in) :: state
  class(StressData),        intent(in) :: stress
  type(DegenerateSubspace), intent(in) :: subspace
  class(SubspaceBasis),     intent(in) :: subspace_basis
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(RealMatrix)                     :: output
  
  type(StressPointer), allocatable :: integrated_stress
  
  integrated_stress = StressPointer(stress)
  call integrated_stress%braket( state,                            &
                               & subspace        = subspace,       &
                               & subspace_basis  = subspace_basis, &
                               & anharmonic_data = anharmonic_data )
  output = integrated_stress%undisplaced_stress()
end function

! Calculates <bra|stress|ket> as a constant, for the potential stress.
recursive function potential_stress_state_state(bra,stress,ket,subspace, &
   & subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceState),     intent(in) :: bra
  class(StressData),        intent(in) :: stress
  class(SubspaceState),     intent(in) :: ket
  type(DegenerateSubspace), intent(in) :: subspace
  class(SubspaceBasis),     intent(in) :: subspace_basis
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(RealMatrix)                     :: output
  
  type(StressPointer), allocatable :: integrated_stress
  
  integrated_stress = StressPointer(stress)
  call integrated_stress%braket( bra,            &
                               & ket,            &
                               & subspace,       &
                               & subspace_basis, &
                               & anharmonic_data )
  output = integrated_stress%undisplaced_stress()
end function
end module
