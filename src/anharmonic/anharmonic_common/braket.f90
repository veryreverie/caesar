! ======================================================================
! Calculates potential energy or potential stress.
! ======================================================================
module braket_module
  use common_module
  
  use basis_state_module
  use subspace_state_module
  use anharmonic_data_module
  use stress_prefactors_module
  use abstract_classes_module
  implicit none
  
  private
  
  public :: potential_energy
  public :: potential_stress
  public :: effective_harmonic_observables
  
  interface potential_energy
    module procedure potential_energy_BasisState
    module procedure potential_energy_BasisState_BasisState
    module procedure potential_energy_SubspaceState
    module procedure potential_energy_SubspaceState_SubspaceState
  end interface
  
  interface potential_stress
    module procedure potential_stress_BasisState
    module procedure potential_stress_BasisState_BasisState
    module procedure potential_stress_SubspaceState
    module procedure potential_stress_SubspaceState_SubspaceState
  end interface
contains

! ----------------------------------------------------------------------
! The bra-ket of an arbitrary potential.
! Integrates across all dimensions and returns a real scalar.
! ----------------------------------------------------------------------
! Calculates <state|V|state> as a constant.
recursive function potential_energy_BasisState(state,potential,subspace, &
   & subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(BasisState),        intent(in) :: state
  class(PotentialData),     intent(in) :: potential
  type(DegenerateSubspace), intent(in) :: subspace
  class(SubspaceBasis),     intent(in) :: subspace_basis
  type(AnharmonicData),     intent(in) :: anharmonic_data
  real(dp)                             :: output
  
  type(PotentialPointer), allocatable :: integrated_potential
  
  integrated_potential = PotentialPointer(potential)
  call integrated_potential%braket( bra             = state,          &
                                  & subspace        = subspace,       &
                                  & subspace_basis  = subspace_basis, &
                                  & anharmonic_data = anharmonic_data )
  output = integrated_potential%undisplaced_energy()
end function

! Calculates <bra|V|ket> as a constant.
recursive function potential_energy_BasisState_BasisState(bra,potential,ket, &
   & subspace,subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(BasisState),        intent(in) :: bra
  class(PotentialData),     intent(in) :: potential
  class(BasisState),        intent(in) :: ket
  type(DegenerateSubspace), intent(in) :: subspace
  class(SubspaceBasis),     intent(in) :: subspace_basis
  type(AnharmonicData),     intent(in) :: anharmonic_data
  real(dp)                             :: output
  
  type(PotentialPointer), allocatable :: integrated_potential
  
  integrated_potential = PotentialPointer(potential)
  call integrated_potential%braket( bra             = bra,            &
                                  & ket             = ket,            &
                                  & subspace        = subspace,       &
                                  & subspace_basis  = subspace_basis, &
                                  & anharmonic_data = anharmonic_data )
  output = integrated_potential%undisplaced_energy()
end function

! ----------------------------------------------------------------------
! The bra-ket of an arbitrary stress.
! Integrates across all dimensions and returns a constant tensor.
! ----------------------------------------------------------------------
! Calculates <state|stress|state> as a constant, for the potential stress.
recursive function potential_stress_BasisState(state,stress,subspace, &
   & subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(BasisState),        intent(in) :: state
  class(StressData),        intent(in) :: stress
  type(DegenerateSubspace), intent(in) :: subspace
  class(SubspaceBasis),     intent(in) :: subspace_basis
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(RealMatrix)                     :: output
  
  type(StressPointer), allocatable :: integrated_stress
  
  integrated_stress = StressPointer(stress)
  call integrated_stress%braket( bra             = state,          &
                               & subspace        = subspace,       &
                               & subspace_basis  = subspace_basis, &
                               & anharmonic_data = anharmonic_data )
  output = integrated_stress%undisplaced_stress()
end function

! Calculates <bra|stress|ket> as a constant, for the potential stress.
recursive function potential_stress_BasisState_BasisState(bra,stress,ket, &
   & subspace,subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(BasisState),        intent(in) :: bra
  class(StressData),        intent(in) :: stress
  class(BasisState),        intent(in) :: ket
  type(DegenerateSubspace), intent(in) :: subspace
  class(SubspaceBasis),     intent(in) :: subspace_basis
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(RealMatrix)                     :: output
  
  type(StressPointer), allocatable :: integrated_stress
  
  integrated_stress = StressPointer(stress)
  call integrated_stress%braket( bra             = bra,            &
                               & ket             = ket,            &
                               & subspace        = subspace,       &
                               & subspace_basis  = subspace_basis, &
                               & anharmonic_data = anharmonic_data )
  output = integrated_stress%undisplaced_stress()
end function

! ----------------------------------------------------------------------
! The bra-ket of an arbitrary potential.
! Integrates across all dimensions and returns a real scalar.
! ----------------------------------------------------------------------
! Calculates <state|V|state> as a constant.
recursive function potential_energy_SubspaceState(state,potential, &
   & anharmonic_data) result(output)
  implicit none
  
  class(SubspaceState), intent(in) :: state
  class(PotentialData), intent(in) :: potential
  type(AnharmonicData), intent(in) :: anharmonic_data
  real(dp)                         :: output
  
  type(PotentialPointer), allocatable :: integrated_potential
  
  integrated_potential = PotentialPointer(potential)
  call integrated_potential%braket( bra             = state,          &
                                  & anharmonic_data = anharmonic_data )
  output = integrated_potential%undisplaced_energy()
end function

! Calculates <bra|V|ket> as a constant.
recursive function potential_energy_SubspaceState_SubspaceState(bra, &
   & potential,ket,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceState), intent(in) :: bra
  class(PotentialData), intent(in) :: potential
  class(SubspaceState), intent(in) :: ket
  type(AnharmonicData), intent(in) :: anharmonic_data
  real(dp)                         :: output
  
  type(PotentialPointer), allocatable :: integrated_potential
  
  integrated_potential = PotentialPointer(potential)
  call integrated_potential%braket( bra             = bra,            &
                                  & ket             = ket,            &
                                  & anharmonic_data = anharmonic_data )
  output = integrated_potential%undisplaced_energy()
end function

! ----------------------------------------------------------------------
! The bra-ket of an arbitrary stress.
! Integrates across all dimensions and returns a constant tensor.
! ----------------------------------------------------------------------
! Calculates <state|stress|state> as a constant, for the potential stress.
recursive function potential_stress_SubspaceState(state,stress, &
   & anharmonic_data) result(output)
  implicit none
  
  class(SubspaceState), intent(in) :: state
  class(StressData),    intent(in) :: stress
  type(AnharmonicData), intent(in) :: anharmonic_data
  type(RealMatrix)                 :: output
  
  type(StressPointer), allocatable :: integrated_stress
  
  integrated_stress = StressPointer(stress)
  call integrated_stress%braket( bra             = state,          &
                               & anharmonic_data = anharmonic_data )
  output = integrated_stress%undisplaced_stress()
end function

! Calculates <bra|stress|ket> as a constant, for the potential stress.
recursive function potential_stress_SubspaceState_SubspaceState(bra,stress, &
   & ket,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceState), intent(in) :: bra
  class(StressData),    intent(in) :: stress
  class(SubspaceState), intent(in) :: ket
  type(AnharmonicData), intent(in) :: anharmonic_data
  type(RealMatrix)                 :: output
  
  type(StressPointer), allocatable :: integrated_stress
  
  integrated_stress = StressPointer(stress)
  call integrated_stress%braket( bra             = bra,            &
                               & ket             = ket,            &
                               & anharmonic_data = anharmonic_data )
  output = integrated_stress%undisplaced_stress()
end function

! Calculate <T+V> for effective harmonic states with
!    effective harmonic weightings.
! N.B. the result is extensive, so will in general need to be normalised
!    to be per unit cell or similar.
impure elemental function effective_harmonic_observables(thermal_energy, &
   & potential,frequency,num_dimensions,anharmonic_data) result(output)
  implicit none
  
  real(dp),             intent(in) :: thermal_energy
  class(PotentialData), intent(in) :: potential
  real(dp),             intent(in) :: frequency
  integer,              intent(in) :: num_dimensions
  type(AnharmonicData), intent(in) :: anharmonic_data
  type(ThermodynamicData)          :: output
  
  real(dp) :: harmonic_expectation
  real(dp) :: potential_expectation
  
  ! Calculate <T+V> for the effective harmonic potential.
  output = ThermodynamicData(thermal_energy, frequency) * num_dimensions
  
  ! Subtract <V> for the effective harmonic potential,
  !    and add <V> for the input potential.
  ! N.B. <V> for the effective harmonic potential is just <T+V>/2.
  ! N.B. the kinetic energy and entropy do not change.
  harmonic_expectation = output%energy/2
  potential_expectation = potential%harmonic_expectation(   &
                      &                   frequency,        &
                      &                   thermal_energy,   &
                      &                   anharmonic_data ) &
                      & * anharmonic_data%anharmonic_supercell%sc_size
  
  output%energy = output%energy        &
              & - harmonic_expectation &
              & + potential_expectation
  output%free_energy = output%free_energy   &
                   & - harmonic_expectation &
                   & + potential_expectation
end function
end module
