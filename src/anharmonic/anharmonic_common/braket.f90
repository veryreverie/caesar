! ======================================================================
! Calculates potential energy or potential stress.
! ======================================================================
module braket_module
  use common_module
  
  use subspace_state_module
  use basis_state_module
  use basis_states_module
  use anharmonic_data_module
  use stress_prefactors_module
  use abstract_classes_module
  use sparse_monomial_module
  implicit none
  
  private
  
  public :: integrate
  public :: potential_energy
  public :: potential_stress
  public :: harmonic_observables
  public :: effective_harmonic_observables
  
  interface integrate
    module procedure integrate_SubspaceState
    module procedure integrate_BasisState
    module procedure integrate_BasisStates
  end interface
  
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

! Integrate a ComplexMonomial.
impure elemental function integrate_SubspaceState(bra,monomial,ket, &
   & anharmonic_data) result(output)
  implicit none
  
  class(SubspaceState),  intent(in)           :: bra
  type(ComplexMonomial), intent(in)           :: monomial
  class(SubspaceState),  intent(in), optional :: ket
  type(AnharmonicData),  intent(in)           :: anharmonic_data
  type(ComplexMonomial)                       :: output
  
  type(ComplexUnivariate), allocatable :: unintegrated_modes(:)
  type(ComplexUnivariate), allocatable :: integrated_modes(:)
  
  unintegrated_modes = monomial%modes(exclude_ids = bra%mode_ids())
  
  output = ComplexMonomial( coefficient = monomial%coefficient, &
                          & modes       = unintegrated_modes    )
  
  integrated_modes = monomial%modes(      &
     & ids        = bra%mode_ids(),       &
     & paired_ids = bra%paired_mode_ids() )
  
  if (any(integrated_modes%total_power()/=0)) then
    output = output                                           &
         & * bra%integrate( SparseMonomial(integrated_modes), &
         &                  ket,                              &
         &                  anharmonic_data                   )
  endif
end function

impure elemental function integrate_BasisState(bra,monomial,ket,subspace, &
   & basis,anharmonic_data) result(output)
  implicit none
  
  class(BasisState),        intent(in)           :: bra
  type(ComplexMonomial),    intent(in)           :: monomial
  class(BasisState),        intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(ComplexMonomial)                          :: output
  
  type(ComplexUnivariate), allocatable :: unintegrated_modes(:)
  type(ComplexUnivariate), allocatable :: integrated_modes(:)
  
  unintegrated_modes = monomial%modes(                        &
     & exclude_ids = basis%mode_ids(subspace,anharmonic_data) )
  
  output = ComplexMonomial( coefficient = monomial%coefficient, &
                          & modes       = unintegrated_modes    )
  
  integrated_modes = monomial%modes(                                &
     & ids        = basis%mode_ids(subspace,anharmonic_data),       &
     & paired_ids = basis%paired_mode_ids(subspace,anharmonic_data) )
  
  if (any(integrated_modes%total_power()/=0)) then
    output = output                                             &
         & * basis%integrate( bra,                              &
         &                    SparseMonomial(integrated_modes), &
         &                    ket,                              &
         &                    subspace,                         &
         &                    anharmonic_data                   )
  endif
end function

impure elemental function integrate_BasisStates(states,thermal_energy, &
   & monomial,subspace,basis,anharmonic_data) result(output)
  implicit none
  
  class(BasisStates),       intent(inout) :: states
  real(dp),                 intent(in)    :: thermal_energy
  type(ComplexMonomial),    intent(in)    :: monomial
  type(DegenerateSubspace), intent(in)    :: subspace
  class(SubspaceBasis),     intent(in)    :: basis
  type(AnharmonicData),     intent(in)    :: anharmonic_data
  type(ComplexMonomial)                   :: output
  
  type(ComplexUnivariate), allocatable :: unintegrated_modes(:)
  type(ComplexUnivariate), allocatable :: integrated_modes(:)
  type(SparseMonomial)                 :: sparse
  
  integer     :: cache_location
  complex(dp) :: expectation
  
  unintegrated_modes = monomial%modes(                        &
     & exclude_ids = basis%mode_ids(subspace,anharmonic_data) )
  
  output = ComplexMonomial( coefficient = monomial%coefficient, &
                          & modes       = unintegrated_modes    )
  
  integrated_modes = monomial%modes(                                &
     & ids        = basis%mode_ids(subspace,anharmonic_data),       &
     & paired_ids = basis%paired_mode_ids(subspace,anharmonic_data) )
  
  if (any(integrated_modes%total_power()/=0)) then
    sparse = SparseMonomial(integrated_modes)
    cache_location = states%expectation_cache%cached_location(sparse)
    if (cache_location==0) then
      expectation = basis%integrate( states,         &
                                   & thermal_energy, &
                                   & sparse,         &
                                   & subspace,       &
                                   & anharmonic_data )
      call states%expectation_cache%cache(sparse, expectation)
    else
      expectation = states%expectation_cache%cached_expectation(cache_location)
    endif
    output = output * expectation
  endif
end function

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

! ----------------------------------------------------------------------
! Thermodynamic quantities derived from arbitrary potentials and stresses.
! ----------------------------------------------------------------------
! Calculate observables for harmonic basis, using a harmonic potential.
! N.B. the result is extensive, so will in general need to be normalised
!    to be per unit cell or similar.
impure elemental function harmonic_observables(thermal_energy,stress, &
   & stress_prefactor,frequency,num_dimensions,supercell_size,        &
   & anharmonic_data) result(output) 
  implicit none
  
  real(dp),             intent(in)           :: thermal_energy
  class(StressData),    intent(in), optional :: stress
  type(RealMatrix),     intent(in), optional :: stress_prefactor
  real(dp),             intent(in)           :: frequency
  integer,              intent(in)           :: num_dimensions
  integer,              intent(in)           :: supercell_size
  type(AnharmonicData), intent(in)           :: anharmonic_data
  type(ThermodynamicData)                    :: output
  
  type(RealMatrix), allocatable :: potential_stress
  real(dp),         allocatable :: volume
  
  if (present(stress).neqv.present(stress_prefactor)) then
    call print_line(CODE_ERROR//': Either both or neither of stress and &
       &stress_prefactor must be given.')
    call err()
  endif
  
  if (present(stress)) then
    potential_stress = stress%harmonic_expectation( frequency,        &
                   &                                thermal_energy,   &
                   &                                supercell_size,   &
                   &                                anharmonic_data ) &
                   & / num_dimensions
    volume = anharmonic_data%structure%volume
  endif
  
  output = ThermodynamicData( thermal_energy,     &
       &                      frequency,          &
       &                      stress_prefactor,   &
       &                      potential_stress,   &
       &                      volume            ) &
       & * num_dimensions
end function

! Calculate observables for harmonic basis, using full potential.
! N.B. the result is extensive, so will in general need to be normalised
!    to be per unit cell or similar.
impure elemental function effective_harmonic_observables(thermal_energy, &
   & potential,stress,stress_prefactor,frequency,num_dimensions,         &
   & supercell_size,anharmonic_data) result(output) 
  implicit none
  
  real(dp),             intent(in)           :: thermal_energy
  class(PotentialData), intent(in)           :: potential
  class(StressData),    intent(in), optional :: stress
  type(RealMatrix),     intent(in), optional :: stress_prefactor
  real(dp),             intent(in)           :: frequency
  integer,              intent(in)           :: num_dimensions
  integer,              intent(in)           :: supercell_size
  type(AnharmonicData), intent(in)           :: anharmonic_data
  type(ThermodynamicData)                    :: output
  
  real(dp) :: harmonic_expectation
  real(dp) :: potential_expectation
  
  real(dp)         :: volume
  type(RealMatrix) :: potential_stress
  type(RealMatrix) :: kinetic_stress
  real(dp)         :: pressure
  
  ! Calculate <T+V> for the effective harmonic potential.
  output = harmonic_observables( thermal_energy,   &
                               & stress,           &
                               & stress_prefactor, &
                               & frequency,        &
                               & num_dimensions,   &
                               & supercell_size,   &
                               & anharmonic_data   )
  
  ! Subtract <V> for the effective harmonic potential,
  !    and add <V> for the input potential.
  ! N.B. <V> for the effective harmonic potential is just <T+V>/2.
  ! N.B. the kinetic energy, entropy and stress do not change.
  harmonic_expectation = output%energy/2
  potential_expectation = potential%harmonic_expectation( &
                                        & frequency,      &
                                        & thermal_energy, &
                                        & supercell_size, &
                                        & anharmonic_data )
  
  output%energy = output%energy        &
              & - harmonic_expectation &
              & + potential_expectation
  output%free_energy = output%free_energy   &
                   & - harmonic_expectation &
                   & + potential_expectation
  if (allocated(output%enthalpy)) then
    output%enthalpy = output%enthalpy &
                  & - harmonic_expectation &
                  & + potential_expectation
    output%gibbs = output%gibbs &
               & - harmonic_expectation &
               & + potential_expectation
  endif
end function
end module
