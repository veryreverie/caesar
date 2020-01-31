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
  public :: integrate_to_constant
  public :: harmonic_observables
  public :: effective_harmonic_observables
  
  interface integrate
    module procedure integrate_SubspaceState
    module procedure integrate_BasisState
    module procedure integrate_BasisStates
  end interface
  
  interface integrate_to_constant
    module procedure integrate_to_constant_SubspaceState_ComplexMonomial
    module procedure integrate_to_constant_BasisState_ComplexMonomial
    module procedure integrate_to_constant_BasisStates_ComplexMonomial
    
    module procedure integrate_to_constant_SubspaceState_ComplexPolynomial
    module procedure integrate_to_constant_BasisState_ComplexPolynomial
    module procedure integrate_to_constant_BasisStates_ComplexPolynomial
  end interface
contains

! Integrate a the parts of a ComplexMonomial in a given subspace.
impure elemental subroutine integrate_SubspaceState(monomial,bra,ket, &
   & anharmonic_data)
  implicit none
  
  type(ComplexMonomial), intent(inout)        :: monomial
  class(SubspaceState),  intent(in)           :: bra
  class(SubspaceState),  intent(in), optional :: ket
  type(AnharmonicData),  intent(in)           :: anharmonic_data
  
  type(ComplexUnivariate), allocatable :: unintegrated_modes(:)
  type(ComplexUnivariate), allocatable :: integrated_modes(:)
  
  unintegrated_modes = monomial%modes(exclude_ids = bra%mode_ids())
  
  integrated_modes = monomial%modes(      &
     & ids        = bra%mode_ids(),       &
     & paired_ids = bra%paired_mode_ids() )
  
  call monomial%set_modes(unintegrated_modes)
  
  if (any(integrated_modes%total_power()/=0)) then
    monomial%coefficient = monomial%coefficient                             &
                       & * bra%integrate( SparseMonomial(integrated_modes), &
                       &                  ket,                              &
                       &                  anharmonic_data                   )
  endif
end subroutine

impure elemental subroutine integrate_BasisState(monomial,bra,ket,subspace, &
   & basis,anharmonic_data)
  implicit none
  
  type(ComplexMonomial),    intent(inout)        :: monomial
  class(BasisState),        intent(in)           :: bra
  class(BasisState),        intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  
  type(ComplexUnivariate), allocatable :: unintegrated_modes(:)
  type(ComplexUnivariate), allocatable :: integrated_modes(:)
  
  unintegrated_modes = monomial%modes(                        &
     & exclude_ids = basis%mode_ids(subspace,anharmonic_data) )
  
  integrated_modes = monomial%modes(                                &
     & ids        = basis%mode_ids(subspace,anharmonic_data),       &
     & paired_ids = basis%paired_mode_ids(subspace,anharmonic_data) )
  
  call monomial%set_modes(unintegrated_modes)
  
  if (any(integrated_modes%total_power()/=0)) then
    monomial%coefficient = monomial%coefficient                               &
                       & * basis%integrate( bra,                              &
                       &                    SparseMonomial(integrated_modes), &
                       &                    ket,                              &
                       &                    subspace,                         &
                       &                    anharmonic_data                   )
  endif
end subroutine

impure elemental subroutine integrate_BasisStates(monomial,states, &
   & thermal_energy,subspace,basis,anharmonic_data) 
  implicit none
  
  type(ComplexMonomial),    intent(inout) :: monomial
  class(BasisStates),       intent(inout) :: states
  real(dp),                 intent(in)    :: thermal_energy
  type(DegenerateSubspace), intent(in)    :: subspace
  class(SubspaceBasis),     intent(in)    :: basis
  type(AnharmonicData),     intent(in)    :: anharmonic_data
  
  type(ComplexUnivariate), allocatable :: unintegrated_modes(:)
  type(ComplexUnivariate), allocatable :: integrated_modes(:)
  type(SparseMonomial)                 :: sparse
  
  integer     :: cache_location
  complex(dp) :: expectation
  
  unintegrated_modes = monomial%modes(                        &
     & exclude_ids = basis%mode_ids(subspace,anharmonic_data) )
  
  integrated_modes = monomial%modes(                                &
     & ids        = basis%mode_ids(subspace,anharmonic_data),       &
     & paired_ids = basis%paired_mode_ids(subspace,anharmonic_data) )
  
  call monomial%set_modes(unintegrated_modes)
  
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
    monomial%coefficient = monomial%coefficient * expectation
  endif
end subroutine

! As integrate, except that the whole ComplexMonomial must be in the given
!    subspace.
! Since the resulting ComplexMonomial would just be a constant, just the
!    constant is returned.
impure elemental function                                                  &
   & integrate_to_constant_SubspaceState_ComplexMonomial(monomial,bra,ket, &
   & anharmonic_data) result(output)
  implicit none
  
  type(ComplexMonomial), intent(in)           :: monomial
  class(SubspaceState),  intent(in)           :: bra
  class(SubspaceState),  intent(in), optional :: ket
  type(AnharmonicData),  intent(in)           :: anharmonic_data
  complex(dp)                                 :: output
  
  type(SparseMonomial) :: sparse_monomial
  
  sparse_monomial = SparseMonomial(monomial%modes( &
              & ids        = bra%mode_ids(),       &
              & paired_ids = bra%paired_mode_ids() ))
  
  output = monomial%coefficient            &
       & * bra%integrate( sparse_monomial, &
       &                  ket,             &
       &                  anharmonic_data  )
end function

impure elemental function integrate_to_constant_BasisState_ComplexMonomial( &
   & monomial,bra,ket,subspace,basis,anharmonic_data) result(output)
  implicit none
  
  type(ComplexMonomial),    intent(in)           :: monomial
  class(BasisState),        intent(in)           :: bra
  class(BasisState),        intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  complex(dp)                                    :: output
  
  type(SparseMonomial) :: sparse_monomial
  
  sparse_monomial = SparseMonomial(monomial%modes(                  &
     & ids        = basis%mode_ids(subspace,anharmonic_data),       &
     & paired_ids = basis%paired_mode_ids(subspace,anharmonic_data) ))
  
  output = monomial%coefficient              &
       & * basis%integrate( bra,             &
       &                    sparse_monomial, &
       &                    ket,             &
       &                    subspace,        &
       &                    anharmonic_data  )
end function

impure elemental function integrate_to_constant_BasisStates_ComplexMonomial( &
   & monomial,states,thermal_energy,subspace,basis,anharmonic_data)          &
   & result(output) 
  implicit none
  
  type(ComplexMonomial),    intent(in)    :: monomial
  class(BasisStates),       intent(inout) :: states
  real(dp),                 intent(in)    :: thermal_energy
  type(DegenerateSubspace), intent(in)    :: subspace
  class(SubspaceBasis),     intent(in)    :: basis
  type(AnharmonicData),     intent(in)    :: anharmonic_data
  complex(dp)                             :: output
  
  type(SparseMonomial) :: sparse_monomial
  
  integer     :: cache_location
  complex(dp) :: expectation
  
  sparse_monomial = SparseMonomial(monomial%modes(                  &
     & ids        = basis%mode_ids(subspace,anharmonic_data),       &
     & paired_ids = basis%paired_mode_ids(subspace,anharmonic_data) ))
  
  cache_location = states%expectation_cache%cached_location(sparse_monomial)
  if (cache_location==0) then
    expectation = basis%integrate( states,          &
                                 & thermal_energy,  &
                                 & sparse_monomial, &
                                 & subspace,        &
                                 & anharmonic_data  )
    call states%expectation_cache%cache(sparse_monomial, expectation)
  else
    expectation = states%expectation_cache%cached_expectation(cache_location)
  endif
  
  output = monomial%coefficient * expectation
end function

impure elemental function                                                  &
   & integrate_to_constant_SubspaceState_ComplexPolynomial(polynomial,bra, &
   & ket,anharmonic_data) result(output)
  implicit none
  
  type(ComplexPolynomial), intent(in)           :: polynomial
  class(SubspaceState),    intent(in)           :: bra
  class(SubspaceState),    intent(in), optional :: ket
  type(AnharmonicData),    intent(in)           :: anharmonic_data
  complex(dp)                                   :: output
  
  output = sum(integrate_to_constant( polynomial%terms, &
                                    & bra,              &
                                    & ket,              &
                                    & anharmonic_data   ))
end function

impure elemental function integrate_to_constant_BasisState_ComplexPolynomial( &
   & polynomial,bra,ket,subspace,basis,anharmonic_data) result(output)
  implicit none
  
  type(ComplexPolynomial),  intent(in)           :: polynomial
  class(BasisState),        intent(in)           :: bra
  class(BasisState),        intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  complex(dp)                                    :: output
  
  output = sum(integrate_to_constant( polynomial%terms, &
                                    & bra,              &
                                    & ket,              &
                                    & subspace,         &
                                    & basis,            &
                                    & anharmonic_data   ))
end function

impure elemental function                                                   &
   & integrate_to_constant_BasisStates_ComplexPolynomial(polynomial,states, &
   & thermal_energy,subspace,basis,anharmonic_data) result(output)
  implicit none
  
  type(ComplexPolynomial),  intent(in)    :: polynomial
  class(BasisStates),       intent(inout) :: states
  real(dp),                 intent(in)    :: thermal_energy
  type(DegenerateSubspace), intent(in)    :: subspace
  class(SubspaceBasis),     intent(in)    :: basis
  type(AnharmonicData),     intent(in)    :: anharmonic_data
  complex(dp)                             :: output
  
  output = sum(integrate_to_constant( polynomial%terms, &
                                    & states,           &
                                    & thermal_energy,   &
                                    & subspace,         &
                                    & basis,            &
                                    & anharmonic_data   ))
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
