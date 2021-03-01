submodule (caesar_braket_module) caesar_braket_submodule
  use caesar_anharmonic_common_module
contains

module procedure integrate_SubspaceBraKet
  type(ComplexUnivariate), allocatable :: unintegrated_modes(:)
  type(ComplexUnivariate), allocatable :: integrated_modes(:)
  
  unintegrated_modes = monomial%modes(exclude_ids = braket%mode_ids)
  
  integrated_modes = monomial%modes(       &
     & ids        = braket%mode_ids,       &
     & paired_ids = braket%paired_mode_ids )
  
  call monomial%set_modes(unintegrated_modes)
  
  if (any(integrated_modes%total_power()/=0)) then
    monomial%coefficient = monomial%coefficient                 &
                       & * braket%integrate(                    &
                       &      SparseMonomial(integrated_modes), &
                       &      anharmonic_data                   )
  endif
end procedure

module procedure integrate_BasisState
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
end procedure

module procedure integrate_BasisStates
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
                                   & sparse,         &
                                   & subspace,       &
                                   & anharmonic_data )
      call states%expectation_cache%cache(sparse, expectation)
    else
      expectation = states%expectation_cache%cached_expectation(cache_location)
    endif
    monomial%coefficient = monomial%coefficient * expectation
  endif
end procedure

module procedure integrate_to_constant_SubspaceBraKet_ComplexMonomial
  type(SparseMonomial) :: sparse_monomial
  
  sparse_monomial%modes = monomial%modes(  &
     & ids        = braket%mode_ids,       &
     & paired_ids = braket%paired_mode_ids )
  
  output = monomial%coefficient               &
       & * braket%integrate( sparse_monomial, &
       &                     anharmonic_data  )
end procedure

module procedure integrate_to_constant_BasisState_ComplexMonomial
  type(SparseMonomial) :: sparse_monomial
  
  sparse_monomial%modes = monomial%modes(                           &
     & ids        = basis%mode_ids(subspace,anharmonic_data),       &
     & paired_ids = basis%paired_mode_ids(subspace,anharmonic_data) )
  
  output = monomial%coefficient              &
       & * basis%integrate( bra,             &
       &                    sparse_monomial, &
       &                    ket,             &
       &                    subspace,        &
       &                    anharmonic_data  )
end procedure

module procedure integrate_to_constant_BasisStates_ComplexMonomial
  type(SparseMonomial) :: sparse_monomial
  
  integer     :: cache_location
  complex(dp) :: expectation
  
  sparse_monomial%modes = monomial%modes(                           &
     & ids        = basis%mode_ids(subspace,anharmonic_data),       &
     & paired_ids = basis%paired_mode_ids(subspace,anharmonic_data) )
  
  cache_location = states%expectation_cache%cached_location(sparse_monomial)
  if (cache_location==0) then
    expectation = basis%integrate( states,          &
                                 & sparse_monomial, &
                                 & subspace,        &
                                 & anharmonic_data  )
    call states%expectation_cache%cache(sparse_monomial, expectation)
  else
    expectation = states%expectation_cache%cached_expectation(cache_location)
  endif
  
  output = monomial%coefficient * expectation
end procedure

module procedure integrate_to_constant_SubspaceBraKet_ComplexPolynomial
  integer :: i
  
  output = 0
  do i=1,size(polynomial%terms)
    output = output + integrate_to_constant( polynomial%terms(i), &
                                           & braket,              &
                                           & anharmonic_data      )
  enddo
end procedure

module procedure integrate_to_constant_BasisState_ComplexPolynomial
  integer :: i
  
  output = 0
  do i=1,size(polynomial%terms)
    output = output + integrate_to_constant( polynomial%terms(i), &
                                           & bra,                 &
                                           & ket,                 &
                                           & subspace,            &
                                           & basis,               &
                                           & anharmonic_data      )
  enddo
end procedure

module procedure integrate_to_constant_BasisStates_ComplexPolynomial
  integer :: i
  
  output = 0
  do i=1,size(polynomial%terms)
    output = output + integrate_to_constant( polynomial%terms(i), &
                                           & states,              &
                                           & subspace,            &
                                           & basis,               &
                                           & anharmonic_data      )
  enddo
end procedure

module procedure harmonic_observables
  type(RealMatrix), allocatable :: potential_stress
  real(dp),         allocatable :: primitive_volume
  
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
    primitive_volume = anharmonic_data%structure%volume
  endif
  
  output = ThermodynamicData( thermal_energy,     &
       &                      frequency,          &
       &                      stress_prefactor,   &
       &                      potential_stress,   &
       &                      primitive_volume  ) &
       & * num_dimensions
end procedure

module procedure effective_harmonic_observables
  real(dp) :: harmonic_expectation
  real(dp) :: potential_expectation
  
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
end procedure
end submodule
