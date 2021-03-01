! ======================================================================
! Calculates potential energy or potential stress.
! ======================================================================
module caesar_braket_module
  use caesar_common_module
  
  use caesar_subspace_state_module
  use caesar_subspace_braket_module
  use caesar_basis_state_module
  use caesar_basis_states_module
  use caesar_anharmonic_data_module
  use caesar_stress_prefactors_module
  use caesar_abstract_classes_module
  use caesar_sparse_monomial_module
  implicit none
  
  private
  
  public :: integrate
  public :: integrate_to_constant
  public :: harmonic_observables
  public :: effective_harmonic_observables
  
  interface integrate
    ! Integrate a the parts of a ComplexMonomial in a given subspace.
    impure elemental module subroutine integrate_SubspaceBraKet(monomial, &
       & braket,anharmonic_data) 
      type(ComplexMonomial), intent(inout) :: monomial
      class(SubspaceBraKet), intent(in)    :: braket
      type(AnharmonicData),  intent(in)    :: anharmonic_data
    end subroutine
  
    impure elemental module subroutine integrate_BasisState(monomial,bra, &
       & ket,subspace,basis,anharmonic_data) 
      type(ComplexMonomial),    intent(inout)        :: monomial
      class(BasisState),        intent(in)           :: bra
      class(BasisState),        intent(in), optional :: ket
      type(DegenerateSubspace), intent(in)           :: subspace
      class(SubspaceBasis),     intent(in)           :: basis
      type(AnharmonicData),     intent(in)           :: anharmonic_data
    end subroutine
  
    impure elemental module subroutine integrate_BasisStates(monomial, &
       & states,subspace,basis,anharmonic_data) 
      type(ComplexMonomial),    intent(inout) :: monomial
      class(BasisStates),       intent(inout) :: states
      type(DegenerateSubspace), intent(in)    :: subspace
      class(SubspaceBasis),     intent(in)    :: basis
      type(AnharmonicData),     intent(in)    :: anharmonic_data
    end subroutine
  end interface
  
  interface integrate_to_constant
    ! As integrate, except that the whole ComplexMonomial must be in the given
    !    subspace.
    ! Since the resulting ComplexMonomial would just be a constant, just the
    !    constant is returned.
    impure elemental module function                                                    integrate_to_constant_SubspaceBraKet_ComplexMonomial(monomial,braket,anharmonic_data) result(output) 
      type(ComplexMonomial), intent(in) :: monomial
      class(SubspaceBraKet), intent(in) :: braket
      type(AnharmonicData),  intent(in) :: anharmonic_data
      complex(dp)                       :: output
    end function
  
    impure elemental module function integrate_to_constant_BasisState_ComplexMonomial(   monomial,bra,ket,subspace,basis,anharmonic_data) result(output) 
      type(ComplexMonomial),    intent(in)           :: monomial
      class(BasisState),        intent(in)           :: bra
      class(BasisState),        intent(in), optional :: ket
      type(DegenerateSubspace), intent(in)           :: subspace
      class(SubspaceBasis),     intent(in)           :: basis
      type(AnharmonicData),     intent(in)           :: anharmonic_data
      complex(dp)                                    :: output
    end function
  
    impure elemental module function integrate_to_constant_BasisStates_ComplexMonomial(   monomial,states,subspace,basis,anharmonic_data) result(output) 
      type(ComplexMonomial),    intent(in)    :: monomial
      class(BasisStates),       intent(inout) :: states
      type(DegenerateSubspace), intent(in)    :: subspace
      class(SubspaceBasis),     intent(in)    :: basis
      type(AnharmonicData),     intent(in)    :: anharmonic_data
      complex(dp)                             :: output
    end function
  
    impure elemental module function                                                 integrate_to_constant_SubspaceBraKet_ComplexPolynomial(polynomial,braket,anharmonic_data) result(output) 
      type(ComplexPolynomial), intent(in) :: polynomial
      class(SubspaceBraKet),   intent(in) :: braket
      type(AnharmonicData),    intent(in) :: anharmonic_data
      complex(dp)                         :: output
    end function
  
    impure elemental module function integrate_to_constant_BasisState_ComplexPolynomial(   polynomial,bra,ket,subspace,basis,anharmonic_data) result(output) 
      type(ComplexPolynomial),  intent(in)           :: polynomial
      class(BasisState),        intent(in)           :: bra
      class(BasisState),        intent(in), optional :: ket
      type(DegenerateSubspace), intent(in)           :: subspace
      class(SubspaceBasis),     intent(in)           :: basis
      type(AnharmonicData),     intent(in)           :: anharmonic_data
      complex(dp)                                    :: output
    end function
  
    impure elemental module function                                                     integrate_to_constant_BasisStates_ComplexPolynomial(polynomial,states,subspace,basis,anharmonic_data) result(output) 
      type(ComplexPolynomial),  intent(in)    :: polynomial
      class(BasisStates),       intent(inout) :: states
      type(DegenerateSubspace), intent(in)    :: subspace
      class(SubspaceBasis),     intent(in)    :: basis
      type(AnharmonicData),     intent(in)    :: anharmonic_data
      complex(dp)                             :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Thermodynamic quantities derived from arbitrary potentials and stresses.
    ! ----------------------------------------------------------------------
    ! Calculate observables for harmonic basis, using a harmonic potential.
    ! N.B. the result is extensive, so will in general need to be normalised
    !    to be per unit cell or similar.
    impure elemental module function harmonic_observables(thermal_energy, &
       & stress,stress_prefactor,frequency,num_dimensions,supercell_size, &
       & anharmonic_data) result(output) 
      real(dp),             intent(in)           :: thermal_energy
      class(StressBase),    intent(in), optional :: stress
      type(RealMatrix),     intent(in), optional :: stress_prefactor
      real(dp),             intent(in)           :: frequency
      integer,              intent(in)           :: num_dimensions
      integer,              intent(in)           :: supercell_size
      type(AnharmonicData), intent(in)           :: anharmonic_data
      type(ThermodynamicData)                    :: output
    end function
  end interface
  
  interface
    ! Calculate observables for harmonic basis, using full potential.
    ! N.B. the result is extensive, so will in general need to be normalised
    !    to be per unit cell or similar.
    impure elemental module function effective_harmonic_observables(thermal_energy,potential,stress,stress_prefactor,frequency,num_dimensions,supercell_size,anharmonic_data) result(output) 
      real(dp),             intent(in)           :: thermal_energy
      class(PotentialBase), intent(in)           :: potential
      class(StressBase),    intent(in), optional :: stress
      type(RealMatrix),     intent(in), optional :: stress_prefactor
      real(dp),             intent(in)           :: frequency
      integer,              intent(in)           :: num_dimensions
      integer,              intent(in)           :: supercell_size
      type(AnharmonicData), intent(in)           :: anharmonic_data
      type(ThermodynamicData)                    :: output
    end function
  end interface
end module
