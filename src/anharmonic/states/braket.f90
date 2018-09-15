! ======================================================================
! Provides functions for objects of class(SubspaceState).
! ======================================================================
! This cannot be acheived with polymorphism since FORTRAN cannot resolve the
!    conflict between f(class(SubspaceState), type(MonomialState)) and
!                     f(type(MonomialState), class(SubspaceState)).
module braket_module
  use common_module
  
  use subspace_state_module
  use monomial_state_module
  use polynomial_state_module
  implicit none
  
  private
  
  public :: braket
  public :: kinetic_energy
  public :: harmonic_potential_energy
  
  interface braket
    module procedure braket_SubspaceStates
    module procedure braket_SubspaceStates_ComplexMonomial
    module procedure braket_SubspaceStates_ComplexPolynomial
  end interface
  
  interface kinetic_energy
    module procedure kinetic_energy_SubspaceStates
  end interface
  
  interface harmonic_potential_energy
    module procedure harmonic_potential_energy_SubspaceStates
  end interface
contains

impure elemental function braket_SubspaceStates(bra,ket) result(output)
  implicit none
  
  class(SubspaceState), intent(in) :: bra
  class(SubspaceState), intent(in) :: ket
  real(dp)                         :: output
  
  select type(bra); type is(MonomialState)
    select type(ket); type is(MonomialState)
      output = braket_MonomialState(bra,ket)
    class default
      call err()
    end select
  type is(PolynomialState)
    select type(ket); type is(PolynomialState)
      output = braket_PolynomialState(bra,ket)
    class default
      call err()
    end select
  class default
    call err()
  end select
end function

impure elemental function braket_SubspaceStates_ComplexMonomial(bra,ket, &
   & monomial,subspace,supercell) result(output)
  implicit none
  
  class(SubspaceState),     intent(in) :: bra
  class(SubspaceState),     intent(in) :: ket
  type(ComplexMonomial),    intent(in) :: monomial
  type(DegenerateSubspace), intent(in) :: subspace
  type(StructureData),      intent(in) :: supercell
  type(ComplexMonomial)                :: output
  
  select type(bra); type is(MonomialState)
    select type(ket); type is(MonomialState)
      output = braket_MonomialState(bra,ket,monomial,subspace,supercell)
    class default
      call err()
    end select
  type is(PolynomialState)
    select type(ket); type is(PolynomialState)
      output = braket_PolynomialState(bra,ket,monomial,subspace,supercell)
    class default
      call err()
    end select
  class default
    call err()
  end select
end function

impure elemental function braket_SubspaceStates_ComplexPolynomial(bra,ket, &
   & polynomial,subspace,supercell) result(output)
  implicit none
  
  class(SubspaceState),     intent(in) :: bra
  class(SubspaceState),     intent(in) :: ket
  type(ComplexPolynomial),  intent(in) :: polynomial
  type(DegenerateSubspace), intent(in) :: subspace
  type(StructureData),      intent(in) :: supercell
  type(ComplexPolynomial)              :: output
  
  output = ComplexPolynomial(braket( bra,              &
                                   & ket,              &
                                   & polynomial%terms, &
                                   & subspace,         &
                                   & supercell         ))
end function

impure elemental function kinetic_energy_SubspaceStates(bra,ket,subspace, &
   & supercell) result(output)
  implicit none
  
  class(SubspaceState),     intent(in) :: bra
  class(SubspaceState),     intent(in) :: ket
  type(DegenerateSubspace), intent(in) :: subspace
  type(StructureData),      intent(in) :: supercell
  real(dp)                             :: output
  
  select type(bra); type is(MonomialState)
    select type(ket); type is(MonomialState)
      output = kinetic_energy_MonomialState(bra,ket,subspace,supercell)
    class default
      call err()
    end select
  type is(PolynomialState)
    select type(ket); type is(PolynomialState)
      output = kinetic_energy_PolynomialState(bra,ket,subspace,supercell)
    class default
      call err()
    end select
  class default
    call err()
  end select
end function

impure elemental function harmonic_potential_energy_SubspaceStates(bra,ket, &
   & subspace,supercell) result(output)
  implicit none
  
  class(SubspaceState),     intent(in) :: bra
  class(SubspaceState),     intent(in) :: ket
  type(DegenerateSubspace), intent(in) :: subspace
  type(StructureData),      intent(in) :: supercell
  real(dp)                             :: output
  
  select type(bra); type is(MonomialState)
    select type(ket); type is(MonomialState)
      output = harmonic_potential_energy_MonomialState( bra,      &
                                                      & ket,      &
                                                      & subspace, &
                                                      & supercell )
    class default
      call err()
    end select
  type is(PolynomialState)
    select type(ket); type is(PolynomialState)
      output = harmonic_potential_energy_PolynomialState( bra,      &
                                                        & ket,      &
                                                        & subspace, &
                                                        & supercell )
    class default
      call err()
    end select
  class default
    call err()
  end select
end function
end module
