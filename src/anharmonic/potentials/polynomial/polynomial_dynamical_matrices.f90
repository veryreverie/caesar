! ======================================================================
! Calculates the contribution to a set of dynamical matrices from
!    a given monomial or polynomial.
! ======================================================================
module caesar_polynomial_dynamical_matrices_module
  use caesar_common_module
  
  use caesar_anharmonic_common_module
  implicit none
  
  private
  
  public :: calculate_dynamical_matrices
  public :: calculate_correction
  
  interface calculate_dynamical_matrices
    ! Calculate the monomial or polynomial's contribution to the effective
    !    dynamical matrix from which the potential can be interpolated in the
    !    large-supercell limit.
    module function calculate_dynamical_matrices_ComplexMonomial(term,    &
       & qpoints,thermal_energy,subspaces,subspace_bases,subspace_states, &
       & subspaces_in_coupling,anharmonic_data) result(output) 
      type(ComplexMonomial),    intent(in)    :: term
      type(QpointData),         intent(in)    :: qpoints(:)
      real(dp),                 intent(in)    :: thermal_energy
      type(DegenerateSubspace), intent(in)    :: subspaces(:)
      class(SubspaceBasis),     intent(in)    :: subspace_bases(:)
      class(BasisStates),       intent(inout) :: subspace_states(:)
      integer,                  intent(in)    :: subspaces_in_coupling(:)
      type(AnharmonicData),     intent(in)    :: anharmonic_data
      type(DynamicalMatrix), allocatable      :: output(:)
    end function
  
    module function calculate_dynamical_matrices_ComplexPolynomial(polynomial,qpoints,thermal_energy,subspaces,subspace_bases,subspace_states,subspaces_in_coupling,anharmonic_data) result(output) 
      type(ComplexPolynomial),  intent(in)    :: polynomial
      type(QpointData),         intent(in)    :: qpoints(:)
      real(dp),                 intent(in)    :: thermal_energy
      type(DegenerateSubspace), intent(in)    :: subspaces(:)
      class(SubspaceBasis),     intent(in)    :: subspace_bases(:)
      class(BasisStates),       intent(inout) :: subspace_states(:)
      integer,                  intent(in)    :: subspaces_in_coupling(:)
      type(AnharmonicData),     intent(in)    :: anharmonic_data
      type(DynamicalMatrix), allocatable      :: output(:)
    end function
  end interface
  
  interface calculate_correction
    ! Calculate the correction due to double-counting when interpolating
    !    polynomials.
    ! It is assumed that a monomial of order n will appear in:
    !    n/2     terms if n is even
    !    (n-1)/2 terms if n is odd.
    ! The double counting of a term is (1-x)<term>, where x is the number of terms
    !    where the term appears.
    module function calculate_correction_ComplexMonomial(monomial,subspaces, &
       & subspace_bases,subspace_states,anharmonic_data) result(output) 
      class(ComplexMonomial),   intent(in)    :: monomial
      type(DegenerateSubspace), intent(in)    :: subspaces(:)
      class(SubspaceBasis),     intent(in)    :: subspace_bases(:)
      class(BasisStates),       intent(inout) :: subspace_states(:)
      type(AnharmonicData),     intent(in)    :: anharmonic_data
      real(dp)                                :: output
    end function
  
    module function calculate_correction_ComplexPolynomial(polynomial, &
       & subspaces,subspace_bases,subspace_states,anharmonic_data)     &
       & result(output) 
      class(ComplexPolynomial), intent(in)    :: polynomial
      type(DegenerateSubspace), intent(in)    :: subspaces(:)
      class(SubspaceBasis),     intent(in)    :: subspace_bases(:)
      class(BasisStates),       intent(inout) :: subspace_states(:)
      type(AnharmonicData),     intent(in)    :: anharmonic_data
      real(dp)                                :: output
    end function
  end interface
end module
