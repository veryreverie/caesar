! ======================================================================
! Basis functions for generating the stress tensor mapping.
! ======================================================================
module caesar_stress_basis_function_module
  use caesar_common_module
  
  use caesar_anharmonic_common_module
  
  use caesar_polynomial_interpolator_module
  use caesar_polynomial_dynamical_matrices_module
  use caesar_polynomial_symmetry_module
  implicit none
  
  private
  
  public :: StressBasisFunction
  public :: generate_stress_basis_functions
  public :: operator(*)
  public :: operator(/)
  
  type, extends(StressBase) :: StressBasisFunction
    ! This will always be 3x3, but is allocatable to avoid stack overflows.
    type(ComplexPolynomial), private, allocatable :: elements_(:,:)
    
    ! A leading coefficient.
    real(dp), private :: coefficient_
  contains
    procedure, public, nopass :: representation => &
                               & representation_StressBasisFunction
    
    procedure, public :: undisplaced_stress => &
                       & undisplaced_stress_StressBasisFunction
    
    procedure, public :: zero_stress => &
                       & zero_stress_StressBasisFunction
    procedure, public :: add_constant => &
                       & add_constant_StressBasisFunction
    
    procedure, public :: stress_RealModeDisplacement => &
                       & stress_RealModeDisplacement_StressBasisFunction
    procedure, public :: stress_ComplexModeDisplacement => &
                       & stress_ComplexModeDisplacement_StressBasisFunction
    
    procedure, public :: braket_SubspaceBraKet => &
                       & braket_SubspaceBraKet_StressBasisFunction
    procedure, public :: braket_BasisState => &
                       & braket_BasisState_StressBasisFunction
    procedure, public :: braket_BasisStates => &
                       & braket_BasisStates_StressBasisFunction
    
    procedure, public :: harmonic_expectation => &
                       & harmonic_expectation_StressBasisFunction
    
    procedure, public :: simplify => simplify_StressBasisFunction
    
    procedure, public :: interpolate_coefficients => &
                       & interpolate_coefficients_StressBasisFunction
    procedure, public :: calculate_dynamical_matrices => &
                       & calculate_dynamical_matrices_StressBasisFunction
    procedure, public :: stress_correction => &
                       & stress_correction_StressBasisFunction
    
    ! I/O.
    procedure, public :: read  => read_StressBasisFunction
    procedure, public :: write => write_StressBasisFunction
  end type
  
  interface StressBasisFunction
    ! Constructor.
    module function new_StressBasisFunction(elements,coefficient) result(this) 
      type(ComplexPolynomial), intent(in)           :: elements(3,3)
      real(dp),                intent(in), optional :: coefficient
      type(StressBasisFunction)                     :: this
    end function
  end interface
  
  interface
    ! Type representation.
    impure elemental module function representation_StressBasisFunction() &
       & result(output) 
      type(String) :: output
    end function
  end interface
  
  interface generate_stress_basis_functions
    ! Generate basis functions.
    module function generate_stress_basis_functions_SubspaceCombination(      &
       & subspace_combination,maximum_coupling_order,structure,complex_modes, &
       & qpoints,subspaces,degenerate_symmetries,vscf_basis_functions_only,   &
       & qpoint_symmetry_groups,subspace_qpoint_stars,logfile) result(output) 
      type(SubspaceCombination), intent(in)    :: subspace_combination
      integer,                   intent(in)    :: maximum_coupling_order
      type(StructureData),       intent(in)    :: structure
      type(ComplexMode),         intent(in)    :: complex_modes(:)
      type(QpointData),          intent(in)    :: qpoints(:)
      type(DegenerateSubspace),  intent(in)    :: subspaces(:)
      type(DegenerateSymmetry),  intent(in)    :: degenerate_symmetries(:)
      logical,                   intent(in)    :: vscf_basis_functions_only
      type(Group),               intent(in)    :: qpoint_symmetry_groups(:)
      type(SubspaceQpointStars), intent(in)    :: subspace_qpoint_stars(:)
      type(OFile),               intent(inout) :: logfile
      type(StressBasisFunction), allocatable   :: output(:)
    end function
  end interface
  
  interface
    ! Given an orthogonal matrix U, s.t. U^n=I, returns the matrix
    !    H = sum(j=0,n-1) U^j / n
    ! H is symmetric, and is a projection matrix which projects out the
    !    eigenvectors of U with eigenvalue 1.
    ! Uses Horner's rule to calculate the sum with minimal matrix multiplication.
    module function projection_matrix(input,order) result(output) 
      type(RealMatrix), intent(in) :: input
      integer,          intent(in) :: order
      type(RealMatrix)             :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Simplify the basis function.
    ! ----------------------------------------------------------------------
    impure elemental module subroutine simplify_StressBasisFunction(this) 
      class(StressBasisFunction), intent(inout) :: this
    end subroutine
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Evaluate the stress due to the basis function.
    ! ----------------------------------------------------------------------
    impure elemental module function stress_RealModeDisplacement_StressBasisFunction(   this,displacement) result(output) 
      class(StressBasisFunction), intent(in) :: this
      type(RealModeDisplacement), intent(in) :: displacement
      type(RealMatrix)                       :: output
    end function
  end interface
  
  interface
    impure elemental module function stress_ComplexModeDisplacement_StressBasisFunction(   this,displacement) result(output) 
      class(StressBasisFunction),    intent(in) :: this
      type(ComplexModeDisplacement), intent(in) :: displacement
      type(ComplexMatrix)                       :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Integrate the basis function between two states.
    ! ----------------------------------------------------------------------
    impure elemental module subroutine braket_SubspaceBraKet_StressBasisFunction(this,braket,whole_subspace,anharmonic_data) 
      class(StressBasisFunction), intent(inout)        :: this
      class(SubspaceBraKet),      intent(in)           :: braket
      logical,                    intent(in), optional :: whole_subspace
      type(AnharmonicData),       intent(in)           :: anharmonic_data
    end subroutine
  end interface
  
  interface
    impure elemental module subroutine braket_BasisState_StressBasisFunction(this,bra,ket,subspace,subspace_basis,whole_subspace,anharmonic_data) 
      class(StressBasisFunction), intent(inout)        :: this
      class(BasisState),          intent(in)           :: bra
      class(BasisState),          intent(in), optional :: ket
      type(DegenerateSubspace),   intent(in)           :: subspace
      class(SubspaceBasis),       intent(in)           :: subspace_basis
      logical,                    intent(in), optional :: whole_subspace
      type(AnharmonicData),       intent(in)           :: anharmonic_data
    end subroutine
  end interface
  
  interface
    impure elemental module subroutine braket_BasisStates_StressBasisFunction(this,states,subspace,subspace_basis,whole_subspace,anharmonic_data) 
      class(StressBasisFunction), intent(inout)        :: this
      class(BasisStates),         intent(inout)        :: states
      type(DegenerateSubspace),   intent(in)           :: subspace
      class(SubspaceBasis),       intent(in)           :: subspace_basis
      logical,                    intent(in), optional :: whole_subspace
      type(AnharmonicData),       intent(in)           :: anharmonic_data
    end subroutine
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Returns the thermal expectation of the basis function.
    ! ----------------------------------------------------------------------
    impure elemental module function harmonic_expectation_StressBasisFunction(this,frequency,thermal_energy,supercell_size,anharmonic_data) result(output) 
      class(StressBasisFunction), intent(in) :: this
      real(dp),                   intent(in) :: frequency
      real(dp),                   intent(in) :: thermal_energy
      integer,                    intent(in) :: supercell_size
      type(AnharmonicData),       intent(in) :: anharmonic_data
      type(RealMatrix)                       :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Returns the stress at zero displacement.
    ! ----------------------------------------------------------------------
    impure elemental module function undisplaced_stress_StressBasisFunction(this) result(output) 
      class(StressBasisFunction), intent(in) :: this
      type(RealMatrix)                       :: output
    end function
  end interface
  
  interface
    impure elemental module subroutine zero_stress_StressBasisFunction(this) 
      class(StressBasisFunction), intent(inout) :: this
    end subroutine
  end interface
  
  interface
    impure elemental module subroutine add_constant_StressBasisFunction(this,input) 
      class(StressBasisFunction), intent(inout) :: this
      type(RealMatrix),           intent(in)    :: input
    end subroutine
  end interface
  
  interface operator(*)
    ! ----------------------------------------------------------------------
    ! Arithmetic.
    ! ----------------------------------------------------------------------
    impure elemental module function multiply_StressBasisFunction_real(this, &
       & that) result(output) 
      type(StressBasisFunction), intent(in) :: this
      real(dp),                  intent(in) :: that
      type(StressBasisFunction)             :: output
    end function
  
    impure elemental module function multiply_real_StressBasisFunction(this, &
       & that) result(output) 
      real(dp),                  intent(in) :: this
      type(StressBasisFunction), intent(in) :: that
      type(StressBasisFunction)             :: output
    end function
  end interface
  
  interface operator(/)
    impure elemental module function divide_StressBasisFunction_real(this, &
       & that) result(output) 
      type(StressBasisFunction), intent(in) :: this
      real(dp),                  intent(in) :: that
      type(StressBasisFunction)             :: output
    end function
  end interface
  
  interface
    ! Calculate the contribution to a given monomial from the interpolation of
    !    this basis function.
    ! The result is given as a cartesian tensor.
    impure elemental module function                                 &
       & interpolate_coefficients_StressBasisFunction(this,monomial, &
       & interpolator) result(output) 
      class(StressBasisFunction),   intent(in) :: this
      type(ComplexMonomial),        intent(in) :: monomial
      type(PolynomialInterpolator), intent(in) :: interpolator
      type(ComplexMatrix)                      :: output
    end function
  end interface
  
  interface
    ! Calculate this basis function's contribution to the effective stress
    !    dynamical matrix from which the stress can be interpolated in the
    !    large-supercell limit.
    module function calculate_dynamical_matrices_StressBasisFunction(this, &
       & qpoints,thermal_energy,subspaces,subspace_bases,subspace_states,  &
       & subspaces_in_coupling,anharmonic_data) result(output) 
      class(StressBasisFunction), intent(in)    :: this
      type(QpointData),           intent(in)    :: qpoints(:)
      real(dp),                   intent(in)    :: thermal_energy
      type(DegenerateSubspace),   intent(in)    :: subspaces(:)
      class(SubspaceBasis),       intent(in)    :: subspace_bases(:)
      class(BasisStates),         intent(inout) :: subspace_states(:)
      integer,                    intent(in)    :: subspaces_in_coupling(:)
      type(AnharmonicData),       intent(in)    :: anharmonic_data
      type(StressDynamicalMatrix), allocatable  :: output(:)
    end function
  end interface
  
  interface
    ! Calculate the correction due to double counting
    !    for the interpolated stress.
    module function stress_correction_StressBasisFunction(this,subspaces, &
       & subspace_bases,subspace_states,anharmonic_data) result(output) 
      class(StressBasisFunction), intent(in)    :: this
      type(DegenerateSubspace),   intent(in)    :: subspaces(:)
      class(SubspaceBasis),       intent(in)    :: subspace_bases(:)
      class(BasisStates),         intent(inout) :: subspace_states(:)
      type(AnharmonicData),       intent(in)    :: anharmonic_data
      type(RealMatrix)                          :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_StressBasisFunction(this,input) 
      class(StressBasisFunction), intent(out) :: this
      type(String),               intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_StressBasisFunction(this) result(output) 
      class(StressBasisFunction), intent(in) :: this
      type(String), allocatable              :: output(:)
    end function
  end interface
  
  interface StressBasisFunction
    module function new_StressBasisFunction_Strings(input) result(this) 
      type(String), intent(in)  :: input(:)
      type(StressBasisFunction) :: this
    end function
  
    impure elemental module function new_StressBasisFunction_StringArray(input) result(this) 
      type(StringArray), intent(in) :: input
      type(StressBasisFunction)     :: this
    end function
  end interface
end module
