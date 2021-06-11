! ======================================================================
! Basis functions for generating the Born-Oppenheimer surface mapping.
! ======================================================================
! N.B. basis functions are initially generated with both a real and complex
!    representation.
! The real representation is used for fitting, and the complex representation
!    is used for calculating overlap integrals.
module caesar_basis_function_module
  use caesar_common_module
  
  use caesar_anharmonic_common_module
  
  use caesar_polynomial_interpolator_module
  use caesar_polynomial_dynamical_matrices_module
  use caesar_polynomial_symmetry_module
  implicit none
  
  private
  
  public :: BasisFunction
  public :: BasisFunctions
  public :: generate_basis_functions
  public :: optimise
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: operator(-)
  
  type, extends(PotentialBase) :: BasisFunction
    ! The basis function in complex co-ordinates.
    type(ComplexPolynomial), private :: complex_representation_
    
    ! A leading coefficient.
    real(dp), private :: coefficient_
  contains
    procedure, public, nopass :: representation => &
                               & representation_BasisFunction
    
    procedure, public :: complex_representation
    
    procedure, public :: undisplaced_energy => undisplaced_energy_BasisFunction
    
    procedure, public :: zero_energy => zero_energy_BasisFunction
    procedure, public :: add_constant => add_constant_BasisFunction
    
    procedure, public :: energy_RealModeDisplacement => &
                       & energy_RealModeDisplacement_BasisFunction
    procedure, public :: energy_ComplexModeDisplacement => &
                       & energy_ComplexModeDisplacement_BasisFunction
    
    procedure, public :: force_RealModeDisplacement => &
                       & force_RealModeDisplacement_BasisFunction
    procedure, public :: force_ComplexModeDisplacement => &
                       & force_ComplexModeDisplacement_BasisFunction
    
    procedure, public :: braket_SubspaceBraKet => &
                       & braket_SubspaceBraKet_BasisFunction
    procedure, public :: braket_BasisState => &
                       & braket_BasisState_BasisFunction
    procedure, public :: braket_BasisStates => &
                       & braket_BasisStates_BasisFunction
    
    procedure, public :: harmonic_expectation => &
                       & harmonic_expectation_BasisFunction
    
    procedure, public :: potential_energy_SubspaceBraKet
    procedure, public :: potential_energy_BasisState
    
    procedure, public :: simplify => simplify_BasisFunction
    
    procedure, public :: power => power_BasisFunction
    procedure, public :: coefficient => coefficient_BasisFunction
    procedure, public :: set_coefficient => set_coefficient_BasisFunction
    
    procedure, public :: interpolate_coefficient => &
                       & interpolate_coefficient_BasisFunction
    procedure, public :: add_interpolated_contribution => &
                       & add_interpolated_contribution_BasisFunction
    procedure, public :: add_harmonic_contribution => &
                       & add_harmonic_contribution_BasisFunction
    procedure, public :: calculate_dynamical_matrices => &
                       & calculate_dynamical_matrices_BasisFunction
    procedure, public :: energy_correction => &
                       & energy_correction_BasisFunction
    
    procedure, public :: terms => terms_BasisFunction
    
    procedure, public :: read  => read_BasisFunction
    procedure, public :: write => write_BasisFunction
  end type
  
  ! An array of type BasisFunction.
  type :: BasisFunctions
    type(BasisFunction), allocatable :: basis_functions(:)
  end type
  
  interface BasisFunction
    ! ----------------------------------------------------------------------
    ! Constructor.
    ! ----------------------------------------------------------------------
    impure elemental module function new_BasisFunction(complex_representation,coefficient) result(this) 
      type(ComplexPolynomial), intent(in)           :: complex_representation
      real(dp),                intent(in), optional :: coefficient
      type(BasisFunction)                           :: this
    end function
  
    recursive module function new_BasisFunction_PotentialBase(input) &
       & result(this) 
      class(PotentialBase), intent(in) :: input
      type(BasisFunction)              :: this
    end function
  end interface
  
  interface
    ! Type representation.
    impure elemental module function representation_BasisFunction() &
       & result(output) 
      type(String) :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Getters.
    ! ----------------------------------------------------------------------
    impure elemental module function complex_representation(this) &
       & result(output) 
      class(BasisFunction), intent(in) :: this
      type(ComplexPolynomial)          :: output
    end function
  end interface
  
  interface generate_basis_functions
    ! ----------------------------------------------------------------------
    ! Generate basis functions.
    ! ----------------------------------------------------------------------
    module function generate_basis_functions_SubspaceCombination(         &
       & subspace_combination,max_subspace_coupling,max_qpoint_coupling,  &
       & structure,complex_modes,qpoints,subspaces,degenerate_symmetries, &
       & vscf_basis_functions_only,qpoint_symmetry_groups,                &
       & subspace_qpoint_stars,logfile) result(output) 
      type(SubspaceCombination), intent(in)    :: subspace_combination
      integer,                   intent(in)    :: max_subspace_coupling
      integer,                   intent(in)    :: max_qpoint_coupling
      type(StructureData),       intent(in)    :: structure
      type(ComplexMode),         intent(in)    :: complex_modes(:)
      type(QpointData),          intent(in)    :: qpoints(:)
      type(DegenerateSubspace),  intent(in)    :: subspaces(:)
      type(DegenerateSymmetry),  intent(in)    :: degenerate_symmetries(:)
      logical,                   intent(in)    :: vscf_basis_functions_only
      type(Group),               intent(in)    :: qpoint_symmetry_groups(:)
      type(SubspaceQpointStars), intent(in)    :: subspace_qpoint_stars(:)
      type(OFile),               intent(inout) :: logfile
      type(BasisFunction), allocatable         :: output(:)
    end function
  end interface
  
  interface
    ! Given a unitary matrix U, s.t. U^n=I, returns the matrix
    !    H = sum(j=0,n-1) U^j / n
    ! H is Hermitian, and is a projection matrix which projects out the
    !    eigenvectors of U with eigenvalue 1.
    ! Uses Horner's rule to calculate the sum with minimal matrix multiplication.
    module function projection_matrix(input,order) result(output) 
      type(ComplexMatrix), intent(in) :: input
      integer,             intent(in) :: order
      type(ComplexMatrix)             :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Simplify the basis function.
    ! ----------------------------------------------------------------------
    impure elemental module subroutine simplify_BasisFunction(this) 
      class(BasisFunction), intent(inout) :: this
    end subroutine
  end interface
  
  interface optimise
    module function optimise_BasisFunctions(input,subspace,subspace_basis, &
       & old_subspace_potential,anharmonic_data) result(output)
      type(BasisFunction),       intent(in)           :: input(:)
      type(DegenerateSubspace),  intent(in)           :: subspace
      class(SubspaceBasis),      intent(in)           :: subspace_basis
      class(PotentialData),      intent(in), optional :: old_subspace_potential
      type(AnharmonicData),      intent(in)           :: anharmonic_data
      type(BasisFunction), allocatable                :: output(:)
    end function
  end interface
  
  interface
    ! Construct the basis functions at a given order.
    module function fit_basis_functions(monomials,basis_functions) &
       & result(output) 
      type(ComplexMonomial),   intent(in) :: monomials(:)
      type(ComplexPolynomial), intent(in) :: basis_functions(:) 
      type(BasisFunction), allocatable    :: output(:)
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Evaluate the energy and forces due to the basis function.
    ! ----------------------------------------------------------------------
    impure elemental module function energy_RealModeDisplacement_BasisFunction(this,displacement) result(output) 
      class(BasisFunction),       intent(in) :: this
      type(RealModeDisplacement), intent(in) :: displacement
      real(dp)                               :: output
    end function
  end interface
  
  interface
    impure elemental module function energy_ComplexModeDisplacement_BasisFunction(this,displacement) result(output) 
      class(BasisFunction),          intent(in) :: this
      type(ComplexModeDisplacement), intent(in) :: displacement
      complex(dp)                               :: output
    end function
  end interface
  
  interface
    impure elemental module function force_RealModeDisplacement_BasisFunction(this,displacement) result(output) 
      class(BasisFunction),       intent(in) :: this
      type(RealModeDisplacement), intent(in) :: displacement
      type(RealModeForce)                    :: output
    end function
  end interface
  
  interface
    impure elemental module function force_ComplexModeDisplacement_BasisFunction(this,displacement) result(output) 
      class(BasisFunction),          intent(in) :: this
      type(ComplexModeDisplacement), intent(in) :: displacement
      type(ComplexModeForce)                    :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Integrate the basis function between two states.
    ! ----------------------------------------------------------------------
    impure elemental module subroutine braket_SubspaceBraKet_BasisFunction(this,braket,whole_subspace,anharmonic_data) 
      class(BasisFunction),  intent(inout)        :: this
      class(SubspaceBraKet), intent(in)           :: braket
      logical,               intent(in), optional :: whole_subspace
      type(AnharmonicData),  intent(in)           :: anharmonic_data
    end subroutine
  end interface
  
  interface
    impure elemental module subroutine braket_BasisState_BasisFunction(this, &
       & bra,ket,subspace,subspace_basis,whole_subspace,anharmonic_data) 
      class(BasisFunction),     intent(inout)        :: this
      class(BasisState),        intent(in)           :: bra
      class(BasisState),        intent(in), optional :: ket
      type(DegenerateSubspace), intent(in)           :: subspace
      class(SubspaceBasis),     intent(in)           :: subspace_basis
      logical,                  intent(in), optional :: whole_subspace
      type(AnharmonicData),     intent(in)           :: anharmonic_data
    end subroutine
  end interface
  
  interface
    impure elemental module subroutine braket_BasisStates_BasisFunction(this,states,subspace,subspace_basis,whole_subspace,anharmonic_data) 
      class(BasisFunction),     intent(inout)        :: this
      class(BasisStates),       intent(inout)        :: states
      type(DegenerateSubspace), intent(in)           :: subspace
      class(SubspaceBasis),     intent(in)           :: subspace_basis
      logical,                  intent(in), optional :: whole_subspace
      type(AnharmonicData),     intent(in)           :: anharmonic_data
    end subroutine
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Returns the thermal expectation of the basis function.
    ! ----------------------------------------------------------------------
    impure elemental module function harmonic_expectation_BasisFunction(this,frequency,thermal_energy,supercell_size,anharmonic_data) result(output) 
      class(BasisFunction), intent(in) :: this
      real(dp),             intent(in) :: frequency
      real(dp),             intent(in) :: thermal_energy
      integer,              intent(in) :: supercell_size
      type(AnharmonicData), intent(in) :: anharmonic_data
      real(dp)                         :: output
    end function
  end interface
  
  interface
    module function potential_energy_SubspaceBraKet(this,braket, &
       & anharmonic_data) result(output) 
      class(BasisFunction),  intent(in) :: this
      class(SubspaceBraKet), intent(in) :: braket
      type(AnharmonicData),  intent(in) :: anharmonic_data
      real(dp)                          :: output
    end function
  end interface
  
  interface
    module function potential_energy_BasisState(this,bra,ket,subspace, &
       & subspace_basis,anharmonic_data) result(output) 
      class(BasisFunction),     intent(in)           :: this
      class(BasisState),        intent(in)           :: bra
      class(BasisState),        intent(in), optional :: ket
      type(DegenerateSubspace), intent(in)           :: subspace
      class(SubspaceBasis),     intent(in)           :: subspace_basis
      type(AnharmonicData),     intent(in)           :: anharmonic_data
      real(dp)                                       :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Returns the energy at zero displacement.
    ! ----------------------------------------------------------------------
    impure elemental module function undisplaced_energy_BasisFunction(this) &
       & result(output) 
      class(BasisFunction), intent(in) :: this
      real(dp)                         :: output
    end function
  end interface
  
  interface
    impure elemental module subroutine zero_energy_BasisFunction(this) 
      class(BasisFunction), intent(inout) :: this
    end subroutine
  end interface
  
  interface
    impure elemental module subroutine add_constant_BasisFunction(this,input) 
      class(BasisFunction), intent(inout) :: this
      real(dp),             intent(in)    :: input
    end subroutine
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Converts the basis function to and from a single coefficient.
    ! These methods allow for simpler linear algebra with basis functions.
    ! ----------------------------------------------------------------------
    impure elemental module function power_BasisFunction(this) result(output) 
      class(BasisFunction), intent(in) :: this
      integer                          :: output
    end function
  end interface
  
  interface
    impure elemental module function coefficient_BasisFunction(this) &
       & result(output) 
      class(BasisFunction), intent(in) :: this
      real(dp)                         :: output
    end function
  end interface
  
  interface
    impure elemental module subroutine set_coefficient_BasisFunction(this, &
       & coefficient) 
      class(BasisFunction), intent(inout) :: this
      real(dp),             intent(in)    :: coefficient
    end subroutine
  end interface
  
  interface operator(*)
    ! ----------------------------------------------------------------------
    ! Arithmetic.
    ! ----------------------------------------------------------------------
    impure elemental module function multiply_BasisFunction_real(this,that) &
       & result(output) 
      type(BasisFunction), intent(in) :: this
      real(dp),            intent(in) :: that
      type(BasisFunction)             :: output
    end function
  
    impure elemental module function multiply_real_BasisFunction(this,that) &
       & result(output) 
      real(dp),            intent(in) :: this
      type(BasisFunction), intent(in) :: that
      type(BasisFunction)             :: output
    end function
  end interface
  
  interface operator(/)
    impure elemental module function divide_BasisFunction_real(this,that) &
       & result(output) 
      type(BasisFunction), intent(in) :: this
      real(dp),            intent(in) :: that
      type(BasisFunction)             :: output
    end function
  end interface
  
  interface operator(+)
    impure elemental module function add_BasisFunction_BasisFunction(this, &
       & that) result(output) 
      type(BasisFunction), intent(in) :: this
      type(BasisFunction), intent(in) :: that
      type(BasisFunction)             :: output
    end function
  end interface
  
  interface operator(-)
    impure elemental module function negative_BasisFunction(this) &
       & result(output) 
      type(BasisFunction), intent(in) :: this
      type(BasisFunction)             :: output
    end function
  
    impure elemental module function subtract_BasisFunction_BasisFunction(this,that) result(output) 
      type(BasisFunction), intent(in) :: this
      type(BasisFunction), intent(in) :: that
      type(BasisFunction)             :: output
    end function
  end interface
  
  interface
    ! Calculate the contribution to a given monomial from the interpolation of
    !    this basis function.
    impure elemental module function interpolate_coefficient_BasisFunction( &
       & this,monomial,interpolator) result(output) 
      class(BasisFunction),         intent(in) :: this
      type(ComplexMonomial),        intent(in) :: monomial
      type(PolynomialInterpolator), intent(in) :: interpolator
      complex(dp)                              :: output
    end function
  end interface
  
  interface
    ! Calculate the contribution to this basis function from
    !    another basis function,
    !    and add this to this basis function's coefficient.
    module subroutine add_interpolated_contribution_BasisFunction(this, &
        basis_function,interpolator) 
      class(BasisFunction),         intent(inout) :: this
      type(BasisFunction),          intent(in)    :: basis_function
      type(PolynomialInterpolator), intent(in)    :: interpolator
    end subroutine
  end interface
  
  interface
    ! Calculate the contribution to this basis function from a set of harmonic
    !    dynamical matrices.
    module subroutine add_harmonic_contribution_BasisFunction(this, &
       & dynamical_matrices,anharmonic_data) 
      class(BasisFunction),  intent(inout) :: this
      type(DynamicalMatrix), intent(in)    :: dynamical_matrices(:)
      type(AnharmonicData),  intent(in)    :: anharmonic_data
    end subroutine
  end interface
  
  interface
    ! Calculate this basis function's contribution to the effective dynamical
    !    matrix from which the potential can be interpolated in the large-supercell
    !    limit.
    module function calculate_dynamical_matrices_BasisFunction(this,qpoints, &
       & thermal_energy,subspaces,subspace_bases,subspace_states,            &
       & subspaces_in_coupling,anharmonic_data) result(output) 
      class(BasisFunction),     intent(in)    :: this
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
  
  interface
    ! Calculate the correction due to double counting
    !    for the interpolated potential.
    module function energy_correction_BasisFunction(this,subspaces, &
       & subspace_bases,subspace_states,anharmonic_data) result(output) 
      class(BasisFunction),     intent(in)    :: this
      type(DegenerateSubspace), intent(in)    :: subspaces(:)
      class(SubspaceBasis),     intent(in)    :: subspace_bases(:)
      class(BasisStates),       intent(inout) :: subspace_states(:)
      type(AnharmonicData),     intent(in)    :: anharmonic_data
      real(dp)                                :: output
    end function
  end interface
  
  interface
    ! Return the monomial terms of this basis function.
    module function terms_BasisFunction(this) result(output) 
      class(BasisFunction), intent(in)   :: this
      type(ComplexMonomial), allocatable :: output(:)
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_BasisFunction(this,input) 
      class(BasisFunction), intent(out) :: this
      type(String),         intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_BasisFunction(this) result(output) 
      class(BasisFunction), intent(in) :: this
      type(String), allocatable        :: output(:)
    end function
  end interface
  
  interface BasisFunction
    module function new_BasisFunction_Strings(input) result(this) 
      type(String), intent(in) :: input(:)
      type(BasisFunction)      :: this
    end function
  
    impure elemental module function new_BasisFunction_StringArray(input) &
       & result(this) 
      type(StringArray), intent(in) :: input
      type(BasisFunction)           :: this
    end function
  end interface
end module
