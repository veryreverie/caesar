! ======================================================================
! A polynomial representation of the stress.
! ======================================================================
module caesar_polynomial_stress_module
  use caesar_common_module
  
  use caesar_anharmonic_common_module
  
  use caesar_polynomial_interpolator_module
  use caesar_stress_basis_function_module
  use caesar_coupling_stress_basis_functions_module
  implicit none
  
  private
  
  public :: PolynomialStress
  
  type, extends(StressData) :: PolynomialStress
    integer,          private :: stress_expansion_order_
    type(RealMatrix), private :: reference_stress_
    type(CouplingStressBasisFunctions), allocatable, private :: &
       & basis_functions_(:)
  contains
    procedure, public, nopass :: representation => &
                               & representation_PolynomialStress
    
    procedure, public :: zero_stress => zero_stress_PolynomialStress
    procedure, public :: add_constant => add_constant_PolynomialStress
    
    procedure, public :: stress_RealModeDisplacement => &
                       & stress_RealModeDisplacement_PolynomialStress
    procedure, public :: stress_ComplexModeDisplacement => &
                       & stress_ComplexModeDisplacement_PolynomialStress
    
    procedure, public :: braket_SubspaceBraKet  => &
                       & braket_SubspaceBraKet_PolynomialStress
    procedure, public :: braket_BasisState  => &
                       & braket_BasisState_PolynomialStress
    procedure, public :: braket_BasisStates => &
                       & braket_BasisStates_PolynomialStress
    
    procedure, public :: simplify => &
                       & simplify_PolynomialStress
    
    procedure, public :: harmonic_expectation => &
                       & harmonic_expectation_PolynomialStress
    
    procedure, public :: can_be_interpolated => &
                       & can_be_interpolated_PolynomialStress
    
    procedure, public :: interpolate_coefficients => &
                       & interpolate_coefficients_PolynomialStress
    
    procedure, public :: calculate_dynamical_matrices => &
                       & calculate_dynamical_matrices_PolynomialStress
    
    procedure, public :: stress_correction => &
                       & stress_correction_PolynomialStress
    
    procedure, public :: expansion_order => expansion_order_PolynomialStress
    
    ! I/O.
    procedure, public :: read  => read_PolynomialStress
    procedure, public :: write => write_PolynomialStress
  end type
  
  interface PolynomialStress
    ! Constructor.
    module function new_PolynomialStress(stress_expansion_order, &
        reference_stress,basis_functions) result(this) 
      integer,                            intent(in) :: stress_expansion_order
      type(RealMatrix),                   intent(in) :: reference_stress
      type(CouplingStressBasisFunctions) ,intent(in) :: basis_functions(:) 
      type(PolynomialStress)                         :: this
    end function
  end interface
  
  interface
    ! Type representation.
    impure elemental module function representation_PolynomialStress() &
       & result(output) 
      type(String) :: output
    end function
  end interface
  
  interface
    ! Set the undisplaced stress to zero.
    impure elemental module subroutine zero_stress_PolynomialStress(this) 
      class(PolynomialStress), intent(inout) :: this
    end subroutine
  end interface
  
  interface
    ! Add a constant to the stress.
    impure elemental module subroutine add_constant_PolynomialStress(this, &
       & input) 
      class(PolynomialStress), intent(inout) :: this
      type(RealMatrix),        intent(in)    :: input
    end subroutine
  end interface
  
  interface
    ! Calculate the stress at a given displacement.
    impure elemental module function stress_RealModeDisplacement_PolynomialStress(this,displacement) result(output) 
      class(PolynomialStress),    intent(in) :: this
      type(RealModeDisplacement), intent(in) :: displacement
      type(RealMatrix)                       :: output
    end function
  end interface
  
  interface
    impure elemental module function stress_ComplexModeDisplacement_PolynomialStress(   this,displacement) result(output) 
      class(PolynomialStress),       intent(in) :: this
      type(ComplexModeDisplacement), intent(in) :: displacement
      type(ComplexMatrix)                       :: output
    end function
  end interface
  
  interface
    ! Integrate the stress between two states.
    impure elemental module subroutine braket_SubspaceBraKet_PolynomialStress(this,braket,whole_subspace,anharmonic_data) 
      class(PolynomialStress), intent(inout)        :: this
      class(SubspaceBraKet),   intent(in)           :: braket
      logical,                 intent(in), optional :: whole_subspace
      type(AnharmonicData),    intent(in)           :: anharmonic_data
    end subroutine
  end interface
  
  interface
    ! Integrate the stress between two states.
    impure elemental module subroutine braket_BasisState_PolynomialStress(this,bra,ket,subspace,subspace_basis,whole_subspace,anharmonic_data) 
      class(PolynomialStress),  intent(inout)        :: this
      class(BasisState),        intent(in)           :: bra
      class(BasisState),        intent(in), optional :: ket
      type(DegenerateSubspace), intent(in)           :: subspace
      class(SubspaceBasis),     intent(in)           :: subspace_basis
      logical,                  intent(in), optional :: whole_subspace
      type(AnharmonicData),     intent(in)           :: anharmonic_data
    end subroutine
  end interface
  
  interface
    impure elemental module subroutine braket_BasisStates_PolynomialStress(this,states,subspace,subspace_basis,whole_subspace,anharmonic_data) 
      class(PolynomialStress),  intent(inout)        :: this
      class(BasisStates),       intent(inout)        :: states
      type(DegenerateSubspace), intent(in)           :: subspace
      class(SubspaceBasis),     intent(in)           :: subspace_basis
      logical,                  intent(in), optional :: whole_subspace
      type(AnharmonicData),     intent(in)           :: anharmonic_data
    end subroutine
  end interface
  
  interface
    ! Identify basis functions which are constant,
    !    add the constant energy to the potential's reference energy,
    !    and remove the term.
    ! Then identify basis functions with the same coupling as a previous coupling,
    !    combine the two and remove the duplicate term.
    module subroutine simplify_PolynomialStress(this) 
      class(PolynomialStress), intent(inout) :: this
    end subroutine
  end interface
  
  interface
    ! Calculate the thermal expectation of the stress, <stress>, for a set of
    !    harmonic states.
    impure elemental module function harmonic_expectation_PolynomialStress(this,frequency,thermal_energy,supercell_size,anharmonic_data) result(output) 
      class(PolynomialStress), intent(in) :: this
      real(dp),                intent(in) :: frequency
      real(dp),                intent(in) :: thermal_energy
      integer,                 intent(in) :: supercell_size
      type(AnharmonicData),    intent(in) :: anharmonic_data
      type(RealMatrix)                    :: output
    end function
  end interface
  
  interface
    ! Calculate the contribution to a given monomial from the interpolation of
    !    this stress.
    ! The result is given as a cartesian tensor.
    module function can_be_interpolated_PolynomialStress(this) result(output) 
      class(PolynomialStress), intent(in) :: this
      logical                             :: output
    end function
  end interface
  
  interface
    ! Constructs the six tensor elements from an array of monomials and a matrix
    !    of coefficients.
    module function generate_stress_elements(monomials,coefficients) &
       & result(output) 
      type(ComplexMonomial), intent(in)      :: monomials(:)
      type(RealMatrix),      intent(in)      :: coefficients
      type(StressBasisFunction), allocatable :: output(:)
    end function
  end interface
  
  interface
    impure elemental module function interpolate_coefficients_PolynomialStress(this,monomial,interpolator) result(output) 
      class(PolynomialStress),      intent(in) :: this
      type(ComplexMonomial),        intent(in) :: monomial
      type(PolynomialInterpolator), intent(in) :: interpolator
      type(ComplexMatrix)                      :: output
    end function
  end interface
  
  interface
    ! Calculate the effective dynamical matrices from which the stress can be
    !    interpolated in the large-supercell limit.
    module function calculate_dynamical_matrices_PolynomialStress(this,   &
       & qpoints,thermal_energy,subspaces,subspace_bases,subspace_states, &
       & anharmonic_data) result(output) 
      class(PolynomialStress),  intent(in)     :: this
      type(QpointData),         intent(in)     :: qpoints(:)
      real(dp),                 intent(in)     :: thermal_energy
      type(DegenerateSubspace), intent(in)     :: subspaces(:)
      class(SubspaceBasis),     intent(in)     :: subspace_bases(:)
      class(BasisStates),       intent(inout)  :: subspace_states(:)
      type(AnharmonicData),     intent(in)     :: anharmonic_data
      type(StressDynamicalMatrix), allocatable :: output(:)
    end function
  end interface
  
  interface
    ! Calculate the correction due to double counting
    !    for the interpolated stress.
    module function stress_correction_PolynomialStress(this,subspaces, &
       & subspace_bases,subspace_states,anharmonic_data) result(output) 
      class(PolynomialStress),  intent(in)    :: this
      type(DegenerateSubspace), intent(in)    :: subspaces(:)
      class(SubspaceBasis),     intent(in)    :: subspace_bases(:)
      class(BasisStates),       intent(inout) :: subspace_states(:)
      type(AnharmonicData),     intent(in)    :: anharmonic_data
      type(RealMatrix)                        :: output
    end function
  end interface
  
  interface
    ! Expansion order.
    impure elemental module function expansion_order_PolynomialStress(this) &
       & result(output) 
      class(PolynomialStress), intent(in) :: this
      integer                             :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_PolynomialStress(this,input) 
      class(PolynomialStress), intent(out) :: this
      type(String),               intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_PolynomialStress(this) result(output) 
      class(PolynomialStress), intent(in) :: this
      type(String), allocatable              :: output(:)
    end function
  end interface
  
  interface PolynomialStress
    module function new_PolynomialStress_Strings(input) result(this) 
      type(String), intent(in) :: input(:)
      type(PolynomialStress)   :: this
    end function
  
    impure elemental module function new_PolynomialStress_StringArray(input) &
       & result(this) 
      type(StringArray), intent(in) :: input
      type(PolynomialStress)        :: this
    end function
  end interface
end module
