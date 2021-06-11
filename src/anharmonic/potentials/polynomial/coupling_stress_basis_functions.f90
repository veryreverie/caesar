! ======================================================================
! A set of stress basis functions spanning a given subspace coupling.
! ======================================================================
module caesar_coupling_stress_basis_functions_module
  use caesar_common_module
  
  use caesar_anharmonic_common_module
  
  use caesar_polynomial_interpolator_module
  use caesar_basis_function_module
  use caesar_stress_basis_function_module
  use caesar_sampling_points_module
  use caesar_sample_result_module
  implicit none
  
  private
  
  public :: CouplingStressBasisFunctions
  public :: size
  public :: generate_stress_basis_functions
  
  type, extends(StressBase) :: CouplingStressBasisFunctions
    type(SubspaceCoupling)                          :: coupling
    type(StressBasisFunction), allocatable, private :: basis_functions_(:)
  contains
    procedure, public, nopass :: representation => &
                               & representation_CouplingStressBasisFunctions
    
    procedure, public :: basis_functions => &
                       & basis_functions_CouplingStressBasisFunctions
    
    procedure, public :: undisplaced_stress => &
                       & undisplaced_stress_CouplingStressBasisFunctions
    
    procedure, public :: zero_stress => &
                       & zero_stress_CouplingStressBasisFunctions
    procedure, public :: add_constant => &
                       & add_constant_CouplingStressBasisFunctions
    
    procedure, public :: stress_RealModeDisplacement => &
       & stress_RealModeDisplacement_CouplingStressBasisFunctions
    procedure, public :: stress_ComplexModeDisplacement => &
       & stress_ComplexModeDisplacement_CouplingStressBasisFunctions
    
    procedure, public :: braket_SubspaceBraKet => &
                       & braket_SubspaceBraKet_CouplingStressBasisFunctions
    procedure, public :: braket_BasisState => &
                       & braket_BasisState_CouplingStressBasisFunctions
    procedure, public :: braket_BasisStates => &
                       & braket_BasisStates_CouplingStressBasisFunctions
    
    procedure, public :: harmonic_expectation => &
                       & harmonic_expectation_CouplingStressBasisFunctions
    
    procedure, public :: append => append_CouplingStressBasisFunctions
    
    procedure, public :: fit_coefficients => &
                       & fit_coefficients_CouplingStressBasisFunctions
    procedure, public :: interpolate_coefficients => &
                       & interpolate_coefficients_CouplingStressBasisFunctions
    procedure, public :: calculate_dynamical_matrices => &
                    & calculate_dynamical_matrices_CouplingStressBasisFunctions
    procedure, public :: stress_correction => &
                       & stress_correction_CouplingStressBasisFunctions
    
    ! I/O.
    procedure, public :: read  => read_CouplingStressBasisFunctions
    procedure, public :: write => write_CouplingStressBasisFunctions
  end type
  
  interface CouplingStressBasisFunctions
    ! Constructor and size function.
    impure elemental module function new_CouplingStressBasisFunctions_empty(coupling) result(this) 
      type(SubspaceCoupling),    intent(in) :: coupling
      type(CouplingStressBasisFunctions)    :: this
    end function
  
    module function new_CouplingStressBasisFunctions(coupling, &
        basis_functions) result(this) 
      type(SubspaceCoupling),    intent(in) :: coupling
      type(StressBasisFunction) ,intent(in) :: basis_functions(:) 
      type(CouplingStressBasisFunctions)    :: this
    end function
  end interface
  
  interface
    impure elemental module function representation_CouplingStressBasisFunctions() result(output) 
      type(String) :: output
    end function
  end interface
  
  interface size
    module function size_CouplingStressBasisFunctions(this) result(output) 
      class(CouplingStressBasisFunctions), intent(in) :: this
      integer                                         :: output
    end function
  end interface
  
  interface
    module function basis_functions_CouplingStressBasisFunctions(this) &
       & result(output) 
      class(CouplingStressBasisFunctions), intent(in) :: this
      type(StressBasisFunction), allocatable          :: output(:)
    end function
  end interface
  
  interface
    impure elemental module function                                             stress_RealModeDisplacement_CouplingStressBasisFunctions(this,displacement) result(output) 
      class(CouplingStressBasisFunctions), intent(in) :: this
      type(RealModeDisplacement),          intent(in) :: displacement
      type(RealMatrix)                                :: output
    end function
  end interface
  
  interface
    impure elemental module function                                                stress_ComplexModeDisplacement_CouplingStressBasisFunctions(this,displacement) result(output) 
      class(CouplingStressBasisFunctions), intent(in) :: this
      type(ComplexModeDisplacement),       intent(in) :: displacement
      type(ComplexMatrix)                             :: output
    end function
  end interface
  
  interface
    impure elemental module subroutine                                            braket_SubspaceBraKet_CouplingStressBasisFunctions(this,braket,whole_subspace,anharmonic_data) 
      class(CouplingStressBasisFunctions), intent(inout)        :: this
      class(SubspaceBraKet),               intent(in)           :: braket
      logical,                             intent(in), optional :: whole_subspace
      type(AnharmonicData),                intent(in)           :: anharmonic_data
    end subroutine
  end interface
  
  interface
    impure elemental module subroutine braket_BasisState_CouplingStressBasisFunctions(   this,bra,ket,subspace,subspace_basis,whole_subspace,anharmonic_data) 
      class(CouplingStressBasisFunctions), intent(inout)        :: this
      class(BasisState),                   intent(in)           :: bra
      class(BasisState),                   intent(in), optional :: ket
      type(DegenerateSubspace),            intent(in)           :: subspace
      class(SubspaceBasis),                intent(in)           :: subspace_basis
      logical,                             intent(in), optional :: whole_subspace
      type(AnharmonicData),                intent(in)           :: anharmonic_data
    end subroutine
  end interface
  
  interface
    impure elemental module subroutine braket_BasisStates_CouplingStressBasisFunctions(   this,states,subspace,subspace_basis,whole_subspace,anharmonic_data) 
      class(CouplingStressBasisFunctions), intent(inout)        :: this
      class(BasisStates),                  intent(inout)        :: states
      type(DegenerateSubspace),            intent(in)           :: subspace
      class(SubspaceBasis),                intent(in)           :: subspace_basis
      logical,                             intent(in), optional :: whole_subspace
      type(AnharmonicData),                intent(in)           :: anharmonic_data
    end subroutine
  end interface
  
  interface
    impure elemental module function harmonic_expectation_CouplingStressBasisFunctions(   this,frequency,thermal_energy,supercell_size,anharmonic_data) result(output) 
      class(CouplingStressBasisFunctions), intent(in) :: this
      real(dp),                            intent(in) :: frequency
      real(dp),                            intent(in) :: thermal_energy
      integer,                             intent(in) :: supercell_size
      type(AnharmonicData),                intent(in) :: anharmonic_data
      type(RealMatrix)                                :: output
    end function
  end interface
  
  interface
    ! Return the stress at zero displacement.
    impure elemental module function undisplaced_stress_CouplingStressBasisFunctions(   this) result(output) 
      class(CouplingStressBasisFunctions), intent(in) :: this
      type(RealMatrix)                                :: output
    end function
  end interface
  
  interface
    impure elemental module subroutine zero_stress_CouplingStressBasisFunctions(this) 
      class(CouplingStressBasisFunctions), intent(inout) :: this
    end subroutine
  end interface
  
  interface
    impure elemental module subroutine add_constant_CouplingStressBasisFunctions(this,input) 
      class(CouplingStressBasisFunctions), intent(inout) :: this
      type(RealMatrix),                    intent(in)    :: input
    end subroutine
  end interface
  
  interface
    ! Append another StressCouplingBasisFunctions to this.
    module subroutine append_CouplingStressBasisFunctions(this,that) 
      class(CouplingStressBasisFunctions), intent(inout) :: this
      class(CouplingStressBasisFunctions), intent(in)    :: that
    end subroutine
  end interface
  
  interface generate_stress_basis_functions
    ! Generate stress basis functions.
    module function generate_stress_basis_functions_SubspaceCoupling(   &
       & coupling,stress_expansion_order,max_subspace_coupling,         &
       & max_qpoint_coupling,structure,complex_modes,qpoints,subspaces, &
       & degenerate_symmetries,vscf_basis_functions_only,               &
       & qpoint_symmetry_groups,subspace_qpoint_stars,logfile) result(output)
      type(SubspaceCoupling),    intent(in)    :: coupling
      integer,                   intent(in)    :: stress_expansion_order
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
      type(CouplingStressBasisFunctions)       :: output
    end function
  end interface
  
  interface
    ! Uses L2 regression to calculate the coefficients of a set of stress basis
    !    functions.
    module subroutine fit_coefficients_CouplingStressBasisFunctions(this, &
       & sampling_points,sample_results,stress) 
      class(CouplingStressBasisFunctions), intent(inout)        :: this
      type(RealModeDisplacement),          intent(in)           :: sampling_points(:)
      type(SampleResult),                  intent(in)           :: sample_results(:)
      class(StressData),                   intent(in), optional :: stress
    end subroutine
  end interface
  
  interface
    ! Calculate the contribution to a given monomial from the interpolation of
    !    this basis function.
    ! The result is given as a cartesian tensor.
    impure elemental module function                                                   interpolate_coefficients_CouplingStressBasisFunctions(this,monomial,interpolator) result(output) 
      class(CouplingStressBasisFunctions), intent(in) :: this
      type(ComplexMonomial),               intent(in) :: monomial
      type(PolynomialInterpolator),        intent(in) :: interpolator
      type(ComplexMatrix)                             :: output
    end function
  end interface
  
  interface
    ! Calculate this basis function's contribution to the effective dynamical
    !    matrix from which the potential can be interpolated in the
    !    large-supercell limit.
    module function calculate_dynamical_matrices_CouplingStressBasisFunctions(this,qpoints,thermal_energy,subspaces,subspace_bases,subspace_states,anharmonic_data) result(output) 
      class(CouplingStressBasisFunctions), intent(in)    :: this
      type(QpointData),                    intent(in)    :: qpoints(:)
      real(dp),                            intent(in)    :: thermal_energy
      type(DegenerateSubspace),            intent(in)    :: subspaces(:)
      class(SubspaceBasis),                intent(in)    :: subspace_bases(:)
      class(BasisStates),                  intent(inout) :: subspace_states(:)
      type(AnharmonicData),                intent(in)    :: anharmonic_data
      type(StressDynamicalMatrix), allocatable           :: output(:)
    end function
  end interface
  
  interface
    ! Calculate the correction due to double counting
    !    for the interpolated stress.
    module function stress_correction_CouplingStressBasisFunctions(this, &
       & subspaces,subspace_bases,subspace_states,anharmonic_data)       &
       & result(output) 
      class(CouplingStressBasisFunctions), intent(in)    :: this
      type(DegenerateSubspace),            intent(in)    :: subspaces(:)
      class(SubspaceBasis),                intent(in)    :: subspace_bases(:)
      class(BasisStates),                  intent(inout) :: subspace_states(:)
      type(AnharmonicData),                intent(in)    :: anharmonic_data
      type(RealMatrix)                                   :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_CouplingStressBasisFunctions(this,input) 
      class(CouplingStressBasisFunctions), intent(out) :: this
      type(String),                        intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_CouplingStressBasisFunctions(this) result(output) 
      class(CouplingStressBasisFunctions), intent(in) :: this
      type(String), allocatable                       :: output(:)
    end function
  end interface
  
  interface CouplingStressBasisFunctions
    module function new_CouplingStressBasisFunctions_Strings(input) &
       & result(this) 
      type(String), intent(in)           :: input(:)
      type(CouplingStressBasisFunctions) :: this
    end function
  end interface
  
  interface CouplingStressBasisFunctions
    impure elemental module function new_CouplingStressBasisFunctions_StringArray(input) result(this) 
      type(StringArray), intent(in)       :: input
      type(CouplingStressBasisFunctions)  :: this
    end function
  end interface
end module
