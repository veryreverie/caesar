! ======================================================================
! A set of basis functions spanning a given subspace coupling.
! ======================================================================
module caesar_coupling_basis_functions_module
  use caesar_common_module
  
  use caesar_anharmonic_common_module
  
  use caesar_polynomial_interpolator_module
  use caesar_basis_function_module
  use caesar_sampling_points_module
  use caesar_sample_result_module
  implicit none
  
  private
  
  public :: CouplingBasisFunctions
  public :: size
  public :: generate_basis_functions
  
  type, extends(PotentialBase) :: CouplingBasisFunctions
    type(SubspaceCoupling)                    :: coupling
    type(BasisFunction), allocatable, private :: basis_functions_(:)
  contains
    procedure, public, nopass :: representation => &
                               & representation_CouplingBasisFunctions
    
    procedure, public :: basis_functions => &
                       & basis_functions_CouplingBasisFunctions
    
    procedure, public :: undisplaced_energy => &
                       & undisplaced_energy_CouplingBasisFunctions
    
    procedure, public :: zero_energy => zero_energy_CouplingBasisFunctions
    procedure, public :: add_constant => add_constant_CouplingBasisFunctions
    
    procedure, public :: energy_RealModeDisplacement => &
                       & energy_RealModeDisplacement_CouplingBasisFunctions
    procedure, public :: energy_ComplexModeDisplacement => &
                       & energy_ComplexModeDisplacement_CouplingBasisFunctions
    
    procedure, public :: force_RealModeDisplacement => &
                       & force_RealModeDisplacement_CouplingBasisFunctions
    procedure, public :: force_ComplexModeDisplacement => &
                       & force_ComplexModeDisplacement_CouplingBasisFunctions
    
    procedure, public :: braket_SubspaceBraKet => &
                       & braket_SubspaceBraKet_CouplingBasisFunctions
    procedure, public :: braket_BasisState => &
                       & braket_BasisState_CouplingBasisFunctions
    procedure, public :: braket_BasisStates => &
                       & braket_BasisStates_CouplingBasisFunctions
    
    procedure, public :: harmonic_expectation => &
                       & harmonic_expectation_CouplingBasisFunctions
    
    procedure, public :: potential_energy_SubspaceBraKet => &
                       & potential_energy_SubspaceBraKet_CouplingBasisFunctions
    procedure, public :: potential_energy_BasisState => &
                       & potential_energy_BasisState_CouplingBasisFunctions
    
    procedure, public :: fit_coefficients => &
                       & fit_coefficients_CouplingBasisFunctions
    procedure, public :: no_coefficients => &
                       & no_coefficients_CouplingBasisFunctions
    procedure, public :: coefficients => &
                       & coefficients_CouplingBasisFunctions
    procedure, public :: set_coefficients => &
                       & set_coefficients_CouplingBasisFunctions
    procedure, public :: zero_coefficients => &
                       & zero_coefficients_CouplingBasisFunctions
    procedure, public :: all_basis_functions => &
                       & all_basis_functions_CouplingBasisFunctions
    procedure, public :: variable_basis_functions => &
                       & variable_basis_functions_CouplingBasisFunctions
    
    procedure, public :: optimise => optimise_CouplingBasisFunctions
    
    procedure, public :: append => append_CouplingBasisFunctions
    
    procedure, public :: interpolate_coefficient => &
                       & interpolate_coefficient_CouplingBasisFunctions
    procedure, public :: add_interpolated_contribution => &
                       & add_interpolated_contribution_CouplingBasisFunctions
    procedure, public :: add_harmonic_contribution => &
                       & add_harmonic_contribution_CouplingBasisFunctions
    procedure, public :: calculate_dynamical_matrices => &
                       & calculate_dynamical_matrices_CouplingBasisFunctions
    procedure, public :: energy_correction => &
                       & energy_correction_CouplingBasisFunctions
    
    ! I/O.
    procedure, public :: read  => read_CouplingBasisFunctions
    procedure, public :: write => write_CouplingBasisFunctions
  end type
  
  interface CouplingBasisFunctions
    impure elemental module function new_CouplingBasisFunctions_empty(coupling) result(this) 
      type(SubspaceCoupling), intent(in) :: coupling
      type(CouplingBasisFunctions)       :: this
    end function
  
    module function new_CouplingBasisFunctions(coupling,basis_functions) &
       & result(this) 
      type(SubspaceCoupling), intent(in) :: coupling
      type(BasisFunction) ,intent(in) :: basis_functions(:) 
      type(CouplingBasisFunctions)       :: this
    end function
  end interface
  
  interface
    impure elemental module function representation_CouplingBasisFunctions() &
       & result(output) 
      type(String) :: output
    end function
  end interface
  
  interface size
    module function size_CouplingBasisFunctions(this) result(output) 
      type(CouplingBasisFunctions), intent(in) :: this
      integer                                  :: output
    end function
  end interface
  
  interface
    module function basis_functions_CouplingBasisFunctions(this) &
       & result(output) 
      class(CouplingBasisFunctions), intent(in) :: this
      type(BasisFunction), allocatable          :: output(:)
    end function
  end interface
  
  interface
    impure elemental module subroutine optimise_CouplingBasisFunctions(this, &
       & subspace,subspace_basis,old_subspace_potential,anharmonic_data) 
      class(CouplingBasisFunctions), intent(inout)        :: this
      type(DegenerateSubspace),      intent(in)           :: subspace
      class(SubspaceBasis),          intent(in)           :: subspace_basis
      class(PotentialData),          intent(in), optional :: old_subspace_potential
      type(AnharmonicData),          intent(in)           :: anharmonic_data
    end subroutine
  end interface
  
  interface
    impure elemental module function energy_RealModeDisplacement_CouplingBasisFunctions(   this,displacement) result(output) 
      class(CouplingBasisFunctions), intent(in) :: this
      type(RealModeDisplacement),    intent(in) :: displacement
      real(dp)                                  :: output
    end function
  end interface
  
  interface
    impure elemental module function                                                       energy_ComplexModeDisplacement_CouplingBasisFunctions(this,displacement) result(output) 
      class(CouplingBasisFunctions), intent(in) :: this
      type(ComplexModeDisplacement), intent(in) :: displacement
      complex(dp)                               :: output
    end function
  end interface
  
  interface
    impure elemental module function force_RealModeDisplacement_CouplingBasisFunctions(   this,displacement) result(output) 
      class(CouplingBasisFunctions), intent(in) :: this
      type(RealModeDisplacement),    intent(in) :: displacement
      type(RealModeForce)                       :: output
    end function
  end interface
  
  interface
    impure elemental module function                                                      force_ComplexModeDisplacement_CouplingBasisFunctions(this,displacement) result(output) 
      class(CouplingBasisFunctions), intent(in) :: this
      type(ComplexModeDisplacement), intent(in) :: displacement
      type(ComplexModeForce)                    :: output
    end function
  end interface
  
  interface
    impure elemental module subroutine braket_SubspaceBraKet_CouplingBasisFunctions(   this,braket,whole_subspace,anharmonic_data) 
      class(CouplingBasisFunctions), intent(inout)        :: this
      class(SubspaceBraKet),         intent(in)           :: braket
      logical,                       intent(in), optional :: whole_subspace
      type(AnharmonicData),          intent(in)           :: anharmonic_data
    end subroutine
  end interface
  
  interface
    impure elemental module subroutine braket_BasisState_CouplingBasisFunctions(this,bra,ket,subspace,subspace_basis,whole_subspace,anharmonic_data) 
      class(CouplingBasisFunctions), intent(inout)        :: this
      class(BasisState),             intent(in)           :: bra
      class(BasisState),             intent(in), optional :: ket
      type(DegenerateSubspace),      intent(in)           :: subspace
      class(SubspaceBasis),          intent(in)           :: subspace_basis
      logical,                       intent(in), optional :: whole_subspace
      type(AnharmonicData),          intent(in)           :: anharmonic_data
    end subroutine
  end interface
  
  interface
    impure elemental module subroutine braket_BasisStates_CouplingBasisFunctions(this,states,subspace,subspace_basis,whole_subspace,anharmonic_data) 
      class(CouplingBasisFunctions), intent(inout)        :: this
      class(BasisStates),            intent(inout)        :: states
      type(DegenerateSubspace),      intent(in)           :: subspace
      class(SubspaceBasis),          intent(in)           :: subspace_basis
      logical,                       intent(in), optional :: whole_subspace
      type(AnharmonicData),          intent(in)           :: anharmonic_data
    end subroutine
  end interface
  
  interface
    impure elemental module function harmonic_expectation_CouplingBasisFunctions(this,frequency,thermal_energy,supercell_size,anharmonic_data) result(output) 
      class(CouplingBasisFunctions), intent(in) :: this
      real(dp),                      intent(in) :: frequency
      real(dp),                      intent(in) :: thermal_energy
      integer,                       intent(in) :: supercell_size
      type(AnharmonicData),          intent(in) :: anharmonic_data
      real(dp)                                  :: output
    end function
  end interface
  
  interface
    module function potential_energy_SubspaceBraKet_CouplingBasisFunctions( &
       & this,braket,anharmonic_data) result(output) 
      class(CouplingBasisFunctions), intent(in) :: this
      class(SubspaceBraKet),         intent(in) :: braket
      type(AnharmonicData),          intent(in) :: anharmonic_data
      real(dp)                                  :: output
    end function
  end interface
  
  interface
    module function potential_energy_BasisState_CouplingBasisFunctions(this, &
       & bra,ket,subspace,subspace_basis,anharmonic_data) result(output) 
      class(CouplingBasisFunctions), intent(in)           :: this
      class(BasisState),             intent(in)           :: bra
      class(BasisState),             intent(in), optional :: ket
      type(DegenerateSubspace),      intent(in)           :: subspace
      class(SubspaceBasis),          intent(in)           :: subspace_basis
      type(AnharmonicData),          intent(in)           :: anharmonic_data
      real(dp)                                            :: output
    end function
  end interface
  
  interface generate_basis_functions
    module function generate_basis_functions_SubspaceCoupling(coupling,    &
       & potential_expansion_order,maximum_coupling_order,structure,       &
       & complex_modes,real_modes,qpoints,subspaces,degenerate_symmetries, &
       & vscf_basis_functions_only,qpoint_symmetry_groups,                 &
       & subspace_qpoint_stars,logfile) result(output) 
      type(SubspaceCoupling),    intent(in)    :: coupling
      integer,                   intent(in)    :: potential_expansion_order
      integer,                   intent(in)    :: maximum_coupling_order
      type(StructureData),       intent(in)    :: structure
      type(ComplexMode),         intent(in)    :: complex_modes(:)
      type(RealMode),            intent(in)    :: real_modes(:)
      type(QpointData),          intent(in)    :: qpoints(:)
      type(DegenerateSubspace),  intent(in)    :: subspaces(:)
      type(DegenerateSymmetry),  intent(in)    :: degenerate_symmetries(:)
      logical,                   intent(in)    :: vscf_basis_functions_only
      type(Group),               intent(in)    :: qpoint_symmetry_groups(:)
      type(SubspaceQpointStars), intent(in)    :: subspace_qpoint_stars(:)
      type(OFile),               intent(inout) :: logfile
      type(CouplingBasisFunctions)             :: output
    end function
  end interface
  
  interface
    ! Return the energy at zero displacement.
    impure elemental module function &
       & undisplaced_energy_CouplingBasisFunctions(this) result(output) 
      class(CouplingBasisFunctions), intent(in) :: this
      real(dp)                                  :: output
    end function
  end interface
  
  interface
    impure elemental module subroutine zero_energy_CouplingBasisFunctions( &
       & this) 
      class(CouplingBasisFunctions), intent(inout) :: this
    end subroutine
  end interface
  
  interface
    impure elemental module subroutine add_constant_CouplingBasisFunctions( &
       & this,input)
      class(CouplingBasisFunctions), intent(inout) :: this
      real(dp),                      intent(in)    :: input
    end subroutine
  end interface
  
  interface
    ! Fit basis function coefficients using L2 regression.
    module subroutine fit_coefficients_CouplingBasisFunctions(this, &
       & sampling_points,sample_results,modes,energy_force_ratio,potential) 
      class(CouplingBasisFunctions), intent(inout)          :: this
      type(RealModeDisplacement),    intent(in)             :: sampling_points(:)
      type(SampleResult),            intent(in)             :: sample_results(:)
      type(RealMode),                intent(in)             :: modes(:)
      real(dp),                      intent(in)             :: energy_force_ratio
      class(PotentialData),          intent(in),   optional :: potential
    end subroutine
  end interface
  
  interface
    ! Get and set basis function coefficients.
    impure elemental module function no_coefficients_CouplingBasisFunctions( &
       & this,maximum_power) result(output) 
      class(CouplingBasisFunctions), intent(in) :: this
      integer,                       intent(in) :: maximum_power
      integer                                   :: output
    end function
  end interface
  
  interface
    module function coefficients_CouplingBasisFunctions(this,maximum_power) &
       & result(output) 
      class(CouplingBasisFunctions), intent(in) :: this
      integer,                       intent(in) :: maximum_power
      real(dp), allocatable                     :: output(:)
    end function
  end interface
  
  interface
    module subroutine set_coefficients_CouplingBasisFunctions(this, &
       & coefficients,maximum_power) 
      class(CouplingBasisFunctions), intent(inout) :: this
      integer,                       intent(in)    :: maximum_power
      real(dp),                      intent(in)    :: coefficients(:)
    end subroutine
  end interface
  
  interface
    module subroutine zero_coefficients_CouplingBasisFunctions(this) 
      class(CouplingBasisFunctions), intent(inout) :: this
    end subroutine
  end interface
  
  interface
    module function all_basis_functions_CouplingBasisFunctions(this) &
       & result(output) 
      class(CouplingBasisFunctions), intent(in) :: this
      type(PotentialBasePointer), allocatable   :: output(:)
    end function
  end interface
  
  interface
    module function variable_basis_functions_CouplingBasisFunctions(this, &
       & maximum_power) result(output) 
      class(CouplingBasisFunctions), intent(in) :: this
      integer,                       intent(in) :: maximum_power
      type(PotentialBasePointer), allocatable   :: output(:)
    end function
  end interface
  
  interface
    ! Append another CouplingBasisFunctions to this.
    module subroutine append_CouplingBasisFunctions(this,that) 
      class(CouplingBasisFunctions), intent(inout) :: this
      type(CouplingBasisFunctions),  intent(in)    :: that
    end subroutine
  end interface
  
  interface
    ! Calculate the contribution to a given monomial from the interpolation of
    !    this basis function.
    impure elemental module function                                   &
       & interpolate_coefficient_CouplingBasisFunctions(this,monomial, &
       & interpolator) result(output) 
      class(CouplingBasisFunctions), intent(in) :: this
      type(ComplexMonomial),         intent(in) :: monomial
      type(PolynomialInterpolator),  intent(in) :: interpolator
      complex(dp)                               :: output
    end function
  end interface
  
  interface
    ! Calculate the contribution to this basis function from another
    !    basis function,
    !    and add this to this basis function's coefficients.
    module subroutine add_interpolated_contribution_CouplingBasisFunctions( &
       & this,basis_function,interpolator) 
      class(CouplingBasisFunctions), intent(inout) :: this
      type(CouplingBasisFunctions),  intent(in)    :: basis_function
      type(PolynomialInterpolator),  intent(in)    :: interpolator
    end subroutine
  end interface
  
  interface
    ! Calculate the contribution to this basis function from a set of harmonic
    !    dynamical matrices.
    module subroutine add_harmonic_contribution_CouplingBasisFunctions(this, &
       & dynamical_matrices,anharmonic_data) 
      class(CouplingBasisFunctions), intent(inout) :: this
      type(DynamicalMatrix),         intent(in)    :: dynamical_matrices(:)
      type(AnharmonicData),          intent(in)    :: anharmonic_data
    end subroutine
  end interface
  
  interface
    ! Calculate this basis function's contribution to the effective dynamical
    !    matrix from which the potential can be interpolated in the
    !    large-supercell limit.
    module function calculate_dynamical_matrices_CouplingBasisFunctions(this, &
       & qpoints,thermal_energy,subspaces,subspace_bases,subspace_states,     &
       & anharmonic_data) result(output) 
      class(CouplingBasisFunctions), intent(in)    :: this
      type(QpointData),              intent(in)    :: qpoints(:)
      real(dp),                      intent(in)    :: thermal_energy
      type(DegenerateSubspace),      intent(in)    :: subspaces(:)
      class(SubspaceBasis),          intent(in)    :: subspace_bases(:)
      class(BasisStates),            intent(inout) :: subspace_states(:)
      type(AnharmonicData),          intent(in)    :: anharmonic_data
      type(DynamicalMatrix), allocatable           :: output(:)
    end function
  end interface
  
  interface
    ! Calculate the correction due to double counting
    !    for the interpolated potential.
    module function energy_correction_CouplingBasisFunctions(this,subspaces, &
       & subspace_bases,subspace_states,anharmonic_data) result(output) 
      class(CouplingBasisFunctions), intent(in)    :: this
      type(DegenerateSubspace),      intent(in)    :: subspaces(:)
      class(SubspaceBasis),          intent(in)    :: subspace_bases(:)
      class(BasisStates),            intent(inout) :: subspace_states(:)
      type(AnharmonicData),          intent(in)    :: anharmonic_data
      real(dp)                                     :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_CouplingBasisFunctions(this,input) 
      class(CouplingBasisFunctions), intent(out) :: this
      type(String),                  intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_CouplingBasisFunctions(this) result(output) 
      class(CouplingBasisFunctions), intent(in) :: this
      type(String), allocatable                 :: output(:)
    end function
  end interface
  
  interface CouplingBasisFunctions
    module function new_CouplingBasisFunctions_Strings(input) result(this) 
      type(String), intent(in)     :: input(:)
      type(CouplingBasisFunctions) :: this
    end function
  
    impure elemental module function new_CouplingBasisFunctions_StringArray( &
      & input) result(this) 
      type(StringArray), intent(in) :: input
      type(CouplingBasisFunctions)  :: this
    end function
  end interface
end module
