! ======================================================================
! A polynomial representation of a potential.
! ======================================================================
module caesar_polynomial_potential_module
  use caesar_common_module
  
  use caesar_states_module
  use caesar_anharmonic_common_module
  
  use caesar_polynomial_interpolator_module
  use caesar_basis_function_module
  use caesar_coupling_basis_functions_module
  use caesar_stress_basis_function_module
  use caesar_coupling_stress_basis_functions_module
  use caesar_vscf_rvectors_module
  use caesar_sampling_points_module
  use caesar_sample_result_module
  use caesar_sample_results_module
  use caesar_polynomial_stress_module
  implicit none
  
  private
  
  public :: startup_polynomial_potential
  
  public :: PolynomialPotential
  
  type, extends(PotentialData) :: PolynomialPotential
    integer,  private :: potential_expansion_order_
    real(dp), private :: reference_energy_
    type(CouplingBasisFunctions), allocatable, private :: basis_functions_(:)
  contains
    procedure, public, nopass :: representation => &
                               & representation_PolynomialPotential
    
    procedure, public :: generate_sampling_points => &
       & generate_sampling_points_PolynomialPotential
    procedure, public :: generate_potential => &
       & generate_potential_PolynomialPotential
    procedure, public :: generate_stress => &
       & generate_stress_PolynomialPotential
    
    procedure, public :: zero_energy => zero_energy_PolynomialPotential
    procedure, public :: add_constant => add_constant_PolynomialPotential
    
    procedure, public :: optimise_subspace_potential => &
                       & optimise_subspace_potential_PolynomialPotential
    
    procedure, public :: energy_RealModeDisplacement => &
                       & energy_RealModeDisplacement_PolynomialPotential
    procedure, public :: energy_ComplexModeDisplacement => &
                       & energy_ComplexModeDisplacement_PolynomialPotential
    procedure, public :: force_RealModeDisplacement => &
                       & force_RealModeDisplacement_PolynomialPotential
    procedure, public :: force_ComplexModeDisplacement => &
                       & force_ComplexModeDisplacement_PolynomialPotential
    
    procedure, public :: braket_SubspaceBraket  => &
                       & braket_SubspaceBraket_PolynomialPotential
    procedure, public :: braket_BasisState  => &
                       & braket_BasisState_PolynomialPotential
    procedure, public :: braket_BasisStates => &
                       & braket_BasisStates_PolynomialPotential
    
    procedure, public :: simplify => &
                       & simplify_PolynomialPotential
    
    procedure, public :: harmonic_expectation => &
                       & harmonic_expectation_PolynomialPotential
    
    procedure, public :: potential_energy_SubspaceBraKet => &
                       & potential_energy_SubspaceBraKet_PolynomialPotential
    procedure, public :: potential_energy_BasisState => &
                       & potential_energy_BasisState_PolynomialPotential
    
    procedure, public :: coefficients => &
                       & coefficients_PolynomialPotential
    procedure, public :: set_coefficients => &
                       & set_coefficients_PolynomialPotential
    procedure, public :: all_basis_functions => &
                       & all_basis_functions_PolynomialPotential
    procedure, public :: variable_basis_functions => &
                       & variable_basis_functions_PolynomialPotential
    
    procedure, public :: can_be_interpolated => &
                       & can_be_interpolated_PolynomialPotential
    procedure, public :: interpolate_potential => &
                       & interpolate_potential_PolynomialPotential
    
    procedure, public :: interpolate_coefficient => &
                       & interpolate_coefficient_PolynomialPotential
    
    procedure, public :: calculate_dynamical_matrices => &
                       & calculate_dynamical_matrices_PolynomialPotential
    
    procedure, public :: energy_correction => &
                       & energy_correction_PolynomialPotential
    
    procedure, public :: expansion_order => expansion_order_PolynomialPotential
    
    ! I/O.
    procedure, public :: read  => read_PolynomialPotential
    procedure, public :: write => write_PolynomialPotential
  end type
  
  interface
    ! Startup procedure.
    module subroutine startup_polynomial_potential() 
    end subroutine
  end interface
  
  interface PolynomialPotential
    ! Constructors.
    module function new_PolynomialPotential(anharmonic_data) result(this) 
      type(AnharmonicData), intent(in) :: anharmonic_data
      type(PolynomialPotential)        :: this
    end function
  
    recursive module function new_PolynomialPotential_PotentialData(input) &
       & result(this) 
      class(PotentialData), intent(in) :: input
      type(PolynomialPotential)        :: this
    end function
  
    module function new_PolynomialPotential_BasisFunctions(potential_expansion_order,reference_energy,basis_functions) result(this) 
      integer,                      intent(in) :: potential_expansion_order
      real(dp),                     intent(in) :: reference_energy
      type(CouplingBasisFunctions) ,intent(in) :: basis_functions(:) 
      type(PolynomialPotential)                :: this
    end function
  end interface
  
  interface
    ! Type representation.
    impure elemental module function representation_PolynomialPotential() &
       & result(output) 
      type(String) :: output
    end function
  end interface
  
  interface
    ! Generate sampling points.
    ! N.B. does not look at use_forces, use_hessians or calculate stress.
    module subroutine generate_sampling_points_PolynomialPotential(this, &
       & anharmonic_data,use_forces,energy_to_force_ratio,use_hessians,  &
       & calculate_stress,sampling_points_dir,calculation_writer,logfile) 
      class(PolynomialPotential), intent(inout) :: this
      type(AnharmonicData),       intent(in)    :: anharmonic_data
      logical,                    intent(in)    :: use_forces
      real(dp),                   intent(in)    :: energy_to_force_ratio
      logical,                    intent(in)    :: use_hessians
      logical,                    intent(in)    :: calculate_stress
      type(String),               intent(in)    :: sampling_points_dir
      type(CalculationWriter),    intent(inout) :: calculation_writer
      type(OFile),                intent(inout) :: logfile
    end subroutine
  end interface
  
  interface
    ! Generate potential.
    module subroutine generate_potential_PolynomialPotential(this,        &
       & anharmonic_data,weighted_energy_force_ratio,sampling_points_dir, &
       & calculation_reader,logfile) 
      class(PolynomialPotential), intent(inout) :: this
      type(AnharmonicData),       intent(in)    :: anharmonic_data
      real(dp),                   intent(in)    :: weighted_energy_force_ratio
      type(String),               intent(in)    :: sampling_points_dir
      type(CalculationReader),    intent(inout) :: calculation_reader
      type(OFile),                intent(inout) :: logfile
    end subroutine
  end interface
  
  interface
    ! Generate stress.
    module function generate_stress_PolynomialPotential(this,        &
       & anharmonic_data,sampling_points_dir,stress_expansion_order, &
        stress_subspace_coupling,vscf_basis_functions_only, &
          & calculation_reader,logfile) result(output) 
      class(PolynomialPotential), intent(in)    :: this
      type(AnharmonicData),       intent(in)    :: anharmonic_data
      type(String),               intent(in)    :: sampling_points_dir
      integer,                    intent(in)    :: stress_expansion_order
      type(SubspaceCoupling),     intent(in)    :: stress_subspace_coupling(:)
      logical,intent(in) :: vscf_basis_functions_only
      type(CalculationReader),    intent(inout) :: calculation_reader
      type(OFile),                intent(inout) :: logfile
      type(StressPointer)                       :: output
    end function
  end interface
  
  interface
    ! Helper functions for generate_potential and generate_stress.
    
    ! Reads sampling points.
    impure elemental module function read_sampling_points(coupling_directory) &
       & result(output) 
      type(String), intent(in) :: coupling_directory
      type(SamplingPoints)     :: output
    end function
  end interface
  
  interface
    ! Generates the sampling point for the equilibrium structure.
    module function generate_equilibrium_sampling_point() result(output) 
      type(RealModeDisplacement) :: output
    end function
  end interface
  
  interface
    ! Reads basis functions.
    impure elemental module function read_basis_functions(coupling_directory) &
       & result(output) 
      type(String), intent(in)     :: coupling_directory
      type(CouplingBasisFunctions) :: output
    end function
  end interface
  
  interface
    ! Generates the basis function which is just a constant.
    module function generate_constant_basis_function() result(output) 
      type(BasisFunction) :: output
    end function
  end interface
  
  interface
    ! Reads sample results.
    impure elemental module function read_sample_results(sampling_points, &
       & coupling_directory,anharmonic_data,calculation_reader) result(output) 
      type(SamplingPoints),    intent(in)    :: sampling_points
      type(String),            intent(in)    :: coupling_directory
      type(AnharmonicData),    intent(in)    :: anharmonic_data
      type(CalculationReader), intent(inout) :: calculation_reader
      type(SampleResults)                    :: output
    end function
  end interface
  
  interface
    ! Reads the sample result corresponding to equilibrium.
    module function read_equilibrium_sample_result(equilibrium_sampling_point,sampling_points_directory,anharmonic_data,calculation_reader) result(output) 
      type(RealModeDisplacement), intent(in)    :: equilibrium_sampling_point
      type(String),               intent(in)    :: sampling_points_directory
      type(AnharmonicData),       intent(in)    :: anharmonic_data
      type(CalculationReader),    intent(inout) :: calculation_reader
      type(SampleResult)                        :: output
    end function
  end interface
  
  interface
    ! Set the undisplaced energy to zero.
    impure elemental module subroutine zero_energy_PolynomialPotential(this) 
      class(PolynomialPotential), intent(inout) :: this
    end subroutine
  end interface
  
  interface
    ! Add a constant to the energy.
    impure elemental module subroutine add_constant_PolynomialPotential(this,input) 
      class(PolynomialPotential), intent(inout) :: this
      real(dp),                   intent(in)    :: input
    end subroutine
  end interface
  
  interface
    ! Finalise a subspace potential.
    ! Re-arranges basis functions to remove duplicates and separate all monomials.
    module subroutine optimise_subspace_potential_PolynomialPotential(this, &
       & subspace,subspace_basis,old_subspace_potential,anharmonic_data) 
      class(PolynomialPotential), intent(inout)        :: this
      type(DegenerateSubspace),   intent(in)           :: subspace
      class(SubspaceBasis),       intent(in)           :: subspace_basis
      class(PotentialData),       intent(in), optional :: old_subspace_potential
      type(AnharmonicData),       intent(in)           :: anharmonic_data
    end subroutine
  end interface
  
  interface
    ! Calculate the energy at a given displacement.
    impure elemental module function energy_RealModeDisplacement_PolynomialPotential(   this,displacement) result(output) 
      class(PolynomialPotential), intent(in) :: this
      type(RealModeDisplacement), intent(in) :: displacement
      real(dp)                               :: output
    end function
  end interface
  
  interface
    impure elemental module function energy_ComplexModeDisplacement_PolynomialPotential(   this,displacement) result(output) 
      class(PolynomialPotential),    intent(in) :: this
      type(ComplexModeDisplacement), intent(in) :: displacement
      complex(dp)                               :: output
    end function
  end interface
  
  interface
    ! Calculate the force at a given displacement.
    impure elemental module function force_RealModeDisplacement_PolynomialPotential(   this,displacement) result(output) 
      class(PolynomialPotential), intent(in) :: this
      type(RealModeDisplacement), intent(in) :: displacement
      type(RealModeForce)                    :: output
    end function
  end interface
  
  interface
    impure elemental module function force_ComplexModeDisplacement_PolynomialPotential(   this,displacement) result(output) 
      class(PolynomialPotential),    intent(in) :: this
      type(ComplexModeDisplacement), intent(in) :: displacement
      type(ComplexModeForce)                    :: output
    end function
  end interface
  
  interface
    ! Integrate the potential between two states.
    impure elemental module subroutine braket_SubspaceBraKet_PolynomialPotential(this,braket,whole_subspace,anharmonic_data) 
      class(PolynomialPotential), intent(inout)        :: this
      class(SubspaceBraKet),      intent(in)           :: braket
      logical,                    intent(in), optional :: whole_subspace
      type(AnharmonicData),       intent(in)           :: anharmonic_data
    end subroutine
  end interface
  
  interface
    ! Integrate the potential between two states.
    impure elemental module subroutine braket_BasisState_PolynomialPotential(this,bra,ket,subspace,subspace_basis,whole_subspace,anharmonic_data) 
      class(PolynomialPotential), intent(inout)        :: this
      class(BasisState),          intent(in)           :: bra
      class(BasisState),          intent(in), optional :: ket
      type(DegenerateSubspace),   intent(in)           :: subspace
      class(SubspaceBasis),       intent(in)           :: subspace_basis
      logical,                    intent(in), optional :: whole_subspace
      type(AnharmonicData),       intent(in)           :: anharmonic_data
    end subroutine
  end interface
  
  interface
    impure elemental module subroutine braket_BasisStates_PolynomialPotential(this,states,subspace,subspace_basis,whole_subspace,anharmonic_data) 
      class(PolynomialPotential), intent(inout)        :: this
      class(BasisStates),         intent(inout)        :: states
      type(DegenerateSubspace),   intent(in)           :: subspace
      class(SubspaceBasis),       intent(in)           :: subspace_basis
      logical,                    intent(in), optional :: whole_subspace
      type(AnharmonicData),       intent(in)           :: anharmonic_data
    end subroutine
  end interface
  
  interface
    ! Identify basis functions which are constant,
    !    add the constant energy to the potential's reference energy,
    !    and remove the term.
    ! Then identify basis functions with the same coupling as a previous coupling,
    !    combine the two and remove the duplicate term.
    module subroutine simplify_PolynomialPotential(this) 
      class(PolynomialPotential), intent(inout) :: this
    end subroutine
  end interface
  
  interface
    ! Calculate the thermal expectation of the potential, <V>, for a set of
    !    harmonic states.
    impure elemental module function harmonic_expectation_PolynomialPotential(this,frequency,thermal_energy,supercell_size,anharmonic_data) result(output) 
      class(PolynomialPotential), intent(in) :: this
      real(dp),                   intent(in) :: frequency
      real(dp),                   intent(in) :: thermal_energy
      integer,                    intent(in) :: supercell_size
      type(AnharmonicData),       intent(in) :: anharmonic_data
      real(dp)                               :: output
    end function
  end interface
  
  interface
    recursive module function potential_energy_SubspaceBraKet_PolynomialPotential(this,braket,anharmonic_data) result(output) 
      class(PolynomialPotential), intent(in) :: this
      class(SubspaceBraKet),      intent(in) :: braket
      type(AnharmonicData),       intent(in) :: anharmonic_data
      real(dp)                               :: output
    end function
  end interface
  
  interface
    recursive module function potential_energy_BasisState_PolynomialPotential(this,bra,ket,subspace,subspace_basis,anharmonic_data) result(output) 
      class(PolynomialPotential), intent(in)           :: this
      class(BasisState),          intent(in)           :: bra
      class(BasisState),          intent(in), optional :: ket
      type(DegenerateSubspace),   intent(in)           :: subspace
      class(SubspaceBasis),       intent(in)           :: subspace_basis
      type(AnharmonicData),       intent(in)           :: anharmonic_data
      real(dp)                                         :: output
    end function
  end interface
  
  interface
    module function coefficients_PolynomialPotential(this,anharmonic_data) &
       & result(output) 
      class(PolynomialPotential), intent(in) :: this
      type(AnharmonicData),       intent(in) :: anharmonic_data
      real(dp), allocatable                  :: output(:)
    end function
  end interface
  
  interface
    module subroutine set_coefficients_PolynomialPotential(this, &
       & coefficients,anharmonic_data) 
      class(PolynomialPotential), intent(inout) :: this
      real(dp),                   intent(in)    :: coefficients(:)
      type(AnharmonicData),       intent(in)    :: anharmonic_data
    end subroutine
  end interface
  
  interface
    module function all_basis_functions_PolynomialPotential(this, &
       & anharmonic_data) result(output) 
      class(PolynomialPotential), intent(in)  :: this
      type(AnharmonicData),       intent(in)  :: anharmonic_data
      type(PotentialBasePointer), allocatable :: output(:)
    end function
  end interface
  
  interface
    module function variable_basis_functions_PolynomialPotential(this, &
       & anharmonic_data) result(output) 
      class(PolynomialPotential), intent(in)  :: this
      type(AnharmonicData),       intent(in)  :: anharmonic_data
      type(PotentialBasePointer), allocatable :: output(:)
    end function
  end interface
  
  interface
    ! Interpolate the potential.
    module function can_be_interpolated_PolynomialPotential(this) &
       & result(output) 
      class(PolynomialPotential), intent(in) :: this
      logical                                :: output
    end function
  end interface
  
  interface
    module subroutine interpolate_potential_PolynomialPotential(this, &
       & anharmonic_min_images,potential,anharmonic_data,             &
       & interpolated_anharmonic_data,difference_dynamical_matrices,logfile) 
      class(PolynomialPotential), intent(inout) :: this
      type(MinImages),            intent(in)    :: anharmonic_min_images(:,:)
      class(PotentialData),       intent(in)    :: potential
      type(AnharmonicData),       intent(in)    :: anharmonic_data
      type(AnharmonicData),       intent(in)    :: interpolated_anharmonic_data
      type(DynamicalMatrix),      intent(in)    :: difference_dynamical_matrices(:)
      type(OFile),                intent(inout) :: logfile
    end subroutine
  end interface
  
  interface
    ! Calculate the contribution to a given monomial from the interpolation of
    !    this potential.
    impure elemental module function interpolate_coefficient_PolynomialPotential(this,monomial,interpolator) result(output) 
      class(PolynomialPotential),   intent(in) :: this
      type(ComplexMonomial),        intent(in) :: monomial
      type(PolynomialInterpolator), intent(in) :: interpolator
      complex(dp)                              :: output
    end function
  end interface
  
  interface
    ! Calculate the effective dynamical matrices from which the potential can be
    !    interpolated in the large-supercell limit.
    module function calculate_dynamical_matrices_PolynomialPotential(this, &
       & qpoints,thermal_energy,subspaces,subspace_bases,subspace_states,  &
       & anharmonic_data) result(output) 
      class(PolynomialPotential), intent(in)    :: this
      type(QpointData),           intent(in)    :: qpoints(:)
      real(dp),                   intent(in)    :: thermal_energy
      type(DegenerateSubspace),   intent(in)    :: subspaces(:)
      class(SubspaceBasis),       intent(in)    :: subspace_bases(:)
      class(BasisStates),         intent(inout) :: subspace_states(:)
      type(AnharmonicData),       intent(in)    :: anharmonic_data
      type(DynamicalMatrix), allocatable        :: output(:)
    end function
  end interface
  
  interface
    ! Calculate the correction due to double counting
    !    for the interpolated potential.
    module function energy_correction_PolynomialPotential(this,subspaces, &
       & subspace_bases,subspace_states,anharmonic_data) result(output) 
      class(PolynomialPotential), intent(in)    :: this
      type(DegenerateSubspace),   intent(in)    :: subspaces(:)
      class(SubspaceBasis),       intent(in)    :: subspace_bases(:)
      class(BasisStates),         intent(inout) :: subspace_states(:)
      type(AnharmonicData),       intent(in)    :: anharmonic_data
      real(dp)                                  :: output
    end function
  end interface
  
  interface
    ! Expansion order.
    impure elemental module function expansion_order_PolynomialPotential(this) result(output) 
      class(PolynomialPotential), intent(in) :: this
      integer                                :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_PolynomialPotential(this,input) 
      class(PolynomialPotential), intent(out) :: this
      type(String),               intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_PolynomialPotential(this) result(output) 
      class(PolynomialPotential), intent(in) :: this
      type(String), allocatable              :: output(:)
    end function
  end interface
  
  interface PolynomialPotential
    module function new_PolynomialPotential_Strings(input) result(this) 
      type(String), intent(in)  :: input(:)
      type(PolynomialPotential) :: this
    end function
  
    impure elemental module function new_PolynomialPotential_StringArray(input) result(this) 
      type(StringArray), intent(in) :: input
      type(PolynomialPotential)     :: this
    end function
  end interface
end module
