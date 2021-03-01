! ======================================================================
! An example module, showing how to use PotentialData and PotentialPointer.
! ======================================================================
module caesar_potential_example_module
  use caesar_common_module
  
  use caesar_anharmonic_common_module
  use caesar_states_module
  
  use caesar_anharmonic_data_module
  implicit none
  
  private
  
  public :: startup_potential_example
  public :: potential_example_subroutine
  public :: PotentialDataExample
  
  type, extends(PotentialData) :: PotentialDataExample
    type(String) :: example_contents
  contains
    procedure, public, nopass :: representation => &
                               & representation_PotentialDataExample
    
    procedure, public :: generate_sampling_points => &
       & generate_sampling_points_PotentialDataExample
    procedure, public :: generate_potential => &
       & generate_potential_PotentialDataExample
    procedure, public :: generate_stress => &
       & generate_stress_PotentialDataExample
    
    procedure, public :: zero_energy => zero_energy_PotentialDataExample
    procedure, public :: add_constant => add_constant_PotentialDataExample
    
    procedure, public :: energy_RealModeDisplacement => &
                       & energy_RealModeDisplacement_PotentialDataExample
    procedure, public :: energy_ComplexModeDisplacement => &
                       & energy_ComplexModeDisplacement_PotentialDataExample
    procedure, public :: force_RealModeDisplacement => &
                       & force_RealModeDisplacement_PotentialDataExample
    procedure, public :: force_ComplexModeDisplacement => &
                       & force_ComplexModeDisplacement_PotentialDataExample
    
    procedure, public :: braket_SubspaceBraKet  => &
                       & braket_SubspaceBraKet_PotentialDataExample
    procedure, public :: braket_BasisState  => &
                       & braket_BasisState_PotentialDataExample
    procedure, public :: braket_BasisStates => &
                       & braket_BasisStates_PotentialDataExample
    
    procedure, public :: harmonic_expectation => &
                       & harmonic_expectation_PotentialDataExample
    
    procedure, public :: coefficients => &
                       & coefficients_PotentialDataExample
    procedure, public :: set_coefficients => &
                       & set_coefficients_PotentialDataExample
    procedure, public :: all_basis_functions => &
                       & all_basis_functions_PotentialDataExample
    procedure, public :: variable_basis_functions => &
                       & variable_basis_functions_PotentialDataExample
    
    procedure, public :: read  => read_PotentialDataExample
    procedure, public :: write => write_PotentialDataExample
  end type
  
  interface
    ! Startup procedure.
    module subroutine startup_potential_example() 
    end subroutine
  end interface
  
  interface PotentialDataExample
    ! --------------------------------------------------
    ! Constructor for example class.
    ! --------------------------------------------------
    ! This is where any PotentialDataExample-specific data is input.
    module function new_PotentialDataExample(example_contents) result(this) 
      type(String), intent(in)   :: example_contents
      type(PotentialDataExample) :: this
    end function
  end interface
  
  interface
    ! Type representation.
    impure elemental module function representation_PotentialDataExample() &
       & result(output) 
      type(String) :: output
    end function
  end interface
  
  interface
    ! --------------------------------------------------
    ! Overloads of PotentialData's methods.
    ! --------------------------------------------------
    module subroutine generate_sampling_points_PotentialDataExample(this, &
       & anharmonic_data,use_forces,energy_to_force_ratio,use_hessians,   &
       & calculate_stress,sampling_points_dir,calculation_writer,logfile) 
      class(PotentialDataExample), intent(inout) :: this
      type(AnharmonicData),        intent(in)    :: anharmonic_data
      logical,                     intent(in)    :: use_forces
      real(dp),                    intent(in)    :: energy_to_force_ratio
      logical,                     intent(in)    :: use_hessians
      logical,                     intent(in)    :: calculate_stress
      type(String),                intent(in)    :: sampling_points_dir
      type(CalculationWriter),     intent(inout) :: calculation_writer
      type(OFile),                 intent(inout) :: logfile
    end subroutine
  end interface
  
  interface
    module subroutine generate_potential_PotentialDataExample(this,       &
       & anharmonic_data,weighted_energy_force_ratio,sampling_points_dir, &
       & calculation_reader,logfile) 
      class(PotentialDataExample), intent(inout) :: this
      type(AnharmonicData),        intent(in)    :: anharmonic_data
      real(dp),                    intent(in)    :: weighted_energy_force_ratio
      type(String),                intent(in)    :: sampling_points_dir
      type(CalculationReader),     intent(inout) :: calculation_reader
      type(OFile),                 intent(inout) :: logfile
    end subroutine
  end interface
  
  interface
    module function generate_stress_PotentialDataExample(this,       &
       & anharmonic_data,sampling_points_dir,stress_expansion_order, &
        stress_subspace_coupling,vscf_basis_functions_only, &
          & calculation_reader,logfile) result(output) 
      class(PotentialDataExample), intent(in)    :: this
      type(AnharmonicData),        intent(in)    :: anharmonic_data
      type(String),                intent(in)    :: sampling_points_dir
      integer,                     intent(in)    :: stress_expansion_order
      type(SubspaceCoupling),      intent(in)    :: stress_subspace_coupling(:)
      logical,intent(in) :: vscf_basis_functions_only
      type(CalculationReader),     intent(inout) :: calculation_reader
      type(OFile),                 intent(inout) :: logfile
      type(StressPointer)                        :: output
    end function
  end interface
  
  interface
    impure elemental module subroutine zero_energy_PotentialDataExample(this) 
      class(PotentialDataExample), intent(inout) :: this
    end subroutine
  end interface
  
  interface
    impure elemental module subroutine add_constant_PotentialDataExample(this,input) 
      class(PotentialDataExample), intent(inout) :: this
      real(dp),                    intent(in)    :: input
    end subroutine
  end interface
  
  interface
    impure elemental module function energy_RealModeDisplacement_PotentialDataExample(   this,displacement) result(output) 
      class(potentialDataExample), intent(in) :: this
      type(RealModeDisplacement),  intent(in) :: displacement
      real(dp)                                :: output
    end function
  end interface
  
  interface
    impure elemental module function energy_ComplexModeDisplacement_PotentialDataExample(  this,displacement) result(output) 
      class(potentialDataExample),    intent(in) :: this
      type(ComplexModeDisplacement),  intent(in) :: displacement
      complex(dp)                                :: output
    end function
  end interface
  
  interface
    impure elemental module function force_RealModeDisplacement_PotentialDataExample(   this,displacement) result(output) 
      class(potentialDataExample), intent(in) :: this
      type(RealModeDisplacement),  intent(in) :: displacement
      type(RealModeForce)                     :: output
    end function
  end interface
  
  interface
    impure elemental module function force_ComplexModeDisplacement_PotentialDataExample(   this,displacement) result(output) 
      class(potentialDataExample),    intent(in) :: this
      type(ComplexModeDisplacement),  intent(in) :: displacement
      type(ComplexModeForce)                     :: output
    end function
  end interface
  
  interface
    impure elemental module subroutine braket_SubspaceBraKet_PotentialDataExample(this,braket,whole_subspace,anharmonic_data) 
      class(PotentialDataExample), intent(inout)        :: this
      class(SubspaceBraKet),       intent(in)           :: braket
      logical,                     intent(in), optional :: whole_subspace
      type(AnharmonicData),        intent(in)           :: anharmonic_data
    end subroutine
  end interface
  
  interface
    impure elemental module subroutine braket_BasisState_PotentialDataExample(this,bra,ket,subspace,subspace_basis,whole_subspace,anharmonic_data) 
      class(PotentialDataExample), intent(inout)        :: this
      class(BasisState),           intent(in)           :: bra
      class(BasisState),           intent(in), optional :: ket
      type(DegenerateSubspace),    intent(in)           :: subspace
      class(SubspaceBasis),        intent(in)           :: subspace_basis
      logical,                     intent(in), optional :: whole_subspace
      type(AnharmonicData),        intent(in)           :: anharmonic_data
    end subroutine
  end interface
  
  interface
    impure elemental module subroutine braket_BasisStates_PotentialDataExample(this,states,subspace,subspace_basis,whole_subspace,anharmonic_data) 
      class(PotentialDataExample), intent(inout)        :: this
      class(BasisStates),          intent(inout)        :: states
      type(DegenerateSubspace),    intent(in)           :: subspace
      class(SubspaceBasis),        intent(in)           :: subspace_basis
      logical,                     intent(in), optional :: whole_subspace
      type(AnharmonicData),        intent(in)           :: anharmonic_data
    end subroutine
  end interface
  
  interface
    impure elemental module function harmonic_expectation_PotentialDataExample(this,frequency,thermal_energy,supercell_size,anharmonic_data) result(output) 
      class(PotentialDataExample), intent(in) :: this
      real(dp),                    intent(in) :: frequency
      real(dp),                    intent(in) :: thermal_energy
      integer,                     intent(in) :: supercell_size
      type(AnharmonicData),        intent(in) :: anharmonic_data
      real(dp)                                :: output
    end function
  end interface
  
  interface
    module function coefficients_PotentialDataExample(this,anharmonic_data) &
       & result(output) 
      class(PotentialDataExample), intent(in) :: this
      type(AnharmonicData),        intent(in) :: anharmonic_data
      real(dp), allocatable                   :: output(:)
    end function
  end interface
  
  interface
    module subroutine set_coefficients_PotentialDataExample(this, &
       & coefficients,anharmonic_data) 
      class(PotentialDataExample), intent(inout) :: this
      real(dp),                    intent(in)    :: coefficients(:)
      type(AnharmonicData),        intent(in)    :: anharmonic_data
    end subroutine
  end interface
  
  interface
    module function all_basis_functions_PotentialDataExample(this, &
       & anharmonic_data) result(output) 
      class(PotentialDataExample), intent(in) :: this
      type(AnharmonicData),        intent(in) :: anharmonic_data
      type(PotentialBasePointer), allocatable :: output(:)
    end function
  end interface
  
  interface
    module function variable_basis_functions_PotentialDataExample(this, &
       & anharmonic_data) result(output) 
      class(PotentialDataExample), intent(in) :: this
      type(AnharmonicData),        intent(in) :: anharmonic_data
      type(PotentialBasePointer), allocatable :: output(:)
    end function
  end interface
  
  interface
    ! --------------------------------------------------
    ! I/O.
    ! --------------------------------------------------
    module subroutine read_PotentialDataExample(this,input) 
      class(PotentialDataExample), intent(out) :: this
      type(String),                intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_PotentialDataExample(this) result(output) 
      class(PotentialDataExample), intent(in) :: this
      type(String), allocatable               :: output(:)
    end function
  end interface
  
  interface PotentialDataExample
    module function new_PotentialDataExample_Strings(input) result(this) 
      type(String), intent(in)   :: input(:)
      type(PotentialDataExample) :: this
    end function
  
    impure elemental module function new_PotentialDataExample_StringArray(input) result(this) 
      type(StringArray), intent(in) :: input
      type(PotentialDataExample)    :: this
    end function
  end interface
  
  interface
    ! --------------------------------------------------
    ! The class in use.
    ! --------------------------------------------------
    module subroutine potential_example_subroutine() 
    end subroutine
  end interface
end module
