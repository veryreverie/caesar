! ======================================================================
! Four abstract classes, and pointers to those classes.
!    - SubspaceBasis defines the basis of states spanning a subspace.
!    - PotentialData defines a potential energy surface.
!    - StressData defines a stress surface.
! ======================================================================
! The intent is that SubspaceBasis is mostly unchanging, and specific sets
!    of states are defined in terms of this basis.
! The pointer types exist because the syntax for polymorphic objects
!    (those defined as class(X) rather than type(X)) is unpleasant,
!    and because compilers do not seem to be able to handle these objects
!    consistently.
! An object of type(XPointer) is effectively an object of class(X),
!    but with better syntax and compiler support.
! See the example module in potential_example.f90 for how to use this module,
!    using the PotentialData and PotentialPointer types as an example.
!
! N.B. PotentialData objects should return energies normalised per primitive
!    cell.
module abstract_classes_module
  use common_module
  
  use subspace_coupling_module
  use anharmonic_data_module
  use subspace_wavefunctions_module
  use stress_prefactors_module
  use subspace_state_module
  use basis_state_module
  use basis_states_module
  use sparse_monomial_module
  implicit none
  
  private
  
  public :: SubspaceBasis
  public :: SubspaceBasisPointer
  
  public :: PotentialData
  public :: PotentialPointer
  
  public :: StressData
  public :: StressPointer
  
  ! ----------------------------------------------------------------------
  ! Abstract and pointer type definitions.
  ! ----------------------------------------------------------------------
  ! N.B. all objects are defined before any type-bound procedures are defined,
  !    to allow the type-bound procedures of all types to use all of the
  !    other types.
  
  type, abstract, extends(Stringsable) :: SubspaceBasis
  contains
    procedure(representation_SubspaceBasis), public, deferred, nopass :: &
       & representation
    procedure, public :: startup => startup_SubspaceBasis
    
    ! Functionality for generating SubspaceStates.
    procedure(initial_states_SubspaceBasis), public, deferred :: initial_states
    procedure(calculate_states_SubspaceBasis), public, deferred :: &
       & calculate_states
    
    procedure, public :: process_subspace_potential => &
                       & process_subspace_potential_SubspaceBasis
    procedure, public :: process_subspace_stress => &
                       & process_subspace_stress_SubspaceBasis
    
    procedure(mode_ids_SubspaceBasis), public, deferred :: mode_ids
    procedure(paired_mode_ids_SubspaceBasis), public, deferred :: &
       & paired_mode_ids
    
    ! If ket is not given, <this|this>, otherwise <this|ket>.
    procedure(inner_product_SubspaceBasis), public, deferred :: inner_product
    
    ! Integrals of the form <i|V|j>
    ! Either <this|V|this> or <this|V|ket>, where V is a ComplexMonomial.
    generic, public :: integrate =>          &
                     & integrate_BasisState, &
                     & integrate_BasisStates
    procedure(integrate_BasisState_SubspaceBasis), public, deferred :: &
       & integrate_BasisState
    procedure(integrate_BasisStates_SubspaceBasis), public, deferred :: &
       & integrate_BasisStates
    
    ! Either <this|T|this> or <this|T|ket>, where T is the kinetic energy.
    procedure(kinetic_energy_SubspaceBasis), public, deferred :: &
       & kinetic_energy
    
    ! Either <this|V|this> or <this|V|ket>, where V is the harmonic potential
    !    energy.
    procedure(harmonic_potential_energy_SubspaceBasis), public, deferred :: &
       & harmonic_potential_energy
    
    ! Either <this|stress|this> or <this|stress|ket>, where stress is the
    !    kinetic stress.
    procedure(kinetic_stress_SubspaceBasis), public, deferred :: &
       & kinetic_stress
    
    ! Functionality involving SubspaceStates.
    procedure(thermodynamic_data_SubspaceBasis), public, deferred :: &
       & thermodynamic_data
    procedure(wavefunctions_SubspaceBasis),      public, deferred :: &
       & wavefunctions
  end type
  
  type, extends(SubspaceBasis) :: SubspaceBasisPointer
    type(String),                      private :: representation_
    class(SubspaceBasis), allocatable, private :: basis_
  contains
    procedure, private :: check => check_SubspaceBasisPointer
    
    procedure, public, nopass :: representation => &
                               & representation_SubspaceBasisPointer
    
    procedure, public :: basis => basis_SubspaceBasisPointer
    
    procedure, public :: initial_states => initial_states_SubspaceBasisPointer
    procedure, public :: calculate_states => &
                       & calculate_states_SubspaceBasisPointer
    
    procedure, public :: process_subspace_potential => &
                       & process_subspace_potential_SubspaceBasisPointer
    procedure, public :: process_subspace_stress => &
                       & process_subspace_stress_SubspaceBasisPointer
    
    procedure, public :: mode_ids => mode_ids_SubspaceBasisPointer
    procedure, public :: paired_mode_ids => &
                       & paired_mode_ids_SubspaceBasisPointer
    
    procedure, public :: inner_product => &
                       & inner_product_SubspaceBasisPointer
    procedure, public :: integrate_BasisState => &
                       & integrate_BasisState_SubspaceBasisPointer
    procedure, public :: kinetic_energy => &
                       & kinetic_energy_SubspaceBasisPointer
    procedure, public :: harmonic_potential_energy => &
                       & harmonic_potential_energy_SubspaceBasisPointer
    procedure, public :: kinetic_stress => &
                       & kinetic_stress_SubspaceBasisPointer
    
    procedure, public :: thermodynamic_data => &
                       & thermodynamic_data_SubspaceBasisPointer
    procedure, public :: wavefunctions => wavefunctions_SubspaceBasisPointer
    procedure, public :: integrate_BasisStates => &
                       & integrate_BasisStates_SubspaceBasisPointer
    
    ! I/O.
    procedure, public :: read  => read_SubspaceBasisPointer
    procedure, public :: write => write_SubspaceBasisPointer
  end type
  
  ! An array of all types which extend SubspaceBasis.
  ! This array will be filled in by startup routines.
  type(SubspaceBasisPointer), allocatable :: TYPES_SubspaceBasis(:)
  
  type, abstract, extends(Stringsable) :: PotentialData
  contains
    procedure(representation_PotentialData), public, deferred, nopass :: &
       & representation
    procedure, public :: startup => startup_PotentialData
    
    ! Given all input data, generate the set of sampling points at which
    !    electronic structure calculations should be run in order to
    !    map out the potential.
    procedure(generate_sampling_points_PotentialData), public, deferred :: &
       & generate_sampling_points
    
    ! Given the same data as above, generate the potential.
    ! All electronic structure calculations generated by
    !    generate_sampling_points will all have been run.
    procedure(generate_potential_PotentialData), public, deferred :: &
       & generate_potential
    
    ! Given the same data as above, and assuming that generate_potential
    !    has already been called, generate a StressData.
    procedure(generate_stress_PotentialData), public, deferred :: &
       & generate_stress
    
    ! Return the energy at zero displacement, or set this energy to zero,
    !    or add a constant to this energy.
    procedure, public :: undisplaced_energy
    procedure(zero_energy_PotentialData), public, deferred :: zero_energy
    procedure(add_constant_PotentialData), public, deferred :: add_constant
    
    ! Finalise the subspace potential.
    procedure, public :: finalise_subspace_potential => &
                       & finalise_subspace_potential_PotentialData
    
    ! Return the energy and force at a given real or complex displacement.
    generic, public :: energy =>                    &
                     & energy_RealModeDisplacement, &
                     & energy_ComplexModeDisplacement
    generic, public :: force  =>                   &
                     & force_RealModeDisplacement, &
                     & force_ComplexModeDisplacement
    procedure(energy_RealModeDisplacement_PotentialData), public, &
       & deferred :: energy_RealModeDisplacement
    procedure(energy_ComplexModeDisplacement_PotentialData), public, &
       & deferred :: energy_ComplexModeDisplacement
    procedure(force_RealModeDisplacement_PotentialData), public, &
       & deferred :: force_RealModeDisplacement
    procedure(force_ComplexModeDisplacement_PotentialData), public, &
       & deferred :: force_ComplexModeDisplacement
    
    ! Evaluate <bra|potential|ket>.
    generic, public :: braket => braket_SubspaceState, &
                               & braket_BasisState,    &
                               & braket_BasisStates
    procedure(braket_SubspaceState_PotentialData), public, deferred :: &
       & braket_SubspaceState
    procedure(braket_BasisState_PotentialData), public, deferred :: &
       & braket_BasisState
    procedure(braket_BasisStates_PotentialData), public, deferred :: &
       & braket_BasisStates
    
    ! Evaluate the thermal expectation of the potential for a set of harmonic
    !    states.
    procedure(harmonic_expectation_PotentialData), public, deferred :: &
       & harmonic_expectation
    
    ! Convert the potential to and from a set of real coefficients.
    ! Required for Pulay scheme.
    procedure(coefficients_PotentialData), public, deferred :: &
       & coefficients
    procedure(set_coefficients_PotentialData), public, deferred :: &
       & set_coefficients
  end type
  
  type, extends(PotentialData) :: PotentialPointer
    type(String),                      private :: representation_
    class(PotentialData), allocatable, private :: potential_
  contains
    procedure, private :: check => check_PotentialPointer
    
    procedure, public, nopass :: representation => &
                               & representation_PotentialPointer
    
    procedure, public :: potential => potential_PotentialPointer
    
    procedure, public :: generate_sampling_points => &
       & generate_sampling_points_PotentialPointer
    procedure, public :: generate_potential => &
       & generate_potential_PotentialPointer
    procedure, public :: generate_stress => &
       & generate_stress_PotentialPointer
    
    procedure, public :: zero_energy => zero_energy_PotentialPointer
    procedure, public :: add_constant => add_constant_PotentialPointer
    
    procedure, public :: finalise_subspace_potential => &
                       & finalise_subspace_potential_PotentialPointer
    
    procedure, public :: energy_RealModeDisplacement => &
                       & energy_RealModeDisplacement_PotentialPointer
    procedure, public :: energy_ComplexModeDisplacement => &
                       & energy_ComplexModeDisplacement_PotentialPointer
    procedure, public :: force_RealModeDisplacement => &
                       & force_RealModeDisplacement_PotentialPointer
    procedure, public :: force_ComplexModeDisplacement => &
                       & force_ComplexModeDisplacement_PotentialPointer
    
    procedure, public :: braket_SubspaceState  => &
                       & braket_SubspaceState_PotentialPointer
    procedure, public :: braket_BasisState  => &
                       & braket_BasisState_PotentialPointer
    procedure, public :: braket_BasisStates => &
                       & braket_BasisStates_PotentialPointer
    
    procedure, public :: harmonic_expectation => &
                       & harmonic_expectation_PotentialPointer
    
    procedure, public :: coefficients => &
                       & coefficients_PotentialPointer
    procedure, public :: set_coefficients => &
                       & set_coefficients_PotentialPointer
    
    ! I/O.
    procedure, public :: read  => read_PotentialPointer
    procedure, public :: write => write_PotentialPointer
  end type
  
  ! An array of all types which extend PotentialData.
  ! This array will be filled in by startup routines.
  type(PotentialPointer), allocatable :: TYPES_PotentialData(:)
  
  type, abstract, extends(Stringsable) :: StressData
  contains
    procedure(representation_StressData), public, deferred, nopass :: &
       & representation
    procedure, public :: startup => startup_StressData
    
    ! Return the stress at zero displacement, or set this stress to zero.
    procedure, public :: undisplaced_stress
    procedure(zero_stress_StressData), public, deferred :: zero_stress
    procedure(add_constant_StressData), public, deferred :: add_constant
    
    ! Return the stress at a given real or complex displacement.
    generic, public :: stress =>                    &
                     & stress_RealModeDisplacement, &
                     & stress_ComplexModeDisplacement
    procedure(stress_RealModeDisplacement_StressData), public, &
       & deferred :: stress_RealModeDisplacement
    procedure(stress_ComplexModeDisplacement_StressData), public, &
       & deferred :: stress_ComplexModeDisplacement
    
    ! Evaluate <bra|stress|ket>.
    generic, public :: braket => braket_SubspaceState, &
                               & braket_BasisState, &
                               & braket_BasisStates
    procedure(braket_SubspaceState_StressData),  public, deferred :: &
       & braket_SubspaceState
    procedure(braket_BasisState_StressData),  public, deferred :: &
       & braket_BasisState
    procedure(braket_BasisStates_StressData), public, deferred :: &
       & braket_BasisStates
    
    ! Evaluate the thermal expectation of the stress for a set of harmonic
    !    states.
    procedure(harmonic_expectation_StressData), public, deferred :: &
       & harmonic_expectation
  end type
  
  type, extends(StressData) :: StressPointer
    type(String),                   private :: representation_
    class(StressData), allocatable, private :: stress_
  contains
    procedure, private :: check => check_StressPointer
    
    procedure, public, nopass :: representation => &
                               & representation_StressPointer
    
    procedure, public :: zero_stress => zero_stress_StressPointer
    procedure, public :: add_constant => add_constant_StressPointer
    
    procedure, public :: stress_RealModeDisplacement => &
                       & stress_RealModeDisplacement_StressPointer
    procedure, public :: stress_ComplexModeDisplacement => &
                       & stress_ComplexModeDisplacement_StressPointer
    
    procedure, public :: braket_SubspaceState  => &
                       & braket_SubspaceState_StressPointer
    procedure, public :: braket_BasisState  => &
                       & braket_BasisState_StressPointer
    procedure, public :: braket_BasisStates => &
                       & braket_BasisStates_StressPointer
    
    procedure, public :: harmonic_expectation => &
                       & harmonic_expectation_StressPointer
    
    ! I/O.
    procedure, public :: read  => read_StressPointer
    procedure, public :: write => write_StressPointer
  end type
  
  ! An array of all types which extend StressData.
  ! This array will be filled in by startup routines.
  type(StressPointer), allocatable :: TYPES_StressData(:)
  
  ! ----------------------------------------------------------------------
  ! Abstract interfaces, defining the functionality of the abstract classes.
  ! ----------------------------------------------------------------------
  abstract interface
    ! SubspaceBasis procedures.
    impure elemental function representation_SubspaceBasis() result(output)
      import String
      implicit none
      
      type(String) :: output
    end function
    
    impure elemental function initial_states_SubspaceBasis(this,subspace, &
       & anharmonic_data) result(output)
      import SubspaceBasis
      import DegenerateSubspace
      import AnharmonicData
      import BasisStatesPointer
      implicit none
      
      class(SubspaceBasis),     intent(in) :: this
      type(DegenerateSubspace), intent(in) :: subspace
      type(AnharmonicData),     intent(in) :: anharmonic_data
      type(BasisStatesPointer)             :: output
    end function
    
    impure elemental function calculate_states_SubspaceBasis(this,subspace,   &
       & subspace_potential,thermal_energy,energy_convergence,                &
       & no_converged_calculations,max_pulay_iterations,pre_pulay_iterations, &
       & pre_pulay_damping,anharmonic_data) result(output)
      import SubspaceBasis
      import DegenerateSubspace
      import PotentialData
      import dp
      import AnharmonicData
      import BasisStatesPointer
      implicit none
      
      class(SubspaceBasis),     intent(in) :: this
      type(DegenerateSubspace), intent(in) :: subspace
      class(PotentialData),     intent(in) :: subspace_potential
      real(dp),                 intent(in) :: thermal_energy
      real(dp),                 intent(in) :: energy_convergence
      integer,                  intent(in) :: no_converged_calculations
      integer,                  intent(in) :: max_pulay_iterations
      integer,                  intent(in) :: pre_pulay_iterations
      real(dp),                 intent(in) :: pre_pulay_damping
      type(AnharmonicData),     intent(in) :: anharmonic_data
      type(BasisStatesPointer)             :: output
    end function
    
    function mode_ids_SubspaceBasis(this,subspace,anharmonic_data) &
       & result(output)
      import SubspaceBasis
      import DegenerateSubspace
      import AnharmonicData
      implicit none
      
      class(SubspaceBasis),     intent(in) :: this
      type(DegenerateSubspace), intent(in) :: subspace
      type(AnharmonicData),     intent(in) :: anharmonic_data
      integer, allocatable                 :: output(:)
    end function
    
    function paired_mode_ids_SubspaceBasis(this,subspace,anharmonic_data) &
       & result(output)
      import SubspaceBasis
      import DegenerateSubspace
      import AnharmonicData
      implicit none
      
      class(SubspaceBasis),     intent(in) :: this
      type(DegenerateSubspace), intent(in) :: subspace
      type(AnharmonicData),     intent(in) :: anharmonic_data
      integer, allocatable                 :: output(:)
    end function
    
    impure elemental function inner_product_SubspaceBasis(this,bra,ket, &
       & subspace,anharmonic_data) result(output)
      import SubspaceBasis
      import BasisState
      import DegenerateSubspace
      import AnharmonicData
      import dp
      implicit none
      
      class(SubspaceBasis),     intent(in)           :: this
      class(BasisState),        intent(in)           :: bra
      class(BasisState),        intent(in), optional :: ket
      type(DegenerateSubspace), intent(in)           :: subspace
      type(AnharmonicData),     intent(in)           :: anharmonic_data
      real(dp)                                       :: output
    end function
    
    impure elemental function integrate_BasisState_SubspaceBasis(this, &
       & bra,monomial,ket,subspace,anharmonic_data) result(output)
      import SubspaceBasis
      import BasisState
      import SparseMonomial
      import DegenerateSubspace
      import AnharmonicData
      import dp
      implicit none
      
      class(SubspaceBasis),     intent(in)           :: this
      class(BasisState),        intent(in)           :: bra
      type(SparseMonomial),     intent(in)           :: monomial
      class(BasisState),        intent(in), optional :: ket
      type(DegenerateSubspace), intent(in)           :: subspace
      type(AnharmonicData),     intent(in)           :: anharmonic_data
      complex(dp)                                    :: output
    end function
    
    impure elemental function kinetic_energy_SubspaceBasis(this,bra,ket, &
       & subspace,anharmonic_data) result(output)
      import SubspaceBasis
      import BasisState
      import DegenerateSubspace
      import AnharmonicData
      import dp
      implicit none
      
      class(SubspaceBasis),     intent(in)           :: this
      class(BasisState),        intent(in)           :: bra
      class(BasisState),        intent(in), optional :: ket
      type(DegenerateSubspace), intent(in)           :: subspace
      type(AnharmonicData),     intent(in)           :: anharmonic_data
      real(dp)                                       :: output
    end function
    
    impure elemental function harmonic_potential_energy_SubspaceBasis(this, &
       & bra,ket,subspace,anharmonic_data) result(output)
      import SubspaceBasis
      import BasisState
      import DegenerateSubspace
      import AnharmonicData
      import dp
      implicit none
      
      class(SubspaceBasis),     intent(in)           :: this
      class(BasisState),        intent(in)           :: bra
      class(BasisState),        intent(in), optional :: ket
      type(DegenerateSubspace), intent(in)           :: subspace
      type(AnharmonicData),     intent(in)           :: anharmonic_data
      real(dp)                                       :: output
    end function
    
    impure elemental function kinetic_stress_SubspaceBasis(this,bra,ket, &
       & subspace,stress_prefactors,anharmonic_data) result(output)
      import SubspaceBasis
      import BasisState
      import DegenerateSubspace
      import StressPrefactors
      import AnharmonicData
      import RealMatrix
      implicit none
      
      class(SubspaceBasis),     intent(in)           :: this
      class(BasisState),        intent(in)           :: bra
      class(BasisState),        intent(in), optional :: ket
      type(DegenerateSubspace), intent(in)           :: subspace
      type(StressPrefactors),   intent(in)           :: stress_prefactors
      type(AnharmonicData),     intent(in)           :: anharmonic_data
      type(RealMatrix)                               :: output
    end function
    
    impure elemental function thermodynamic_data_SubspaceBasis(this, &
       & thermal_energy,states,subspace,subspace_potential,          &
       & subspace_stress,stress_prefactors,anharmonic_data) result(output)
      import SubspaceBasis
      import dp
      import BasisStates
      import DegenerateSubspace
      import PotentialData
      import StressData
      import StressPrefactors
      import AnharmonicData
      import ThermodynamicData
      implicit none
      
      class(SubspaceBasis),     intent(in)           :: this
      real(dp),                 intent(in)           :: thermal_energy
      class(BasisStates),       intent(in)           :: states
      type(DegenerateSubspace), intent(in)           :: subspace
      class(PotentialData),     intent(in)           :: subspace_potential
      class(StressData),        intent(in), optional :: subspace_stress
      type(StressPrefactors),   intent(in), optional :: stress_prefactors
      type(AnharmonicData),     intent(in)           :: anharmonic_data
      type(ThermodynamicData)                        :: output
    end function
    
    impure elemental function wavefunctions_SubspaceBasis(this,states, &
       & subspace,anharmonic_data) result(output)
      import SubspaceBasis
      import BasisStates
      import DegenerateSubspace
      import AnharmonicData
      import SubspaceWavefunctionsPointer
      implicit none
      
      class(SubspaceBasis),     intent(in) :: this
      class(BasisStates),       intent(in) :: states
      type(DegenerateSubspace), intent(in) :: subspace
      type(AnharmonicData),     intent(in) :: anharmonic_data
      type(SubspaceWavefunctionsPointer)   :: output
    end function
    
    impure elemental function integrate_BasisStates_SubspaceBasis(this, &
       & states,thermal_energy,monomial,subspace,anharmonic_data)       &
       & result(output)
      import SubspaceBasis
      import BasisStates
      import dp
      import SparseMonomial
      import DegenerateSubspace
      import QpointData
      import AnharmonicData
      implicit none
      
      class(SubspaceBasis),     intent(in) :: this
      class(BasisStates),       intent(in) :: states
      real(dp),                 intent(in) :: thermal_energy
      type(SparseMonomial),     intent(in) :: monomial
      type(DegenerateSubspace), intent(in) :: subspace
      type(AnharmonicData),     intent(in) :: anharmonic_data
      complex(dp)                          :: output
    end function
    
    ! PotentialData procedures.
    impure elemental function representation_PotentialData() result(output)
      import String
      implicit none
      
      type(String) :: output
    end function
    
    subroutine generate_sampling_points_PotentialData(this,anharmonic_data, &
       & use_forces,use_hessians,calculate_stress,sampling_points_dir,      &
       & calculation_writer,logfile)
      import PotentialData
      import AnharmonicData
      import String
      import CalculationWriter
      import OFile
      implicit none
      
      class(PotentialData),    intent(inout) :: this
      type(AnharmonicData),    intent(in)    :: anharmonic_data
      logical,                 intent(in)    :: use_forces
      logical,                 intent(in)    :: use_hessians
      logical,                 intent(in)    :: calculate_stress
      type(String),            intent(in)    :: sampling_points_dir
      type(CalculationWriter), intent(inout) :: calculation_writer
      type(OFile),             intent(inout) :: logfile
    end subroutine
    
    subroutine generate_potential_PotentialData(this,anharmonic_data,        &
       & weighted_energy_force_ratio,sampling_points_dir,calculation_reader, &
       & logfile)
      import dp
      import PotentialData
      import AnharmonicData
      import String
      import CalculationReader
      import OFile
      implicit none
      
      class(PotentialData),    intent(inout) :: this
      type(AnharmonicData),    intent(in)    :: anharmonic_data
      real(dp),                intent(in)    :: weighted_energy_force_ratio
      type(String),            intent(in)    :: sampling_points_dir
      type(CalculationReader), intent(inout) :: calculation_reader
      type(OFile),             intent(inout) :: logfile
    end subroutine
    
    function generate_stress_PotentialData(this,anharmonic_data, &
       & sampling_points_dir,stress_expansion_order,             &
       & stress_subspace_coupling,vscf_basis_functions_only,     &
       & calculation_reader,logfile) result(output)
      import dp
      import PotentialData
      import AnharmonicData
      import String
      import SubspaceCoupling
      import CalculationReader
      import OFile
      import StressPointer
      implicit none
      
      class(PotentialData),    intent(in)    :: this
      type(AnharmonicData),    intent(in)    :: anharmonic_data
      type(String),            intent(in)    :: sampling_points_dir
      integer,                 intent(in)    :: stress_expansion_order
      type(SubspaceCoupling),  intent(in)    :: stress_subspace_coupling(:)
      logical,                 intent(in)    :: vscf_basis_functions_only
      type(CalculationReader), intent(inout) :: calculation_reader
      type(OFile),             intent(inout) :: logfile
      type(StressPointer)                    :: output
    end function
    
    impure elemental subroutine zero_energy_PotentialData(this)
      import PotentialData
      implicit none
      
      class(PotentialData), intent(inout) :: this
    end subroutine
    
    impure elemental subroutine add_constant_PotentialData(this,input)
      import PotentialData
      import dp
      implicit none
      
      class(PotentialData), intent(inout) :: this
      real(dp),             intent(in)    :: input
    end subroutine
    
    impure elemental function energy_RealModeDisplacement_PotentialData(this, &
       & displacement) result(output)
      import dp
      import PotentialData
      import RealModeDisplacement
      implicit none
      
      class(PotentialData),       intent(in) :: this
      type(RealModeDisplacement), intent(in) :: displacement
      real(dp)                               :: output
    end function
    
    impure elemental function energy_ComplexModeDisplacement_PotentialData( &
       & this,displacement) result(output)
      import dp
      import PotentialData
      import ComplexModeDisplacement
      implicit none
      
      class(PotentialData),          intent(in) :: this
      type(ComplexModeDisplacement), intent(in) :: displacement
      complex(dp)                               :: output
    end function
    
    impure elemental function force_RealModeDisplacement_PotentialData(this, &
       & displacement) result(output)
      import PotentialData
      import RealModeDisplacement
      import RealModeForce
      implicit none
      
      class(PotentialData),       intent(in) :: this
      type(RealModeDisplacement), intent(in) :: displacement
      type(RealModeForce)                    :: output
    end function
    
    impure elemental function force_ComplexModeDisplacement_PotentialData( &
       & this,displacement) result(output)
      import PotentialData
      import ComplexModeDisplacement
      import ComplexModeForce
      implicit none
      
      class(PotentialData),          intent(in) :: this
      type(ComplexModeDisplacement), intent(in) :: displacement
      type(ComplexModeForce)                    :: output
    end function
    
    subroutine braket_SubspaceState_PotentialData(this,bra,ket, &
       & whole_subspace,anharmonic_data)
      import PotentialData
      import SubspaceState
      import AnharmonicData
      implicit none
      
      class(PotentialData), intent(inout)        :: this
      class(SubspaceState), intent(in)           :: bra
      class(SubspaceState), intent(in), optional :: ket
      logical,              intent(in), optional :: whole_subspace
      type(AnharmonicData), intent(in)           :: anharmonic_data
    end subroutine
    
    subroutine braket_BasisState_PotentialData(this,bra,ket,subspace, &
       & subspace_basis,whole_subspace,anharmonic_data)
      import PotentialData
      import BasisState
      import DegenerateSubspace
      import SubspaceBasis
      import AnharmonicData
      import QpointData
      implicit none
      
      class(PotentialData),     intent(inout)        :: this
      class(BasisState),        intent(in)           :: bra
      class(BasisState),        intent(in), optional :: ket
      type(DegenerateSubspace), intent(in)           :: subspace
      class(SubspaceBasis),     intent(in)           :: subspace_basis
      logical,                  intent(in), optional :: whole_subspace
      type(AnharmonicData),     intent(in)           :: anharmonic_data
    end subroutine
    
    subroutine braket_BasisStates_PotentialData(this,states,thermal_energy, &
       & subspace,subspace_basis,whole_subspace,anharmonic_data)
      import PotentialData
      import BasisStates
      import dp
      import DegenerateSubspace
      import SubspaceBasis
      import AnharmonicData
      implicit none
      
      class(PotentialData),     intent(inout)        :: this
      class(BasisStates),       intent(in)           :: states
      real(dp),                 intent(in)           :: thermal_energy
      type(DegenerateSubspace), intent(in)           :: subspace
      class(SubspaceBasis),     intent(in)           :: subspace_basis
      logical,                  intent(in), optional :: whole_subspace
      type(AnharmonicData),     intent(in)           :: anharmonic_data
    end subroutine
    
    function harmonic_expectation_PotentialData(this,frequency, &
       & thermal_energy,anharmonic_data) result(output)
      import PotentialData
      import dp
      import AnharmonicData
      implicit none
      
      class(PotentialData), intent(in) :: this
      real(dp),             intent(in) :: frequency
      real(dp),             intent(in) :: thermal_energy
      type(AnharmonicData), intent(in) :: anharmonic_data
      real(dp)                         :: output
    end function
    
    function coefficients_PotentialData(this,frequency,anharmonic_data) &
       & result(output)
      import PotentialData
      import dp
      import AnharmonicData
      implicit none
      
      class(potentialData), intent(in) :: this
      real(dp),             intent(in) :: frequency
      type(AnharmonicData), intent(in) :: anharmonic_data
      real(dp), allocatable            :: output(:)
    end function
    
    subroutine set_coefficients_PotentialData(this,coefficients,frequency, &
       & anharmonic_data)
      import PotentialData
      import dp
      import AnharmonicData
      implicit none
      
      class(PotentialData), intent(inout) :: this
      real(dp),             intent(in)    :: coefficients(:)
      real(dp),             intent(in)    :: frequency
      type(AnharmonicData), intent(in)    :: anharmonic_data
    end subroutine
    
    ! StressData procedures.
    impure elemental function representation_StressData() result(output)
      import String
      implicit none
      
      type(String) :: output
    end function
    
    impure elemental subroutine zero_stress_StressData(this)
      import StressData
      implicit none
      
      class(StressData), intent(inout) :: this
    end subroutine
    
    impure elemental subroutine add_constant_StressData(this,input)
      import StressData
      import RealMatrix
      implicit none
      
      class(StressData), intent(inout) :: this
      type(RealMatrix),  intent(in)    :: input
    end subroutine
    
    impure elemental function stress_RealModeDisplacement_StressData(this, &
       & displacement) result(output)
      import StressData
      import RealModeDisplacement
      import RealMatrix
      implicit none
      
      class(StressData),          intent(in) :: this
      type(RealModeDisplacement), intent(in) :: displacement
      type(RealMatrix)                       :: output
    end function
    
    impure elemental function stress_ComplexModeDisplacement_StressData( &
       & this,displacement) result(output)
      import StressData
      import ComplexModeDisplacement
      import ComplexMatrix
      implicit none
      
      class(StressData),             intent(in) :: this
      type(ComplexModeDisplacement), intent(in) :: displacement
      type(ComplexMatrix)                       :: output
    end function
    
    subroutine braket_SubspaceState_StressData(this,bra,ket,whole_subspace, &
       & anharmonic_data)
      import StressData
      import SubspaceState
      import AnharmonicData
      implicit none
      
      class(StressData),    intent(inout)        :: this
      class(SubspaceState), intent(in)           :: bra
      class(SubspaceState), intent(in), optional :: ket
      logical,              intent(in), optional :: whole_subspace
      type(AnharmonicData), intent(in)           :: anharmonic_data
    end subroutine
    
    subroutine braket_BasisState_StressData(this,bra,ket,subspace, &
       & subspace_basis,whole_subspace,anharmonic_data)
      import StressData
      import BasisState
      import DegenerateSubspace
      import SubspaceBasis
      import AnharmonicData
      implicit none
      
      class(StressData),        intent(inout)        :: this
      class(BasisState),        intent(in)           :: bra
      class(BasisState),        intent(in), optional :: ket
      type(DegenerateSubspace), intent(in)           :: subspace
      class(SubspaceBasis),     intent(in)           :: subspace_basis
      logical,                  intent(in), optional :: whole_subspace
      type(AnharmonicData),     intent(in)           :: anharmonic_data
    end subroutine
    
    subroutine braket_BasisStates_StressData(this,states,thermal_energy, &
       & subspace,subspace_basis,whole_subspace,anharmonic_data)
      import StressData
      import BasisStates
      import dp
      import DegenerateSubspace
      import SubspaceBasis
      import AnharmonicData
      implicit none
      
      class(StressData),        intent(inout)        :: this
      class(BasisStates),       intent(in)           :: states
      real(dp),                 intent(in)           :: thermal_energy
      type(DegenerateSubspace), intent(in)           :: subspace
      class(SubspaceBasis),     intent(in)           :: subspace_basis
      logical,                  intent(in), optional :: whole_subspace
      type(AnharmonicData),     intent(in)           :: anharmonic_data
    end subroutine
    
    function harmonic_expectation_StressData(this,frequency, &
       & thermal_energy,anharmonic_data) result(output)
      import StressData
      import dp
      import AnharmonicData
      import RealMatrix
      implicit none
      
      class(StressData),    intent(in) :: this
      real(dp),             intent(in) :: frequency
      real(dp),             intent(in) :: thermal_energy
      type(AnharmonicData), intent(in) :: anharmonic_data
      type(RealMatrix)                 :: output
    end function
  end interface
  
  ! --------------------------------------------------
  ! Pointer constructor interfaces.
  ! --------------------------------------------------
  interface SubspaceBasisPointer
    module procedure new_SubspaceBasisPointer
    module procedure new_SubspaceBasisPointer_Strings
    module procedure new_SubspaceBasisPointer_StringArray
  end interface
  
  interface PotentialPointer
    module procedure new_PotentialPointer
    module procedure new_PotentialPointer_Strings
    module procedure new_PotentialPointer_StringArray
  end interface
  
  interface StressPointer
    module procedure new_StressPointer
    module procedure new_StressPointer_Strings
    module procedure new_StressPointer_StringArray
  end interface
contains

! ----------------------------------------------------------------------
! Startup methods.
! ----------------------------------------------------------------------
subroutine startup_SubspaceBasis(this)
  implicit none
  
  class(SubspaceBasis), intent(in) :: this
  
  integer :: i
  
  if (.not.allocated(TYPES_SubspaceBasis)) then
    TYPES_SubspaceBasis = [SubspaceBasisPointer(this)]
  elseif (.not.any([(                                        &
     & this%representation()                                 &
     &    == TYPES_SubspaceBasis(i)%basis_%representation(), &
     & i=1,                                                  &
     & size(TYPES_SubspaceBasis)                             )])) then
    TYPES_SubspaceBasis = [TYPES_SubspaceBasis, SubspaceBasisPointer(this)]
  endif
end subroutine

subroutine startup_PotentialData(this)
  implicit none
  
  class(PotentialData), intent(in) :: this
  
  integer :: i
  
  if (.not.allocated(TYPES_PotentialData)) then
    TYPES_PotentialData = [PotentialPointer(this)]
  elseif (.not.any([(                                            &
     & this%representation()                                     &
     &    == TYPES_PotentialData(i)%potential_%representation(), &
     & i=1,                                                      &
     & size(TYPES_PotentialData)                                 )])) then
    TYPES_PotentialData = [TYPES_PotentialData, PotentialPointer(this)]
  endif
end subroutine

subroutine startup_StressData(this)
  implicit none
  
  class(StressData), intent(in) :: this
  
  integer :: i
  
  if (.not.allocated(TYPES_StressData)) then
    TYPES_StressData = [StressPointer(this)]
  elseif (.not.any([(                                      &
     & this%representation()                               &
     &    == TYPES_StressData(i)%stress_%representation(), &
     & i=1,                                                &
     & size(TYPES_StressData)                              )])) then
    TYPES_StressData = [TYPES_StressData, StressPointer(this)]
  endif
end subroutine

! ----------------------------------------------------------------------
! SubspaceBasisPointer methods.
! ----------------------------------------------------------------------
! Construct a SubspaceBasisPointer from any type which extends SubspaceBasis.
impure elemental function new_SubspaceBasisPointer(basis) result(this)
  implicit none
  
  class(SubspaceBasis), intent(in) :: basis
  type(SubspaceBasisPointer)       :: this
  
  integer :: ialloc
  
  select type(basis); type is(SubspaceBasisPointer)
    this = basis
  class default
    this%representation_ = basis%representation()
    allocate( this%basis_, source=basis, &
            & stat=ialloc); call err(ialloc)
  end select
end function

! Checks that the pointer has been allocated before it is used.
subroutine check_SubspaceBasisPointer(this)
  implicit none
  
  class(SubspaceBasisPointer), intent(in) :: this
  
  if (.not. allocated(this%basis_)) then
    call print_line(CODE_ERROR//': Trying to use a SubspaceBasisPointer &
       &before it has been allocated.')
    call err()
  endif
end subroutine

! Type representation.
impure elemental function representation_SubspaceBasisPointer() result(output)
  implicit none
  
  type(String) :: output
  
  output = 'pointer'
end function

! SubspaceBasis methods.
function basis_SubspaceBasisPointer(this) result(output)
  implicit none
  
  class(SubspaceBasisPointer), intent(in) :: this
  class(SubspaceBasis), allocatable       :: output
  
  output = this%basis_
end function

impure elemental function initial_states_SubspaceBasisPointer(this,subspace, &
   & anharmonic_data) result(output)
  implicit none
  
  class(SubspaceBasisPointer), intent(in) :: this
  type(DegenerateSubspace),    intent(in) :: subspace
  type(AnharmonicData),        intent(in) :: anharmonic_data
  type(BasisStatesPointer)                :: output
  
  call this%check()
  
  output = this%basis_%initial_states(subspace, anharmonic_data)
end function

impure elemental function calculate_states_SubspaceBasisPointer(this,     &
   & subspace,subspace_potential,thermal_energy,energy_convergence,       &
   & no_converged_calculations,max_pulay_iterations,pre_pulay_iterations, &
   & pre_pulay_damping,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceBasisPointer), intent(in) :: this
  type(DegenerateSubspace),    intent(in) :: subspace
  class(PotentialData),        intent(in) :: subspace_potential
  real(dp),                    intent(in) :: thermal_energy
  real(dp),                    intent(in) :: energy_convergence
  integer,                     intent(in) :: no_converged_calculations
  integer,                     intent(in) :: max_pulay_iterations
  integer,                     intent(in) :: pre_pulay_iterations
  real(dp),                    intent(in) :: pre_pulay_damping
  type(AnharmonicData),        intent(in) :: anharmonic_data
  type(BasisStatesPointer)                :: output
  
  call this%check()
  
  output = this%basis_%calculate_states( subspace,                  &
                                       & subspace_potential,        &
                                       & thermal_energy,            &
                                       & energy_convergence,        &
                                       & no_converged_calculations, &
                                       & max_pulay_iterations,      &
                                       & pre_pulay_iterations,      &
                                       & pre_pulay_damping,         &
                                       & anharmonic_data            )
end function

impure elemental function process_subspace_potential_SubspaceBasisPointer( &
   & this,potential,states,subspace,thermal_energy,anharmonic_data)        &
   & result(output)
  implicit none
  
  class(SubspaceBasisPointer), intent(in) :: this
  class(PotentialData),        intent(in) :: potential
  class(BasisStates),          intent(in) :: states
  type(DegenerateSubspace),    intent(in) :: subspace
  real(dp),                    intent(in) :: thermal_energy
  type(AnharmonicData),        intent(in) :: anharmonic_data
  type(PotentialPointer)                  :: output
  
  call this%check()
  
  output = this%basis_%process_subspace_potential( potential,      &
                                                 & states,         &
                                                 & subspace,       &
                                                 & thermal_energy, &
                                                 & anharmonic_data )
end function

impure elemental function process_subspace_stress_SubspaceBasisPointer(this, &
   & stress,states,subspace,thermal_energy,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceBasisPointer), intent(in) :: this
  class(stressData),           intent(in) :: stress
  class(BasisStates),          intent(in) :: states
  type(DegenerateSubspace),    intent(in) :: subspace
  real(dp),                    intent(in) :: thermal_energy
  type(AnharmonicData),        intent(in) :: anharmonic_data
  type(stressPointer)                     :: output
  
  call this%check()
  
  output = this%basis_%process_subspace_stress( stress,         &
                                              & states,         &
                                              & subspace,       &
                                              & thermal_energy, &
                                              & anharmonic_data )
end function

function mode_ids_SubspaceBasisPointer(this,subspace,anharmonic_data) &
   & result(output)
  implicit none
  
  class(SubspaceBasisPointer), intent(in) :: this
  type(DegenerateSubspace),    intent(in) :: subspace
  type(AnharmonicData),        intent(in) :: anharmonic_data
  integer, allocatable                    :: output(:)
  
  call this%check()
  
  output = this%basis_%mode_ids(subspace,anharmonic_data)
end function

function paired_mode_ids_SubspaceBasisPointer(this,subspace,anharmonic_data) &
   & result(output)
  implicit none
  
  class(SubspaceBasisPointer), intent(in) :: this
  type(DegenerateSubspace),    intent(in) :: subspace
  type(AnharmonicData),        intent(in) :: anharmonic_data
  integer, allocatable                    :: output(:)
  
  call this%check()
  
  output = this%basis_%paired_mode_ids(subspace,anharmonic_data)
end function

impure elemental function inner_product_SubspaceBasisPointer(this,bra,ket, &
   & subspace,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceBasisPointer), intent(in)           :: this
  class(BasisState),           intent(in)           :: bra
  class(BasisState),           intent(in), optional :: ket
  type(DegenerateSubspace),    intent(in)           :: subspace
  type(AnharmonicData),        intent(in)           :: anharmonic_data
  real(dp)                                          :: output
  
  call this%check()
  
  output = this%basis_%inner_product( bra,            &
                                    & ket,            &
                                    & subspace,       &
                                    & anharmonic_data )
end function

impure elemental function integrate_BasisState_SubspaceBasisPointer( &
   & this,bra,monomial,ket,subspace,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceBasisPointer), intent(in)           :: this
  class(BasisState),           intent(in)           :: bra
  type(SparseMonomial),        intent(in)           :: monomial
  class(BasisState),           intent(in), optional :: ket
  type(DegenerateSubspace),    intent(in)           :: subspace
  type(AnharmonicData),        intent(in)           :: anharmonic_data
  complex(dp)                                       :: output
  
  call this%check()
  
  output = this%basis_%integrate( bra,            &
                                & monomial,       &
                                & ket,            &
                                & subspace,       &
                                & anharmonic_data )
end function

impure elemental function kinetic_energy_SubspaceBasisPointer(this,bra,ket, &
   & subspace,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceBasisPointer), intent(in)           :: this
  class(BasisState),           intent(in)           :: bra
  class(BasisState),           intent(in), optional :: ket
  type(DegenerateSubspace),    intent(in)           :: subspace
  type(AnharmonicData),        intent(in)           :: anharmonic_data
  real(dp)                                          :: output
  
  call this%check()
  
  output = this%basis_%kinetic_energy( bra,            &
                                     & ket,            &
                                     & subspace,       &
                                     & anharmonic_data )
end function

impure elemental function harmonic_potential_energy_SubspaceBasisPointer( &
   & this,bra,ket,subspace,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceBasisPointer), intent(in)           :: this
  class(BasisState),           intent(in)           :: bra
  class(BasisState),           intent(in), optional :: ket
  type(DegenerateSubspace),    intent(in)           :: subspace
  type(AnharmonicData),        intent(in)           :: anharmonic_data
  real(dp)                                          :: output
  
  call this%check()
  
  output = this%basis_%harmonic_potential_energy( bra,            &
                                                & ket,            &
                                                & subspace,       &
                                                & anharmonic_data )
end function

impure elemental function kinetic_stress_SubspaceBasisPointer(this,bra,ket, &
   & subspace,stress_prefactors,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceBasisPointer), intent(in)           :: this
  class(BasisState),           intent(in)           :: bra
  class(BasisState),           intent(in), optional :: ket
  type(DegenerateSubspace),    intent(in)           :: subspace
  type(StressPrefactors),      intent(in)           :: stress_prefactors
  type(AnharmonicData),        intent(in)           :: anharmonic_data
  type(RealMatrix)                                  :: output
  
  call this%check()
  
  output = this%basis_%kinetic_stress( bra,               &
                                     & ket,               &
                                     & subspace,          &
                                     & stress_prefactors, &
                                     & anharmonic_data    )
end function

impure elemental function thermodynamic_data_SubspaceBasisPointer(this, &
   & thermal_energy,states,subspace,subspace_potential,subspace_stress, &
   & stress_prefactors,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceBasisPointer), intent(in)           :: this
  real(dp),                    intent(in)           :: thermal_energy
  class(BasisStates),          intent(in)           :: states
  type(DegenerateSubspace),    intent(in)           :: subspace
  class(PotentialData),        intent(in)           :: subspace_potential
  class(StressData),           intent(in), optional :: subspace_stress
  type(StressPrefactors),      intent(in), optional :: stress_prefactors
  type(AnharmonicData),        intent(in)           :: anharmonic_data
  type(ThermodynamicData)                           :: output
  
  call this%check()
  
  output = this%basis_%thermodynamic_data( thermal_energy,     &
                                         & states,             &
                                         & subspace,           &
                                         & subspace_potential, &
                                         & subspace_stress,    &
                                         & stress_prefactors,  &
                                         & anharmonic_data     )
end function

impure elemental function wavefunctions_SubspaceBasisPointer(this,states, &
   & subspace,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceBasisPointer), intent(in) :: this
  class(BasisStates),          intent(in) :: states
  type(DegenerateSubspace),    intent(in) :: subspace
  type(AnharmonicData),        intent(in) :: anharmonic_data
  type(SubspaceWavefunctionsPointer)      :: output
  
  call this%check()
  
  output = this%basis_%wavefunctions( states,         &
                                    & subspace,       &
                                    & anharmonic_data )
end function

impure elemental function integrate_BasisStates_SubspaceBasisPointer( &
   & this,states,thermal_energy,monomial,subspace,anharmonic_data)    &
   & result(output)
  implicit none
  
  class(SubspaceBasisPointer), intent(in) :: this
  class(BasisStates),          intent(in) :: states
  real(dp),                    intent(in) :: thermal_energy
  type(SparseMonomial),        intent(in) :: monomial
  type(DegenerateSubspace),    intent(in) :: subspace
  type(AnharmonicData),        intent(in) :: anharmonic_data
  complex(dp)                             :: output
  
  call this%check()
  
  output = this%basis_%integrate( states,         &
                                & thermal_energy, &
                                & monomial,       &
                                & subspace,       &
                                & anharmonic_data )
end function

! I/O.
subroutine read_SubspaceBasisPointer(this,input)
  implicit none
  
  class(SubspaceBasisPointer), intent(out) :: this
  type(String),                intent(in)  :: input(:)
  
  type(String), allocatable :: line(:)
  
  type(String) :: representation
  
  integer :: i
  
  select type(this); type is(SubspaceBasisPointer)
    line = split_line(input(1))
    representation = line(3)
    
    ! Identify which type corresponds to the representation.
    i = first([(                                                         &
       & TYPES_SubspaceBasis(i)%basis_%representation()==representation, &
       & i=1,                                                            &
       & size(TYPES_SubspaceBasis)                                       )])
    
    ! Read the input into the element of the correct type,
    !    and copy that element into the output.
    call TYPES_SubspaceBasis(i)%basis_%read(input(2:))
    this = SubspaceBasisPointer(TYPES_SubspaceBasis(i))
  class default
    call err()
  end select
end subroutine

function write_SubspaceBasisPointer(this) result(output)
  implicit none
  
  class(SubspaceBasisPointer), intent(in) :: this
  type(String), allocatable               :: output(:)
  
  select type(this); type is(SubspaceBasisPointer)
    output = [ 'SubspaceBasis representation: '//this%representation_, &
             & str(this%basis_)                                        ]
  end select
end function

function new_SubspaceBasisPointer_Strings(input) result(this)
  implicit none
  
  type(String), intent(in)   :: input(:)
  type(SubspaceBasisPointer) :: this
  
  call this%read(input)
end function

impure elemental function new_SubspaceBasisPointer_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(SubspaceBasisPointer)    :: this
  
  this = SubspaceBasisPointer(str(input))
end function

! ----------------------------------------------------------------------
! PotentialPointer methods.
! ----------------------------------------------------------------------
! Construct a PotentialPointer from any type which extends PotentialData.
impure elemental function new_PotentialPointer(potential) result(this)
  implicit none
  
  class(PotentialData), intent(in) :: potential
  type(PotentialPointer)           :: this
  
  integer :: ialloc
  
  select type(potential); type is(PotentialPointer)
    this = potential
  class default
    this%representation_ = potential%representation()
    allocate( this%potential_, source=potential, &
            & stat=ialloc); call err(ialloc)
  end select
end function

! Checks that the pointer has been allocated before it is used.
impure elemental subroutine check_PotentialPointer(this)
  implicit none
  
  class(PotentialPointer), intent(in) :: this
  
  if (.not. allocated(this%potential_)) then
    call print_line(CODE_ERROR//': Trying to use a PotentialPointer before &
       &it has been allocated.')
    call err()
  endif
end subroutine

! Type representation.
impure elemental function representation_PotentialPointer() result(output)
  implicit none
  
  type(String) :: output
  
  output = 'pointer'
end function

! Return the pointed-to potential.
function potential_PotentialPointer(this) result(output)
  implicit none
  
  class(PotentialPointer), intent(in) :: this
  class(PotentialData), allocatable   :: output
  
  output = this%potential_
end function

! Wrappers for all of PotentialData's methods.
subroutine generate_sampling_points_PotentialPointer(this,anharmonic_data, &
   & use_forces,use_hessians,calculate_stress,sampling_points_dir,         &
   & calculation_writer,logfile)
  implicit none
  
  class(PotentialPointer), intent(inout) :: this
  type(AnharmonicData),    intent(in)    :: anharmonic_data
  logical,                 intent(in)    :: use_forces
  logical,                 intent(in)    :: use_hessians
  logical,                 intent(in)    :: calculate_stress
  type(String),            intent(in)    :: sampling_points_dir
  type(CalculationWriter), intent(inout) :: calculation_writer
  type(OFile),             intent(inout) :: logfile
  
  call this%check()
  
  call this%potential_%generate_sampling_points( anharmonic_data,     &
                                               & use_forces,          &
                                               & use_hessians,        &
                                               & calculate_stress,    &
                                               & sampling_points_dir, &
                                               & calculation_writer,  &
                                               & logfile              )
end subroutine

subroutine generate_potential_PotentialPointer(this,anharmonic_data,     &
   & weighted_energy_force_ratio,sampling_points_dir,calculation_reader, &
   & logfile)
  implicit none
  
  class(PotentialPointer), intent(inout) :: this
  type(AnharmonicData),    intent(in)    :: anharmonic_data
  real(dp),                intent(in)    :: weighted_energy_force_ratio
  type(String),            intent(in)    :: sampling_points_dir
  type(CalculationReader), intent(inout) :: calculation_reader
  type(OFile),             intent(inout) :: logfile
  
  call this%check()
  
  call this%potential_%generate_potential( anharmonic_data,              &
                                         & weighted_energy_force_ratio,  &
                                         & sampling_points_dir,          &
                                         & calculation_reader,           &
                                         & logfile                       )
end subroutine

function generate_stress_PotentialPointer(this,anharmonic_data,           &
   & sampling_points_dir,stress_expansion_order,stress_subspace_coupling, &
   & vscf_basis_functions_only,calculation_reader,logfile) result(output)
  implicit none
  
  class(PotentialPointer), intent(in)    :: this
  type(AnharmonicData),    intent(in)    :: anharmonic_data
  type(String),            intent(in)    :: sampling_points_dir
  integer,                 intent(in)    :: stress_expansion_order
  type(SubspaceCoupling),  intent(in)    :: stress_subspace_coupling(:)
  logical,                 intent(in)    :: vscf_basis_functions_only
  type(CalculationReader), intent(inout) :: calculation_reader
  type(OFile),             intent(inout) :: logfile
  type(StressPointer)                    :: output
  
  call this%check()
  
  output = this%potential_%generate_stress( anharmonic_data,           &
                                          & sampling_points_dir,       &
                                          & stress_expansion_order,    &
                                          & stress_subspace_coupling,  &
                                          & vscf_basis_functions_only, &
                                          & calculation_reader,        &
                                          & logfile                    )
end function

impure elemental subroutine zero_energy_PotentialPointer(this)
  implicit none
  
  class(PotentialPointer), intent(inout) :: this
  
  call this%check()
  
  call this%potential_%zero_energy()
end subroutine

impure elemental subroutine add_constant_PotentialPointer(this,input)
  implicit none
  
  class(PotentialPointer), intent(inout) :: this
  real(dp),                intent(in)    :: input
  
  call this%check()
  
  call this%potential_%add_constant(input)
end subroutine

impure elemental subroutine finalise_subspace_potential_PotentialPointer( &
   & this,subspace,anharmonic_data)
  implicit none
  
  class(PotentialPointer),  intent(inout) :: this
  type(DegenerateSubspace), intent(in)    :: subspace
  type(AnharmonicData),     intent(in)    :: anharmonic_data
  
  call this%check()
  
  call this%potential_%finalise_subspace_potential(subspace,anharmonic_data)
end subroutine

impure elemental function energy_RealModeDisplacement_PotentialPointer(this, &
   & displacement) result(output)
  implicit none
  
  class(PotentialPointer),    intent(in) :: this
  type(RealModeDisplacement), intent(in) :: displacement
  real(dp)                               :: output
  
  call this%check()
  
  output = this%potential_%energy(displacement)
end function

impure elemental function energy_ComplexModeDisplacement_PotentialPointer( &
   & this,displacement) result(output)
  implicit none
  
  class(PotentialPointer),       intent(in) :: this
  type(ComplexModeDisplacement), intent(in) :: displacement
  complex(dp)                               :: output
  
  call this%check()
  
  output = this%potential_%energy(displacement)
end function

impure elemental function force_RealModeDisplacement_PotentialPointer(this, &
   & displacement) result(output)
  implicit none
  
  class(PotentialPointer),    intent(in) :: this
  type(RealModeDisplacement), intent(in) :: displacement
  type(RealModeForce)                    :: output
  
  call this%check()
  
  output = this%potential_%force(displacement)
end function

impure elemental function force_ComplexModeDisplacement_PotentialPointer( &
   & this,displacement) result(output)
  implicit none
  
  class(PotentialPointer),       intent(in) :: this
  type(ComplexModeDisplacement), intent(in) :: displacement
  type(ComplexModeForce)                    :: output
  
  call this%check()
  
  output = this%potential_%force(displacement)
end function

subroutine braket_SubspaceState_PotentialPointer(this,bra,ket, &
   & whole_subspace,anharmonic_data)
  implicit none
  
  class(PotentialPointer), intent(inout)        :: this
  class(SubspaceState),    intent(in)           :: bra
  class(SubspaceState),    intent(in), optional :: ket
  logical,                 intent(in), optional :: whole_subspace
  type(AnharmonicData),    intent(in)           :: anharmonic_data
  
  call this%check()
  
  call this%potential_%braket(bra,ket,whole_subspace,anharmonic_data)
end subroutine

subroutine braket_BasisState_PotentialPointer(this,bra,ket,subspace, &
   & subspace_basis,whole_subspace,anharmonic_data)
  implicit none
  
  class(PotentialPointer),  intent(inout)        :: this
  class(BasisState),        intent(in)           :: bra
  class(BasisState),        intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  logical,                  intent(in), optional :: whole_subspace
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  
  call this%check()
  
  call this%potential_%braket( bra,            &
                             & ket,            &
                             & subspace,       &
                             & subspace_basis, &
                             & whole_subspace, &
                             & anharmonic_data )
end subroutine

subroutine braket_BasisStates_PotentialPointer(this,states,thermal_energy, &
   & subspace,subspace_basis,whole_subspace,anharmonic_data)
  implicit none
  
  class(PotentialPointer),  intent(inout)        :: this
  class(BasisStates),       intent(in)           :: states
  real(dp),                 intent(in)           :: thermal_energy
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  logical,                  intent(in), optional :: whole_subspace
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  
  call this%check()
  
  call this%potential_%braket( states,         &
                             & thermal_energy, &
                             & subspace,       &
                             & subspace_basis, &
                             & whole_subspace, &
                             & anharmonic_data )
end subroutine

function harmonic_expectation_PotentialPointer(this,frequency, &
   & thermal_energy,anharmonic_data) result(output)
  implicit none
  
  class(PotentialPointer), intent(in) :: this
  real(dp),                intent(in) :: frequency
  real(dp),                intent(in) :: thermal_energy
  type(AnharmonicData),    intent(in) :: anharmonic_data
  real(dp)                            :: output
  
  call this%check()
  
  output = this%potential_%harmonic_expectation( frequency,      &
                                               & thermal_energy, &
                                               & anharmonic_data )
end function

! Gets or sets the potential from a set of real coefficients.
function coefficients_PotentialPointer(this,frequency,anharmonic_data) &
   & result(output)
  implicit none
  
  class(PotentialPointer), intent(in) :: this
  real(dp),                intent(in) :: frequency
  type(AnharmonicData),    intent(in) :: anharmonic_data
  real(dp), allocatable               :: output(:)
  
  call this%check()
  
  output = this%potential_%coefficients(frequency, anharmonic_data)
end function

subroutine set_coefficients_PotentialPointer(this,coefficients,frequency, &
   & anharmonic_data)
  implicit none
  
  class(PotentialPointer), intent(inout) :: this
  real(dp),                intent(in)    :: coefficients(:)
  real(dp),                intent(in)    :: frequency
  type(AnharmonicData),    intent(in)    :: anharmonic_data
  
  call this%check()
  
  call this%potential_%set_coefficients( coefficients,   &
                                       & frequency,      &
                                       & anharmonic_data )
end subroutine

! I/O.
subroutine read_PotentialPointer(this,input)
  implicit none
  
  class(PotentialPointer), intent(out) :: this
  type(String),            intent(in)  :: input(:)
  
  type(String), allocatable :: line(:)
  
  type(String) :: representation
  
  integer :: i
  
  select type(this); type is(PotentialPointer)
    line = split_line(input(1))
    representation = line(3)
    
    ! Identify which type corresponds to the representation.
    i = first([(                                                             &
       & TYPES_PotentialData(i)%potential_%representation()==representation, &
       & i=1,                                                                &
       & size(TYPES_PotentialData)                                          )])
    
    ! Read the input into the element of the correct type,
    !    and copy that element into the output.
    call TYPES_PotentialData(i)%potential_%read(input(2:))
    this = PotentialPointer(TYPES_PotentialData(i))
  class default
    call err()
  end select
end subroutine

function write_PotentialPointer(this) result(output)
  implicit none
  
  class(PotentialPointer), intent(in) :: this
  type(String), allocatable           :: output(:)
  
  select type(this); type is(PotentialPointer)
    output = [ 'Potential representation: '//this%representation_, &
             & str(this%potential_)                                ]
  end select
end function

function new_PotentialPointer_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(PotentialPointer)   :: this
  
  call this%read(input)
end function

impure elemental function new_PotentialPointer_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(PotentialPointer)        :: this
  
  this = PotentialPointer(str(input))
end function

! ----------------------------------------------------------------------
! StressPointer methods.
! ----------------------------------------------------------------------
! Construct a StressPointer from any type which extends StressData.
impure elemental function new_StressPointer(stress) result(this)
  implicit none
  
  class(StressData), intent(in) :: stress
  type(StressPointer)           :: this
  
  integer :: ialloc
  
  select type(stress); type is(StressPointer)
    this = stress
  class default
    this%representation_ = stress%representation()
    allocate( this%stress_, source=stress, &
            & stat=ialloc); call err(ialloc)
  end select
end function

! Checks that the pointer has been allocated before it is used.
impure elemental subroutine check_StressPointer(this)
  implicit none
  
  class(StressPointer), intent(in) :: this
  
  if (.not. allocated(this%stress_)) then
    call print_line(CODE_ERROR//': Trying to use a StressPointer before &
       &it has been allocated.')
    call err()
  endif
end subroutine

! Type representation.
impure elemental function representation_StressPointer() result(output)
  implicit none
  
  type(String) :: output
  
  output = 'pointer'
end function

impure elemental subroutine zero_stress_StressPointer(this)
  implicit none
  
  class(StressPointer), intent(inout) :: this
  
  call this%check()
  
  call this%stress_%zero_stress()
end subroutine

impure elemental subroutine add_constant_StressPointer(this,input)
  implicit none
  
  class(StressPointer), intent(inout) :: this
  type(RealMatrix),     intent(in)    :: input
  
  call this%check()
  
  call this%stress_%add_constant(input)
end subroutine

impure elemental function stress_RealModeDisplacement_StressPointer(this, &
   & displacement) result(output)
  implicit none
  
  class(StressPointer),       intent(in) :: this
  type(RealModeDisplacement), intent(in) :: displacement
  type(RealMatrix)                       :: output
  
  call this%check()
  
  output = this%stress_%stress(displacement)
end function

impure elemental function stress_ComplexModeDisplacement_StressPointer( &
   & this,displacement) result(output)
  implicit none
  
  class(StressPointer),          intent(in) :: this
  type(ComplexModeDisplacement), intent(in) :: displacement
  type(ComplexMatrix)                       :: output
  
  call this%check()
  
  output = this%stress_%stress(displacement)
end function

subroutine braket_SubspaceState_StressPointer(this,bra,ket,whole_subspace, &
   & anharmonic_data)
  implicit none
  
  class(StressPointer),     intent(inout)        :: this
  class(SubspaceState),     intent(in)           :: bra
  class(SubspaceState),     intent(in), optional :: ket
  logical,                  intent(in), optional :: whole_subspace
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  
  call this%check()
  
  call this%stress_%braket(bra,ket,whole_subspace,anharmonic_data)
end subroutine

subroutine braket_BasisState_StressPointer(this,bra,ket,subspace, &
   & subspace_basis,whole_subspace,anharmonic_data)
  implicit none
  
  class(StressPointer),     intent(inout)        :: this
  class(BasisState),        intent(in)           :: bra
  class(BasisState),        intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  logical,                  intent(in), optional :: whole_subspace
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  
  call this%check()
  
  call this%stress_%braket( bra,            &
                          & ket,            &
                          & subspace,       &
                          & subspace_basis, &
                          & whole_subspace, &
                          & anharmonic_data )
end subroutine

subroutine braket_BasisStates_StressPointer(this,states,thermal_energy, &
   & subspace,subspace_basis,whole_subspace,anharmonic_data)
  implicit none
  
  class(StressPointer),     intent(inout)        :: this
  class(BasisStates),       intent(in)           :: states
  real(dp),                 intent(in)           :: thermal_energy
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  logical,                  intent(in), optional :: whole_subspace
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  
  call this%check()
  
  call this%stress_%braket( states,         &
                          & thermal_energy, &
                          & subspace,       &
                          & subspace_basis, &
                          & whole_subspace, &
                          & anharmonic_data )
end subroutine

function harmonic_expectation_StressPointer(this,frequency, &
   & thermal_energy,anharmonic_data) result(output)
  implicit none
  
  class(StressPointer), intent(in) :: this
  real(dp),             intent(in) :: frequency
  real(dp),             intent(in) :: thermal_energy
  type(AnharmonicData), intent(in) :: anharmonic_data
  type(RealMatrix)                 :: output
  
  call this%check()
  
  output = this%stress_%harmonic_expectation( frequency,      &
                                            & thermal_energy, &
                                            & anharmonic_data )
end function

! I/O.
subroutine read_StressPointer(this,input)
  implicit none
  
  class(StressPointer), intent(out) :: this
  type(String),         intent(in)  :: input(:)
  
  type(String), allocatable :: line(:)
  
  type(String) :: representation
  
  integer :: i
  
  select type(this); type is(StressPointer)
    line = split_line(input(1))
    representation = line(3)
    
    ! Identify which type corresponds to the representation.
    i = first([(                                                       &
       & TYPES_StressData(i)%stress_%representation()==representation, &
       & i=1,                                                          &
       & size(TYPES_StressData)                                        )])
    
    ! Read the input into the element of the correct type,
    !    and copy that element into the output.
    call TYPES_StressData(i)%stress_%read(input(2:))
    this = StressPointer(TYPES_StressData(i))
  class default
    call err()
  end select
end subroutine

function write_StressPointer(this) result(output)
  implicit none
  
  class(StressPointer), intent(in) :: this
  type(String), allocatable        :: output(:)
  
  select type(this); type is(StressPointer)
    output = [ 'Stress representation: '//this%representation_, &
             & str(this%stress_)                                ]
  end select
end function

function new_StressPointer_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(StressPointer)      :: this
  
  call this%read(input)
end function

impure elemental function new_StressPointer_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(StressPointer)           :: this
  
  this = StressPointer(str(input))
end function

! ----------------------------------------------------------------------
! Concrete methods of abstract classes.
! ----------------------------------------------------------------------
! SubspaceBasis methods.
impure elemental function process_subspace_potential_SubspaceBasis(this, &
   & potential,states,subspace,thermal_energy,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceBasis),     intent(in) :: this
  class(PotentialData),     intent(in) :: potential
  class(BasisStates),       intent(in) :: states
  type(DegenerateSubspace), intent(in) :: subspace
  real(dp),                 intent(in) :: thermal_energy
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(PotentialPointer)               :: output
  
  output = PotentialPointer(potential)
end function

impure elemental function process_subspace_stress_SubspaceBasis(this,stress, &
   & states,subspace,thermal_energy,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceBasis),     intent(in) :: this
  class(stressData),        intent(in) :: stress
  class(BasisStates),       intent(in) :: states
  type(DegenerateSubspace), intent(in) :: subspace
  real(dp),                 intent(in) :: thermal_energy
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(stressPointer)                  :: output
  
  output = stressPointer(stress)
end function

! PotentialData methods.
function undisplaced_energy(this) result(output)
  implicit none
  
  class(PotentialData), intent(in) :: this
  real(dp)                         :: output
  
  type(RealSingleDisplacement) :: zero_displacement(0)
  
  output = this%energy(RealModeDisplacement(zero_displacement))
end function

impure elemental subroutine finalise_subspace_potential_PotentialData(this, &
   & subspace,anharmonic_data)
  implicit none
  
  class(PotentialData),     intent(inout) :: this
  type(DegenerateSubspace), intent(in)    :: subspace
  type(AnharmonicData),     intent(in)    :: anharmonic_data
  
  ! By default this doesn't do anything.
end subroutine

! StressData methods.
function undisplaced_stress(this) result(output)
  implicit none
  
  class(StressData), intent(in) :: this
  type(RealMatrix)              :: output
  
  type(RealSingleDisplacement) :: zero_displacement(0)
  
  output = this%stress(RealModeDisplacement(zero_displacement))
end function
end module
