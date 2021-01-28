! ======================================================================
! Four abstract classes, and pointers to those classes.
!    - SubspaceBasis defines the basis of states spanning a subspace.
!    - PotentialBase defines a potential energy surface.
!    - StressBase defines a stress surface.
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
!
! N.B. the Basis and State methods take other Basis and State class() arguments
!    with the TARGET attribute, so that they can be cast to their concrete type
!    pointers.
! The %state_pointer() and %states_pointer() methods should
!    not be called on objects which do not have the TARGET attribute,
!    as this will silently lead to undefined behaviour.
module caesar_abstract_classes_module
  use caesar_common_module
  
  use caesar_subspace_coupling_module
  use caesar_anharmonic_data_module
  use caesar_degenerate_symmetry_module
  use caesar_subspace_wavefunctions_module
  use caesar_stress_prefactors_module
  use caesar_subspace_state_module
  use caesar_subspace_braket_module
  use caesar_basis_state_module
  use caesar_basis_states_module
  use caesar_sparse_monomial_module
  use caesar_pulay_module
  implicit none
  
  private
  
  public :: SubspaceBasis
  public :: SubspaceBasisPointer
  
  public :: PotentialBase
  public :: PotentialBasePointer
  
  public :: StressBase
  public :: StressBasePointer
  
  ! ----------------------------------------------------------------------
  ! Abstract and pointer type definitions.
  ! ----------------------------------------------------------------------
  ! N.B. all objects are defined before any type-bound procedures are defined,
  !    to allow the type-bound procedures of all types to use all of the
  !    other types.
  
  type, abstract, extends(Stringsable) :: SubspaceBasis
    real(dp) :: frequency
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
    
    ! Calculate the derivative of the free energy with respect to a set of
    !    basis function coefficients.
    procedure(free_energy_gradient_SubspaceBasis), public, deferred :: &
       & free_energy_gradient
    
    ! Select the symmetries which are conserved under this basis.
    ! By default, this will simply return all of them.
    procedure, public :: select_symmetries => select_symmetries_SubspaceBasis
  end type
  
  type, extends(SubspaceBasis) :: SubspaceBasisPointer
    type(String),                      private :: representation_
    class(SubspaceBasis), allocatable, private :: basis_
  contains
    procedure, private :: check => check_SubspaceBasisPointer
    
    procedure, public, nopass :: representation => &
                               & representation_SubspaceBasisPointer
    
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
    
    procedure, public :: free_energy_gradient => &
                       & free_energy_gradient_SubspaceBasisPointer
    
    procedure, public :: select_symmetries => &
                       & select_symmetries_SubspaceBasisPointer
    
    ! I/O.
    procedure, public :: read  => read_SubspaceBasisPointer
    procedure, public :: write => write_SubspaceBasisPointer
  end type
  
  ! An array of all types which extend SubspaceBasis.
  ! This array will be filled in by startup routines.
  type(SubspaceBasisPointer), allocatable :: TYPES_SubspaceBasis(:)
  
  type, abstract, extends(Stringsable) :: PotentialBase
  contains
    procedure(representation_PotentialBase), public, deferred, nopass :: &
       & representation
    procedure, public :: startup => startup_PotentialBase
    
    ! Return the energy at zero displacement
    procedure, public :: undisplaced_energy
    
    ! Set the energy at zero displacement to zero,
    !    or add a constant to this energy.
    procedure(zero_energy_PotentialBase),  public, deferred :: zero_energy
    procedure(add_constant_PotentialBase), public, deferred :: add_constant
    
    ! Return the energy and force at a given real or complex displacement.
    generic, public :: energy =>                    &
                     & energy_RealModeDisplacement, &
                     & energy_ComplexModeDisplacement
    generic, public :: force  =>                   &
                     & force_RealModeDisplacement, &
                     & force_ComplexModeDisplacement
    procedure(energy_RealModeDisplacement_PotentialBase), public, &
       & deferred :: energy_RealModeDisplacement
    procedure(energy_ComplexModeDisplacement_PotentialBase), public, &
       & deferred :: energy_ComplexModeDisplacement
    procedure(force_RealModeDisplacement_PotentialBase), public, &
       & deferred :: force_RealModeDisplacement
    procedure(force_ComplexModeDisplacement_PotentialBase), public, &
       & deferred :: force_ComplexModeDisplacement
    
    ! Evaluate <bra|potential|ket>.
    generic, public :: braket =>              &
                     & braket_SubspaceBraKet, &
                     & braket_BasisState,     &
                     & braket_BasisStates
    procedure(braket_SubspaceBraKet_PotentialBase), public, deferred :: &
       & braket_SubspaceBraKet
    procedure(braket_BasisState_PotentialBase), public, deferred :: &
       & braket_BasisState
    procedure(braket_BasisStates_PotentialBase), public, deferred :: &
       & braket_BasisStates
    
    ! Evaluate the thermal expectation of the potential for a set of harmonic
    !    states.
    procedure(harmonic_expectation_PotentialBase), public, deferred :: &
       & harmonic_expectation
    
    ! Calculate the potential energy of the potential w/r/t a given state.
    generic,   public :: potential_energy =>              &
                       & potential_energy_SubspaceBraKet, &
                       & potential_energy_BasisState
    procedure, public :: potential_energy_SubspaceBraKet => &
                       & potential_energy_SubspaceBraKet_PotentialBase
    procedure, public :: potential_energy_BasisState => &
                       & potential_energy_BasisState_PotentialBase
  end type
  
  type, extends(PotentialBase) :: PotentialBasePointer
    type(String),                      private :: representation_
    class(PotentialBase), allocatable, private :: potential_
  contains
    procedure, private :: check => check_PotentialBasePointer
    
    procedure, public, nopass :: representation => &
                               & representation_PotentialBasePointer
    
    procedure, public :: potential => potential_PotentialBasePointer
    
    procedure, public :: zero_energy => zero_energy_PotentialBasePointer
    procedure, public :: add_constant => add_constant_PotentialBasePointer
    
    procedure, public :: energy_RealModeDisplacement => &
                       & energy_RealModeDisplacement_PotentialBasePointer
    procedure, public :: energy_ComplexModeDisplacement => &
                       & energy_ComplexModeDisplacement_PotentialBasePointer
    procedure, public :: force_RealModeDisplacement => &
                       & force_RealModeDisplacement_PotentialBasePointer
    procedure, public :: force_ComplexModeDisplacement => &
                       & force_ComplexModeDisplacement_PotentialBasePointer
    
    procedure, public :: braket_SubspaceBraKet  => &
                       & braket_SubspaceBraKet_PotentialBasePointer
    procedure, public :: braket_BasisState  => &
                       & braket_BasisState_PotentialBasePointer
    procedure, public :: braket_BasisStates => &
                       & braket_BasisStates_PotentialBasePointer
    
    procedure, public :: harmonic_expectation => &
                       & harmonic_expectation_PotentialBasePointer
    
    procedure, public :: potential_energy_SubspaceBraKet => &
                       & potential_energy_SubspaceBraKet_PotentialBasePointer
    procedure, public :: potential_energy_BasisState => &
                       & potential_energy_BasisState_PotentialBasePointer
    
    ! I/O.
    procedure, public :: read  => read_PotentialBasePointer
    procedure, public :: write => write_PotentialBasePointer
  end type
  
  ! An array of all types which extend PotentialBase.
  ! This array will be filled in by startup routines.
  type(PotentialBasePointer), allocatable :: TYPES_PotentialBase(:)
  
  type, abstract, extends(Stringsable) :: StressBase
  contains
    procedure(representation_StressBase), public, deferred, nopass :: &
       & representation
    procedure, public :: startup => startup_StressBase
    
    ! Return the stress at zero displacement.
    procedure, public :: undisplaced_stress
    
    ! Set the stress at zero displacement to zero,
    !    or add a constant to this stress.
    procedure(zero_stress_StressBase), public, deferred :: zero_stress
    procedure(add_constant_StressBase), public, deferred :: add_constant
    
    ! Return the stress at a given real or complex displacement.
    generic, public :: stress =>                    &
                     & stress_RealModeDisplacement, &
                     & stress_ComplexModeDisplacement
    procedure(stress_RealModeDisplacement_StressBase), public, &
       & deferred :: stress_RealModeDisplacement
    procedure(stress_ComplexModeDisplacement_StressBase), public, &
       & deferred :: stress_ComplexModeDisplacement
    
    ! Evaluate <bra|stress|ket>.
    generic, public :: braket =>              &
                     & braket_SubspaceBraKet, &
                     & braket_BasisState,     &
                     & braket_BasisStates
    procedure(braket_SubspaceBraKet_StressBase),  public, deferred :: &
       & braket_SubspaceBraKet
    procedure(braket_BasisState_StressBase),  public, deferred :: &
       & braket_BasisState
    procedure(braket_BasisStates_StressBase), public, deferred :: &
       & braket_BasisStates
    
    ! Evaluate the thermal expectation of the stress for a set of harmonic
    !    states.
    procedure(harmonic_expectation_StressBase), public, deferred :: &
       & harmonic_expectation
    
    ! Calculate the potential stress of the stress w/r/t a given state.
    generic,   public :: potential_stress =>              &
                       & potential_stress_SubspaceBraKet, &
                       & potential_stress_BasisState
    procedure, public :: potential_stress_SubspaceBraKet => &
                       & potential_stress_SubspaceBraKet_StressBase
    procedure, public :: potential_stress_BasisState => &
                       & potential_stress_BasisState_StressBase
  end type
  
  type, extends(StressBase) :: StressBasePointer
    type(String),                   private :: representation_
    class(StressBase), allocatable, private :: stress_
  contains
    procedure, private :: check => check_StressBasePointer
    
    procedure, public, nopass :: representation => &
                               & representation_StressBasePointer
    
    procedure, public :: zero_stress => zero_stress_StressBasePointer
    procedure, public :: add_constant => add_constant_StressBasePointer
    
    procedure, public :: stress_RealModeDisplacement => &
                       & stress_RealModeDisplacement_StressBasePointer
    procedure, public :: stress_ComplexModeDisplacement => &
                       & stress_ComplexModeDisplacement_StressBasePointer
    
    procedure, public :: braket_SubspaceBraKet  => &
                       & braket_SubspaceBraKet_StressBasePointer
    procedure, public :: braket_BasisState  => &
                       & braket_BasisState_StressBasePointer
    procedure, public :: braket_BasisStates => &
                       & braket_BasisStates_StressBasePointer
    
    procedure, public :: harmonic_expectation => &
                       & harmonic_expectation_StressBasePointer
    
    procedure, public :: potential_stress_SubspaceBraKet => &
                       & potential_stress_SubspaceBraKet_StressBasePointer
    procedure, public :: potential_stress_BasisState => &
                       & potential_stress_BasisState_StressBasePointer
    
    ! I/O.
    procedure, public :: read  => read_StressBasePointer
    procedure, public :: write => write_StressBasePointer
  end type
  
  ! An array of all types which extend StressBase.
  ! This array will be filled in by startup routines.
  type(StressBasePointer), allocatable :: TYPES_StressBase(:)
  
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
       & thermal_energy,anharmonic_data) result(output)
      import SubspaceBasis
      import DegenerateSubspace
      import dp
      import AnharmonicData
      import BasisStatesPointer
      implicit none
      
      class(SubspaceBasis),     intent(in) :: this
      type(DegenerateSubspace), intent(in) :: subspace
      real(dp),                 intent(in) :: thermal_energy
      type(AnharmonicData),     intent(in) :: anharmonic_data
      type(BasisStatesPointer)             :: output
    end function
    
    impure elemental function calculate_states_SubspaceBasis(this,subspace, &
       & subspace_potential,thermal_energy,state_energy_cutoff,             &
       & convergence_data,anharmonic_data) result(output) 
      import SubspaceBasis
      import DegenerateSubspace
      import PotentialBase
      import dp
      import ConvergenceData
      import AnharmonicData
      import BasisStatesPointer
      implicit none
      
      class(SubspaceBasis),     intent(in) :: this
      type(DegenerateSubspace), intent(in) :: subspace
      class(PotentialBase),     intent(in) :: subspace_potential
      real(dp),                 intent(in) :: thermal_energy
      real(dp),                 intent(in) :: state_energy_cutoff
      type(ConvergenceData),    intent(in) :: convergence_data
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
      
      class(SubspaceBasis),     intent(in)                   :: this
      class(BasisState),        intent(in),           target :: bra
      class(BasisState),        intent(in), optional, target :: ket
      type(DegenerateSubspace), intent(in)                   :: subspace
      type(AnharmonicData),     intent(in)                   :: anharmonic_data
      real(dp)                                               :: output
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
      
      class(SubspaceBasis),     intent(in)                   :: this
      class(BasisState),        intent(in),           target :: bra
      type(SparseMonomial),     intent(in)                   :: monomial
      class(BasisState),        intent(in), optional, target :: ket
      type(DegenerateSubspace), intent(in)                   :: subspace
      type(AnharmonicData),     intent(in)                   :: anharmonic_data
      complex(dp)                                            :: output
    end function
    
    impure elemental function kinetic_energy_SubspaceBasis(this,bra,ket, &
       & subspace,anharmonic_data) result(output)
      import SubspaceBasis
      import BasisState
      import DegenerateSubspace
      import AnharmonicData
      import dp
      implicit none
      
      class(SubspaceBasis),     intent(in)                   :: this
      class(BasisState),        intent(in),           target :: bra
      class(BasisState),        intent(in), optional, target :: ket
      type(DegenerateSubspace), intent(in)                   :: subspace
      type(AnharmonicData),     intent(in)                   :: anharmonic_data
      real(dp)                                               :: output
    end function
    
    impure elemental function harmonic_potential_energy_SubspaceBasis(this, &
       & bra,ket,subspace,anharmonic_data) result(output)
      import SubspaceBasis
      import BasisState
      import DegenerateSubspace
      import AnharmonicData
      import dp
      implicit none
      
      class(SubspaceBasis),     intent(in)                   :: this
      class(BasisState),        intent(in),           target :: bra
      class(BasisState),        intent(in), optional, target :: ket
      type(DegenerateSubspace), intent(in)                   :: subspace
      type(AnharmonicData),     intent(in)                   :: anharmonic_data
      real(dp)                                               :: output
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
      
      class(SubspaceBasis),     intent(in)                   :: this
      class(BasisState),        intent(in),           target :: bra
      class(BasisState),        intent(in), optional, target :: ket
      type(DegenerateSubspace), intent(in)                   :: subspace
      type(StressPrefactors),   intent(in)                   :: stress_prefactors
      type(AnharmonicData),     intent(in)                   :: anharmonic_data
      type(RealMatrix)                                       :: output
    end function
    
    impure elemental function thermodynamic_data_SubspaceBasis(this, &
       & thermal_energy,states,subspace,subspace_potential,          &
       & subspace_stress,stress_prefactors,anharmonic_data) result(output)
      import SubspaceBasis
      import dp
      import BasisStates
      import DegenerateSubspace
      import PotentialBase
      import StressBase
      import StressPrefactors
      import AnharmonicData
      import ThermodynamicData
      implicit none
      
      class(SubspaceBasis),     intent(in)                  :: this
      real(dp),                 intent(in)                  :: thermal_energy
      class(BasisStates),       intent(in),          target :: states
      type(DegenerateSubspace), intent(in)                  :: subspace
      class(PotentialBase),     intent(in)                  :: subspace_potential
      class(StressBase),        intent(in), optional        :: subspace_stress
      type(StressPrefactors),   intent(in), optional        :: stress_prefactors
      type(AnharmonicData),     intent(in)                  :: anharmonic_data
      type(ThermodynamicData)                               :: output
    end function
    
    impure elemental function wavefunctions_SubspaceBasis(this,states, &
       & subspace,anharmonic_data) result(output)
      import SubspaceBasis
      import BasisStates
      import DegenerateSubspace
      import AnharmonicData
      import SubspaceWavefunctionsPointer
      implicit none
      
      class(SubspaceBasis),     intent(in)         :: this
      class(BasisStates),       intent(in), target :: states
      type(DegenerateSubspace), intent(in)         :: subspace
      type(AnharmonicData),     intent(in)         :: anharmonic_data
      type(SubspaceWavefunctionsPointer)           :: output
    end function
    
    impure elemental function integrate_BasisStates_SubspaceBasis(this, &
       & states,monomial,subspace,anharmonic_data) result(output)
      import SubspaceBasis
      import BasisStates
      import SparseMonomial
      import DegenerateSubspace
      import QpointData
      import AnharmonicData
      import dp
      implicit none
      
      class(SubspaceBasis),     intent(in)         :: this
      class(BasisStates),       intent(in), target :: states
      type(SparseMonomial),     intent(in)         :: monomial
      type(DegenerateSubspace), intent(in)         :: subspace
      type(AnharmonicData),     intent(in)         :: anharmonic_data
      complex(dp)                                  :: output
    end function
    
    function free_energy_gradient_SubspaceBasis(this,subspace_potential,     &
       & basis_functions,subspace,states,thermal_energy,state_energy_cutoff, &
       & anharmonic_data) result(output)
      import SubspaceBasis
      import PotentialBase
      import DegenerateSubspace
      import BasisStates
      import dp
      import AnharmonicData
      implicit none
      
      class(SubspaceBasis),     intent(in) :: this
      class(PotentialBase),     intent(in) :: subspace_potential
      class(PotentialBase),     intent(in) :: basis_functions(:)
      type(DegenerateSubspace), intent(in) :: subspace
      class(BasisStates),       intent(in) :: states
      real(dp),                 intent(in) :: thermal_energy
      real(dp),                 intent(in) :: state_energy_cutoff
      type(AnharmonicData),     intent(in) :: anharmonic_data
      real(dp), allocatable                :: output(:)
    end function
    
    ! PotentialBase procedures.
    impure elemental function representation_PotentialBase() result(output)
      import String
      implicit none
      
      type(String) :: output
    end function
    
    impure elemental function energy_RealModeDisplacement_PotentialBase(this, &
       & displacement) result(output)
      import dp
      import PotentialBase
      import RealModeDisplacement
      implicit none
      
      class(PotentialBase),       intent(in) :: this
      type(RealModeDisplacement), intent(in) :: displacement
      real(dp)                               :: output
    end function
    
    impure elemental function energy_ComplexModeDisplacement_PotentialBase( &
       & this,displacement) result(output)
      import dp
      import PotentialBase
      import ComplexModeDisplacement
      implicit none
      
      class(PotentialBase),          intent(in) :: this
      type(ComplexModeDisplacement), intent(in) :: displacement
      complex(dp)                               :: output
    end function
    
    impure elemental function force_RealModeDisplacement_PotentialBase(this, &
       & displacement) result(output)
      import PotentialBase
      import RealModeDisplacement
      import RealModeForce
      implicit none
      
      class(PotentialBase),       intent(in) :: this
      type(RealModeDisplacement), intent(in) :: displacement
      type(RealModeForce)                    :: output
    end function
    
    impure elemental function force_ComplexModeDisplacement_PotentialBase( &
       & this,displacement) result(output)
      import PotentialBase
      import ComplexModeDisplacement
      import ComplexModeForce
      implicit none
      
      class(PotentialBase),          intent(in) :: this
      type(ComplexModeDisplacement), intent(in) :: displacement
      type(ComplexModeForce)                    :: output
    end function
    
    impure elemental subroutine braket_SubspaceBraKet_PotentialBase(this, &
       & braket,whole_subspace,anharmonic_data)
      import PotentialBase
      import SubspaceBraKet
      import AnharmonicData
      implicit none
      
      class(PotentialBase),  intent(inout)        :: this
      class(SubspaceBraKet), intent(in)           :: braket
      logical,               intent(in), optional :: whole_subspace
      type(AnharmonicData),  intent(in)           :: anharmonic_data
    end subroutine
    
    impure elemental subroutine braket_BasisState_PotentialBase(this,bra,ket, &
       & subspace,subspace_basis,whole_subspace,anharmonic_data)
      import PotentialBase
      import BasisState
      import DegenerateSubspace
      import SubspaceBasis
      import AnharmonicData
      import QpointData
      implicit none
      
      class(PotentialBase),     intent(inout)        :: this
      class(BasisState),        intent(in)           :: bra
      class(BasisState),        intent(in), optional :: ket
      type(DegenerateSubspace), intent(in)           :: subspace
      class(SubspaceBasis),     intent(in)           :: subspace_basis
      logical,                  intent(in), optional :: whole_subspace
      type(AnharmonicData),     intent(in)           :: anharmonic_data
    end subroutine
    
    impure elemental subroutine braket_BasisStates_PotentialBase(this,states, &
       & subspace,subspace_basis,whole_subspace,anharmonic_data)
      import PotentialBase
      import BasisStates
      import DegenerateSubspace
      import SubspaceBasis
      import AnharmonicData
      implicit none
      
      class(PotentialBase),     intent(inout)        :: this
      class(BasisStates),       intent(inout)        :: states
      type(DegenerateSubspace), intent(in)           :: subspace
      class(SubspaceBasis),     intent(in)           :: subspace_basis
      logical,                  intent(in), optional :: whole_subspace
      type(AnharmonicData),     intent(in)           :: anharmonic_data
    end subroutine
    
    impure elemental function harmonic_expectation_PotentialBase(this, &
       & frequency,thermal_energy,supercell_size,anharmonic_data)      &
       & result(output)
      import PotentialBase
      import dp
      import AnharmonicData
      implicit none
      
      class(PotentialBase), intent(in) :: this
      real(dp),             intent(in) :: frequency
      real(dp),             intent(in) :: thermal_energy
      integer,              intent(in) :: supercell_size
      type(AnharmonicData), intent(in) :: anharmonic_data
      real(dp)                         :: output
    end function
    
    impure elemental subroutine zero_energy_PotentialBase(this)
      import PotentialBase
      implicit none
      
      class(PotentialBase), intent(inout) :: this
    end subroutine
    
    impure elemental subroutine add_constant_PotentialBase(this,input)
      import PotentialBase
      import dp
      implicit none
      
      class(PotentialBase), intent(inout) :: this
      real(dp),             intent(in)    :: input
    end subroutine
    
    ! StressBase procedures.
    impure elemental function representation_StressBase() result(output)
      import String
      implicit none
      
      type(String) :: output
    end function
    
    impure elemental subroutine zero_stress_StressBase(this)
      import StressBase
      implicit none
      
      class(StressBase), intent(inout) :: this
    end subroutine
    
    impure elemental subroutine add_constant_StressBase(this,input)
      import StressBase
      import RealMatrix
      implicit none
      
      class(StressBase), intent(inout) :: this
      type(RealMatrix),  intent(in)    :: input
    end subroutine
    
    impure elemental function stress_RealModeDisplacement_StressBase(this, &
       & displacement) result(output)
      import StressBase
      import RealModeDisplacement
      import RealMatrix
      implicit none
      
      class(StressBase),          intent(in) :: this
      type(RealModeDisplacement), intent(in) :: displacement
      type(RealMatrix)                       :: output
    end function
    
    impure elemental function stress_ComplexModeDisplacement_StressBase( &
       & this,displacement) result(output)
      import StressBase
      import ComplexModeDisplacement
      import ComplexMatrix
      implicit none
      
      class(StressBase),             intent(in) :: this
      type(ComplexModeDisplacement), intent(in) :: displacement
      type(ComplexMatrix)                       :: output
    end function
    
    impure elemental subroutine braket_SubspaceBraKet_StressBase(this,braket, &
       & whole_subspace,anharmonic_data)
      import StressBase
      import SubspaceBraKet
      import AnharmonicData
      implicit none
      
      class(StressBase),     intent(inout)        :: this
      class(SubspaceBraKet), intent(in)           :: braket
      logical,               intent(in), optional :: whole_subspace
      type(AnharmonicData),  intent(in)           :: anharmonic_data
    end subroutine
    
    impure elemental subroutine braket_BasisState_StressBase(this,bra,ket, &
       & subspace,subspace_basis,whole_subspace,anharmonic_data)
      import StressBase
      import BasisState
      import DegenerateSubspace
      import SubspaceBasis
      import AnharmonicData
      implicit none
      
      class(StressBase),        intent(inout)        :: this
      class(BasisState),        intent(in)           :: bra
      class(BasisState),        intent(in), optional :: ket
      type(DegenerateSubspace), intent(in)           :: subspace
      class(SubspaceBasis),     intent(in)           :: subspace_basis
      logical,                  intent(in), optional :: whole_subspace
      type(AnharmonicData),     intent(in)           :: anharmonic_data
    end subroutine
    
    impure elemental subroutine braket_BasisStates_StressBase(this,states, &
       & subspace,subspace_basis,whole_subspace,anharmonic_data)
      import StressBase
      import BasisStates
      import DegenerateSubspace
      import SubspaceBasis
      import AnharmonicData
      implicit none
      
      class(StressBase),        intent(inout)        :: this
      class(BasisStates),       intent(inout)        :: states
      type(DegenerateSubspace), intent(in)           :: subspace
      class(SubspaceBasis),     intent(in)           :: subspace_basis
      logical,                  intent(in), optional :: whole_subspace
      type(AnharmonicData),     intent(in)           :: anharmonic_data
    end subroutine
    
    impure elemental function harmonic_expectation_StressBase(this,frequency, &
       & thermal_energy,supercell_size,anharmonic_data) result(output)
      import StressBase
      import dp
      import AnharmonicData
      import RealMatrix
      implicit none
      
      class(StressBase),    intent(in) :: this
      real(dp),             intent(in) :: frequency
      real(dp),             intent(in) :: thermal_energy
      integer,              intent(in) :: supercell_size
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
  
  interface PotentialBasePointer
    module procedure new_PotentialBasePointer
    module procedure new_PotentialBasePointer_Strings
    module procedure new_PotentialBasePointer_StringArray
  end interface
  
  interface StressBasePointer
    module procedure new_StressBasePointer
    module procedure new_StressBasePointer_Strings
    module procedure new_StressBasePointer_StringArray
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

subroutine startup_PotentialBase(this)
  implicit none
  
  class(PotentialBase), intent(in) :: this
  
  integer :: i
  
  if (.not.allocated(TYPES_PotentialBase)) then
    TYPES_PotentialBase = [PotentialBasePointer(this)]
  elseif (.not.any([(                                            &
     & this%representation()                                     &
     &    == TYPES_PotentialBase(i)%potential_%representation(), &
     & i=1,                                                      &
     & size(TYPES_PotentialBase)                                 )])) then
    TYPES_PotentialBase = [TYPES_PotentialBase, PotentialBasePointer(this)]
  endif
end subroutine

subroutine startup_StressBase(this)
  implicit none
  
  class(StressBase), intent(in) :: this
  
  integer :: i
  
  if (.not.allocated(TYPES_StressBase)) then
    TYPES_StressBase = [StressBasePointer(this)]
  elseif (.not.any([(                                      &
     & this%representation()                               &
     &    == TYPES_StressBase(i)%stress_%representation(), &
     & i=1,                                                &
     & size(TYPES_StressBase)                              )])) then
    TYPES_StressBase = [TYPES_StressBase, StressBasePointer(this)]
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
    this%frequency = this%basis_%frequency
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
impure elemental function initial_states_SubspaceBasisPointer(this,subspace, &
   & thermal_energy,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceBasisPointer), intent(in) :: this
  type(DegenerateSubspace),    intent(in) :: subspace
  real(dp),                    intent(in) :: thermal_energy
  type(AnharmonicData),        intent(in) :: anharmonic_data
  type(BasisStatesPointer)                :: output
  
  call this%check()
  
  output = this%basis_%initial_states( subspace,       &
                                     & thermal_energy, &
                                     & anharmonic_data )
end function

impure elemental function calculate_states_SubspaceBasisPointer(this, &
   & subspace,subspace_potential,thermal_energy,state_energy_cutoff,  &
   & convergence_data,anharmonic_data) result(output) 
  implicit none
  
  class(SubspaceBasisPointer), intent(in) :: this
  type(DegenerateSubspace),    intent(in) :: subspace
  class(PotentialBase),        intent(in) :: subspace_potential
  real(dp),                    intent(in) :: thermal_energy
  real(dp),                    intent(in) :: state_energy_cutoff
  type(ConvergenceData),       intent(in) :: convergence_data
  type(AnharmonicData),        intent(in) :: anharmonic_data
  type(BasisStatesPointer)                :: output
  
  call this%check()
  
  output = this%basis_%calculate_states( subspace,                &
                                       & subspace_potential,      &
                                       & thermal_energy,          &
                                       & state_energy_cutoff,     &
                                       & convergence_data,        &
                                       & anharmonic_data          )
end function

impure elemental subroutine process_subspace_potential_SubspaceBasisPointer( &
   & this,potential,states,subspace,anharmonic_data)
  implicit none
  
  class(SubspaceBasisPointer), intent(in)            :: this
  class(PotentialBase),        intent(inout)         :: potential
  class(BasisStates),          intent(inout), target :: states
  type(DegenerateSubspace),    intent(in)            :: subspace
  type(AnharmonicData),        intent(in)            :: anharmonic_data
  
  call this%check()
  
  call this%basis_%process_subspace_potential( potential,      &
                                             & states,         &
                                             & subspace,       &
                                             & anharmonic_data )
end subroutine

impure elemental subroutine process_subspace_stress_SubspaceBasisPointer( &
   & this,stress,states,subspace,anharmonic_data)
  implicit none
  
  class(SubspaceBasisPointer), intent(in)            :: this
  class(stressBase),           intent(inout)         :: stress
  class(BasisStates),          intent(inout), target :: states
  type(DegenerateSubspace),    intent(in)            :: subspace
  type(AnharmonicData),        intent(in)            :: anharmonic_data
  
  call this%check()
  
  call this%basis_%process_subspace_stress( stress,         &
                                          & states,         &
                                          & subspace,       &
                                          & anharmonic_data )
end subroutine

function select_symmetries_SubspaceBasisPointer(this,symmetries, &
   & anharmonic_data) result(output)
  implicit none
  
  class(SubspaceBasisPointer), intent(in) :: this
  type(DegenerateSymmetry),    intent(in) :: symmetries(:)
  type(AnharmonicData),        intent(in) :: anharmonic_data
  type(DegenerateSymmetry), allocatable   :: output(:)
  
  call this%check()
  
  output = this%basis_%select_symmetries(symmetries, anharmonic_data)
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
  
  class(SubspaceBasisPointer), intent(in)                   :: this
  class(BasisState),           intent(in),           target :: bra
  class(BasisState),           intent(in), optional, target :: ket
  type(DegenerateSubspace),    intent(in)                   :: subspace
  type(AnharmonicData),        intent(in)                   :: anharmonic_data
  real(dp)                                                  :: output
  
  call this%check()
  
  output = this%basis_%inner_product( bra,            &
                                    & ket,            &
                                    & subspace,       &
                                    & anharmonic_data )
end function

impure elemental function integrate_BasisState_SubspaceBasisPointer( &
   & this,bra,monomial,ket,subspace,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceBasisPointer), intent(in)                   :: this
  class(BasisState),           intent(in),           target :: bra
  type(SparseMonomial),        intent(in)                   :: monomial
  class(BasisState),           intent(in), optional, target :: ket
  type(DegenerateSubspace),    intent(in)                   :: subspace
  type(AnharmonicData),        intent(in)                   :: anharmonic_data
  complex(dp)                                               :: output
  
  call this%check()
  
  output = this%basis_%integrate( bra,            &
                                & monomial,       &
                                & ket,            &
                                & subspace,       &
                                & anharmonic_data )
end function

function free_energy_gradient_SubspaceBasisPointer(this,subspace_potential, &
   & basis_functions,subspace,states,thermal_energy,state_energy_cutoff,    &
   & anharmonic_data) result(output) 
  implicit none
  
  class(SubspaceBasisPointer), intent(in) :: this
  class(PotentialBase),        intent(in) :: subspace_potential
  class(PotentialBase),        intent(in) :: basis_functions(:)
  type(DegenerateSubspace),    intent(in) :: subspace
  class(BasisStates),          intent(in) :: states
  real(dp),                    intent(in) :: thermal_energy
  real(dp),                    intent(in) :: state_energy_cutoff
  type(AnharmonicData),        intent(in) :: anharmonic_data
  real(dp), allocatable                   :: output(:)
  
  call this%check()
  
  output = this%basis_%free_energy_gradient( subspace_potential,  &
                                           & basis_functions,     &
                                           & subspace,            &
                                           & states,              &
                                           & thermal_energy,      &
                                           & state_energy_cutoff, &
                                           & anharmonic_data      )
end function

impure elemental function kinetic_energy_SubspaceBasisPointer(this,bra,ket, &
   & subspace,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceBasisPointer), intent(in)                   :: this
  class(BasisState),           intent(in),           target :: bra
  class(BasisState),           intent(in), optional, target :: ket
  type(DegenerateSubspace),    intent(in)                   :: subspace
  type(AnharmonicData),        intent(in)                   :: anharmonic_data
  real(dp)                                                  :: output
  
  call this%check()
  
  output = this%basis_%kinetic_energy( bra,            &
                                     & ket,            &
                                     & subspace,       &
                                     & anharmonic_data )
end function

impure elemental function harmonic_potential_energy_SubspaceBasisPointer( &
   & this,bra,ket,subspace,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceBasisPointer), intent(in)                   :: this
  class(BasisState),           intent(in),           target :: bra
  class(BasisState),           intent(in), optional, target :: ket
  type(DegenerateSubspace),    intent(in)                   :: subspace
  type(AnharmonicData),        intent(in)                   :: anharmonic_data
  real(dp)                                                  :: output
  
  call this%check()
  
  output = this%basis_%harmonic_potential_energy( bra,            &
                                                & ket,            &
                                                & subspace,       &
                                                & anharmonic_data )
end function

impure elemental function kinetic_stress_SubspaceBasisPointer(this,bra,ket, &
   & subspace,stress_prefactors,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceBasisPointer), intent(in)                   :: this
  class(BasisState),           intent(in),           target :: bra
  class(BasisState),           intent(in), optional, target :: ket
  type(DegenerateSubspace),    intent(in)                   :: subspace
  type(StressPrefactors),      intent(in)                   :: stress_prefactors
  type(AnharmonicData),        intent(in)                   :: anharmonic_data
  type(RealMatrix)                                          :: output
  
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
  
  class(SubspaceBasisPointer), intent(in)                  :: this
  real(dp),                    intent(in)                  :: thermal_energy
  class(BasisStates),          intent(in),          target :: states
  type(DegenerateSubspace),    intent(in)                  :: subspace
  class(PotentialBase),        intent(in)                  :: subspace_potential
  class(StressBase),           intent(in), optional        :: subspace_stress
  type(StressPrefactors),      intent(in), optional        :: stress_prefactors
  type(AnharmonicData),        intent(in)                  :: anharmonic_data
  type(ThermodynamicData)                                  :: output
  
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
  
  class(SubspaceBasisPointer), intent(in)         :: this
  class(BasisStates),          intent(in), target :: states
  type(DegenerateSubspace),    intent(in)         :: subspace
  type(AnharmonicData),        intent(in)         :: anharmonic_data
  type(SubspaceWavefunctionsPointer)              :: output
  
  call this%check()
  
  output = this%basis_%wavefunctions( states,         &
                                    & subspace,       &
                                    & anharmonic_data )
end function

impure elemental function integrate_BasisStates_SubspaceBasisPointer( &
   & this,states,monomial,subspace,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceBasisPointer), intent(in)         :: this
  class(BasisStates),          intent(in), target :: states
  type(SparseMonomial),        intent(in)         :: monomial
  type(DegenerateSubspace),    intent(in)         :: subspace
  type(AnharmonicData),        intent(in)         :: anharmonic_data
  complex(dp)                                     :: output
  
  call this%check()
  
  output = this%basis_%integrate( states,         &
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
! PotentialBasePointer methods.
! ----------------------------------------------------------------------
! Construct a PotentialBasePointer from any type which extends PotentialBase.
impure elemental function new_PotentialBasePointer(potential) result(this)
  implicit none
  
  class(PotentialBase), intent(in) :: potential
  type(PotentialBasePointer)       :: this
  
  integer :: ialloc
  
  select type(potential); type is(PotentialBasePointer)
    this = potential
  class default
    this%representation_ = potential%representation()
    allocate( this%potential_, source=potential, &
            & stat=ialloc); call err(ialloc)
  end select
end function

! Checks that the pointer has been allocated before it is used.
impure elemental subroutine check_PotentialBasePointer(this)
  implicit none
  
  class(PotentialBasePointer), intent(in) :: this
  
  if (.not. allocated(this%potential_)) then
    call print_line(CODE_ERROR//': Trying to use a PotentialBasePointer &
       &before it has been allocated.')
    call err()
  endif
end subroutine

! Type representation.
impure elemental function representation_PotentialBasePointer() result(output)
  implicit none
  
  type(String) :: output
  
  output = 'pointer'
end function

! Return the stored potential.
impure function potential_PotentialBasePointer(this) result(output)
  implicit none
  
  class(PotentialBasePointer), intent(in) :: this
  class(PotentialBase), allocatable       :: output
  
  call this%check()
  
  output = this%potential_
end function

! Wrappers for all of PotentialBase's methods.
impure elemental subroutine zero_energy_PotentialBasePointer(this)
  implicit none
  
  class(PotentialBasePointer), intent(inout) :: this
  
  call this%check()
  
  call this%potential_%zero_energy()
end subroutine

impure elemental subroutine add_constant_PotentialBasePointer(this,input)
  implicit none
  
  class(PotentialBasePointer), intent(inout) :: this
  real(dp),                    intent(in)    :: input
  
  call this%check()
  
  call this%potential_%add_constant(input)
end subroutine

impure elemental function energy_RealModeDisplacement_PotentialBasePointer(this, &
   & displacement) result(output)
  implicit none
  
  class(PotentialBasePointer),    intent(in) :: this
  type(RealModeDisplacement), intent(in) :: displacement
  real(dp)                               :: output
  
  call this%check()
  
  output = this%potential_%energy(displacement)
end function

impure elemental function energy_ComplexModeDisplacement_PotentialBasePointer( &
   & this,displacement) result(output)
  implicit none
  
  class(PotentialBasePointer),       intent(in) :: this
  type(ComplexModeDisplacement), intent(in) :: displacement
  complex(dp)                               :: output
  
  call this%check()
  
  output = this%potential_%energy(displacement)
end function

impure elemental function force_RealModeDisplacement_PotentialBasePointer(this, &
   & displacement) result(output)
  implicit none
  
  class(PotentialBasePointer),    intent(in) :: this
  type(RealModeDisplacement), intent(in) :: displacement
  type(RealModeForce)                    :: output
  
  call this%check()
  
  output = this%potential_%force(displacement)
end function

impure elemental function force_ComplexModeDisplacement_PotentialBasePointer( &
   & this,displacement) result(output)
  implicit none
  
  class(PotentialBasePointer),       intent(in) :: this
  type(ComplexModeDisplacement), intent(in) :: displacement
  type(ComplexModeForce)                    :: output
  
  call this%check()
  
  output = this%potential_%force(displacement)
end function

impure elemental subroutine braket_SubspaceBraKet_PotentialBasePointer(this, &
   & braket,whole_subspace,anharmonic_data) 
  implicit none
  
  class(PotentialBasePointer), intent(inout)        :: this
  class(SubspaceBraKet),       intent(in)           :: braket
  logical,                     intent(in), optional :: whole_subspace
  type(AnharmonicData),        intent(in)           :: anharmonic_data
  
  call this%check()
  
  call this%potential_%braket(braket,whole_subspace,anharmonic_data)
end subroutine

impure elemental subroutine braket_BasisState_PotentialBasePointer(this,bra, &
   & ket,subspace,subspace_basis,whole_subspace,anharmonic_data) 
  implicit none
  
  class(PotentialBasePointer), intent(inout)        :: this
  class(BasisState),           intent(in)           :: bra
  class(BasisState),           intent(in), optional :: ket
  type(DegenerateSubspace),    intent(in)           :: subspace
  class(SubspaceBasis),        intent(in)           :: subspace_basis
  logical,                     intent(in), optional :: whole_subspace
  type(AnharmonicData),        intent(in)           :: anharmonic_data
  
  call this%check()
  
  call this%potential_%braket( bra,            &
                             & ket,            &
                             & subspace,       &
                             & subspace_basis, &
                             & whole_subspace, &
                             & anharmonic_data )
end subroutine

impure elemental subroutine braket_BasisStates_PotentialBasePointer(this, &
   & states,subspace,subspace_basis,whole_subspace,anharmonic_data) 
  implicit none
  
  class(PotentialBasePointer), intent(inout)        :: this
  class(BasisStates),          intent(inout)        :: states
  type(DegenerateSubspace),    intent(in)           :: subspace
  class(SubspaceBasis),        intent(in)           :: subspace_basis
  logical,                     intent(in), optional :: whole_subspace
  type(AnharmonicData),        intent(in)           :: anharmonic_data
  
  call this%check()
  
  call this%potential_%braket( states,         &
                             & subspace,       &
                             & subspace_basis, &
                             & whole_subspace, &
                             & anharmonic_data )
end subroutine

impure elemental function harmonic_expectation_PotentialBasePointer(this, &
   & frequency,thermal_energy,supercell_size,anharmonic_data) result(output)
  implicit none
  
  class(PotentialBasePointer), intent(in) :: this
  real(dp),                intent(in) :: frequency
  real(dp),                intent(in) :: thermal_energy
  integer,                 intent(in) :: supercell_size
  type(AnharmonicData),    intent(in) :: anharmonic_data
  real(dp)                            :: output
  
  call this%check()
  
  output = this%potential_%harmonic_expectation( frequency,      &
                                               & thermal_energy, &
                                               & supercell_size, &
                                               & anharmonic_data )
end function

recursive function potential_energy_SubspaceBraKet_PotentialBasePointer(this, &
   & braket,anharmonic_data) result(output) 
  implicit none
  
  class(PotentialBasePointer), intent(in) :: this
  class(SubspaceBraKet),       intent(in) :: braket
  type(AnharmonicData),        intent(in) :: anharmonic_data
  real(dp)                                :: output
  
  call this%check()
  
  output = this%potential_%potential_energy(braket, anharmonic_data)
end function

recursive function potential_energy_BasisState_PotentialBasePointer(this, &
   & bra,ket,subspace,subspace_basis,anharmonic_data) result(output) 
  implicit none
  
  class(PotentialBasePointer), intent(in)           :: this
  class(BasisState),           intent(in)           :: bra
  class(BasisState),           intent(in), optional :: ket
  type(DegenerateSubspace),    intent(in)           :: subspace
  class(SubspaceBasis),        intent(in)           :: subspace_basis
  type(AnharmonicData),        intent(in)           :: anharmonic_data
  real(dp)                                          :: output
  
  call this%check()
  
  output = this%potential_%potential_energy( bra,            &
                                           & ket,            &
                                           & subspace,       &
                                           & subspace_basis, &
                                           & anharmonic_data )
end function

! I/O.
subroutine read_PotentialBasePointer(this,input)
  implicit none
  
  class(PotentialBasePointer), intent(out) :: this
  type(String),            intent(in)  :: input(:)
  
  type(String), allocatable :: line(:)
  
  type(String) :: representation
  
  integer :: i
  
  select type(this); type is(PotentialBasePointer)
    line = split_line(input(1))
    representation = line(3)
    
    ! Identify which type corresponds to the representation.
    i = first([(                                                             &
       & TYPES_PotentialBase(i)%potential_%representation()==representation, &
       & i=1,                                                                &
       & size(TYPES_PotentialBase)                                          )])
    
    ! Read the input into the element of the correct type,
    !    and copy that element into the output.
    call TYPES_PotentialBase(i)%potential_%read(input(2:))
    this = PotentialBasePointer(TYPES_PotentialBase(i))
  class default
    call err()
  end select
end subroutine

function write_PotentialBasePointer(this) result(output)
  implicit none
  
  class(PotentialBasePointer), intent(in) :: this
  type(String), allocatable               :: output(:)
  
  select type(this); type is(PotentialBasePointer)
    output = [ 'Potential representation: '//this%representation_, &
             & str(this%potential_)                                ]
  end select
end function

function new_PotentialBasePointer_Strings(input) result(this)
  implicit none
  
  type(String), intent(in)   :: input(:)
  type(PotentialBasePointer) :: this
  
  call this%read(input)
end function

impure elemental function new_PotentialBasePointer_StringArray(input) &
   & result(this) 
  implicit none
  
  type(StringArray), intent(in) :: input
  type(PotentialBasePointer)        :: this
  
  this = PotentialBasePointer(str(input))
end function

! ----------------------------------------------------------------------
! StressBasePointer methods.
! ----------------------------------------------------------------------
! Construct a StressBasePointer from any type which extends StressBase.
impure elemental function new_StressBasePointer(stress) result(this)
  implicit none
  
  class(StressBase), intent(in) :: stress
  type(StressBasePointer)       :: this
  
  integer :: ialloc
  
  select type(stress); type is(StressBasePointer)
    this = stress
  class default
    this%representation_ = stress%representation()
    allocate( this%stress_, source=stress, &
            & stat=ialloc); call err(ialloc)
  end select
end function

! Checks that the pointer has been allocated before it is used.
impure elemental subroutine check_StressBasePointer(this)
  implicit none
  
  class(StressBasePointer), intent(in) :: this
  
  if (.not. allocated(this%stress_)) then
    call print_line(CODE_ERROR//': Trying to use a StressBasePointer before &
       &it has been allocated.')
    call err()
  endif
end subroutine

! Type representation.
impure elemental function representation_StressBasePointer() result(output)
  implicit none
  
  type(String) :: output
  
  output = 'pointer'
end function

impure elemental subroutine zero_stress_StressBasePointer(this)
  implicit none
  
  class(StressBasePointer), intent(inout) :: this
  
  call this%check()
  
  call this%stress_%zero_stress()
end subroutine

impure elemental subroutine add_constant_StressBasePointer(this,input)
  implicit none
  
  class(StressBasePointer), intent(inout) :: this
  type(RealMatrix),         intent(in)    :: input
  
  call this%check()
  
  call this%stress_%add_constant(input)
end subroutine

impure elemental function stress_RealModeDisplacement_StressBasePointer(this, &
   & displacement) result(output)
  implicit none
  
  class(StressBasePointer),       intent(in) :: this
  type(RealModeDisplacement), intent(in) :: displacement
  type(RealMatrix)                       :: output
  
  call this%check()
  
  output = this%stress_%stress(displacement)
end function

impure elemental function stress_ComplexModeDisplacement_StressBasePointer( &
   & this,displacement) result(output)
  implicit none
  
  class(StressBasePointer),          intent(in) :: this
  type(ComplexModeDisplacement), intent(in) :: displacement
  type(ComplexMatrix)                       :: output
  
  call this%check()
  
  output = this%stress_%stress(displacement)
end function

impure elemental subroutine braket_SubspaceBraKet_StressBasePointer(this, &
   & braket,whole_subspace,anharmonic_data) 
  implicit none
  
  class(StressBasePointer), intent(inout)        :: this
  class(SubspaceBraKet),    intent(in)           :: braket
  logical,                  intent(in), optional :: whole_subspace
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  
  call this%check()
  
  call this%stress_%braket(braket,whole_subspace,anharmonic_data)
end subroutine

impure elemental subroutine braket_BasisState_StressBasePointer(this,bra, &
   & ket,subspace,subspace_basis,whole_subspace,anharmonic_data) 
  implicit none
  
  class(StressBasePointer), intent(inout)        :: this
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

impure elemental subroutine braket_BasisStates_StressBasePointer(this, &
   & states,subspace,subspace_basis,whole_subspace,anharmonic_data) 
  implicit none
  
  class(StressBasePointer), intent(inout)        :: this
  class(BasisStates),       intent(inout)        :: states
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  logical,                  intent(in), optional :: whole_subspace
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  
  call this%check()
  
  call this%stress_%braket( states,         &
                          & subspace,       &
                          & subspace_basis, &
                          & whole_subspace, &
                          & anharmonic_data )
end subroutine

impure elemental function harmonic_expectation_StressBasePointer(this,frequency, &
   & thermal_energy,supercell_size,anharmonic_data) result(output)
  implicit none
  
  class(StressBasePointer), intent(in) :: this
  real(dp),             intent(in) :: frequency
  real(dp),             intent(in) :: thermal_energy
  integer,              intent(in) :: supercell_size
  type(AnharmonicData), intent(in) :: anharmonic_data
  type(RealMatrix)                 :: output
  
  call this%check()
  
  output = this%stress_%harmonic_expectation( frequency,      &
                                            & thermal_energy, &
                                            & supercell_size, &
                                            & anharmonic_data )
end function

recursive function potential_stress_SubspaceBraKet_StressBasePointer(this,braket, &
   & anharmonic_data) result(output) 
  implicit none
  
  class(StressBasePointer),  intent(in) :: this
  class(SubspaceBraKet), intent(in) :: braket
  type(AnharmonicData),  intent(in) :: anharmonic_data
  type(RealMatrix)                  :: output
  
  call this%check()
  
  output = this%stress_%potential_stress(braket, anharmonic_data)
end function

recursive function potential_stress_BasisState_StressBasePointer(this,bra,ket, &
   & subspace,subspace_basis,anharmonic_data) result(output) 
  implicit none
  
  class(StressBasePointer),     intent(in)           :: this
  class(BasisState),        intent(in)           :: bra
  class(BasisState),        intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(RealMatrix)                               :: output
  
  call this%check()
  
  output = this%stress_%potential_stress( bra,            &
                                        & ket,            &
                                        & subspace,       &
                                        & subspace_basis, &
                                        & anharmonic_data )
end function

! I/O.
subroutine read_StressBasePointer(this,input)
  implicit none
  
  class(StressBasePointer), intent(out) :: this
  type(String),             intent(in)  :: input(:)
  
  type(String), allocatable :: line(:)
  
  type(String) :: representation
  
  integer :: i
  
  select type(this); type is(StressBasePointer)
    line = split_line(input(1))
    representation = line(3)
    
    ! Identify which type corresponds to the representation.
    i = first([(                                                       &
       & TYPES_StressBase(i)%stress_%representation()==representation, &
       & i=1,                                                          &
       & size(TYPES_StressBase)                                        )])
    
    ! Read the input into the element of the correct type,
    !    and copy that element into the output.
    call TYPES_StressBase(i)%stress_%read(input(2:))
    this = StressBasePointer(TYPES_StressBase(i))
  class default
    call err()
  end select
end subroutine

function write_StressBasePointer(this) result(output)
  implicit none
  
  class(StressBasePointer), intent(in) :: this
  type(String), allocatable            :: output(:)
  
  select type(this); type is(StressBasePointer)
    output = [ 'Stress representation: '//this%representation_, &
             & str(this%stress_)                                ]
  end select
end function

function new_StressBasePointer_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(StressBasePointer)  :: this
  
  call this%read(input)
end function

impure elemental function new_StressBasePointer_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(StressBasePointer)       :: this
  
  this = StressBasePointer(str(input))
end function

! ----------------------------------------------------------------------
! Concrete methods of abstract classes.
! ----------------------------------------------------------------------
! SubspaceBasis methods.
impure elemental subroutine process_subspace_potential_SubspaceBasis(this, &
   & potential,states,subspace,anharmonic_data)
  implicit none
  
  class(SubspaceBasis),     intent(in)            :: this
  class(PotentialBase),     intent(inout)         :: potential
  class(BasisStates),       intent(inout), target :: states
  type(DegenerateSubspace), intent(in)            :: subspace
  type(AnharmonicData),     intent(in)            :: anharmonic_data
  
  ! This does nothing by default.
end subroutine

impure elemental subroutine process_subspace_stress_SubspaceBasis(this,stress, &
   & states,subspace,anharmonic_data)
  implicit none
  
  class(SubspaceBasis),     intent(in)            :: this
  class(stressBase),        intent(inout)         :: stress
  class(BasisStates),       intent(inout), target :: states
  type(DegenerateSubspace), intent(in)            :: subspace
  type(AnharmonicData),     intent(in)            :: anharmonic_data
  
  ! This does nothing by default.
end subroutine

function select_symmetries_SubspaceBasis(this,symmetries,anharmonic_data) &
   & result(output) 
  implicit none
  
  class(SubspaceBasis),     intent(in)  :: this
  type(DegenerateSymmetry), intent(in)  :: symmetries(:)
  type(AnharmonicData),     intent(in)  :: anharmonic_data
  type(DegenerateSymmetry), allocatable :: output(:)
  
  ! By default, this just returns the whole list.
  output = symmetries
end function

! Concrete PotentialBase methods.
impure elemental function undisplaced_energy(this) result(output)
  implicit none
  
  class(PotentialBase), intent(in) :: this
  real(dp)                         :: output
  
  type(RealSingleDisplacement) :: zero_displacement(0)
  
  output = this%energy(RealModeDisplacement(zero_displacement))
end function

recursive function potential_energy_SubspaceBraKet_PotentialBase(this,braket, &
   & anharmonic_data) result(output) 
  implicit none
  
  class(PotentialBase),  intent(in) :: this
  class(SubspaceBraKet), intent(in) :: braket
  type(AnharmonicData),  intent(in) :: anharmonic_data
  real(dp)                          :: output
  
  type(PotentialBasePointer), allocatable :: integrated_potential
  
  integrated_potential = PotentialBasePointer(this)
  call integrated_potential%braket( braket          = braket,         &
                                  & anharmonic_data = anharmonic_data )
  output = integrated_potential%undisplaced_energy()
end function

recursive function potential_energy_BasisState_PotentialBase(this,bra,ket, &
   & subspace,subspace_basis,anharmonic_data) result(output) 
  implicit none
  
  class(PotentialBase),     intent(in)           :: this
  class(BasisState),        intent(in)           :: bra
  class(BasisState),        intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  real(dp)                                       :: output
  
  type(PotentialBasePointer), allocatable :: integrated_potential
  
  integrated_potential = PotentialBasePointer(this)
  call integrated_potential%braket( bra             = bra,            &
                                  & ket             = ket,            &
                                  & subspace        = subspace,       &
                                  & subspace_basis  = subspace_basis, &
                                  & anharmonic_data = anharmonic_data )
  output = integrated_potential%undisplaced_energy()
end function

! StressBase methods.
impure elemental function undisplaced_stress(this) result(output)
  implicit none
  
  class(StressBase), intent(in) :: this
  type(RealMatrix)              :: output
  
  type(RealSingleDisplacement) :: zero_displacement(0)
  
  output = this%stress(RealModeDisplacement(zero_displacement))
end function

recursive function potential_stress_SubspaceBraKet_StressBase(this,braket, &
   & anharmonic_data) result(output) 
  implicit none
  
  class(StressBase),     intent(in) :: this
  class(SubspaceBraKet), intent(in) :: braket
  type(AnharmonicData),  intent(in) :: anharmonic_data
  type(RealMatrix)                  :: output
  
  type(StressBasePointer), allocatable :: integrated_stress
  
  integrated_stress = StressBasePointer(this)
  call integrated_stress%braket( braket          = braket,         &
                               & anharmonic_data = anharmonic_data )
  output = integrated_stress%undisplaced_stress()
end function

recursive function potential_stress_BasisState_StressBase(this,bra,ket, &
   & subspace,subspace_basis,anharmonic_data) result(output) 
  implicit none
  
  class(StressBase),        intent(in)           :: this
  class(BasisState),        intent(in)           :: bra
  class(BasisState),        intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(RealMatrix)                               :: output
  
  type(StressBasePointer), allocatable :: integrated_stress
  
  integrated_stress = StressBasePointer(this)
  call integrated_stress%braket( bra             = bra,            &
                               & ket             = ket,            &
                               & subspace        = subspace,       &
                               & subspace_basis  = subspace_basis, &
                               & anharmonic_data = anharmonic_data )
  output = integrated_stress%undisplaced_stress()
end function
end module
