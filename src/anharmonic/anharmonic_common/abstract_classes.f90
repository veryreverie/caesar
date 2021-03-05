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
  
  type, abstract, extends(Stringsable) :: PotentialBase
  contains
    procedure(representation_PotentialBase), public, deferred, nopass :: &
       & representation
    
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
  
  type, abstract, extends(Stringsable) :: StressBase
  contains
    procedure(representation_StressBase), public, deferred, nopass :: &
       & representation
    
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
  
  interface
    module function types_SubspaceBasis() result(output)
      type(SubspaceBasisPointer), allocatable :: output(:)
    end function
    
    module function types_PotentialBase() result(output)
      type(PotentialBasePointer), allocatable :: output(:)
    end function
    
    module function types_StressBase() result(output)
      type(StressBasePointer), allocatable :: output(:)
    end function
  end interface
  
  interface SubspaceBasisPointer
    ! ----------------------------------------------------------------------
    ! SubspaceBasisPointer methods.
    ! ----------------------------------------------------------------------
    ! Construct a SubspaceBasisPointer from any type which extends SubspaceBasis.
    impure elemental module function new_SubspaceBasisPointer(basis) &
       & result(this) 
      class(SubspaceBasis), intent(in) :: basis
      type(SubspaceBasisPointer)       :: this
    end function
  end interface
  
  interface
    ! Checks that the pointer has been allocated before it is used.
    module subroutine check_SubspaceBasisPointer(this) 
      class(SubspaceBasisPointer), intent(in) :: this
    end subroutine
  end interface
  
  interface
    ! Type representation.
    impure elemental module function representation_SubspaceBasisPointer() &
       & result(output) 
      type(String) :: output
    end function
  end interface
  
  interface
    ! SubspaceBasis methods.
    impure elemental module function initial_states_SubspaceBasisPointer(this,subspace,thermal_energy,anharmonic_data) result(output) 
      class(SubspaceBasisPointer), intent(in) :: this
      type(DegenerateSubspace),    intent(in) :: subspace
      real(dp),                    intent(in) :: thermal_energy
      type(AnharmonicData),        intent(in) :: anharmonic_data
      type(BasisStatesPointer)                :: output
    end function
  end interface
  
  interface
    impure elemental module function calculate_states_SubspaceBasisPointer(this,subspace,subspace_potential,thermal_energy,state_energy_cutoff,convergence_data,anharmonic_data) result(output) 
      class(SubspaceBasisPointer), intent(in) :: this
      type(DegenerateSubspace),    intent(in) :: subspace
      class(PotentialBase),        intent(in) :: subspace_potential
      real(dp),                    intent(in) :: thermal_energy
      real(dp),                    intent(in) :: state_energy_cutoff
      type(ConvergenceData),       intent(in) :: convergence_data
      type(AnharmonicData),        intent(in) :: anharmonic_data
      type(BasisStatesPointer)                :: output
    end function
  end interface
  
  interface
    impure elemental module subroutine process_subspace_potential_SubspaceBasisPointer(   this,potential,states,subspace,anharmonic_data) 
      class(SubspaceBasisPointer), intent(in)            :: this
      class(PotentialBase),        intent(inout)         :: potential
      class(BasisStates),          intent(inout), target :: states
      type(DegenerateSubspace),    intent(in)            :: subspace
      type(AnharmonicData),        intent(in)            :: anharmonic_data
    end subroutine
  end interface
  
  interface
    impure elemental module subroutine process_subspace_stress_SubspaceBasisPointer(   this,stress,states,subspace,anharmonic_data) 
      class(SubspaceBasisPointer), intent(in)            :: this
      class(stressBase),           intent(inout)         :: stress
      class(BasisStates),          intent(inout), target :: states
      type(DegenerateSubspace),    intent(in)            :: subspace
      type(AnharmonicData),        intent(in)            :: anharmonic_data
    end subroutine
  end interface
  
  interface
    module function select_symmetries_SubspaceBasisPointer(this,symmetries, &
       & anharmonic_data) result(output) 
      class(SubspaceBasisPointer), intent(in) :: this
      type(DegenerateSymmetry),    intent(in) :: symmetries(:)
      type(AnharmonicData),        intent(in) :: anharmonic_data
      type(DegenerateSymmetry), allocatable   :: output(:)
    end function
  end interface
  
  interface
    module function mode_ids_SubspaceBasisPointer(this,subspace, &
       & anharmonic_data) result(output) 
      class(SubspaceBasisPointer), intent(in) :: this
      type(DegenerateSubspace),    intent(in) :: subspace
      type(AnharmonicData),        intent(in) :: anharmonic_data
      integer, allocatable                    :: output(:)
    end function
  end interface
  
  interface
    module function paired_mode_ids_SubspaceBasisPointer(this,subspace, &
       & anharmonic_data) result(output) 
      class(SubspaceBasisPointer), intent(in) :: this
      type(DegenerateSubspace),    intent(in) :: subspace
      type(AnharmonicData),        intent(in) :: anharmonic_data
      integer, allocatable                    :: output(:)
    end function
  end interface
  
  interface
    impure elemental module function inner_product_SubspaceBasisPointer(this,bra,ket,subspace,anharmonic_data) result(output) 
      class(SubspaceBasisPointer), intent(in)                   :: this
      class(BasisState),           intent(in),           target :: bra
      class(BasisState),           intent(in), optional, target :: ket
      type(DegenerateSubspace),    intent(in)                   :: subspace
      type(AnharmonicData),        intent(in)                   :: anharmonic_data
      real(dp)                                                  :: output
    end function
  end interface
  
  interface
    impure elemental module function integrate_BasisState_SubspaceBasisPointer(   this,bra,monomial,ket,subspace,anharmonic_data) result(output) 
      class(SubspaceBasisPointer), intent(in)                   :: this
      class(BasisState),           intent(in),           target :: bra
      type(SparseMonomial),        intent(in)                   :: monomial
      class(BasisState),           intent(in), optional, target :: ket
      type(DegenerateSubspace),    intent(in)                   :: subspace
      type(AnharmonicData),        intent(in)                   :: anharmonic_data
      complex(dp)                                               :: output
    end function
  end interface
  
  interface
    module function free_energy_gradient_SubspaceBasisPointer(this,         &
        subspace_potential,basis_functions,subspace,states,thermal_energy, &
          & state_energy_cutoff,anharmonic_data) result(output) 
      class(SubspaceBasisPointer), intent(in) :: this
      class(PotentialBase),        intent(in) :: subspace_potential
      class(PotentialBase) ,intent(in) :: basis_functions(:) 
      type(DegenerateSubspace),    intent(in) :: subspace
      class(BasisStates),          intent(in) :: states
      real(dp),                    intent(in) :: thermal_energy
      real(dp),                    intent(in) :: state_energy_cutoff
      type(AnharmonicData),        intent(in) :: anharmonic_data
      real(dp), allocatable                   :: output(:)
    end function
  end interface
  
  interface
    impure elemental module function kinetic_energy_SubspaceBasisPointer(this,bra,ket,subspace,anharmonic_data) result(output) 
      class(SubspaceBasisPointer), intent(in)                   :: this
      class(BasisState),           intent(in),           target :: bra
      class(BasisState),           intent(in), optional, target :: ket
      type(DegenerateSubspace),    intent(in)                   :: subspace
      type(AnharmonicData),        intent(in)                   :: anharmonic_data
      real(dp)                                                  :: output
    end function
  end interface
  
  interface
    impure elemental module function harmonic_potential_energy_SubspaceBasisPointer(   this,bra,ket,subspace,anharmonic_data) result(output) 
      class(SubspaceBasisPointer), intent(in)                   :: this
      class(BasisState),           intent(in),           target :: bra
      class(BasisState),           intent(in), optional, target :: ket
      type(DegenerateSubspace),    intent(in)                   :: subspace
      type(AnharmonicData),        intent(in)                   :: anharmonic_data
      real(dp)                                                  :: output
    end function
  end interface
  
  interface
    impure elemental module function kinetic_stress_SubspaceBasisPointer(this,bra,ket,subspace,stress_prefactors,anharmonic_data) result(output) 
      class(SubspaceBasisPointer), intent(in)                   :: this
      class(BasisState),           intent(in),           target :: bra
      class(BasisState),           intent(in), optional, target :: ket
      type(DegenerateSubspace),    intent(in)                   :: subspace
      type(StressPrefactors),      intent(in)                   :: stress_prefactors
      type(AnharmonicData),        intent(in)                   :: anharmonic_data
      type(RealMatrix)                                          :: output
    end function
  end interface
  
  interface
    impure elemental module function thermodynamic_data_SubspaceBasisPointer(this,thermal_energy,states,subspace,subspace_potential,subspace_stress,stress_prefactors,anharmonic_data) result(output) 
      class(SubspaceBasisPointer), intent(in)                  :: this
      real(dp),                    intent(in)                  :: thermal_energy
      class(BasisStates),          intent(in),          target :: states
      type(DegenerateSubspace),    intent(in)                  :: subspace
      class(PotentialBase),        intent(in)                  :: subspace_potential
      class(StressBase),           intent(in), optional        :: subspace_stress
      type(StressPrefactors),      intent(in), optional        :: stress_prefactors
      type(AnharmonicData),        intent(in)                  :: anharmonic_data
      type(ThermodynamicData)                                  :: output
    end function
  end interface
  
  interface
    impure elemental module function wavefunctions_SubspaceBasisPointer(this,states,subspace,anharmonic_data) result(output) 
      class(SubspaceBasisPointer), intent(in)         :: this
      class(BasisStates),          intent(in), target :: states
      type(DegenerateSubspace),    intent(in)         :: subspace
      type(AnharmonicData),        intent(in)         :: anharmonic_data
      type(SubspaceWavefunctionsPointer) :: output
    end function
  end interface
  
  interface
    impure elemental module function integrate_BasisStates_SubspaceBasisPointer(   this,states,monomial,subspace,anharmonic_data) result(output) 
      class(SubspaceBasisPointer), intent(in)         :: this
      class(BasisStates),          intent(in), target :: states
      type(SparseMonomial),        intent(in)         :: monomial
      type(DegenerateSubspace),    intent(in)         :: subspace
      type(AnharmonicData),        intent(in)         :: anharmonic_data
      complex(dp)                                     :: output
    end function
  end interface
  
  interface
    ! I/O.
    module subroutine read_SubspaceBasisPointer(this,input) 
      class(SubspaceBasisPointer), intent(out) :: this
      type(String),                intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_SubspaceBasisPointer(this) result(output) 
      class(SubspaceBasisPointer), intent(in) :: this
      type(String), allocatable               :: output(:)
    end function
  end interface
  
  interface SubspaceBasisPointer
    module function new_SubspaceBasisPointer_Strings(input) result(this) 
      type(String), intent(in)   :: input(:)
      type(SubspaceBasisPointer) :: this
    end function
  
    impure elemental module function new_SubspaceBasisPointer_StringArray(input) result(this) 
      type(StringArray), intent(in) :: input
      type(SubspaceBasisPointer)    :: this
    end function
  end interface
  
  interface PotentialBasePointer
    ! ----------------------------------------------------------------------
    ! PotentialBasePointer methods.
    ! ----------------------------------------------------------------------
    ! Construct a PotentialBasePointer from any type which extends PotentialBase.
    impure elemental module function new_PotentialBasePointer(potential) &
       & result(this) 
      class(PotentialBase), intent(in) :: potential
      type(PotentialBasePointer)       :: this
    end function
  end interface
  
  interface
    ! Checks that the pointer has been allocated before it is used.
    impure elemental module subroutine check_PotentialBasePointer(this) 
      class(PotentialBasePointer), intent(in) :: this
    end subroutine
  end interface
  
  interface
    ! Type representation.
    impure elemental module function representation_PotentialBasePointer() &
       & result(output) 
      type(String) :: output
    end function
  end interface
  
  interface
    ! Return the stored potential.
    impure module function potential_PotentialBasePointer(this) result(output) 
      class(PotentialBasePointer), intent(in) :: this
      class(PotentialBase), allocatable       :: output
    end function
  end interface
  
  interface
    ! Wrappers for all of PotentialBase's methods.
    impure elemental module subroutine zero_energy_PotentialBasePointer(this) 
      class(PotentialBasePointer), intent(inout) :: this
    end subroutine
  end interface
  
  interface
    impure elemental module subroutine add_constant_PotentialBasePointer(this,input) 
      class(PotentialBasePointer), intent(inout) :: this
      real(dp),                    intent(in)    :: input
    end subroutine
  end interface
  
  interface
    impure elemental module function energy_RealModeDisplacement_PotentialBasePointer(this,displacement) result(output) 
      class(PotentialBasePointer),    intent(in) :: this
      type(RealModeDisplacement), intent(in) :: displacement
      real(dp)                               :: output
    end function
  end interface
  
  interface
    impure elemental module function energy_ComplexModeDisplacement_PotentialBasePointer(   this,displacement) result(output) 
      class(PotentialBasePointer),       intent(in) :: this
      type(ComplexModeDisplacement), intent(in) :: displacement
      complex(dp)                               :: output
    end function
  end interface
  
  interface
    impure elemental module function force_RealModeDisplacement_PotentialBasePointer(this,displacement) result(output) 
      class(PotentialBasePointer),    intent(in) :: this
      type(RealModeDisplacement), intent(in) :: displacement
      type(RealModeForce)                    :: output
    end function
  end interface
  
  interface
    impure elemental module function force_ComplexModeDisplacement_PotentialBasePointer(   this,displacement) result(output) 
      class(PotentialBasePointer),       intent(in) :: this
      type(ComplexModeDisplacement), intent(in) :: displacement
      type(ComplexModeForce)                    :: output
    end function
  end interface
  
  interface
    impure elemental module subroutine braket_SubspaceBraKet_PotentialBasePointer(this,braket,whole_subspace,anharmonic_data) 
      class(PotentialBasePointer), intent(inout)        :: this
      class(SubspaceBraKet),       intent(in)           :: braket
      logical,                     intent(in), optional :: whole_subspace
      type(AnharmonicData),        intent(in)           :: anharmonic_data
    end subroutine
  end interface
  
  interface
    impure elemental module subroutine braket_BasisState_PotentialBasePointer(this,bra,ket,subspace,subspace_basis,whole_subspace,anharmonic_data) 
      class(PotentialBasePointer), intent(inout)        :: this
      class(BasisState),           intent(in)           :: bra
      class(BasisState),           intent(in), optional :: ket
      type(DegenerateSubspace),    intent(in)           :: subspace
      class(SubspaceBasis),        intent(in)           :: subspace_basis
      logical,                     intent(in), optional :: whole_subspace
      type(AnharmonicData),        intent(in)           :: anharmonic_data
    end subroutine
  end interface
  
  interface
    impure elemental module subroutine braket_BasisStates_PotentialBasePointer(this,states,subspace,subspace_basis,whole_subspace,anharmonic_data) 
      class(PotentialBasePointer), intent(inout)        :: this
      class(BasisStates),          intent(inout)        :: states
      type(DegenerateSubspace),    intent(in)           :: subspace
      class(SubspaceBasis),        intent(in)           :: subspace_basis
      logical,                     intent(in), optional :: whole_subspace
      type(AnharmonicData),        intent(in)           :: anharmonic_data
    end subroutine
  end interface
  
  interface
    impure elemental module function harmonic_expectation_PotentialBasePointer(this,frequency,thermal_energy,supercell_size,anharmonic_data) result(output) 
      class(PotentialBasePointer), intent(in) :: this
      real(dp),                intent(in) :: frequency
      real(dp),                intent(in) :: thermal_energy
      integer,                 intent(in) :: supercell_size
      type(AnharmonicData),    intent(in) :: anharmonic_data
      real(dp)                            :: output
    end function
  end interface
  
  interface
    recursive module function potential_energy_SubspaceBraKet_PotentialBasePointer(this,braket,anharmonic_data) result(output) 
      class(PotentialBasePointer), intent(in) :: this
      class(SubspaceBraKet),       intent(in) :: braket
      type(AnharmonicData),        intent(in) :: anharmonic_data
      real(dp)                                :: output
    end function
  end interface
  
  interface
    recursive module function potential_energy_BasisState_PotentialBasePointer(this,bra,ket,subspace,subspace_basis,anharmonic_data) result(output) 
      class(PotentialBasePointer), intent(in)           :: this
      class(BasisState),           intent(in)           :: bra
      class(BasisState),           intent(in), optional :: ket
      type(DegenerateSubspace),    intent(in)           :: subspace
      class(SubspaceBasis),        intent(in)           :: subspace_basis
      type(AnharmonicData),        intent(in)           :: anharmonic_data
      real(dp)                                          :: output
    end function
  end interface
  
  interface
    ! I/O.
    module subroutine read_PotentialBasePointer(this,input) 
      class(PotentialBasePointer), intent(out) :: this
      type(String),            intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_PotentialBasePointer(this) result(output) 
      class(PotentialBasePointer), intent(in) :: this
      type(String), allocatable               :: output(:)
    end function
  end interface
  
  interface PotentialBasePointer
    module function new_PotentialBasePointer_Strings(input) result(this) 
      type(String), intent(in)   :: input(:)
      type(PotentialBasePointer) :: this
    end function
  
    impure elemental module function new_PotentialBasePointer_StringArray(input) result(this) 
      type(StringArray), intent(in) :: input
      type(PotentialBasePointer)        :: this
    end function
  end interface
  
  interface StressBasePointer
    ! ----------------------------------------------------------------------
    ! StressBasePointer methods.
    ! ----------------------------------------------------------------------
    ! Construct a StressBasePointer from any type which extends StressBase.
    impure elemental module function new_StressBasePointer(stress) &
       & result(this) 
      class(StressBase), intent(in) :: stress
      type(StressBasePointer)       :: this
    end function
  end interface
  
  interface
    ! Checks that the pointer has been allocated before it is used.
    impure elemental module subroutine check_StressBasePointer(this) 
      class(StressBasePointer), intent(in) :: this
    end subroutine
  end interface
  
  interface
    ! Type representation.
    impure elemental module function representation_StressBasePointer() &
       & result(output) 
      type(String) :: output
    end function
  end interface
  
  interface
    impure elemental module subroutine zero_stress_StressBasePointer(this) 
      class(StressBasePointer), intent(inout) :: this
    end subroutine
  end interface
  
  interface
    impure elemental module subroutine add_constant_StressBasePointer(this, &
       & input) 
      class(StressBasePointer), intent(inout) :: this
      type(RealMatrix),         intent(in)    :: input
    end subroutine
  end interface
  
  interface
    impure elemental module function stress_RealModeDisplacement_StressBasePointer(this,displacement) result(output) 
      class(StressBasePointer),       intent(in) :: this
      type(RealModeDisplacement), intent(in) :: displacement
      type(RealMatrix)                       :: output
    end function
  end interface
  
  interface
    impure elemental module function stress_ComplexModeDisplacement_StressBasePointer(   this,displacement) result(output) 
      class(StressBasePointer),          intent(in) :: this
      type(ComplexModeDisplacement), intent(in) :: displacement
      type(ComplexMatrix)                       :: output
    end function
  end interface
  
  interface
    impure elemental module subroutine braket_SubspaceBraKet_StressBasePointer(this,braket,whole_subspace,anharmonic_data) 
      class(StressBasePointer), intent(inout)        :: this
      class(SubspaceBraKet),    intent(in)           :: braket
      logical,                  intent(in), optional :: whole_subspace
      type(AnharmonicData),     intent(in)           :: anharmonic_data
    end subroutine
  end interface
  
  interface
    impure elemental module subroutine braket_BasisState_StressBasePointer(this,bra,ket,subspace,subspace_basis,whole_subspace,anharmonic_data) 
      class(StressBasePointer), intent(inout)        :: this
      class(BasisState),        intent(in)           :: bra
      class(BasisState),        intent(in), optional :: ket
      type(DegenerateSubspace), intent(in)           :: subspace
      class(SubspaceBasis),     intent(in)           :: subspace_basis
      logical,                  intent(in), optional :: whole_subspace
      type(AnharmonicData),     intent(in)           :: anharmonic_data
    end subroutine
  end interface
  
  interface
    impure elemental module subroutine braket_BasisStates_StressBasePointer(this,states,subspace,subspace_basis,whole_subspace,anharmonic_data) 
      class(StressBasePointer), intent(inout)        :: this
      class(BasisStates),       intent(inout)        :: states
      type(DegenerateSubspace), intent(in)           :: subspace
      class(SubspaceBasis),     intent(in)           :: subspace_basis
      logical,                  intent(in), optional :: whole_subspace
      type(AnharmonicData),     intent(in)           :: anharmonic_data
    end subroutine
  end interface
  
  interface
    impure elemental module function harmonic_expectation_StressBasePointer(this,frequency,thermal_energy,supercell_size,anharmonic_data) result(output) 
      class(StressBasePointer), intent(in) :: this
      real(dp),             intent(in) :: frequency
      real(dp),             intent(in) :: thermal_energy
      integer,              intent(in) :: supercell_size
      type(AnharmonicData), intent(in) :: anharmonic_data
      type(RealMatrix)                 :: output
    end function
  end interface
  
  interface
    recursive module function potential_stress_SubspaceBraKet_StressBasePointer(this,braket,anharmonic_data) result(output) 
      class(StressBasePointer),  intent(in) :: this
      class(SubspaceBraKet), intent(in) :: braket
      type(AnharmonicData),  intent(in) :: anharmonic_data
      type(RealMatrix)                  :: output
    end function
  end interface
  
  interface
    recursive module function potential_stress_BasisState_StressBasePointer(this,bra,ket,subspace,subspace_basis,anharmonic_data) result(output) 
      class(StressBasePointer),     intent(in)           :: this
      class(BasisState),        intent(in)           :: bra
      class(BasisState),        intent(in), optional :: ket
      type(DegenerateSubspace), intent(in)           :: subspace
      class(SubspaceBasis),     intent(in)           :: subspace_basis
      type(AnharmonicData),     intent(in)           :: anharmonic_data
      type(RealMatrix)                               :: output
    end function
  end interface
  
  interface
    ! I/O.
    module subroutine read_StressBasePointer(this,input) 
      class(StressBasePointer), intent(out) :: this
      type(String),             intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_StressBasePointer(this) result(output) 
      class(StressBasePointer), intent(in) :: this
      type(String), allocatable            :: output(:)
    end function
  end interface
  
  interface StressBasePointer
    module function new_StressBasePointer_Strings(input) result(this) 
      type(String), intent(in) :: input(:)
      type(StressBasePointer)  :: this
    end function
  
    impure elemental module function new_StressBasePointer_StringArray(input) &
       & result(this) 
      type(StringArray), intent(in) :: input
      type(StressBasePointer)       :: this
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Concrete methods of abstract classes.
    ! ----------------------------------------------------------------------
    ! SubspaceBasis methods.
    impure elemental module subroutine process_subspace_potential_SubspaceBasis(this,potential,states,subspace,anharmonic_data) 
      class(SubspaceBasis),     intent(in)            :: this
      class(PotentialBase),     intent(inout)         :: potential
      class(BasisStates),       intent(inout), target :: states
      type(DegenerateSubspace), intent(in)            :: subspace
      type(AnharmonicData),     intent(in)            :: anharmonic_data
    end subroutine
  end interface
  
  interface
    impure elemental module subroutine process_subspace_stress_SubspaceBasis(this,stress,states,subspace,anharmonic_data) 
      class(SubspaceBasis),     intent(in)            :: this
      class(stressBase),        intent(inout)         :: stress
      class(BasisStates),       intent(inout), target :: states
      type(DegenerateSubspace), intent(in)            :: subspace
      type(AnharmonicData),     intent(in)            :: anharmonic_data
    end subroutine
  end interface
  
  interface
    module function select_symmetries_SubspaceBasis(this,symmetries, &
       & anharmonic_data) result(output) 
      class(SubspaceBasis),     intent(in)  :: this
      type(DegenerateSymmetry), intent(in)  :: symmetries(:)
      type(AnharmonicData),     intent(in)  :: anharmonic_data
      type(DegenerateSymmetry), allocatable :: output(:)
    end function
  end interface
  
  interface
    ! Concrete PotentialBase methods.
    impure elemental module function undisplaced_energy(this) result(output) 
      class(PotentialBase), intent(in) :: this
      real(dp)                         :: output
    end function
  end interface
  
  interface
    recursive module function potential_energy_SubspaceBraKet_PotentialBase(this,braket,anharmonic_data) result(output) 
      class(PotentialBase),  intent(in) :: this
      class(SubspaceBraKet), intent(in) :: braket
      type(AnharmonicData),  intent(in) :: anharmonic_data
      real(dp)                          :: output
    end function
  end interface
  
  interface
    recursive module function potential_energy_BasisState_PotentialBase(this,bra,ket,subspace,subspace_basis,anharmonic_data) result(output) 
      class(PotentialBase),     intent(in)           :: this
      class(BasisState),        intent(in)           :: bra
      class(BasisState),        intent(in), optional :: ket
      type(DegenerateSubspace), intent(in)           :: subspace
      class(SubspaceBasis),     intent(in)           :: subspace_basis
      type(AnharmonicData),     intent(in)           :: anharmonic_data
      real(dp)                                       :: output
    end function
  end interface
  
  interface
    ! StressBase methods.
    impure elemental module function undisplaced_stress(this) result(output) 
      class(StressBase), intent(in) :: this
      type(RealMatrix)              :: output
    end function
  end interface
  
  interface
    recursive module function potential_stress_SubspaceBraKet_StressBase(this,braket,anharmonic_data) result(output) 
      class(StressBase),     intent(in) :: this
      class(SubspaceBraKet), intent(in) :: braket
      type(AnharmonicData),  intent(in) :: anharmonic_data
      type(RealMatrix)                  :: output
    end function
  end interface
  
  interface
    recursive module function potential_stress_BasisState_StressBase(this, &
       & bra,ket,subspace,subspace_basis,anharmonic_data) result(output) 
      class(StressBase),        intent(in)           :: this
      class(BasisState),        intent(in)           :: bra
      class(BasisState),        intent(in), optional :: ket
      type(DegenerateSubspace), intent(in)           :: subspace
      class(SubspaceBasis),     intent(in)           :: subspace_basis
      type(AnharmonicData),     intent(in)           :: anharmonic_data
      type(RealMatrix)                               :: output
    end function
  end interface
end module
