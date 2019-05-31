! ======================================================================
! Four abstract classes, and pointers to those classes.
!    - SubspaceBasis defines the basis of states spanning a subspace.
!    - BasisState defines a specific state, in terms of the basis.
!    - BasisStates defines a set of states, again in terms of the basis.
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
module abstract_classes_module
  use common_module
  
  use subspace_coupling_module
  use anharmonic_data_module
  use energy_spectrum_module
  use subspace_wavefunctions_module
  use stress_prefactors_module
  use subspace_state_module
  use basis_state_module
  use basis_states_module
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
    
    procedure(modes_SubspaceBasis), public, deferred :: modes
    
    ! If ket is not given, <this|this>, otherwise <this|ket>.
    procedure(inner_product_BasisState), public, deferred :: inner_product
    
    ! Integrals of the form <i|V|j>
    generic, public :: braket =>                 &
                     & braket_ComplexMonomial,   &
                     & braket_ComplexPolynomial
    ! Either <this|V|this> or <this|V|ket>, where V is a ComplexMonomial.
    procedure(braket_ComplexMonomial_BasisState), public, deferred :: &
       & braket_ComplexMonomial
    ! Either <this|V|this> or <this|V|ket>, where V is a ComplexPolynomial.
    procedure, public :: braket_ComplexPolynomial => &
                       & braket_ComplexPolynomial_BasisState
    
    ! Either <this|T|this> or <this|T|ket>, where T is the kinetic energy.
    procedure(kinetic_energy_BasisState), public, deferred :: &
       & kinetic_energy
    
    ! Either <this|V|this> or <this|V|ket>, where V is the harmonic potential
    !    energy.
    procedure(harmonic_potential_energy_BasisState), public, deferred :: &
       & harmonic_potential_energy
    
    ! Either <this|stress|this> or <this|stress|ket>, where stress is the
    !    kinetic stress.
    procedure(kinetic_stress_BasisState), public, deferred :: &
       & kinetic_stress
    
    ! Functionality involving SubspaceStates.
    procedure(spectra_BasisStates),       public, deferred :: spectra
    procedure(wavefunctions_BasisStates), public, deferred :: wavefunctions
    
    generic, public :: integrate =>               &
                     & integrate_ComplexMonomial, &
                     & integrate_ComplexPolynomial
    procedure(integrate_ComplexMonomial_BasisStates), public, deferred :: &
       & integrate_ComplexMonomial
    procedure, public :: integrate_ComplexPolynomial => &
                       & integrate_ComplexPolynomial_BasisStates
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
    
    procedure, public :: modes => modes_SubspaceBasisPointer
    
    procedure, public :: inner_product => &
                       & inner_product_BasisStatePointer
    procedure, public :: braket_ComplexMonomial => &
                       & braket_ComplexMonomial_BasisStatePointer
    procedure, public :: kinetic_energy => &
                       & kinetic_energy_BasisStatePointer
    procedure, public :: harmonic_potential_energy => &
                       & harmonic_potential_energy_BasisStatePointer
    procedure, public :: kinetic_stress => &
                       & kinetic_stress_BasisStatePointer
    
    procedure, public :: spectra => spectra_BasisStatesPointer
    procedure, public :: wavefunctions => wavefunctions_BasisStatesPointer
    procedure, public :: integrate_ComplexMonomial => &
                       & integrate_ComplexMonomial_BasisStatesPointer
    
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
    
    ! Return the energy at zero displacement, or set this energy to zero.
    procedure, public :: undisplaced_energy
    procedure(zero_energy_PotentialData), public, deferred :: zero_energy
    
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
    
    ! Update the potential from previous iterations, either using damped
    !    iteration or a Pulay scheme.
    procedure(iterate_damped_PotentialData), public, deferred :: &
       & iterate_damped
    procedure(iterate_pulay_PotentialData), public, deferred :: &
       & iterate_pulay
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
    
    procedure, public :: iterate_damped => &
                       & iterate_damped_PotentialPointer
    procedure, public :: iterate_pulay => &
                       & iterate_pulay_PotentialPointer
    
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
    
    impure elemental function calculate_states_SubspaceBasis(this,subspace, &
       & subspace_potential,energy_convergence,no_converged_calculations,   &
       & max_pulay_iterations,pre_pulay_iterations,pre_pulay_damping,       &
       & anharmonic_data) result(output)
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
      real(dp),                 intent(in) :: energy_convergence
      integer,                  intent(in) :: no_converged_calculations
      integer,                  intent(in) :: max_pulay_iterations
      integer,                  intent(in) :: pre_pulay_iterations
      real(dp),                 intent(in) :: pre_pulay_damping
      type(AnharmonicData),     intent(in) :: anharmonic_data
      type(BasisStatesPointer)             :: output
    end function
    
    function modes_SubspaceBasis(this,subspace,anharmonic_data) result(output)
      import SubspaceBasis
      import DegenerateSubspace
      import AnharmonicData
      implicit none
      
      class(SubspaceBasis),     intent(in) :: this
      type(DegenerateSubspace), intent(in) :: subspace
      type(AnharmonicData),     intent(in) :: anharmonic_data
      integer, allocatable                 :: output(:)
    end function
    
    impure elemental function inner_product_BasisState(this,bra,ket, &
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
    
    impure elemental function braket_ComplexMonomial_BasisState(this, &
       & bra,monomial,ket,subspace,anharmonic_data) result(output)
      import SubspaceBasis
      import BasisState
      import ComplexMonomial
      import DegenerateSubspace
      import AnharmonicData
      implicit none
      
      class(SubspaceBasis),     intent(in)           :: this
      class(BasisState),        intent(in)           :: bra
      type(ComplexMonomial),    intent(in)           :: monomial
      class(BasisState),        intent(in), optional :: ket
      type(DegenerateSubspace), intent(in)           :: subspace
      type(AnharmonicData),     intent(in)           :: anharmonic_data
      type(ComplexMonomial)                          :: output
    end function
    
    impure elemental function kinetic_energy_BasisState(this,bra,ket, &
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
    
    impure elemental function harmonic_potential_energy_BasisState(this, &
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
    
    impure elemental function kinetic_stress_BasisState(this,bra,ket, &
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
    
    impure elemental function spectra_BasisStates(this,states,subspace, &
       & subspace_potential,subspace_stress,stress_prefactors,          &
       & anharmonic_data) result(output)
      import SubspaceBasis
      import BasisStates
      import DegenerateSubspace
      import PotentialData
      import StressData
      import StressPrefactors
      import AnharmonicData
      import EnergySpectra
      implicit none
      
      class(SubspaceBasis),     intent(in)           :: this
      class(BasisStates),       intent(in)           :: states
      type(DegenerateSubspace), intent(in)           :: subspace
      class(PotentialData),     intent(in)           :: subspace_potential
      class(StressData),        intent(in), optional :: subspace_stress
      type(StressPrefactors),   intent(in), optional :: stress_prefactors
      type(AnharmonicData),     intent(in)           :: anharmonic_data
      type(EnergySpectra)                            :: output
    end function
    
    impure elemental function wavefunctions_BasisStates(this,states, &
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
    
    impure elemental function integrate_ComplexMonomial_BasisStates(this, &
       & states,monomial,subspace,anharmonic_data) result(output)
      import SubspaceBasis
      import BasisStates
      import ComplexMonomial
      import DegenerateSubspace
      import QpointData
      import AnharmonicData
      implicit none
      
      class(SubspaceBasis),     intent(in) :: this
      class(BasisStates),       intent(in) :: states
      type(ComplexMonomial),    intent(in) :: monomial
      type(DegenerateSubspace), intent(in) :: subspace
      type(AnharmonicData),     intent(in) :: anharmonic_data
      type(ComplexMonomial)                :: output
    end function
    
    ! PotentialData procedures.
    impure elemental function representation_PotentialData() result(output)
      import String
      implicit none
      
      type(String) :: output
    end function
    
    subroutine generate_sampling_points_PotentialData(this,anharmonic_data, &
       & sampling_points_dir,calculation_writer,logfile)
      import PotentialData
      import AnharmonicData
      import String
      import CalculationWriter
      import OFile
      implicit none
      
      class(PotentialData),    intent(inout) :: this
      type(AnharmonicData),    intent(in)    :: anharmonic_data
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
       & anharmonic_data)
      import PotentialData
      import SubspaceState
      import AnharmonicData
      implicit none
      
      class(PotentialData),     intent(inout)        :: this
      class(SubspaceState),     intent(in)           :: bra
      class(SubspaceState),     intent(in), optional :: ket
      type(AnharmonicData),     intent(in)           :: anharmonic_data
    end subroutine
    
    subroutine braket_BasisState_PotentialData(this,bra,ket,subspace, &
       & subspace_basis,anharmonic_data)
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
      type(AnharmonicData),     intent(in)           :: anharmonic_data
    end subroutine
    
    subroutine braket_BasisStates_PotentialData(this,states,subspace, &
       & subspace_basis,anharmonic_data)
      import PotentialData
      import BasisStates
      import DegenerateSubspace
      import SubspaceBasis
      import AnharmonicData
      implicit none
      
      class(PotentialData),     intent(inout) :: this
      class(BasisStates),       intent(in)    :: states
      type(DegenerateSubspace), intent(in)    :: subspace
      class(SubspaceBasis),     intent(in)    :: subspace_basis
      type(AnharmonicData),     intent(in)    :: anharmonic_data
    end subroutine
    
    function harmonic_expectation_PotentialData(this,frequency, &
       & thermal_energy,subspace,anharmonic_data) result(output)
      import PotentialData
      import dp
      import DegenerateSubspace
      import AnharmonicData
      implicit none
      
      class(PotentialData),     intent(in) :: this
      real(dp),                 intent(in) :: frequency
      real(dp),                 intent(in) :: thermal_energy
      type(DegenerateSubspace), intent(in) :: subspace
      type(AnharmonicData),     intent(in) :: anharmonic_data
      real(dp)                             :: output
    end function
    
    impure elemental function iterate_damped_PotentialData(this, &
       & new_potential,damping,anharmonic_data) result(output)
      import PotentialData
      import dp
      import AnharmonicData
      import PotentialPointer
      implicit none
      
      class(PotentialData), intent(in)  :: this
      class(PotentialData), intent(in)  :: new_potential
      real(dp),             intent(in)  :: damping
      type(AnharmonicData), intent(in)  :: anharmonic_data
      type(PotentialPointer)            :: output
    end function
    
    function iterate_pulay_PotentialData(this,input_potentials, &
       & output_potentials,anharmonic_data) result(output)
      import PotentialData
      import PotentialPointer
      import AnharmonicData
      implicit none
      
      class(PotentialData),   intent(in)  :: this
      type(PotentialPointer), intent(in)  :: input_potentials(:)
      type(PotentialPointer), intent(in)  :: output_potentials(:)
      type(AnharmonicData),   intent(in)  :: anharmonic_data
      type(PotentialPointer)              :: output
    end function
    
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
    
    subroutine braket_SubspaceState_StressData(this,bra,ket,anharmonic_data)
      import StressData
      import SubspaceState
      import AnharmonicData
      implicit none
      
      class(StressData),        intent(inout)        :: this
      class(SubspaceState),     intent(in)           :: bra
      class(SubspaceState),     intent(in), optional :: ket
      type(AnharmonicData),     intent(in)           :: anharmonic_data
    end subroutine
    
    subroutine braket_BasisState_StressData(this,bra,ket,subspace, &
       & subspace_basis,anharmonic_data)
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
      type(AnharmonicData),     intent(in)           :: anharmonic_data
    end subroutine
    
    subroutine braket_BasisStates_StressData(this,states,subspace, &
       & subspace_basis,anharmonic_data)
      import StressData
      import BasisStates
      import DegenerateSubspace
      import SubspaceBasis
      import AnharmonicData
      implicit none
      
      class(StressData),        intent(inout) :: this
      class(BasisStates),       intent(in)    :: states
      type(DegenerateSubspace), intent(in)    :: subspace
      class(SubspaceBasis),     intent(in)    :: subspace_basis
      type(AnharmonicData),     intent(in)    :: anharmonic_data
    end subroutine
    
    function harmonic_expectation_StressData(this,frequency, &
       & thermal_energy,subspace,anharmonic_data) result(output)
      import StressData
      import dp
      import DegenerateSubspace
      import AnharmonicData
      import RealMatrix
      implicit none
      
      class(StressData),        intent(in) :: this
      real(dp),                 intent(in) :: frequency
      real(dp),                 intent(in) :: thermal_energy
      type(DegenerateSubspace), intent(in) :: subspace
      type(AnharmonicData),     intent(in) :: anharmonic_data
      type(RealMatrix)                     :: output
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
   & subspace,subspace_potential,energy_convergence,                      &
   & no_converged_calculations,max_pulay_iterations,pre_pulay_iterations, &
   & pre_pulay_damping,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceBasisPointer), intent(in) :: this
  type(DegenerateSubspace),    intent(in) :: subspace
  class(PotentialData),        intent(in) :: subspace_potential
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
                                       & energy_convergence,        &
                                       & no_converged_calculations, &
                                       & max_pulay_iterations,      &
                                       & pre_pulay_iterations,      &
                                       & pre_pulay_damping,         &
                                       & anharmonic_data            )
end function

function modes_SubspaceBasisPointer(this,subspace,anharmonic_data) &
   & result(output)
  implicit none
  
  class(SubspaceBasisPointer), intent(in) :: this
  type(DegenerateSubspace),    intent(in) :: subspace
  type(AnharmonicData),        intent(in) :: anharmonic_data
  integer, allocatable                    :: output(:)
  
  call this%check()
  
  output = this%basis_%modes(subspace,anharmonic_data)
end function

impure elemental function inner_product_BasisStatePointer(this,bra,ket, &
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

impure elemental function braket_ComplexMonomial_BasisStatePointer(this,bra, &
   & monomial,ket,subspace,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceBasisPointer), intent(in)           :: this
  class(BasisState),           intent(in)           :: bra
  type(ComplexMonomial),       intent(in)           :: monomial
  class(BasisState),           intent(in), optional :: ket
  type(DegenerateSubspace),    intent(in)           :: subspace
  type(AnharmonicData),        intent(in)           :: anharmonic_data
  type(ComplexMonomial)                             :: output
  
  call this%check()
  
  output = this%basis_%braket( bra,            &
                             & monomial,       &
                             & ket,            &
                             & subspace,       &
                             & anharmonic_data )
end function

impure elemental function kinetic_energy_BasisStatePointer(this,bra,ket, &
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

impure elemental function harmonic_potential_energy_BasisStatePointer( &
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

impure elemental function kinetic_stress_BasisStatePointer(this,bra,ket, &
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

impure elemental function spectra_BasisStatesPointer(this,states,subspace, &
   & subspace_potential,subspace_stress,stress_prefactors,anharmonic_data) &
   & result(output)
  implicit none
  
  class(SubspaceBasisPointer), intent(in)           :: this
  class(BasisStates),          intent(in)           :: states
  type(DegenerateSubspace),    intent(in)           :: subspace
  class(PotentialData),        intent(in)           :: subspace_potential
  class(StressData),           intent(in), optional :: subspace_stress
  type(StressPrefactors),      intent(in), optional :: stress_prefactors
  type(AnharmonicData),        intent(in)           :: anharmonic_data
  type(EnergySpectra)                               :: output
  
  call this%check()
  
  output = this%basis_%spectra( states,             &
                              & subspace,           &
                              & subspace_potential, &
                              & subspace_stress,    &
                              & stress_prefactors,  &
                              & anharmonic_data     )
end function

impure elemental function wavefunctions_BasisStatesPointer(this,states, &
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

impure elemental function integrate_ComplexMonomial_BasisStatesPointer(this, &
   & states,monomial,subspace,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceBasisPointer), intent(in) :: this
  class(BasisStates),          intent(in) :: states
  type(ComplexMonomial),       intent(in) :: monomial
  type(DegenerateSubspace),    intent(in) :: subspace
  type(AnharmonicData),        intent(in) :: anharmonic_data
  type(ComplexMonomial)                   :: output
  
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
   & sampling_points_dir,calculation_writer,logfile)
  implicit none
  
  class(PotentialPointer), intent(inout) :: this
  type(AnharmonicData),    intent(in)    :: anharmonic_data
  type(String),            intent(in)    :: sampling_points_dir
  type(CalculationWriter), intent(inout) :: calculation_writer
  type(OFile),             intent(inout) :: logfile
  
  call this%check()
  
  call this%potential_%generate_sampling_points( anharmonic_data,     &
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
   & anharmonic_data)
  implicit none
  
  class(PotentialPointer), intent(inout)        :: this
  class(SubspaceState),    intent(in)           :: bra
  class(SubspaceState),    intent(in), optional :: ket
  type(AnharmonicData),    intent(in)           :: anharmonic_data
  
  call this%check()
  
  call this%potential_%braket(bra,ket,anharmonic_data)
end subroutine

subroutine braket_BasisState_PotentialPointer(this,bra,ket,subspace, &
   & subspace_basis,anharmonic_data)
  implicit none
  
  class(PotentialPointer),  intent(inout)        :: this
  class(BasisState),        intent(in)           :: bra
  class(BasisState),        intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  
  call this%check()
  
  call this%potential_%braket( bra,            &
                             & ket,            &
                             & subspace,       &
                             & subspace_basis, &
                             & anharmonic_data )
end subroutine

subroutine braket_BasisStates_PotentialPointer(this,states,subspace, &
   & subspace_basis,anharmonic_data)
  implicit none
  
  class(PotentialPointer),  intent(inout) :: this
  class(BasisStates),       intent(in)    :: states
  type(DegenerateSubspace), intent(in)    :: subspace
  class(SubspaceBasis),     intent(in)    :: subspace_basis
  type(AnharmonicData),     intent(in)    :: anharmonic_data
  
  call this%check()
  
  call this%potential_%braket( states,         &
                             & subspace,       &
                             & subspace_basis, &
                             & anharmonic_data )
end subroutine

function harmonic_expectation_PotentialPointer(this,frequency, &
   & thermal_energy,subspace,anharmonic_data) result(output)
  implicit none
  
  class(PotentialPointer),  intent(in) :: this
  real(dp),                 intent(in) :: frequency
  real(dp),                 intent(in) :: thermal_energy
  type(DegenerateSubspace), intent(in) :: subspace
  type(AnharmonicData),     intent(in) :: anharmonic_data
  real(dp)                             :: output
  
  call this%check()
  
  output = this%potential_%harmonic_expectation( frequency,      &
                                               & thermal_energy, &
                                               & subspace,       &
                                               & anharmonic_data )
end function

! Generates the next iteration of the potential, either following a damped
!    iterative scheme or a pulay scheme.
impure elemental function iterate_damped_PotentialPointer(this,new_potential, &
   & damping,anharmonic_data) result(output)
  implicit none
  
  class(PotentialPointer), intent(in) :: this
  class(PotentialData),    intent(in) :: new_potential
  real(dp),                intent(in) :: damping
  type(AnharmonicData),    intent(in) :: anharmonic_data
  type(PotentialPointer)              :: output
  
  class(PotentialData), allocatable :: potential
  
  call this%check()
  
  select type(new_potential); type is(PotentialPointer)
    call new_potential%check()
    potential = new_potential%potential_
  class default
    potential = new_potential
  end select
  
  output = this%potential_%iterate_damped( potential,      &
                                         & damping,        &
                                         & anharmonic_data )
end function

function iterate_pulay_PotentialPointer(this,input_potentials, &
   & output_potentials,anharmonic_data) result(output)
  implicit none
  
  class(PotentialPointer), intent(in) :: this
  type(PotentialPointer),  intent(in) :: input_potentials(:)
  type(PotentialPointer),  intent(in) :: output_potentials(:)
  type(AnharmonicData),    intent(in) :: anharmonic_data
  type(PotentialPointer)              :: output
  
  call this%check()
  
  output = this%potential_%iterate_pulay( input_potentials,  &
                                        & output_potentials, &
                                        & anharmonic_data    )
end function

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

subroutine braket_SubspaceState_StressPointer(this,bra,ket,anharmonic_data)
  implicit none
  
  class(StressPointer),     intent(inout)        :: this
  class(SubspaceState),     intent(in)           :: bra
  class(SubspaceState),     intent(in), optional :: ket
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  
  call this%check()
  
  call this%stress_%braket(bra,ket,anharmonic_data)
end subroutine

subroutine braket_BasisState_StressPointer(this,bra,ket,subspace, &
   & subspace_basis,anharmonic_data)
  implicit none
  
  class(StressPointer),     intent(inout)        :: this
  class(BasisState),        intent(in)           :: bra
  class(BasisState),        intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  
  call this%check()
  
  call this%stress_%braket( bra,            &
                          & ket,            &
                          & subspace,       &
                          & subspace_basis, &
                          & anharmonic_data )
end subroutine

subroutine braket_BasisStates_StressPointer(this,states,subspace, &
   & subspace_basis,anharmonic_data)
  implicit none
  
  class(StressPointer),     intent(inout) :: this
  class(BasisStates),       intent(in)    :: states
  type(DegenerateSubspace), intent(in)    :: subspace
  class(SubspaceBasis),     intent(in)    :: subspace_basis
  type(AnharmonicData),     intent(in)    :: anharmonic_data
  
  call this%check()
  
  call this%stress_%braket( states,         &
                          & subspace,       &
                          & subspace_basis, &
                          & anharmonic_data )
end subroutine

function harmonic_expectation_StressPointer(this,frequency, &
   & thermal_energy,subspace,anharmonic_data) result(output)
  implicit none
  
  class(StressPointer),     intent(in) :: this
  real(dp),                 intent(in) :: frequency
  real(dp),                 intent(in) :: thermal_energy
  type(DegenerateSubspace), intent(in) :: subspace
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(RealMatrix)                     :: output
  
  call this%check()
  
  output = this%stress_%harmonic_expectation( frequency,      &
                                            & thermal_energy, &
                                            & subspace,       &
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
! BasisState methods.
impure elemental function braket_ComplexPolynomial_BasisState(this, &
   & bra,polynomial,ket,subspace,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceBasis),     intent(in)           :: this
  class(BasisState),        intent(in)           :: bra
  type(ComplexPolynomial),  intent(in)           :: polynomial
  class(BasisState),        intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(ComplexPolynomial)                        :: output
  
  type(ComplexMonomial), allocatable :: monomials(:)
  
  monomials = this%braket( bra,              &
                         & polynomial%terms, &
                         & ket,              &
                         & subspace,         &
                         & anharmonic_data   )
  output = ComplexPolynomial(monomials)
end function

! BasisStates methods.
impure elemental function integrate_ComplexPolynomial_BasisStates(this, &
   & states,polynomial,subspace,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceBasis),     intent(in) :: this
  class(BasisStates),       intent(in) :: states
  type(ComplexPolynomial),  intent(in) :: polynomial
  type(DegenerateSubspace), intent(in) :: subspace
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(ComplexPolynomial)              :: output
  
  type(ComplexMonomial), allocatable :: monomials(:)
  
  monomials = this%integrate( states,           &
                            & polynomial%terms, &
                            & subspace,         &
                            & anharmonic_data   )
  output = ComplexPolynomial(monomials)
end function

! PotentialData methods.
function undisplaced_energy(this) result(output)
  implicit none
  
  class(PotentialData), intent(in) :: this
  real(dp)                         :: output
  
  type(RealModeDisplacement) :: zero_displacement
  
  zero_displacement = RealModeDisplacement([RealSingleDisplacement::])
  
  output = this%energy(zero_displacement)
end function

! ----------------------------------------------------------------------
! StressData methods.
! ----------------------------------------------------------------------
function undisplaced_stress(this) result(output)
  implicit none
  
  class(StressData), intent(in) :: this
  type(RealMatrix)              :: output
  
  type(RealModeDisplacement) :: zero_displacement
  
  zero_displacement = RealModeDisplacement([RealSingleDisplacement::])
  
  output = this%stress(zero_displacement)
end function
end module
