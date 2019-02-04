! ======================================================================
! Four abstract classes, and pointers to those classes.
!    - SubspaceBasis defines the basis of states spanning a subspace.
!    - SubspaceState defines a specific state, in terms of the basis.
!    - SubspaceStates defines a set of states, again in terms of the basis.
!    - PotentialData defines a potential.
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
  
  use anharmonic_data_module
  use energy_spectrum_module
  use subspace_wavefunctions_module
  use subspace_wavefunctions_pointer_module
  implicit none
  
  private
  
  public :: SubspaceBasis
  public :: SubspaceBasisPointer
  
  public :: SubspaceState
  public :: SubspaceStatePointer
  
  public :: SubspaceStates
  public :: SubspaceStatesPointer
  
  public :: PotentialData
  public :: PotentialPointer
  
  public :: generate_subspace_potentials
  public :: braket
  public :: kinetic_energy
  public :: harmonic_potential_energy
  public :: potential_energy
  
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
    
    procedure(initial_states_SubspaceBasis), public, deferred :: initial_states
    procedure(calculate_states_SubspaceBasis), public, deferred :: &
       & calculate_states
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
    
    ! I/O.
    procedure, public :: read  => read_SubspaceBasisPointer
    procedure, public :: write => write_SubspaceBasisPointer
  end type
  
  ! An array of all types which extend SubspaceBasisPointer.
  ! This array will be filled in by startup routines.
  type(SubspaceBasisPointer), allocatable :: TYPES_SubspaceBasis(:)
  
  type, abstract, extends(Stringsable) :: SubspaceState
    integer :: subspace_id
  contains
    procedure(representation_SubspaceState), public, deferred, nopass :: &
       & representation
    procedure, public :: startup => startup_SubspaceState
    
    ! Integrals of the form <i|V|j>
    generic, public :: braket =>                 &
                     & braket_SubspaceState,     &
                     & braket_ComplexUnivariate, &
                     & braket_ComplexMonomial,   &
                     & braket_ComplexPolynomial
    ! If ket is not given, <this|this>, otherwise <this|ket>.
    procedure(braket_SubspaceState_SubspaceState), public, deferred :: &
       & braket_SubspaceState
    ! Either <this|V|this> or <this|V|ket>, where V is a ComplexUnivariate.
    procedure(braket_ComplexUnivariate_SubspaceState), public, deferred :: &
       & braket_ComplexUnivariate
    ! Either <this|V|this> or <this|V|ket>, where V is a ComplexMonomial.
    procedure(braket_ComplexMonomial_SubspaceState), public, deferred :: &
       & braket_ComplexMonomial
    ! Either <this|V|this> or <this|V|ket>, where V is a ComplexPolynomial.
    procedure, public :: braket_ComplexPolynomial => &
                       & braket_ComplexPolynomial_SubspaceState
    
    ! Either <this|T|this> or <this|T|ket>, where T is the kinetic energy.
    procedure(kinetic_energy_SubspaceState), public, deferred :: &
       & kinetic_energy
    
    ! Either <this|V|this> or <this|V|ket>, where V is the harmonic potential
    !    energy.
    procedure(harmonic_potential_energy_SubspaceState), public, deferred :: &
       & harmonic_potential_energy
  end type
  
  type, extends(SubspaceState) :: SubspaceStatePointer
    type(String),                      private :: representation_
    class(SubspaceState), allocatable, private :: state_
  contains
    procedure, private :: check => check_SubspaceStatePointer
    
    procedure, public, nopass :: representation => &
                               & representation_SubspaceStatePointer
    
    procedure, public :: braket_SubspaceState => &
                       & braket_SubspaceState_SubspaceStatePointer
    procedure, public :: braket_ComplexUnivariate => &
                       & braket_ComplexUnivariate_SubspaceStatePointer
    procedure, public :: braket_ComplexMonomial => &
                       & braket_ComplexMonomial_SubspaceStatePointer
    procedure, public :: kinetic_energy => &
                       & kinetic_energy_SubspaceStatePointer
    procedure, public :: harmonic_potential_energy => &
                       & harmonic_potential_energy_SubspaceStatePointer
    
    ! I/O.
    procedure, public :: read  => read_SubspaceStatePointer
    procedure, public :: write => write_SubspaceStatePointer
  end type
  
  ! An array of all types which extend SubspaceStatePointer.
  ! This array will be filled in by startup routines.
  type(SubspaceStatePointer), allocatable :: TYPES_SubspaceState(:)
  
  type, abstract, extends(Stringsable) :: SubspaceStates
  contains
    procedure(representation_SubspaceStates), public, deferred, nopass :: &
       & representation
    procedure, public :: startup => startup_SubspaceStates
    
    procedure(states_SubspaceStates),        public, deferred :: states
    procedure(spectra_SubspaceStates),       public, deferred :: spectra
    procedure(wavefunctions_SubspaceStates), public, deferred :: wavefunctions
    procedure(integrate_SubspaceStates),     public, deferred :: integrate
  end type
  
  type, extends(SubspaceStates) :: SubspaceStatesPointer
    type(String),                       private :: representation_
    class(SubspaceStates), allocatable, private :: states_
  contains
    procedure, private :: check => check_SubspaceStatesPointer
    
    procedure, public, nopass :: representation => &
                               & representation_SubspaceStatesPointer
    
    procedure, public :: states => states_SubspaceStatesPointer
    procedure, public :: spectra => spectra_SubspaceStatesPointer
    procedure, public :: wavefunctions => wavefunctions_SubspaceStatesPointer
    procedure, public :: integrate => integrate_SubspaceStatesPointer
    
    ! I/O.
    procedure, public :: read  => read_SubspaceStatesPointer
    procedure, public :: write => write_SubspaceStatesPointer
  end type
  
  ! An array of all types which extend SubspaceStatesPointer.
  ! This array will be filled in by startup routines.
  type(SubspaceStatesPointer), allocatable :: TYPES_SubspaceStates(:)
  
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
    procedure(braket_PotentialData), public, deferred :: braket
    
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
    
    procedure, public :: zero_energy => zero_energy_PotentialPointer
    
    procedure, public :: energy_RealModeDisplacement => &
                       & energy_RealModeDisplacement_PotentialPointer
    procedure, public :: energy_ComplexModeDisplacement => &
                       & energy_ComplexModeDisplacement_PotentialPointer
    procedure, public :: force_RealModeDisplacement => &
                       & force_RealModeDisplacement_PotentialPointer
    procedure, public :: force_ComplexModeDisplacement => &
                       & force_ComplexModeDisplacement_PotentialPointer
    
    procedure, public :: braket => braket_PotentialPointer
    
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
      import SubspaceStatesPointer
      implicit none
      
      class(SubspaceBasis),     intent(in) :: this
      type(DegenerateSubspace), intent(in) :: subspace
      type(AnharmonicData),     intent(in) :: anharmonic_data
      type(SubspaceStatesPointer)          :: output
    end function
    
    impure elemental function calculate_states_SubspaceBasis(this,subspace, &
       & subspace_potential,anharmonic_data) result(output)
      import SubspaceBasis
      import DegenerateSubspace
      import PotentialData
      import AnharmonicData
      import SubspaceStatesPointer
      implicit none
      
      class(SubspaceBasis),     intent(in) :: this
      type(DegenerateSubspace), intent(in) :: subspace
      class(PotentialData),     intent(in) :: subspace_potential
      type(AnharmonicData),     intent(in) :: anharmonic_data
      type(SubspaceStatesPointer)          :: output
    end function
    
    ! SubspaceState procedures.
    impure elemental function representation_SubspaceState() result(output)
      import String
      implicit none
      
      type(String) :: output
    end function
    
    impure elemental function braket_SubspaceState_SubspaceState(this,ket, &
       & subspace,subspace_basis,anharmonic_data) result(output)
      import SubspaceState
      import DegenerateSubspace
      import SubspaceBasis
      import AnharmonicData
      import dp
      implicit none
      
      class(SubspaceState),     intent(in)           :: this
      class(SubspaceState),     intent(in), optional :: ket
      type(DegenerateSubspace), intent(in)           :: subspace
      class(SubspaceBasis),     intent(in)           :: subspace_basis
      type(AnharmonicData),     intent(in)           :: anharmonic_data
      real(dp)                                       :: output
    end function
    
    impure elemental function braket_ComplexUnivariate_SubspaceState(this, &
       & univariate,ket,subspace,subspace_basis,anharmonic_data)           &
       & result(output)
      import SubspaceState
      import ComplexUnivariate
      import DegenerateSubspace
      import SubspaceBasis
      import AnharmonicData
      import ComplexMonomial
      implicit none
      
      class(SubspaceState),     intent(in)           :: this
      type(ComplexUnivariate),  intent(in)           :: univariate
      class(SubspaceState),     intent(in), optional :: ket
      type(DegenerateSubspace), intent(in)           :: subspace
      class(SubspaceBasis),     intent(in)           :: subspace_basis
      type(AnharmonicData),     intent(in)           :: anharmonic_data
      type(ComplexMonomial)                          :: output
    end function
    
    impure elemental function braket_ComplexMonomial_SubspaceState(this, &
       & monomial,ket,subspace,subspace_basis,anharmonic_data) result(output)
      import SubspaceState
      import ComplexMonomial
      import DegenerateSubspace
      import SubspaceBasis
      import AnharmonicData
      implicit none
      
      class(SubspaceState),     intent(in)           :: this
      type(ComplexMonomial),    intent(in)           :: monomial
      class(SubspaceState),     intent(in), optional :: ket
      type(DegenerateSubspace), intent(in)           :: subspace
      class(SubspaceBasis),     intent(in)           :: subspace_basis
      type(AnharmonicData),     intent(in)           :: anharmonic_data
      type(ComplexMonomial)                          :: output
    end function
    
    impure elemental function kinetic_energy_SubspaceState(this,ket, &
       & subspace,subspace_basis,anharmonic_data) result(output)
      import SubspaceState
      import DegenerateSubspace
      import SubspaceBasis
      import AnharmonicData
      import dp
      implicit none
      
      class(SubspaceState),     intent(in)           :: this
      class(SubspaceState),     intent(in), optional :: ket
      type(DegenerateSubspace), intent(in)           :: subspace
      class(SubspaceBasis),     intent(in)           :: subspace_basis
      type(AnharmonicData),     intent(in)           :: anharmonic_data
      real(dp)                                       :: output
    end function
    
    impure elemental function harmonic_potential_energy_SubspaceState(this, &
       & ket,subspace,subspace_basis,anharmonic_data) result(output)
      import SubspaceState
      import DegenerateSubspace
      import SubspaceBasis
      import AnharmonicData
      import dp
      implicit none
      
      class(SubspaceState),     intent(in)           :: this
      class(SubspaceState),     intent(in), optional :: ket
      type(DegenerateSubspace), intent(in)           :: subspace
      class(SubspaceBasis),     intent(in)           :: subspace_basis
      type(AnharmonicData),     intent(in)           :: anharmonic_data
      real(dp)                                       :: output
    end function
    
    ! SubspaceStates procedures.
    impure elemental function representation_SubspaceStates() result(output)
      import String
      implicit none
      
      type(String) :: output
    end function
    
    function states_SubspaceStates(this,subspace,subspace_basis, &
       & anharmonic_data) result(output)
      import SubspaceStates
      import DegenerateSubspace
      import SubspaceBasis
      import AnharmonicData
      import SubspaceStatePointer
      implicit none
      
      class(SubspaceStates),    intent(in)    :: this
      type(DegenerateSubspace), intent(in)    :: subspace
      class(SubspaceBasis),     intent(in)    :: subspace_basis
      type(AnharmonicData),     intent(in)    :: anharmonic_data
      type(SubspaceStatePointer), allocatable :: output(:)
    end function
    
    impure elemental function spectra_SubspaceStates(this,subspace, &
       & subspace_basis,anharmonic_data) result(output)
      import SubspaceStates
      import DegenerateSubspace
      import SubspaceBasis
      import AnharmonicData
      import EnergySpectra
      implicit none
      
      class(SubspaceStates),    intent(in) :: this
      type(DegenerateSubspace), intent(in) :: subspace
      class(SubspaceBasis),     intent(in) :: subspace_basis
      type(AnharmonicData),     intent(in) :: anharmonic_data
      type(EnergySpectra)                  :: output
    end function
    
    impure elemental function wavefunctions_SubspaceStates(this,subspace, &
       & subspace_basis,anharmonic_data) result(output)
      import SubspaceStates
      import DegenerateSubspace
      import SubspaceBasis
      import AnharmonicData
      import SubspaceWavefunctionsPointer
      implicit none
      
      class(SubspaceStates),    intent(in) :: this
      type(DegenerateSubspace), intent(in) :: subspace
      class(SubspaceBasis),     intent(in) :: subspace_basis
      type(AnharmonicData),     intent(in) :: anharmonic_data
      type(SubspaceWavefunctionsPointer)   :: output
    end function
    
    impure elemental function integrate_SubspaceStates(this,potential, &
       & subspace,subspace_basis,anharmonic_data) result(output)
      import SubspaceStates
      import PotentialData
      import DegenerateSubspace
      import SubspaceBasis
      import AnharmonicData
      import PotentialPointer
      implicit none
      
      class(SubspaceStates),    intent(in) :: this
      class(PotentialData),     intent(in) :: potential
      type(DegenerateSubspace), intent(in) :: subspace
      class(SubspaceBasis),     intent(in) :: subspace_basis
      type(AnharmonicData),     intent(in) :: anharmonic_data
      type(PotentialPointer)               :: output
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
    
    subroutine generate_potential_PotentialData(this,anharmonic_data,      &
       & weighted_energy_force_ratio,calculate_stress,sampling_points_dir, &
       & calculation_reader,logfile)
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
      logical,                 intent(in)    :: calculate_stress
      type(String),            intent(in)    :: sampling_points_dir
      type(CalculationReader), intent(inout) :: calculation_reader
      type(OFile),             intent(inout) :: logfile
    end subroutine
    
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
    
    function braket_PotentialData(this,bra,ket,subspace,subspace_basis, &
       & anharmonic_data) result(output)
      import PotentialData
      import SubspaceState
      import DegenerateSubspace
      import SubspaceBasis
      import AnharmonicData
      import PotentialPointer
      implicit none
      
      class(PotentialData),     intent(in)           :: this
      class(SubspaceState),     intent(in)           :: bra
      class(SubspaceState),     intent(in), optional :: ket
      type(DegenerateSubspace), intent(in)           :: subspace
      class(SubspaceBasis),     intent(in)           :: subspace_basis
      type(AnharmonicData),     intent(in)           :: anharmonic_data
      type(PotentialPointer)                         :: output
    end function
    
    function harmonic_expectation_PotentialData(this,frequency, &
       & thermal_energy,no_states,subspace,anharmonic_data) result(output)
      import PotentialData
      import dp
      import DegenerateSubspace
      import AnharmonicData
      implicit none
      
      class(PotentialData),     intent(in) :: this
      real(dp),                 intent(in) :: frequency
      real(dp),                 intent(in) :: thermal_energy
      integer,                  intent(in) :: no_states
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
  end interface
  
  ! --------------------------------------------------
  ! Pointer constructor interfaces.
  ! --------------------------------------------------
  interface SubspaceBasisPointer
    module procedure new_SubspaceBasisPointer
    module procedure new_SubspaceBasisPointer_Strings
    module procedure new_SubspaceBasisPointer_StringArray
  end interface
  
  interface SubspaceStatePointer
    module procedure new_SubspaceStatePointer
    module procedure new_SubspaceStatePointer_Strings
    module procedure new_SubspaceStatePointer_StringArray
  end interface
  
  interface SubspaceStatesPointer
    module procedure new_SubspaceStatesPointer
    module procedure new_SubspaceStatesPointer_Strings
    module procedure new_SubspaceStatesPointer_StringArray
  end interface
  
  interface PotentialPointer
    module procedure new_PotentialPointer
    module procedure new_PotentialPointer_Strings
    module procedure new_PotentialPointer_StringArray
  end interface
  
  ! --------------------------------------------------
  ! The braket method, which calculates <i|j> or <i|V|j>,
  !    the kinetic energy method, which calculates <i|T|j>,
  !    and the harmonic potential energy method which calculates <i|Vh|j>.
  ! --------------------------------------------------
  interface braket
    module procedure braket_state
    module procedure braket_state_state
    module procedure braket_state_ComplexUnivariate
    module procedure braket_state_ComplexUnivariate_state
    module procedure braket_state_ComplexMonomial
    module procedure braket_state_ComplexMonomial_state
    module procedure braket_state_ComplexPolynomial
    module procedure braket_state_ComplexPolynomial_state
    module procedure braket_state_potential
    module procedure braket_state_potential_state
  end interface
  
  interface kinetic_energy
    module procedure kinetic_energy_state
    module procedure kinetic_energy_state_state
  end interface
  
  interface harmonic_potential_energy
    module procedure harmonic_potential_energy_state
    module procedure harmonic_potential_energy_state_state
  end interface
  
  interface potential_energy
    module procedure potential_energy_state
    module procedure potential_energy_state_state
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

subroutine startup_SubspaceState(this)
  implicit none
  
  class(SubspaceState), intent(in) :: this
  
  integer :: i
  
  if (.not.allocated(TYPES_SubspaceState)) then
    TYPES_SubspaceState = [SubspaceStatePointer(this)]
  elseif (.not.any([(                                        &
     & this%representation()                                 &
     &    == TYPES_SubspaceState(i)%state_%representation(), &
     & i=1,                                                  &
     & size(TYPES_SubspaceState)                             )])) then
    TYPES_SubspaceState = [TYPES_SubspaceState, SubspaceStatePointer(this)]
  endif
end subroutine

subroutine startup_SubspaceStates(this)
  implicit none
  
  class(SubspaceStates), intent(in) :: this
  
  integer :: i
  
  if (.not.allocated(TYPES_SubspaceStates)) then
    TYPES_SubspaceStates = [SubspaceStatesPointer(this)]
  elseif (.not.any([(                                          &
     & this%representation()                                   &
     &    == TYPES_SubspaceStates(i)%states_%representation(), &
     & i=1,                                                    &
     & size(TYPES_SubspaceStates)                              )])) then
    TYPES_SubspaceStates = [TYPES_SubspaceStates, SubspaceStatesPointer(this)]
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
impure elemental function initial_states_SubspaceBasisPointer(this,subspace, &
   & anharmonic_data) result(output)
  implicit none
  
  class(SubspaceBasisPointer), intent(in) :: this
  type(DegenerateSubspace),    intent(in) :: subspace
  type(AnharmonicData),        intent(in) :: anharmonic_data
  type(SubspaceStatesPointer)             :: output
  
  call this%check()
  
  output = this%basis_%initial_states(subspace, anharmonic_data)
end function

impure elemental function calculate_states_SubspaceBasisPointer(this, &
   & subspace,subspace_potential,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceBasisPointer), intent(in) :: this
  type(DegenerateSubspace),    intent(in) :: subspace
  class(PotentialData),        intent(in) :: subspace_potential
  type(AnharmonicData),        intent(in) :: anharmonic_data
  type(SubspaceStatesPointer)             :: output
  
  call this%check()
  
  output = this%basis_%calculate_states( subspace,           &
                                       & subspace_potential, &
                                       & anharmonic_data     )
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
       & TYPES_Subspacebasis(i)%basis_%representation()==representation, &
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
! SubspaceStatePointer methods.
! ----------------------------------------------------------------------
! Construct a SubspaceStatePointer from any type which extends SubspaceState.
impure elemental function new_SubspaceStatePointer(state) result(this)
  implicit none
  
  class(SubspaceState), intent(in) :: state
  type(SubspaceStatePointer)       :: this
  
  integer :: ialloc
  
  select type(state); type is(SubspaceStatePointer)
    this = state
  class default
    this%representation_ = state%representation()
    allocate( this%state_, source=state, &
            & stat=ialloc); call err(ialloc)
  end select
end function

! Checks that the pointer has been allocated before it is used.
subroutine check_SubspaceStatePointer(this)
  implicit none
  
  class(SubspaceStatePointer), intent(in) :: this
  
  if (.not. allocated(this%state_)) then
    call print_line(CODE_ERROR//': Trying to use a SubspaceStatePointer &
       &before it has been allocated.')
    call err()
  endif
end subroutine

! Type representation.
impure elemental function representation_SubspaceStatePointer() &
   & result(output)
  implicit none
  
  type(String) :: output
  
  output = 'pointer'
end function

! SubspaceState methods.
impure elemental function braket_SubspaceState_SubspaceStatePointer(this, &
   & ket,subspace,subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceStatePointer), intent(in)           :: this
  class(SubspaceState),        intent(in), optional :: ket
  type(DegenerateSubspace),    intent(in)           :: subspace
  class(SubspaceBasis),        intent(in)           :: subspace_basis
  type(AnharmonicData),        intent(in)           :: anharmonic_data
  real(dp)                                          :: output
  
  call this%check()
  
  if (present(ket)) then
    select type(ket); type is(SubspaceStatePointer)
      call ket%check()
      output = this%state_%braket( ket%state_,     &
                                 & subspace,       &
                                 & subspace_basis, &
                                 & anharmonic_data )
    class default
      output = this%state_%braket( ket,            &
                                 & subspace,       &
                                 & subspace_basis, &
                                 & anharmonic_data )
    end select
  else
    output = this%state_%braket( subspace        = subspace,       &
                               & subspace_basis  = subspace_basis, &
                               & anharmonic_data = anharmonic_data )
  endif
end function

impure elemental function braket_ComplexUnivariate_SubspaceStatePointer( &
   & this,univariate,ket,subspace,subspace_basis,anharmonic_data)        &
   & result(output)
  implicit none
  
  class(SubspaceStatePointer), intent(in)           :: this
  type(ComplexUnivariate),     intent(in)           :: univariate
  class(SubspaceState),        intent(in), optional :: ket
  type(DegenerateSubspace),    intent(in)           :: subspace
  class(SubspaceBasis),        intent(in)           :: subspace_basis
  type(AnharmonicData),        intent(in)           :: anharmonic_data
  type(ComplexMonomial)                             :: output
  
  call this%check()
  
  if (present(ket)) then
    select type(ket); type is(SubspaceStatePointer)
      call ket%check()
      output = this%state_%braket( univariate,     &
                                 & ket%state_,     &
                                 & subspace,       &
                                 & subspace_basis, &
                                 & anharmonic_data )
    class default
      output = this%state_%braket( univariate,     &
                                 & ket,            &
                                 & subspace,       &
                                 & subspace_basis, &
                                 & anharmonic_data )
    end select
  else
    output = this%state_%braket( univariate      = univariate,     &
                               & subspace        = subspace,       &
                               & subspace_basis  = subspace_basis, &
                               & anharmonic_data = anharmonic_data )
  endif
end function

impure elemental function braket_ComplexMonomial_SubspaceStatePointer(this, &
   & monomial,ket,subspace,subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceStatePointer), intent(in)           :: this
  type(ComplexMonomial),       intent(in)           :: monomial
  class(SubspaceState),        intent(in), optional :: ket
  type(DegenerateSubspace),    intent(in)           :: subspace
  class(SubspaceBasis),        intent(in)           :: subspace_basis
  type(AnharmonicData),        intent(in)           :: anharmonic_data
  type(ComplexMonomial)                             :: output
  
  call this%check()
  
  if (present(ket)) then
    select type(ket); type is(SubspaceStatePointer)
      call ket%check()
      output = this%state_%braket( monomial,       &
                                 & ket%state_,     &
                                 & subspace,       &
                                 & subspace_basis, &
                                 & anharmonic_data )
    class default
      output = this%state_%braket( monomial,       &
                                 & ket,            &
                                 & subspace,       &
                                 & subspace_basis, &
                                 & anharmonic_data )
    end select
  else
    output = this%state_%braket( monomial        = monomial,       &
                               & subspace        = subspace,       &
                               & subspace_basis  = subspace_basis, &
                               & anharmonic_data = anharmonic_data )
  endif
end function

impure elemental function kinetic_energy_SubspaceStatePointer(this,ket, &
   & subspace,subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceStatePointer), intent(in)           :: this
  class(SubspaceState),        intent(in), optional :: ket
  type(DegenerateSubspace),    intent(in)           :: subspace
  class(SubspaceBasis),        intent(in)           :: subspace_basis
  type(AnharmonicData),        intent(in)           :: anharmonic_data
  real(dp)                                          :: output
  
  call this%check()
  
  if (present(ket)) then
    select type(ket); type is(SubspaceStatePointer)
      call ket%check()
      output = this%state_%kinetic_energy( ket%state_,     &
                                         & subspace,       &
                                         & subspace_basis, &
                                         & anharmonic_data )
    class default
      output = this%state_%kinetic_energy( ket,            &
                                         & subspace,       &
                                         & subspace_basis, &
                                         & anharmonic_data )
    end select
  else
    output = this%state_%kinetic_energy( subspace        = subspace,       &
                                       & subspace_basis  = subspace_basis, &
                                       & anharmonic_data = anharmonic_data )
  endif
end function

impure elemental function harmonic_potential_energy_SubspaceStatePointer( &
   & this,ket,subspace,subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceStatePointer), intent(in)           :: this
  class(SubspaceState),        intent(in), optional :: ket
  type(DegenerateSubspace),    intent(in)           :: subspace
  class(SubspaceBasis),        intent(in)           :: subspace_basis
  type(AnharmonicData),        intent(in)           :: anharmonic_data
  real(dp)                                          :: output
  
  call this%check()
  
  if (present(ket)) then
    select type(ket); type is(SubspaceStatePointer)
      call ket%check()
      output = this%state_%harmonic_potential_energy( ket%state_,     &
                                                    & subspace,       &
                                                    & subspace_basis, &
                                                    & anharmonic_data )
    class default
      output = this%state_%harmonic_potential_energy( ket,            &
                                                    & subspace,       &
                                                    & subspace_basis, &
                                                    & anharmonic_data )
    end select
  else
    output = this%state_%harmonic_potential_energy( &
                & subspace        = subspace,       &
                & subspace_basis  = subspace_basis, &
                & anharmonic_data = anharmonic_data )
  endif
end function

! I/O.
subroutine read_SubspaceStatePointer(this,input)
  implicit none
  
  class(SubspaceStatePointer), intent(out) :: this
  type(String),                intent(in)  :: input(:)
  
  type(String), allocatable :: line(:)
  
  type(String) :: representation
  
  integer :: i
  
  select type(this); type is(SubspaceStatePointer)
    line = split_line(input(1))
    representation = line(3)
    
    ! Identify which type corresponds to the representation.
    i = first([(                                                         &
       & TYPES_SubspaceState(i)%state_%representation()==representation, &
       & i=1,                                                            &
       & size(TYPES_SubspaceState)                                       )])
    
    ! Read the input into the element of the correct type,
    !    and copy that element into the output.
    call TYPES_SubspaceState(i)%state_%read(input(2:))
    this = SubspaceStatePointer(TYPES_SubspaceState(i))
  class default
    call err()
  end select
end subroutine

function write_SubspaceStatePointer(this) result(output)
  implicit none
  
  class(SubspaceStatePointer), intent(in) :: this
  type(String), allocatable               :: output(:)
  
  select type(this); type is(SubspaceStatePointer)
    output = [ 'SubspaceState representation: '//this%representation_, &
             & str(this%state_)                                        ]
  end select
end function

function new_SubspaceStatePointer_Strings(input) result(this)
  implicit none
  
  type(String), intent(in)   :: input(:)
  type(SubspaceStatePointer) :: this
  
  call this%read(input)
end function

impure elemental function new_SubspaceStatePointer_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(SubspaceStatePointer)    :: this
  
  this = SubspaceStatePointer(str(input))
end function

! ----------------------------------------------------------------------
! SubspaceStatesPointer methods.
! ----------------------------------------------------------------------
! Construct a SubspaceStatesPointer from any type which extends SubspaceStates.
impure elemental function new_SubspaceStatesPointer(states) result(this)
  implicit none
  
  class(SubspaceStates), intent(in) :: states
  type(SubspaceStatesPointer)       :: this
  
  integer :: ialloc
  
  select type(states); type is(SubspaceStatesPointer)
    this = states
  class default
    this%representation_ = states%representation()
    allocate( this%states_, source=states, &
            & stat=ialloc); call err(ialloc)
  end select
end function

! Checks that the pointer has been allocated before it is used.
subroutine check_SubspaceStatesPointer(this)
  implicit none
  
  class(SubspaceStatesPointer), intent(in) :: this
  
  if (.not. allocated(this%states_)) then
    call print_line(CODE_ERROR//': Trying to use a SubspaceStatesPointer &
       &before it has been allocated.')
    call err()
  endif
end subroutine

! Type representation.
impure elemental function representation_SubspaceStatesPointer() result(output)
  implicit none
  
  type(String) :: output
  
  output = 'pointer'
end function

! SubspaceStates methods.
function states_SubspaceStatesPointer(this,subspace,subspace_basis, &
   & anharmonic_data) result(output)
  implicit none
  
  class(SubspaceStatesPointer), intent(in) :: this
  type(DegenerateSubspace),     intent(in) :: subspace
  class(SubspaceBasis),         intent(in) :: subspace_basis
  type(AnharmonicData),         intent(in) :: anharmonic_data
  type(SubspaceStatePointer), allocatable  :: output(:)
  
  call this%check()
  
 output = this%states_%states( subspace,       &
                             & subspace_basis, &
                             & anharmonic_data )
end function

impure elemental function spectra_SubspaceStatesPointer(this,subspace, &
   & subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceStatesPointer), intent(in) :: this
  type(DegenerateSubspace),     intent(in) :: subspace
  class(SubspaceBasis),         intent(in) :: subspace_basis
  type(AnharmonicData),         intent(in) :: anharmonic_data
  type(EnergySpectra)                      :: output
  
  call this%check()
  
  output = this%states_%spectra(subspace, subspace_basis, anharmonic_data)
end function

impure elemental function wavefunctions_SubspaceStatesPointer(this,subspace, &
   & subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceStatesPointer), intent(in) :: this
  type(DegenerateSubspace),     intent(in) :: subspace
  class(SubspaceBasis),         intent(in) :: subspace_basis
  type(AnharmonicData),         intent(in) :: anharmonic_data
  type(SubspaceWavefunctionsPointer)       :: output
  
  call this%check()
  
  output = this%states_%wavefunctions( subspace,       &
                                     & subspace_basis, &
                                     & anharmonic_data )
end function

impure elemental function integrate_SubspaceStatesPointer(this,potential, &
   & subspace,subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceStatesPointer), intent(in) :: this
  class(PotentialData),         intent(in) :: potential
  type(DegenerateSubspace),     intent(in) :: subspace
  class(SubspaceBasis),         intent(in) :: subspace_basis
  type(AnharmonicData),         intent(in) :: anharmonic_data
  type(PotentialPointer)                   :: output
  
  call this%check()
  
  output = this%states_%integrate( potential,      &
                                 & subspace,       &
                                 & subspace_basis, &
                                 & anharmonic_data )
end function

! I/O.
subroutine read_SubspaceStatesPointer(this,input)
  implicit none
  
  class(SubspaceStatesPointer), intent(out) :: this
  type(String),                 intent(in)  :: input(:)
  
  type(String), allocatable :: line(:)
  
  type(String) :: representation
  
  integer :: i
  
  select type(this); type is(SubspaceStatesPointer)
    line = split_line(input(1))
    representation = line(3)
    
    ! Identify which type corresponds to the representation.
    i = first([(                                                           &
       & TYPES_SubspaceStates(i)%states_%representation()==representation, &
       & i=1,                                                              &
       & size(TYPES_SubspaceStates)                                        )])
    
    ! Read the input into the element of the correct type,
    !    and copy that element into the output.
    call TYPES_SubspaceStates(i)%states_%read(input(2:))
    this = SubspaceStatesPointer(TYPES_SubspaceStates(i))
  class default
    call err()
  end select
end subroutine

function write_SubspaceStatesPointer(this) result(output)
  implicit none
  
  class(SubspaceStatesPointer), intent(in) :: this
  type(String), allocatable                :: output(:)
  
  select type(this); type is(SubspaceStatesPointer)
    output = [ 'SubspaceStates representation: '//this%representation_, &
             & str(this%states_)                                        ]
  end select
end function

function new_SubspaceStatesPointer_Strings(input) result(this)
  implicit none
  
  type(String), intent(in)    :: input(:)
  type(SubspaceStatesPointer) :: this
  
  call this%read(input)
end function

impure elemental function new_SubspaceStatesPointer_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(SubspaceStatesPointer)   :: this
  
  this = SubspaceStatesPointer(str(input))
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
impure elemental function potential_PotentialPointer(this) result(output)
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

subroutine generate_potential_PotentialPointer(this,anharmonic_data,   &
   & weighted_energy_force_ratio,calculate_stress,sampling_points_dir, &
   & calculation_reader,logfile)
  implicit none
  
  class(PotentialPointer), intent(inout) :: this
  type(AnharmonicData),    intent(in)    :: anharmonic_data
  real(dp),                intent(in)    :: weighted_energy_force_ratio
  logical,                 intent(in)    :: calculate_stress
  type(String),            intent(in)    :: sampling_points_dir
  type(CalculationReader), intent(inout) :: calculation_reader
  type(OFile),             intent(inout) :: logfile
  
  call this%check()
  
  call this%potential_%generate_potential( anharmonic_data,              &
                                         & weighted_energy_force_ratio,  &
                                         & calculate_stress,             &
                                         & sampling_points_dir,          &
                                         & calculation_reader,           &
                                         & logfile                       )
end subroutine

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

function braket_PotentialPointer(this,bra,ket,subspace,subspace_basis, &
   & anharmonic_data) result(output)
  implicit none
  
  class(PotentialPointer),  intent(in)           :: this
  class(SubspaceState),     intent(in)           :: bra
  class(SubspaceState),     intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(PotentialPointer)                         :: output
  
  call this%check()
  
  output = this%potential_%braket( bra,            &
                                 & ket,            &
                                 & subspace,       &
                                 & subspace_basis, &
                                 & anharmonic_data )
end function

function harmonic_expectation_PotentialPointer(this,frequency, &
   & thermal_energy,no_states,subspace,anharmonic_data) result(output)
  implicit none
  
  class(PotentialPointer),  intent(in) :: this
  real(dp),                 intent(in) :: frequency
  real(dp),                 intent(in) :: thermal_energy
  integer,                  intent(in) :: no_states
  type(DegenerateSubspace), intent(in) :: subspace
  type(AnharmonicData),     intent(in) :: anharmonic_data
  real(dp)                             :: output
  
  call this%check()
  
  output = this%potential_%harmonic_expectation( frequency,      &
                                               & thermal_energy, &
                                               & no_states,      &
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
  
  integer :: i
  
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
! Concrete methods of abstract classes.
! ----------------------------------------------------------------------
! SubspaceState methods.
impure elemental function braket_ComplexPolynomial_SubspaceState(this, &
   & polynomial,ket,subspace,subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceState),     intent(in)           :: this
  type(ComplexPolynomial),  intent(in)           :: polynomial
  class(SubspaceState),     intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(ComplexPolynomial)                        :: output
  
  type(ComplexMonomial), allocatable :: monomials(:)
  
  monomials = this%braket( polynomial%terms, &
                         & ket,              &
                         & subspace,         &
                         & subspace_basis,   &
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
! Takes a potential V and an array of subspace states {|i>}, and generates the
!    set of single-subspace potentials {V_i}, defined by
!    V_i = (prod_{j/=i}<j|)V(prod_{j/=i}|j>).
! ----------------------------------------------------------------------
! The naive method of calculating {V_i} for n subspaces takes
!    n(n-1) operations.
! This can be accelerated using a bisection method, outlined below.
! 
! V0(1) = V, the input potential.
! 
! The first iteration splits the states into two intervals,
!    [1,s-1] and [s,n], where s=n/2, and two potentials are calculated:
!    - V1(1) = (<s|<s+1|...<n|)V0(1)(|s>|s+1>...|n>)
!    - V1(s) = (<1|<2|...<s-1|)V0(1)(|1>|2>...|s-1>)
! These intervals are recorded in terms of their min and max values:
!   mins = [1  , s]
!   maxs = [s-1, n]
!
! The next iteration splits each of the intervals into two intervals,
!    copies the potential to both intervals, and integrates the potential
!    corresponding to each interval over the states in the other interval.
! This method takes O(n.log(n)) operations.
function generate_subspace_potentials(potential,subspaces,subspace_bases, &
   & subspace_states,anharmonic_data) result(output)
  implicit none
  
  class(PotentialData),     intent(in) :: potential
  type(DegenerateSubspace), intent(in) :: subspaces(:)
  class(SubspaceBasis),     intent(in) :: subspace_bases(:)
  class(SubspaceStates),    intent(in) :: subspace_states(:)
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(PotentialPointer), allocatable  :: output(:)
  
  ! The minimum and maximum indices in each interval.
  integer, allocatable :: mins_in(:)
  integer, allocatable :: maxs_in(:)
  integer, allocatable :: mins_out(:)
  integer, allocatable :: maxs_out(:)
  
  ! The half-way point of each interval.
  integer :: s
  
  integer :: i,j,ialloc
  
  if (size(subspace_states)==0) then
    call print_line(ERROR//': No states to integrate over.')
    call err()
  endif
  
  ! Initialise to a single interval spanning all states,
  !    and the un-integrated potential.
  mins_in = [1]
  maxs_in = [size(subspace_states)]
  allocate(output(size(subspace_states)), stat=ialloc); call err(ialloc)
  output(1) = PotentialPointer(potential)
  
  ! Loop over iteration until every interval contains exactly one subspace.
  do while (any(mins_in/=maxs_in))
    mins_out = [integer::]
    maxs_out = [integer::]
    
    ! Loop over intervals.
    do i=1,size(mins_in)
      if (mins_in(i)==maxs_in(i)) then
        ! The interval contains only one subspace; nothing needs doing.
        mins_out = [mins_out, mins_in(i)]
        maxs_out = [maxs_out, maxs_in(i)]
      else
        ! The interval contains more than one subspace;
        !    and generate two integrated potentials.
        
        ! Bisect the interval (or split it as evenly as possible).
        s = mins_in(i) + (maxs_in(i)-mins_in(i)+1)/2
        mins_out = [mins_out, mins_in(i), s     ]
        maxs_out = [maxs_out, s-1   , maxs_in(i)]
        
        ! Copy the potential from the first interval to the second.
        output(s) = output(mins_in(i))
        
        ! Integrate the first potential over all states in the second interval,
        !    and the second potential over all states in the first interval.
        do j=s,maxs_in(i)
          output(mins_in(i)) = PotentialPointer(                   &
             & subspace_states(j)%integrate( output(mins_in(i)),   &
             &                               subspaces(j),         &
             &                               subspace_bases(j),    &
             &                               anharmonic_data     ) )
        enddo
        do j=mins_in(i),s-1
          output(s) = PotentialPointer(                           &
             & subspace_states(j)%integrate( output(s),           &
             &                               subspaces(j),        &
             &                               subspace_bases(j),   &
             &                               anharmonic_data    ) )
        enddo
      endif
    enddo
    
    mins_in = mins_out
    maxs_in = maxs_out
  enddo
  
  ! Set the constant term to zero.
  do i=1,size(output)
    call output(i)%zero_energy()
  enddo
end function

! ----------------------------------------------------------------------
! The braket function.
! ----------------------------------------------------------------------
! Calculates <state|state>.
recursive function braket_state(state,subspace,subspace_basis, &
   & anharmonic_data) result(output)
  implicit none
  
  class(SubspaceState),     intent(in) :: state
  type(DegenerateSubspace), intent(in) :: subspace
  class(SubspaceBasis),     intent(in) :: subspace_basis
  type(AnharmonicData),     intent(in) :: anharmonic_data
  real(dp)                             :: output
  
  output = state%braket( subspace        = subspace,       &
                       & subspace_basis  = subspace_basis, &
                       & anharmonic_data = anharmonic_data )
end function

! Calculates <bra|ket>.
recursive function braket_state_state(bra,ket,subspace, &
   & subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceState),     intent(in) :: bra
  class(SubspaceState),     intent(in) :: ket
  type(DegenerateSubspace), intent(in) :: subspace
  class(SubspaceBasis),     intent(in) :: subspace_basis
  type(AnharmonicData),     intent(in) :: anharmonic_data
  real(dp)                             :: output
  
  output = bra%braket( ket             = ket,            &
                     & subspace        = subspace,       &
                     & subspace_basis  = subspace_basis, &
                     & anharmonic_data = anharmonic_data )
end function

! Calculates <state|V|state>.
recursive function braket_state_ComplexUnivariate(state,potential, &
   & subspace,subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceState),     intent(in) :: state
  type(ComplexUnivariate),  intent(in) :: potential
  type(DegenerateSubspace), intent(in) :: subspace
  class(SubspaceBasis),     intent(in) :: subspace_basis
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(ComplexMonomial)                :: output
  
  output = state%braket( univariate      = potential,      &
                       & subspace        = subspace,       &
                       & subspace_basis  = subspace_basis, &
                       & anharmonic_data = anharmonic_data )
end function

! Calculates <bra|V|ket>.
recursive function braket_state_ComplexUnivariate_state(bra, &
   & potential,ket,subspace,subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceState),     intent(in) :: bra
  type(ComplexUnivariate),  intent(in) :: potential
  class(SubspaceState),     intent(in) :: ket
  type(DegenerateSubspace), intent(in) :: subspace
  class(SubspaceBasis),     intent(in) :: subspace_basis
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(ComplexMonomial)                :: output
  
  output = bra%braket( univariate      = potential,      &
                     & ket             = ket,            &
                     & subspace        = subspace,       &
                     & subspace_basis  = subspace_basis, &
                     & anharmonic_data = anharmonic_data )
end function

! Calculates <state|V|state>.
recursive function braket_state_ComplexMonomial(state,potential, &
   & subspace,subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceState),     intent(in) :: state
  type(ComplexMonomial),    intent(in) :: potential
  type(DegenerateSubspace), intent(in) :: subspace
  class(SubspaceBasis),     intent(in) :: subspace_basis
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(ComplexMonomial)                :: output
  
  output = state%braket( monomial        = potential,      &
                       & subspace        = subspace,       &
                       & subspace_basis  = subspace_basis, &
                       & anharmonic_data = anharmonic_data )
end function

! Calculates <bra|V|ket>.
recursive function braket_state_ComplexMonomial_state(bra, &
   & potential,ket,subspace,subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceState),     intent(in) :: bra
  type(ComplexMonomial),    intent(in) :: potential
  class(SubspaceState),     intent(in) :: ket
  type(DegenerateSubspace), intent(in) :: subspace
  class(SubspaceBasis),     intent(in) :: subspace_basis
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(ComplexMonomial)                :: output
  
  output = bra%braket( monomial        = potential,      &
                     & ket             = ket,            &
                     & subspace        = subspace,       &
                     & subspace_basis  = subspace_basis, &
                     & anharmonic_data = anharmonic_data )
end function

! Calculates <state|V|state>.
recursive function braket_state_ComplexPolynomial(state,potential, &
   & subspace,subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceState),     intent(in) :: state
  type(ComplexPolynomial),  intent(in) :: potential
  type(DegenerateSubspace), intent(in) :: subspace
  class(SubspaceBasis),     intent(in) :: subspace_basis
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(ComplexPolynomial)              :: output
  
  output = state%braket( polynomial      = potential,      &
                       & subspace        = subspace,       &
                       & subspace_basis  = subspace_basis, &
                       & anharmonic_data = anharmonic_data )
end function

! Calculates <bra|V|ket>.
recursive function braket_state_ComplexPolynomial_state(bra, &
   & potential,ket,subspace,subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceState),     intent(in) :: bra
  type(ComplexPolynomial),  intent(in) :: potential
  class(SubspaceState),     intent(in) :: ket
  type(DegenerateSubspace), intent(in) :: subspace
  class(SubspaceBasis),     intent(in) :: subspace_basis
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(ComplexPolynomial)              :: output
  
  output = bra%braket( polynomial      = potential,      &
                     & ket             = ket,            &
                     & subspace        = subspace,       &
                     & subspace_basis  = subspace_basis, &
                     & anharmonic_data = anharmonic_data )
end function

! Calculates <state|V|state>.
recursive function braket_state_potential(state,potential,subspace, &
   & subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceState),     intent(in) :: state
  class(PotentialData),     intent(in) :: potential
  type(DegenerateSubspace), intent(in) :: subspace
  class(SubspaceBasis),     intent(in) :: subspace_basis
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(PotentialPointer)               :: output
  
  output = potential%braket( bra             = state,          &
                           & subspace        = subspace,       &
                           & subspace_basis  = subspace_basis, &
                           & anharmonic_data = anharmonic_data )
end function

! Calculates <bra|V|ket>.
recursive function braket_state_potential_state(bra,potential,ket, &
   & subspace,subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceState),     intent(in) :: bra
  class(PotentialData),     intent(in) :: potential
  class(SubspaceState),     intent(in) :: ket
  type(DegenerateSubspace), intent(in) :: subspace
  class(SubspaceBasis),     intent(in) :: subspace_basis
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(PotentialPointer)               :: output
  
  output = potential%braket( bra             = bra,            &
                           & ket             = ket,            &
                           & subspace        = subspace,       &
                           & subspace_basis  = subspace_basis, &
                           & anharmonic_data = anharmonic_data )
end function

! Calculates <state|T|state>.
recursive function kinetic_energy_state(state,subspace, &
   & subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceState),     intent(in) :: state
  type(DegenerateSubspace), intent(in) :: subspace
  class(SubspaceBasis),     intent(in) :: subspace_basis
  type(AnharmonicData),     intent(in) :: anharmonic_data
  real(dp)                             :: output
  
  output = state%kinetic_energy( subspace        = subspace,       &
                               & subspace_basis  = subspace_basis, &
                               & anharmonic_data = anharmonic_data )
end function

! Calculates <bra|T|ket>.
recursive function kinetic_energy_state_state(bra,ket,subspace, &
   & subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceState),     intent(in) :: bra
  class(SubspaceState),     intent(in) :: ket
  type(DegenerateSubspace), intent(in) :: subspace
  class(SubspaceBasis),     intent(in) :: subspace_basis
  type(AnharmonicData),     intent(in) :: anharmonic_data
  real(dp)                             :: output
  
  output = bra%kinetic_energy( ket             = ket,            &
                             & subspace        = subspace,       &
                             & subspace_basis  = subspace_basis, &
                             & anharmonic_data = anharmonic_data )
end function

! Calculates <state|Vh|state>, where Vh is the harmonic potential.
recursive function harmonic_potential_energy_state(state,subspace, &
   & subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceState),     intent(in) :: state
  type(DegenerateSubspace), intent(in) :: subspace
  class(SubspaceBasis),     intent(in) :: subspace_basis
  type(AnharmonicData),     intent(in) :: anharmonic_data
  real(dp)                             :: output
  
  output = state%harmonic_potential_energy( &
        & subspace        = subspace,       &
        & subspace_basis  = subspace_basis, &
        & anharmonic_data = anharmonic_data )
end function

! Calculates <bra|Vh|ket>, where Vh is the harmonic potential.
recursive function harmonic_potential_energy_state_state(bra,ket, &
   & subspace,subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceState),     intent(in) :: bra
  class(SubspaceState),     intent(in) :: ket
  type(DegenerateSubspace), intent(in) :: subspace
  class(SubspaceBasis),     intent(in) :: subspace_basis
  type(AnharmonicData),     intent(in) :: anharmonic_data
  real(dp)                             :: output
  
  output = bra%harmonic_potential_energy( ket             = ket,            &
                                        & subspace        = subspace,       &
                                        & subspace_basis  = subspace_basis, &
                                        & anharmonic_data = anharmonic_data )
end function

! Calculates <state|V|state> as a constant.
recursive function potential_energy_state(state,potential,subspace, &
   & subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceState),     intent(in) :: state
  class(PotentialData),     intent(in) :: potential
  type(DegenerateSubspace), intent(in) :: subspace
  class(SubspaceBasis),     intent(in) :: subspace_basis
  type(AnharmonicData),     intent(in) :: anharmonic_data
  real(dp)                             :: output
  
  type(PotentialPointer), allocatable :: integrated_potential
  
  integrated_potential = braket( state,          &
                               & potential,      &
                               & subspace,       &
                               & subspace_basis, &
                               & anharmonic_data )
  output = integrated_potential%undisplaced_energy()
end function

! Calcualtes <bra|V|ket> as a constant.
recursive function potential_energy_state_state(bra,potential,ket,subspace, &
   & subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceState),     intent(in) :: bra
  class(PotentialData),     intent(in) :: potential
  class(SubspaceState),     intent(in) :: ket
  type(DegenerateSubspace), intent(in) :: subspace
  class(SubspaceBasis),     intent(in) :: subspace_basis
  type(AnharmonicData),     intent(in) :: anharmonic_data
  real(dp)                             :: output
  
  type(PotentialPointer), allocatable :: integrated_potential
  
  integrated_potential = braket( bra,            &
                               & potential,      &
                               & ket,            &
                               & subspace,       &
                               & subspace_basis, &
                               & anharmonic_data )
  output = integrated_potential%undisplaced_energy()
end function
end module
