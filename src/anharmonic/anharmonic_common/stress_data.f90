! ======================================================================
! Defines the StressData type and the StressPointer.
! This type extends StressBase with methods for constructing
!    and interpolating the stress.
! ======================================================================
module caesar_stress_data_module
  use caesar_common_module
  
  use caesar_subspace_coupling_module
  use caesar_anharmonic_data_module
  use caesar_subspace_state_module
  use caesar_subspace_braket_module
  use caesar_basis_state_module
  use caesar_basis_states_module
  use caesar_abstract_classes_module
  implicit none
  
  public :: StressData
  public :: StressPointer
  
  type, abstract, extends(StressBase) :: StressData
  contains
    procedure(representation_StressData), public, deferred, nopass :: &
       & representation
    procedure, public :: startup => startup_StressData
    
    ! Interpolation of the stress.
    procedure, public :: can_be_interpolated => &
                       & can_be_interpolated_StressData
    procedure, public :: calculate_dynamical_matrices => &
                       & calculate_dynamical_matrices_StressData
    procedure, public :: stress_correction => &
                       & stress_correction_StressData
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
    
    procedure, public :: braket_SubspaceBraKet  => &
                       & braket_SubspaceBraKet_StressPointer
    procedure, public :: braket_BasisState  => &
                       & braket_BasisState_StressPointer
    procedure, public :: braket_BasisStates => &
                       & braket_BasisStates_StressPointer
    
    procedure, public :: harmonic_expectation => &
                       & harmonic_expectation_StressPointer
    
    procedure, public :: potential_stress_SubspaceBraKet => &
                       & potential_stress_SubspaceBraKet_StressPointer
    procedure, public :: potential_stress_BasisState => &
                       & potential_stress_BasisState_StressPointer
    
    procedure, public :: can_be_interpolated => &
                       & can_be_interpolated_StressPointer
    procedure, public :: calculate_dynamical_matrices => &
                       & calculate_dynamical_matrices_StressPointer
    procedure, public :: stress_correction => &
                       & stress_correction_StressPointer
    
    ! I/O.
    procedure, public :: read  => read_StressPointer
    procedure, public :: write => write_StressPointer
  end type
  
  abstract interface
    ! StressData procedures.
    impure elemental function representation_StressData() result(output)
      import String
      implicit none
      
      type(String) :: output
    end function
  end interface
  
  interface
    module subroutine startup_StressData(this) 
      class(StressData), intent(in) :: this
    end subroutine
  end interface
  
  interface StressPointer
    ! Construct a StressPointer from any type which extends StressData.
    impure elemental module function new_StressPointer(stress) result(this) 
      class(StressData), intent(in) :: stress
      type(StressPointer)           :: this
    end function
  end interface
  
  interface
    ! Checks that the pointer has been allocated before it is used.
    impure elemental module subroutine check_StressPointer(this) 
      class(StressPointer), intent(in) :: this
    end subroutine
  end interface
  
  interface
    ! Type representation.
    impure elemental module function representation_StressPointer() &
       & result(output) 
      type(String) :: output
    end function
  end interface
  
  interface
    impure elemental module subroutine zero_stress_StressPointer(this) 
      class(StressPointer), intent(inout) :: this
    end subroutine
  end interface
  
  interface
    impure elemental module subroutine add_constant_StressPointer(this,input) 
      class(StressPointer), intent(inout) :: this
      type(RealMatrix),     intent(in)    :: input
    end subroutine
  end interface
  
  interface
    impure elemental module function stress_RealModeDisplacement_StressPointer(this,displacement) result(output) 
      class(StressPointer),       intent(in) :: this
      type(RealModeDisplacement), intent(in) :: displacement
      type(RealMatrix)                       :: output
    end function
  end interface
  
  interface
    impure elemental module function stress_ComplexModeDisplacement_StressPointer(   this,displacement) result(output) 
      class(StressPointer),          intent(in) :: this
      type(ComplexModeDisplacement), intent(in) :: displacement
      type(ComplexMatrix)                       :: output
    end function
  end interface
  
  interface
    impure elemental module subroutine braket_SubspaceBraKet_StressPointer(this,braket,whole_subspace,anharmonic_data) 
      class(StressPointer),  intent(inout)        :: this
      class(SubspaceBraKet), intent(in)           :: braket
      logical,               intent(in), optional :: whole_subspace
      type(AnharmonicData),  intent(in)           :: anharmonic_data
    end subroutine
  end interface
  
  interface
    impure elemental module subroutine braket_BasisState_StressPointer(this, &
       & bra,ket,subspace,subspace_basis,whole_subspace,anharmonic_data) 
      class(StressPointer),     intent(inout)        :: this
      class(BasisState),        intent(in)           :: bra
      class(BasisState),        intent(in), optional :: ket
      type(DegenerateSubspace), intent(in)           :: subspace
      class(SubspaceBasis),     intent(in)           :: subspace_basis
      logical,                  intent(in), optional :: whole_subspace
      type(AnharmonicData),     intent(in)           :: anharmonic_data
    end subroutine
  end interface
  
  interface
    impure elemental module subroutine braket_BasisStates_StressPointer(this,states,subspace,subspace_basis,whole_subspace,anharmonic_data) 
      class(StressPointer),     intent(inout)        :: this
      class(BasisStates),       intent(inout)        :: states
      type(DegenerateSubspace), intent(in)           :: subspace
      class(SubspaceBasis),     intent(in)           :: subspace_basis
      logical,                  intent(in), optional :: whole_subspace
      type(AnharmonicData),     intent(in)           :: anharmonic_data
    end subroutine
  end interface
  
  interface
    impure elemental module function harmonic_expectation_StressPointer(this,frequency,thermal_energy,supercell_size,anharmonic_data) result(output) 
      class(StressPointer), intent(in) :: this
      real(dp),             intent(in) :: frequency
      real(dp),             intent(in) :: thermal_energy
      integer,              intent(in) :: supercell_size
      type(AnharmonicData), intent(in) :: anharmonic_data
      type(RealMatrix)                 :: output
    end function
  end interface
  
  interface
    recursive module function potential_stress_SubspaceBraKet_StressPointer(this,braket,anharmonic_data) result(output) 
      class(StressPointer),  intent(in) :: this
      class(SubspaceBraKet), intent(in) :: braket
      type(AnharmonicData),  intent(in) :: anharmonic_data
      type(RealMatrix)                  :: output
    end function
  end interface
  
  interface
    recursive module function potential_stress_BasisState_StressPointer(this,bra,ket,subspace,subspace_basis,anharmonic_data) result(output) 
      class(StressPointer),     intent(in)           :: this
      class(BasisState),        intent(in)           :: bra
      class(BasisState),        intent(in), optional :: ket
      type(DegenerateSubspace), intent(in)           :: subspace
      class(SubspaceBasis),     intent(in)           :: subspace_basis
      type(AnharmonicData),     intent(in)           :: anharmonic_data
      type(RealMatrix)                               :: output
    end function
  end interface
  
  interface
    module function can_be_interpolated_StressPointer(this) result(output) 
      class(StressPointer), intent(in) :: this
      logical                          :: output
    end function
  end interface
  
  interface
    module function calculate_dynamical_matrices_StressPointer(this,qpoints, &
       & thermal_energy,subspaces,subspace_bases,subspace_states,            &
       & anharmonic_data) result(output) 
      class(StressPointer),     intent(in)     :: this
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
    module function stress_correction_StressPointer(this,subspaces, &
       & subspace_bases,subspace_states,anharmonic_data) result(output) 
      class(StressPointer),     intent(in)    :: this
      type(DegenerateSubspace), intent(in)    :: subspaces(:)
      class(SubspaceBasis),     intent(in)    :: subspace_bases(:)
      class(BasisStates),       intent(inout) :: subspace_states(:)
      type(AnharmonicData),     intent(in)    :: anharmonic_data
      type(RealMatrix)                        :: output
    end function
  end interface
  
  interface
    module function can_be_interpolated_StressData(this) result(output) 
      class(StressData), intent(in) :: this
      logical                       :: output
    end function
  end interface
  
  interface
    module function calculate_dynamical_matrices_StressData(this,qpoints, &
       & thermal_energy,subspaces,subspace_bases,subspace_states,         &
       & anharmonic_data) result(output) 
      class(StressData),        intent(in)     :: this
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
    module function stress_correction_StressData(this,subspaces, &
       & subspace_bases,subspace_states,anharmonic_data) result(output) 
      class(StressData),        intent(in)    :: this
      type(DegenerateSubspace), intent(in)    :: subspaces(:)
      class(SubspaceBasis),     intent(in)    :: subspace_bases(:)
      class(BasisStates),       intent(inout) :: subspace_states(:)
      type(AnharmonicData),     intent(in)    :: anharmonic_data
      type(RealMatrix)                        :: output
    end function
  end interface
  
  interface
    ! I/O.
    module subroutine read_StressPointer(this,input) 
      class(StressPointer), intent(out) :: this
      type(String),         intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_StressPointer(this) result(output) 
      class(StressPointer), intent(in) :: this
      type(String), allocatable        :: output(:)
    end function
  end interface
  
  interface StressPointer
    module function new_StressPointer_Strings(input) result(this) 
      type(String), intent(in) :: input(:)
      type(StressPointer)      :: this
    end function
  
    impure elemental module function new_StressPointer_StringArray(input) &
       & result(this) 
      type(StringArray), intent(in) :: input
      type(StressPointer)           :: this
    end function
  end interface
end module
