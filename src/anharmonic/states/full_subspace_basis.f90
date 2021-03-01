! ======================================================================
! A basis of states which spans a full subspace.
! ======================================================================
module caesar_full_subspace_basis_module
  use caesar_common_module
  
  use caesar_anharmonic_common_module
  
  use caesar_wavevector_state_module
  use caesar_wavevector_states_module
  use caesar_wavevector_basis_module
  use caesar_full_subspace_wavefunctions_module
  use caesar_core_shell_thermodynamics_module
  implicit none
  
  private
  
  public :: startup_full_subspace_basis
  
  public :: FullSubspaceBasis
  
  ! All states spanning the subspace.
  type, extends(SubspaceBasis) :: FullSubspaceBasis
    ! The number of primitive cells in the supercell.
    integer :: supercell_size
    ! The maximum power of the monomial states.
    ! This is also the maximum occupation of the harmonic basis states.
    integer  :: maximum_power
    ! The expansion order of the potential.
    ! This is also the limit on coupling between states.
    integer :: expansion_order
    ! The ID of the subspace.
    integer  :: subspace_id
    ! The states, wavevector by wavevector.
    ! N.B. this only includes one wavevector from each symmetry-related set.
    type(WavevectorBasis), allocatable :: wavevectors(:)
  contains
    procedure, public, nopass :: representation => &
                               & representation_FullSubspaceBasis
    
    ! Print the harmonic ground-state wavefunction of the basis.
    procedure, public :: ground_state_wavefunction
    
    ! Generate the first guess at states.
    procedure, public :: initial_states => initial_states_FullSubspaceBasis
    
    ! Generate the eigenstates of a single-subspace potential.
    procedure, public :: calculate_states => calculate_states_FullSubspaceBasis
    
    ! Return the modes spanned by this basis.
    procedure, public :: mode_ids => mode_ids_FullSubspaceBasis
    procedure, public :: paired_mode_ids => paired_mode_ids_FullSubspaceBasis
    
    ! Procedures involving individual states.
    procedure, public :: inner_product => &
                       & inner_product_FullSubspaceBasis
    procedure, public :: integrate_BasisState => &
                       & integrate_BasisState_FullSubspaceBasis
    procedure, public :: kinetic_energy => &
                       & kinetic_energy_FullSubspaceBasis
    procedure, public :: harmonic_potential_energy => &
                       & harmonic_potential_energy_FullSubspaceBasis
    procedure, public :: kinetic_stress => &
                       & kinetic_stress_FullSubspaceBasis
    procedure, public :: wavefunction => wavefunction_FullSubspaceBasis
    
    ! Procedures involving sets of states.
    procedure, public :: thermodynamic_data => &
                       & thermodynamic_data_FullSubspaceBasis
    procedure, public :: wavefunctions => wavefunctions_FullSubspaceBasis
    procedure, public :: integrate_BasisStates => &
                       & integrate_BasisStates_FullSubspaceBasis
    
    ! The derivative of the free energy.
    procedure, public :: free_energy_gradient => &
                       & free_energy_gradient_FullSubspaceBasis
    
    ! I/O.
    procedure, public :: read  => read_FullSubspaceBasis
    procedure, public :: write => write_FullSubspaceBasis
  end type
  
  interface
    ! Startup procedure.
    module subroutine startup_full_subspace_basis() 
    end subroutine
  end interface
  
  interface FullSubspaceBasis
    ! ----------------------------------------------------------------------
    ! FullSubspaceBasis methods.
    ! ----------------------------------------------------------------------
    ! Constructors.
    module function new_FullSubspaceBasis(supercell_size,maximum_power, &
       & expansion_order,subspace_id,frequency,wavevectors) result(this) 
      integer,               intent(in) :: supercell_size
      integer,               intent(in) :: maximum_power
      integer,               intent(in) :: expansion_order
      integer,               intent(in) :: subspace_id
      real(dp),              intent(in) :: frequency
      type(WavevectorBasis), intent(in) :: wavevectors(:)
      type(FullSubspaceBasis)           :: this
    end function
  end interface
  
  interface
    ! Type representation.
    impure elemental module function representation_FullSubspaceBasis() &
       & result(output) 
      type(String) :: output
    end function
  end interface
  
  interface FullSubspaceBasis
    ! Generates states up to a given power, spanning the whole subspace.
    module function new_FullSubspaceBasis_subspace(subspace,frequency,modes, &
       & qpoints,supercell,maximum_power,potential_expansion_order,          &
       & symmetries) result(output) 
      type(DegenerateSubspace), intent(in) :: subspace
      real(dp),                 intent(in) :: frequency
      type(ComplexMode),        intent(in) :: modes(:)
      type(QpointData),         intent(in) :: qpoints(:)
      type(StructureData),      intent(in) :: supercell
      integer,                  intent(in) :: maximum_power
      integer,                  intent(in) :: potential_expansion_order
      type(SymmetryOperator),   intent(in) :: symmetries(:)
      type(FullSubspaceBasis)              :: output
    end function
  end interface
  
  interface
    ! Returns the harmonic ground-state wavefunction for the basis.
    module function ground_state_wavefunction(this,subspace,supercell) &
       & result(output) 
      class(FullSubspaceBasis), intent(in) :: this
      type(DegenerateSubspace), intent(in) :: subspace
      type(StructureData),      intent(in) :: supercell
      type(String)                         :: output
    end function
  end interface
  
  interface
    ! Generate an initial guess at states.
    impure elemental module function initial_states_FullSubspaceBasis(this, &
       & subspace,thermal_energy,anharmonic_data) result(output) 
      class(FullSubspaceBasis), intent(in) :: this
      type(DegenerateSubspace), intent(in) :: subspace
      real(dp),                 intent(in) :: thermal_energy
      type(AnharmonicData),     intent(in) :: anharmonic_data
      type(BasisStatesPointer)             :: output
    end function
  end interface
  
  interface
    ! Calculate the eigenstates of a single-subspace potential.
    impure elemental module function calculate_states_FullSubspaceBasis(this,subspace,subspace_potential,thermal_energy,state_energy_cutoff,convergence_data,anharmonic_data) result(output) 
      class(FullSubspaceBasis), intent(in) :: this
      type(DegenerateSubspace), intent(in) :: subspace
      class(PotentialBase),     intent(in) :: subspace_potential
      real(dp),                 intent(in) :: thermal_energy
      real(dp),                 intent(in) :: state_energy_cutoff
      type(ConvergenceData),    intent(in) :: convergence_data
      type(AnharmonicData),     intent(in) :: anharmonic_data
      type(BasisStatesPointer)             :: output
    end function
  end interface
  
  interface
    impure elemental module function wavefunction_FullSubspaceBasis(this, &
       & state,supercell) result(output) 
      class(FullSubspaceBasis), intent(in) :: this
      type(WavevectorState),    intent(in) :: state
      type(StructureData),      intent(in) :: supercell
      type(String)                         :: output
    end function
  end interface
  
  interface
    module function mode_ids_FullSubspaceBasis(this,subspace,anharmonic_data) &
       & result(output) 
      class(FullSubspaceBasis), intent(in) :: this
      type(DegenerateSubspace), intent(in) :: subspace
      type(AnharmonicData),     intent(in) :: anharmonic_data
      integer, allocatable                 :: output(:)
    end function
  end interface
  
  interface
    module function paired_mode_ids_FullSubspaceBasis(this,subspace, &
       & anharmonic_data) result(output) 
      class(FullSubspaceBasis), intent(in) :: this
      type(DegenerateSubspace), intent(in) :: subspace
      type(AnharmonicData),     intent(in) :: anharmonic_data
      integer, allocatable                 :: output(:)
    end function
  end interface
  
  interface
    impure elemental module function inner_product_FullSubspaceBasis(this, &
       & bra,ket,subspace,anharmonic_data) result(output) 
      class(FullSubspaceBasis), intent(in)                   :: this
      class(BasisState),        intent(in),           target :: bra
      class(BasisState),        intent(in), optional, target :: ket
      type(DegenerateSubspace), intent(in)                   :: subspace
      type(AnharmonicData),     intent(in)                   :: anharmonic_data
      real(dp)                                               :: output
    end function
  end interface
  
  interface
    impure elemental module function integrate_BasisState_FullSubspaceBasis(this,bra,monomial,ket,subspace,anharmonic_data) result(output) 
      class(FullSubspaceBasis), intent(in)                   :: this
      class(BasisState),        intent(in),           target :: bra
      type(SparseMonomial),     intent(in)                   :: monomial
      class(BasisState),        intent(in), optional, target :: ket
      type(DegenerateSubspace), intent(in)                   :: subspace
      type(AnharmonicData),     intent(in)                   :: anharmonic_data
      complex(dp)                                            :: output
    end function
  end interface
  
  interface
    impure elemental module function kinetic_energy_FullSubspaceBasis(this, &
       & bra,ket,subspace,anharmonic_data) result(output) 
      class(FullSubspaceBasis), intent(in)                   :: this
      class(BasisState),        intent(in),           target :: bra
      class(BasisState),        intent(in), optional, target :: ket
      type(DegenerateSubspace), intent(in)                   :: subspace
      type(AnharmonicData),     intent(in)                   :: anharmonic_data
      real(dp)                                               :: output
    end function
  end interface
  
  interface
    impure elemental module function harmonic_potential_energy_FullSubspaceBasis(   this,bra,ket,subspace,anharmonic_data) result(output) 
      class(FullSubspaceBasis), intent(in)                   :: this
      class(BasisState),        intent(in),           target :: bra
      class(BasisState),        intent(in), optional, target :: ket
      type(DegenerateSubspace), intent(in)                   :: subspace
      type(AnharmonicData),     intent(in)                   :: anharmonic_data
      real(dp)                                               :: output
    end function
  end interface
  
  interface
    impure elemental module function kinetic_stress_FullSubspaceBasis(this, &
       & bra,ket,subspace,stress_prefactors,anharmonic_data) result(output) 
      class(FullSubspaceBasis), intent(in)                   :: this
      class(BasisState),        intent(in),           target :: bra
      class(BasisState),        intent(in), optional, target :: ket
      type(DegenerateSubspace), intent(in)                   :: subspace
      type(StressPrefactors),   intent(in)                   :: stress_prefactors
      type(AnharmonicData),     intent(in)                   :: anharmonic_data
      type(RealMatrix)                                       :: output
    end function
  end interface
  
  interface
    ! Thermodynamic data. Energy, entropy, free energy etc.
    impure elemental module function thermodynamic_data_FullSubspaceBasis(this,thermal_energy,states,subspace,subspace_potential,subspace_stress,stress_prefactors,anharmonic_data) result(output) 
      class(FullSubspaceBasis), intent(in)                  :: this
      real(dp),                 intent(in)                  :: thermal_energy
      class(BasisStates),       intent(in),          target :: states
      type(DegenerateSubspace), intent(in)                  :: subspace
      class(PotentialBase),     intent(in)                  :: subspace_potential
      class(StressBase),        intent(in), optional        :: subspace_stress
      type(StressPrefactors),   intent(in), optional        :: stress_prefactors
      type(AnharmonicData),     intent(in)                  :: anharmonic_data
      type(ThermodynamicData)                               :: output
    end function
  end interface
  
  interface
    ! Wavefunctions.
    impure elemental module function wavefunctions_FullSubspaceBasis(this, &
       & states,subspace,anharmonic_data) result(output) 
      class(FullSubspaceBasis),  intent(in)         :: this
      class(BasisStates),        intent(in), target :: states
      type(DegenerateSubspace),  intent(in)         :: subspace
      type(AnharmonicData),      intent(in)         :: anharmonic_data
      type(SubspaceWavefunctionsPointer) :: output
    end function
  end interface
  
  interface
    ! Integrate a monomial.
    impure elemental module function integrate_BasisStates_FullSubspaceBasis(this,states,monomial,subspace,anharmonic_data) result(output) 
      class(FullSubspaceBasis), intent(in)         :: this
      class(BasisStates),       intent(in), target :: states
      type(SparseMonomial),     intent(in)         :: monomial
      type(DegenerateSubspace), intent(in)         :: subspace
      type(AnharmonicData),     intent(in)         :: anharmonic_data
      complex(dp)                                  :: output
    end function
  end interface
  
  interface
    ! Calculate the derivative of the free energy.
    module function free_energy_gradient_FullSubspaceBasis(this,            &
        subspace_potential,basis_functions,subspace,states,thermal_energy, &
          & state_energy_cutoff,anharmonic_data) result(output) 
      class(FullSubspaceBasis), intent(in) :: this
      class(PotentialBase),     intent(in) :: subspace_potential
      class(PotentialBase) ,intent(in) :: basis_functions(:) 
      type(DegenerateSubspace), intent(in) :: subspace
      class(BasisStates),       intent(in) :: states
      real(dp),                 intent(in) :: thermal_energy
      real(dp),                 intent(in) :: state_energy_cutoff
      type(AnharmonicData),     intent(in) :: anharmonic_data
      real(dp), allocatable                :: output(:)
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_FullSubspaceBasis(this,input) 
      class(FullSubspaceBasis), intent(out) :: this
      type(String),             intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_FullSubspaceBasis(this) result(output) 
      class(FullSubspaceBasis), intent(in) :: this
      type(String), allocatable            :: output(:)
    end function
  end interface
  
  interface FullSubspaceBasis
    module function new_FullSubspaceBasis_Strings(input) result(this) 
      type(String), intent(in) :: input(:)
      type(FullSubspaceBasis)  :: this
    end function
  
    impure elemental module function new_FullSubspaceBasis_StringArray(input) &
       & result(this) 
      type(StringArray), intent(in) :: input
      type(FullSubspaceBasis)       :: this
    end function
  end interface
end module
