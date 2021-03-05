! ======================================================================
! A basis of states which treats q-points separately.
! ======================================================================
! N.B. only states at a single q-point are calculated and stored.
! All other q-points are included by symmetry.
module caesar_split_qpoints_basis_module
  use caesar_common_module
  
  use caesar_anharmonic_common_module
  
  use caesar_wavevector_state_module
  use caesar_wavevector_states_module
  use caesar_wavevector_basis_module
  use caesar_split_qpoints_wavefunctions_module
  use caesar_core_shell_thermodynamics_module
  implicit none
  
  private
  
  public :: QpointModeIDs
  public :: SplitQpointsBasis
  
  type, extends(Stringable) :: QpointModeIDs
    integer, allocatable :: mode_ids(:)
    integer, allocatable :: paired_mode_ids(:)
  contains
    procedure, public :: read  => read_QpointModeIDs
    procedure, public :: write => write_QpointModeIDs
  end type
  
  ! All states spanning the subspace.
  type, extends(SubspaceBasis) :: SplitQpointsBasis
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
    ! The q-points, and the modes at each q-point.
    type(QpointModeIDs), allocatable :: qpoints(:)
    
    ! The q-points which should be integrated.
    integer, allocatable, private :: qpoints_to_integrate_(:)
  contains
    procedure, public, nopass :: representation => &
                               & representation_SplitQpointsBasis
    
    ! Print the harmonic ground-state wavefunction of the basis.
    procedure, public :: ground_state_wavefunction
    
    ! Generate the first guess at states.
    procedure, public :: initial_states => initial_states_SplitQpointsBasis
    
    ! Generate the eigenstates of a single-subspace potential.
    procedure, public :: calculate_states => calculate_states_SplitQpointsBasis
    
    ! Process a subspace potential or stress.
    ! For a SplitQpointsBasis, this converts the single-subspace object
    !    into a single-q-point object.
    procedure, public :: process_subspace_potential => &
                       & process_subspace_potential_SplitQpointsBasis
    procedure, public :: process_subspace_stress => &
                       & process_subspace_stress_SplitQpointsBasis
    
    ! Return the modes spanned by this basis.
    procedure, public :: mode_ids => mode_ids_SplitQpointsBasis
    procedure, public :: paired_mode_ids => paired_mode_ids_SplitQpointsBasis
    
    ! Procedures involving individual states.
    procedure, public :: inner_product => &
                       & inner_product_SplitQpointsBasis
    procedure, public :: integrate_BasisState => &
                       & integrate_BasisState_SplitQpointsBasis
    procedure, public :: kinetic_energy => &
                       & kinetic_energy_SplitQpointsBasis
    procedure, public :: harmonic_potential_energy => &
                       & harmonic_potential_energy_SplitQpointsBasis
    procedure, public :: kinetic_stress => &
                       & kinetic_stress_SplitQpointsBasis
    procedure, public :: wavefunction => wavefunction_SplitQpointsBasis
    
    ! Procedures involving sets of states.
    procedure, public :: thermodynamic_data => &
                       & thermodynamic_data_SplitQpointsBasis
    procedure, public :: wavefunctions => wavefunctions_SplitQpointsBasis
    procedure, public :: integrate_BasisStates => &
                       & integrate_BasisStates_SplitQpointsBasis
    
    ! The derivative of the free energy.
    procedure, public :: free_energy_gradient => &
                       & free_energy_gradient_SplitQpointsBasis
    
    ! Select symmetries which map the principle q-point to itself or its pair.
    procedure, public :: select_symmetries => &
                       & select_symmetries_SplitQpointsBasis
    
    ! I/O.
    procedure, public :: read  => read_SplitQpointsBasis
    procedure, public :: write => write_SplitQpointsBasis
    
    ! Private helper functions.
    !procedure, private :: integrate_potential
  end type
  
  interface QpointModeIDs
    ! QpointModeIDs methods.
    module function new_QpointModeIDs(mode_ids,paired_mode_ids) result(this) 
      integer, intent(in) :: mode_ids(:)
      integer, intent(in) :: paired_mode_ids(:)
      type(QpointModeIDs) :: this
    end function
  end interface
  
  interface SplitQpointsBasis
    ! ----------------------------------------------------------------------
    ! SplitQpointsBasis methods.
    ! ----------------------------------------------------------------------
    ! Constructors.
    module function new_SplitQpointsBasis(supercell_size,maximum_power, &
       & expansion_order,subspace_id,frequency,wavevectors,qpoints)     &
       & result(this) 
      integer,               intent(in) :: supercell_size
      integer,               intent(in) :: maximum_power
      integer,               intent(in) :: expansion_order
      integer,               intent(in) :: subspace_id
      real(dp),              intent(in) :: frequency
      type(WavevectorBasis), intent(in) :: wavevectors(:)
      type(QpointModeIDs),   intent(in) :: qpoints(:)
      type(SplitQpointsBasis)           :: this
    end function
  end interface
  
  interface
    ! Type representation.
    impure elemental module function representation_SplitQpointsBasis() &
       & result(output) 
      type(String) :: output
    end function
  end interface
  
  interface SplitQpointsBasis
    ! Generates states up to a given power, spanning the whole subspace.
    module function new_SplitQpointsBasis_subspace(subspace,frequency,modes, &
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
      type(SplitQpointsBasis)              :: output
    end function
  end interface
  
  interface
    ! Returns the harmonic ground-state wavefunction for the basis.
    module function ground_state_wavefunction(this,subspace,supercell) &
       & result(output) 
      class(SplitQpointsBasis), intent(in) :: this
      type(DegenerateSubspace), intent(in) :: subspace
      type(StructureData),      intent(in) :: supercell
      type(String)                         :: output
    end function
  end interface
  
  interface
    ! Generate an initial guess at states.
    impure elemental module function initial_states_SplitQpointsBasis(this, &
       & subspace,thermal_energy,anharmonic_data) result(output) 
      class(SplitQpointsBasis), intent(in) :: this
      type(DegenerateSubspace), intent(in) :: subspace
      real(dp),                 intent(in) :: thermal_energy
      type(AnharmonicData),     intent(in) :: anharmonic_data
      type(BasisStatesPointer)             :: output
    end function
  end interface
  
  interface
    ! Calculate the eigenstates of a single-subspace potential.
    impure elemental module function calculate_states_SplitQpointsBasis(this,subspace,subspace_potential,thermal_energy,state_energy_cutoff,convergence_data,anharmonic_data) result(output) 
      class(SplitQpointsBasis), intent(in) :: this
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
    ! Integrates the potential over all but the first q-point in the subspace.
    impure elemental module subroutine process_subspace_potential_SplitQpointsBasis(   this,potential,states,subspace,anharmonic_data) 
      class(SplitQpointsBasis), intent(in)            :: this
      class(PotentialBase),     intent(inout)         :: potential
      class(BasisStates),       intent(inout), target :: states
      type(DegenerateSubspace), intent(in)            :: subspace
      type(AnharmonicData),     intent(in)            :: anharmonic_data
    end subroutine
  end interface
  
  interface
    impure elemental module subroutine process_subspace_stress_SplitQpointsBasis(this,stress,states,subspace,anharmonic_data) 
      class(SplitQpointsBasis), intent(in)            :: this
      class(StressBase),        intent(inout)         :: stress
      class(BasisStates),       intent(inout), target :: states
      type(DegenerateSubspace), intent(in)            :: subspace
      type(AnharmonicData),     intent(in)            :: anharmonic_data
    end subroutine
  end interface
  
  interface
    !function integrate_potential(this,potential,states,thermal_energy,subspace,
    !   & anharmonic_data) result(output)
    !  implicit none
    !  
    !  class(SplitQpointsBasis), intent(in)    :: this
    !  class(PotentialBase),     intent(in)    :: potential
    !  type(WavevectorStates),   intent(inout) :: states
    !  real(dp),                 intent(in)    :: thermal_energy
    !  type(DegenerateSubspace), intent(in)    :: subspace
    !  type(AnharmonicData),     intent(in)    :: anharmonic_data
    !  type(PotentialBasePointer)              :: output
    !  
    !  type(SplitQpointsBasis) :: basis
    !  
    !  type(PotentialBasePointer) :: integrated_potential
    !  real(dp)                   :: correction_energy
    !  
    !  integer :: i
    !  
    !  basis = this
    !  
    !  ! Calculate (prod_{j/=i}<j|)V(prod_{j/=i}|j>).
    !  output = PotentialBasePointer(potential)
    !  
    !  if (size(this%qpoints)==1) then
    !    return
    !  endif
    !  
    !  basis%qpoints_to_integrate_ = [(i,i=2,size(this%qpoints))]
    !  
    !  call output%braket( states,                           &
    !                    & subspace        = subspace,       &
    !                    & thermal_energy  = thermal_energy, &
    !                    & subspace_basis  = basis,          &
    !                    & whole_subspace  = .false.,        &
    !                    & anharmonic_data = anharmonic_data )
    !  
    !  ! Calculate (prod_i<i|)V(prod_i|i>), and use this to correct for
    !  !    over-counting. See main vscf method for details.
    !  integrated_potential = output
    !  
    !  basis%qpoints_to_integrate_ = [1]
    !  
    !  call integrated_potential%braket( states,                           &
    !                                  & subspace        = subspace,       &
    !                                  & thermal_energy  = thermal_energy, &
    !                                  & subspace_basis  = basis,          &
    !                                  & anharmonic_data = anharmonic_data )
    !  
    !  correction_energy = integrated_potential%undisplaced_energy() &
    !                  & * (1.0_dp-size(this%qpoints))/size(this%qpoints)
    !  
    !  call output%add_constant(correction_energy)
    !end function
    
    impure elemental module function wavefunction_SplitQpointsBasis(this, &
       & state,supercell) result(output) 
      class(SplitQpointsBasis), intent(in) :: this
      type(WavevectorState),    intent(in) :: state
      type(StructureData),      intent(in) :: supercell
      type(String)                         :: output
    end function
  end interface
  
  interface
    module function mode_ids_SplitQpointsBasis(this,subspace,anharmonic_data) &
       & result(output) 
      class(SplitQpointsBasis), intent(in) :: this
      type(DegenerateSubspace), intent(in) :: subspace
      type(AnharmonicData),     intent(in) :: anharmonic_data
      integer, allocatable                 :: output(:)
    end function
  end interface
  
  interface
    module function paired_mode_ids_SplitQpointsBasis(this,subspace, &
       & anharmonic_data) result(output) 
      class(SplitQpointsBasis), intent(in) :: this
      type(DegenerateSubspace), intent(in) :: subspace
      type(AnharmonicData),     intent(in) :: anharmonic_data
      integer, allocatable                 :: output(:)
    end function
  end interface
  
  interface
    impure elemental module function inner_product_SplitQpointsBasis(this, &
       & bra,ket,subspace,anharmonic_data) result(output) 
      class(SplitQpointsBasis), intent(in)                   :: this
      class(BasisState),        intent(in),           target :: bra
      class(BasisState),        intent(in), optional, target :: ket
      type(DegenerateSubspace), intent(in)                   :: subspace
      type(AnharmonicData),     intent(in)                   :: anharmonic_data
      real(dp)                                               :: output
    end function
  end interface
  
  interface
    impure elemental module function integrate_BasisState_SplitQpointsBasis(this,bra,monomial,ket,subspace,anharmonic_data) result(output) 
      class(SplitQpointsBasis), intent(in)                   :: this
      class(BasisState),        intent(in),           target :: bra
      type(SparseMonomial),     intent(in)                   :: monomial
      class(BasisState),        intent(in), optional, target :: ket
      type(DegenerateSubspace), intent(in)                   :: subspace
      type(AnharmonicData),     intent(in)                   :: anharmonic_data
      complex(dp)                                            :: output
    end function
  end interface
  
  interface
    ! Helper function. Takes a monomial, and splits into an array of monomials,
    !    where the output monomials correspond to each of this%qpoints in turn.
    ! Transforms each monomial so it appears to be at the first q-point.
    ! Integrating this monomial at the first q-point is equivalent to
    !    integrating the input monomial in the basis of states at the
    !    selected q-point, but is much more efficient.
    module function split_monomial(this,monomial) result(output) 
      class(SplitQpointsBasis), intent(in) :: this
      type(SparseMonomial),     intent(in) :: monomial
      type(SparseMonomial), allocatable    :: output(:)
    end function
  end interface
  
  interface
    impure elemental module function kinetic_energy_SplitQpointsBasis(this, &
       & bra,ket,subspace,anharmonic_data) result(output) 
      class(SplitQpointsBasis), intent(in)                   :: this
      class(BasisState),        intent(in),           target :: bra
      class(BasisState),        intent(in), optional, target :: ket
      type(DegenerateSubspace), intent(in)                   :: subspace
      type(AnharmonicData),     intent(in)                   :: anharmonic_data
      real(dp)                                               :: output
    end function
  end interface
  
  interface
    impure elemental module function harmonic_potential_energy_SplitQpointsBasis(this,bra,ket,subspace,anharmonic_data) result(output) 
      class(SplitQpointsBasis), intent(in)                   :: this
      class(BasisState),        intent(in),           target :: bra
      class(BasisState),        intent(in), optional, target :: ket
      type(DegenerateSubspace), intent(in)                   :: subspace
      type(AnharmonicData),     intent(in)                   :: anharmonic_data
      real(dp)                                               :: output
    end function
  end interface
  
  interface
    impure elemental module function kinetic_stress_SplitQpointsBasis(this, &
       & bra,ket,subspace,stress_prefactors,anharmonic_data) result(output) 
      class(SplitQpointsBasis), intent(in)                   :: this
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
    impure elemental module function thermodynamic_data_SplitQpointsBasis(this,thermal_energy,states,subspace,subspace_potential,subspace_stress,stress_prefactors,anharmonic_data) result(output) 
      class(SplitQpointsBasis), intent(in)                  :: this
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
    impure elemental module function wavefunctions_SplitQpointsBasis(this, &
       & states,subspace,anharmonic_data) result(output) 
      class(SplitQpointsBasis),  intent(in)         :: this
      class(BasisStates),        intent(in), target :: states
      type(DegenerateSubspace),  intent(in)         :: subspace
      type(AnharmonicData),      intent(in)         :: anharmonic_data
      type(SubspaceWavefunctionsPointer) :: output
    end function
  end interface
  
  interface
    ! Integrate a monomial.
    impure elemental module function integrate_BasisStates_SplitQpointsBasis(this,states,monomial,subspace,anharmonic_data) result(output) 
      class(SplitQpointsBasis),  intent(in)         :: this
      class(BasisStates),        intent(in), target :: states
      type(SparseMonomial),      intent(in)         :: monomial
      type(DegenerateSubspace),  intent(in)         :: subspace
      type(AnharmonicData),      intent(in)         :: anharmonic_data
      complex(dp)                                   :: output
    end function
  end interface
  
  interface
    ! Calculate the derivative of the free energy.
    module function free_energy_gradient_SplitQpointsBasis(this,            &
        subspace_potential,basis_functions,subspace,states,thermal_energy, &
          & state_energy_cutoff,anharmonic_data) result(output) 
      class(SplitQpointsBasis), intent(in) :: this
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
    ! Select the symmetries which map the modes at the principle q-point (and its
    !    pair) onto other such modes.
    module function select_symmetries_SplitQpointsBasis(this,symmetries, &
       & anharmonic_data) result(output) 
      class(SplitQpointsBasis), intent(in)  :: this
      type(DegenerateSymmetry), intent(in)  :: symmetries(:)
      type(AnharmonicData),     intent(in)  :: anharmonic_data
      type(DegenerateSymmetry), allocatable :: output(:)
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_QpointModeIDs(this,input) 
      class(QpointModeIDs), intent(out) :: this
      type(String),         intent(in)  :: input
    end subroutine
  end interface
  
  interface
    module function write_QpointModeIDs(this) result(output) 
      class(QpointModeIDs), intent(in) :: this
      type(String)                         :: output
    end function
  end interface
  
  interface QpointModeIDs
    impure elemental module function new_QpointModeIDs_String(input) &
       & result(this) 
      type(String), intent(in) :: input
      type(QpointModeIDs)      :: this
    end function
  end interface
  
  interface
    module subroutine read_SplitQpointsBasis(this,input) 
      class(SplitQpointsBasis), intent(out) :: this
      type(String),             intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_SplitQpointsBasis(this) result(output) 
      class(SplitQpointsBasis), intent(in) :: this
      type(String), allocatable            :: output(:)
    end function
  end interface
  
  interface SplitQpointsBasis
    module function new_SplitQpointsBasis_Strings(input) result(this) 
      type(String), intent(in) :: input(:)
      type(SplitQpointsBasis)  :: this
    end function
  
    impure elemental module function new_SplitQpointsBasis_StringArray(input) &
       & result(this) 
      type(StringArray), intent(in) :: input
      type(SplitQpointsBasis)       :: this
    end function
  end interface
end module
