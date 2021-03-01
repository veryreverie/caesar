! ======================================================================
! A basis harmonic states for a given subspace at a given wavevector.
! N.B. the wavevector in this context is the wavevector of the state,
!    not the wavevector of the mode.
! e.g. the state |p> along a mode at q-point q has a wavevector of p*q.
! ======================================================================
! N.B. The WavevectorBasis has a polymorphic allocatable component,
!   harmonic_states_. Since compiler support for such things is not fantastic,
!   care should be taken when modifying the code referring to this component.
module caesar_wavevector_basis_module
  use caesar_common_module
  
  use caesar_anharmonic_common_module
  
  use caesar_harmonic_state_1d_module
  use caesar_harmonic_state_2d_module
  use caesar_harmonic_state_real_module
  use caesar_harmonic_state_complex_module
  use caesar_harmonic_braket_real_module
  use caesar_harmonic_braket_complex_module
  use caesar_coupled_states_module
  use caesar_density_matrix_module
  use caesar_wavevector_state_module
  use caesar_wavevector_states_module
  use caesar_calculate_weights_module
  implicit none
  
  private
  
  public :: startup_wavevector_basis
  
  public :: WavevectorBasis
  
  public :: calculate_states
  public :: core_harmonic_observables
  public :: core_effective_harmonic_observables
  public :: core_vci_observables
  public :: free_energy_gradient
  
  public :: HarmonicBraKetReal
  public :: HarmonicBraKetComplex
  
  type, extends(SubspaceBasis) :: WavevectorBasis
    integer                                    :: maximum_power
    integer                                    :: expansion_order
    integer                                    :: subspace_id
    type(FractionVector)                       :: wavevector
    class(SubspaceState), allocatable, private :: harmonic_states_(:)
    type(CoupledStates),  allocatable, private :: harmonic_couplings_(:)
  contains
    procedure, public, nopass :: representation => &
                               & representation_WavevectorBasis
    
    ! State procedures.
    procedure, public :: inner_product => &
                       & inner_product_WavevectorBasis
    procedure, public :: integrate_BasisState => &
                       & integrate_BasisState_WavevectorBasis
    procedure, public :: kinetic_energy => &
                       & kinetic_energy_WavevectorBasis
    procedure, public :: harmonic_potential_energy => &
                       & harmonic_potential_energy_WavevectorBasis
    procedure, public :: kinetic_stress => &
                       & kinetic_stress_WavevectorBasis
    
    ! Basis procedures.
    procedure, public :: initial_ground_state => &
                       & initial_ground_state_WavevectorBasis
    procedure, public :: initial_states => &
                       & initial_states_WavevectorBasis
    procedure, public :: calculate_states => &
                       & calculate_states_WavevectorBasis
    
    ! Mode IDs.
    procedure, public :: mode_ids => mode_ids_WavevectorBasis
    procedure, public :: paired_mode_ids => paired_mode_ids_WavevectorBasis
    
    ! I/O.
    procedure, public :: read  => read_WavevectorBasis
    procedure, public :: write => write_WavevectorBasis
    
    ! Inherited procedures which do not apply to this type.
    procedure, public :: integrate_BasisStates => &
                       & integrate_BasisStates_WavevectorBasis
    procedure, public :: thermodynamic_data => &
                       & thermodynamic_data_WavevectorBasis
    procedure, public :: wavefunctions => &
                       & wavefunctions_WavevectorBasis
    procedure, public :: free_energy_gradient => &
                       & free_energy_gradient_WavevectorBasis
  end type
  
  type, extends(NoDefaultConstructor) :: SelectedStatesHamiltonian
    integer, allocatable   :: selected_states(:)
    type(SparseRealMatrix) :: hamiltonian
  end type
  
  interface
    ! Startup procedure.
    module subroutine startup_wavevector_basis() 
    end subroutine
  end interface
  
  interface
    ! Type representation.
    impure elemental module function representation_WavevectorBasis() &
       & result(output) 
      type(String) :: output
    end function
  end interface
  
  interface WavevectorBasis
    ! Constructor and size module function.
    module function new_WavevectorBasis(maximum_power,expansion_order,        &
       & subspace_id,frequency,wavevector,harmonic_states,harmonic_couplings) &
       & result(this) 
      integer,              intent(in) :: maximum_power
      integer,              intent(in) :: expansion_order
      integer,              intent(in) :: subspace_id
      real(dp),             intent(in) :: frequency
      type(FractionVector), intent(in) :: wavevector
      class(SubspaceState), intent(in) :: harmonic_states(:)
      type(CoupledStates),  intent(in) :: harmonic_couplings(:)
      type(WavevectorBasis)            :: this
    end function
  end interface
  
  interface HarmonicBraKetReal
    ! Construct a BraKet from a basis.
    module function new_HarmonicBraKetReal_WavevectorBasis(basis,state, &
       & anharmonic_data) result(this) 
      type(WavevectorBasis),   intent(in) :: basis
      type(HarmonicStateReal), intent(in) :: state
      type(AnharmonicData),    intent(in) :: anharmonic_data
      type(HarmonicBraKetReal)            :: this
    end function
  end interface
  
  interface HarmonicBraKetComplex
    module function new_HarmonicBraKetComplex_WavevectorBasis(basis,state, &
       & anharmonic_data) result(this) 
      type(WavevectorBasis),      intent(in) :: basis
      type(HarmonicStateComplex), intent(in) :: state
      type(AnharmonicData),       intent(in) :: anharmonic_data
      type(HarmonicBraKetComplex)            :: this
    end function
  end interface
  
  interface WavevectorBasis
    ! ----------------------------------------------------------------------
    ! Generates states in a given subspace, up to a given power.
    ! ----------------------------------------------------------------------
    ! If qpoint is specified, only modes at that q-point are included.
    module function new_WavevectorBasis_subspace(subspace,frequency,modes, &
       & qpoints,maximum_power,potential_expansion_order,supercell_size,   &
       & symmetries,qpoint) result(output) 
      type(DegenerateSubspace), intent(in)           :: subspace
      real(dp),                 intent(in)           :: frequency
      type(ComplexMode),        intent(in)           :: modes(:)
      type(QpointData),         intent(in)           :: qpoints(:)
      integer,                  intent(in)           :: maximum_power
      integer,                  intent(in)           :: potential_expansion_order
      integer,                  intent(in)           :: supercell_size
      type(SymmetryOperator),   intent(in)           :: symmetries(:)
      type(QpointData),         intent(in), optional :: qpoint
      type(WavevectorBasis), allocatable             :: output(:)
    end function
  end interface
  
  interface
    module function generate_mode_basis(subspace,frequency,mode, &
       & maximum_power,potential_expansion_order,supercell_size) &
       & result(output) 
      type(DegenerateSubspace), intent(in) :: subspace
      real(dp),                 intent(in) :: frequency
      type(ComplexMode),        intent(in) :: mode
      integer,                  intent(in) :: maximum_power
      integer,                  intent(in) :: potential_expansion_order
      integer,                  intent(in) :: supercell_size
      type(WavevectorBasis)                :: output
    end function
  end interface
  
  interface
    module function generate_mode_basis_1d(subspace,frequency,mode, &
       & maximum_power,potential_expansion_order,supercell_size)    &
       & result(output) 
      type(DegenerateSubspace), intent(in) :: subspace
      real(dp),                 intent(in) :: frequency
      type(ComplexMode),        intent(in) :: mode
      integer,                  intent(in) :: maximum_power
      integer,                  intent(in) :: potential_expansion_order
      integer,                  intent(in) :: supercell_size
      type(WavevectorBasis)                :: output
    end function
  end interface
  
  interface
    module function generate_mode_basis_2d(subspace,frequency,mode, &
       & maximum_power,potential_expansion_order,supercell_size)    &
       & result(output) 
      type(DegenerateSubspace), intent(in) :: subspace
      real(dp),                 intent(in) :: frequency
      type(ComplexMode),        intent(in) :: mode
      integer,                  intent(in) :: maximum_power
      integer,                  intent(in) :: potential_expansion_order
      integer,                  intent(in) :: supercell_size
      type(WavevectorBasis)                :: output
    end function
  end interface
  
  interface
    ! Takes the basis over one set of modes in the subspace, and the basis over a
    !    disjoint set of modes in the same subspace, and constructs the basis
    !    over the union of the two sets, up to the maximum power.
    ! e.g. if 'this has states |0>, |u1> and |u1^2>,
    !    and 'that' has states |0>, |u2> and |u2^2>, then the output will have
    !    states |0>, |u1>, |u2>, |u1^2>, |u1u2> and |u2^2>.
    ! N.B. the states of 'that' must be sorted in ascending order of total power.
    module function tensor_product(this,that) result(output) 
      type(WavevectorBasis), intent(in) :: this
      type(WavevectorBasis), intent(in) :: that
      type(WavevectorBasis)             :: output
    end function
  end interface
  
  interface
    ! Splits up a WavevectorBasis by wavevector.
    module function split_by_wavevector(input,modes,qpoints,symmetries) &
       & result(output) 
      type(WavevectorBasis),  intent(in) :: input
      type(ComplexMode),      intent(in) :: modes(:)
      type(QpointData),       intent(in) :: qpoints(:)
      type(SymmetryOperator), intent(in) :: symmetries(:)
      type(WavevectorBasis), allocatable :: output(:)
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Operations involving WavevectorStates.
    ! ----------------------------------------------------------------------
    impure elemental module function inner_product_WavevectorBasis(this,bra, &
       & ket,subspace,anharmonic_data) result(output) 
      class(WavevectorBasis),   intent(in)                   :: this
      class(BasisState),        intent(in),           target :: bra
      class(BasisState),        intent(in), optional, target :: ket
      type(DegenerateSubspace), intent(in)                   :: subspace
      type(AnharmonicData),     intent(in)                   :: anharmonic_data
      real(dp)                                               :: output
    end function
  end interface
  
  interface
    impure elemental module function integrate_BasisState_WavevectorBasis(this,bra,monomial,ket,subspace,anharmonic_data) result(output) 
      class(WavevectorBasis),   intent(in)                   :: this
      class(BasisState),        intent(in),           target :: bra
      type(SparseMonomial),     intent(in)                   :: monomial
      class(BasisState),        intent(in), optional, target :: ket
      type(DegenerateSubspace), intent(in)                   :: subspace
      type(AnharmonicData),     intent(in)                   :: anharmonic_data
      complex(dp)                                            :: output
    end function
  end interface
  
  interface
    module function calculate_integral_real(basis,harmonic_states,bra,ket, &
       & monomial,anharmonic_data) result(output) 
      type(WavevectorBasis),   intent(in)         :: basis
      type(HarmonicStateReal), intent(in), target :: harmonic_states(:)
      type(WavevectorState),   intent(in)         :: bra
      type(WavevectorState),   intent(in)         :: ket
      type(SparseMonomial),    intent(in)         :: monomial
      type(AnharmonicData),    intent(in)         :: anharmonic_data
      complex(dp)                                 :: output
    end function
  end interface
  
  interface
    module function calculate_integral_complex(basis,harmonic_states,bra, &
       & ket,monomial,anharmonic_data) result(output) 
      type(WavevectorBasis),      intent(in)         :: basis
      type(HarmonicStateComplex), intent(in), target :: harmonic_states(:)
      type(WavevectorState),      intent(in)         :: bra
      type(WavevectorState),      intent(in)         :: ket
      type(SparseMonomial),       intent(in)         :: monomial
      type(AnharmonicData),       intent(in)         :: anharmonic_data
      complex(dp)                                    :: output
    end function
  end interface
  
  interface
    ! Integrate a monomial between sets of states.
    impure elemental module function integrate_BasisStates_WavevectorBasis(this,states,monomial,subspace,anharmonic_data) result(output) 
      class(WavevectorBasis),   intent(in)         :: this
      class(BasisStates),       intent(in), target :: states
      type(SparseMonomial),     intent(in)         :: monomial
      type(DegenerateSubspace), intent(in)         :: subspace
      type(AnharmonicData),     intent(in)         :: anharmonic_data
      complex(dp)                                  :: output
    end function
  end interface
  
  interface
    module function integrate_real(basis,harmonic_states,density_matrix, &
       & monomial,anharmonic_data) result(output) 
      type(WavevectorBasis),   intent(in)         :: basis
      type(HarmonicStateReal), intent(in), target :: harmonic_states(:)
      type(DensityMatrix),     intent(in)         :: density_matrix
      type(SparseMonomial),    intent(in)         :: monomial
      type(AnharmonicData),    intent(in)         :: anharmonic_data
      complex(dp)                                 :: output
    end function
  end interface
  
  interface
    module function integrate_complex(basis,harmonic_states,density_matrix, &
       & monomial,anharmonic_data) result(output) 
      type(WavevectorBasis),      intent(in)         :: basis
      type(HarmonicStateComplex), intent(in), target :: harmonic_states(:)
      type(DensityMatrix),        intent(in)         :: density_matrix
      type(SparseMonomial),       intent(in)         :: monomial
      type(AnharmonicData),       intent(in)         :: anharmonic_data
      complex(dp)                                    :: output
    end function
  end interface
  
  interface
    impure elemental module function kinetic_energy_WavevectorBasis(this, &
       & bra,ket,subspace,anharmonic_data) result(output) 
      class(WavevectorBasis),   intent(in)                   :: this
      class(BasisState),        intent(in),           target :: bra
      class(BasisState),        intent(in), optional, target :: ket
      type(DegenerateSubspace), intent(in)                   :: subspace
      type(AnharmonicData),     intent(in)                   :: anharmonic_data
      real(dp)                                               :: output
    end function
  end interface
  
  interface
    module function kinetic_energy_real(basis,harmonic_states,bra,ket, &
       & anharmonic_data) result(output) 
      type(WavevectorBasis),   intent(in)         :: basis
      type(HarmonicStateReal), intent(in), target :: harmonic_states(:)
      type(WavevectorState),   intent(in)         :: bra
      type(WavevectorState),   intent(in)         :: ket
      type(AnharmonicData),    intent(in)         :: anharmonic_data
      real(dp)                                    :: output
    end function
  end interface
  
  interface
    module function kinetic_energy_complex(basis,harmonic_states,bra,ket, &
       & anharmonic_data) result(output) 
      type(WavevectorBasis),      intent(in)         :: basis
      type(HarmonicStateComplex), intent(in), target :: harmonic_states(:)
      type(WavevectorState),      intent(in)         :: bra
      type(WavevectorState),      intent(in)         :: ket
      type(AnharmonicData),       intent(in)         :: anharmonic_data
      real(dp)                                       :: output
    end function
  end interface
  
  interface
    impure elemental module function harmonic_potential_energy_WavevectorBasis(this,bra,ket,subspace,anharmonic_data) result(output) 
      class(WavevectorBasis),   intent(in)                   :: this
      class(BasisState),        intent(in),           target :: bra
      class(BasisState),        intent(in), optional, target :: ket
      type(DegenerateSubspace), intent(in)                   :: subspace
      type(AnharmonicData),     intent(in)                   :: anharmonic_data
      real(dp)                                               :: output
    end function
  end interface
  
  interface
    module function harmonic_potential_energy_real(basis,harmonic_states, &
       & bra,ket,anharmonic_data) result(output) 
      type(WavevectorBasis),   intent(in)         :: basis
      type(HarmonicStateReal), intent(in), target :: harmonic_states(:)
      type(WavevectorState),   intent(in)         :: bra
      type(WavevectorState),   intent(in)         :: ket
      type(AnharmonicData),    intent(in)         :: anharmonic_data
      real(dp)                                    :: output
    end function
  end interface
  
  interface
    module function harmonic_potential_energy_complex(basis,harmonic_states, &
       & bra,ket,anharmonic_data) result(output) 
      type(WavevectorBasis),      intent(in)         :: basis
      type(HarmonicStateComplex), intent(in), target :: harmonic_states(:)
      type(WavevectorState),      intent(in)         :: bra
      type(WavevectorState),      intent(in)         :: ket
      type(AnharmonicData),       intent(in)         :: anharmonic_data
      real(dp)                                       :: output
    end function
  end interface
  
  interface
    impure elemental module function kinetic_stress_WavevectorBasis(this, &
       & bra,ket,subspace,stress_prefactors,anharmonic_data) result(output) 
      class(WavevectorBasis),   intent(in)                   :: this
      class(BasisState),        intent(in),           target :: bra
      class(BasisState),        intent(in), optional, target :: ket
      type(DegenerateSubspace), intent(in)                   :: subspace
      type(StressPrefactors),   intent(in)                   :: stress_prefactors
      type(AnharmonicData),     intent(in)                   :: anharmonic_data
      type(RealMatrix)                                       :: output
    end function
  end interface
  
  interface
    module function kinetic_stress_real(basis,harmonic_states,bra,ket, &
       & stress_prefactors,anharmonic_data) result(output) 
      type(WavevectorBasis),   intent(in)         :: basis
      type(HarmonicStateReal), intent(in), target :: harmonic_states(:)
      type(WavevectorState),   intent(in)         :: bra
      type(WavevectorState),   intent(in)         :: ket
      type(StressPrefactors),  intent(in)         :: stress_prefactors
      type(AnharmonicData),    intent(in)         :: anharmonic_data
      type(RealMatrix)                            :: output
    end function
  end interface
  
  interface
    module function kinetic_stress_complex(basis,harmonic_states,bra,ket, &
       & stress_prefactors,anharmonic_data) result(output) 
      type(WavevectorBasis),      intent(in)         :: basis
      type(HarmonicStateComplex), intent(in), target :: harmonic_states(:)
      type(WavevectorState),      intent(in)         :: bra
      type(WavevectorState),      intent(in)         :: ket
      type(StressPrefactors),     intent(in)         :: stress_prefactors
      type(AnharmonicData),       intent(in)         :: anharmonic_data
      type(RealMatrix)                               :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Operations which generate WavevectorStates.
    ! ----------------------------------------------------------------------
    ! Generate initial guess. This is simply the basis state |0>, i.e. the
    !    ground state of the effective harmonic potential from which the basis
    !    states were generated.
    impure elemental module function initial_ground_state_WavevectorBasis(this) result(output) 
      class(WavevectorBasis), intent(in) :: this
      type(WavevectorState)              :: output
    end function
  end interface
  
  interface
    impure elemental module function initial_states_WavevectorBasis(this, &
       & subspace,thermal_energy,anharmonic_data) result(output) 
      class(WavevectorBasis),   intent(in) :: this
      type(DegenerateSubspace), intent(in) :: subspace
      real(dp),                 intent(in) :: thermal_energy
      type(AnharmonicData),     intent(in) :: anharmonic_data
      type(BasisStatesPointer)             :: output
    end function
  end interface
  
  interface
    impure elemental module function calculate_states_WavevectorBasis(this, &
       & subspace,subspace_potential,thermal_energy,state_energy_cutoff,    &
       & convergence_data,anharmonic_data) result(output) 
      class(WavevectorBasis),   intent(in) :: this
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
    module function calculate_hamiltonian_real(basis,harmonic_states,   &
       & state_energy_cutoff,subspace_potential,include_kinetic_energy, &
       & anharmonic_data) result(output) 
      type(WavevectorBasis),   intent(in)           :: basis
      type(HarmonicStateReal), intent(in), target   :: harmonic_states(:)
      real(dp),                intent(in)           :: state_energy_cutoff
      class(PotentialBase),    intent(in)           :: subspace_potential
      logical,                 intent(in), optional :: include_kinetic_energy
      type(AnharmonicData),    intent(in)           :: anharmonic_data
      type(SelectedStatesHamiltonian)               :: output
    end function
  end interface
  
  interface
    module function calculate_hamiltonian_complex(basis,harmonic_states, &
       & state_energy_cutoff,subspace_potential,include_kinetic_energy,  &
       & anharmonic_data) result(output) 
      type(WavevectorBasis),      intent(in)           :: basis
      type(HarmonicStateComplex), intent(in), target   :: harmonic_states(:)
      real(dp),                   intent(in)           :: state_energy_cutoff
      class(PotentialBase),       intent(in)           :: subspace_potential
      logical,                    intent(in), optional :: include_kinetic_energy
      type(AnharmonicData),       intent(in)           :: anharmonic_data
      type(SelectedStatesHamiltonian)                  :: output
    end function
  end interface
  
  interface
    ! Calculate the eigenstates of a wavevector basis.
    module function calculate_states(basis,subspace,subspace_potential,       &
       & thermal_energy,state_energy_cutoff,convergence_data,anharmonic_data) &
       & result(output) 
      type(WavevectorBasis),    intent(in) :: basis(:)
      type(DegenerateSubspace), intent(in) :: subspace
      class(PotentialBase),     intent(in) :: subspace_potential
      real(dp),                 intent(in) :: thermal_energy
      real(dp),                 intent(in) :: state_energy_cutoff
      type(ConvergenceData),    intent(in) :: convergence_data
      type(AnharmonicData),     intent(in) :: anharmonic_data
      type(BasisStatesPointer)             :: output
    end function
  end interface
  
  interface core_harmonic_observables
    ! ----------------------------------------------------------------------
    ! Operations which calculate thermodynamic data.
    ! ----------------------------------------------------------------------
    ! Calculate the harmonic expectation of the harmonic potential,
    !    but only using the states in the basis.
    ! N.B. this normalises the sum of thermal weights of the included states
    !    to be one,
    !    rather than normalising to the sum across all states to be one.
    module function core_harmonic_observables_WavevectorBasis(bases, &
       & thermal_energy,stress,stress_prefactors,anharmonic_data)    &
       & result(output) 
      type(WavevectorBasis),  intent(in), target   :: bases(:)
      real(dp),               intent(in)           :: thermal_energy
      class(StressBase),      intent(in), optional :: stress
      type(StressPrefactors), intent(in), optional :: stress_prefactors
      type(AnharmonicData),   intent(in)           :: anharmonic_data
      type(ThermodynamicData)                      :: output
    end function
  end interface
  
  interface core_effective_harmonic_observables
    ! Calculate the harmonic expectation of a potential,
    !    but only using the states in the basis.
    ! N.B. this normalises the sum of thermal weights of the included states
    !    to be one,
    !    rather than normalising to the sum across all states to be one.
    module function core_effective_harmonic_observables_WavevectorBasis(bases,thermal_energy,potential,stress,stress_prefactors,anharmonic_data) result(output) 
      type(WavevectorBasis),  intent(in), target   :: bases(:)
      real(dp),               intent(in)           :: thermal_energy
      class(PotentialBase),   intent(in)           :: potential
      class(StressBase),      intent(in), optional :: stress
      type(StressPrefactors), intent(in), optional :: stress_prefactors
      type(AnharmonicData),   intent(in)           :: anharmonic_data
      type(ThermodynamicData)                      :: output
    end function
  end interface
  
  interface core_vci_observables
    ! Calculate thermodynamic data for the core VCI states.
    module function core_vci_observables_WavevectorBasis(bases,             &
       & thermal_energy,states,subspace,subspace_potential,subspace_stress, &
       & stress_prefactors,anharmonic_data) result(output) 
      type(WavevectorBasis),    intent(in)                  :: bases(:)
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
    module function internal_energy_real(basis,harmonic_states, &
       & density_matrix,potential,anharmonic_data) result(output) 
      type(WavevectorBasis),   intent(in)         :: basis
      type(HarmonicStateReal), intent(in), target :: harmonic_states(:)
      type(DensityMatrix),     intent(in)         :: density_matrix
      class(PotentialBase),    intent(in)         :: potential
      type(AnharmonicData),    intent(in)         :: anharmonic_data
      real(dp)                                    :: output
    end function
  end interface
  
  interface
    module function internal_energy_complex(basis,harmonic_states, &
       & density_matrix,potential,anharmonic_data) result(output) 
      type(WavevectorBasis),      intent(in)         :: basis
      type(HarmonicStateComplex), intent(in), target :: harmonic_states(:)
      type(DensityMatrix),        intent(in)         :: density_matrix
      class(PotentialBase),       intent(in)         :: potential
      type(AnharmonicData),       intent(in)         :: anharmonic_data
      real(dp)                                       :: output
    end function
  end interface
  
  interface
    module function stress_real(basis,harmonic_states,density_matrix,stress, &
       & stress_prefactors,anharmonic_data) result(output) 
      type(WavevectorBasis),   intent(in)         :: basis
      type(HarmonicStateReal), intent(in), target :: harmonic_states(:)
      type(DensityMatrix),     intent(in)         :: density_matrix
      class(StressBase),       intent(in)         :: stress
      type(StressPrefactors),  intent(in)         :: stress_prefactors
      type(AnharmonicData),    intent(in)         :: anharmonic_data
      type(RealMatrix)                            :: output
    end function
  end interface
  
  interface
    module function stress_complex(basis,harmonic_states,density_matrix, &
       & stress,stress_prefactors,anharmonic_data) result(output) 
      type(WavevectorBasis),      intent(in)         :: basis
      type(HarmonicStateComplex), intent(in), target :: harmonic_states(:)
      type(DensityMatrix),        intent(in)         :: density_matrix
      class(StressBase),          intent(in)         :: stress
      type(StressPrefactors),     intent(in)         :: stress_prefactors
      type(AnharmonicData),       intent(in)         :: anharmonic_data
      type(RealMatrix)                               :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Mode IDs.
    ! ----------------------------------------------------------------------
    module function mode_ids_WavevectorBasis(this,subspace,anharmonic_data) &
       & result(output) 
      class(WavevectorBasis),   intent(in) :: this
      type(DegenerateSubspace), intent(in) :: subspace
      type(AnharmonicData),     intent(in) :: anharmonic_data
      integer, allocatable                 :: output(:)
    end function
  end interface
  
  interface
    module function paired_mode_ids_WavevectorBasis(this,subspace, &
       & anharmonic_data) result(output) 
      class(WavevectorBasis),   intent(in) :: this
      type(DegenerateSubspace), intent(in) :: subspace
      type(AnharmonicData),     intent(in) :: anharmonic_data
      integer, allocatable                 :: output(:)
    end function
  end interface
  
  interface
    ! Calculate the derivative of the free energy.
    module function free_energy_gradient_WavevectorBasis(this,              &
        subspace_potential,basis_functions,subspace,states,thermal_energy, &
          & state_energy_cutoff,anharmonic_data) result(output) 
      class(WavevectorBasis),   intent(in) :: this
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
  
  interface free_energy_gradient
    module function free_energy_gradient_WavevectorBases(bases,             &
        subspace_potential,basis_functions,subspace,states,thermal_energy, &
          & state_energy_cutoff,anharmonic_data) result(output) 
      type(WavevectorBasis),    intent(in) :: bases(:)
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
    module subroutine read_WavevectorBasis(this,input) 
      class(WavevectorBasis), intent(out) :: this
      type(String),           intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_WavevectorBasis(this) result(output) 
      class(WavevectorBasis), intent(in) :: this
      type(String), allocatable          :: output(:)
    end function
  end interface
  
  interface WavevectorBasis
    module function new_WavevectorBasis_Strings(input) result(this) 
      type(String), intent(in) :: input(:)
      type(WavevectorBasis)    :: this
    end function
  
    impure elemental module function new_WavevectorBasis_StringArray(input) &
       & result(this) 
      type(StringArray), intent(in) :: input
      type(WavevectorBasis)         :: this
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Inherited procedures which do not apply to this type.
    ! ----------------------------------------------------------------------
    impure elemental module function thermodynamic_data_WavevectorBasis(this,thermal_energy,states,subspace,subspace_potential,subspace_stress,stress_prefactors,anharmonic_data) result(output) 
      class(WavevectorBasis),   intent(in)                  :: this
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
    impure elemental module function wavefunctions_WavevectorBasis(this, &
       & states,subspace,anharmonic_data) result(output) 
      class(WavevectorBasis),   intent(in)         :: this
      class(BasisStates),       intent(in), target :: states
      type(DegenerateSubspace), intent(in)         :: subspace
      type(AnharmonicData),     intent(in)         :: anharmonic_data
      type(SubspaceWavefunctionsPointer) :: output
    end function
  end interface
end module
