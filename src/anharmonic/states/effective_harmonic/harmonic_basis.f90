! ======================================================================
! A basis of harmonic states.
! ======================================================================
module caesar_harmonic_basis_module
  use caesar_common_module
  
  use caesar_anharmonic_common_module
  
  use caesar_harmonic_states_module
  implicit none
  
  private
  
  public :: startup_harmonic_basis
  
  public :: HarmonicBasis
  
  type, extends(SubspaceBasis) :: HarmonicBasis
    integer :: subspace_id
    integer :: supercell_size
  contains
    procedure, public, nopass :: representation => representation_HarmonicBasis
    
    procedure, public :: initial_states => initial_states_HarmonicBasis
    procedure, public :: calculate_states => calculate_states_HarmonicBasis
    procedure, public :: mode_ids => mode_ids_HarmonicBasis
    procedure, public :: paired_mode_ids => paired_mode_ids_HarmonicBasis
    
    ! Procedures involving individual states.
    ! N.B. these are all left blank, as individual harmonic states are
    !    currently treated under a different framework.
    procedure, public :: inner_product => &
                       & inner_product_HarmonicBasis
    procedure, public :: integrate_BasisState => &
                       & integrate_BasisState_HarmonicBasis
    procedure, public :: kinetic_energy => &
                       & kinetic_energy_HarmonicBasis
    procedure, public :: harmonic_potential_energy => &
                       & harmonic_potential_energy_HarmonicBasis
    procedure, public :: kinetic_stress => &
                       & kinetic_stress_HarmonicBasis
    
    ! Procedures involving sets of states.
    procedure, public :: thermodynamic_data => &
                       & thermodynamic_data_HarmonicBasis
    procedure, public :: wavefunctions => &
                       & wavefunctions_HarmonicBasis
    procedure, public :: integrate_BasisStates => &
                       & integrate_BasisStates_HarmonicBasis
    
    ! Calculating the derivative of the free energy with respect to
    !    basis function coefficients.
    procedure, public :: free_energy_gradient => &
                       & free_energy_gradient_HarmonicBasis
    
    ! I/O.
    procedure, public :: read  => read_HarmonicBasis
    procedure, public :: write => write_HarmonicBasis
  end type
  
  interface
    ! Startup procedure and type representation.
    module subroutine startup_harmonic_basis() 
    end subroutine
  end interface
  
  interface
    impure elemental module function representation_HarmonicBasis() &
       & result(output) 
      type(String) :: output
    end function
  end interface
  
  interface HarmonicBasis
    ! Constructor.
    impure elemental module function new_HarmonicBasis(subspace_id, &
       & frequency,supercell_size) result(this) 
      integer,  intent(in) :: subspace_id
      real(dp), intent(in) :: frequency
      integer,  intent(in) :: supercell_size
      type(HarmonicBasis)  :: this
    end function
  end interface
  
  interface
    ! Calculate states.
    impure elemental module function initial_states_HarmonicBasis(this, &
       & subspace,thermal_energy,anharmonic_data) result(output) 
      class(HarmonicBasis),     intent(in) :: this
      type(DegenerateSubspace), intent(in) :: subspace
      real(dp),                 intent(in) :: thermal_energy
      type(AnharmonicData),     intent(in) :: anharmonic_data
      type(BasisStatesPointer)             :: output
    end function
  end interface
  
  interface
    impure elemental module function calculate_states_HarmonicBasis(this, &
       & subspace,subspace_potential,thermal_energy,state_energy_cutoff,  &
       & convergence_data,anharmonic_data) result(output) 
      class(HarmonicBasis),     intent(in) :: this
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
    module function mode_ids_HarmonicBasis(this,subspace,anharmonic_data) &
       & result(output) 
      class(HarmonicBasis),     intent(in) :: this
      type(DegenerateSubspace), intent(in) :: subspace
      type(AnharmonicData),     intent(in) :: anharmonic_data
      integer, allocatable                 :: output(:)
    end function
  end interface
  
  interface
    module function paired_mode_ids_HarmonicBasis(this,subspace, &
       & anharmonic_data) result(output) 
      class(HarmonicBasis),     intent(in) :: this
      type(DegenerateSubspace), intent(in) :: subspace
      type(AnharmonicData),     intent(in) :: anharmonic_data
      integer, allocatable                 :: output(:)
    end function
  end interface
  
  interface
    ! Procedures involving individual states.
    ! N.B. these are all left blank, as individual harmonic states are
    !    currently treated under a different framework.
    impure elemental module function inner_product_HarmonicBasis(this,bra, &
       & ket,subspace,anharmonic_data) result(output) 
      class(HarmonicBasis),     intent(in)                   :: this
      class(BasisState),        intent(in),           target :: bra
      class(BasisState),        intent(in), optional, target :: ket
      type(DegenerateSubspace), intent(in)                   :: subspace
      type(AnharmonicData),     intent(in)                   :: anharmonic_data
      real(dp)                                               :: output
    end function
  end interface
  
  interface
    impure elemental module function integrate_BasisState_HarmonicBasis(this,bra,monomial,ket,subspace,anharmonic_data) result(output) 
      class(HarmonicBasis),     intent(in)                   :: this
      class(BasisState),        intent(in),           target :: bra
      type(SparseMonomial),     intent(in)                   :: monomial
      class(BasisState),        intent(in), optional, target :: ket
      type(DegenerateSubspace), intent(in)                   :: subspace
      type(AnharmonicData),     intent(in)                   :: anharmonic_data
      complex(dp)                                            :: output
    end function
  end interface
  
  interface
    impure elemental module function kinetic_energy_HarmonicBasis(this,bra, &
       & ket,subspace,anharmonic_data) result(output) 
      class(HarmonicBasis),     intent(in)                   :: this
      class(BasisState),        intent(in),           target :: bra
      class(BasisState),        intent(in), optional, target :: ket
      type(DegenerateSubspace), intent(in)                   :: subspace
      type(AnharmonicData),     intent(in)                   :: anharmonic_data
      real(dp)                                               :: output
    end function
  end interface
  
  interface
    impure elemental module function harmonic_potential_energy_HarmonicBasis(   this,bra,ket,subspace,anharmonic_data) result(output) 
      class(HarmonicBasis),     intent(in)                   :: this
      class(BasisState),        intent(in),           target :: bra
      class(BasisState),        intent(in), optional, target :: ket
      type(DegenerateSubspace), intent(in)                   :: subspace
      type(AnharmonicData),     intent(in)                   :: anharmonic_data
      real(dp)                                               :: output
    end function
  end interface
  
  interface
    impure elemental module function kinetic_stress_HarmonicBasis(this,bra, &
       & ket,subspace,stress_prefactors,anharmonic_data) result(output) 
      class(HarmonicBasis),     intent(in)                   :: this
      class(BasisState),        intent(in),           target :: bra
      class(BasisState),        intent(in), optional, target :: ket
      type(DegenerateSubspace), intent(in)                   :: subspace
      type(StressPrefactors),   intent(in)                   :: stress_prefactors
      type(AnharmonicData),     intent(in)                   :: anharmonic_data
      type(RealMatrix)                                       :: output
    end function
  end interface
  
  interface
    ! Procedures involving sets of states.
    impure elemental module function thermodynamic_data_HarmonicBasis(this, &
       & thermal_energy,states,subspace,subspace_potential,subspace_stress, &
       & stress_prefactors,anharmonic_data) result(output) 
      class(HarmonicBasis),     intent(in)                  :: this
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
    impure elemental module function wavefunctions_HarmonicBasis(this, &
       & states,subspace,anharmonic_data) result(output) 
      class(HarmonicBasis),      intent(in)         :: this
      class(BasisStates),        intent(in), target :: states
      type(DegenerateSubspace),  intent(in)         :: subspace
      type(AnharmonicData),      intent(in)         :: anharmonic_data
      type(SubspaceWavefunctionsPointer) :: output
    end function
  end interface
  
  interface
    impure elemental module function integrate_BasisStates_HarmonicBasis(this,states,monomial,subspace,anharmonic_data) result(output) 
      class(HarmonicBasis),     intent(in)         :: this
      class(BasisStates),       intent(in), target :: states
      type(SparseMonomial),     intent(in)         :: monomial
      type(DegenerateSubspace), intent(in)         :: subspace
      type(AnharmonicData),     intent(in)         :: anharmonic_data
      complex(dp)                                  :: output
    end function
  end interface
  
  interface
    ! Calculate the derivative of the free energy.
    module function free_energy_gradient_HarmonicBasis(this,                &
        subspace_potential,basis_functions,subspace,states,thermal_energy, &
          & state_energy_cutoff,anharmonic_data) result(output) 
      class(HarmonicBasis),     intent(in) :: this
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
    ! I/O.
    module subroutine read_HarmonicBasis(this,input) 
      class(HarmonicBasis), intent(out) :: this
      type(String),         intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_HarmonicBasis(this) result(output) 
      class(HarmonicBasis), intent(in) :: this
      type(String), allocatable        :: output(:)
    end function
  end interface
  
  interface HarmonicBasis
    module function new_HarmonicBasis_Strings(input) result(this) 
      type(String), intent(in) :: input(:)
      type(HarmonicBasis)      :: this
    end function
  
    impure elemental module function new_HarmonicBasis_StringArray(input) &
       & result(this) 
      type(StringArray), intent(in) :: input
      type(HarmonicBasis)           :: this
    end function
  end interface
end module
