! ======================================================================
! A state along each complex mode in a degenerate subspace,
!    at a set of q-points for which 2q=G.
! ======================================================================
module caesar_monomial_state_real_module
  use caesar_common_module
  
  use caesar_anharmonic_common_module
  
  use caesar_monomial_state_1d_module
  implicit none
  
  private
  
  public :: MonomialStateReal
  
  public :: monomial_state_real_pointer
  
  public :: generate_monomial_states
  
  public :: finite_overlap
  
  type, extends(SubspaceState) :: MonomialStateReal
    integer                                     :: supercell_size
    real(dp)                                    :: frequency
    real(dp),                           private :: log_2nw_
    type(MonomialState1D), allocatable, private :: modes_(:)
  contains
    procedure, public, nopass :: representation => &
                               & representation_MonomialStateReal
    
    procedure, public :: mode_ids => mode_ids_MonomialStateReal
    procedure, public :: paired_mode_ids => paired_mode_ids_MonomialStateReal
    
    procedure, public :: occupation => occupation_MonomialStateReal
    
    procedure, public :: change_modes => change_modes_MonomialStateReal
    
    procedure, public :: wavevector => wavevector_MonomialStateReal
    
    procedure, public :: wavefunction => wavefunction_MonomialStateReal
    
    procedure, public :: inner_product => &
                       & inner_product_MonomialStateReal
    procedure, public :: integrate => integrate_MonomialStateReal
    procedure, public :: kinetic_energy => &
                       & kinetic_energy_MonomialStateReal
    procedure, public :: harmonic_potential_energy => &
                       & harmonic_potential_energy_MonomialStateReal
    procedure, public :: kinetic_stress => &
                       & kinetic_stress_MonomialStateReal
    
    ! I/O.
    procedure, public :: read  => read_MonomialStateReal
    procedure, public :: write => write_MonomialStateReal
  end type
  
  interface MonomialStateReal
    ! ----------------------------------------------------------------------
    ! Constructor.
    ! ----------------------------------------------------------------------
    module function new_MonomialStateReal(supercell_size,frequency,modes) &
       & result(this) 
      integer,               intent(in) :: supercell_size
      real(dp),              intent(in) :: frequency
      type(MonomialState1D), intent(in) :: modes(:)
      type(MonomialStateReal)           :: this
    end function
  
    recursive module function new_MonomialStateReal_SubspaceState(input) &
       & result(this) 
      class(SubspaceState), intent(in) :: input
      type(MonomialStateReal)          :: this
    end function
  end interface
  
  interface
    ! Cast a class(SubspaceState) to a pointer of type(MonomialStateReal).
    ! N.B. this must only be called on inputs with the TARGET attribute.
    recursive module function monomial_state_real_pointer(input) result(this) 
      class(SubspaceState), intent(in), target :: input
      type(MonomialStateReal), pointer         :: this
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Type representation.
    ! ----------------------------------------------------------------------
    impure elemental module function representation_MonomialStateReal() &
       & result(output) 
      type(String) :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Returns the modes spanned by the state.
    ! ----------------------------------------------------------------------
    module function mode_ids_MonomialStateReal(this) result(output) 
      class(MonomialStateReal), intent(in) :: this
      integer, allocatable                 :: output(:)
    end function
  end interface
  
  interface
    module function paired_mode_ids_MonomialStateReal(this) result(output) 
      class(MonomialStateReal), intent(in) :: this
      integer, allocatable                 :: output(:)
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Generates all monomial states in a subspace up to a given power.
    ! ----------------------------------------------------------------------
    module function generate_monomial_states(subspace,supercell_size, &
       & frequency,modes,maximum_power) result(output) 
      type(DegenerateSubspace), intent(in) :: subspace
      integer,                  intent(in) :: supercell_size
      real(dp),                 intent(in) :: frequency
      type(ComplexMode),        intent(in) :: modes(:)
      integer,                  intent(in) :: maximum_power
      type(MonomialStateReal), allocatable :: output(:)
    end function
  end interface
  
  interface
    recursive module function generate_monomial_states_helper(ids,power, &
       & state) result(output) 
      integer,                 intent(in)  :: ids(:)
      integer,                 intent(in)  :: power
      type(MonomialStateReal), intent(in)  :: state
      type(MonomialStateReal), allocatable :: output(:)
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Returns the total power of a given state.
    ! ----------------------------------------------------------------------
    ! The total power of the state product_{q,i} |(u_{q,i})^(n_{q,i})> is equal to
    !    sum_{q,i} n_{q,i}.
    impure elemental module function occupation_MonomialStateReal(this) &
       & result(output) 
      class(MonomialStateReal), intent(in) :: this
      integer                              :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Returns the wavevector of a given state.
    ! ----------------------------------------------------------------------
    ! The wavevector of the state product_{q,i} |(u_{q,i})^(n_{q,i})> is equal to
    !    sum_{q,i} n_{q,i}q.
    module function wavevector_MonomialStateReal(this,modes,qpoints) &
       & result(output) 
      class(MonomialStateReal), intent(in) :: this
      type(ComplexMode),        intent(in) :: modes(:)
      type(QpointData),         intent(in) :: qpoints(:)
      type(FractionVector)                 :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Returns the wavefunction of the state,
    !    with all coefficients accounted for.
    ! ----------------------------------------------------------------------
    impure elemental module function wavefunction_MonomialStateReal(this, &
       & frequency,supercell) result(output) 
      class(MonomialStateReal), intent(in) :: this
      real(dp),                 intent(in) :: frequency
      type(StructureData),      intent(in) :: supercell
      type(String)                         :: output
    end function
  end interface
  
  interface finite_overlap
    ! ----------------------------------------------------------------------
    ! SubspaceState methods.
    ! ----------------------------------------------------------------------
    ! Returns whether or not braket(bra,ket) is non-zero.
    impure elemental module function finite_overlap_MonomialStateReals(bra, &
       & ket,anharmonic_data) result(output) 
      type(MonomialStateReal), intent(in) :: bra
      type(MonomialStateReal), intent(in) :: ket
      type(AnharmonicData),    intent(in) :: anharmonic_data
      logical                             :: output
    end function
  end interface
  
  interface
    impure elemental module function inner_product_MonomialStateReal(this, &
       & ket,anharmonic_data) result(output) 
      class(MonomialStateReal), intent(in)                   :: this
      class(SubspaceState),     intent(in), optional, target :: ket
      type(AnharmonicData),     intent(in)                   :: anharmonic_data
      real(dp)                                               :: output
    end function
  end interface
  
  interface
    impure elemental module function integrate_MonomialStateReal(this, &
       & monomial,ket,anharmonic_data) result(output) 
      class(MonomialStateReal), intent(in)                   :: this
      type(SparseMonomial),     intent(in)                   :: monomial
      class(SubspaceState),     intent(in), optional, target :: ket
      type(AnharmonicData),     intent(in)                   :: anharmonic_data
      complex(dp)                                            :: output
    end function
  end interface
  
  interface
    impure elemental module function kinetic_energy_MonomialStateReal(this, &
       & ket,anharmonic_data) result(output) 
      class(MonomialStateReal), intent(in)                   :: this
      class(SubspaceState),     intent(in), optional, target :: ket
      type(AnharmonicData),     intent(in)                   :: anharmonic_data
      real(dp)                                               :: output
    end function
  end interface
  
  interface
    impure elemental module function harmonic_potential_energy_MonomialStateReal(   this,ket,anharmonic_data) result(output) 
      class(MonomialStateReal), intent(in)                   :: this
      class(SubspaceState),     intent(in), optional, target :: ket
      type(AnharmonicData),     intent(in)                   :: anharmonic_data
      real(dp)                                               :: output
    end function
  end interface
  
  interface
    impure elemental module function kinetic_stress_MonomialStateReal(this, &
       & ket,stress_prefactors,anharmonic_data) result(output) 
      class(MonomialStateReal), intent(in)                   :: this
      class(SubspaceState),     intent(in), optional, target :: ket
      type(StressPrefactors),   intent(in)                   :: stress_prefactors
      type(AnharmonicData),     intent(in)                   :: anharmonic_data
      type(RealMatrix)                                       :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Change the modes of the state by the specified group.
    ! ----------------------------------------------------------------------
    impure elemental module function change_modes_MonomialStateReal(this, &
       & mode_group) result(output) 
      class(MonomialStateReal), intent(in) :: this
      type(Group),          intent(in) :: mode_group
      type(MonomialStateReal)              :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_MonomialStateReal(this,input) 
      class(MonomialStateReal), intent(out) :: this
      type(String),             intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_MonomialStateReal(this) result(output) 
      class(MonomialStateReal), intent(in) :: this
      type(String), allocatable        :: output(:)
    end function
  end interface
  
  interface MonomialStateReal
    module function new_MonomialStateReal_Strings(input) result(this) 
      type(String), intent(in) :: input(:)
      type(MonomialStateReal)      :: this
    end function
  
    impure elemental module function new_MonomialStateReal_StringArray(input) &
       & result(this) 
      type(StringArray), intent(in) :: input
      type(MonomialStateReal)           :: this
    end function
  end interface
end module
