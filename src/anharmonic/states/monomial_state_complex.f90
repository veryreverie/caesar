! ======================================================================
! A state along each complex mode in a degenerate subspace,
!    at a set of q-points for which 2q/=G.
! ======================================================================
module caesar_monomial_state_complex_module
  use caesar_common_module
  
  use caesar_anharmonic_common_module
  
  use caesar_monomial_state_2d_module
  implicit none
  
  private
  
  public :: startup_monomial_state_complex
  
  public :: MonomialStateComplex
  
  public :: monomial_state_complex_pointer
  
  public :: finite_overlap
  
  type, extends(SubspaceState) :: MonomialStateComplex
    integer                                     :: supercell_size
    real(dp)                                    :: frequency
    real(dp),                           private :: log_2nw_
    type(MonomialState2D), allocatable, private :: modes_(:)
  contains
    procedure, public, nopass :: representation => &
                               & representation_MonomialStateComplex
    
    procedure, public :: mode_ids => mode_ids_MonomialStateComplex
    procedure, public :: paired_mode_ids => &
                       & paired_mode_ids_MonomialStateComplex
    
    procedure, public :: occupation => occupation_MonomialStateComplex
    
    procedure, public :: change_modes => change_modes_MonomialStateComplex
    
    procedure, public :: wavevector => wavevector_MonomialStateComplex
    
    procedure, public :: wavefunction => wavefunction_MonomialStateComplex
    
    procedure, public :: inner_product => &
                       & inner_product_MonomialStateComplex
    procedure, public :: integrate => integrate_MonomialStateComplex
    procedure, public :: kinetic_energy => &
                       & kinetic_energy_MonomialStateComplex
    procedure, public :: harmonic_potential_energy => &
                       & harmonic_potential_energy_MonomialStateComplex
    procedure, public :: kinetic_stress => &
                       & kinetic_stress_MonomialStateComplex
    
    ! I/O.
    procedure, public :: read  => read_MonomialStateComplex
    procedure, public :: write => write_MonomialStateComplex
  end type
  
  interface
    ! Startup procedure.
    module subroutine startup_monomial_state_complex() 
    end subroutine
  end interface
  
  interface MonomialStateComplex
    ! ----------------------------------------------------------------------
    ! Constructor.
    ! ----------------------------------------------------------------------
    module function new_MonomialStateComplex(supercell_size,frequency,modes) &
       & result(this) 
      integer,               intent(in) :: supercell_size
      real(dp),              intent(in) :: frequency
      type(MonomialState2D), intent(in) :: modes(:)
      type(MonomialStateComplex)        :: this
    end function
  
    recursive module function new_MonomialStateComplex_SubspaceState(input) &
       & result(this) 
      class(SubspaceState), intent(in) :: input
      type(MonomialStateComplex)       :: this
    end function
  end interface
  
  interface
    ! Cast a class(SubspaceState) to a pointer of type(MonomialStateComplex).
    ! N.B. this must only be called on inputs with the TARGET attribute.
    recursive module function monomial_state_complex_pointer(input) &
       & result(this) 
      class(SubspaceState), intent(in), target :: input
      type(MonomialStateComplex), pointer      :: this
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Type representation.
    ! ----------------------------------------------------------------------
    impure elemental module function representation_MonomialStateComplex() &
       & result(output) 
      type(String) :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Returns the modes spanned by the state.
    ! ----------------------------------------------------------------------
    module function mode_ids_MonomialStateComplex(this) result(output) 
      class(MonomialStateComplex), intent(in) :: this
      integer, allocatable                    :: output(:)
    end function
  end interface
  
  interface
    module function paired_mode_ids_MonomialStateComplex(this) result(output) 
      class(MonomialStateComplex), intent(in) :: this
      integer, allocatable                    :: output(:)
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Returns the total power of a given state.
    ! ----------------------------------------------------------------------
    ! The total power of the state product_{q,i} |(u_{q,i})^(n_{q,i})> is equal to
    !    sum_{q,i} n_{q,i}.
    impure elemental module function occupation_MonomialStateComplex(this) &
       & result(output) 
      class(MonomialStateComplex), intent(in) :: this
      integer                                 :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Returns the wavevector of a given state.
    ! ----------------------------------------------------------------------
    ! The wavevector of the state product_{q,i} |(u_{q,i})^(n_{q,i})> is equal to
    !    sum_{q,i} n_{q,i}q.
    module function wavevector_MonomialStateComplex(this,modes,qpoints) &
       & result(output) 
      class(MonomialStateComplex), intent(in) :: this
      type(ComplexMode),           intent(in) :: modes(:)
      type(QpointData),            intent(in) :: qpoints(:)
      type(FractionVector)                    :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Returns the wavefunction of the state,
    !    with all coefficients accounted for.
    ! ----------------------------------------------------------------------
    impure elemental module function wavefunction_MonomialStateComplex(this, &
       & frequency,supercell) result(output) 
      class(MonomialStateComplex), intent(in) :: this
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
    impure elemental module function finite_overlap_MonomialStateComplexs(bra,ket,anharmonic_data) result(output) 
      type(MonomialStateComplex), intent(in) :: bra
      type(MonomialStateComplex), intent(in) :: ket
      type(AnharmonicData),       intent(in) :: anharmonic_data
      logical                                :: output
    end function
  end interface
  
  interface
    impure elemental module function inner_product_MonomialStateComplex(this,ket,anharmonic_data) result(output) 
      class(MonomialStateComplex), intent(in)                   :: this
      class(SubspaceState),        intent(in), optional, target :: ket
      type(AnharmonicData),        intent(in)                   :: anharmonic_data
      real(dp)                                                  :: output
    end function
  end interface
  
  interface
    impure elemental module function integrate_MonomialStateComplex(this, &
       & monomial,ket,anharmonic_data) result(output) 
      class(MonomialStateComplex), intent(in)                   :: this
      type(SparseMonomial),        intent(in)                   :: monomial
      class(SubspaceState),        intent(in), optional, target :: ket
      type(AnharmonicData),        intent(in)                   :: anharmonic_data
      complex(dp)                                               :: output
    end function
  end interface
  
  interface
    impure elemental module function kinetic_energy_MonomialStateComplex(this,ket,anharmonic_data) result(output) 
      class(MonomialStateComplex), intent(in)                   :: this
      class(SubspaceState),        intent(in), optional, target :: ket
      type(AnharmonicData),        intent(in)                   :: anharmonic_data
      real(dp)                                                  :: output
    end function
  end interface
  
  interface
    impure elemental module function harmonic_potential_energy_MonomialStateComplex(   this,ket,anharmonic_data) result(output) 
      class(MonomialStateComplex), intent(in)                   :: this
      class(SubspaceState),        intent(in), optional, target :: ket
      type(AnharmonicData),        intent(in)                   :: anharmonic_data
      real(dp)                                                  :: output
    end function
  end interface
  
  interface
    impure elemental module function kinetic_stress_MonomialStateComplex(this,ket,stress_prefactors,anharmonic_data) result(output) 
      class(MonomialStateComplex), intent(in)                   :: this
      class(SubspaceState),        intent(in), optional, target :: ket
      type(StressPrefactors),      intent(in)                   :: stress_prefactors
      type(AnharmonicData),        intent(in)                   :: anharmonic_data
      type(RealMatrix)                                          :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Change the modes of the state by the specified group.
    ! ----------------------------------------------------------------------
    impure elemental module function change_modes_MonomialStateComplex(this, &
       & mode_group) result(output) 
      class(MonomialStateComplex), intent(in) :: this
      type(Group),                 intent(in) :: mode_group
      type(MonomialStateComplex)              :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_MonomialStateComplex(this,input) 
      class(MonomialStateComplex), intent(out) :: this
      type(String),                intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_MonomialStateComplex(this) result(output) 
      class(MonomialStateComplex), intent(in) :: this
      type(String), allocatable        :: output(:)
    end function
  end interface
  
  interface MonomialStateComplex
    module function new_MonomialStateComplex_Strings(input) result(this) 
      type(String), intent(in)   :: input(:)
      type(MonomialStateComplex) :: this
    end function
  
    impure elemental module function new_MonomialStateComplex_StringArray(input) result(this) 
      type(StringArray), intent(in) :: input
      type(MonomialStateComplex)           :: this
    end function
  end interface
end module
