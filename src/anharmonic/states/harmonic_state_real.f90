! ======================================================================
! A state along each complex mode in a degenerate subspace,
!    at a set of q-points for which 2q=G.
! ======================================================================
module caesar_harmonic_state_real_module
  use caesar_common_module
  
  use caesar_anharmonic_common_module
  
  use caesar_harmonic_state_1d_module
  implicit none
  
  private
  
  public :: HarmonicStateReal
  
  public :: harmonic_state_real_pointer
  
  public :: prod_real
  
  type, extends(SubspaceState) :: HarmonicStateReal
    type(HarmonicState1D), allocatable :: modes_(:)
  contains
    procedure, public, nopass :: representation => &
                               & representation_HarmonicStateReal
    
    procedure, public :: mode_ids => mode_ids_HarmonicStateReal
    procedure, public :: paired_mode_ids => paired_mode_ids_HarmonicStateReal
    
    procedure, public :: change_modes => change_modes_HarmonicStateReal
    
    procedure, public :: occupation => occupation_HarmonicStateReal
    
    procedure, public :: wavevector => wavevector_HarmonicStateReal
    
    procedure, public :: wavefunction => wavefunction_HarmonicStateReal
    
    ! I/O.
    procedure, public :: read  => read_HarmonicStateReal
    procedure, public :: write => write_HarmonicStateReal
  end type
  
  interface
    impure elemental module function prod_real(lhs,rhs) result(output) 
      type(HarmonicStateReal), intent(in) :: lhs
      type(HarmonicStateReal), intent(in) :: rhs
      type(HarmonicStateReal)             :: output
    end function
  end interface
  
  interface HarmonicStateReal
    ! ----------------------------------------------------------------------
    ! Constructor.
    ! ----------------------------------------------------------------------
    module function new_HarmonicStateReal(modes) result(this) 
      type(HarmonicState1D), intent(in) :: modes(:)
      type(HarmonicStateReal)           :: this
    end function
  
    recursive module function new_HarmonicStateReal_SubspaceState(input) &
       & result(this) 
      class(SubspaceState), intent(in) :: input
      type(HarmonicStateReal)          :: this
    end function
  end interface
  
  interface
    ! Cast a class(SubspaceState) to a pointer of type(HarmonicStateReal).
    ! N.B. this must only be called on inputs with the TARGET attribute.
    recursive module function harmonic_state_real_pointer(input) result(this) 
      class(SubspaceState), intent(in), target :: input
      type(HarmonicStateReal), pointer         :: this
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Type representation.
    ! ----------------------------------------------------------------------
    impure elemental module function representation_HarmonicStateReal() &
       & result(output) 
      type(String) :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Returns the modes spanned by the state.
    ! ----------------------------------------------------------------------
    module function mode_ids_HarmonicStateReal(this) result(output) 
      class(HarmonicStateReal), intent(in) :: this
      integer, allocatable                 :: output(:)
    end function
  end interface
  
  interface
    module function paired_mode_ids_HarmonicStateReal(this) result(output) 
      class(HarmonicStateReal), intent(in) :: this
      integer, allocatable                 :: output(:)
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Returns the total occupation of a given state.
    ! ----------------------------------------------------------------------
    ! The total occupation of the state product_{q,i} |n_{q,i}> is equal to
    !    sum_{q,i} n_{q,i}.
    impure elemental module function occupation_HarmonicStateReal(this) &
       & result(output) 
      class(HarmonicStateReal), intent(in) :: this
      integer                              :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Returns the wavevector of a given state.
    ! ----------------------------------------------------------------------
    ! The wavevector of the state product_{q,i} |n_{q,i}> is equal to
    !    sum_{q,i} n_{q,i}q.
    module function wavevector_HarmonicStateReal(this,modes,qpoints) &
       & result(output) 
      class(HarmonicStateReal), intent(in) :: this
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
    impure elemental module function wavefunction_HarmonicStateReal(this, &
       & frequency,supercell) result(output) 
      class(HarmonicStateReal), intent(in) :: this
      real(dp),                 intent(in) :: frequency
      type(StructureData),      intent(in) :: supercell
      type(String)                         :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Change the modes of the state by the specified group.
    ! ----------------------------------------------------------------------
    impure elemental module function change_modes_HarmonicStateReal(this, &
       & mode_group) result(output) 
      class(HarmonicStateReal), intent(in) :: this
      type(Group),              intent(in) :: mode_group
      type(HarmonicStateReal)              :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_HarmonicStateReal(this,input) 
      class(HarmonicStateReal), intent(out) :: this
      type(String),             intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_HarmonicStateReal(this) result(output) 
      class(HarmonicStateReal), intent(in) :: this
      type(String), allocatable        :: output(:)
    end function
  end interface
  
  interface HarmonicStateReal
    module function new_HarmonicStateReal_Strings(input) result(this) 
      type(String), intent(in) :: input(:)
      type(HarmonicStateReal)  :: this
    end function
  
    impure elemental module function new_HarmonicStateReal_StringArray(input) &
       & result(this) 
      type(StringArray), intent(in) :: input
      type(HarmonicStateReal)       :: this
    end function
  end interface
end module
