! ======================================================================
! A state along each complex mode in a degenerate subspace,
!    at a set of q-points for which 2q/=G.
! ======================================================================
module caesar_harmonic_state_complex_module
  use caesar_common_module
  
  use caesar_anharmonic_common_module
  
  use caesar_harmonic_state_2d_module
  implicit none
  
  private
  
  public :: HarmonicStateComplex
  
  public :: harmonic_state_complex_pointer
  
  public :: prod_complex
  
  type, extends(SubspaceState) :: HarmonicStateComplex
    type(HarmonicState2D), allocatable :: modes_(:)
  contains
    procedure, public, nopass :: representation => &
                               & representation_HarmonicStateComplex
    
    procedure, public :: mode_ids => &
                       & mode_ids_HarmonicStateComplex
    procedure, public :: paired_mode_ids => &
                       & paired_mode_ids_HarmonicStateComplex
    
    procedure, public :: occupation => occupation_HarmonicStateComplex
    
    procedure, public :: change_modes => change_modes_HarmonicStateComplex
    
    procedure, public :: wavevector => wavevector_HarmonicStateComplex
    
    procedure, public :: wavefunction => wavefunction_HarmonicStateComplex
    
    ! I/O.
    procedure, public :: read  => read_HarmonicStateComplex
    procedure, public :: write => write_HarmonicStateComplex
  end type
  
  interface
    impure elemental module function prod_complex(lhs,rhs) result(output) 
      type(HarmonicStateComplex), intent(in) :: lhs
      type(HarmonicStateComplex), intent(in) :: rhs
      type(HarmonicStateComplex)             :: output
    end function
  end interface
  
  interface HarmonicStateComplex
    ! ----------------------------------------------------------------------
    ! Constructor.
    ! ----------------------------------------------------------------------
    module function new_HarmonicStateComplex(modes) result(this) 
      type(HarmonicState2D), intent(in) :: modes(:)
      type(HarmonicStateComplex)        :: this
    end function
  
    recursive module function new_HarmonicStateComplex_SubspaceState(input) &
       & result(this) 
      class(SubspaceState), intent(in) :: input
      type(HarmonicStateComplex)       :: this
    end function
  end interface
  
  interface
    ! Cast a class(SubspaceState) to a pointer of type(HarmonicStateComplex).
    ! N.B. this must only be called on inputs with the TARGET attribute.
    recursive module function harmonic_state_complex_pointer(input) &
       & result(this) 
      class(SubspaceState), intent(in), target :: input
      type(HarmonicStateComplex), pointer      :: this
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Type representation.
    ! ----------------------------------------------------------------------
    impure elemental module function representation_HarmonicStateComplex() &
       & result(output) 
      type(String) :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Returns the modes spanned by the state.
    ! ----------------------------------------------------------------------
    module function mode_ids_HarmonicStateComplex(this) result(output) 
      class(HarmonicStateComplex), intent(in) :: this
      integer, allocatable                    :: output(:)
    end function
  end interface
  
  interface
    module function paired_mode_ids_HarmonicStateComplex(this) result(output) 
      class(HarmonicStateComplex), intent(in) :: this
      integer, allocatable                    :: output(:)
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Returns the total occupation of a given state.
    ! ----------------------------------------------------------------------
    ! The total occupation of the state product_{q,i}|n_{q,i}> is equal to
    !    sum_{q,i} n_{q,i}.
    impure elemental module function occupation_HarmonicStateComplex(this) &
       & result(output) 
      class(HarmonicStateComplex), intent(in) :: this
      integer                                 :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Returns the wavevector of a given state.
    ! ----------------------------------------------------------------------
    ! The wavevector of the state product_{q,i} |n_{q,i}> is equal to
    !    sum_{q,i} n_{q,i}q.
    module function wavevector_HarmonicStateComplex(this,modes,qpoints) &
       & result(output) 
      class(HarmonicStateComplex), intent(in) :: this
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
    impure elemental module function wavefunction_HarmonicStateComplex(this, &
       & frequency,supercell) result(output) 
      class(HarmonicStateComplex), intent(in) :: this
      real(dp),                    intent(in) :: frequency
      type(StructureData),         intent(in) :: supercell
      type(String)                            :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Change the modes of the state by the specified group.
    ! ----------------------------------------------------------------------
    impure elemental module function change_modes_HarmonicStateComplex(this, &
       & mode_group) result(output) 
      class(HarmonicStateComplex), intent(in) :: this
      type(Group),                 intent(in) :: mode_group
      type(HarmonicStateComplex)              :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_HarmonicStateComplex(this,input) 
      class(HarmonicStateComplex), intent(out) :: this
      type(String),                intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_HarmonicStateComplex(this) result(output) 
      class(HarmonicStateComplex), intent(in) :: this
      type(String), allocatable               :: output(:)
    end function
  end interface
  
  interface HarmonicStateComplex
    module function new_HarmonicStateComplex_Strings(input) result(this) 
      type(String), intent(in)   :: input(:)
      type(HarmonicStateComplex) :: this
    end function
  
    impure elemental module function new_HarmonicStateComplex_StringArray(input) result(this) 
      type(StringArray), intent(in) :: input
      type(HarmonicStateComplex)    :: this
    end function
  end interface
end module
