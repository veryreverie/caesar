! ======================================================================
! A set of states in a wavevector basis.
! See wavevector_basis.f90 for more details.
! ======================================================================
module caesar_wavevector_states_module
  use caesar_common_module
  
  use caesar_anharmonic_common_module
  
  use caesar_wavevector_state_module
  use caesar_density_matrix_module
  implicit none
  
  private
  
  public :: WavevectorStates
  
  public :: wavevector_states_pointer
  
  type, extends(BasisStates) :: WavevectorStates
    type(WavevectorState), allocatable :: states(:)
    real(dp),              allocatable :: energies(:)
    real(dp),              allocatable :: weights(:)
    type(DensityMatrix),   allocatable :: density_matrices(:)
  contains
    procedure, public, nopass :: representation => &
                               & representation_WavevectorStates
    
    ! I/O.
    procedure, public :: read  => read_WavevectorStates
    procedure, public :: write => write_WavevectorStates
  end type
  
  interface WavevectorStates
    ! ----------------------------------------------------------------------
    ! WavevectorStates methods.
    ! ----------------------------------------------------------------------
    ! Constructors.
    module function new_WavevectorStates(subspace_id,states,energies,weights) &
       & result(this) 
      integer,               intent(in)           :: subspace_id
      type(WavevectorState), intent(in)           :: states(:)
      real(dp),              intent(in)           :: energies(:)
      real(dp),              intent(in), optional :: weights(:)
      type(WavevectorStates)                      :: this
    end function
  
    recursive module function new_WavevectorStates_BasisStates(input) &
       & result(this) 
      class(BasisStates), intent(in) :: input
      type(WavevectorStates)         :: this
    end function
  end interface
  
  interface
    ! Cast a class(BasisStates) to a pointer of type(WavevectorStates).
    ! N.B. this must only be called on inputs with the TARGET attribute.
    recursive module function wavevector_states_pointer(input) result(this) 
      class(BasisStates), intent(in), target :: input
      type(WavevectorStates), pointer        :: this
    end function
  end interface
  
  interface
    ! Type representation.
    impure elemental module function representation_WavevectorStates() &
       & result(output) 
      type(String) :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_WavevectorStates(this,input) 
      class(WavevectorStates), intent(out) :: this
      type(String),            intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_WavevectorStates(this) result(output) 
      class(WavevectorStates), intent(in) :: this
      type(String), allocatable           :: output(:)
    end function
  end interface
  
  interface WavevectorStates
    module function new_WavevectorStates_Strings(input) result(this) 
      type(String), intent(in) :: input(:)
      type(WavevectorStates)   :: this
    end function
  
    impure elemental module function new_WavevectorStates_StringArray(input) &
       & result(this) 
      type(StringArray), intent(in) :: input
      type(WavevectorStates)        :: this
    end function
  end interface
end module
