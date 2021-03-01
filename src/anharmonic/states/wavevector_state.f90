! ======================================================================
! A state in a wavevector basis.
! See wavevector_basis.f90 for more information.
! ======================================================================
module caesar_wavevector_state_module
  use caesar_common_module
  
  use caesar_anharmonic_common_module
  implicit none
  
  private
  
  public :: startup_wavevector_state
  
  public :: WavevectorState
  
  public :: wavevector_state_pointer
  
  type, extends(BasisState) :: WavevectorState
    type(FractionVector)  :: wavevector
    integer,  allocatable :: state_ids(:)
    real(dp), allocatable :: coefficients(:)
  contains
    procedure, public, nopass :: representation => &
                               & representation_WavevectorState
    ! I/O.
    procedure, public :: read  => read_WavevectorState
    procedure, public :: write => write_WavevectorState
  end type
  
  interface
    ! Startup procedure.
    module module subroutine startup_wavevector_state() 
    end subroutine
  end interface
  
  interface
    ! Type representation.
    impure elemental module module function representation_WavevectorState() &
       & result(output) 
      type(String) :: output
    end function
  end interface
  
  interface WavevectorState
    ! Constructors.
    module module function new_WavevectorState(subspace_id,wavevector, &
       & state_ids,coefficients) result(this) 
      integer,              intent(in) :: subspace_id
      type(FractionVector), intent(in) :: wavevector
      integer,              intent(in) :: state_ids(:)
      real(dp),             intent(in) :: coefficients(:)
      type(WavevectorState)            :: this
    end function
  
    recursive module module function new_WavevectorState_BasisState(input) &
       & result(this) 
      class(BasisState), intent(in) :: input
      type(WavevectorState)         :: this
    end function
  end interface
  
  interface
    ! Cast a class(BasisState) to a pointer of type(WavevectorState).
    ! N.B. this must only be called on inputs with the TARGET attribute.
    recursive module module function wavevector_state_pointer(input) &
       & result(this) 
      class(BasisState), intent(in), target :: input
      type(WavevectorState), pointer        :: this
    end function
  end interface
  
  interface
    ! I/O.
    module module subroutine read_WavevectorState(this,input) 
      class(WavevectorState), intent(out) :: this
      type(String),           intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module module function write_WavevectorState(this) result(output) 
      class(WavevectorState), intent(in) :: this
      type(String), allocatable          :: output(:)
    end function
  end interface
  
  interface WavevectorState
    module module function new_WavevectorState_Strings(input) result(this) 
      type(String), intent(in) :: input(:)
      type(WavevectorState)    :: this
    end function
  
    impure elemental module module function new_WavevectorState_StringArray(input) result(this) 
      type(StringArray), intent(in) :: input
      type(WavevectorState)         :: this
    end function
  end interface
end module
