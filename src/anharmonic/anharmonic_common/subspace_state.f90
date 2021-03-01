! ======================================================================
! The SubspaceState abstract class, defining states which do not need to
!    reference a basis.
! ======================================================================
! N.B. the Basis and State methods take other Basis and State class() arguments
!    with the TARGET attribute, so that they can be cast to their concrete type
!    pointers.
! The %state_pointer() and %states_pointer() methods should
!    not be called on objects which do not have the TARGET attribute,
!    as this will silently lead to undefined behaviour.
module caesar_subspace_state_module
  use caesar_common_module
  
  use caesar_stress_prefactors_module
  use caesar_anharmonic_data_module
  use caesar_sparse_monomial_module
  implicit none
  
  private
  
  public :: SubspaceState
  public :: SubspaceStatePointer
  
  type, abstract, extends(Stringsable) :: SubspaceState
  contains
    procedure(representation_SubspaceState), public, deferred, nopass :: &
       & representation
    procedure, public :: startup => startup_SubspaceState
    
    ! Return the id of the modes across which the state is defined.
    procedure(mode_ids_SubspaceState), public, deferred :: mode_ids
    procedure(paired_mode_ids_SubspaceState), public, deferred :: &
       & paired_mode_ids
    
    ! Returns the occupation of the state.
    procedure(occupation_SubspaceState), public, deferred :: occupation
    
    ! Returns the wavevector of the state.
    procedure(wavevector_SubspaceState), public, deferred :: wavevector
  end type
  
  type, extends(SubspaceState) :: SubspaceStatePointer
    type(String),                      private :: representation_
    ! N.B. state_ should never be modified.
    ! It is public for performance reasons only.
    class(SubspaceState), allocatable, public :: state_
  contains
    procedure, private :: check => check_SubspaceStatePointer
    
    procedure, public, nopass :: representation => &
                               & representation_SubspaceStatePointer
    
    procedure, public :: state         => state_SubspaceStatePointer
    procedure, public :: state_pointer => state_pointer_SubspaceStatePointer
    
    procedure, public :: mode_ids => &
                       & mode_ids_SubspaceStatePointer
    procedure, public :: paired_mode_ids => &
                       & paired_mode_ids_SubspaceStatePointer
    
    procedure, public :: occupation => occupation_SubspaceStatePointer
    
    procedure, public :: wavevector => wavevector_SubspaceStatePointer
    
    ! I/O.
    procedure, public :: read  => read_SubspaceStatePointer
    procedure, public :: write => write_SubspaceStatePointer
  end type
  
  ! Abstract interface for SubspaceState functionality.
  abstract interface
    impure elemental function representation_SubspaceState() result(output)
      import String
      implicit none
      
      type(String) :: output
    end function
    
    function mode_ids_SubspaceState(this) result(output)
      import SubspaceState
      implicit none
      
      class(SubspaceState), intent(in) :: this
      integer, allocatable             :: output(:)
    end function
    
    function paired_mode_ids_SubspaceState(this) result(output)
      import SubspaceState
      implicit none
      
      class(SubspaceState), intent(in) :: this
      integer, allocatable             :: output(:)
    end function
    
    impure elemental function occupation_SubspaceState(this) result(output)
      import SubspaceState
      implicit none
      
      class(SubspaceState), intent(in) :: this
      integer                          :: output
    end function
    
    function wavevector_SubspaceState(this,modes,qpoints) &
       & result(output)
      import SubspaceState
      import ComplexMode
      import QpointData
      import FractionVector
      implicit none
      
      class(SubspaceState), intent(in) :: this
      type(ComplexMode),    intent(in) :: modes(:)
      type(QpointData),     intent(in) :: qpoints(:)
      type(FractionVector)             :: output
    end function
  end interface
  
  interface
    ! Startup method.
    module subroutine startup_SubspaceState(this) 
      class(SubspaceState), intent(in) :: this
    end subroutine
  end interface
  
  interface SubspaceStatePointer
    ! --------------------------------------------------
    ! SubspaceStatePointer methods.
    ! --------------------------------------------------
    ! Construct a SubspaceStatePointer from any type which extends SubspaceState.
    impure elemental module function new_SubspaceStatePointer(state) &
       & result(this) 
      class(SubspaceState), intent(in) :: state
      type(SubspaceStatePointer)       :: this
    end function
  end interface
  
  interface
    ! Checks that the pointer has been allocated before it is used.
    module subroutine check_SubspaceStatePointer(this) 
      class(SubspaceStatePointer), intent(in) :: this
    end subroutine
  end interface
  
  interface
    ! Type representation.
    impure elemental module function representation_SubspaceStatePointer() &
       & result(output) 
      type(String) :: output
    end function
  end interface
  
  interface
    ! SubspaceState methods.
    module function state_SubspaceStatePointer(this) result(output) 
      class(SubspaceStatePointer), intent(in) :: this
      class(SubspaceState), allocatable       :: output
    end function
  end interface
  
  interface
    module function state_pointer_SubspaceStatePointer(this) result(output) 
      class(SubspaceStatePointer), intent(in), target :: this
      class(SubspaceState), pointer                   :: output
    end function
  end interface
  
  interface
    ! --------------------------------------------------
    ! SubspaceStatePointer wrappers for SubspaceState methods.
    ! --------------------------------------------------
    module function mode_ids_SubspaceStatePointer(this) result(output) 
      class(SubspaceStatePointer), intent(in) :: this
      integer, allocatable                    :: output(:)
    end function
  end interface
  
  interface
    module function paired_mode_ids_SubspaceStatePointer(this) result(output) 
      class(SubspaceStatePointer), intent(in) :: this
      integer, allocatable                    :: output(:)
    end function
  end interface
  
  interface
    impure elemental module function occupation_SubspaceStatePointer(this) &
       & result(output) 
      class(SubspaceStatePointer), intent(in) :: this
      integer                                 :: output
    end function
  end interface
  
  interface
    module function wavevector_SubspaceStatePointer(this,modes,qpoints) &
       & result(output) 
      class(SubspaceStatePointer), intent(in) :: this
      type(ComplexMode),           intent(in) :: modes(:)
      type(QpointData),            intent(in) :: qpoints(:)
      type(FractionVector)                    :: output
    end function
  end interface
  
  interface
    ! --------------------------------------------------
    ! SubspaceStatePointer I/O.
    ! --------------------------------------------------
    module subroutine read_SubspaceStatePointer(this,input) 
      class(SubspaceStatePointer), intent(out) :: this
      type(String),                intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_SubspaceStatePointer(this) result(output) 
      class(SubspaceStatePointer), intent(in) :: this
      type(String), allocatable               :: output(:)
    end function
  end interface
  
  interface SubspaceStatePointer
    module function new_SubspaceStatePointer_Strings(input) result(this) 
      type(String), intent(in)   :: input(:)
      type(SubspaceStatePointer) :: this
    end function
  end interface
  
  interface SubspaceStatePointer
    impure elemental module function new_SubspaceStatePointer_StringArray(input) result(this) 
      type(StringArray), intent(in) :: input
      type(SubspaceStatePointer)    :: this
    end function
  end interface
end module
