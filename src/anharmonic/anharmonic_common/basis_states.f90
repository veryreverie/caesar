! ======================================================================
! A set of basis states, defined in terms of a SubspaceBasis.
! ======================================================================
module caesar_basis_states_module
  use caesar_common_module
  use caesar_expectation_cache_module
  implicit none
  
  private
  
  public :: BasisStates
  public :: BasisStatesPointer
  
  type, abstract, extends(Stringsable) :: BasisStates
    integer :: subspace_id
    
    type(ExpectationCache) :: expectation_cache
  contains
    procedure(representation_BasisStates), public, deferred, nopass :: &
       & representation
  end type
  
  type, extends(BasisStates) :: BasisStatesPointer
    type(String),                    private :: representation_
    class(BasisStates), allocatable, private :: states_
  contains
    procedure, private :: check => check_BasisStatesPointer
    
    procedure, public, nopass :: representation => &
                               & representation_BasisStatesPointer
    
    procedure, public :: states => states_BasisStatesPointer
    procedure, public :: states_pointer => states_pointer_BasisStatesPointer
    
    ! I/O.
    procedure, public :: read  => read_BasisStatesPointer
    procedure, public :: write => write_BasisStatesPointer
  end type
  
  abstract interface
    impure elemental function representation_BasisStates() result(output)
      import String
      implicit none
      
      type(String) :: output
    end function
  end interface
  
  interface
    module function types_BasisStates() result(output)
      type(BasisStatesPointer), allocatable :: output(:)
    end function
  end interface
  
  interface BasisStatesPointer
    ! ----------------------------------------------------------------------
    ! BasisStatesPointer methods.
    ! ----------------------------------------------------------------------
    ! Construct a BasisStatesPointer from any type which extends BasisStates.
    impure elemental module function new_BasisStatesPointer(states) &
       & result(this) 
      class(BasisStates), intent(in) :: states
      type(BasisStatesPointer)       :: this
    end function
  end interface
  
  interface
    ! Checks that the pointer has been allocated before it is used.
    module subroutine check_BasisStatesPointer(this) 
      class(BasisStatesPointer), intent(in) :: this
    end subroutine
  end interface
  
  interface
    ! Type representation.
    impure elemental module function representation_BasisStatesPointer() &
       & result(output) 
      type(String) :: output
    end function
  end interface
  
  interface
    ! BasisStates methods.
    module function states_BasisStatesPointer(this) result(output) 
      class(BasisStatesPointer), intent(in) :: this
      class(BasisStates), allocatable       :: output
    end function
  end interface
  
  interface
    module function states_pointer_BasisStatesPointer(this) result(output) 
      class(BasisStatesPointer), intent(in), target :: this
      class(BasisStates), pointer                   :: output
    end function
  end interface
  
  interface
    ! I/O.
    module subroutine read_BasisStatesPointer(this,input) 
      class(BasisStatesPointer), intent(out) :: this
      type(String),              intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_BasisStatesPointer(this) result(output) 
      class(BasisStatesPointer), intent(in) :: this
      type(String), allocatable             :: output(:)
    end function
  end interface
  
  interface BasisStatesPointer
    module function new_BasisStatesPointer_Strings(input) result(this) 
      type(String), intent(in) :: input(:)
      type(BasisStatesPointer) :: this
    end function
  
    impure elemental module function new_BasisStatesPointer_StringArray(input) result(this) 
      type(StringArray), intent(in) :: input
      type(BasisStatesPointer)      :: this
    end function
  end interface
end module
