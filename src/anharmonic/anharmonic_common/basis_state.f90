! ======================================================================
! A basis state, defined in terms of a SubspaceBasis.
! ======================================================================
module caesar_basis_state_module
  use caesar_common_module
  implicit none
  
  private
  
  public :: BasisState
  public :: BasisStatePointer
  
  type, abstract, extends(Stringsable) :: BasisState
    integer :: subspace_id
  contains
    procedure(representation_BasisState), public, deferred, nopass :: &
       & representation
  end type
  
  type, extends(BasisState) :: BasisStatePointer
    type(String),                   private :: representation_
    class(BasisState), allocatable, private :: state_
  contains
    procedure, private :: check => check_BasisStatePointer
    
    procedure, public, nopass :: representation => &
                               & representation_BasisStatePointer
    
    procedure, public :: state => state_BasisStatePointer
    procedure, public :: state_pointer => state_pointer_BasisStatePointer
    
    ! I/O.
    procedure, public :: read  => read_BasisStatePointer
    procedure, public :: write => write_BasisStatePointer
  end type
  
  abstract interface
    impure elemental function representation_BasisState() result(output)
      import String
      implicit none
      
      type(String) :: output
    end function
  end interface
  
  interface
    module function types_BasisState() result(output)
      type(BasisStatePointer), allocatable :: output(:)
    end function
  end interface
  
  interface BasisStatePointer
    ! ----------------------------------------------------------------------
    ! BasisStatePointer methods.
    ! ----------------------------------------------------------------------
    ! Construct a BasisStatePointer from any type which extends BasisState.
    impure elemental module function new_BasisStatePointer(state) result(this) 
      class(BasisState), intent(in) :: state
      type(BasisStatePointer)       :: this
    end function
  end interface
  
  interface
    ! Checks that the pointer has been allocated before it is used.
    module subroutine check_BasisStatePointer(this) 
      class(BasisStatePointer), intent(in) :: this
    end subroutine
  end interface
  
  interface
    ! Type representation.
    impure elemental module function representation_BasisStatePointer() &
       & result(output) 
      type(String) :: output
    end function
  end interface
  
  interface
    ! BasisState methods.
    module function state_BasisStatePointer(this) result(output) 
      class(BasisStatePointer), intent(in) :: this
      class(BasisState), allocatable       :: output
    end function
  end interface
  
  interface
    module function state_pointer_BasisStatePointer(this) result(output) 
      class(BasisStatePointer), intent(in), target :: this
      class(BasisState), pointer                   :: output
    end function
  end interface
  
  interface
    ! I/O.
    module subroutine read_BasisStatePointer(this,input) 
      class(BasisStatePointer), intent(out) :: this
      type(String),             intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_BasisStatePointer(this) result(output) 
      class(BasisStatePointer), intent(in) :: this
      type(String), allocatable            :: output(:)
    end function
  end interface
  
  interface BasisStatePointer
    module function new_BasisStatePointer_Strings(input) result(this) 
      type(String), intent(in) :: input(:)
      type(BasisStatePointer)  :: this
    end function
  
    impure elemental module function new_BasisStatePointer_StringArray(input) &
       & result(this) 
      type(StringArray), intent(in) :: input
      type(BasisStatePointer)       :: this
    end function
  end interface
end module
