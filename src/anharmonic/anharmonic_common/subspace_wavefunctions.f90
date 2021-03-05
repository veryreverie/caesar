! ======================================================================
! An abstract class to hold a printable representation of the wavefunctions
!    spanning a subspace.
! ======================================================================
module caesar_subspace_wavefunctions_module
  use caesar_common_module
  implicit none
  
  private
  
  public :: SubspaceWavefunctions
  public :: SubspaceWavefunctionsPointer
  
  type, abstract, extends(Stringsable) :: SubspaceWavefunctions
  contains
    procedure(representation_SubspaceWavefunctions), public, deferred, &
       & nopass :: representation
  end type
  
  type, extends(SubspaceWavefunctions) :: SubspaceWavefunctionsPointer
    type(String),                              private :: representation_
    class(SubspaceWavefunctions), allocatable, private :: wavefunctions_
  contains
    procedure, public :: check => check_SubspaceWavefunctionsPointer
    
    procedure, public, nopass :: representation => &
                               & representation_SubspaceWavefunctionsPointer
    
    ! I/O.
    procedure, public :: read  => read_SubspaceWavefunctionsPointer
    procedure, public :: write => write_SubspaceWavefunctionsPointer
  end type
  
  ! Abstract SubspaceWavfunction routines.
  abstract interface
    impure elemental function representation_SubspaceWavefunctions() &
       & result(output)
      import String
      implicit none
      
      type(String) :: output
    end function
  end interface
  
  interface
    module function types_SubspaceWavefunctions() result(output)
      type(SubspaceWavefunctionsPointer), allocatable :: output(:)
    end function
  end interface
  
  interface SubspaceWavefunctionsPointer
    ! Construct a SubspaceWavefunctionsPointer from any type which extends
    !    SubspaceWavefunctions.
    impure elemental module function new_SubspaceWavefunctionsPointer( &
       & wavefunctions) result(this) 
      class(SubspaceWavefunctions), intent(in) :: wavefunctions
      type(SubspaceWavefunctionsPointer)       :: this
    end function
  end interface
  
  interface
    ! Checks that the pointer has been allocated before it is used.
    module subroutine check_SubspaceWavefunctionsPointer(this) 
      class(SubspaceWavefunctionsPointer) ,intent(in) :: this
    end subroutine
  end interface
  
  interface
    ! Type representation.
    impure elemental module function &
       & representation_SubspaceWavefunctionsPointer() result(output) 
      type(String) :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_SubspaceWavefunctionsPointer(this,input) 
      class(SubspaceWavefunctionsPointer) ,intent(out) :: this
      type(String),                        intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_SubspaceWavefunctionsPointer(this) result(output) 
      class(SubspaceWavefunctionsPointer) ,intent(in) :: this
      type(String), allocatable                       :: output(:)
    end function
  end interface
  
  interface SubspaceWavefunctionsPointer
    module function new_SubspaceWavefunctionsPointer_Strings(input) &
       & result(this) 
      type(String), intent(in)           :: input(:)
      type(SubspaceWavefunctionsPointer) :: this
    end function
  
    impure elemental module function &
       & new_SubspaceWavefunctionsPointer_StringArray(input) result(this) 
      type(StringArray), intent(in)      :: input
      type(SubspaceWavefunctionsPointer) :: this
    end function
  end interface
end module
