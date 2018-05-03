! ======================================================================
! An abstract type combining StringReadable and StringWriteable.
! ======================================================================
! Any type which extends Stringable has all the properties of StringReadable
!    and StringWriteable.
! Any type which extends Stringable must overload %write() and %read(string).
! See example module below for how to extend this type.
!
! N.B.: Stringable extends StringWriteable, but duplicates the functionality 
!    of StringReadable since Fortran does not support multiple inheritance.
module stringable_submodule
  use string_submodule
  use string_array_submodule
  use string_readable_submodule
  use string_writeable_submodule
  implicit none
  
  private
  
  public :: Stringable
  public :: assignment(=)
  
  type, abstract, extends(StringWriteable) :: Stringable
  contains
    procedure(read_Stringable), deferred :: read
  end type
  
  abstract interface
    recursive subroutine read_Stringable(this,input)
      import String
      import Stringable
      implicit none
      
      class(Stringable), intent(out) :: this
      type(String),      intent(in)  :: input
    end subroutine
  end interface
  
  interface assignment(=)
    module procedure assign_Stringable_String
    module procedure assign_Stringable_character
    module procedure assign_Stringables_Strings
    module procedure assign_Stringables_StringArray
  end interface
contains

! ----------------------------------------------------------------------
! Assign a Stringable type from a String type.
! ----------------------------------------------------------------------
recursive subroutine assign_Stringable_String(output,input)
  implicit none
  
  class(Stringable), intent(out) :: output
  type(String),      intent(in)  :: input
  
  call output%read(input)
end subroutine

recursive subroutine assign_Stringable_character(output,input)
  implicit none
  
  class(Stringable), intent(out) :: output
  character(*),      intent(in)  :: input
  
  output = str(input)
end subroutine

recursive subroutine assign_Stringables_Strings(output,input)
  implicit none
  
  class(Stringable), allocatable, intent(out) :: output(:)
  type(String),                   intent(in)  :: input(:)
  
  integer :: i,ialloc
  
  allocate(output(size(input)), mold=output, stat=ialloc); call err(ialloc)
  do i=1,size(input)
    output(i) = input(i)
  enddo
end subroutine

recursive subroutine assign_Stringables_StringArray(output,input)
  implicit none
  
  class(Stringable), allocatable, intent(out) :: output(:)
  type(StringArray),              intent(in)  :: input
  
  output = input%strings
end subroutine
end module

! ======================================================================
! An example module showing how to extend Stringable.
! ======================================================================
module stringable_example_submodule
  use string_submodule
  use stringable_submodule
  use error_submodule
  implicit none
  
  private
  
  public :: StringableExample
  
  type, extends(Stringable) :: StringableExample
    type(String) :: contents
  contains
    procedure, public :: read  => read_StringableExample
    procedure, public :: write => write_StringableExample
  end type
contains

subroutine read_StringableExample(this,input)
  implicit none
  
  class(StringableExample), intent(out) :: this
  type(String),             intent(in)  :: input
  
  ! Select type needed to call non-polymorphic procedures, and to ensure that
  !    read() is overloaded by any type which extends StringableExample.
  select type(this); type is(StringableExample)
    this%contents = input
  class default
    call print_line(CODE_ERROR//': Called the StringableExample version &
       &of read() from a type other than StringableExample.')
    call err()
  end select
end subroutine

function write_StringableExample(this) result(output)
  implicit none
  
  class(StringableExample), intent(in) :: this
  type(String)                         :: output
  
  ! Select type needed to call non-polymorphic procedures, and to ensure that
  !    write() is overloaded by any type which extends StringableExample.
  select type(this); type is(StringableExample)
    output = str(this%contents)
  class default
    call print_line(CODE_ERROR//': Called the StringableExample version &
       &of write() from a type other than StringableExample.')
    call err()
  end select
end function
end module
