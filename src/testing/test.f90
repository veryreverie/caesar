! ======================================================================
! A test space, for temporary checking of misc. parts of the code.
! ======================================================================
module test_module
  use common_module
  implicit none
  
  private
  
  public :: test_keywords
  public :: test_mode
  public :: test
  
  type, abstract, extends(NoDefaultConstructor) :: Base
  contains
    procedure(shout_Base), public, deferred :: shout
  end type
  
  abstract interface
    subroutine shout_Base(this)
      import Base
      implicit none
      
      class(Base), intent(in) :: this
    end subroutine
  end interface
  
  type, extends(Base) :: Concrete1
    integer :: contents
  contains
    procedure, public :: shout => shout_Concrete1
  end type
  
  interface Concrete1
    module procedure new_Concrete1
  end interface
  
  type, extends(Base) :: Concrete2
    integer :: contents
  contains
    procedure, public :: shout => shout_Concrete2
  end type
  
  interface Concrete2
    module procedure new_Concrete2
  end interface
  
  type, extends(Base) :: BasePointer
    class(Base), allocatable :: contents
  contains
    procedure, public :: shout => shout_BasePointer
  end type
  
  interface BasePointer
    module procedure new_BasePointer
  end interface
  
  interface assignment(=)
    module procedure assign_BasePointer_Base
  end interface
contains

function new_Concrete1(contents) result(this)
  implicit none
  
  integer, intent(in) :: contents
  type(Concrete1)     :: this
  
  this%contents = contents
end function

function new_Concrete2(contents) result(this)
  implicit none
  
  integer, intent(in) :: contents
  type(Concrete2)     :: this
  
  this%contents = contents
end function

function new_BasePointer(contents) result(this)
  implicit none
  
  class(Base), intent(in) :: contents
  type(BasePointer)       :: this
  
  integer :: ialloc
  
  allocate(this%contents, source=contents, stat=ialloc); call err(ialloc)
end function

subroutine assign_BasePointer_Base(output,input)
  implicit none
  
  type(BasePointer), intent(out) :: output
  class(Base),       intent(in)  :: input
  
  integer :: ialloc
  
  select type(input); class is(BasePointer)
    allocate( output%contents, source=input%contents, &
            & stat=ialloc); call err(ialloc)
  class default
    allocate( output%contents, source=input, &
            & stat=ialloc); call err(ialloc)
  end select
end subroutine

subroutine shout_Concrete1(this)
  implicit none
  
  class(Concrete1), intent(in) :: this
  
  call print_line('Concrete1. Value = '//this%contents)
end subroutine

subroutine shout_Concrete2(this)
  implicit none
  
  class(Concrete2), intent(in) :: this
  
  call print_line('Concrete2. Value = '//this%contents)
end subroutine

subroutine shout_BasePointer(this)
  implicit none
  
  class(BasePointer), intent(in) :: this
  
  call this%contents%shout()
end subroutine

function choose(concrete,contents) result(output)
  implicit none
  
  integer, intent(in)      :: concrete
  integer, intent(in)      :: contents
  type(BasePointer)        :: output
  
  if (concrete==1) then
    output = Concrete1(contents)
  elseif (concrete==2) then
    output = Concrete2(contents)
  else
    call err()
  endif
end function

! ----------------------------------------------------------------------
! Generates keywords and helptext.
! ----------------------------------------------------------------------
function test_keywords() result(keywords)
  implicit none
  
  type(KeywordData), allocatable :: keywords(:)
  
  keywords = [KeywordData::]
end function

function test_mode() result(output)
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'test'
  output%description = 'Runs temporary code for testing purposes.'
  output%keywords = test_keywords()
  output%main_subroutine => test
  output%suppress_from_helptext = .true.
end function

! ----------------------------------------------------------------------
! Main function.
! ----------------------------------------------------------------------
subroutine test(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  type(String) :: wd
  
  type(BasePointer) :: thing
  
  wd = arguments%value('working_directory')
  
  thing = Concrete1(12)
  call thing%shout()
  
  thing = Concrete2(43)
  call thing%shout()
  
  thing = choose(1,42)
  call thing%shout()
  
  thing = choose(2,39)
  call thing%shout()
end subroutine
end module
