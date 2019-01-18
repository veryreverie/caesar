! ======================================================================
! A test space, for temporary checking of misc. parts of the code.
! ======================================================================
module TAbs_module
  use common_module
  implicit none
  
  private
  
  public :: TAbs
  public :: TPtr
  
  type, abstract, extends(Stringable) :: TAbs
  contains
    procedure(representation_TAbs), public, deferred, nopass :: representation
    
    procedure, public :: update
  end type
  
  type, extends(TAbs) :: TPtr
    class(TAbs), allocatable :: pointer_
    type(String)             :: representation_
  contains
    procedure, public, nopass :: representation => representation_TPtr
    
    procedure, public :: read  => read_TPtr
    procedure, public :: write => write_TPtr
  end type
  
  type(TPtr), allocatable :: TYPES_TAbs(:)
  
  abstract interface
    impure elemental function representation_TAbs() result(output)
      import String
      implicit none
      
      type(String) :: output
    end function
  end interface
  
  interface TPtr
    module procedure new_TPtr
    module procedure new_TPtr_String
  end interface
contains

subroutine update(this)
  implicit none
  
  class(TAbs), intent(in) :: this
  
  integer :: i
  
  if (.not.allocated(TYPES_TAbs)) then
    TYPES_TAbs = [TPtr(this)]
  elseif (.not.any(this%representation()==[(TYPES_TAbs(i)%pointer_%representation(),i=1,size(TYPES_TAbs))])) then
    TYPES_TAbs = [TYPES_TAbs, TPtr(this)]
  endif
end subroutine

impure elemental function new_TPtr(input) result(output)
  implicit none
  
  class(TAbs), intent(in) :: input
  type(TPtr)              :: output
  
  select type(input); type is(TPtr)
    output = input
  class default
    allocate(output%pointer_, source=input)
    output%representation_ = input%representation()
  end select
end function

impure elemental function representation_TPtr() result(output)
  implicit none
  
  type(String) :: output
  
  output = 'TPtr'
end function

subroutine read_TPtr(this,input)
  implicit none
  
  class(TPtr),  intent(out) :: this
  type(String), intent(in)  :: input
  
  type(String), allocatable :: line(:)
  
  type(String) :: representation
  type(String) :: contents
  
  integer :: i
  
  select type(this); type is(TPtr)
    line = split_line(input)
    
    representation = line(1)
    contents = line(2)
    
    i = first(representation==[(TYPES_TAbs(i)%pointer_%representation(),i=1,size(TYPES_TAbs))])
    call TYPES_TAbs(i)%pointer_%read(contents)
    this = TPtr(TYPES_TAbs(i)%pointer_)
  end select
end subroutine

function write_TPtr(this) result(output)
  implicit none
  
  class(TPtr), intent(in) :: this
  type(String)            :: output
  
  select type(this); type is(TPtr)
    output = this%representation_//' '//str(this%pointer_)
  end select
end function

impure elemental function new_TPtr_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(TPtr)               :: this
  
  call this%read(input)
end function
end module

module TCon_module
  use common_module
  
  use TAbs_module
  implicit none
  
  private
  
  public :: TCon1
  public :: TCon2
  
  type, extends(TAbs) :: TCon1
    integer :: contents_
  contains
    procedure, public, nopass :: representation => representation_TCon1
    
    procedure, public :: read  => read_TCon1
    procedure, public :: write => write_TCon1
  end type
  
  interface TCon1
    module procedure new_TCon1
    module procedure new_TCon1_String
  end interface
  
  type, extends(TAbs) :: TCon2
    integer :: contents_
  contains
    procedure, public, nopass :: representation => representation_TCon2
    
    procedure, public :: read  => read_TCon2
    procedure, public :: write => write_TCon2
  end type
  
  interface TCon2
    module procedure new_TCon2
    module procedure new_TCon2_String
  end interface
contains

impure elemental function new_TCon1(contents) result(this)
  implicit none
  
  integer, intent(in) :: contents
  type(TCon1)         :: this
  
  this%contents_ = contents
end function

impure elemental function new_TCon2(contents) result(this)
  implicit none
  
  integer, intent(in) :: contents
  type(TCon2)         :: this
  
  this%contents_ = contents
end function

impure elemental function representation_TCon1() result(output)
  implicit none
  
  type(String) :: output
  
  output = 'TCon1'
end function

impure elemental function representation_TCon2() result(output)
  implicit none
  
  type(String) :: output
  
  output = 'TCon2'
end function

subroutine read_TCon1(this,input)
  implicit none
  
  class(TCon1), intent(out) :: this
  type(String), intent(in)  :: input
  
  integer :: contents
  
  select type(this); type is(TCon1)
    contents = int(input)
    this = TCon1(contents)
  end select
end subroutine

subroutine read_TCon2(this,input)
  implicit none
  
  class(TCon2), intent(out) :: this
  type(String), intent(in)  :: input
  
  integer :: contents
  
  select type(this); type is(TCon2)
    contents = int(input)
    this = TCon2(contents)
  end select
end subroutine

function write_TCon1(this) result(output)
  implicit none
  
  class(TCon1), intent(in) :: this
  type(String)             :: output
  
  select type(this); type is(TCon1)
    output = str(this%contents_)
  end select
end function

function write_TCon2(this) result(output)
  implicit none
  
  class(TCon2), intent(in) :: this
  type(String)             :: output
  
  select type(this); type is(TCon2)
    output = str(this%contents_)
  end select
end function

impure elemental function new_TCon1_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(TCon1)              :: this
  
  call this%read(input)
end function

impure elemental function new_TCon2_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(TCon2)              :: this
  
  call this%read(input)
end function
end module

module test_module
  use common_module
  
  use TAbs_module
  use TCon_module
  implicit none
  
  private
  
  public :: test
contains

! ----------------------------------------------------------------------
! Generates keywords and helptext.
! ----------------------------------------------------------------------
function test() result(output)
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'test'
  output%description = 'Runs temporary code for testing purposes.'
  output%keywords = [KeywordData::]
  output%main_subroutine => test_subroutine
  output%suppress_from_helptext = .true.
end function

! ----------------------------------------------------------------------
! Main function.
! ----------------------------------------------------------------------
subroutine test_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  integer :: i
  
  class(TAbs), allocatable :: a_
  
  type(TPtr) :: c_
  
  class(TAbs), allocatable :: as_(:)
  class(TPtr), allocatable :: cs_(:)
  
  type(String), allocatable :: strings(:)
  
  type(TCon1) :: t1
  type(TCon2) :: t2
  
  call t1%update()
  call t2%update()
  
  call print_line('')
  allocate(a_, source=TCon1(1))
  call print_line(a_)
  
  call print_line('')
  select type(a_); type is(TCon1)
    call print_line('Type is TCon1')
  type is(TCon2)
    call print_line('Type is TCon2')
  end select
  
  call print_line('')
  c_ = TPtr(a_)
  call print_line(c_)
  
  call print_line('')
  allocate(as_(4), source=TCon1(1))
  do i=1,size(as_)
    call print_line(as_(i))
  enddo
  
  call print_line('')
  deallocate(as_)
  allocate(as_, source=[TCon1(1),TCon1(2)])
  do i=1,size(as_)
    call print_line(as_(i))
  enddo
  
  call print_line('')
  cs_ = [TPtr(TCon1(1)), TPtr(TCon1(2)), TPtr(TCon2(3))]
  do i=1,size(cs_)
    call print_line(cs_(i))
  enddo
  
  call print_line('')
  strings = str(cs_)
  call print_lines(strings)
  
  call print_line('')
  cs_ = TPtr(strings)
  call print_lines(cs_)
end subroutine
end module
