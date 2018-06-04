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
  
  type, extends(Stringable) :: TestType
    integer, allocatable :: contents(:)
  contains
    procedure, public :: read  => read_TestType
    procedure, public :: write => write_TestType
  end type
  
  interface TestType
    module procedure new_TestType
    module procedure new_TestType_String
  end interface
contains

function new_TestType(input) result(this)
  implicit none
  
  integer, intent(in) :: input(:)
  type(TestType)      :: this
  
  this%contents = input
end function

subroutine read_TestType(this,input)
  implicit none
  
  class(TestType), intent(out) :: this
  type(String),    intent(in)  :: input
  
  select type(this); type is(TestType)
    this = TestType(int(split_line(input)))
  end select
end subroutine

function write_TestType(this) result(output)
  implicit none
  
  class(TestType), intent(in) :: this
  type(String)                :: output
  
  select type(this); type is(TestType)
    output = join(this%contents)
  end select
end function

impure elemental function new_TestType_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(TestType)           :: this
  
  this = input
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

subroutine write_test(filename,test1)
  implicit none
  
  type(String),   intent(in) :: filename
  type(TestType), intent(in) :: test1(:)
  
  type(OFile) :: o_file
  
  o_file = OFile(filename)
  call o_file%print_lines(test1)
end subroutine

! ----------------------------------------------------------------------
! Main function.
! ----------------------------------------------------------------------
subroutine test(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  type(String) :: wd
  
  type(TestType), allocatable :: test1(:)
  type(OFile)                 :: o_file
  type(IFile)                 :: i_file
  type(TestType), allocatable :: test2(:)
  
  integer :: a,b,c
  
  wd = arguments%value('working_directory')
  
  a = 3
  b = 4
  c = 7
  
  test1 = [ TestType([a,b,c]),   &
          & TestType([c,b,b,b]), &
          & TestType([a,c,a])    ]
  
  call write_test(wd//'/file.dat', test1)
  
  i_file = IFile(wd//'/file.dat')
  test2 = TestType(i_file%lines())
  call print_lines(test2)
end subroutine
end module
