! ======================================================================
! Output file handling. (Second of two such modules.)
! ======================================================================
! Multiple OFiles can point to a single OFileTarget, and will all write to
!    the same file.
! Uses a SharedCounter to keep track of how many OFiles point to each
!    OFileTarget, and calls close() when the last is deallocated.
module ofile_submodule
  use precision_module
  use io_basic_module
  use abstract_module
  
  use ofile_target_submodule
  use string_writeable_submodule
  use strings_writeable_submodule
  implicit none
  
  private
  
  public :: OFile
  public :: assignment(=)
  
  type :: OFile
    type(OFileTarget),   pointer     :: ofile_target
    type(SharedCounter), allocatable :: counter
  contains
    final :: final_OFile
    
    procedure, private :: check_associated
    
    procedure, public :: make_stdout
    
    generic,   public  :: print_line =>               &
                        & print_line_character,       &
                        & print_line_String,          &
                        & print_line_StringWriteable, &
                        & print_line_integer,         &
                        & print_line_real,            &
                        & print_line_logical,         &
                        & print_line_complex,         &
                        & print_line_integers,        &
                        & print_line_reals,           &
                        & print_line_logicals,        &
                        & print_line_complexes
    procedure, private :: print_line_character
    procedure, private :: print_line_String
    procedure, private :: print_line_StringWriteable
    procedure, private :: print_line_integer
    procedure, private :: print_line_real
    procedure, private :: print_line_logical
    procedure, private :: print_line_complex
    procedure, private :: print_line_integers
    procedure, private :: print_line_reals
    procedure, private :: print_line_logicals
    procedure, private :: print_line_complexes
    
    generic,   public  :: print_lines =>                          &
                        & print_lines_Strings_character,          &
                        & print_lines_Strings_String,             &
                        & print_lines_StringWriteables_character, &
                        & print_lines_StringWriteables_String,    &
                        & print_lines_StringsWriteable,           &
                        & print_lines_StringsWriteables_character,&
                        & print_lines_StringsWriteables_String
    procedure, private :: print_lines_Strings_character
    procedure, private :: print_lines_Strings_String
    procedure, private :: print_lines_StringWriteables_character
    procedure, private :: print_lines_StringWriteables_String
    procedure, private :: print_lines_StringsWriteable
    procedure, private :: print_lines_StringsWriteables_character
    procedure, private :: print_lines_StringsWriteables_String
  end type
  
  interface OFile
    module procedure new_OFile_character
    module procedure new_OFile_String
  end interface
  
  interface assignment(=)
    module procedure assign_OFile_OFile
  end interface
contains

! Constructor, assignment and finalization.
function new_OFile_character(filename) result(this)
  implicit none
  
  character(*), intent(in) :: filename
  type(OFile)              :: this
  
  integer :: ialloc
  
  type(OFileTarget) :: ofile_target
  
  ofile_target = OFileTarget(filename)
  
  allocate(this%ofile_target, stat=ialloc); call err(ialloc)
  this%ofile_target = ofile_target
  
  allocate(this%counter, stat=ialloc); call err(ialloc)
  this%counter = SharedCounter()
end function

function new_OFile_String(filename) result(this)
  implicit none
  
  type(String), intent(in) :: filename
  type(OFile)              :: this
  
  this = OFile(char(filename))
end function

subroutine assign_OFile_OFile(output,input)
  implicit none
  
  type(OFile), intent(out) :: output
  type(OFile), intent(in)  :: input
  
  integer :: ialloc
  
  output%ofile_target => input%ofile_target
  
  allocate(output%counter, stat=ialloc); call err(ialloc)
  output%counter = input%counter
end subroutine

subroutine final_OFile(this)
  implicit none
  
  type(OFile), intent(inout) :: this
  
  integer :: ialloc
  
  if (allocated(this%counter)) then
    if (this%counter%is_only_pointer()) then
      call this%ofile_target%close()
      deallocate(this%ofile_target, stat=ialloc); call err(ialloc)
    endif
    deallocate(this%counter, stat=ialloc); call err(ialloc)
  endif
end subroutine

! Checks that this OFile has been associated with a file.
subroutine check_associated(this)
  implicit none
  
  class(OFile), intent(in) :: this
  
  if (.not. associated(this%ofile_target)) then
    call print_line(CODE_ERROR//': Attempting to call output file operations &
       &on OFile which has not been associated.')
    call err()
  endif
end subroutine

! Calls functionality of the associated file.
subroutine make_stdout(this)
  implicit none
  
  class(OFile), intent(inout) :: this
  
  call this%check_associated()
  call this%ofile_target%make_stdout()
end subroutine

subroutine print_line_character(this,input)
  implicit none
  
  class(OFile), intent(inout) :: this
  character(*), intent(in)    :: input
  
  call this%check_associated()
  call this%ofile_target%print_line(input)
end subroutine

subroutine print_line_String(this,input)
  implicit none
  
  class(OFile), intent(inout) :: this
  type(String), intent(in)    :: input
  
  call this%print_line(char(input))
end subroutine

subroutine print_line_StringWriteable(this,input)
  implicit none
  
  class(OFile),           intent(inout) :: this
  class(StringWriteable), intent(in)    :: input
  
  call this%print_line(str(input))
end subroutine

subroutine print_line_integer(this,input)
  implicit none
  
  class(OFile), intent(inout) :: this
  integer,      intent(in)    :: input
  
  call this%print_line(str(input))
end subroutine

subroutine print_line_real(this,input)
  implicit none
  
  class(OFile), intent(inout) :: this
  real(dp),     intent(in)    :: input
  
  call this%print_line(str(input))
end subroutine

subroutine print_line_logical(this,input)
  implicit none
  
  class(OFile), intent(inout) :: this
  logical,      intent(in)    :: input
  
  call this%print_line(str(input))
end subroutine

subroutine print_line_complex(this,input)
  implicit none
  
  class(OFile), intent(inout) :: this
  complex(dp),  intent(in)    :: input
  
  call this%print_line(str(input))
end subroutine

subroutine print_line_integers(this,input)
  implicit none
  
  class(OFile), intent(inout) :: this
  integer,      intent(in)    :: input(:)
  
  call this%print_line(join(input))
end subroutine

subroutine print_line_reals(this,input)
  implicit none
  
  class(OFile), intent(inout) :: this
  real(dp),     intent(in)    :: input(:)
  
  call this%print_line(join(input))
end subroutine

subroutine print_line_logicals(this,input)
  implicit none
  
  class(OFile), intent(inout) :: this
  logical,      intent(in)    :: input(:)
  
  call this%print_line(join(input))
end subroutine

subroutine print_line_complexes(this,input)
  implicit none
  
  class(OFile), intent(inout) :: this
  complex(dp),  intent(in)    :: input(:)
  
  call this%print_line(join(input))
end subroutine

subroutine print_lines_Strings_character(this,input,separating_line)
  implicit none
  
  class(OFile), intent(inout)        :: this
  type(String), intent(in)           :: input(:)
  character(*), intent(in), optional :: separating_line
  
  integer :: i
  
  do i=1,size(input)
    call this%print_line(input(i))
    if (present(separating_line) .and. i<size(input)) then
      call print_line(separating_line)
    endif
  enddo
end subroutine

subroutine print_lines_Strings_String(this,input,separating_line)
  implicit none
  
  class(OFile), intent(inout) :: this
  type(String), intent(in)    :: input(:)
  type(String), intent(in)    :: separating_line
  
  call this%print_lines(input,char(separating_line))
end subroutine

subroutine print_lines_StringWriteables_character(this,input,separating_line)
  implicit none
  
  class(OFile),           intent(inout)        :: this
  class(StringWriteable), intent(in)           :: input(:)
  character(*),           intent(in), optional :: separating_line
  
  call this%print_lines(str(input,separating_line))
end subroutine

subroutine print_lines_StringWriteables_String(this,input,separating_line)
  implicit none
  
  class(OFile),           intent(inout) :: this
  class(StringWriteable), intent(in)    :: input(:)
  type(String),           intent(in)    :: separating_line
  
  call this%print_lines(str(input,separating_line))
end subroutine

subroutine print_lines_StringsWriteable(this,input)
  implicit none
  
  class(OFile),            intent(inout) :: this
  class(StringsWriteable), intent(in)    :: input
  
  call this%print_lines(str(input))
end subroutine

subroutine print_lines_StringsWriteables_character(this,input,separating_line)
  implicit none
  
  class(OFile),            intent(inout)        :: this
  class(StringsWriteable), intent(in)           :: input(:)
  character(*),            intent(in), optional :: separating_line
  
  call this%print_lines(str(input,separating_line))
end subroutine

subroutine print_lines_StringsWriteables_String(this,input,separating_line)
  implicit none
  
  class(OFile),            intent(inout) :: this
  class(StringsWriteable), intent(in)    :: input(:)
  type(String),            intent(in)    :: separating_line
  
  call this%print_lines(str(input,separating_line))
end subroutine
end module
