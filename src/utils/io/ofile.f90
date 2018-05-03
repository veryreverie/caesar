! ======================================================================
! Output file.
! ======================================================================
module ofile_submodule
  use precision_module
  
  use error_submodule
  use string_submodule
  use print_submodule
  use intrinsics_submodule
  use string_writeable_submodule
  use printable_submodule
  use io_submodule
  use file_submodule
  implicit none
  
  private
  
  public :: OFile
  
  type :: OFile
    logical,      private :: open_ = .false.
    type(String), private :: filename_
    integer,      private :: file_unit_
    logical,      private :: is_stdout_ = .false.
  contains
    procedure, public :: make_stdout
    
    generic, public :: print_line => print_line_character,       &
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
    procedure, private ::            print_line_character
    procedure, private ::            print_line_String
    procedure, private ::            print_line_StringWriteable
    procedure, private ::            print_line_integer
    procedure, private ::            print_line_real
    procedure, private ::            print_line_logical
    procedure, private ::            print_line_complex
    procedure, private ::            print_line_integers
    procedure, private ::            print_line_reals
    procedure, private ::            print_line_logicals
    procedure, private ::            print_line_complexes
    
    generic,   public  :: print_lines => print_lines_Strings_character,     &
                                       & print_lines_Strings_String,        &
                                       & print_lines_StringWriteables_character, &
                                       & print_lines_StringWriteables_String,    &
                                       & print_lines_Printable,             &
                                       & print_lines_Printables_character,  &
                                       & print_lines_Printables_String
    procedure, private ::                print_lines_Strings_character
    procedure, private ::                print_lines_Strings_String
    procedure, private ::                print_lines_StringWriteables_character
    procedure, private ::                print_lines_StringWriteables_String
    procedure, private ::                print_lines_Printable
    procedure, private ::                print_lines_Printables_character
    procedure, private ::                print_lines_Printables_String
    
    final :: finalizer
  end type
  
  interface OFile
    module procedure new_OFile_character
    module procedure new_OFile_String
  end interface
contains

! Opens a file for writing to, from its filename.
function new_OFile_character(filename) result(this)
  implicit none
  
  character(*), intent(in) :: filename
  type(OFile)              :: this
  
  this%open_      = .true.
  this%filename_  = filename
  this%file_unit_ = open_write_file(filename)
end function

function new_OFile_String(filename) result(this)
  implicit none
  
  type(String), intent(in) :: filename
  type(OFile)              :: this
  
  this = OFile(char(filename))
end function

! Makes this file be stdout.
subroutine make_stdout(this)
  implicit none
  
  class(OFile), intent(inout) :: this
  
  if (.not. this%open_) then
    call print_line('Code Error: attempted to point stdout to a file which &
       &has either not been opened or has already been closed.')
    call err()
  endif
  
  call set_output_unit(this%file_unit_)
  this%is_stdout_ = .true.
end subroutine

! Writes a line to the file.
subroutine print_line_character(this,input)
  implicit none
  
  class(OFile), intent(inout) :: this
  character(*), intent(in)    :: input
  
  integer :: ierr
  
  if (.not. this%open_) then
    call print_line('Code Error: attempted to write to a file which has &
       &either not been opened or has already been closed.')
    call err()
  endif
  
  write(this%file_unit_,'(a)',iostat=ierr) input
  
  if (ierr/=0) then
    call print_line('Error in OFile::print_line.')
    call err()
  endif
  
  flush(this%file_unit_,iostat=ierr)
  
  if (ierr/=0) then
    call print_line('Error in OFile::print_line.')
    call err()
  endif
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

! Writes multiple lines to a file.
subroutine print_lines_Strings_character(this,input,separating_line)
  implicit none
  
  class(OFile), intent(inout)        :: this
  type(String), intent(in)           :: input(:)
  character(*), intent(in), optional :: separating_line
  
  integer :: i
  
  do i=1,size(input)
    call this%print_line(input(i))
    if (present(separating_line)) then
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
  
  integer :: i
  
  do i=1,size(input)
    call this%print_line(input(i))
    if (present(separating_line)) then
      call print_line(separating_line)
    endif
  enddo
end subroutine

subroutine print_lines_StringWriteables_String(this,input,separating_line)
  implicit none
  
  class(OFile),           intent(inout) :: this
  class(StringWriteable), intent(in)    :: input(:)
  type(String),           intent(in)    :: separating_line
  
  call this%print_lines(input,char(separating_line))
end subroutine

subroutine print_lines_Printable(this,input)
  implicit none
  
  class(OFile),     intent(inout) :: this
  class(Printable), intent(in)    :: input
  
  type(String), allocatable :: lines(:)
  integer                   :: i
  
  lines = str(input)
  do i=1,size(lines)
    call this%print_line(lines(i))
  enddo
end subroutine

subroutine print_lines_Printables_character(this,input,separating_line)
  implicit none
  
  class(OFile),     intent(inout)        :: this
  class(Printable), intent(in)           :: input(:)
  character(*),     intent(in), optional :: separating_line
  
  integer :: i
  
  do i=1,size(input)
    call this%print_lines(input(i))
    if (present(separating_line)) then
      call print_line(separating_line)
    endif
  enddo
end subroutine

subroutine print_lines_Printables_String(this,input,separating_line)
  implicit none
  
  class(OFile),     intent(inout) :: this
  class(Printable), intent(in)    :: input(:)
  type(String),     intent(in)    :: separating_line
  
  call this%print_lines(input,char(separating_line))
end subroutine

! ----------------------------------------------------------------------
! Handles closing the file (if open) when the OFile is finalized.
! Redirects stdout if relevant.
! ----------------------------------------------------------------------
subroutine finalizer(this)
  implicit none
  
  type(OFile), intent(inout) :: this
  
  integer :: ierr
  
  ! Reset stdout.
  if (this%is_stdout_) then
    call unset_output_unit()
    this%is_stdout_ = .false.
  endif
  
  ! Close file (if open).
  if (this%open_) then
    close(this%file_unit_,iostat=ierr)
    
    if (ierr/=0) then
      call print_line('Error: could not close file.')
      call err()
    endif
  endif
end subroutine
end module
