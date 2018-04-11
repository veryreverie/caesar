! ======================================================================
! Output file.
! ======================================================================
module ofile_submodule
  use precision_module
  use string_submodule
  use stringable_submodule
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
    generic, public :: assignment(=) => open_character, &
                                      & open_String
    procedure, private ::               open_character
    procedure, private ::               open_String
    
    procedure, public :: make_stdout
    
    generic, public :: print_line => print_line_character,  &
                                   & print_line_String,     &
                                   & print_line_Stringable, &
                                   & print_line_Printable, &
                                   & print_line_integer,    &
                                   & print_line_real,       &
                                   & print_line_logical,    &
                                   & print_line_complex,    &
                                   & print_line_integers,   &
                                   & print_line_reals,      &
                                   & print_line_logicals,   &
                                   & print_line_complexes
    procedure, private ::            print_line_character
    procedure, private ::            print_line_String
    procedure, private ::            print_line_Stringable
    procedure, private ::            print_line_Printable
    procedure, private ::            print_line_integer
    procedure, private ::            print_line_real
    procedure, private ::            print_line_logical
    procedure, private ::            print_line_complex
    procedure, private ::            print_line_integers
    procedure, private ::            print_line_reals
    procedure, private ::            print_line_logicals
    procedure, private ::            print_line_complexes
    
    final :: finalizer
  end type
contains

! Opens a file for writing to, from its filename.
subroutine open_character(this,filename)
  implicit none
  
  class(OFile), intent(out) :: this
  character(*), intent(in)  :: filename
  
  this%open_ = .true.
  this%filename_ = filename
  this%file_unit_ = open_write_file(filename)
end subroutine

subroutine open_String(this,filename)
  implicit none
  
  class(OFile), intent(out) :: this
  type(String), intent(in)  :: filename
  
  this = char(filename)
end subroutine

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
subroutine print_line_character(this,line)
  implicit none
  
  class(OFile), intent(inout) :: this
  character(*), intent(in)    :: line
  
  integer :: ierr
  
  if (.not. this%open_) then
    call print_line('Code Error: attempted to write to a file which has &
       &either not been opened or has already been closed.')
    call err()
  endif
  
  write(this%file_unit_,'(a)',iostat=ierr) line
  
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

subroutine print_line_String(this,line)
  implicit none
  
  class(OFile), intent(inout) :: this
  type(String), intent(in)    :: line
  
  call this%print_line(char(line))
end subroutine

subroutine print_line_Stringable(this,line)
  implicit none
  
  class(OFile),      intent(inout) :: this
  class(Stringable), intent(in)    :: line
  
  call this%print_line(str(line))
end subroutine

subroutine print_line_Printable(this,line)
  use printable_submodule
  implicit none
  
  class(OFile),     intent(inout) :: this
  class(Printable), intent(in)    :: line
  
  type(String), allocatable :: lines(:)
  integer                   :: i
  
  lines = str(line)
  do i=1,size(lines)
    call this%print_line(lines(i))
  enddo
end subroutine

subroutine print_line_integer(this,line)
  implicit none
  
  class(OFile), intent(inout) :: this
  integer,      intent(in)    :: line
  
  call this%print_line(str(line))
end subroutine

subroutine print_line_real(this,line)
  implicit none
  
  class(OFile), intent(inout) :: this
  real(dp),     intent(in)    :: line
  
  call this%print_line(str(line))
end subroutine

subroutine print_line_logical(this,line)
  implicit none
  
  class(OFile), intent(inout) :: this
  logical,      intent(in)    :: line
  
  call this%print_line(str(line))
end subroutine

subroutine print_line_complex(this,line)
  implicit none
  
  class(OFile), intent(inout) :: this
  complex(dp),  intent(in)    :: line
  
  call this%print_line(str(line))
end subroutine

subroutine print_line_integers(this,line)
  implicit none
  
  class(OFile), intent(inout) :: this
  integer,      intent(in)    :: line(:)
  
  call this%print_line(join(line))
end subroutine

subroutine print_line_reals(this,line)
  implicit none
  
  class(OFile), intent(inout) :: this
  real(dp),     intent(in)    :: line(:)
  
  call this%print_line(join(line))
end subroutine

subroutine print_line_logicals(this,line)
  implicit none
  
  class(OFile), intent(inout) :: this
  logical,      intent(in)    :: line(:)
  
  call this%print_line(join(line))
end subroutine

subroutine print_line_complexes(this,line)
  implicit none
  
  class(OFile), intent(inout) :: this
  complex(dp),  intent(in)    :: line(:)
  
  call this%print_line(join(line))
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
