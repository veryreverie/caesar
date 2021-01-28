! ======================================================================
! Output file handling. (First of two such modules.)
! ======================================================================
! OFileTarget keeps track of a single file.
! Multiple OFiles may point to a single OFileTarget.
module caesar_ofile_target_module
  use caesar_precision_module
  use caesar_io_basic_module
  
  use caesar_string_writeable_module
  use caesar_strings_writeable_module
  use caesar_file_module
  implicit none
  
  private
  
  public :: OFileTarget
  
  type :: OFileTarget
    logical,      private :: open_ = .false.
    type(String), private :: filename_
    integer,      private :: file_unit_
    logical,      private :: is_stdout_
  contains
    procedure, public :: close => close_OFileTarget
    
    procedure, public :: is_open
    
    procedure, public :: make_stdout
    
    generic,   public  :: print_line => &
                        & print_line_character
    procedure, private :: print_line_character
  end type
  
  interface OFileTarget
    module procedure new_OFileTarget
  end interface
contains

! Constructor and finalizer. Opens a file on construction, closes the file
!    (and redirects stdout, if relevant) when finalised.
function new_OFileTarget(filename) result(this)
  implicit none
  
  character(*), intent(in) :: filename
  type(OFileTarget)        :: this
  
  this%open_      = .true.
  this%filename_  = filename
  this%file_unit_ = open_write_file(filename)
  this%is_stdout_ = .false.
end function

subroutine close_OFileTarget(this)
  implicit none
  
  class(OFileTarget), intent(inout) :: this
  
  integer :: ierr
  
  ! Check file is open.
  if (.not. this%is_open()) then
    call print_line(ERROR//': Trying to close file which is not open.')
    call err()
  endif
  
  ! Close file
  close(this%file_unit_,iostat=ierr)
  if (ierr/=0) then
    call print_line(ERROR//': could not close file. Error code '//ierr)
    call err()
  endif
  
  ! Reset stdout if necessary.
  if (this%is_stdout_) then
    call unset_output_unit()
    this%is_stdout_ = .false.
  endif
end subroutine

! Returns whether or not this file is pointed to by any OFiles.
function is_open(this) result(output)
  implicit none
  
  class(OFileTarget), intent(in) :: this
  logical                        :: output
  
  output = this%open_
end function

! Makes this file be stdout.
subroutine make_stdout(this)
  implicit none
  
  class(OFileTarget), intent(inout) :: this
  
  if (.not. this%is_open()) then
    call print_line(CODE_ERROR//': attempted to point stdout to a file which &
       &has either not been opened or has already been closed.')
    call err()
  endif
  
  call set_output_unit(this%file_unit_)
  this%is_stdout_ = .true.
end subroutine

! Writes a line to the file.
subroutine print_line_character(this,input,settings)
  implicit none
  
  class(OFileTarget),  intent(inout)        :: this
  type(PrintSettings), intent(in), optional :: settings
  character(*),        intent(in)           :: input
  
  integer :: ierr
  
  if (.not. this%is_open()) then
    call print_line(CODE_ERROR//': attempted to write to a file which has &
       &either not been opened or has already been closed.')
    call err()
  endif
  
  write(this%file_unit_,'(a)',iostat=ierr) input
  
  if (ierr/=0) then
    call print_line(ERROR//': writing to file failed. Error code '//ierr)
    call err()
  endif
  
  flush(this%file_unit_,iostat=ierr)
  
  if (ierr/=0) then
    call print_line(ERROR//': flushing file failed. Error code '//ierr)
    call err()
  endif
end subroutine
end module
