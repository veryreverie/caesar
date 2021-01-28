! ======================================================================
! Input file.
! ======================================================================
module caesar_ifile_module
  use caesar_precision_module
  use caesar_io_basic_module
  
  use caesar_file_module
  use caesar_string_array_module
  implicit none
  
  private
  
  public :: IFile
  public :: size
  
  type :: IFile
    type(String), private              :: filename_
    type(String), private, allocatable :: lines_(:)
  contains
    procedure, public :: line
    
    generic,   public  :: lines =>   &
                        & lines_all, &
                        & lines_slice
    procedure, private :: lines_all
    procedure, private :: lines_slice
    
    generic,   public  :: sections =>         &
                        & sections_character, &
                        & sections_String
    procedure, private :: sections_character
    procedure, private :: sections_String
  end type
  
  interface IFile
    module procedure new_IFile_character
    module procedure new_IFile_String
  end interface
  
  interface size
    module procedure size_IFile
  end interface
contains

! Read a file from its filename.
function new_IFile_character(filename) result(this)
  implicit none
  
  character(*), intent(in) :: filename
  type(IFile)              :: this
  
  integer          :: file_length
  integer          :: file_unit
  character(10000) :: line
  
  integer :: i,ierr,ialloc
  
  if (.not. file_exists(filename)) then
    call print_line(ERROR//': file does not exist: '//filename)
    call err()
  endif
  
  this%filename_ = filename
  
  file_length = count_lines(filename)
  
  allocate(this%lines_(file_length),stat=ialloc); call err(ialloc)
  
  file_unit = open_read_file(filename)
  do i=1,file_length
    read(file_unit,'(a)',iostat=ierr) line
    if (ierr/=0) then
      call print_line(ERROR//': failed to read file: '//filename)
      call err()
    endif
    this%lines_(i) = trim(line)
  enddo
  close(file_unit)
end function

function new_IFile_String(filename) result(this)
  implicit none
  
  type(String), intent(in) :: filename
  type(IFile)              :: this
  
  this = IFile(char(filename))
end function

! The number of lines in the file.
function size_IFile(this) result(output)
  implicit none
  
  type(IFile), intent(in) :: this
  integer                 :: output
  
  output = size(this%lines_)
end function

! Returns a line from the file.
function line(this,line_number) result(output)
  implicit none
  
  class(IFile), intent(in) :: this
  integer,      intent(in) :: line_number
  type(String)             :: output
  
  output = this%lines_(line_number)
end function

! Returns an array of lines from the file.
function lines_all(this) result(output)
  implicit none
  
  class(IFile), intent(in)  :: this
  type(String), allocatable :: output(:)
  
  output = this%lines_
end function

function lines_slice(this,first_line_number,last_line_number) result(output)
  implicit none
  
  class(IFile), intent(in)  :: this
  integer,      intent(in)  :: first_line_number
  integer,      intent(in)  :: last_line_number
  type(String), allocatable :: output(:)
  
  output = this%lines_(first_line_number:last_line_number)
end function

! Equivalent to split_into_sections(this%lines(),separating_line).
function sections_character(this,separating_line) result(output)
  implicit none
  
  class(IFile), intent(in)           :: this
  character(*), intent(in), optional :: separating_line
  type(StringArray), allocatable     :: output(:)
  
  type(String), allocatable :: lines(:)
  
  lines = this%lines()
  
  output = split_into_sections(lines, separating_line)
end function

function sections_String(this,separating_line) result(output)
  implicit none
  
  class(IFile), intent(in)       :: this
  type(String), intent(in)       :: separating_line
  type(StringArray), allocatable :: output(:)
  
  output = this%sections(char(separating_line))
end function

! Returns the number of lines in a file.
function count_lines(filename) result(output)
  implicit none
  
  character(*), intent(in) :: filename
  integer                  :: output
  
  integer      :: file_unit
  integer      :: iostat
  
  file_unit = open_read_file(filename)
  output = 0
  iostat = 0
  do while (iostat==0)
    read(file_unit, '(a)', iostat=iostat)
    if (iostat==0) then
      output = output+1
    elseif (iostat>0) then
      call print_line('Error counting lines of '//filename)
      call err()
    endif
  enddo
  close(file_unit)
end function
end module  
