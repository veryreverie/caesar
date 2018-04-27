! ======================================================================
! Input file.
! ======================================================================
module ifile_submodule
  use precision_module
  
  use error_submodule
  use string_submodule
  use io_submodule
  use file_submodule
  use string_array_submodule
  implicit none
  
  private
  
  public :: IFile
  public :: size
  
  type :: IFile
    type(String), private              :: filename_
    type(String), private, allocatable :: lines_(:)
  contains
    generic, public :: assignment(=) => read_file_character, &
                                      & read_file_String
    procedure, private ::               read_file_character
    procedure, private ::               read_file_String
    
    procedure, public :: line
    procedure, public :: split_line
    procedure, public :: lines
    procedure, public :: split_by_blank_lines
  end type
  
  interface size
    module procedure size_IFile
  end interface
contains

! Read a file from its filename.
subroutine read_file_character(this,filename)
  implicit none
  
  class(IFile), intent(out) :: this
  character(*), intent(in)  :: filename
  
  integer         :: file_length
  integer         :: file_unit
  character(1000) :: line
  
  integer :: i,ierr,ialloc
  
  this%filename_ = filename
  
  file_length = count_lines(filename)
  
  allocate(this%lines_(file_length),stat=ialloc); call err(ialloc)
  
  file_unit = open_read_file(filename)
  do i=1,file_length
    read(file_unit,'(a)',iostat=ierr) line
    if (ierr/=0) then
      call print_line('Error reading from '//filename)
      call err()
    endif
    this%lines_(i) = trim(line)
  enddo
  close(file_unit)
end subroutine

subroutine read_file_String(this,filename)
  implicit none
  
  class(IFile), intent(out) :: this
  type(String), intent(in)  :: filename
  
  this = char(filename)
end subroutine

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

! Returns a line from the file, and splits it.
function split_line(this,line_number,delimiter) result(output)
  implicit none
  
  class(IFile), intent(in)  :: this
  integer,      intent(in)  :: line_number
  character,    optional    :: delimiter
  type(String), allocatable :: output(:)
  
  type(String) :: line
  
  line = this%line(line_number)
  
  if (present(delimiter)) then
    output = split(line,delimiter)
  else
    output = split(line)
  endif
end function

! Returns an array of lines from the file.
function lines(this,first_line_number,last_line_number) result(output)
  implicit none
  
  class(IFile), intent(in)  :: this
  integer,      intent(in)  :: first_line_number
  integer,      intent(in)  :: last_line_number
  type(String), allocatable :: output(:)
  
  output = this%lines_(first_line_number:last_line_number)
end function

! Returns the number of lines in a file.
function count_lines(filename) result(output)
  implicit none
  
  character(*), intent(in) :: filename
  integer                  :: output
  
  integer      :: file_unit
  integer      :: iostat
  character(1) :: line
  
  file_unit = open_read_file(filename)
  output = 0
  iostat = 0
  do while (iostat==0)
    read(file_unit, '(a)', iostat=iostat) line
    if (iostat==0) then
      output = output+1
    elseif (iostat>0) then
      call print_line('Error counting lines of '//filename)
      call err()
    endif
  enddo
  close(file_unit)
end function

! Splits the file into sections, split by one or more blank lines.
function split_by_blank_lines(this) result(output)
  implicit none
  
  class(IFile), intent(in)       :: this
  type(StringArray), allocatable :: output(:)
  
  logical, allocatable :: blank_line(:)
  
  integer :: no_sections
  logical :: reading_section
  
  integer, allocatable :: first_lines(:)
  integer, allocatable :: last_lines(:)
  
  integer :: i,ialloc
  
  allocate( first_lines(size(this)), &
          & last_lines(size(this)),  &
          & stat=ialloc); call err(ialloc)
  no_sections = 0
  reading_section = .false.
  do i=1,size(this)
    if (len(this%line(i))==0) then
      ! This line is blank.
      ! If reading a section, then the end of that section is the line above.
      if (reading_section) then
        last_lines(no_sections) = i-1
        reading_section = .false.
      endif
    else
      ! This line is not blank.
      ! If not reading a section, then this line is the start of a new section.
      if (.not. reading_section) then
        no_sections = no_sections+1
        first_lines(no_sections) = i
        reading_section = .true.
      endif
    endif
  enddo
  
  ! If a section is still being read, then that section ends on the last line
  !    of the file.
  if (reading_section) then
    last_lines(no_sections) = size(this)
  endif
  
  allocate(output(no_sections), stat=ialloc); call err(ialloc)
  do i=1,no_sections
    output(i) = StringArray(this%lines_(first_lines(i):last_lines(i)))
  enddo
end function
end module  
