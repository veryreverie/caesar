! ======================================================================
! Provides functionallity for writing to the terminal.
! ======================================================================
! print_line acts as the intrinsic write(*,*), but with the ability to be
!    redirected to a file. Provides error checking and formatting dependent on
!    whether the output is the terminal or a file.
! print_lines acts as print_line, but for objects which are printed onto more
!    than one line.
! colour adds terminal escape characters so that a string is coloured when
!    printed to the terminal.
module print_module
  use iso_fortran_env, only : OUTPUT_UNIT
  use string_module
  use terminal_module
  use error_module
  use print_settings_module
  implicit none
  
  private
  
  public :: print_line
  public :: print_lines
  public :: colour
  public :: set_output_unit
  public :: unset_output_unit
  public :: set_terminal_width
  
  ! Terminal settings.
  integer :: TERMINAL_WIDTH = 79
  integer :: OUTPUT_FILE_UNIT = OUTPUT_UNIT
  
  interface print_line
    module procedure print_line_character
    module procedure print_line_String
  end interface
  
  interface print_lines
    module procedure print_lines_Strings_character
    module procedure print_lines_Strings_String
  end interface
  
  interface colour
    module procedure colour_character_character
    module procedure colour_String_character
    module procedure colour_character_String
    module procedure colour_String_String
  end interface
  
  ! C tput cols interface.
  interface
    function get_terminal_width_c(width) bind(c) result(success)
      use, intrinsic :: iso_c_binding
      implicit none
      
      integer(kind=c_int), intent(out) :: width
      logical(kind=c_bool)             :: success
    end function
  end interface
contains

! ----------------------------------------------------------------------
! Prints a character or string variable, either to the terminal or
!    to the output file specified by OUTPUT_FILE_UNIT.
! ----------------------------------------------------------------------
subroutine print_line_character(line,settings)
  implicit none
  
  character(*),        intent(in)           :: line
  type(PrintSettings), intent(in), optional :: settings
  
  type(PrintSettings) :: print_settings
  
  integer                   :: last_space
  integer                   :: first_none_space
  logical                   :: reading_special
  integer                   :: current_terminal_width
  character(:), allocatable :: indent_spaces
  character(:), allocatable :: overhang_spaces
  
  logical,      allocatable :: is_space(:)
  logical,      allocatable :: is_special(:)
  integer,      allocatable :: positions(:)
  integer                   :: line_start
  type(String), allocatable :: lines(:)
  
  integer :: i,ierr,ialloc
  
  if (present(settings)) then
    print_settings = settings
  else
    print_settings = PrintSettings()
  endif
  
  indent_spaces = spaces(print_settings%indent)
  overhang_spaces = spaces(print_settings%overhang)
  
  ! Break the line up into multiple lines such that line breaks happen at
  !    spaces only.
  if (OUTPUT_FILE_UNIT/=OUTPUT_UNIT) then
    ! If not writing to the terminal, ignore line-breaking.
    lines = [str(line)]
  elseif (len(line)==0) then
    lines = [str(line)]
  else
    ! Identify which characters are spaces.
    is_space = [(line(i:i)==' ', i=1, len(line))]
    
    ! Identify which characters are part of invisible escape strings.
    allocate(is_special(len(line)), stat=ialloc); call err(ialloc)
    reading_special = .false.
    do i=1,len(line)
      if (line(i:i)==ESC) then
        reading_special = .true.
        is_special(i) = .true.
      elseif (reading_special) then
        if (line(i:i)=='m') then
          reading_special = .false.
        endif
        is_special(i) = .true.
      else
        is_special(i) = .false.
      endif
    enddo
    
    ! Identify the position of each character in the visible string,
    !    ignoring invisible escape strings.
    allocate(positions(len(line)), stat=ialloc); call err(ialloc)
    if (is_special(1)) then
      positions(1) = 0
    else
      positions(1) = 1
    endif
    do i=2,size(positions)
      if (is_special(i)) then
        positions(i) = positions(i-1)
      else
        positions(i) = positions(i-1)+1
      endif
    enddo
    
    ! Identify lines.
    current_terminal_width = TERMINAL_WIDTH - print_settings%indent
    line_start = 1
    lines = [String::]
    do
      ! Find the last space which fits onto one terminal line,
      !    and the first none-space character which does not fit.
      last_space = 0
      first_none_space = 0
      do i=1,size(positions)
        if (positions(i)<=current_terminal_width+1 .and. is_space(i)) then
          last_space = i
        elseif (positions(i)>current_terminal_width .and..not.is_space(i)) then
          first_none_space = i
          exit
        endif
      enddo
      
      ! If the line does not need to be split, or cannot be split,
      !    add it in its entirety.
      if (last_space==0 .or. first_none_space==0) then
        lines = [lines, str(line(line_start:len(line)))]
        exit
      endif
      
      ! Find the first none-space character on the next line.
      do i=first_none_space,last_space,-1
        if (is_space(i)) then
          first_none_space = i+1
          exit
        endif
      enddo
      
      ! Add the line, and update variables to process next line.
      lines = [lines, str(line(line_start:line_start+last_space-2))]
      line_start = line_start + first_none_space-1
      is_space = is_space(first_none_space:)
      positions = positions(first_none_space:) - positions(first_none_space-1)
      
      ! Account for the overhang in the terminal width.
      current_terminal_width = TERMINAL_WIDTH        &
                           & - print_settings%indent &
                           & - print_settings%overhang
    enddo
  endif
  
  ! Write lines.
  do i=1,size(lines)
    if (i==1) then
      write(OUTPUT_FILE_UNIT,'(a)',iostat=ierr) &
         & indent_spaces//char(lines(i))
    else
      write(OUTPUT_FILE_UNIT,'(a)',iostat=ierr) &
         & indent_spaces//overhang_spaces//char(lines(i))
    endif
    
    ! Check for errors and flush the write buffer.
    if (ierr /= 0) then
      write(OUTPUT_FILE_UNIT,'(a)') 'Error in print_line.'
      call err()
    endif
    
    flush(OUTPUT_FILE_UNIT,iostat=ierr)
    
    if (ierr /= 0) then
      write(OUTPUT_FILE_UNIT,'(a)') 'Error in print_line.'
      call err()
    endif
  enddo
end subroutine

subroutine print_line_String(line,settings)
  implicit none
  
  type(String),        intent(in)           :: line
  type(PrintSettings), intent(in), optional :: settings
  
  call print_line(char(line), settings)
end subroutine

! ----------------------------------------------------------------------
! As print_line, but for printing over multiple lines.
! ----------------------------------------------------------------------
subroutine print_lines_Strings_character(lines,separating_line,settings)
  implicit none
  
  type(String),        intent(in)           :: lines(:)
  character(*),        intent(in), optional :: separating_line
  type(PrintSettings), intent(in), optional :: settings
  
  integer :: i
  
  do i=1,size(lines)
    call print_line(lines(i), settings)
    if (present(separating_line)) then
      call print_line(separating_line, settings)
    endif
  enddo
end subroutine

subroutine print_lines_Strings_String(lines,separating_line,settings)
  implicit none
  
  type(String),        intent(in)           :: lines(:)
  type(String),        intent(in)           :: separating_line
  type(PrintSettings), intent(in), optional :: settings
  
  call print_lines(lines, char(separating_line), settings)
end subroutine

! ----------------------------------------------------------------------
! Adds terminal colour to a string.
! Does nothing if the output is going to a file.
! ----------------------------------------------------------------------
! N.B. can't detect "caesar > file", only supresses colour if "caesar -o file".
function colour_character_character(input,colour_name) result(output)
  implicit none
  
  character(*), intent(in) :: input
  character(*), intent(in) :: colour_name
  type(String)             :: output
  
  type(String) :: lower_case_name
  character(2) :: colour_code
  
  ! Supress colour if outputting to file.
  if (OUTPUT_FILE_UNIT/=OUTPUT_UNIT) then
    output = input
    return
  endif
  
  lower_case_name = lower_case(colour_name)
  
  if (lower_case_name=='black') then
    colour_code = '30'
  elseif (lower_case_name=='red') then
    colour_code = '31'
  elseif (lower_case_name=='green') then
    colour_code = '32'
  elseif (lower_case_name=='yellow') then
    colour_code = '33'
  elseif (lower_case_name=='blue') then
    colour_code = '34'
  elseif (lower_case_name=='magenta') then
    colour_code = '35'
  elseif (lower_case_name=='cyan') then
    colour_code = '36'
  elseif (lower_case_name=='light gray') then
    colour_code = '37'
  elseif (lower_case_name=='dark gray') then
    colour_code = '90'
  elseif (lower_case_name=='light red') then
    colour_code = '91'
  elseif (lower_case_name=='light green') then
    colour_code = '92'
  elseif (lower_case_name=='light yellow') then
    colour_code = '93'
  elseif (lower_case_name=='light blue') then
    colour_code = '94'
  elseif (lower_case_name=='light magenta') then
    colour_code = '95'
  elseif (lower_case_name=='light cyan') then
    colour_code = '96'
  elseif (lower_case_name=='white') then
    colour_code = '97'
  else
    call print_line('Error: '//colour_name//' is not an accepted colour.')
    call err()
  endif
  
  output = ESC//'['//colour_code//'m'//input//ESC//'[0m'
end function

function colour_String_character(input,colour_name) result(output)
  implicit none
  
  character(*), intent(in) :: input
  type(String), intent(in) :: colour_name
  type(String)             :: output
  
  output = colour(input,char(colour_name))
end function

function colour_character_String(input,colour_name) result(output)
  implicit none
  
  type(String), intent(in) :: input
  character(*), intent(in) :: colour_name
  type(String)             :: output
  
  output = colour(char(input),colour_name)
end function

impure elemental function colour_String_String(input,colour_name) &
   & result(output)
  implicit none
  
  type(String), intent(in) :: input
  type(String), intent(in) :: colour_name
  type(String)             :: output
  
  output = colour(char(input),char(colour_name))
end function

! ----------------------------------------------------------------------
! Set or unset OUTPUT_FILE_UNIT as a file unit.
! ----------------------------------------------------------------------
subroutine set_output_unit(file_unit)
  implicit none
  
  integer, intent(in) :: file_unit
  
  if (OUTPUT_FILE_UNIT/=OUTPUT_UNIT) then
    call print_line('Code Error: attempted to redirect stdout when it was &
       &already redirected.')
    call err()
  endif
  
  OUTPUT_FILE_UNIT = file_unit
  call set_error_strings_uncoloured()
end subroutine

subroutine unset_output_unit()
  implicit none
  
  OUTPUT_FILE_UNIT = OUTPUT_UNIT
  call set_error_strings_coloured()
end subroutine

subroutine set_terminal_width()
  implicit none
  
  integer :: width
  logical :: success
  
  success = get_terminal_width_c(width)
  if (success) then
    TERMINAL_WIDTH = width
  else
    call print_line( 'Failed to get terminal width. &
                     &Reverting to default of 79 characters.')
    TERMINAL_WIDTH = 79
  endif
end subroutine
end module
