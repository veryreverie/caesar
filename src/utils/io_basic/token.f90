! ======================================================================
! Provides the token function, which splits a string by a given delimiter
!    and returns the token at the requested index, and the tokens function,
!    which does the same but returning an array of tokens rather than a
!    single token, as specified by a start index and end index.
! ======================================================================
module token_module
  use error_module
  use string_module
  use print_module
  use intrinsics_module
  implicit none
  
  private
  
  public :: token
  public :: tokens
  
  interface token
    module procedure token_character
    module procedure token_String
  end interface
  
  interface tokens
    module procedure tokens_character
    module procedure tokens_String
  end interface
contains

! Split a string, and return the token at the given index.
function token_String(input,index,delimiter,delimiters) result(output)
  implicit none
  
  type(String), intent(in)           :: input
  integer,      intent(in)           :: index
  character(1), intent(in), optional :: delimiter
  character(1), intent(in), optional :: delimiters(:)
  type(String)                       :: output
  
  type(String), allocatable :: line(:)
  
  type(String) :: delimiters_string
  
  integer :: i
  
  if (index<1) then
    call print_line(ERROR//': index<1.')
    call err()
  endif
  
  line = split_line(input, delimiter, delimiters)
  
  if (index>size(line)) then
    call print_line(ERROR//': String does not contain enough tokens.')
    call print_line('Input string    : '//input)
    call print_line('Requested index : '//index)
    if (present(delimiter)) then
      if (delimiter=='"') then
        call print_line("Delimiter       : '"//delimiter//"'")
      else
        call print_line('Delimiter       : "'//delimiter//'"')
      endif
    elseif (present(delimiters)) then
      if (size(delimiters)==0) then
        call print_line('No Delimiters given.')
      else
        delimiters_string = str("")
        do i=1,size(delimiters)
          if (delimiters(i)=='"') then
            delimiters_string = delimiters_string//"'"//delimiters(i)//"'"
          else
            delimiters_string = delimiters_string//'"'//delimiters(i)//'"'
          endif
          if (i<size(delimiters)) then
            delimiters_string = delimiters_string//', '
          endif
        enddo
        call print_line('Delimiters      : '//delimiters_string)
      endif
    endif
    call err()
  endif
  
  output = line(index)
end function

function token_character(input,index,delimiter,delimiters) result(output)
  implicit none
  
  character(*), intent(in)           :: input
  integer,      intent(in)           :: index
  character(1), intent(in), optional :: delimiter
  character(1), intent(in), optional :: delimiters(:)
  type(String)                       :: output
  
  output = token(str(input), index, delimiter)
end function

! Split a string, and return the tokens in the given range of indices.
! first defaults to 1.
! last defaults to the number of tokens.
function tokens_String(input,first,last,delimiter,delimiters) result(output)
  implicit none
  
  type(String), intent(in)           :: input
  integer,      intent(in), optional :: first
  integer,      intent(in), optional :: last
  character(1), intent(in), optional :: delimiter
  character(1), intent(in), optional :: delimiters(:)
  type(String), allocatable          :: output(:)
  
  type(String), allocatable :: line(:)
  
  line = split_line(input, delimiter, delimiters)
  
  ! If first is present, check that 1 <= first <= len(line).
  if (present(first)) then
    if (first<1) then
      call print_line(ERROR//': first<1.')
      call err()
    endif
    
    if (first>size(line)) then
      call print_line(ERROR//': String does not contain enough tokens.')
      call print_line('Input string    : '//input)
      call print_line('Requested index : '//first)
      if (present(delimiter)) then
        call print_line('Delimiter       : "'//delimiter//'"')
      endif
      call err()
    endif
  endif
    
  ! If last is present, check that 1 <= last <= len(line).
  if (present(last)) then
    if (last<1) then
      call print_line(ERROR//': first<1.')
      call err()
    endif
    
    if (last>size(line)) then
      call print_line(ERROR//': String does not contain enough tokens.')
      call print_line('Input string    : '//input)
      call print_line('Requested index : '//last)
      if (present(delimiter)) then
        call print_line('Delimiter       : "'//delimiter//'"')
      endif
      call err()
    endif
  endif
  
  ! Return the requested slice.
  if (present(first) .and. present(last)) then
    if (last<first) then
      call print_line(ERROR//': last<first.')
    endif
    output = line(first:last)
  elseif (present(first)) then
    output = line(first:)
  elseif (present(last)) then
    output = line(:last)
  else
    output = line
  endif
end function

function tokens_character(input,first,last,delimiter,delimiters) result(output)
  implicit none
  
  character(*), intent(in)           :: input
  integer,      intent(in), optional :: first
  integer,      intent(in), optional :: last
  character(1), intent(in), optional :: delimiter
  character(1), intent(in), optional :: delimiters(:)
  type(String), allocatable          :: output(:)
  
  output = tokens(str(input), first, last, delimiter, delimiters)
end function
end module
