submodule (caesar_string_module) caesar_string_submodule
contains

module procedure str_character
  output = this
end procedure

module procedure len_String
  output = len(char(this))
end procedure

module procedure join_String
  type(String) :: delimiter_character
  integer      :: i
  
  if (present(delimiter)) then
    delimiter_character = delimiter
  else
    delimiter_character = ' '
  endif
  
  if (size(this)==0) then
    output = ''
  else
    output = this(1)
    do i=2,size(this)
      output = output//delimiter_character//this(i)
    enddo
  endif
end procedure

module procedure concatenate_String_character_
  output = char(this)//that
end procedure

module procedure concatenate_character_String_
  output = this//char(that)
end procedure

module procedure concatenate_String_String_
  output = char(this)//char(that)
end procedure

module procedure lower_case_character
  character(*), parameter :: lower_chars = "abcdefghijklmnopqrstuvwxyz"
  character(*), parameter :: upper_chars = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
  
  integer :: i,j
  
  output = input
  
  do i=1,len(output)
    j = index(upper_chars, output(i:i))
    if (j/=0) then
      output(i:i) = lower_chars(j:j)
    endif
  enddo
end procedure

module procedure lower_case_String
  output = lower_case(char(this))
end procedure

module procedure spaces
  integer :: i,ialloc
  
  allocate(character(no_spaces) :: output, stat=ialloc); call err(ialloc)
  
  do i=1,no_spaces
    output(i:i) = ' '
  enddo
end procedure

module procedure split_line_character
  ! The position of a delimiter.
  integer :: first
  ! The position of the next delimiter after 'first'.
  integer :: second
  ! Temporary variables.
  integer :: ialloc
  
  first = 0
  second = 0
  allocate(output(0), stat=ialloc); call err(ialloc)
  do
    ! Search after previously found delimiter.
    first = second
    ! Exit if entire word split.
    if (first == len(input)+1) then
      exit
    endif
    ! Find the next delimiter.
    if (present(delimiter)) then
      second = first + index(input(first+1:), delimiter)
    elseif (present(delimiters)) then
      second = first + first_index(input(first+1:), delimiters)
    else
      ! If delimiter is not set, find the next space or tab character.
      ! N.B. the delimiter array is [space, tab].
      second = first + first_index(input(first+1:), [' ', '	'])
    endif
    ! If second==first there is no next delimiter. Parse the final token.
    if (second == first) then
      second = len(input)+1
    endif
    ! If second==first+1, there are multiple delimiters in a row.
    ! They are treated as a single delimiter.
    if (second == first+1) then
      cycle
    endif
    ! Append the token to the output.
    output = [output, str(input(first+1:second-1))]
  enddo
end procedure

module procedure split_line_String
  output = split_line(char(input),delimiter,delimiters)
end procedure

module procedure first_index
  integer :: i,j
  
  output = 0
  
  do i=1,size(delimiters)
    j = index(input,delimiters(i))
    if (j/=0) then
      if (output==0) then
        output = j
      else
        output = min(output,j)
      endif
    endif
  enddo
end procedure

module procedure slice_character
  output = input(first:last)
end procedure

module procedure slice_String
  output = slice(char(input),first,last)
end procedure

module procedure trim_String
  output = trim(adjustl(char(input)))
end procedure

module procedure replace_character
  integer :: i
  
  output = ''
  do i=1,len(input)
    if (input(i:i)==from) then
      output = output//to
    else
      output = output//input(i:i)
    endif
  enddo
end procedure

module procedure replace_String
  output = replace(char(input),from,to)
end procedure
end submodule
