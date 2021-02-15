submodule (caesar_string_array_module) caesar_string_array_submodule
  implicit none
contains
module procedure new_StringArray_Strings
  this%strings = input
end procedure

module procedure size_StringArray
  output = size(this%strings)
end procedure

module procedure str_StringArray
  output = this%strings
end procedure

module procedure str_StringArrays_String
  output = str(join(this, separating_line=separating_line))
end procedure

module procedure str_StringArrays_character
  output = str(join(this, separating_line=str(separating_line)))
end procedure

module procedure split_into_sections_Strings_String
  type(String) :: separating_line_string
  
  integer :: no_sections
  logical :: reading_section
  
  integer, allocatable :: first_lines(:)
  integer, allocatable :: last_lines(:)
  
  integer :: i,ialloc
  
  if (present(separating_line)) then
    separating_line_string = separating_line
  else
    separating_line_string = ''
  endif
  
  allocate( first_lines(size(this)), &
          & last_lines(size(this)),  &
          & stat=ialloc); call err(ialloc)
  no_sections = 0
  reading_section = .false.
  do i=1,size(this)
    if (this(i)==separating_line_string) then
      ! This line is a separating_line string.
      ! If reading a section, then the end of that section is the line above.
      if (reading_section) then
        last_lines(no_sections) = i-1
        reading_section = .false.
      endif
    else
      ! This line is not a separating_line string.
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
    output(i) = StringArray(this(first_lines(i):last_lines(i)))
  enddo
end procedure

module procedure split_into_sections_Strings_character
  output = split_into_sections(this,str(separating_line))
end procedure

module procedure split_into_sections_StringArray_String
  output = split_into_sections(this%strings,separating_line)
end procedure

module procedure split_into_sections_StringArray_character
  output = split_into_sections(this%strings,separating_line)
end procedure

module procedure concatenate_StringArray_StringArray
  output = StringArray([this%strings, that%strings])
end procedure

module procedure concatenate_StringArray_String
  output = StringArray([this%strings, that])
end procedure

module procedure concatenate_StringArray_character
  output = StringArray([this%strings, str(that)])
end procedure

module procedure concatenate_StringArray_Strings
  output = StringArray([this%strings, that])
end procedure

module procedure concatenate_String_StringArray
  output = StringArray([this, that%strings])
end procedure

module procedure concatenate_character_StringArray
  output = StringArray([str(this), that%strings])
end procedure

module procedure concatenate_Strings_StringArray
  output = StringArray([this, that%strings])
end procedure

module procedure join_StringArrays_String
  type(String), allocatable :: strings(:)
  
  integer :: i,ialloc
  
  allocate(strings(0), stat=ialloc); call err(ialloc)
  do i=1,size(input)
    strings = [strings, input(i)%strings]
    if (i/=size(input) .and. present(separating_line)) then
      strings = [strings, separating_line]
    endif
  enddo
  output = StringArray(strings)
end procedure

module procedure join_StringArrays_character
  output = join(input,str(separating_line))
end procedure
end submodule
