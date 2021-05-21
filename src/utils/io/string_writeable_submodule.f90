submodule (caesar_string_writeable_module) caesar_string_writeable_submodule
  use caesar_io_module
contains
module procedure str_StringWriteable
  if (present(settings)) then
    call set_print_settings(settings)
  endif
  
  output = this%write()
  
  if (present(settings)) then
    call unset_print_settings()
  endif
end procedure

module procedure str_StringWriteables_String
  integer :: i,ialloc
  
  if (present(settings)) then
    call set_print_settings(settings)
  endif
  
  if (present(separating_line)) then
    allocate(output(max(2*size(this)-1,0)), stat=ialloc); call err(ialloc)
    do i=1,size(this)
      output(2*i-1) = str(this(i))
      if (i<size(this)) then
        output(2*i) = separating_line
      endif
    enddo
  else
    allocate(output(size(this)), stat=ialloc); call err(ialloc)
    do i=1,size(this)
      output(i) = str(this(i))
    enddo
  endif
  
  if (present(settings)) then
    call unset_print_settings()
  endif
end procedure

module procedure str_StringWriteables_character
  output = str(this, str(separating_line), settings)
end procedure

module procedure concatenate_StringWriteable_character
  output = str(this)//that
end procedure

module procedure concatenate_character_StringWriteable
  output = this//str(that)
end procedure

module procedure concatenate_StringWriteable_String
  output = str(this)//that
end procedure

module procedure concatenate_String_StringWriteable
  output = this//str(that)
end procedure

module procedure join_StringWriteable
  if (present(settings)) then
    call set_print_settings(settings)
  endif
  
  output = join(str(this), delimiter)
  
  if (present(settings)) then
    call unset_print_settings()
  endif
end procedure

module procedure print_line_StringWriteable
  if (present(settings)) then
    call set_print_settings(settings)
  endif
  
  call print_line(str(this))
  
  if (present(settings)) then
    call unset_print_settings()
  endif
end procedure

module procedure print_lines_StringWriteables_String
  if (present(settings)) then
    call set_print_settings(settings)
  endif
  
  call print_lines(str(this,separating_line))
  
  if (present(settings)) then
    call unset_print_settings()
  endif
end procedure

module procedure print_lines_StringWriteables_character
  if (present(settings)) then
    call set_print_settings(settings)
  endif
  
  call print_lines(str(this,separating_line))
  
  if (present(settings)) then
    call unset_print_settings()
  endif
end procedure
end submodule
