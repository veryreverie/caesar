submodule (caesar_strings_writeable_module) caesar_strings_writeable_submodule
  use caesar_io_module
contains
module procedure str_StringsWriteable
  if (present(settings)) then
    call set_print_settings(settings)
  endif
  
  output = this%write()
  
  if (present(settings)) then
    call unset_print_settings()
  endif
end procedure

module procedure str_StringsWriteables_String
  type(StringArray), allocatable :: sections(:)
  
  integer :: i,ialloc
  
  if (present(settings)) then
    call set_print_settings(settings)
  endif
  
  allocate(sections(size(this)), stat=ialloc); call err(ialloc)
  do i=1,size(this)
    sections(i) = StringArray(str(this(i)))
  enddo
  
  output = str(join(sections,separating_line))
  
  if (present(settings)) then
    call unset_print_settings()
  endif
end procedure

module procedure str_StringsWriteables_character
  output = str(this,str(separating_line))
end procedure

module procedure print_lines_StringsWriteable
  if (present(settings)) then
    call set_print_settings(settings)
  endif
  
  call print_lines(str(this))
  
  if (present(settings)) then
    call unset_print_settings()
  endif
end procedure

module procedure print_lines_StringsWriteables_String
  if (present(settings)) then
    call set_print_settings(settings)
  endif
  
  call print_lines(str(this,separating_line))
  
  if (present(settings)) then
    call unset_print_settings()
  endif
end procedure

module procedure print_lines_StringsWriteables_character
  if (present(settings)) then
    call set_print_settings(settings)
  endif
  
  call print_lines(str(this,separating_line))
  
  if (present(settings)) then
    call unset_print_settings()
  endif
end procedure
end submodule
