submodule (caesar_ofile_module) caesar_ofile_submodule
  use caesar_io_module
contains
module procedure new_OFile_character
  integer :: ialloc
  
  type(OFileTarget) :: ofile_target
  
  ofile_target = OFileTarget(filename)
  
  allocate(this%ofile_target, stat=ialloc); call err(ialloc)
  this%ofile_target = ofile_target
  
  allocate(this%counter, stat=ialloc); call err(ialloc)
  this%counter = SharedCounter()
end procedure

module procedure new_OFile_String
  this = OFile(char(filename))
end procedure

module procedure assign_OFile_OFile
  integer :: ialloc
  
  output%ofile_target => input%ofile_target
  
  allocate(output%counter, stat=ialloc); call err(ialloc)
  output%counter = input%counter
end procedure

module procedure final_OFile
  integer :: ialloc
  
  if (allocated(this%counter)) then
    if (this%counter%is_only_copy()) then
      call this%ofile_target%close()
      deallocate(this%ofile_target, stat=ialloc); call err(ialloc)
    endif
    deallocate(this%counter, stat=ialloc); call err(ialloc)
  endif
end procedure

module procedure check_associated
  if (.not. associated(this%ofile_target)) then
    call print_line(CODE_ERROR//': Attempting to call output file operations &
       &on OFile which has not been associated.')
    call err()
  endif
end procedure

module procedure make_stdout
  call this%check_associated()
  call this%ofile_target%make_stdout()
end procedure

module procedure print_line_character
  call this%check_associated()
  call this%ofile_target%print_line(input,settings)
end procedure

module procedure print_line_String
  call this%print_line(char(input),settings)
end procedure

module procedure print_line_StringWriteable
  call this%print_line(str(input,settings))
end procedure

module procedure print_line_logical
  call this%print_line(str(input,settings))
end procedure

module procedure print_line_integer
  call this%print_line(str(input,settings))
end procedure

module procedure print_line_real
  call this%print_line(str(input,settings))
end procedure

module procedure print_line_complex
  call this%print_line(str(input, settings=settings))
end procedure

module procedure print_line_logicals
  call this%print_line(join(input, settings=settings))
end procedure

module procedure print_line_integers
  call this%print_line(join(input, settings=settings))
end procedure

module procedure print_line_reals
  call this%print_line(join(input, settings=settings))
end procedure

module procedure print_line_complexes
  call this%print_line(join(input, settings=settings))
end procedure

module procedure print_lines_Strings_String
  integer :: i
  
  do i=1,size(input)
    call this%print_line(input(i),settings)
    if (present(separating_line) .and. i<size(input)) then
      call print_line(separating_line,settings)
    endif
  enddo
end procedure

module procedure print_lines_Strings_character
  call this%print_lines(input,str(separating_line),settings)
end procedure

module procedure print_lines_StringWriteables_String
  call this%print_lines(str(input,separating_line,settings))
end procedure

module procedure print_lines_StringWriteables_character
  call this%print_lines(str(input,separating_line,settings))
end procedure

module procedure print_lines_StringsWriteable
  call this%print_lines(str(input,settings))
end procedure

module procedure print_lines_StringsWriteables_String
  call this%print_lines(str(input,separating_line,settings))
end procedure

module procedure print_lines_StringsWriteables_character
  if (lbound(input,1)/=1) then
    call print_line(CODE_ERROR//': The lower bound of an array is not 1. This &
       &is probably a compiler bug.')
    call err()
  endif
  
  call this%print_lines(str(input,separating_line,settings))
end procedure

module procedure print_lines_logicals_String
  call this%print_lines(str(input,separating_line,settings))
end procedure

module procedure print_lines_logicals_character
  call this%print_lines(str(input,separating_line,settings))
end procedure

module procedure print_lines_integers_String
  call this%print_lines(str(input,separating_line,settings))
end procedure

module procedure print_lines_integers_character
  call this%print_lines(str(input,separating_line,settings))
end procedure

module procedure print_lines_reals_String
  call this%print_lines(str(input,separating_line,settings))
end procedure

module procedure print_lines_reals_character
  call this%print_lines(str(input,separating_line,settings))
end procedure

module procedure print_lines_complexes_String
  call this%print_lines(str(input,separating_line,settings))
end procedure

module procedure print_lines_complexes_character
  call this%print_lines(str(input,separating_line,settings))
end procedure
end submodule
