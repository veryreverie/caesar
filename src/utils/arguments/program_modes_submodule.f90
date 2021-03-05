submodule (caesar_program_modes_module) caesar_program_modes_submodule
contains

module procedure new_ProgramModes_ProgramMode
  this%modes_ = modes
end procedure

module procedure print_help_ProgramModes
  integer :: i
  
  do i=1,size(this%modes_)
    call this%modes_(i)%print_help()
  enddo
end procedure

module procedure mode_String
  type(String) :: lower_case_mode_name
  
  integer :: i
  
  lower_case_mode_name = lower_case(mode_name)
  
  i = first(this%modes_%mode_name==lower_case_mode_name, default=0)
  
  if (i==0) then
    call print_line(ERROR//': Unrecognised mode: '//mode_name)
    call print_line('Call '//colour('caesar -h','white')//' for help.')
    call quit()
  endif
  
  output = this%modes_(i)
end procedure

module procedure mode_character
  output = this%mode(str(mode_name))
end procedure
end submodule
