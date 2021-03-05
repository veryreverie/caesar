submodule (caesar_c_string_module) caesar_c_string_submodule
  use caesar_print_module
  use caesar_error_module
contains

module procedure parse_c_string
  integer :: i
  
  output = ''
  do i=1,len(input)
    if (input(i:i)==char(0)) then
      output = input(:i-1)
      exit
    endif
  enddo
  
  if (output=='') then
    call print_line(ERROR//': C string string not nul-terminated: '// &
       & input)
    call err()
  endif
end procedure
end submodule
