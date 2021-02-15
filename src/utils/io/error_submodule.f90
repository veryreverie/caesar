submodule (caesar_error_module) caesar_error_submodule
contains
module procedure quit
  call quit_implementation()
end procedure

module procedure err_none
  write(*,'(a)') ''
  write(*,'(a)') LIGHT_MAGENTA_ESC//'Intentionally aborting with &
                 &stacktrace.'//RESET_ESC
  write(*,'(a)') ''
  write(*,'(a)') LIGHT_MAGENTA_ESC//'compiler_specific.f90 and error.f90 &
                 &can be ignored in stacktrace.'//RESET_ESC
  write(*,'(a)') ''
  
  call err_implementation()
end procedure

module procedure err_allocate_flag
  if (this/=0) then
    write(*,'(a)') ''
    write(*,'(a)') ERROR//': Allocation error.'
    call err()
  endif
end procedure

module procedure set_error_strings_coloured
  ERROR = RED_ESC//'Error'//RESET_ESC
  CODE_ERROR = RED_ESC//'Code Error'//RESET_ESC
  WARNING = LIGHT_MAGENTA_ESC//'Warning'//RESET_ESC
end procedure

module procedure set_error_strings_uncoloured
  ERROR = 'Error'
  CODE_ERROR = 'Code Error'
  WARNING = 'Warning'
end procedure
end submodule
