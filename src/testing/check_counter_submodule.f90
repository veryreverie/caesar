submodule (caesar_check_counter_module) caesar_check_counter_submodule
  use caesar_testing_module
contains

module procedure startup_check_counter
  type(CaesarMode) :: mode
  
  integer :: ialloc
  
  mode%mode_name = 'check_counter'
  mode%description = 'Checks for a gfortran bug with shared counters.'
  allocate(mode%keywords(0), stat=ialloc); call err(ialloc)
  mode%main_subroutine => check_counter_subroutine
  
  call add_mode(mode)
end procedure

module procedure check_counter_subroutine
  type(SharedCounter) :: a
  
  call print_line('Initialising counter A.')
  a = SharedCounter()
  call print_line('Does counter A believe itself to be unique? '// &
                & colour_check(a%is_only_copy(),.true.))
  call print_line('Passing counter A into module procedure.')
  call check_counter_subroutine_2(a)
  call print_line('Does counter A believe itself to be unique? '// &
                & colour_check(a%is_only_copy(),.true.))
end procedure

module procedure check_counter_subroutine_2
  type(SharedCounter) :: b
  
  call print_line('Does counter A believe itself to be unique? '// &
                & colour_check(a%is_only_copy(),.true.))
  call print_line('Initialising counter B to counter A.')
  b = a
  call print_line('Does counter A believe itself to be unique? '// &
                & colour_check(a%is_only_copy(),.false.))
  call print_line('Does counter B believe itself to be unique? '// &
                & colour_check(b%is_only_copy(),.false.))
  call print_line('Exiting module procedure.')
end procedure

module procedure colour_check
  if (input.eqv.expected) then
    output = colour(str(input),'green')//' (Should be '//expected//')'
  else
    output = colour(str(input),'red')//' (Should be '//expected//')'
  endif
end procedure
end submodule
