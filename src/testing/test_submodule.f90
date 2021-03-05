submodule (caesar_test_module) caesar_test_submodule
  use caesar_testing_module
contains

module procedure test_mode
  integer :: ialloc
  
  output%mode_name = 'test'
  output%description = 'Runs temporary code for testing purposes.'
  allocate(output%keywords(0), stat=ialloc); call err(ialloc)
  output%main_subroutine => test_subroutine
  output%suppress_from_helptext = .true.
end procedure

module procedure test_subroutine
end procedure
end submodule
