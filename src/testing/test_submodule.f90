submodule (caesar_test_module) caesar_test_submodule
  use caesar_testing_module
contains

module procedure startup_test
  type(CaesarMode) :: mode
  
  integer :: ialloc
  
  mode%mode_name = 'test'
  mode%description = 'Runs temporary code for testing purposes.'
  allocate(mode%keywords(0), stat=ialloc); call err(ialloc)
  mode%main_subroutine => test_subroutine
  mode%suppress_from_helptext = .true.
  
  call add_mode(mode)
end procedure

module procedure test_subroutine
end procedure
end submodule
