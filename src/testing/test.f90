! ======================================================================
! A test space, for temporary checking of misc. parts of the code.
! ======================================================================
module test_module
  use common_module
  
  implicit none
  
  private
  
  public :: startup_test
contains

! ----------------------------------------------------------------------
! Generates keywords and helptext.
! ----------------------------------------------------------------------
subroutine startup_test()
  implicit none
  
  type(CaesarMode) :: mode
  
  mode%mode_name = 'test'
  mode%description = 'Runs temporary code for testing purposes.'
  mode%keywords = [ &
     & KeywordData('a','a',exclusive_with=[str('b')]),          &
     & KeywordData('b','b',exclusive_with=[str('a')]),          &
     & KeywordData('c','c',exclusive_with=[str('d'),str('e'),str('f')]), &
     & KeywordData('d','d',exclusive_with=[str('e'),str('f'),str('c')]), &
     & KeywordData('e','e',exclusive_with=[str('c'),str('d'),str('f')]), &
     & KeywordData('f','f',exclusive_with=[str('c'),str('d'),str('e')]), &
     & KeywordData('g','g')                                     &
     & ]
  mode%main_subroutine => test_subroutine
  mode%suppress_from_helptext = .true.
  
  call add_mode(mode)
end subroutine

! ----------------------------------------------------------------------
! Main function.
! ----------------------------------------------------------------------
subroutine test_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  if (arguments%is_set('a')) then
    call print_line('a: '//arguments%value('a'))
  elseif (arguments%is_set('b')) then
    call print_line('b: '//arguments%value('b'))
  else
    call print_line('a and b unset.')
  endif
  
  if (arguments%is_set('c')) then
    call print_line('c: '//arguments%value('c'))
  elseif (arguments%is_set('d')) then
    call print_line('d: '//arguments%value('d'))
  elseif (arguments%is_set('e')) then
    call print_line('e: '//arguments%value('e'))
  elseif (arguments%is_set('f')) then
    call print_line('f: '//arguments%value('f'))
  else
    call print_line('c,d,e and f unset.')
  endif
end subroutine
end module
