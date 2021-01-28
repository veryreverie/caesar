! ======================================================================
! Checks SharedCounter to see whether a gfortran bug is present or not.
! ======================================================================
module caesar_check_counter_module
  use caesar_common_module
  implicit none
  
  private
  
  public :: startup_check_counter
contains

! ----------------------------------------------------------------------
! Generates keywords and helptext.
! ----------------------------------------------------------------------
subroutine startup_check_counter()
  implicit none
  
  type(CaesarMode) :: mode
  
  integer :: ialloc
  
  mode%mode_name = 'check_counter'
  mode%description = 'Checks for a gfortran bug with shared counters.'
  allocate(mode%keywords(0), stat=ialloc); call err(ialloc)
  mode%main_subroutine => check_counter_subroutine
  
  call add_mode(mode)
end subroutine

! ----------------------------------------------------------------------
! Main function.
! ----------------------------------------------------------------------
subroutine check_counter_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  type(SharedCounter) :: a
  
  call print_line('Initialising counter A.')
  a = SharedCounter()
  call print_line('Does counter A believe itself to be unique? '// &
                & colour_check(a%is_only_pointer(),.true.))
  call print_line('Passing counter A into subroutine.')
  call check_counter_subroutine_2(a)
  call print_line('Does counter A believe itself to be unique? '// &
                & colour_check(a%is_only_pointer(),.true.))
end subroutine

subroutine check_counter_subroutine_2(a)
  implicit none
  
  type(SharedCounter), intent(inout) :: a
  
  type(SharedCounter) :: b
  
  call print_line('Does counter A believe itself to be unique? '// &
                & colour_check(a%is_only_pointer(),.true.))
  call print_line('Initialising counter B to counter A.')
  b = a
  call print_line('Does counter A believe itself to be unique? '// &
                & colour_check(a%is_only_pointer(),.false.))
  call print_line('Does counter B believe itself to be unique? '// &
                & colour_check(b%is_only_pointer(),.false.))
  call print_line('Exiting subroutine.')
end subroutine

function colour_check(input,expected) result(output)
  implicit none
  
  logical, intent(in) :: input
  logical, intent(in) :: expected
  type(String)        :: output
  
  if (input.eqv.expected) then
    output = colour(str(input),'green')//' (Should be '//expected//')'
  else
    output = colour(str(input),'red')//' (Should be '//expected//')'
  endif
end function
end module
