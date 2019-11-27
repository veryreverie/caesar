! ======================================================================
! A test space, for temporary checking of misc. parts of the code.
! ======================================================================
module test_module
  use common_module
  
  use anharmonic_module
  
  implicit none
  
  private
  
  public :: startup_test
  
  type :: T1
    integer :: contents
  contains
    procedure, public :: print => print_T1
  end type
  
  type, extends(T1) :: T2
    type(T1) :: t_1
  contains
    procedure, public :: print => print_T2
  end type
  
  interface T1
    module procedure new_T1
    module procedure new_T1_T2
  end interface
  
  interface T2
    module procedure new_T2
  end interface
contains

function new_T1(contents) result(this)
  implicit none
  
  integer, intent(in) :: contents
  type(T1)            :: this
  
  this%contents = contents
end function

function new_T2(t_1) result(this)
  implicit none
  
  type(T1), intent(in) :: t_1
  type(T2)             :: this
  
  this%t_1 = t_1
  this%contents = t_1%contents
end function

function new_T1_T2(this) result(output)
  implicit none
  
  class(T2), intent(in) :: this
  type(T1), pointer     :: output
  
  output = this%t_1
end function

subroutine print_T1(this)
  implicit none
  
  class(T1), intent(in) :: this
  
  call print_line('T1: '//this%contents)
end subroutine

subroutine print_T2(this)
  implicit none
  
  class(T2), intent(in) :: this
  
  call print_line('T2: '//this%contents//' = '//this%t_1%contents)
end subroutine

! ----------------------------------------------------------------------
! Generates keywords and helptext.
! ----------------------------------------------------------------------
subroutine startup_test()
  implicit none
  
  type(CaesarMode) :: mode
  
  integer :: ialloc
  
  mode%mode_name = 'test'
  mode%description = 'Runs temporary code for testing purposes.'
  allocate(mode%keywords(0), stat=ialloc); call err(ialloc)
  mode%main_subroutine => test_subroutine
  mode%suppress_from_helptext = .true.
  
  call add_mode(mode)
end subroutine

! ----------------------------------------------------------------------
! Main function.
! ----------------------------------------------------------------------
subroutine test_subroutine(arguments)
  use interpolation_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  call heaps_algorithm([1,2,3,3])
  
  call print_line('')
  
  call no_repeats([1,2,1,2,1,2])
  call print_line('')
  call no_repeats([1,2,3,2,1,2])
  call print_line('')
  call no_repeats([1,2,3,4])
  call print_line('')
  call no_repeats([1,2,3,4,5,6])
  call print_line('')
  call no_repeats([1,2,3,4,5,6,7,8])
  call print_line('')
  call no_repeats([1,2,3,4,5,6,7,8,9,10])
  
  call print_line('')
  call print_line('')
  call perm([1,2,3],[1,2,3])
  call print_line('')
  call print_line('')
  call perm([1,2,3],[1,2,2])
  call print_line('')
  call print_line('')
  call perm([1,2,3,4],[1,2,2,1])
end subroutine
end module
