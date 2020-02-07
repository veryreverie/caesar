! ======================================================================
! A test space, for temporary checking of misc. parts of the code.
! ======================================================================
module test__module
  use common_module
  implicit none
  
  type, extends(NoDefaultConstructor), abstract :: Parent
  contains
    procedure(f1_Parent), public, deferred :: f1
    procedure(f2_Parent), public, deferred :: f2
  end type
  
  abstract interface
    subroutine f1_Parent(this)
      import Parent
      implicit none
      
      class(Parent), intent(in) :: this
    end subroutine
    
    subroutine f2_Parent(this,that)
      import Parent
      implicit none
      
      class(Parent),         intent(in) :: this
      class(Parent), target, intent(in) :: that
    end subroutine
  end interface
  
  type, extends(Parent) :: ParentPointer
    class(parent), allocatable :: parent_
  contains
    procedure, public :: f1 => f1_Pointer
    procedure, public :: f2 => f2_Pointer
  end type
  
  type, extends(Parent) :: Child
    real(dp), allocatable :: big_array(:)
  contains
    procedure, public :: f1 => f1_Child
    procedure, public :: f2 => f2_Child
  end type
  
  interface Child
    module procedure new_Child
    module procedure new_Child_Parent
  end interface
  
  type :: Holder
    type(Child) :: c1
    type(Child) :: c2
  contains
    procedure, public :: f => f_Holder
  end type
contains

subroutine f1_Pointer(this)
  implicit none
  
  class(ParentPointer), intent(in) :: this
  
  call this%parent_%f1()
end subroutine

subroutine f2_Pointer(this,that)
  implicit none
  
  class(ParentPointer),  intent(in) :: this
  class(Parent), target, intent(in) :: that
  
  call this%parent_%f2(that)
end subroutine

subroutine f1_Child(this)
  implicit none
  
  class(Child), intent(in) :: this
  
  call print_line(this%big_array(1))
end subroutine

subroutine f2_Child(this,that)
  implicit none
  
  class(Child),          intent(in) :: this
  class(Parent), target, intent(in) :: that
  
  type(Child), pointer :: child_
  
  real(dp) :: thing
  
  child_ => Child(that)
  thing = dot_product(this%big_array, child_%big_array)
  call print_line(thing)
end subroutine

function new_Child() result(this)
  implicit none
  
  type(Child) :: this
  
  integer :: i
  
  allocate(this%big_array(1000000))
  do i=1,size(this%big_array)
    call random_number(this%big_array(i))
  enddo
end function

recursive function new_Child_Parent(input) result(this)
  implicit none
  
  class(Parent), target, intent(in) :: input
  type(Child),   pointer            :: this
  
  select type(input); type is(ParentPointer)
    this => new_Child_Parent(input%parent_)
  type is(Child)
    this => input
  end select
end function

subroutine f_Holder(this)
  implicit none
  
  class(Holder), target, intent(in) :: this
  
  call this%c1%f2(this%c2)
end subroutine
end module

module test_module
  use common_module
  
  use anharmonic_module
  
  use test__module
  
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
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  type(Child), target :: child_
  
  type(Holder) :: holder_
  
  integer :: i
  
  class(Parent), pointer :: p_
  
  child_ = Child()
  
  do i=1,1000
    call child_%f2(child_)
  enddo
  
  holder_ = Holder(Child(),Child())
  do i=1,1000
    call holder_%f()
  enddo
  
  p_ => child_
  select type(p_); type is(Child)
    call g(p_)
  end select
end subroutine

subroutine g(input)
  class(Child), intent(in) :: input
end subroutine
end module
