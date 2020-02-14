! ======================================================================
! A test space, for temporary checking of misc. parts of the code.
! ======================================================================
module test__module
  use common_module
  implicit none
  
  type, extends(NoDefaultConstructor), abstract :: Parent
  contains
    procedure(f1_Parent), public, deferred                         :: f1
    procedure(f2_Parent), public, deferred                         :: f2
    !procedure(f3_Parent), public, deferred, pass(this), pass(that) :: f3
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
      
      class(Parent), intent(in) :: this
      class(Parent), intent(in) :: that
    end subroutine
    
    !subroutine f3_Parent(this,that)
    !  import Parent
    !  implicit none
    !  
    !  class(Parent), intent(in) :: this
    !  class(Parent), intent(in) :: that
    !end subroutine
  end interface
  
  type, extends(Parent) :: ParentPointer
    class(parent), allocatable :: parent_
  contains
    procedure, public                  :: f1 => f1_Pointer
    procedure, public                  :: f2 => f2_Pointer
    !procedure, public, pass(this,that) :: f3 => f3_Pointer
  end type
  
  type, extends(Parent) :: Child
    real(dp), allocatable :: big_array(:)
  contains
    procedure, public                  :: f1 => f1_Child
    procedure, public                  :: f2 => f2_Child
    !procedure, public, pass(this,that) :: f3 => f3_Child
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
  class(Parent),         intent(in) :: that
  
  call this%parent_%f2(that)
end subroutine

!subroutine f3_Pointer(this,that)
!  implicit none
!  
!  class(ParentPointer), intent(in) :: this
!  class(ParentPointer), intent(in) :: that
!  
!  ! TODO
!end subroutine

subroutine f1_Child(this)
  implicit none
  
  class(Child), intent(in) :: this
  
  call print_line(this%big_array(1))
end subroutine

subroutine f2_Child(this,that)
  implicit none
  
  class(Child),  intent(in) :: this
  class(Parent), intent(in) :: that
  
  type(Child), pointer :: child_
  
  real(dp) :: thing
  
  child_ => Child(that)
  thing = dot_product(this%big_array, child_%big_array)
  call print_line(thing)
end subroutine

!subroutine f3_Child(this,that)
!  implicit none
!  
!  class(Child), intent(in) :: this
!  class(Child), intent(in) :: that
!  
!  ! TODO
!end subroutine

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

module a_module
  type :: T1
  end type
  
  type :: T2
    type(T1), pointer :: a
    type(T1), pointer :: b
  end type
end module

module test_module
  use common_module
  
  use anharmonic_module
  
  use test__module
  
  use a_module
  
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
  
  type(Child) :: a
  
  class(Parent), allocatable :: b(:)
  
  integer :: i
  
  type(T1), target :: c
  type(T2) :: d
  
  a = Child()
  
  b = [(Child(), i=1, 100)]
  
  d%a=>c
  call print_line(associated(d%a))
  call print_line(associated(d%b))
  call x(d)
end subroutine

subroutine x(this)
  implicit none
  
  type(T2), intent(in) :: this
  
  call print_line(associated(this%a))
  call print_line(associated(this%b))
end subroutine
end module
