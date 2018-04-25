! ======================================================================
! An abstract type, which allows extended types to be turned into strings.
! ======================================================================
! An abstract type, which allows extended types to be turned into strings.
! Any type which extends Stringable can be:
!    - converted to String, using string=this or str(this).
!    - concatenated, using string//this or character//this.
!    - printed to stdout, using print_line(this).
!    - printed to file, using file%print_line(this).
! Any type which extends Stringable must overload %to_String(). See the example
!    module below.
module stringable_submodule
  use string_submodule
  use error_submodule
  implicit none
  
  private
  
  public :: Stringable
  
  public :: assignment(=)
  public :: operator(//)
  public :: str
  
  type, abstract :: Stringable
  contains
    procedure(to_String_Stringable), deferred :: to_String
  end type
  
  abstract interface
    recursive function to_String_Stringable(this) result(output)
      import String
      import Stringable
      implicit none
      
      class(Stringable), intent(in) :: this
      type(String)                  :: output
    end function
  end interface
  
  ! String = this.
  interface assignment(=)
    module procedure assign_String_Stringable
  end interface
    
  ! String = char//this or this//char or String//this or this//String.
  interface operator(//)
    module procedure concatenate_Stringable_character
    module procedure concatenate_character_Stringable
    module procedure concatenate_Stringable_String
    module procedure concatenate_String_Stringable
  end interface
  
  ! String = str(this).
  interface str
    module procedure str_Stringable_0d
    module procedure str_Stringable_1d
    module procedure str_Stringable_2d
  end interface
contains

! ----------------------------------------------------------------------
! String = Stringable.
! ----------------------------------------------------------------------
recursive subroutine assign_String_Stringable(output,input)
  implicit none
  
  type(String),      intent(out) :: output
  class(Stringable), intent(in)  :: input
  
  output = input%to_String()
end subroutine

! ----------------------------------------------------------------------
! Concatenation of string types and Stringable types.
! ----------------------------------------------------------------------
! String = Stringable//character
recursive function concatenate_Stringable_character(this,that) result(output)
  implicit none
  
  class(Stringable), intent(in) :: this
  character(*),      intent(in) :: that
  type(String)                  :: output
  
  output = this%to_String()//that
end function

! String = character//Stringable
recursive function concatenate_character_Stringable(this,that) result(output)
  implicit none
  
  character(*),      intent(in) :: this
  class(Stringable), intent(in) :: that
  type(String)                  :: output
  
  output = this//that%to_String()
end function

! String = Stringable//String
recursive function concatenate_Stringable_String(this,that) result(output)
  implicit none
  
  class(Stringable), intent(in) :: this
  type(String),      intent(in) :: that
  type(String)                  :: output
  
  output = this%to_String()//that
end function

! String = String//Stringable
recursive function concatenate_String_Stringable(this,that) result(output)
  implicit none
  
  type(String),      intent(in) :: this
  class(Stringable), intent(in) :: that
  type(String)                  :: output
  
  output = this//that%to_String()
end function

! ----------------------------------------------------------------------
! String = str(Stringable).
! ----------------------------------------------------------------------
! N.B. can't use impure elemental because this must be recursive.
! N.B. can't use err() because that relies on this module.
recursive function str_Stringable_0d(this) result(output)
  implicit none
  
  class(Stringable), intent(in) :: this
  type(String)                  :: output
  
  output = this
end function

recursive function str_Stringable_1d(this) result(output)
  implicit none
  
  class(Stringable), intent(in) :: this(:)
  type(String), allocatable     :: output(:)
  
  integer :: i,ialloc
  
  allocate(output(size(this)), stat=ialloc); call err(ialloc)
  do i=1,size(this)
    output(i) = str(this(i))
  enddo
end function

recursive function str_Stringable_2d(this) result(output)
  implicit none
  
  class(Stringable), intent(in) :: this(:,:)
  type(String), allocatable     :: output(:,:)
  
  integer :: i,ialloc
  
  allocate(output(size(this,1),size(this,2)), stat=ialloc); call err(ialloc)
  do i=1,size(this,2)
    output(:,i) = str(this(:,i))
  enddo
end function
end module

! ======================================================================
! An example module to demonstrate the use of Stringable.
! ======================================================================
module stringable_example_submodule
  use stringable_submodule
  use string_submodule
  implicit none
  
  type, extends(Stringable) :: StringableExample
    character(1) :: contents
  contains
    procedure :: to_String => to_String_StringableExample
  end type
contains

recursive function to_String_StringableExample(this) result(output)
  implicit none
  
  class(StringableExample), intent(in) :: this
  type(String)                         :: output
  
  output = str(this%contents)
end function
end module
