! ======================================================================
! A set of modes which are coupled to one another.
! Does not specify whether the modes are real or complex.
! ======================================================================
module mode_coupling_module
  use common_module
  implicit none
  
  private
  
  public :: ModeCoupling
  public :: size
  public :: operator(==)
  public :: operator(/=)
  
  type, extends(Stringable) :: ModeCoupling
    ! The ids of the modes which are coupled together.
    ! ids is assumed to always be sorted in ascending order.
    integer, allocatable :: ids(:)
  contains
    ! I/O.
    procedure, public :: read  => read_ModeCoupling
    procedure, public :: write => write_ModeCoupling
  end type
  
  interface ModeCoupling
    module procedure new_ModeCoupling
    module procedure new_ModeCoupling_RealMonomial
    module procedure new_ModeCoupling_ComplexMonomial
    module procedure new_ModeCoupling_String
  end interface
  
  interface size
    module procedure size_ModeCoupling
  end interface
  
  interface operator(==)
    module procedure equality_ModeCoupling_ModeCoupling
  end interface
  
  interface operator(/=)
    module procedure non_equality_ModeCoupling_ModeCoupling
  end interface
contains

! ----------------------------------------------------------------------
! Constructors.
! ----------------------------------------------------------------------
function new_ModeCoupling(ids) result(this)
  implicit none
  
  integer, intent(in) :: ids(:)
  type(ModeCoupling)  :: this
  
  this%ids = ids
end function

impure elemental function new_ModeCoupling_RealMonomial(input) result(this)
  implicit none
  
  type(RealMonomial), intent(in) :: input
  type(ModeCoupling)             :: this
  
  integer, allocatable :: ids(:)
  
  ids = input%modes%id
  ids = sort(ids)
  this = ModeCoupling(ids)
end function

impure elemental function new_ModeCoupling_ComplexMonomial(input) result(this)
  implicit none
  
  type(ComplexMonomial), intent(in) :: input
  type(ModeCoupling)                :: this
  
  integer, allocatable :: ids(:)
  
  ids = input%modes%id
  ids = sort(ids)
  this = ModeCoupling(ids)
end function

! ----------------------------------------------------------------------
! Basis type functionality: size() and comparison operators.
! ----------------------------------------------------------------------
function size_ModeCoupling(this) result(output)
  implicit none
  
  type(ModeCoupling), intent(in) :: this
  integer                        :: output
  
  output = size(this%ids)
end function

impure elemental function equality_ModeCoupling_ModeCoupling(this,that) &
   & result(output)
  implicit none
  
  type(ModeCoupling), intent(in) :: this
  type(ModeCoupling), intent(in) :: that
  logical                        :: output
  
  if (size(this)/=size(that)) then
    output = .false.
  else
    if (all(this%ids==that%ids)) then
      output = .true.
    else
      output = .false.
    endif
  endif
end function

impure elemental function non_equality_ModeCoupling_ModeCoupling(this,that) &
   & result(output)
  implicit none
  
  type(ModeCoupling), intent(in) :: this
  type(ModeCoupling), intent(in) :: that
  logical                        :: output
  
  output = .not. this==that
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_ModeCoupling(this,input)
  implicit none
  
  class(ModeCoupling), intent(out) :: this
  type(String),        intent(in)  :: input
  
  select type(this); type is(ModeCoupling)
    this = ModeCoupling(int(split_line(input)))
  end select
end subroutine

function write_ModeCoupling(this) result(output)
  implicit none
  
  class(ModeCoupling), intent(in) :: this
  type(String)                    :: output
  
  select type(this); type is(ModeCoupling)
    output = join(this%ids)
  end select
end function

impure elemental function new_ModeCoupling_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(ModeCoupling)       :: this
  
  this = input
end function
end module
