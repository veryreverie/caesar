! ======================================================================
! A set of modes which are coupled to one another.
! Does not specify whether the modes are real or complex.
! ======================================================================
module coupled_modes_module
  use common_module
  implicit none
  
  private
  
  public :: CoupledModes
  public :: size
  public :: operator(==)
  public :: operator(/=)
  
  type, extends(Stringable) :: CoupledModes
    ! The ids of the modes which are coupled together.
    ! ids is assumed to always be sorted in ascending order.
    integer, allocatable :: ids(:)
  contains
    ! I/O.
    procedure, public :: read  => read_CoupledModes
    procedure, public :: write => write_CoupledModes
  end type
  
  interface CoupledModes
    module procedure new_CoupledModes
  end interface
  
  interface size
    module procedure size_CoupledModes
  end interface
  
  interface operator(==)
    module procedure equality_CoupledModes_CoupledModes
  end interface
  
  interface operator(/=)
    module procedure non_equality_CoupledModes_CoupledModes
  end interface
contains

! ----------------------------------------------------------------------
! Constructor.
! ----------------------------------------------------------------------
function new_CoupledModes(ids) result(this)
  implicit none
  
  integer, intent(in) :: ids(:)
  type(CoupledModes)  :: this
  
  this%ids = ids
end function

! ----------------------------------------------------------------------
! Basis type functionality: size() and comparison operators.
! ----------------------------------------------------------------------
function size_CoupledModes(this) result(output)
  implicit none
  
  type(CoupledModes), intent(in) :: this
  integer                        :: output
  
  output = size(this%ids)
end function

impure elemental function equality_CoupledModes_CoupledModes(this,that) &
   & result(output)
  implicit none
  
  type(CoupledModes), intent(in) :: this
  type(CoupledModes), intent(in) :: that
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

impure elemental function non_equality_CoupledModes_CoupledModes(this,that) &
   & result(output)
  implicit none
  
  type(CoupledModes), intent(in) :: this
  type(CoupledModes), intent(in) :: that
  logical                        :: output
  
  output = .not. this==that
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_CoupledModes(this,input)
  implicit none
  
  class(CoupledModes), intent(out) :: this
  type(String),        intent(in)  :: input
  
  select type(this); type is(CoupledModes)
    this = CoupledModes(int(split_line(input)))
  end select
end subroutine

function write_CoupledModes(this) result(output)
  implicit none
  
  class(CoupledModes), intent(in) :: this
  type(String)                    :: output
  
  select type(this); type is(CoupledModes)
    output = join(this%ids)
  end select
end function
end module
