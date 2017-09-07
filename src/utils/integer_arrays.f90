! ======================================================================
! A heterogeneous multi-dimensional integer array.
! ======================================================================
module integer_arrays_module
  use constants_module, only : dp
  use string_module
  use io_module
  
  private
  
  public :: IntArray1D
  public :: IntArray2D
  
  public :: array
  public :: size
  
  ! An array of integers.
  type :: IntArray1D
    integer, public, allocatable :: i(:)
  contains
    generic, public :: assignment (= ) => assign_IntArray1D_integers
    generic, public :: operator   (//) => concatenate_IntArray1D_integers, &
                                        & concatenate_integers_IntArray1D, &
                                        & concatenate_IntArray1D_IntArray1D
    
    procedure, private             :: assign_IntArray1D_integers
    procedure, private             :: concatenate_IntArray1D_integers
    procedure, private, pass(that) :: concatenate_integers_IntArray1D
    procedure, private             :: concatenate_IntArray1D_IntArray1D
  end type
  
  ! A array of IntArrays.
  type :: IntArray2D
    type(IntArray1D), allocatable :: i(:)
  contains
    generic, public :: assignment (= ) => assign_IntArray2D_IntArray1Ds
    generic, public :: operator   (//) => concatenate_IntArray2D_IntArray1Ds, &
                                        & concatenate_IntArray1Ds_IntArray2D, &
                                        & concatenate_IntArray2D_IntArray2D
    
    procedure, private             :: assign_IntArray2D_IntArray1Ds
    procedure, private             :: concatenate_IntArray2D_IntArray1Ds
    procedure, private, pass(that) :: concatenate_IntArray1Ds_IntArray2D
    procedure, private             :: concatenate_IntArray2D_IntArray2D
  end type
  
  ! Type conversion to array.
  interface array
    module procedure array_IntArray1D_integers
    module procedure array_IntArray2D_IntArray1Ds
  end interface
  
  ! Size.
  interface size
    module procedure size_IntArray1D
    module procedure size_IntArray2D
  end interface
  
contains

! ----------------------------------------------------------------------
! Conversions from integer(:) to IntArray1D and IntArray1D(:) to IntArray2D
! ----------------------------------------------------------------------
subroutine assign_IntArray1D_integers(this,that)
  implicit none
  
  class(IntArray1D), intent(out) :: this
  integer,           intent(in)  :: that(:)
  
  this%i = that
end subroutine

subroutine assign_IntArray2D_IntArray1Ds(this,that)
  implicit none
  
  class(IntArray2D), intent(out) :: this
  type(IntArray1D),  intent(in)  :: that(:)
  
  this%i = that
end subroutine

function array_IntArray1D_integers(input) result(output)
  implicit none
  
  integer, intent(in) :: input(:)
  type(IntArray1D)    :: output
  
  output = input
end function

function array_IntArray2D_IntArray1Ds(input) result(output)
  implicit none
  
  type(IntArray1D), intent(in) :: input(:)
  type(IntArray2D)             :: output
  
  output = input
end function

! ----------------------------------------------------------------------
! Concatenation operations
! ----------------------------------------------------------------------

! IntArray1D = IntArray1D // integer(:)
function concatenate_IntArray1D_integers(this,that) result(output)
  implicit none
  
  class(IntArray1D), intent(in) :: this
  integer,           intent(in) :: that(:)
  type(IntArray1D)              :: output
  
  integer :: ialloc
  
  allocate(output%i(size(this)+size(that)), stat=ialloc); call err(ialloc)
  output%i(            :size(this)) = this%i
  output%i(size(this)+1:          ) = that
end function

! IntArray1D = integer(:) // IntArray1D
function concatenate_integers_IntArray1D(this,that) result(output)
  implicit none
  
  integer,           intent(in) :: this(:)
  class(IntArray1D), intent(in) :: that
  type(IntArray1D)              :: output
  
  integer :: ialloc
  
  allocate(output%i(size(this)+size(that)), stat=ialloc); call err(ialloc)
  output%i(            :size(this)) = this
  output%i(size(this)+1:          ) = that%i
end function

! IntArray1D = IntArray1D // IntArray1D
function concatenate_IntArray1D_IntArray1D(this,that) result(output)
  implicit none
  
  class(IntArray1D), intent(in) :: this
  class(IntArray1D), intent(in) :: that
  type(IntArray1D)              :: output
  
  integer :: ialloc
  
  allocate(output%i(size(this)+size(that)), stat=ialloc); call err(ialloc)
  output%i(            :size(this)) = this%i
  output%i(size(this)+1:          ) = that%i
end function

! IntArray2D = IntArray2D // IntArray1D(:)
function concatenate_IntArray2D_IntArray1Ds(this,that) result(output)
  implicit none
  
  class(IntArray2D), intent(in) :: this
  type(IntArray1D),  intent(in) :: that(:)
  type(IntArray2D)              :: output
  
  integer :: ialloc
  
  allocate(output%i(size(this)+size(that)), stat=ialloc); call err(ialloc)
  output%i(            :size(this)) = this%i
  output%i(size(this)+1:          ) = that
end function

! IntArray2D = IntArray1D(:) // IntArray2D
function concatenate_IntArray1Ds_IntArray2D(this,that) result(output)
  implicit none
  
  type(IntArray1D),  intent(in) :: this(:)
  class(IntArray2D), intent(in) :: that
  type(IntArray2D)              :: output
  
  integer :: ialloc
  
  allocate(output%i(size(this)+size(that)), stat=ialloc); call err(ialloc)
  output%i(            :size(this)) = this
  output%i(size(this)+1:          ) = that%i
end function

! IntArray2D = IntArray2D // IntArray2D
function concatenate_IntArray2D_IntArray2D(this,that) result(output)
  implicit none
  
  class(IntArray2D), intent(in) :: this
  class(IntArray2D), intent(in) :: that
  type(IntArray2D)              :: output
  
  integer :: ialloc
  
  allocate(output%i(size(this)+size(that)), stat=ialloc); call err(ialloc)
  output%i(            :size(this)) = this%i
  output%i(size(this)+1:          ) = that%i
end function

! ----------------------------------------------------------------------
! Size function.
! ----------------------------------------------------------------------
function size_IntArray1D(input) result(output)
  implicit none
  
  type(IntArray1D), intent(in) :: input
  integer                      :: output
  
  output = size(input%i)
end function

function size_IntArray2D(input) result(output)
  implicit none
  
  type(IntArray2D), intent(in) :: input
  integer                      :: output
  
  output = size(input%i)
end function
end module
