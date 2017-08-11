! ======================================================================
! A heterogeneous multi-dimensional integer array.
! ======================================================================
module integer_arrays_module
  use constants_module, only : dp
  use string_module
  use io_module
  
  ! An array of integers.
  type :: IntArray1D
    integer, allocatable :: i(:)
  end type
  
  ! A array of IntArrays.
  type :: IntArray2D
    type(IntArray1D), allocatable :: i(:)
  end type
  
  ! Conversions from arrays to IntArray types.
  interface assignment(=)
    module procedure assign_IntArray1D_integers
    module procedure assign_IntArray2D_IntArray1Ds
  end interface
  
  interface array
    module procedure array_IntArray1D_integers
    module procedure array_IntArray2D_IntArray1Ds
  end interface
  
  ! Concatenation.
  interface operator(//)
    module procedure concatenate_IntArray1D_integers
    module procedure concatenate_integers_IntArray1D
    module procedure concatenate_IntArray1D_IntArray1D
    
    module procedure concatenate_IntArray2D_IntArray1Ds
    module procedure concatenate_IntArray1Ds_IntArray2D
    module procedure concatenate_IntArray2D_IntArray2D
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
subroutine assign_IntArray1D_integers(output,input)
  implicit none
  
  integer,          intent(in)  :: input(:)
  type(IntArray1D), intent(out) :: output
  
  output%i = input
end subroutine

subroutine assign_IntArray2D_IntArray1Ds(output,input)
  implicit none
  
  type(IntArray1D), intent(in)  :: input(:)
  type(IntArray2D), intent(out) :: output
  
  output%i = input
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
function concatenate_IntArray1D_integers(a,b) result(output)
  implicit none
  
  type(IntArray1D), intent(in) :: a
  integer,          intent(in) :: b(:)
  type(IntArray1D)             :: output
  
  integer :: ialloc
  
  allocate(output%i(size(a)+size(b)), stat=ialloc); call err(ialloc)
  output%i(:size(a)) = a%i
  output%i(size(a)+1:) = b
end function

! IntArray1D = integer(:) // IntArray1D
function concatenate_integers_IntArray1D(a,b) result(output)
  implicit none
  
  integer,          intent(in) :: a(:)
  type(IntArray1D), intent(in) :: b
  type(IntArray1D)             :: output
  
  integer :: ialloc
  
  allocate(output%i(size(a)+size(b)), stat=ialloc); call err(ialloc)
  output%i(:size(a)) = a
  output%i(size(a)+1:) = b%i
end function

! IntArray1D = IntArray1D // IntArray1D
function concatenate_IntArray1D_IntArray1D(a,b) result(output)
  implicit none
  
  type(IntArray1D), intent(in) :: a
  type(IntArray1D), intent(in) :: b
  type(IntArray1D)             :: output
  
  integer :: ialloc
  
  allocate(output%i(size(a)+size(b)), stat=ialloc); call err(ialloc)
  output%i(:size(a)) = a%i
  output%i(size(a)+1:) = b%i
end function

! IntArray2D = IntArray2D // IntArray1D(:)
function concatenate_IntArray2D_IntArray1Ds(a,b) result(output)
  implicit none
  
  type(IntArray2D), intent(in) :: a
  type(IntArray1D), intent(in) :: b(:)
  type(IntArray2D)             :: output
  
  integer :: ialloc
  
  allocate(output%i(size(a)+size(b)), stat=ialloc); call err(ialloc)
  output%i(:size(a)) = a%i
  output%i(size(a)+1:) = b
end function

! IntArray2D = IntArray1D(:) // IntArray2D
function concatenate_IntArray1Ds_IntArray2D(a,b) result(output)
  implicit none
  
  type(IntArray1D), intent(in) :: a(:)
  type(IntArray2D), intent(in) :: b
  type(IntArray2D)             :: output
  
  integer :: ialloc
  
  allocate(output%i(size(a)+size(b)), stat=ialloc); call err(ialloc)
  output%i(:size(a)) = a
  output%i(size(a)+1:) = b%i
end function

! IntArray2D = IntArray2D // IntArray2D
function concatenate_IntArray2D_IntArray2D(a,b) result(output)
  implicit none
  
  type(IntArray2D), intent(in) :: a
  type(IntArray2D), intent(in) :: b
  type(IntArray2D)             :: output
  
  integer :: ialloc
  
  allocate(output%i(size(a)+size(b)), stat=ialloc); call err(ialloc)
  output%i(:size(a)) = a%i
  output%i(size(a)+1:) = b%i
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
