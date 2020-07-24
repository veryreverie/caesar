! ======================================================================
! A type which can be a scalar, vector or tensor.
! Allows for polymorphic integration procedures.
! ======================================================================
module integrated_module
  use common_module
  implicit none
  
  public :: Integrated
  public :: operator(+)
  public :: sum
  
  type :: Integrated
    real(dp),         allocatable :: scalar_
    type(RealVector), allocatable :: vector_
    type(RealMatrix), allocatable :: tensor_
    contains
    procedure, public :: scalar => scalar_Integrated
    procedure, public :: vector => vector_Integrated
    procedure, public :: tensor => tensor_Integrated
  end type
  
  interface Integrated
    module procedure new_Integrated_scalar
    module procedure new_Integrated_vector
    module procedure new_Integrated_tensor
  end interface
  
  interface operator(+)
    module procedure add_Integrated_Integrated
  end interface
  
  interface sum
    module procedure sum_Integrateds
  end interface
contains
impure elemental function new_Integrated_scalar(input) result(this)
  implicit none
  
  real(dp), intent(in) :: input
  type(Integrated)     :: this
  
  this%scalar_ = input
end function

impure elemental function new_Integrated_vector(input) result(this)
  implicit none
  
  type(RealVector), intent(in) :: input
  type(Integrated)             :: this
  
  this%vector_ = input
end function

impure elemental function new_Integrated_tensor(input) result(this)
  implicit none
  
  type(RealMatrix), intent(in) :: input
  type(Integrated)             :: this
  
  this%tensor_ = input
end function

impure elemental function scalar_Integrated(this) result(output)
  implicit none
  
  class(Integrated), intent(in) :: this
  real(dp)                      :: output
  
  if (allocated(this%scalar_)) then
    output = this%scalar_
  else
    call err()
  endif
end function

impure elemental function vector_Integrated(this) result(output)
  implicit none
  
  class(Integrated), intent(in) :: this
  type(RealVector)             :: output
  
  if (allocated(this%vector_)) then
    output = this%vector_
  else
    call err()
  endif
end function

impure elemental function tensor_Integrated(this) result(output)
  implicit none
  
  class(Integrated), intent(in) :: this
  type(RealMatrix)              :: output
  
  if (allocated(this%tensor_)) then
    output = this%tensor_
  else
    call err()
  endif
end function

impure elemental function add_Integrated_Integrated(this,that) result(output)
  implicit none
  
  class(Integrated), intent(in) :: this
  class(Integrated), intent(in) :: that
  type(Integrated)              :: output
  
  if (allocated(this%scalar_)) then
    output = Integrated(this%scalar_+that%scalar_)
  elseif (allocated(this%vector_)) then
    output = Integrated(this%vector_+that%vector_)
  elseif (allocated(this%tensor_)) then
    output = Integrated(this%tensor_+that%tensor_)
  else
    call err()
  endif
end function

function sum_Integrateds(input) result(output)
  implicit none
  
  class(Integrated), intent(in) :: input(:)
  type(Integrated)              :: output
  
  integer :: i
  
  if (size(input)>0) then
    output = input(1)
    do i=2,size(input)
      output = output + input(i)
    enddo
  else
    call err()
  endif
end function
end module
