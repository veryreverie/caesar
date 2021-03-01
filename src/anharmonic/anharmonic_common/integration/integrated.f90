! ======================================================================
! A type which can be a scalar, vector or tensor.
! Allows for polymorphic integration procedures.
! ======================================================================
module caesar_integrated_module
  use caesar_common_module
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
    impure elemental module function new_Integrated_scalar(input) result(this) 
      real(dp), intent(in) :: input
      type(Integrated)     :: this
    end function
  
    impure elemental module function new_Integrated_vector(input) result(this) 
      type(RealVector), intent(in) :: input
      type(Integrated)             :: this
    end function
  
    impure elemental module function new_Integrated_tensor(input) result(this) 
      type(RealMatrix), intent(in) :: input
      type(Integrated)             :: this
    end function
  end interface
  
  interface
    impure elemental module function scalar_Integrated(this) result(output) 
      class(Integrated), intent(in) :: this
      real(dp)                      :: output
    end function
  end interface
  
  interface
    impure elemental module function vector_Integrated(this) result(output) 
      class(Integrated), intent(in) :: this
      type(RealVector)             :: output
    end function
  end interface
  
  interface
    impure elemental module function tensor_Integrated(this) result(output) 
      class(Integrated), intent(in) :: this
      type(RealMatrix)              :: output
    end function
  end interface
  
  interface operator(+)
    impure elemental module function add_Integrated_Integrated(this,that) &
       & result(output) 
      class(Integrated), intent(in) :: this
      class(Integrated), intent(in) :: that
      type(Integrated)              :: output
    end function
  end interface
  
  interface sum
    module function sum_Integrateds(input) result(output) 
      class(Integrated), intent(in) :: input(:)
      type(Integrated)              :: output
    end function
  end interface
end module
