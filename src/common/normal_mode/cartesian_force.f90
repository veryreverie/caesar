! ======================================================================
! A force in cartesian co-ordinates.
! ======================================================================
module caesar_cartesian_force_module
  use caesar_utils_module
  
  use caesar_structure_module
  implicit none
  
  private
  
  public :: CartesianForce
  public :: size
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: operator(-)
  public :: sum
  
  type, extends(Stringsable) :: CartesianForce
    type(RealVector), allocatable :: vectors(:)
  contains
    procedure, public :: read  => read_CartesianForce
    procedure, public :: write => write_CartesianForce
  end type
  
  interface CartesianForce
    ! ----------------------------------------------------------------------
    ! Constructor and size() module function.
    ! ----------------------------------------------------------------------
    module function new_CartesianForce(forces) result(this) 
      type(RealVector), intent(in) :: forces(:)
      type(CartesianForce)         :: this
    end function
  end interface
  
  interface size
    module function size_CartesianForce(this) result(output) 
      class(CartesianForce), intent(in) :: this
      integer                           :: output
    end function
  end interface
  
  interface CartesianForce
    ! ----------------------------------------------------------------------
    ! Construct a zero force.
    ! ----------------------------------------------------------------------
    impure elemental module function new_CartesianForce_zero(structure) &
       & result(this) 
      type(StructureData), intent(in) :: structure
      type(CartesianForce)            :: this
    end function
  end interface
  
  interface operator(*)
    ! ----------------------------------------------------------------------
    ! Algebra.
    ! ----------------------------------------------------------------------
    impure elemental module function multiply_real_CartesianForce(this,that) &
       & result(output) 
      real(dp),             intent(in) :: this
      type(CartesianForce), intent(in) :: that
      type(CartesianForce)             :: output
    end function
  
    impure elemental module function multiply_CartesianForce_real(this,that) &
       & result(output) 
      type(CartesianForce), intent(in) :: this
      real(dp),             intent(in) :: that
      type(CartesianForce)             :: output
    end function
  end interface
  
  interface operator(/)
    impure elemental module function divide_CartesianForce_real(this,that) &
       & result(output) 
      type(CartesianForce), intent(in) :: this
      real(dp),             intent(in) :: that
      type(CartesianForce)             :: output
    end function
  end interface
  
  interface operator(+)
    impure elemental module function add_CartesianForce_CartesianForce(this, &
       & that) result(output) 
      type(CartesianForce), intent(in) :: this
      type(CartesianForce), intent(in) :: that
      type(CartesianForce)             :: output
    end function
  end interface
  
  interface sum
    module function sum_CartesianForces(input) result(output) 
      type(CartesianForce), intent(in) :: input(:)
      type(CartesianForce)             :: output
    end function
  end interface
  
  interface operator(-)
    impure elemental module function negative_CartesianForce(this) &
       & result(output) 
      type(CartesianForce), intent(in) :: this
      type(CartesianForce)             :: output
    end function
  
    impure elemental module function subtract_CartesianForce_CartesianForce(this,that) result(output) 
      type(CartesianForce), intent(in) :: this
      type(CartesianForce), intent(in) :: that
      type(CartesianForce)             :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_CartesianForce(this,input) 
      class(CartesianForce), intent(out) :: this
      type(String),          intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_CartesianForce(this) result(output) 
      class(CartesianForce), intent(in) :: this
      type(String), allocatable         :: output(:)
    end function
  end interface
  
  interface CartesianForce
    module function new_CartesianForce_Strings(input) result(this) 
      type(String), intent(in) :: input(:)
      type(CartesianForce)     :: this
    end function
  
    impure elemental module function new_CartesianForce_StringArray(input) &
       & result(this) 
      type(StringArray), intent(in) :: input
      type(CartesianForce)          :: this
    end function
  end interface
end module
