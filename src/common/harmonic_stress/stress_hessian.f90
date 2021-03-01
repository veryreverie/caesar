! ======================================================================
! The harmonic approximation to the stress, in cartesian co-ordinates.
! ======================================================================
module caesar_stress_hessian_module
  use caesar_utils_module
  
  use caesar_structure_module
  use caesar_normal_mode_module
  implicit none
  
  private
  
  public :: StressHessian
  
  type, extends(Stringsable) :: StressHessian
    type(CartesianHessian), allocatable :: elements(:,:)
  contains
    ! I/O.
    procedure, public :: read  => read_StressHessian
    procedure, public :: write => write_StressHessian
  end type
  
  interface StressHessian
    module function new_StressHessian(elements) result(this) 
      type(CartesianHessian), intent(in) :: elements(:,:)
      type(StressHessian)                :: this
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_StressHessian(this,input) 
      class(StressHessian), intent(out) :: this
      type(String),         intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_StressHessian(this) result(output) 
      class(StressHessian), intent(in) :: this
      type(String), allocatable        :: output(:)
    end function
  end interface
  
  interface StressHessian
    module function new_StressHessian_Strings(input) result(this) 
      type(String), intent(in) :: input(:)
      type(StressHessian)      :: this
    end function
  
    impure elemental module function new_StressHessian_StringArray(input) &
       & result(this) 
      type(StringArray), intent(in) :: input
      type(StressHessian)           :: this
    end function
  end interface
end module
