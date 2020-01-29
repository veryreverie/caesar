! ======================================================================
! The harmonic approximation to the stress, in cartesian co-ordinates.
! ======================================================================
module stress_hessian_module
  use utils_module
  
  use structure_module
  use normal_mode_module
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
    module procedure new_StressHessian
    module procedure new_StressHessian_Strings
    module procedure new_StressHessian_StringArray
  end interface
contains

function new_StressHessian(elements) result(this)
  implicit none
  
  type(CartesianHessian), intent(in) :: elements(:,:)
  type(StressHessian)                :: this
  
  this%elements = elements
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_StressHessian(this,input)
  implicit none
  
  class(StressHessian), intent(out) :: this
  type(String),         intent(in)  :: input(:)
  
  type(CartesianHessian), allocatable :: elements(:,:)
  
  select type(this); type is(StressHessian)
    elements = reshape(                            &
       & CartesianHessian(split_into_sections(     &
       &        input,                             &
       &        separating_line=repeat('-',50) )), &
       & [3,3]                                     )
    this = StressHessian(elements)
  class default
    call err()
  end select
end subroutine

function write_StressHessian(this) result(output)
  implicit none
  
  class(StressHessian), intent(in) :: this
  type(String), allocatable        :: output(:)
  
  select type(this); type is(StressHessian)
    output = str( [ this%elements(:,1),          &
                &   this%elements(:,2),          &
                &   this%elements(:,3)  ],       &
                & separating_line=repeat('-',50) )
  class default
    call err()
  end select
end function

function new_StressHessian_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(StressHessian)      :: this
  
  call this%read(input)
end function

impure elemental function new_StressHessian_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(StressHessian)           :: this
  
  this = StressHessian(str(input))
end function
end module
