! ======================================================================
! A term of the form (u_i)**n_i.
! ======================================================================
module univariate_module
  use constants_module, only : dp
  use string_module
  use io_module
  
  use stringable_module
  implicit none
  
  private
  
  public :: Univariate
  
  type, extends(Stringable) :: Univariate
    integer :: id
    integer :: power
  contains
    procedure :: evaluate   => evaluate_Univariate
    procedure :: str        => str_Univariate
  end type
contains

! Evaluate the univariate at a give displacement.
! Checks that the displacement and univariate are in the same mode.
function evaluate_Univariate(this,displacement) result(output)
  use single_mode_displacement_module
  implicit none
  
  class(Univariate),            intent(in) :: this
  type(SingleModeDisplacement), intent(in) :: displacement
  complex(dp)                              :: output
  
  if (this%id/=displacement%id) then
    call print_line(CODE_ERROR//': Trying to evaluate a univariate at an &
       &incompatible displacement.')
    call err()
  endif
  
  output = displacement%displacement**this%power
end function

! I/O.
function str_Univariate(this) result(output)
  implicit none
  
  class(Univariate), intent(in) :: this
  type(String)                  :: output
  
  output = 'u'//this%id//'^'//this%power
end function
end module
