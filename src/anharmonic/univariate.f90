! ======================================================================
! A term of the form (u_i)**n_i.
! ======================================================================
module univariate_module
  use common_module
  
  use single_mode_displacement_module
  implicit none
  
  private
  
  public :: Univariate
  public :: operator(==)
  public :: operator(/=)
  
  type, extends(Stringable) :: Univariate
    integer :: id
    integer :: power
  contains
    procedure :: evaluate   => evaluate_Univariate
    procedure :: str        => str_Univariate
  end type
  
  interface operator(==)
    module procedure equality_Univariate_Univariate
  end interface
  
  interface operator(/=)
    module procedure non_equality_Univariate_Univariate
  end interface
contains

! Evaluate the univariate at a give displacement.
! Checks that the displacement and univariate are in the same mode.
function evaluate_Univariate(this,displacement) result(output)
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

! Comparison between univariates.
impure elemental function equality_Univariate_Univariate(this,that) &
   & result(output)
  implicit none
  
  type(Univariate), intent(in) :: this
  type(Univariate), intent(in) :: that
  logical                      :: output
  
  output = this%id==that%id .and. this%power==that%power
end function

impure elemental function non_equality_Univariate_Univariate(this,that) &
   & result(output)
  implicit none
  
  type(Univariate), intent(in) :: this
  type(Univariate), intent(in) :: that
  logical                      :: output
  
  output = .not. this==that
end function

! I/O.
function str_Univariate(this) result(output)
  implicit none
  
  class(Univariate), intent(in) :: this
  type(String)                  :: output
  
  output = 'u'//this%id//'^'//this%power
end function
end module
