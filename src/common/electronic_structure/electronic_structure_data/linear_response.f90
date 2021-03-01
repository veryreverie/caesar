! ======================================================================
! Linear response data, including permittivity and Born effective charges.
! ======================================================================
module caesar_linear_response_module
  use caesar_utils_module
  
  use caesar_structure_module
  use caesar_normal_mode_module
  implicit none
  
  private
  
  public :: LinearResponse
  
  type, extends(Stringsable) :: LinearResponse
    type(realMatrix)              :: permittivity
    type(RealMatrix), allocatable :: born_charges(:)
  contains
    procedure, public :: read  => read_LinearResponse
    procedure, public :: write => write_LinearResponse
  end type
  
  interface LinearResponse
    module function new_LinearResponse(permittivity,born_charges) result(this) 
      type(RealMatrix), intent(in) :: permittivity
      type(RealMatrix), intent(in) :: born_charges(:)
      type(LinearResponse)         :: this
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_LinearResponse(this,input) 
      class(LinearResponse), intent(out) :: this
      type(String),          intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_LinearResponse(this) result(output) 
      class(LinearResponse), intent(in) :: this
      type(String), allocatable         :: output(:)
    end function
  end interface
  
  interface LinearResponse
    module function new_LinearResponse_Strings(input) result(this) 
      type(String), intent(in) :: input(:)
      type(LinearResponse)     :: this
    end function
  
    impure elemental module function new_LinearResponse_StringArray(input) &
       & result(this) 
      type(StringArray), intent(in) :: input
      type(LinearResponse)          :: this
    end function
  end interface
end module
