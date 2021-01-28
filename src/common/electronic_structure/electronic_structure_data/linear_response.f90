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
    module procedure new_LinearResponse
    module procedure new_LinearResponse_Strings
    module procedure new_LinearResponse_StringArray
  end interface
contains

function new_LinearResponse(permittivity,born_charges) result(this)
  implicit none
  
  type(RealMatrix), intent(in) :: permittivity
  type(RealMatrix), intent(in) :: born_charges(:)
  type(LinearResponse)         :: this
  
  this%permittivity = permittivity
  this%born_charges = born_charges
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_LinearResponse(this,input)
  implicit none
  
  class(LinearResponse), intent(out) :: this
  type(String),          intent(in)  :: input(:)
  
  type(RealMatrix)              :: permittivity
  type(RealMatrix), allocatable :: born_charges(:)
  
  select type(this); type is(LinearResponse)
    permittivity = RealMatrix(input(2:4))
    born_charges = RealMatrix(split_into_sections( &
                              & input(6:),         &
                              & separating_line='' ))
    
    this = LinearResponse(permittivity,born_charges)
  class default
    call err()
  end select
end subroutine

function write_LinearResponse(this) result(output)
  implicit none
  
  class(LinearResponse), intent(in) :: this
  type(String), allocatable         :: output(:)
  
  select type(this); type is(LinearResponse)
    output = [ str('Permittivity'),                      &
             & str(this%permittivity),                   &
             & str('Born Effective Charges'),            &
             & str(this%born_charges,separating_line='') ]
  class default
    call err()
  end select
end function

function new_LinearResponse_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(LinearResponse)     :: this
  
  call this%read(input)
end function

impure elemental function new_LinearResponse_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(LinearResponse)          :: this
  
  this = LinearResponse(str(input))
end function
end module
