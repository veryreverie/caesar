! ======================================================================
! The harmonic approximation to the stress, in normal-mode co-ordinates.
! ======================================================================
module stress_dynamical_matrix_module
  use utils_module
  
  use structure_module
  use dynamical_matrices_module
  implicit none
  
  private
  
  public :: StressDynamicalMatrix
  
  type, extends(Stringsable) :: StressDynamicalMatrix
    type(DynamicalMatrix), allocatable :: elements(:,:)
  contains
    ! I/O.
    procedure, public :: read  => read_StressDynamicalMatrix
    procedure, public :: write => write_StressDynamicalMatrix
  end type
  
  interface StressDynamicalMatrix
    module procedure new_StressDynamicalMatrix
    module procedure new_StressDynamicalMatrix_Strings
    module procedure new_StressDynamicalMatrix_StringArray
  end interface
contains

function new_StressDynamicalMatrix(elements) result(this)
  implicit none
  
  type(DynamicalMatrix), intent(in) :: elements(:,:)
  type(StressDynamicalMatrix)       :: this
  
  this%elements = elements
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_StressDynamicalMatrix(this,input)
  implicit none
  
  class(StressDynamicalMatrix), intent(out) :: this
  type(String),                 intent(in)  :: input(:)
  
  type(DynamicalMatrix), allocatable :: elements(:,:)
  
  select type(this); type is(StressDynamicalMatrix)
    elements = reshape(                        &
       & DynamicalMatrix(split_into_sections(  &
       &    input,                             &
       &    separating_line=repeat('-',50) )), &
       & [3,3]                                 )
    this = StressDynamicalMatrix(elements)
  class default
    call err()
  end select
end subroutine

function write_StressDynamicalMatrix(this) result(output)
  implicit none
  
  class(StressDynamicalMatrix), intent(in) :: this
  type(String), allocatable                :: output(:)
  
  select type(this); type is(StressDynamicalMatrix)
    output = str( [ this%elements(:,1),          &
                &   this%elements(:,2),          &
                &   this%elements(:,3)  ],       &
                & separating_line=repeat('-',50) )
  class default
    call err()
  end select
end function

function new_StressDynamicalMatrix_Strings(input) result(this)
  implicit none
  
  type(String), intent(in)    :: input(:)
  type(StressDynamicalMatrix) :: this
  
  call this%read(input)
end function

impure elemental function new_StressDynamicalMatrix_StringArray(input) &
   & result(this) 
  implicit none
  
  type(StringArray), intent(in) :: input
  type(StressDynamicalMatrix)   :: this
  
  this = StressDynamicalMatrix(str(input))
end function
end module
