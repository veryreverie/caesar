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
  public :: conjg
  
  public :: operator(+)
  public :: operator(-)
  public :: operator(*)
  public :: operator(/)
  
  type, extends(Stringsable) :: StressDynamicalMatrix
    type(DynamicalMatrix), allocatable :: elements(:,:)
  contains
    ! I/O.
    procedure, public :: read  => read_StressDynamicalMatrix
    procedure, public :: write => write_StressDynamicalMatrix
  end type
  
  interface StressDynamicalMatrix
    module procedure new_StressDynamicalMatrix
    module procedure new_StressDynamicalMatrix_zeroes
    module procedure new_StressDynamicalMatrix_Strings
    module procedure new_StressDynamicalMatrix_StringArray
  end interface
  
  interface operator(+)
    module procedure add_StressDynamicalMatrix_StressDynamicalMatrix
  end interface
  
  interface operator(-)
    module procedure negative_StressDynamicalMatrix
    module procedure subtract_StressDynamicalMatrix_StressDynamicalMatrix
  end interface
  
  interface operator(*)
    module procedure multiply_StressDynamicalMatrix_real
    module procedure multiply_real_StressDynamicalMatrix
    module procedure multiply_StressDynamicalMatrix_complex
    module procedure multiply_complex_StressDynamicalMatrix
  end interface
  
  interface operator(/)
    module procedure divide_StressDynamicalMatrix_real
    module procedure divide_StressDynamicalMatrix_complex
  end interface
  
  interface conjg
    module procedure conjg_StressDynamicalMatrix
  end interface
contains

! Constructors.
function new_StressDynamicalMatrix(elements) result(this)
  implicit none
  
  type(DynamicalMatrix), intent(in) :: elements(:,:)
  type(StressDynamicalMatrix)       :: this
  
  this%elements = elements
end function

function new_StressDynamicalMatrix_zeroes(no_atoms) result(this)
  implicit none
  
  integer, intent(in)         :: no_atoms
  type(StressDynamicalMatrix) :: this
  
  integer :: ialloc
  
  allocate(this%elements(3,3), stat=ialloc); call err(ialloc)
  this%elements = DynamicalMatrix(no_atoms)
end function

! Algebra.
impure elemental function add_StressDynamicalMatrix_StressDynamicalMatrix( &
   & this,that) result(output) 
  implicit none
  
  type(StressDynamicalMatrix), intent(in) :: this
  type(StressDynamicalMatrix), intent(in) :: that
  type(StressDynamicalMatrix)             :: output
  
  output = StressDynamicalMatrix(this%elements + that%elements)
end function

impure elemental function negative_StressDynamicalMatrix(this) result(output)
  implicit none
  
  type(StressDynamicalMatrix), intent(in) :: this
  type(StressDynamicalMatrix)             :: output
  
  output = StressDynamicalMatrix(-this%elements)
end function

impure elemental function                                            &
   & subtract_StressDynamicalMatrix_StressDynamicalMatrix(this,that) &
   & result(output)
  implicit none
  
  type(StressDynamicalMatrix), intent(in) :: this
  type(StressDynamicalMatrix), intent(in) :: that
  type(StressDynamicalMatrix)             :: output
  
  output = StressDynamicalMatrix(this%elements - that%elements)
end function

impure elemental function multiply_StressDynamicalMatrix_real(this,that) &
   & result(output) 
  implicit none
  
  type(StressDynamicalMatrix), intent(in) :: this
  real(dp),                    intent(in) :: that
  type(StressDynamicalMatrix)             :: output
  
  output = StressDynamicalMatrix(this%elements*that)
end function

impure elemental function multiply_real_StressDynamicalMatrix(this,that) &
   & result(output) 
  implicit none
  
  real(dp),                    intent(in) :: this
  type(StressDynamicalMatrix), intent(in) :: that
  type(StressDynamicalMatrix)             :: output
  
  output = StressDynamicalMatrix(this*that%elements)
end function

impure elemental function multiply_StressDynamicalMatrix_complex(this,that) &
   & result(output) 
  implicit none
  
  type(StressDynamicalMatrix), intent(in) :: this
  complex(dp),                 intent(in) :: that
  type(StressDynamicalMatrix)             :: output
  
  output = StressDynamicalMatrix(this%elements*that)
end function

impure elemental function multiply_complex_StressDynamicalMatrix(this,that) &
   & result(output) 
  implicit none
  
  complex(dp),                 intent(in) :: this
  type(StressDynamicalMatrix), intent(in) :: that
  type(StressDynamicalMatrix)             :: output
  
  output = StressDynamicalMatrix(this*that%elements)
end function

impure elemental function divide_StressDynamicalMatrix_real(this,that) &
   & result(output) 
  implicit none
  
  type(StressDynamicalMatrix), intent(in) :: this
  real(dp),                    intent(in) :: that
  type(StressDynamicalMatrix)             :: output
  
  output = StressDynamicalMatrix(this%elements/that)
end function

impure elemental function divide_StressDynamicalMatrix_complex(this,that) &
   & result(output) 
  implicit none
  
  type(StressDynamicalMatrix), intent(in) :: this
  complex(dp),                 intent(in) :: that
  type(StressDynamicalMatrix)             :: output
  
  output = StressDynamicalMatrix(this%elements/that)
end function

impure elemental function conjg_StressDynamicalMatrix(input) result(output)
  implicit none
  
  type(StressDynamicalMatrix), intent(in) :: input
  type(StressDynamicalMatrix)             :: output
  
  output = StressDynamicalMatrix(conjg(input%elements))
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
