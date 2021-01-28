! ======================================================================
! Vectors and matrices in sparse representation.
! ======================================================================
module caesar_sparse_algebra_module
  use caesar_precision_module
  use caesar_io_module
  
  use caesar_linear_algebra_module
  implicit none
  
  public :: IntElement
  public :: RealElement
  public :: ComplexElement
  
  public :: SparseIntMatrix
  public :: SparseRealMatrix
  public :: SparseComplexMatrix
  
  public :: mat
  
  type, extends(Stringable) :: IntElement
    integer :: i
    integer :: j
    integer :: element
  contains
    procedure, public :: read  => read_IntElement
    procedure, public :: write => write_IntElement
  end type
  
  interface IntElement
    module procedure new_IntElement
    module procedure new_IntElement_String
  end interface
  
  type, extends(Stringable) :: RealElement
    integer  :: i
    integer  :: j
    real(dp) :: element
  contains
    procedure, public :: read  => read_RealElement
    procedure, public :: write => write_RealElement
  end type
  
  interface RealElement
    module procedure new_RealElement
    module procedure new_RealElement_String
  end interface
  
  type, extends(Stringable) :: ComplexElement
    integer     :: i
    integer     :: j
    complex(dp) :: element
  contains
    procedure, public :: read  => read_ComplexElement
    procedure, public :: write => write_ComplexElement
  end type
  
  interface ComplexElement
    module procedure new_ComplexElement
    module procedure new_ComplexElement_String
  end interface
  
  type, extends(Stringsable) :: SparseIntMatrix
    integer,                       private :: no_rows_
    integer,                       private :: no_cols_
    integer,                       private :: no_elements_
    type(IntElement), allocatable, private :: elements_(:)
  contains
    procedure, public :: elements => elements_SparseIntMatrix
    
    procedure, public :: add_element => add_element_SparseIntMatrix
    procedure, public :: add_elements => add_elements_SparseIntMatrix
    
    procedure, public :: read  => read_SparseIntMatrix
    procedure, public :: write => write_SparseIntMatrix
  end type
  
  interface SparseIntMatrix
    module procedure new_SparseIntMatrix
    module procedure new_SparseIntMatrix_Strings
    module procedure new_SparseIntMatrix_StringArray
  end interface
  
  type, extends(Stringsable) :: SparseRealMatrix
    integer,                        private :: no_rows_
    integer,                        private :: no_cols_
    integer,                        private :: no_elements_
    type(RealElement), allocatable, private :: elements_(:)
  contains
    procedure, public :: elements => elements_SparseRealMatrix
    
    procedure, public :: add_element => add_element_SparseRealMatrix
    procedure, public :: add_elements => add_elements_SparseRealMatrix
    
    procedure, public :: read  => read_SparseRealMatrix
    procedure, public :: write => write_SparseRealMatrix
  end type
  
  interface SparseRealMatrix
    module procedure new_SparseRealMatrix
    module procedure new_SparseRealMatrix_Strings
    module procedure new_SparseRealMatrix_StringArray
  end interface
  
  type, extends(Stringsable) :: SparseComplexMatrix
    integer,                           private :: no_rows_
    integer,                           private :: no_cols_
    integer,                           private :: no_elements_
    type(ComplexElement), allocatable, private :: elements_(:)
  contains
    procedure, public :: elements => elements_SparseComplexMatrix
    
    procedure, public :: add_element => add_element_SparseComplexMatrix
    procedure, public :: add_elements => add_elements_SparseComplexMatrix
    
    procedure, public :: read  => read_SparseComplexMatrix
    procedure, public :: write => write_SparseComplexMatrix
  end type
  
  interface SparseComplexMatrix
    module procedure new_SparseComplexMatrix
    module procedure new_SparseComplexMatrix_Strings
    module procedure new_SparseComplexMatrix_StringArray
  end interface
  
  ! Convert to dense representation.
  interface mat
    module procedure mat_SparseIntMatrix
    module procedure mat_SparseRealMatrix
    module procedure mat_SparseComplexMatrix
  end interface

! Include preprocessed procedure interfaces.
#include "sparse_algebra_includes.fpp"

contains

! Include preprocessed procedure bodies.
#define MACRO_BODY
#include "sparse_algebra_includes.fpp"

! Constructors.
impure elemental function new_IntElement(i,j,element) result(this)
  implicit none
  
  integer, intent(in) :: i
  integer, intent(in) :: j
  integer, intent(in) :: element
  type(IntElement)    :: this
  
  this%i = i
  this%j = j
  this%element = element
end function

impure elemental function new_RealElement(i,j,element) result(this)
  implicit none
  
  integer,  intent(in) :: i
  integer,  intent(in) :: j
  real(dp), intent(in) :: element
  type(RealElement)    :: this
  
  this%i = i
  this%j = j
  this%element = element
end function

impure elemental function new_ComplexElement(i,j,element) result(this)
  implicit none
  
  integer,     intent(in) :: i
  integer,     intent(in) :: j
  complex(dp), intent(in) :: element
  type(ComplexElement)    :: this
  
  this%i = i
  this%j = j
  this%element = element
end function

function new_SparseIntMatrix(no_rows,no_cols,elements) result(this)
  implicit none
  
  integer,          intent(in)           :: no_rows
  integer,          intent(in)           :: no_cols
  type(IntElement), intent(in), optional :: elements(:)
  type(SparseIntMatrix)                  :: this
  
  integer :: ialloc
  
  this%no_rows_ = no_rows
  this%no_cols_ = no_cols
  if (present(elements)) then
    this%no_elements_ = size(elements)
    this%elements_ = elements
  else
    this%no_elements_ = 0
    allocate(this%elements_(0), stat=ialloc); call err(ialloc)
  endif
end function

function new_SparseRealMatrix(no_rows,no_cols,elements) result(this)
  implicit none
  
  integer,           intent(in)           :: no_rows
  integer,           intent(in)           :: no_cols
  type(RealElement), intent(in), optional :: elements(:)
  type(SparseRealMatrix)                  :: this
  
  integer :: ialloc
  
  this%no_rows_ = no_rows
  this%no_cols_ = no_cols
  if (present(elements)) then
    this%no_elements_ = size(elements)
    this%elements_ = elements
  else
    this%no_elements_ = 0
    allocate(this%elements_(0), stat=ialloc); call err(ialloc)
  endif
end function

function new_SparseComplexMatrix(no_rows,no_cols,elements) result(this)
  implicit none
  
  integer,              intent(in)           :: no_rows
  integer,              intent(in)           :: no_cols
  type(ComplexElement), intent(in), optional :: elements(:)
  type(SparseComplexMatrix)                  :: this
  
  integer :: ialloc
  
  this%no_rows_ = no_rows
  this%no_cols_ = no_cols
  if (present(elements)) then
    this%no_elements_ = size(elements)
    this%elements_ = elements
  else
    this%no_elements_ = 0
    allocate(this%elements_(0), stat=ialloc); call err(ialloc)
  endif
end function

! Getters and setters for elements.
function elements_SparseIntMatrix(this) result(output)
  implicit none
  
  class(SparseIntMatrix), intent(in) :: this
  type(IntElement), allocatable      :: output(:)
  
  output = this%elements_(:this%no_elements_)
end function

subroutine add_element_SparseIntMatrix(this,element)
  implicit none
  
  class(SparseIntMatrix), intent(inout) :: this
  type(IntElement),       intent(in)    :: element
  
  call this%add_elements([element])
end subroutine

subroutine add_elements_SparseIntMatrix(this,elements)
  implicit none
  
  class(SparseIntMatrix), intent(inout) :: this
  type(IntElement),       intent(in)    :: elements(:)
  
  type(IntElement), allocatable :: temp(:)
  
  integer :: ialloc
  
  if (this%no_elements_+size(elements)>size(this%elements_)) then
    allocate( temp(2*(this%no_elements_+size(elements))), &
            & stat=ialloc); call err(ialloc)
    temp(:this%no_elements_) = this%elements_(:this%no_elements_)
    this%elements_ = temp
  endif
  
  temp(this%no_elements_+1:this%no_elements_+size(elements)) = elements
  this%no_elements_ = this%no_elements_+size(elements)
end subroutine

function elements_SparseRealMatrix(this) result(output)
  implicit none
  
  class(SparseRealMatrix), intent(in) :: this
  type(RealElement), allocatable      :: output(:)
  
  output = this%elements_(:this%no_elements_)
end function

subroutine add_element_SparseRealMatrix(this,element)
  implicit none
  
  class(SparseRealMatrix), intent(inout) :: this
  type(RealElement),       intent(in)    :: element
  
  call this%add_elements([element])
end subroutine

subroutine add_elements_SparseRealMatrix(this,elements)
  implicit none
  
  class(SparseRealMatrix), intent(inout) :: this
  type(RealElement),       intent(in)    :: elements(:)
  
  type(RealElement), allocatable :: temp(:)
  
  integer :: ialloc
  
  if (this%no_elements_+size(elements)>size(this%elements_)) then
    allocate( temp(2*(this%no_elements_+size(elements))), &
            & stat=ialloc); call err(ialloc)
    temp(:this%no_elements_) = this%elements_(:this%no_elements_)
    this%elements_ = temp
  endif
  
  temp(this%no_elements_+1:this%no_elements_+size(elements)) = elements
  this%no_elements_ = this%no_elements_+size(elements)
end subroutine

function elements_SparseComplexMatrix(this) result(output)
  implicit none
  
  class(SparseComplexMatrix), intent(in) :: this
  type(ComplexElement), allocatable      :: output(:)
  
  output = this%elements_(:this%no_elements_)
end function

subroutine add_element_SparseComplexMatrix(this,element)
  implicit none
  
  class(SparseComplexMatrix), intent(inout) :: this
  type(ComplexElement),       intent(in)    :: element
  
  call this%add_elements([element])
end subroutine

subroutine add_elements_SparseComplexMatrix(this,elements)
  implicit none
  
  class(SparseComplexMatrix), intent(inout) :: this
  type(ComplexElement),       intent(in)    :: elements(:)
  
  type(ComplexElement), allocatable :: temp(:)
  
  integer :: ialloc
  
  if (this%no_elements_+size(elements)>size(this%elements_)) then
    allocate( temp(2*(this%no_elements_+size(elements))), &
            & stat=ialloc); call err(ialloc)
    temp(:this%no_elements_) = this%elements_(:this%no_elements_)
    this%elements_ = temp
  endif
  
  temp(this%no_elements_+1:this%no_elements_+size(elements)) = elements
  this%no_elements_ = this%no_elements_+size(elements)
end subroutine

! Conversion to dense representation.
impure elemental function mat_SparseIntMatrix(this) result(output)
  implicit none
  
  class(SparseIntMatrix), intent(in) :: this
  type(IntMatrix)                    :: output
  
  integer, allocatable :: elements(:,:)
  
  integer :: i,ialloc
  
  allocate( elements(this%no_rows_, this%no_cols_), &
          & stat=ialloc); call err(ialloc)
  elements = 0
  do i=1,this%no_elements_
    associate(element=>this%elements_(i))
      elements(element%i,element%j) = element%element
    end associate
  enddo
  
  output = mat(elements)
end function

impure elemental function mat_SparseRealMatrix(this) result(output)
  implicit none
  
  class(SparseRealMatrix), intent(in) :: this
  type(RealMatrix)                    :: output
  
  real(dp), allocatable :: elements(:,:)
  
  integer :: i,ialloc
  
  allocate( elements(this%no_rows_, this%no_cols_), &
          & stat=ialloc); call err(ialloc)
  elements = 0
  do i=1,this%no_elements_
    associate(element=>this%elements_(i))
      elements(element%i,element%j) = element%element
    end associate
  enddo
  
  output = mat(elements)
end function

impure elemental function mat_SparseComplexMatrix(this) result(output)
  implicit none
  
  class(SparseComplexMatrix), intent(in) :: this
  type(ComplexMatrix)                    :: output
  
  complex(dp), allocatable :: elements(:,:)
  
  integer :: i,ialloc
  
  allocate( elements(this%no_rows_, this%no_cols_), &
          & stat=ialloc); call err(ialloc)
  elements = 0
  do i=1,this%no_elements_
    associate(element=>this%elements_(i))
      elements(element%i,element%j) = element%element
    end associate
  enddo
  
  output = mat(elements)
end function

! I/O.
subroutine read_IntElement(this,input)
  implicit none
  
  class(IntElement), intent(out) :: this
  type(String),      intent(in)  :: input
  
  integer :: i
  integer :: j
  integer :: element
  
  select type(this); type is(IntElement)
    i = int(token(input,1))
    j = int(token(input,2))
    element = int(token(input,3))
    
    this = IntElement(i,j,element)
  class default
    call err()
  end select
end subroutine

function write_IntElement(this) result(output)
  implicit none
  
  class(IntElement), intent(in) :: this
  type(String)                  :: output
  
  select type(this); type is(IntElement)
    output = this%i//' '//this%j//' '//this%element
  class default
    call err()
  end select
end function

impure elemental function new_IntElement_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(IntElement)         :: this
  
  call this%read(input)
end function

subroutine read_RealElement(this,input)
  implicit none
  
  class(RealElement), intent(out) :: this
  type(String),       intent(in)  :: input
  
  integer  :: i
  integer  :: j
  real(dp) :: element
  
  select type(this); type is(RealElement)
    i = int(token(input,1))
    j = int(token(input,2))
    element = int(token(input,3))
    
    this = RealElement(i,j,element)
  class default
    call err()
  end select
end subroutine

function write_RealElement(this) result(output)
  implicit none
  
  class(RealElement), intent(in) :: this
  type(String)                   :: output
  
  select type(this); type is(RealElement)
    output = this%i//' '//this%j//' '//this%element
  class default
    call err()
  end select
end function

impure elemental function new_RealElement_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(RealElement)        :: this
  
  call this%read(input)
end function

subroutine read_ComplexElement(this,input)
  implicit none
  
  class(ComplexElement), intent(out) :: this
  type(String),          intent(in)  :: input
  
  integer     :: i
  integer     :: j
  complex(dp) :: element
  
  select type(this); type is(ComplexElement)
    i = int(token(input,1))
    j = int(token(input,2))
    element = int(token(input,3))
    
    this = ComplexElement(i,j,element)
  class default
    call err()
  end select
end subroutine

function write_ComplexElement(this) result(output)
  implicit none
  
  class(ComplexElement), intent(in) :: this
  type(String)                      :: output
  
  select type(this); type is(ComplexElement)
    output = this%i//' '//this%j//' '//this%element
  class default
    call err()
  end select
end function

impure elemental function new_ComplexElement_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(ComplexElement)     :: this
  
  call this%read(input)
end function

subroutine read_SparseIntMatrix(this,input)
  implicit none
  
  class(SparseIntMatrix), intent(out) :: this
  type(String),           intent(in)  :: input(:)
  
  integer                       :: no_rows
  integer                       :: no_cols
  type(IntElement), allocatable :: elements(:)
  
  select type(this); type is(SparseIntMatrix)
    no_rows = int(token(input(2),1))
    no_cols = int(token(input(2),4))
    elements = IntElement(input(4:))
    this = SparseIntMatrix(no_rows,no_cols,elements)
  class default
    call err()
  end select
end subroutine

function write_SparseIntMatrix(this) result(output)
  implicit none
  
  class(SparseIntMatrix), intent(in) :: this
  type(String), allocatable          :: output(:)
  
  select type(this); type is(SparseIntMatrix)
    output = [ str('Sparse Matrix'),                                  &
             & this%no_rows_//' rows by '//this%no_cols_//' columns', &
             & str('i | j | element(i,j)'),                           &
             & str(this%elements_(:this%no_elements_))                ]
  class default
    call err()
  end select
end function

function new_SparseIntMatrix_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(SparseIntMatrix)    :: this
  
  call this%read(input)
end function

impure elemental function new_SparseIntMatrix_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(SparseIntMatrix)         :: this
  
  this = SparseIntMatrix(str(input))
end function

subroutine read_SparseRealMatrix(this,input)
  implicit none
  
  class(SparseRealMatrix), intent(out) :: this
  type(String),            intent(in)  :: input(:)
  
  integer                        :: no_rows
  integer                        :: no_cols
  type(RealElement), allocatable :: elements(:)
  
  select type(this); type is(SparseRealMatrix)
    no_rows = int(token(input(2),1))
    no_cols = int(token(input(2),4))
    elements = RealElement(input(4:))
    this = SparseRealMatrix(no_rows,no_cols,elements)
  class default
    call err()
  end select
end subroutine

function write_SparseRealMatrix(this) result(output)
  implicit none
  
  class(SparseRealMatrix), intent(in) :: this
  type(String), allocatable           :: output(:)
  
  select type(this); type is(SparseRealMatrix)
    output = [ str('Sparse Matrix'),                                  &
             & this%no_rows_//' rows by '//this%no_cols_//' columns', &
             & str('i | j | element(i,j)'),                           &
             & str(this%elements_(:this%no_elements_))                ]
  class default
    call err()
  end select
end function

function new_SparseRealMatrix_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(SparseRealMatrix)   :: this
  
  call this%read(input)
end function

impure elemental function new_SparseRealMatrix_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(SparseRealMatrix)        :: this
  
  this = SparseRealMatrix(str(input))
end function

subroutine read_SparseComplexMatrix(this,input)
  implicit none
  
  class(SparseComplexMatrix), intent(out) :: this
  type(String),               intent(in)  :: input(:)
  
  integer                           :: no_rows
  integer                           :: no_cols
  type(ComplexElement), allocatable :: elements(:)
  
  select type(this); type is(SparseComplexMatrix)
    no_rows = int(token(input(2),1))
    no_cols = int(token(input(2),4))
    elements = ComplexElement(input(4:))
    this = SparseComplexMatrix(no_rows,no_cols,elements)
  class default
    call err()
  end select
end subroutine

function write_SparseComplexMatrix(this) result(output)
  implicit none
  
  class(SparseComplexMatrix), intent(in) :: this
  type(String), allocatable              :: output(:)
  
  select type(this); type is(SparseComplexMatrix)
    output = [ str('Sparse Matrix'),                                  &
             & this%no_rows_//' rows by '//this%no_cols_//' columns', &
             & str('i | j | element(i,j)'),                           &
             & str(this%elements_(:this%no_elements_))                ]
  class default
    call err()
  end select
end function

function new_SparseComplexMatrix_Strings(input) result(this)
  implicit none
  
  type(String), intent(in)  :: input(:)
  type(SparseComplexMatrix) :: this
  
  call this%read(input)
end function

impure elemental function new_SparseComplexMatrix_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(SparseComplexMatrix)     :: this
  
  this = SparseComplexMatrix(str(input))
end function
end module
