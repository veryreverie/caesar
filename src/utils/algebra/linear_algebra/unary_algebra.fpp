! ======================================================================
! Templates for unary operations involving linear algebra types.
! When this file is included, the macro MACRO_TYPE_NAME must be defined.
! ======================================================================
! N.B. for each type, this file should be included twice:
!    - Before "contains", with MACRO_BODY not defined.
!    - After "contains", with MACRO_BODY defined.
!
! Type-bound functions:
!    - check, which checks that contents_ has been allocated.
!    - element, which returns a given element.
!    - row, which returns a given row (only applicable to matrices.)
!    - column, which returns a given column (only applicable to matrices.)
!
! Functions:
!    - int/dble/cmplx which return contents_ as an array.
!    - vec, which converts a 1-D array to a vector.
!    - mat, which converts a 2-D array to a matrix.
!    - row_matrix
!    - column_matrix
!    - size
!    - trace
!    - sum
!    - negative
!
! Also defines read, write and the standard I/O constructors.
! ======================================================================

! Include useful macros.
#include "macros.fpp"

! Define MACRO_TYPE_DECLARATION, MACRO_TYPE_VEC_NAME
!    and MACRO_TYPE_MAT_NAME.
#include "unary_type_declaration.fpp"

! Define function names.
#define MACRO_CHECK_VEC           MACRO_CAT2(check_,MACRO_TYPE_VEC_NAME)
#define MACRO_CHECK_MAT           MACRO_CAT2(check_,MACRO_TYPE_MAT_NAME)
#define MACRO_ELEMENT_VEC         MACRO_CAT2(element_,MACRO_TYPE_VEC_NAME)
#define MACRO_ELEMENT_MAT         MACRO_CAT2(element_,MACRO_TYPE_MAT_NAME)
#define MACRO_ROW_MAT             MACRO_CAT2(row_,MACRO_TYPE_MAT_NAME)
#define MACRO_COLUMN_MAT          MACRO_CAT2(column_,MACRO_TYPE_MAT_NAME)

#define MACRO_CONVERT_VEC         MACRO_CAT3(MACRO_TYPE_CONVERSION,_,MACRO_TYPE_VEC_NAME)
#define MACRO_CONVERT_MAT         MACRO_CAT3(MACRO_TYPE_CONVERSION,_,MACRO_TYPE_MAT_NAME)
#define MACRO_VEC_VEC             MACRO_CAT3(vec_,MACRO_TYPE_NAME,s)
#define MACRO_MAT_MAT             MACRO_CAT3(mat_,MACRO_TYPE_NAME,s)
#define MACRO_MAT_MAT_SHAPE       MACRO_CAT3(mat_,MACRO_TYPE_NAME,s_shape)
#define MACRO_ROW_MATRIX_VECS     MACRO_CAT3(row_matrix_,MACRO_TYPE_VEC_NAME,s)
#define MACRO_COLUMN_MATRIX_VECS  MACRO_CAT3(column_matrix_,MACRO_TYPE_VEC_NAME,s)
#define MACRO_SIZE_VEC            MACRO_CAT2(size_,MACRO_TYPE_VEC_NAME)
#define MACRO_SIZE_MAT            MACRO_CAT2(size_,MACRO_TYPE_MAT_NAME)
#define MACRO_TRACE_MAT           MACRO_CAT2(trace_,MACRO_TYPE_MAT_NAME)
#define MACRO_SUM_VECS            MACRO_CAT2(sum_,MACRO_TYPE_VEC_NAME)
#define MACRO_SUM_MATS            MACRO_CAT2(sum_,MACRO_TYPE_MAT_NAME)
#define MACRO_NEGATIVE_VEC        MACRO_CAT2(negative_,MACRO_TYPE_VEC_NAME)
#define MACRO_NEGATIVE_MAT        MACRO_CAT2(negative_,MACRO_TYPE_MAT_NAME)

#define MACRO_READ_VEC            MACRO_CAT2(read_,MACRO_TYPE_VEC_NAME)
#define MACRO_READ_MAT            MACRO_CAT2(read_,MACRO_TYPE_MAT_NAME)
#define MACRO_WRITE_VEC           MACRO_CAT2(write_,MACRO_TYPE_VEC_NAME)
#define MACRO_WRITE_MAT           MACRO_CAT2(write_,MACRO_TYPE_MAT_NAME)
#define MACRO_NEW_VEC_STRING      MACRO_CAT3(new_,MACRO_TYPE_VEC_NAME,_String)
#define MACRO_NEW_MAT_STRINGS     MACRO_CAT3(new_,MACRO_TYPE_MAT_NAME,_Strings)
#define MACRO_NEW_MAT_STRINGARRAY MACRO_CAT3(new_,MACRO_TYPE_MAT_NAME,_StringArray)

#ifndef MACRO_BODY
! ----------------------------------------------------------------------
! Code to appear in the header, before "contains".
! ----------------------------------------------------------------------

interface MACRO_TYPE_CONVERSION
  module procedure MACRO_CONVERT_VEC
  module procedure MACRO_CONVERT_MAT
end interface

interface vec
  module procedure MACRO_VEC_VEC
end interface

interface mat
  module procedure MACRO_MAT_MAT
  module procedure MACRO_MAT_MAT_SHAPE
end interface

interface row_matrix
  module procedure MACRO_ROW_MATRIX_VECS
end interface

interface column_matrix
  module procedure MACRO_COLUMN_MATRIX_VECS
end interface

interface size
  module procedure MACRO_SIZE_VEC
  module procedure MACRO_SIZE_MAT
end interface

interface trace
  module procedure MACRO_TRACE_MAT
end interface

interface sum
  module procedure MACRO_SUM_VECS
  module procedure MACRO_SUM_MATS
end interface

interface operator(-)
  module procedure MACRO_NEGATIVE_VEC
  module procedure MACRO_NEGATIVE_MAT
end interface

interface MACRO_TYPE_VEC_NAME
  module procedure MACRO_NEW_VEC_STRING
end interface

interface MACRO_TYPE_MAT_NAME
  module procedure MACRO_NEW_MAT_STRINGS
  module procedure MACRO_NEW_MAT_STRINGARRAY
end interface

#else
! ----------------------------------------------------------------------
! Code to appear in the body, after "contains".
! ----------------------------------------------------------------------

! Check that contents_ is allocated.
impure elemental subroutine MACRO_CHECK_VEC(this)
  implicit none
  
  class(MACRO_TYPE_VEC_NAME), intent(in) :: this
  
  if (.not.allocated(this%contents_)) then
    call print_line(CODE_ERROR//': Trying to use the contents of a vector &
       &before it has been allocated.')
    call err()
  endif
end subroutine

impure elemental subroutine MACRO_CHECK_MAT(this)
  implicit none
  
  class(MACRO_TYPE_MAT_NAME), intent(in) :: this
  
  if (.not.allocated(this%contents_)) then
    call print_line(CODE_ERROR//': Trying to use the contents of a matrix &
       &before it has been allocated.')
    call err()
  endif
end subroutine

! Return an element, row or column.
impure elemental function MACRO_ELEMENT_VEC(this,i) result(output)
  implicit none
  
  class(MACRO_TYPE_VEC_NAME), intent(in) :: this
  integer,                    intent(in) :: i
  MACRO_TYPE_DECLARATION                 :: output
  
  if (i>0 .and. i<=size(this)) then
    output = this%contents_(i)
  else
    call print_line(CODE_ERROR//': Trying to access an element outside of &
       &the vector.')
    call err()
  endif
end function

impure elemental function MACRO_ELEMENT_MAT(this,i,j) result(output)
  implicit none
  
  class(MACRO_TYPE_MAT_NAME), intent(in) :: this
  integer,                    intent(in) :: i
  integer,                    intent(in) :: j
  MACRO_TYPE_DECLARATION                 :: output
  
  if (i>0 .and. i<=size(this,1) .and. j>0 .and. j<=size(this,2)) then
    output = this%contents_(i,j)
  else
    call print_line(CODE_ERROR//': Trying to access an element outside of &
       &the matrix.')
    call err()
  endif
end function

impure elemental function MACRO_ROW_MAT(this,i) result(output)
  implicit none
  
  class(MACRO_TYPE_MAT_NAME), intent(in) :: this
  integer,                    intent(in) :: i
  type(MACRO_TYPE_VEC_NAME)              :: output
  
  if (i>0 .and. i<=size(this,1)) then
    output = vec(this%contents_(i,:))
  else
    call print_line(CODE_ERROR//': Trying to access an element outside of &
       &the matrix.')
    call err()
  endif
end function

impure elemental function MACRO_COLUMN_MAT(this,j) result(output)
  implicit none
  
  class(MACRO_TYPE_MAT_NAME), intent(in) :: this
  integer,                    intent(in) :: j
  type(MACRO_TYPE_VEC_NAME)              :: output
  
  if (j>0 .and. j<=size(this,2)) then
    output = vec(this%contents_(:,j))
  else
    call print_line(CODE_ERROR//': Trying to access an element outside of &
       &the matrix.')
    call err()
  endif
end function

! Convert to an array of the base type.
function MACRO_CONVERT_VEC(input) result(output)
  implicit none
  
  type(MACRO_TYPE_VEC_NAME), intent(in) :: input
  MACRO_TYPE_DECLARATION, allocatable   :: output(:)
  
  call input%check()
  
  output = input%contents_
end function

function MACRO_CONVERT_MAT(input) result(output)
  implicit none
  
  type(MACRO_TYPE_MAT_NAME), intent(in) :: input
  MACRO_TYPE_DECLARATION, allocatable   :: output(:,:)
  
  call input%check()
  
  output = input%contents_
end function

! Convert from an array of the base type.
function MACRO_VEC_VEC(input) result(output)
  implicit none
  
  MACRO_TYPE_DECLARATION, intent(in) :: input(:)
  type(MACRO_TYPE_VEC_NAME)          :: output
  
  output%contents_ = input
end function

function MACRO_MAT_MAT(input) result(output)
  implicit none
  
  MACRO_TYPE_DECLARATION, intent(in) :: input(:,:)
  type(MACRO_TYPE_MAT_NAME)          :: output
  
  output%contents_ = input
end function

function MACRO_MAT_MAT_SHAPE(input,m,n) result(output)
  implicit none
  
  MACRO_TYPE_DECLARATION, intent(in) :: input(:)
  integer,                intent(in) :: m
  integer,                intent(in) :: n
  type(MACRO_TYPE_MAT_NAME)          :: output
  
  output%contents_ = transpose(reshape(input, [m,n]))
end function

! Construct a matrix from vectors.
function MACRO_ROW_MATRIX_VECS(input) result(output)
  implicit none
  
  type(MACRO_TYPE_VEC_NAME), intent(in) :: input(:)
  MACRO_TYPE_DECLARATION, allocatable   :: output(:,:)
  
  integer :: i,ialloc
  
  if (size(input)==0) then
    call print_line(CODE_ERROR//': Trying to make row matrix from empty &
       &array.')
    call err()
  endif
  
  do i=2,size(input)
    if (size(input(i))/=size(input(1))) then
      call print_line(CODE_ERROR//': Trying to make row matrix from &
         &inconsistent vectors.')
      call err()
    endif
  enddo
  
  allocate(output(size(input), size(input(1))), stat=ialloc); call err(ialloc)
  do i=1,size(input)
    output(i,:) = input(i)%contents_
  enddo
end function

function MACRO_COLUMN_MATRIX_VECS(input) result(output)
  implicit none
  
  type(MACRO_TYPE_VEC_NAME), intent(in) :: input(:)
  MACRO_TYPE_DECLARATION, allocatable   :: output(:,:)
  
  integer :: i,ialloc
  
  if (size(input)==0) then
    call print_line(CODE_ERROR//': Trying to make column matrix from empty &
       &array.')
    call err()
  endif
  
  do i=2,size(input)
    if (size(input(i))/=size(input(1))) then
      call print_line(CODE_ERROR//': Trying to make column matrix from &
         &inconsistent vectors.')
      call err()
    endif
  enddo
  
  allocate(output(size(input(1)), size(input)), stat=ialloc); call err(ialloc)
  do i=1,size(input)
    output(:,i) = input(i)%contents_
  enddo
end function

! The size() function.
function MACRO_SIZE_VEC(input) result(output)
  implicit none
  
  type(MACRO_TYPE_VEC_NAME), intent(in) :: input
  integer                               :: output
  
  call input%check()
  
  output = size(input%contents_)
end function

function MACRO_SIZE_MAT(input,dim) result(output)
  implicit none
  
  type(MACRO_TYPE_MAT_NAME), intent(in) :: input
  integer,                   intent(in) :: dim
  integer                               :: output
  
  call input%check()
  
  output = size(input%contents_, dim)
end function

! Trace.
function MACRO_TRACE_MAT(input) result(output)
  implicit none
  
  type(MACRO_TYPE_MAT_NAME), intent(in) :: input
  MACRO_TYPE_DECLARATION                :: output
  
  integer :: i
  
  if (size(input%contents_,1)/=size(input%contents_,2)) then
    call print_line(CODE_ERROR//': Trying to take the trace of a non-square &
       &matrix.')
    call err()
  endif
  
  output = 0
  do i=1,size(input%contents_,1)
    output = output + input%contents_(i,i)
  enddo
end function

! Sum an array of vectors or matrices.
function MACRO_SUM_VECS(input) result(output)
  implicit none
  
  type(MACRO_TYPE_VEC_NAME), intent(in) :: input(:)
  type(MACRO_TYPE_VEC_NAME)             :: output
  
  integer :: i
  
  if (size(input)==0) then
    call print_line(CODE_ERROR//': Trying to take the sum of an empty array.')
    call err()
  endif
  
  call input(1)%check()
  output%contents_ = input(1)%contents_
  do i=2,size(input)
    if (size(input(i))/=size(input(1))) then
      call print_line(CODE_ERROR//': Trying to sum vectors of different &
         &lengths.')
    endif
    output%contents_ = output%contents_ + input(i)%contents_
  enddo
end function

function MACRO_SUM_MATS(input) result(output)
  implicit none
  
  type(MACRO_TYPE_MAT_NAME), intent(in) :: input(:)
  type(MACRO_TYPE_MAT_NAME)             :: output
  
  integer :: i
  
  if (size(input)==0) then
    call print_line(CODE_ERROR//': Trying to take the sum of an empty array.')
    call err()
  endif
  
  call input(1)%check()
  output%contents_ = input(1)%contents_
  do i=2,size(input)
    if ( size(input(i),1)/=size(input(1),1) .or. &
       & size(input(i),2)/=size(input(1),2)      ) then
      call print_line(CODE_ERROR//': Trying to sum vectors of different &
         &lengths.')
    endif
    output%contents_ = output%contents_ + input(i)%contents_
  enddo
end function

! The negative.
impure elemental function MACRO_NEGATIVE_VEC(input) result(output)
  implicit none
  
  type(MACRO_TYPE_VEC_NAME), intent(in) :: input
  type(MACRO_TYPE_VEC_NAME)             :: output
  
  call input%check()
  output%contents_ = -input%contents_
end function

impure elemental function MACRO_NEGATIVE_MAT(input) result(output)
  implicit none
  
  type(MACRO_TYPE_MAT_NAME), intent(in) :: input
  type(MACRO_TYPE_MAT_NAME)             :: output
  
  call input%check()
  output%contents_ = -input%contents_
end function

! I/O.
subroutine MACRO_READ_VEC(this,input)
  implicit none
  
  class(MACRO_TYPE_VEC_NAME), intent(out) :: this
  type(String),               intent(in)  :: input
  
  select type(this); type is(MACRO_TYPE_VEC_NAME)
    this = vec(MACRO_TYPE_CONVERSION(tokens(input)))
  class default
    call err()
  end select
end subroutine

function MACRO_WRITE_VEC(this) result(output)
  implicit none
  
  class(MACRO_TYPE_VEC_NAME), intent(in) :: this
  type(String)                           :: output
  
  select type(this); type is(MACRO_TYPE_VEC_NAME)
    call this%check()
    output = join(MACRO_TYPE_CONVERSION(this))
  class default
    call err()
  end select
end function

impure elemental function MACRO_NEW_VEC_STRING(input) result(this)
  implicit none
  
  type(String), intent(in)  :: input
  type(MACRO_TYPE_VEC_NAME) :: this
  
  call this%read(input)
end function

subroutine MACRO_READ_MAT(this,input)
  implicit none
  
  class(MACRO_TYPE_MAT_NAME), intent(out) :: this
  type(String),               intent(in)  :: input(:)
  
  type(String), allocatable :: line(:)
  
  integer :: m,n
  
  integer :: i,ialloc
  
  select type(this); type is(MACRO_TYPE_MAT_NAME)
    if (size(input)==0) then
      allocate(this%contents_(0,0), stat=ialloc); call err(ialloc)
    else
      m = size(input)
      n = size(tokens(input(1)))
      allocate(this%contents_(m,n), stat=ialloc); call err(ialloc)
      do i=1,size(input)
        line = tokens(input(i))
        if (size(line)/=n) then
          call print_line(ERROR//': Reading matrix: rows of different &
             &lengths.')
          call err()
        endif
        this%contents_(i,:) = MACRO_TYPE_CONVERSION(line)
      enddo
    endif
  class default
    call err()
  end select
end subroutine

function MACRO_WRITE_MAT(this) result(output)
  implicit none
  
  Class(MACRO_TYPE_MAT_NAME), intent(in) :: this
  type(String), allocatable              :: output(:)
  
  integer :: i,ialloc
  
  call this%check()
  
  select type(this); type is(MACRO_TYPE_MAT_NAME)
    allocate(output(size(this%contents_,1)), stat=ialloc); call err(ialloc)
    do i=1,size(this%contents_,1)
      output(i) = join(this%contents_(i,:))
    enddo
  class default
    call err()
  end select
end function

function MACRO_NEW_MAT_STRINGS(input) result(this)
  implicit none
  
  type(String), intent(in)  :: input(:)
  type(MACRO_TYPE_MAT_NAME) :: this
  
  call this%read(input)
end function

impure elemental function MACRO_NEW_MAT_STRINGARRAY(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(MACRO_TYPE_MAT_NAME)     :: this
  
  this = MACRO_TYPE_MAT_NAME(str(input))
end function

#endif

! Undefine function names.
#undef MACRO_CHECK_VEC
#undef MACRO_CHECK_MAT
#undef MACRO_ELEMENT_VEC
#undef MACRO_ELEMENT_MAT
#undef MACRO_ROW_MAT
#undef MACRO_COLUMN_MAT

#undef MACRO_CONVERT_VEC
#undef MACRO_CONVERT_MAT
#undef MACRO_VEC_VEC
#undef MACRO_MAT_MAT
#undef MACRO_MAT_MAT_SHAPE
#undef MACRO_ROW_MATRIX_VECS
#undef MACRO_COLUMN_MATRIX_VECS
#undef MACRO_SIZE_VEC
#undef MACRO_SIZE_MAT
#undef MACRO_TRACE_MAT
#undef MACRO_SUM_VECS
#undef MACRO_SUM_MATS
#undef MACRO_NEGATIVE_VEC
#undef MACRO_NEGATIVE_MAT

#undef MACRO_READ_VEC
#undef MACRO_READ_MAT
#undef MACRO_WRITE_VEC
#undef MACRO_WRITE_MAT
#undef MACRO_NEW_VEC_STRING
#undef MACRO_NEW_MAT_STRINGS
#undef MACRO_NEW_MAT_STRINGARRAY
