! ======================================================================
! Templates for binary operations involving linear algebra types.
! When this file is included, the macros MACRO_TYPE1_NAME and MACRO_TYPE2_NAME
!    must be defined.
! ======================================================================
! N.B. for each set of types, this file should be included twice:
!    - Before "contains", with MACRO_BODY not defined.
!    - After "contains", with MACRO_BODY defined.
!
! Functions:
!    - +,- : addition of vectors / matrices.
!    - *,/ : multiplication and division of vectors / matrices by scalars.
!    - *   : dot-products of vectors.
!    - *   : matrix multiplication between matrices and matrices / vectors.
!    - cross_product.
!    - outer_product.
!    - commutator.
!    - anticommutator.
!  
! ======================================================================

! Include useful macros.
#include "macros.fpp"

! Define MACRO_TYPE3_NAME from MACRO_TYPE1_NAME and MACRO_TYPE2_NAME.
! Define MACRO_TYPE*_DECLARATION, MACRO_TYPE*_VEC_NAME
!    and MACRO_TYPE*_MAT_NAME for * in [1,2,3].
#include "binary_type_declaration.fpp"

! Define function names.
#define MACRO_ADD_VEC        MACRO_CAT4(add_,MACRO_TYPE1_VEC_NAME,_,MACRO_TYPE2_VEC_NAME)
#define MACRO_ADD_MAT        MACRO_CAT4(add_,MACRO_TYPE1_MAT_NAME,_,MACRO_TYPE2_MAT_NAME)
#define MACRO_SUB_VEC        MACRO_CAT4(subtract_,MACRO_TYPE1_VEC_NAME,_,MACRO_TYPE2_VEC_NAME)
#define MACRO_SUB_MAT        MACRO_CAT4(subtract_,MACRO_TYPE1_MAT_NAME,_,MACRO_TYPE2_MAT_NAME)
#define MACRO_MUL_VEC_SCALAR MACRO_CAT4(multiply_,MACRO_TYPE1_VEC_NAME,_,MACRO_TYPE2_NAME)
#define MACRO_MUL_SCALAR_VEC MACRO_CAT4(multiply_,MACRO_TYPE1_NAME,_,MACRO_TYPE2_VEC_NAME)
#define MACRO_MUL_MAT_SCALAR MACRO_CAT4(multiply_,MACRO_TYPE1_MAT_NAME,_,MACRO_TYPE2_NAME)
#define MACRO_MUL_SCALAR_MAT MACRO_CAT4(multiply_,MACRO_TYPE1_NAME,_,MACRO_TYPE2_MAT_NAME)
#define MACRO_DIV_VEC_SCALAR MACRO_CAT4(divide_,MACRO_TYPE1_VEC_NAME,_,MACRO_TYPE2_NAME)
#define MACRO_DIV_MAT_SCALAR MACRO_CAT4(divide_,MACRO_TYPE1_MAT_NAME,_,MACRO_TYPE2_NAME)
#define MACRO_DOT_VEC_VEC    MACRO_CAT4(dot_,MACRO_TYPE1_VEC_NAME,_,MACRO_TYPE2_VEC_NAME)
#define MACRO_DOT_VEC_MAT    MACRO_CAT4(dot_,MACRO_TYPE1_VEC_NAME,_,MACRO_TYPE2_MAT_NAME)
#define MACRO_DOT_MAT_VEC    MACRO_CAT4(dot_,MACRO_TYPE1_MAT_NAME,_,MACRO_TYPE2_VEC_NAME)
#define MACRO_DOT_MAT_MAT    MACRO_CAT4(dot_,MACRO_TYPE1_MAT_NAME,_,MACRO_TYPE2_MAT_NAME)
#define MACRO_CROSS_PRODUCT  MACRO_CAT4(cross_product_,MACRO_TYPE1_VEC_NAME,_,MACRO_TYPE2_VEC_NAME)
#define MACRO_OUTER_PRODUCT  MACRO_CAT4(outer_product_,MACRO_TYPE1_VEC_NAME,_,MACRO_TYPE2_VEC_NAME)
#define MACRO_COMMUTATOR     MACRO_CAT4(commutator_,MACRO_TYPE1_MAT_NAME,_,MACRO_TYPE2_MAT_NAME)
#define MACRO_ANTICOMMUTATOR MACRO_CAT4(anticommutator_,MACRO_TYPE1_MAT_NAME,_,MACRO_TYPE2_MAT_NAME)

#ifndef MACRO_BODY
! ----------------------------------------------------------------------
! Code to appear in the header, before "contains".
! ----------------------------------------------------------------------
interface operator(+)
  module procedure MACRO_ADD_VEC
  module procedure MACRO_ADD_MAT
end interface

interface operator(-)
  module procedure MACRO_SUB_VEC
  module procedure MACRO_SUB_MAT
end interface

interface operator(*)
  module procedure MACRO_MUL_VEC_SCALAR
  module procedure MACRO_MUL_SCALAR_VEC
  module procedure MACRO_MUL_MAT_SCALAR
  module procedure MACRO_MUL_SCALAR_MAT
  
  module procedure MACRO_DOT_VEC_VEC
  module procedure MACRO_DOT_VEC_MAT
  module procedure MACRO_DOT_MAT_VEC
  module procedure MACRO_DOT_MAT_MAT
end interface

interface operator(/)
  module procedure MACRO_DIV_VEC_SCALAR
  module procedure MACRO_DIV_MAT_SCALAR
end interface

interface cross_product
  module procedure MACRO_CROSS_PRODUCT
end interface

interface outer_product
  module procedure MACRO_OUTER_PRODUCT
end interface

interface commutator
  module procedure MACRO_COMMUTATOR
end interface

interface anticommutator
  module procedure MACRO_ANTICOMMUTATOR
end interface

#else
! ----------------------------------------------------------------------
! Code to appear in the body, after "contains".
! ----------------------------------------------------------------------

! Addition and subtraction.
impure elemental function MACRO_ADD_VEC(a,b) result(output)
  implicit none
  
  type(MACRO_TYPE1_VEC_NAME), intent(in) :: a
  type(MACRO_TYPE2_VEC_NAME), intent(in) :: b
  type(MACRO_TYPE3_VEC_NAME)             :: output
  
  if (size(a)/=size(b)) then
    call print_line(CODE_ERROR//': Adding vectors of different lengths.')
    call err()
  endif
  
  output%contents_ = a%contents_ + b%contents_
end function

impure elemental function MACRO_ADD_MAT(a,b) result(output)
  implicit none
  
  type(MACRO_TYPE1_MAT_NAME), intent(in) :: a
  type(MACRO_TYPE2_MAT_NAME), intent(in) :: b
  type(MACRO_TYPE3_MAT_NAME)             :: output
  
  if (size(a,1)/=size(b,1) .or. size(a,2)/=size(b,2)) then
    call print_line(CODE_ERROR//': Adding vectors of different lengths.')
    call err()
  endif
  
  output%contents_ = a%contents_ + b%contents_
end function

impure elemental function MACRO_SUB_VEC(a,b) result(output)
  implicit none
  
  type(MACRO_TYPE1_VEC_NAME), intent(in) :: a
  type(MACRO_TYPE2_VEC_NAME), intent(in) :: b
  type(MACRO_TYPE3_VEC_NAME)             :: output
  
  if (size(a)/=size(b)) then
    call print_line(CODE_ERROR//': Subtracting vectors of different lengths.')
    call err()
  endif
  
  output%contents_ = a%contents_ - b%contents_
end function

impure elemental function MACRO_SUB_MAT(a,b) result(output)
  implicit none
  
  type(MACRO_TYPE1_MAT_NAME), intent(in) :: a
  type(MACRO_TYPE2_MAT_NAME), intent(in) :: b
  type(MACRO_TYPE3_MAT_NAME)             :: output
  
  if (size(a,1)/=size(b,1) .or. size(a,2)/=size(b,2)) then
    call print_line(CODE_ERROR//': Subtracting vectors of different lengths.')
    call err()
  endif
  
  output%contents_ = a%contents_ - b%contents_
end function

! Multiplication and division by scalar.
impure elemental function MACRO_MUL_VEC_SCALAR(a,b) result(output)
  implicit none
  
  type(MACRO_TYPE1_VEC_NAME), intent(in) :: a
  MACRO_TYPE2_DECLARATION,    intent(in) :: b
  type(MACRO_TYPE3_VEC_NAME)             :: output
  
  call a%check()
  
  output%contents_ = a%contents_ * b
end function

impure elemental function MACRO_MUL_SCALAR_VEC(a,b) result(output)
  implicit none
  
  MACRO_TYPE1_DECLARATION,    intent(in) :: a
  type(MACRO_TYPE2_VEC_NAME), intent(in) :: b
  type(MACRO_TYPE3_VEC_NAME)             :: output
  
  call b%check()
  
  output%contents_ = a * b%contents_
end function

impure elemental function MACRO_MUL_MAT_SCALAR(a,b) result(output)
  implicit none
  
  type(MACRO_TYPE1_MAT_NAME), intent(in) :: a
  MACRO_TYPE2_DECLARATION,    intent(in) :: b
  type(MACRO_TYPE3_MAT_NAME)             :: output
  
  call a%check()
  
  output%contents_ = a%contents_ * b
end function

impure elemental function MACRO_MUL_SCALAR_MAT(a,b) result(output)
  implicit none
  
  MACRO_TYPE1_DECLARATION,    intent(in) :: a
  type(MACRO_TYPE2_MAT_NAME), intent(in) :: b
  type(MACRO_TYPE3_MAT_NAME)             :: output
  
  call b%check()
  
  output%contents_ = a * b%contents_
end function

impure elemental function MACRO_DIV_VEC_SCALAR(a,b) result(output)
  implicit none
  
  type(MACRO_TYPE1_VEC_NAME), intent(in) :: a
  MACRO_TYPE2_DECLARATION,    intent(in) :: b
  type(MACRO_TYPE3_VEC_NAME)             :: output
  
  call a%check()
  
  output%contents_ = a%contents_ / b
end function

impure elemental function MACRO_DIV_MAT_SCALAR(a,b) result(output)
  implicit none
  
  type(MACRO_TYPE1_MAT_NAME), intent(in) :: a
  MACRO_TYPE2_DECLARATION,    intent(in) :: b
  type(MACRO_TYPE3_MAT_NAME)             :: output
  
  call a%check()
  
  output%contents_ = a%contents_ / b
end function

! Dot products and matrix multiplication.
impure elemental function MACRO_DOT_VEC_VEC(a,b) result(output)
  implicit none
  
  type(MACRO_TYPE1_VEC_NAME), intent(in) :: a
  type(MACRO_TYPE2_VEC_NAME), intent(in) :: b
  MACRO_TYPE3_DECLARATION                :: output
  
  integer :: i
  
  if (size(a)/=size(b)) then
    call print_line(CODE_ERROR//': Taking dot product of vectors of different &
       &lengths.')
    call err()
  endif
  
  output = 0
  do i=1,size(a)
    output = output + a%contents_(i)*b%contents_(i)
  enddo
end function

impure elemental function MACRO_DOT_VEC_MAT(a,b) result(output)
  implicit none
  
  type(MACRO_TYPE1_VEC_NAME), intent(in) :: a
  type(MACRO_TYPE2_MAT_NAME), intent(in) :: b
  type(MACRO_TYPE3_VEC_NAME)             :: output
  
  if (size(a)/=size(b,1)) then
    call print_line(CODE_ERROR//': Taking dot product of incompatible vector &
       &and matrix.')
    call err()
  endif
  
  output%contents_ = matmul(a%contents_, b%contents_)
end function

impure elemental function MACRO_DOT_MAT_VEC(a,b) result(output)
  implicit none
  
  type(MACRO_TYPE1_MAT_NAME), intent(in) :: a
  type(MACRO_TYPE2_VEC_NAME), intent(in) :: b
  type(MACRO_TYPE3_VEC_NAME)             :: output
  
  if (size(a,2)/=size(b)) then
    call print_line(CODE_ERROR//': Taking dot product of incompatible matrix &
       &and vector.')
    call err()
  endif
  
  output%contents_ = matmul(a%contents_, b%contents_)
end function

impure elemental function MACRO_DOT_MAT_MAT(a,b) result(output)
  implicit none
  
  type(MACRO_TYPE1_MAT_NAME), intent(in) :: a
  type(MACRO_TYPE2_MAT_NAME), intent(in) :: b
  type(MACRO_TYPE3_MAT_NAME)             :: output
  
  if (size(a,2)/=size(b,1)) then
    call print_line(CODE_ERROR//': Taking dot product of incompatible matrix &
       &and matrix.')
    call err()
  endif
  
  output%contents_ = matmul(a%contents_, b%contents_)
end function

! Cross product, outer product, commutator and anti-commutator.
impure elemental function MACRO_CROSS_PRODUCT(a,b) result(output)
  implicit none
  
  type(MACRO_TYPE1_VEC_NAME), intent(in) :: a
  type(MACRO_TYPE2_VEC_NAME), intent(in) :: b
  type(MACRO_TYPE3_VEC_NAME)             :: output
  
  if (size(a)/=3 .or. size(b)/=3) then
    call print_line(CODE_ERROR//': Trying to take a cross product involving &
       &a vector which has other than three components.')
    call err()
  endif
  
  output%contents_ = [                                                &
     & a%contents_(2)*b%contents_(3) - b%contents_(2)*a%contents_(3), &
     & a%contents_(3)*b%contents_(1) - b%contents_(3)*a%contents_(1), &
     & a%contents_(1)*b%contents_(2) - b%contents_(1)*a%contents_(2)  ]
end function

impure elemental function MACRO_OUTER_PRODUCT(a,b) result(output)
  implicit none
  
  type(MACRO_TYPE1_VEC_NAME), intent(in) :: a
  type(MACRO_TYPE2_VEC_NAME), intent(in) :: b
  type(MACRO_TYPE3_MAT_NAME)             :: output
  
  integer :: i,j,ialloc
  
  allocate(output%contents_(size(a),size(b)), stat=ialloc); call err(ialloc)
  do i=1,size(b)
    do j=1,size(a)
      output%contents_(j,i) = a%contents_(j)*b%contents_(i)
    enddo
  enddo
end function

impure elemental function MACRO_COMMUTATOR(a,b) result(output)
  implicit none
  
  type(MACRO_TYPE1_MAT_NAME), intent(in) :: a
  type(MACRO_TYPE2_MAT_NAME), intent(in) :: b
  type(MACRO_TYPE3_MAT_NAME)             :: output
  
  if ( size(a,1)/=size(a,2) .or. &
     & size(a,1)/=size(b,1) .or. &
     & size(a,1)/=size(b,2)      ) then
    call print_line(CODE_ERROR//': Taking the commutator of incompatible &
       &matrices.')
    call err()
  endif
  
  output%contents_ = matmul(a%contents_,b%contents_) &
                 & - matmul(b%contents_,a%contents_)
end function

impure elemental function MACRO_ANTICOMMUTATOR(a,b) result(output)
  implicit none
  
  type(MACRO_TYPE1_MAT_NAME), intent(in) :: a
  type(MACRO_TYPE2_MAT_NAME), intent(in) :: b
  type(MACRO_TYPE3_MAT_NAME)             :: output
  
  if ( size(a,1)/=size(a,2) .or. &
     & size(a,1)/=size(b,1) .or. &
     & size(a,1)/=size(b,2)      ) then
    call print_line(CODE_ERROR//': Taking the anti-commutator of incompatible &
       &matrices.')
    call err()
  endif
  
  output%contents_ = matmul(a%contents_,b%contents_) &
                 & + matmul(b%contents_,a%contents_)
end function

#endif

! Undefine function names.
#undef MACRO_ADD_VEC
#undef MACRO_ADD_MAT
#undef MACRO_SUB_VEC
#undef MACRO_SUB_MAT
#undef MACRO_MUL_VEC_SCALAR
#undef MACRO_MUL_SCALAR_VEC
#undef MACRO_MUL_MAT_SCALAR
#undef MACRO_MUL_SCALAR_MAT
#undef MACRO_DIV_VEC_SCALAR
#undef MACRO_DIV_MAT_SCALAR
#undef MACRO_DOT_VEC_VEC
#undef MACRO_DOT_VEC_MAT
#undef MACRO_DOT_MAT_VEC
#undef MACRO_DOT_MAT_MAT
#undef MACRO_CROSS_PRODUCT
#undef MACRO_OUTER_PRODUCT
#undef MACRO_COMMUTATOR
#undef MACRO_ANTICOMMUTATOR
