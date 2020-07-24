! ======================================================================
! Templates for binary operations involving linear algebra types.
! When this file is included, the macros MACRO_TYPE1_NAME and MACRO_TYPE2_NAME
!    must be defined.
! ======================================================================
! N.B. for each set of types, this file should be included twice:
!    - Before "contains", with MACRO_BODY not defined.
!    - After "contains", with MACRO_BODY defined.

! Include useful macros.
#include "macros.fpp"

! Define MACRO_TYPE3_NAME from MACRO_TYPE1_NAME and MACRO_TYPE2_NAME.
! Define MACRO_TYPE*_DECLARATION, MACRO_TYPE*_VEC_NAME
!    and MACRO_TYPE*_MAT_NAME for * in [1,2,3].
#include "binary_type_declaration.fpp"

! Define function names.
#define MACRO_DOT_SPARSEMAT_VEC MACRO_CAT4(dot_Sparse,MACRO_TYPE1_MAT_NAME,_,MACRO_TYPE2_VEC_NAME)
#define MACRO_DOT_VEC_SPARSEMAT MACRO_CAT4(dot_,MACRO_TYPE1_VEC_NAME,_Sparse,MACRO_TYPE2_MAT_NAME)
#define MACRO_DOT_SPARSEMAT_MAT MACRO_CAT4(dot_Sparse,MACRO_TYPE1_MAT_NAME,_,MACRO_TYPE2_MAT_NAME)
#define MACRO_DOT_MAT_SPARSEMAT MACRO_CAT4(dot_,MACRO_TYPE1_MAT_NAME,_Sparse,MACRO_TYPE2_MAT_NAME)

! Define sparse type names.
#define MACRO_TYPE1_SPARSEMAT_NAME MACRO_CAT2(Sparse,MACRO_TYPE1_MAT_NAME)
#define MACRO_TYPE2_SPARSEMAT_NAME MACRO_CAT2(Sparse,MACRO_TYPE2_MAT_NAME)
#define MACRO_TYPE3_SPARSEMAT_NAME MACRO_CAT2(Sparse,MACRO_TYPE3_MAT_NAME)

#ifndef MACRO_BODY
! ----------------------------------------------------------------------
! Code to appear in the header, before "contains".
! ----------------------------------------------------------------------
interface operator(*)
  module procedure MACRO_DOT_SPARSEMAT_VEC
  module procedure MACRO_DOT_VEC_SPARSEMAT
  module procedure MACRO_DOT_SPARSEMAT_MAT
  module procedure MACRO_DOT_MAT_SPARSEMAT
end interface

#else
! ----------------------------------------------------------------------
! Code to appear in the body, after "contains".
! ----------------------------------------------------------------------

impure elemental function MACRO_DOT_SPARSEMAT_VEC(a,b) result(output)
  implicit none
  
  type(MACRO_TYPE1_SPARSEMAT_NAME), intent(in) :: a
  type(MACRO_TYPE2_VEC_NAME),       intent(in) :: b
  type(MACRO_TYPE3_VEC_NAME)                   :: output
  
  MACRO_TYPE2_DECLARATION, allocatable :: b_elements(:)
  MACRO_TYPE3_DECLARATION, allocatable :: elements(:)
  
  integer :: i,ialloc
  
  b_elements = MACRO_TYPE2_CONVERT(b)
  
  allocate(elements(a%no_rows_), stat=ialloc); call err(ialloc)
  elements = 0
  do i=1,a%no_elements_
    elements(a%elements_(i)%i) = elements(a%elements_(i)%i) &
                             & + a%elements_(i)%element     &
                             & * b_elements(a%elements_(i)%j)
  enddo
  output = vec(elements)
end function

impure elemental function MACRO_DOT_VEC_SPARSEMAT(a,b) result(output)
  implicit none
  
  type(MACRO_TYPE1_VEC_NAME),       intent(in) :: a
  type(MACRO_TYPE2_SPARSEMAT_NAME), intent(in) :: b
  type(MACRO_TYPE3_VEC_NAME)                   :: output
  
  MACRO_TYPE1_DECLARATION, allocatable :: a_elements(:)
  MACRO_TYPE3_DECLARATION, allocatable :: elements(:)
  
  integer :: i,ialloc
  
  a_elements = MACRO_TYPE1_CONVERT(a)
  
  allocate(elements(b%no_cols_), stat=ialloc); call err(ialloc)
  elements = 0
  do i=1,b%no_elements_
    elements(b%elements_(i)%j) = elements(b%elements_(i)%j)   &
                             & + a_elements(b%elements_(i)%i) &
                             & * b%elements_(i)%element
  enddo
  output = vec(elements)
end function

impure elemental function MACRO_DOT_SPARSEMAT_MAT(a,b) result(output)
  implicit none
  
  type(MACRO_TYPE1_SPARSEMAT_NAME), intent(in) :: a
  type(MACRO_TYPE2_MAT_NAME),       intent(in) :: b
  type(MACRO_TYPE3_MAT_NAME)                   :: output
  
  MACRO_TYPE2_DECLARATION, allocatable :: b_elements(:,:)
  MACRO_TYPE3_DECLARATION, allocatable :: elements(:,:)
  
  integer :: i,ialloc
  
  b_elements = MACRO_TYPE2_CONVERT(b)
  
  allocate(elements(a%no_rows_, size(b,2)), stat=ialloc); call err(ialloc)
  elements = 0
  do i=1,a%no_elements_
    elements(a%elements_(i)%i,:) = elements(a%elements_(i)%i,:) &
                               & + a%elements_(i)%element &
                               & * b_elements(a%elements_(i)%j,:)
  enddo
  output = mat(elements)
end function

impure elemental function MACRO_DOT_MAT_SPARSEMAT(a,b) result(output)
  implicit none
  
  type(MACRO_TYPE1_MAT_NAME),       intent(in) :: a
  type(MACRO_TYPE2_SPARSEMAT_NAME), intent(in) :: b
  type(MACRO_TYPE3_MAT_NAME)                   :: output
  
  MACRO_TYPE1_DECLARATION, allocatable :: a_elements(:,:)
  MACRO_TYPE3_DECLARATION, allocatable :: elements(:,:)
  
  integer :: i,ialloc
  
  a_elements = MACRO_TYPE1_CONVERT(a)
  
  allocate(elements(size(a,1), b%no_cols_), stat=ialloc); call err(ialloc)
  elements = 0
  do i=1,b%no_elements_
    elements(:,b%elements_(i)%j) = elements(:,b%elements_(i)%j)   &
                               & + a_elements(:,b%elements_(i)%i) &
                               & * b%elements_(i)%element
  enddo
  output = mat(elements)
end function
#endif
