! ======================================================================
! Templates for binary operations involving linear algebra types.
! When this file is included, the macros MACRO_TYPE*_NAME for * in [1,2,3]
!    must be defined.
! These types must satisfy x1+x2=x3 if x* is of type MACRO_TYPE*_NAME.
! ======================================================================
! N.B. for each set of types, this file should be included twice:
!    - Before "contains", with MACRO_BODY not defined.
!    - After "contains", with MACRO_BODY defined.

! Include useful macros.
#include "macros.fpp"

! Define MACRO_TYPE*_DECLARATION, MACRO_TYPE*_VEC_NAME
!    and MACRO_TYPE*_MAT_NAME for * in [1,2,3].
#include "type1_algebra_declaration.fpp"
#include "type2_algebra_declaration.fpp"
#include "type3_algebra_declaration.fpp"

! Define function names.
#define MACRO_BLOB_VEC_VEC MACRO_CAT4(blob_,MACRO_TYPE1_VEC_NAME,_,MACRO_TYPE2_VEC_NAME)

#ifndef MACRO_BODY

interface blob
  module procedure MACRO_BLOB_VEC_VEC
end interface

#else

function MACRO_BLOB_VEC_VEC(this,that) result(output)
  implicit none
  
  type(MACRO_TYPE1_VEC_NAME), intent(in) :: this
  type(MACRO_TYPE2_VEC_NAME), intent(in) :: that
  type(MACRO_TYPE3_VEC_NAME)             :: output
  
  output%contents_ = this%contents_+that%contents_
end function

#endif

! Undefine function names.
#undef MACRO_BLOB_VEC_VEC
