! ======================================================================
! Templates for unary operations involving linear algebra types.
! When this file is included, the macro MACRO_TYPE_NAME must be defined.
! ======================================================================
! N.B. for each type, this file should be included twice:
!    - Before "contains", with MACRO_BODY not defined.
!    - After "contains", with MACRO_BODY defined.

! Include useful macros.
#include "macros.fpp"

! Define MACRO_TYPE_DECLARATION, MACRO_TYPE_VEC_NAME
!    and MACRO_TYPE_MAT_NAME.
#include "type_algebra_declaration.fpp"

! Define function names.
#define MACRO_NEW_VEC             MACRO_CAT2(new_,MACRO_TYPE_VEC_NAME)
#define MACRO_NEW_MAT             MACRO_CAT2(new_,MACRO_TYPE_MAT_NAME)
#define MACRO_CHECK_VEC           MACRO_CAT2(check_,MACRO_TYPE_VEC_NAME)
#define MACRO_CHECK_MAT           MACRO_CAT2(check_,MACRO_TYPE_MAT_NAME)
#define MACRO_ELEMENT_VEC         MACRO_CAT2(element_,MACRO_TYPE_VEC_NAME)
#define MACRO_ELEMENT_MAT         MACRO_CAT2(element_,MACRO_TYPE_MAT_NAME)
#define MACRO_ROW_MAT             MACRO_CAT2(row_,MACRO_TYPE_MAT_NAME)
#define MACRO_COLUMN_MAT          MACRO_CAT2(column_,MACRO_TYPE_MAT_NAME)
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

!interface MACRO_TYPE_VEC_NAME
!  module procedure MACRO_NEW_VEC
!end interface

! ----------------------------------------------------------------------
#else
! ----------------------------------------------------------------------
! Code to appear in the body, after "contains".
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
#endif

! Undefine function names.
#undef MACRO_NEW_VEC
#undef MACRO_NEW_MAT
#undef MACRO_CHECK_VEC
#undef MACRO_CHECK_MAT
#undef MACRO_ELEMENT_VEC
#undef MACRO_ELEMENT_MAT
#undef MACRO_ROW_MAT
#undef MACRO_COLUMN_MAT
#undef MACRO_READ_VEC
#undef MACRO_READ_MAT
#undef MACRO_WRITE_VEC
#undef MACRO_WRITE_MAT
#undef MACRO_NEW_VEC_STRING
#undef MACRO_NEW_MAT_STRINGS
#undef MACRO_NEW_MAT_STRINGARRAY
