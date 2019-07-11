! ======================================================================
! Defines various macros based on MACRO_TYPE_NAME.
!
! MACRO_TYPE_NAME        | integer    real        complex
! -----------------------+-------------------------------------
! MACRO_TYPE_DECLARATION | integer    real(dp)    complex(dp)
! MACRO_TYPE_VEC_NAME    | IntVector  RealVector  ComplexVector
! MACRO_TYPE_MAT_NAME    | IntMatrix  RealMatrix  ComplexMatrix
! ======================================================================

! Define MACRO_TYPE_DECLARATION
#include "type_declaration.fpp"

! Undefine MACRO_TYPE_VEC_NAME and MACRO_TYPE_MAT_NAME
!    if they are already defined
#ifdef MACRO_TYPE_VEC_NAME
#undef MACRO_TYPE_VEC_NAME
#endif

#ifdef MACRO_TYPE_MAT_NAME
#undef MACRO_TYPE_MAT_NAME
#endif

! Identify MACRO_TYPE_NAME from integer, real and complex.
! Define the relevant MACRO_*.
#define integer 1
#if MACRO_TYPE_NAME==1
#define MACRO_INTEGER
#endif
#undef integer

#define real 1
#if MACRO_TYPE_NAME==1
#define MACRO_REAL
#endif
#undef real

#define complex 1
#if MACRO_TYPE_NAME==1
#define MACRO_COMPLEX
#endif
#undef complex

! Define MACRO_TYPE_VEC_NAME and MACRO_TYPE2_VEC_NAME.
! Also undefine the relevant MACRO_*.
#if defined MACRO_INTEGER

#define MACRO_TYPE_VEC_NAME IntVector
#define MACRO_TYPE_MAT_NAME IntMatrix
#undef MACRO_INTEGER

#elif defined MACRO_REAL

#define MACRO_TYPE_VEC_NAME RealVector
#define MACRO_TYPE_MAT_NAME RealMatrix
#undef MACRO_REAL

#elif defined MACRO_COMPLEX

#define MACRO_TYPE_VEC_NAME ComplexVector
#define MACRO_TYPE_MAT_NAME ComplexMatrix
#undef MACRO_COMPLEX

#else

#error "MACRO_TYPE_NAME is not integer, real or complex."

#endif
