! ======================================================================
! Defines various macros based on MACRO_TYPE1_NAME and MACRO_TYPE2_NAME.
!
! Defines MACRO_TYPE3_NAME from MACRO_TYPE1_NAME and MACRO_TYPE2_NAME,
!    such that the following function is valid:
! 
! function add(x,y) result(z)
!    MACRO_TYPE1_NAME :: x
!    MACRO_TYPE2_NAME :: y
!    MACRO_TYPE3_NAME :: z
!    z = x+y
! end function
!
! Also defines, for * in [1,2,3]:
!
! MACRO_TYPE*_NAME        | integer    real        complex
! ------------------------+-------------------------------------
! MACRO_TYPE*_DECLARATION | integer    real(dp)    complex(dp)
! MACRO_TYPE*_VEC_NAME    | IntVector  RealVector  ComplexVector
! MACRO_TYPE*_MAT_NAME    | IntMatrix  RealMatrix  ComplexMatrix
! ======================================================================

! Undefine macros if they are already defined
#ifdef MACRO_TYPE1_VEC_NAME
#undef MACRO_TYPE1_VEC_NAME
#endif

#ifdef MACRO_TYPE1_MAT_NAME
#undef MACRO_TYPE1_MAT_NAME
#endif

#ifdef MACRO_TYPE2_VEC_NAME
#undef MACRO_TYPE2_VEC_NAME
#endif

#ifdef MACRO_TYPE2_MAT_NAME
#undef MACRO_TYPE2_MAT_NAME
#endif

#ifdef MACRO_TYPE3_NAME
#undef MACRO_TYPE3_NAME
#endif

#ifdef MACRO_TYPE3_VEC_NAME
#undef MACRO_TYPE3_VEC_NAME
#endif

#ifdef MACRO_TYPE3_MAT_NAME
#undef MACRO_TYPE3_MAT_NAME
#endif

! Define MACRO_TYPE3_NAME from MACRO_TYPE1_NAME and MACRO_TYPE2_NAME.
#define integer 1
#define real 3
#define complex 4

#if MACRO_TYPE1_NAME>MACRO_TYPE2_NAME
#define MACRO_TYPE3_NAME MACRO_TYPE1_NAME
#else
#define MACRO_TYPE3_NAME MACRO_TYPE2_NAME
#endif

#undef integer
#undef real
#undef complex

! Identify each MACRO_TYPE*_NAME, and define the relevant MACRO_*_*.
#define integer 1
#if MACRO_TYPE1_NAME==1
#define MACRO_1_INTEGER
#endif
#if MACRO_TYPE2_NAME==1
#define MACRO_2_INTEGER
#endif
#if MACRO_TYPE3_NAME==1
#define MACRO_3_INTEGER
#endif
#undef integer

#define real 1
#if MACRO_TYPE1_NAME==1
#define MACRO_1_REAL
#endif
#if MACRO_TYPE2_NAME==1
#define MACRO_2_REAL
#endif
#if MACRO_TYPE3_NAME==1
#define MACRO_3_REAL
#endif
#undef real

#define complex 1
#if MACRO_TYPE1_NAME==1
#define MACRO_1_COMPLEX
#endif
#if MACRO_TYPE2_NAME==1
#define MACRO_2_COMPLEX
#endif
#if MACRO_TYPE3_NAME==1
#define MACRO_3_COMPLEX
#endif
#undef complex

! Define MACRO_TYPE1_VEC_NAME and MACRO_TYPE2_VEC_NAME.
! Also undefine the relevant MACRO_1_*.
#if defined MACRO_1_INTEGER

#define MACRO_TYPE1_VEC_NAME IntVector
#define MACRO_TYPE1_MAT_NAME IntMatrix
#undef MACRO_1_INTEGER

#elif defined MACRO_1_REAL

#define MACRO_TYPE1_VEC_NAME RealVector
#define MACRO_TYPE1_MAT_NAME RealMatrix
#undef MACRO_1_REAL

#elif defined MACRO_1_COMPLEX

#define MACRO_TYPE1_VEC_NAME ComplexVector
#define MACRO_TYPE1_MAT_NAME ComplexMatrix
#undef MACRO_1_COMPLEX

#else

#error "MACRO_TYPE1_NAME is not integer, real or complex."

#endif

! Define MACRO_TYPE2_VEC_NAME and MACRO_TYPE2_VEC_NAME.
! Also undefine the relevant MACRO_2_*.
#if defined MACRO_2_INTEGER

#define MACRO_TYPE2_VEC_NAME IntVector
#define MACRO_TYPE2_MAT_NAME IntMatrix
#undef MACRO_2_INTEGER

#elif defined MACRO_2_REAL

#define MACRO_TYPE2_VEC_NAME RealVector
#define MACRO_TYPE2_MAT_NAME RealMatrix
#undef MACRO_2_REAL

#elif defined MACRO_2_COMPLEX

#define MACRO_TYPE2_VEC_NAME ComplexVector
#define MACRO_TYPE2_MAT_NAME ComplexMatrix
#undef MACRO_2_COMPLEX

#else

#error "MACRO_TYPE2_NAME is not integer, real or complex."

#endif

! Define MACRO_TYPE3_VEC_NAME and MACRO_TYPE3_VEC_NAME.
! Also undefine the relevant MACRO_3_*.
#if defined MACRO_3_INTEGER

#define MACRO_TYPE3_VEC_NAME IntVector
#define MACRO_TYPE3_MAT_NAME IntMatrix
#undef MACRO_3_INTEGER

#elif defined MACRO_3_REAL

#define MACRO_TYPE3_VEC_NAME RealVector
#define MACRO_TYPE3_MAT_NAME RealMatrix
#undef MACRO_3_REAL

#elif defined MACRO_3_COMPLEX

#define MACRO_TYPE3_VEC_NAME ComplexVector
#define MACRO_TYPE3_MAT_NAME ComplexMatrix
#undef MACRO_3_COMPLEX

#else

#error "MACRO_TYPE3_NAME is not integer, real or complex."

#endif

! Define MACRO_TYPE*_DECLARATION
#include "type1_declaration.fpp"
#include "type2_declaration.fpp"
#include "type3_declaration.fpp"
