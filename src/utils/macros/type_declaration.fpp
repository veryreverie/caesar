! ======================================================================
! Defines MACRO_TYPE_DECLARATION based on MACRO_TYPE_NAME.
!
! MACRO_TYPE_NAME | MACRO_TYPE_DECLARATION
! ----------------+----------------------------
! character       | character(*)
! integer         | integer
! real            | real(dp)
! complex         | complex(dp)
!                 |
! String          | type(String)
!
! Other user-defined typenames follow the same pattern as String.
! ======================================================================

! Check that MACRO_TYPE_NAME is set.
#ifndef MACRO_TYPE_NAME
#error "Error: MACRO_TYPE_NAME is not defined."
#endif

! Undefine MACRO_TYPE_DECLARATION if it is already defined.
#ifdef MACRO_TYPE_DECLARATION
#undef MACRO_TYPE_DECLARATION
#endif

! Test to see if MACRO_TYPE_NAME is an intrinsic type,
!    and define the relevant MACRO_* if so.
#define character 1
#if MACRO_TYPE_NAME==1
#define MACRO_CHARACTER
#endif
#undef character

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

! Use the MACRO_* macros to define MACRO_TYPE_DECLARATION.
! Also undefines whichever MACRO_* macro is defined.
#if defined MACRO_CHARACTER
#define MACRO_TYPE_DECLARATION character(*)
#undef MACRO_CHARACTER
#elif defined MACRO_INTEGER
#define MACRO_TYPE_DECLARATION integer
#undef MACRO_INTEGER
#elif defined MACRO_REAL
#define MACRO_TYPE_DECLARATION real(dp)
#undef MACRO_REAL
#elif defined MACRO_COMPLEX
#define MACRO_TYPE_DECLARATION complex(dp)
#undef MACRO_COMPLEX
#else
#define MACRO_TYPE_DECLARATION type(MACRO_TYPE_NAME)
#endif
