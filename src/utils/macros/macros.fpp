! ======================================================================
! Various preprocessor macros.
! ======================================================================
#ifndef MACROS_FPP
#define MACROS_FPP

! The MACRO_CAT macro, which expands and concatenates two macros.
! If M1=a and M2=b then MACRO_CAT(M1,M2)=ab.
#if defined __GFORTRAN__ || defined NAGFOR
#define MACRO_PASTE(a) a
#define MACRO_CAT(a,b) MACRO_PASTE(a)b
#else
#define MACRO_PASTE(a,b) a ## b
#define MACRO_CAT(a,b) MACRO_PASTE(a,b)
#endif

! Define concatenation macros for multiple arguments.
#define MACRO_CAT2(a,b) MACRO_CAT(a,b)
#define MACRO_CAT3(a,b,c) MACRO_CAT(a,MACRO_CAT(b,c))
#define MACRO_CAT4(a,b,c,d) MACRO_CAT(a,MACRO_CAT3(b,c,d))

#endif
