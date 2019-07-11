! ======================================================================
! Includes unary_algebra for all types,
!    and binary_algebra for all combinations of types.
! ======================================================================

! Include unary algebra for all types.
#ifdef MACRO_TYPE_NAME
#undef MACRO_TYPE_NAME
#endif
#define MACRO_TYPE_NAME integer
#include "unary_algebra.fpp"
#undef MACRO_TYPE_NAME
#define MACRO_TYPE_NAME real
#include "unary_algebra.fpp"
#undef MACRO_TYPE_NAME
#define MACRO_TYPE_NAME complex
#include "unary_algebra.fpp"
#undef MACRO_TYPE_NAME

! Include binary_algebra for all combinations of types.
#ifdef MACRO_TYPE1_NAME
#undef MACRO_TYPE1_NAME
#endif
#ifdef MACRO_TYPE2_NAME
#undef MACRO_TYPE2_NAME
#endif
#ifdef MACRO_TYPE3_NAME
#undef MACRO_TYPE3_NAME
#endif
#define MACRO_TYPE1_NAME integer
#define MACRO_TYPE2_NAME integer
#define MACRO_TYPE3_NAME integer
#include "binary_algebra.fpp"
#undef MACRO_TYPE1_NAME
#undef MACRO_TYPE2_NAME
#undef MACRO_TYPE3_NAME
#define MACRO_TYPE1_NAME integer
#define MACRO_TYPE2_NAME real
#define MACRO_TYPE3_NAME real
#include "binary_algebra.fpp"
#undef MACRO_TYPE1_NAME
#undef MACRO_TYPE2_NAME
#undef MACRO_TYPE3_NAME
#define MACRO_TYPE1_NAME integer
#define MACRO_TYPE2_NAME complex
#define MACRO_TYPE3_NAME complex
#include "binary_algebra.fpp"
#undef MACRO_TYPE1_NAME
#undef MACRO_TYPE2_NAME
#undef MACRO_TYPE3_NAME
#define MACRO_TYPE1_NAME real
#define MACRO_TYPE2_NAME integer
#define MACRO_TYPE3_NAME real
#include "binary_algebra.fpp"
#undef MACRO_TYPE1_NAME
#undef MACRO_TYPE2_NAME
#undef MACRO_TYPE3_NAME
#define MACRO_TYPE1_NAME real
#define MACRO_TYPE2_NAME real
#define MACRO_TYPE3_NAME real
#include "binary_algebra.fpp"
#undef MACRO_TYPE1_NAME
#undef MACRO_TYPE2_NAME
#undef MACRO_TYPE3_NAME
#define MACRO_TYPE1_NAME real
#define MACRO_TYPE2_NAME complex
#define MACRO_TYPE3_NAME complex
#include "binary_algebra.fpp"
#undef MACRO_TYPE1_NAME
#undef MACRO_TYPE2_NAME
#undef MACRO_TYPE3_NAME
#define MACRO_TYPE1_NAME complex
#define MACRO_TYPE2_NAME integer
#define MACRO_TYPE3_NAME complex
#include "binary_algebra.fpp"
#undef MACRO_TYPE1_NAME
#undef MACRO_TYPE2_NAME
#undef MACRO_TYPE3_NAME
#define MACRO_TYPE1_NAME complex
#define MACRO_TYPE2_NAME real
#define MACRO_TYPE3_NAME complex
#include "binary_algebra.fpp"
#undef MACRO_TYPE1_NAME
#undef MACRO_TYPE2_NAME
#undef MACRO_TYPE3_NAME
#define MACRO_TYPE1_NAME complex
#define MACRO_TYPE2_NAME complex
#define MACRO_TYPE3_NAME complex
#include "binary_algebra.fpp"
#undef MACRO_TYPE1_NAME
#undef MACRO_TYPE2_NAME
#undef MACRO_TYPE3_NAME
