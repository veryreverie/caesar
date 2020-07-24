! ======================================================================
! Includes binary_sparse for all combinations of types.
! ======================================================================

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
#include "binary_sparse.fpp"
#undef MACRO_TYPE1_NAME
#undef MACRO_TYPE2_NAME
#undef MACRO_TYPE3_NAME
#define MACRO_TYPE1_NAME integer
#define MACRO_TYPE2_NAME real
#define MACRO_TYPE3_NAME real
#include "binary_sparse.fpp"
#undef MACRO_TYPE1_NAME
#undef MACRO_TYPE2_NAME
#undef MACRO_TYPE3_NAME
#define MACRO_TYPE1_NAME integer
#define MACRO_TYPE2_NAME complex
#define MACRO_TYPE3_NAME complex
#include "binary_sparse.fpp"
#undef MACRO_TYPE1_NAME
#undef MACRO_TYPE2_NAME
#undef MACRO_TYPE3_NAME
#define MACRO_TYPE1_NAME real
#define MACRO_TYPE2_NAME integer
#define MACRO_TYPE3_NAME real
#include "binary_sparse.fpp"
#undef MACRO_TYPE1_NAME
#undef MACRO_TYPE2_NAME
#undef MACRO_TYPE3_NAME
#define MACRO_TYPE1_NAME real
#define MACRO_TYPE2_NAME real
#define MACRO_TYPE3_NAME real
#include "binary_sparse.fpp"
#undef MACRO_TYPE1_NAME
#undef MACRO_TYPE2_NAME
#undef MACRO_TYPE3_NAME
#define MACRO_TYPE1_NAME real
#define MACRO_TYPE2_NAME complex
#define MACRO_TYPE3_NAME complex
#include "binary_sparse.fpp"
#undef MACRO_TYPE1_NAME
#undef MACRO_TYPE2_NAME
#undef MACRO_TYPE3_NAME
#define MACRO_TYPE1_NAME complex
#define MACRO_TYPE2_NAME integer
#define MACRO_TYPE3_NAME complex
#include "binary_sparse.fpp"
#undef MACRO_TYPE1_NAME
#undef MACRO_TYPE2_NAME
#undef MACRO_TYPE3_NAME
#define MACRO_TYPE1_NAME complex
#define MACRO_TYPE2_NAME real
#define MACRO_TYPE3_NAME complex
#include "binary_sparse.fpp"
#undef MACRO_TYPE1_NAME
#undef MACRO_TYPE2_NAME
#undef MACRO_TYPE3_NAME
#define MACRO_TYPE1_NAME complex
#define MACRO_TYPE2_NAME complex
#define MACRO_TYPE3_NAME complex
#include "binary_sparse.fpp"
#undef MACRO_TYPE1_NAME
#undef MACRO_TYPE2_NAME
#undef MACRO_TYPE3_NAME
