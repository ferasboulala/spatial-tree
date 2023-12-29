#pragma once

#if ((defined(__GNUC__) || defined(__GNUG__)) && !defined(__clang__))
#define GNU_COMPILER
#elif (defined(__clang__))
#define CLANG_COMPILER
#endif

#if defined(GNU_COMPILER)
#define UNROLL_2 _Pragma("GCC unroll 2")
#define UNROLL_3 _Pragma("GCC unroll 3")
#define UNROLL_4 _Pragma("GCC unroll 4")
#define UNROLL_8 _Pragma("GCC unroll 8")
#elif defined(CLANG_COMPILER)
#define UNROLL_2 _Pragma("clang loop unroll_count(2)")
#define UNROLL_3 _Pragma("clang loop unroll_count(3)")
#define UNROLL_4 _Pragma("clang loop unroll_count(4)")
#define UNROLL_8 _Pragma("clang loop unroll_count(8)")
#define UNROLL(factor) DO_PRAGMA("clang loop unroll_count(" #factor ")")
#endif

#ifndef __always_inline
#define __always_inline __inline __attribute__ ((__always_inline__))
#endif