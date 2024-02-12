#pragma once

#if ((defined(__GNUC__) || defined(__GNUG__)) && !defined(__clang__))
#define __gnu_compiler
#elif (defined(__clang__))
#define __clang_compiler
#endif

#if defined(__gnu_compiler)
#define __unroll_2 _Pragma("GCC unroll 2")
#define __unroll_3 _Pragma("GCC unroll 3")
#define __unroll_4 _Pragma("GCC unroll 4")
#define __unroll_8 _Pragma("GCC unroll 8")
#elif defined(__clang_compiler)
#define __unroll_2 _Pragma("clang loop unroll_count(2)")
#define __unroll_3 _Pragma("clang loop unroll_count(3)")
#define __unroll_4 _Pragma("clang loop unroll_count(4)")
#define __unroll_8 _Pragma("clang loop unroll_count(8)")
#endif

#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <functional>
#include <initializer_list>
#include <limits>
#include <numeric>
#include <tuple>
#include <vector>

namespace st {



}  // namespace st

#ifdef __gnu_compiler
#undef __gnu_compiler
#endif

#ifdef __clang_compiler
#undef __clang_compiler
#endif

#ifdef __unroll_2
#undef __unroll_2
#undef __unroll_3
#undef __unroll_4
#undef __unroll_8
#endif