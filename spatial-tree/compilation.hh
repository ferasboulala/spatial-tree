#pragma once

#if ((defined(__GNUC__) || defined(__GNUG__)) && !defined(__clang__))
#define GNU_COMPILER
#elif (defined(__clang__))
#define CLANG_COMPILER
#endif

// TODO: Avoid the use of __always_inline and other aliases. Make that platform independent.
