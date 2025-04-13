/* Needed for Apple clang version 17 and above. Turn off certain
   aspects of floating point optimisation that interfere with the
   exact maths routines in triangle.c and prevent reproduction of some
   triangulations on the Mac M1 platform. */
#pragma STDC FENV_ACCESS ON
#pragma STDC FP_CONTRACT OFF

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <stdint.h>
#define printf Rprintf
#define TRILIBRARY
/* Needed to ensure test test-triangulate/"triangulate can triangulate
   an example that has crashed on Win i386" passes on Windows 32
   bit. See triangle.c for explanation of WIN32. */
#ifdef _WIN32
#define CPU86
#endif /* _WIN32 */
