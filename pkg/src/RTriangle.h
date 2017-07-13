#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <stdint.h>
#undef PI
#define printf Rprintf
#define TRILIBRARY
/* Needed to ensure test test-triangulate/"triangulate can triangulate
   an example that has crashed on Win i386" passes on Windows 32
   bit. See triangle.c for explanation of WIN32. */
#ifdef _WIN32
#define CPU86
#endif /* _WIN32 */
