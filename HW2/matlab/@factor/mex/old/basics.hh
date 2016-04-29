#ifndef __BASICS_H
#define __BASICS_H

#ifdef MEX
#include "mex.h"
#endif

#include <stdint.h>

#ifdef MEX
// VarID		: uint32 by m-files ; used to index into # of dimensions, etc
// VarStates	: not in m-files ; mwSize or mwIndex (=size_t) by mex ; used to index into factor table
// Value        : double ; double
typedef uint32_t vindex;             // define "variable index" type (long)
typedef mwSize   vsize;              // define "variable size" (dimensions) type (long)
typedef double   value;              // define "factor value" type (double)
#else
typedef size_t   vindex;             // 
typedef size_t   vsize;              // 
typedef double   value;              // 
#endif

#endif
