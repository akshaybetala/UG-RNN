//////////////////////////////////////////////////////////////////////////
// changeCell.cpp : Mex file for changing cell array entries by reference
//
// changeCell(C,i,A);  // C{i}=A, by reference
//
// WARNING! WARNING! WARNING!  This MEX file is not for the faint of heart.
//  It attempts to circumvent Matlab's memory-handling functions in order
//  to give more efficient code, and uses undocumented functions to do so.
//  It comes with no guarantees whatsoever.  It may not work on other
//  versions of Matlab, or other platforms, or even at all.  It may result
//  in memory leaks, segmentation faults, destroy data, set your computer on 
//  fire, or suck it into a black hole, for all I know.  If you do not accept 
//  these possible risks, do not use this MEX file.
//
//////////////////////////////////////////////////////////////////////////

#include "mex.h"

extern "C" mxArray *mxCreateSharedCopy(const mxArray *pr);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

void* tmp;
mxArray *cells, *data;

tmp = (void*) prhs[0]; cells = (mxArray*) tmp;
int ind = (int) mxGetScalar(prhs[1])-1;
if (ind<0) return;
tmp = (void*) prhs[2]; data = (mxArray*) tmp;

mxDestroyArray(mxGetCell(cells,ind)); 		    // Remove old cell contents
mxSetCell(cells,ind, mxCreateSharedCopy(prhs[2])); // Link to new contents

}

