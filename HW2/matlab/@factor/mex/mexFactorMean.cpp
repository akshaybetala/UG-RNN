//
// Matlab MEX interface for Factors (test functions)
//
// Written by Alex Ihler
// Copyright (C) 2003 Alexander Ihler; distributable under GPL -- see README.txt
//

/*
addpath ~/code/matlab/shared/toolbox/
Fa=factor([1,2],rand(2,3)); Fb=factor([1],rand(2,1));
mex FactorSumOp.cpp -DSUM
Ftest = FactorSumOp(Fa,[2])

*/

#define MEX
#include <mex.h>
#include "VarSet.h"
#include "Factor.h"
#include "../../cpp/src/Factor.cpp"
#include <cstdio>
#include <iostream>
#include <cfloat>
#include <cmath>
#include <limits>
#include <stdint.h>

using namespace std;

// Create MEX versions of mean-like operators ("mean", "geomean")

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  // check for the right number of arguments
  if(nrhs < 1)
    mexErrMsgTxt("Takes >= 1 input arguments");
  if(nlhs > 1)
    mexErrMsgTxt("Outputs only one result"); 

	mex::Factor F1; F1.mxCopy(prhs[0]);

  for (int i=1;i<nrhs;i++) {
    mex::Factor F2; F2.mxSet((mxArray*)prhs[i]);
#ifdef MEAN
    F1 += F2;
#endif
#ifdef GEOMEAN
    F1 *= F2;
#endif
  }

#ifdef MEAN
  F1 /= (double) nrhs;
#endif
#ifdef GEOMEAN
  F1 ^= (1.0/nrhs);
#endif

  plhs[0]=F1.mxGet();						// return the pointer to the new copy

}
