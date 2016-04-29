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

// Create MEX versions of distance functions:
//   Entropy, KL, norms, etc.

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  // check for the right number of arguments
#ifdef ENTROPY
  if(nrhs !=1)
    mexErrMsgTxt("Takes 1 input arguments");
#endif
#ifdef NORM
  if(nrhs < 1 || nrhs > 2)
    mexErrMsgTxt("Takes 1-2 input arguments");
#endif
#ifdef DIST
  if(nrhs > 3 || nrhs < 2)
    mexErrMsgTxt("Takes 2-3 input arguments");
#endif
  if(nlhs > 1)
    mexErrMsgTxt("Outputs only one result"); 

	mex::Factor F1; F1.mxSet((mxArray*)prhs[0]);	// get 1st argument (a factor)
  double retval;
  int type;
#ifdef DIST
	mex::Factor F2; F2.mxSet((mxArray*)prhs[1]);	// distance functions need 2nd factor
  if (nrhs > 2)
	type = (int) mxGetScalar(prhs[2]);	// get method type
  else type = (int) mex::Factor::DIST_L2;	// default to L2-distance
#endif 
#ifdef NORM
  if (nrhs > 1)
	type = (int) mxGetScalar(prhs[1]);	// get method type
  else type = (int) mex::Factor::DIST_L2;	// default to L2-distance
#endif 
  
#ifdef ENTROPY
  retval = F1.entropy();
#endif
#ifdef DIST
  retval = F1.distance(F2, (Factor::DistType) type);
#endif
#ifdef NORM
  retval = F1.norm( (Factor::DistType) type);
#endif

  //if (nlhs)
    plhs[0] = mxCreateDoubleScalar( retval );

}
