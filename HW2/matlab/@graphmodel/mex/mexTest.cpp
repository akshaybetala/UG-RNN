//
// Written by Alex Ihler
// Copyright (C) 2015 Alexander Ihler; distributable under GPL -- see README.txt
//

/*
*/

#define MEX
#include <mex.h>
#include "VarSet.h"
#include "Factor.h"
#include "wmbe.h"
#include "../../cpp/src/Factor.cpp"
#include "../../cpp/src/graphmodel.cpp"
#include "../../cpp/src/wmbe.cpp"
#include <cstdio>
#include <iostream>
#include <cfloat>
#include <cmath>
#include <limits>
#include <stdint.h>

//using namespace std;

using mex::midx;
using mex::vector;
using mex::Var;
using mex::VarSet;
using mex::Factor;
using mex::wmbe;

// Create MEX versions of WMB operators:
//   wmbFwd, wmbBwd, wmbImpSample

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

/*
  // check for the right number of arguments
  int minrhs,maxrhs;
  if(nrhs > 2 || nrhs < 1)				// note -- only checking for one!
    mexErrMsgTxt("Takes 2 input arguments");

  if(nlhs > 2)
    mexErrMsgTxt("Outputs 1-2 results"); 

  size_t nBatch = 1;
  if (nrhs == 2) nBatch = mxGetScalar(prhs[1]);
*/

	wmbe gm;
	gm.mxSet( (mxArray*) prhs[0] );
	return;

}
