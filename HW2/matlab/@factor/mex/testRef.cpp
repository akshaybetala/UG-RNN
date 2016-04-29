//
// Matlab MEX interface for Factors (test functions)
//
// Written by Alex Ihler
// Copyright (C) 2003 Alexander Ihler; distributable under GPL -- see README.txt
//

/*
addpath ~/code/matlab/shared/toolbox/
Fa=factor([1,2],rand(2,3)); Fb=factor([1],rand(2,1));
mex FactorBinOp.cpp -DPLUS
Ftest = FactorBinOp(Fa,Fb)

*/

#define MEX
#include "mex.h"
#include "../../include/Variables.h"
#include "../../include/Factor.h"
#include <cstdio>
#include <iostream>
#include <cfloat>
#include <cmath>
#include <limits>
#include <stdint.h>

using namespace std;

#ifdef PLUS
#define BINOP +=
#endif

#ifdef TIMES
#define BINOP *=
#endif

#ifdef RDIV
#define BINOP /=
#endif

#ifdef MINUS
#define BINOP -=
#endif


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  // check for the right number of arguments
  if(nrhs !=2)
    mexErrMsgTxt("Takes 2 input arguments");
  if(nlhs > 1)
    mexErrMsgTxt("Outputs only one result"); 

	mex::Factor F1,F2; 

  if (mxIsNumeric(prhs[0])) { 				// if the first argument is a numeric
    F1.mxSet((mxArray*)prhs[1]); 					// the second must be the factor class variable; copy it
		if (!mxIsEmpty(prhs[0]))
      F1 BINOP mxGetScalar(prhs[0]);			// and act directly
  } else {
    F1.mxSet((mxArray*)prhs[0]);						// otherwise, it *is* a factor class variable (assumed) so copy it,
    if (mxIsNumeric(prhs[1])) { 				//  if the 2nd arg is a scalar,
		  if (!mxIsEmpty(prhs[1]))
        F1 BINOP mxGetScalar(prhs[1]);		//  use that;
		} else {
      F2.mxSet((mxArray*)prhs[1]);			// otherwise, wrap it with a c++ class
      F1 BINOP F2;							//  and overwrite the new copy
	  }
  }
  plhs[0]=F1.mxGet();						// return the pointer to the new copy

}
