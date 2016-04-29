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
#include "set.h"
#include "../../cpp/src/Factor.cpp"
#include <cstdio>
#include <iostream>
#include <cfloat>
#include <cmath>
#include <limits>
#include <stdint.h>
#include <string.h>

using mex::Factor;
using mex::vector;

// Create MEX versions of mean-like operators ("mean", "geomean")

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  // check for the right number of arguments
  if(nrhs != 3)
    mexErrMsgTxt("Takes 3 input arguments");
  if(nlhs > 0)
    mexErrMsgTxt("Outputs no results (write by reference!)"); 

	Factor F1; F1.mxCopy(prhs[0]);                                // Copy the original factor
  Factor msgOut; msgOut.mxSet((mxArray*)prhs[1]);               // Write to out message by reference
  vector<Factor> msgIn; msgIn.mxSet((mxArray*)prhs[2]);    // Get in messages by reference

  for (size_t i=0;i<msgIn.size();++i) F1*=msgIn[i];
  F1.marginalInto( msgOut.vars(), msgOut );
  msgOut /= msgOut.max();

}

/*

addpath ~/code/matlab/shared/toolbox/
f=factor([1,2], rand(3,2));          
F1=decompProd(f, {uint32(1),uint32(2)} , 'l2')
mex mexFactorDecomp.cpp 
!mv mexFactorDecomp.mexa64 ../decompProdMX.mexa64
F2=decompProdMX(f, {uint32(1),uint32(2)} , 'l2')

f=factor([1,2,3], rand(3,2,2));          
F1=decompSum(f, {uint32(1),uint32([2 3])} , 'l2'); F1{1},
mex mexFactorDecomp.cpp 
!mv mexFactorDecomp.mexa64 ../decompProdMX.mexa64
F2=decompProdMX(f, {uint32(1),uint32([2 3])} , 'l2'); F2{1},

*/

