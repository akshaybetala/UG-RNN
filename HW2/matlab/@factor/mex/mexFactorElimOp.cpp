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

// Create MEX versions of elimination operators:
//   Sum, Marginal
//   LogSumExp
//   SumPower
//   Max, MaxMarginal
//   Min, MinMarginal

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  // check for the right number of arguments
  int minrhs,maxrhs;
#if defined(SUMPOWER) || defined(CONDITION)
  if(nrhs !=3)
    mexErrMsgTxt("Takes 3 input arguments");
#else
  if(nrhs > 2 || nrhs < 1)				// note -- only checking for one!
    mexErrMsgTxt("Takes 2 input arguments");
#endif
#if defined(MARGINAL) || defined(MAXMARGINAL) || defined(MINMARGINAL)
  if (nrhs < 2) mexErrMsgTxt("Takes 2 input arguments");
#endif
#ifdef CONDITION
  if (mxGetN(prhs[2])!=mxGetN(prhs[1])) mexErrMsgTxt("The number of variables should match the tuple length");
  if (mxGetM(prhs[2])!=1) mexErrMsgTxt("The number of tuples should be 1");
#endif

  if(nlhs > 1)
    mexErrMsgTxt("Outputs only one result"); 

	mex::Factor F1,F2; F1.mxSet((mxArray*)prhs[0]); 
  F2.mxGet();
  //F2.mxCopy(prhs[0]);						// maybe faster than creating a new one? depends on size...

	mex::VarSet v;
  if (nrhs == 1) {
    v = F1.vars();
  } else {
    if (mxGetClassID(prhs[1])==mxUINT32_CLASS) {		// check that our variable list is directly convertable
      v.mxSet((mxArray*)prhs[1]);
    } else {
      mxArray* conv;									// otherwise callback convert to uint32
      mexCallMATLAB(1,&conv,1,(mxArray**) &(prhs[1]),"uint32");
      v.mxSet(conv);
    }
  }


#ifdef SUM
  F2 = F1.sum( v );							// compute sum over desired variables
#endif
#ifdef LOGSUMEXP
  F2 = F1.logsumexp( v );
#endif
#ifdef SUMPOWER
  double pow = mxGetScalar(prhs[2]);
  F2 = F1.sumPower( v , pow );
#endif
#ifdef MAX
  F2 = F1.max( v );							//    or max, or min
#endif
#ifdef MIN
  F2 = F1.min( v );
#endif
#ifdef MARGINAL
  F2 = F1.marginal( F1.variables() & v );	// need to fix up dimensions (required?) & compute marginal
#endif
#ifdef MAXMARGINAL
  F2 = F1.maxmarginal( F1.variables() & v );//    or max, or min
#endif
#ifdef MINMARGINAL
  F2 = F1.minmarginal( F1.variables() & v );
#endif
#ifdef CONDITION
  double* vidx; size_t m,n;
  m = mxGetM((mxArray*)prhs[2]); n = mxGetN((mxArray*)prhs[2]);
  //if (n != v.size()) mexErrMsgTxt("The number of variables should match the tuple length");
  //if (m != 1) mexErrMsgTxt("The number of tuples should be 1");
  if (mxGetClassID(prhs[2])==mxDOUBLE_CLASS) {  // check that our variable list is directly convertable
    vidx = (double*) mxGetData(prhs[2]);
  } else {
    mxArray* conv;                              // otherwise callback convert to double
    mexCallMATLAB(1,&conv,1,(mxArray**) &(prhs[2]),"double");
    vidx = (double*) mxGetData(conv);
  }

  size_t maxIdx = 0;
  if (v.nvar() > 0) maxIdx = v.last().label()+1;
  mex::vector<size_t> sub(maxIdx);        // allocate storage for subscripts
  for (size_t vi=0;vi<v.nvar();++vi) sub[ v[vi] ] = ((size_t)vidx[vi])-1;  // TODO: better 0/1 indexing solution?
  mex::VarSet vs = F1.vars() & v;
  for (size_t vi=0;vi<vs.nvar();++vi) if (sub[vs[vi]] >= vs[vi].states()) {
    //mexErrMsgTxt("Invalid tuple; value outside variable range");   // TODO: can't exit with error on range error?
    sub[vs[vi]] = vs[vi].states()-1;
  }

  F2 = F1.condition(vs, sub2ind(vs,sub));  // Need "&" to ensure dimension information
  //F2 = F1.condition(v, sub2ind(F1.vars()&v,sub));  // Need "&" to ensure dimension information
#endif

  if (nrhs > 1) 
    plhs[0]=F2.mxGet();						// return the pointer to the new copy
  else
    plhs[0]=mxCreateDoubleScalar( F2[0] );  // or, return scalar if no vars specified

}
