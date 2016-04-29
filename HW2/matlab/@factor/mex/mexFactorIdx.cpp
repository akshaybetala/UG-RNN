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

// Create MEX versions of indexing and subscripting operations:
//   ind2sub
//   subv2ind
//   Condition? Sample? Others?

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

#ifdef IND2SUBV
// ind2sub(F, idx) = [v1 .. vk]
  if(nrhs > 2 || nrhs < 2) 
    mexErrMsgTxt("Takes 2 input arguments");
  if(nlhs > 1)
    mexErrMsgTxt("Outputs only one result");
  mex::Factor F1; F1.mxSet((mxArray*)prhs[0]); 

  // get idx = vector of doubles
  double* idx;
  size_t m,n; m = mxGetM((mxArray*)prhs[1]); if (m==1) m=mxGetN((mxArray*)prhs[1]);
  //mex::vector<double> idx;
  if (mxGetClassID(prhs[1])==mxDOUBLE_CLASS) {  // check that our variable list is directly convertable
    //idx.mxSet((mxArray*)prhs[1]);
    idx = (double*) mxGetData(prhs[1]);
  } else {
    mxArray* conv;                              // otherwise callback convert to uint32
    mexCallMATLAB(1,&conv,1,(mxArray**) &(prhs[1]),"double");
    //idx.mxSet(conv);
    idx = (double*) mxGetData(conv);
  }
  //size_t idx = mxGetScalar((mxArray*)prhs[1]);
  
  // manually build len(idx) x nvar() matrix of doubles & get pointer to memory
  //plhs[0] = mxCreateNumericMatrix( idx.size(), F1.nvar(), mxDOUBLE_CLASS, mxREAL);
  plhs[0] = mxCreateNumericMatrix( m, F1.nvar(), mxDOUBLE_CLASS, mxREAL);
  double* ret = (double*) mxGetData(plhs[0]);
  //mex::vector<double> ret(F1.nvar()); ret.mxGet();

  size_t maxIdx = 0;
  if (F1.vars().nvar() > 0) maxIdx = F1.vars().last().label()+1;
  mex::vector<size_t> sub(maxIdx);        // allocate storage for subscripts
  //for (size_t i=0; i<idx.size(); ++i) {
  for (size_t i=0; i<m; ++i) {
    mex::ind2sub(F1.vars(),idx[i]-1,sub);                              // TODO:  vvv vvv
    for (size_t v=0;v<F1.nvar();++v) ret[v*m+i] = sub[ F1.vars()[v] ]+1; // TODO: better 0/1 indexing sol'n
  }
  //plhs[0] = ret.mxGet();

//TODO: FIX: 0- vs 1-indexing subscripts (midx<uint32_t>)
//TODO: FIX: cell array for ind2sub, multiple indices (idx is a vector)
#endif

#ifdef IND2SUB
// // Reuse ind2subv, but make vector<vector<double> > to return a cell array?
// vector<vector<double> > ret(F1.nvar());
// for (size_t v=0; v<F1.nvar(); ++v) ret[v].resize(idx.size());
// ...
// ret[v][i] = ...
#endif

#ifdef SUBV2IND
// sub2ind(F, [v1..vk]) = idx
  if(nrhs > 2 || nrhs < 2)
    mexErrMsgTxt("Takes 2 input arguments");
  if(nlhs > 1)
    mexErrMsgTxt("Outputs only one result");
  mex::Factor F1; F1.mxSet((mxArray*)prhs[0]);
  //mex::vector<double> vidx;
  double* vidx; size_t m,n;
  m = mxGetM((mxArray*)prhs[1]); n = mxGetN((mxArray*)prhs[1]);
  if (n != F1.nvar()) mexErrMsgTxt("The number of variables in F should match the tuple length");
  plhs[0] = mxCreateNumericMatrix( m, 1, mxDOUBLE_CLASS, mxREAL);
  double* idx = (double*) mxGetData(plhs[0]);
  if (mxGetClassID(prhs[1])==mxDOUBLE_CLASS) {  // check that our variable list is directly convertable
    //vidx.mxSet((mxArray*)prhs[1]);
    vidx = (double*) mxGetData((mxArray*)prhs[1]);
  } else {
    mxArray* conv;                              // otherwise callback convert to uint32
    mexCallMATLAB(1,&conv,1,(mxArray**) &(prhs[1]),"double");
    //vidx.mxSet(conv);
    vidx = (double*) mxGetData(conv);
  }
  size_t maxIdx = 0;
  if (F1.vars().nvar() > 0) maxIdx = F1.vars().last().label()+1;
  mex::vector<size_t> sub(maxIdx);        // allocate storage for subscripts
  for (size_t i=0;i<m;++i) { 
    //for (size_t v=0;v<F1.nvar();++v) sub[ F1.vars()[v] ] = vidx[v]-1;  // TODO: better 0/1 indexing solution?
    for (size_t v=0;v<F1.nvar();++v) sub[ F1.vars()[v] ] = vidx[v*m+i]-1;  // TODO: better 0/1 indexing solution?
    idx[i] = (double) (sub2ind(F1.vars(), sub)+1);
  }
  //size_t idx = sub2ind(F1.vars(),sub)+1; 
  //plhs[0]=mxCreateDoubleScalar( idx ); 
#endif

#ifdef VALUE
  // take in matrix of tuples (subv) or vector (ind)
  // output column vector of values from lookup
  if(nrhs > 2 || nrhs < 2)
    mexErrMsgTxt("Takes 2 input arguments");
  if(nlhs > 1)
    mexErrMsgTxt("Outputs only one result");
  mex::Factor F1; F1.mxSet((mxArray*)prhs[0]);
  //mex::vector<double> vidx;
  double* vidx; size_t m,n;
  m = mxGetM((mxArray*)prhs[1]); n = mxGetN((mxArray*)prhs[1]);
  if (n != F1.nvar()) mexErrMsgTxt("The number of variables in F should match the tuple length");
  plhs[0] = mxCreateNumericMatrix( m, 1, mxDOUBLE_CLASS, mxREAL);
  double* idx = (double*) mxGetData(plhs[0]);
  if (mxGetClassID(prhs[1])==mxDOUBLE_CLASS) {  // check that our variable list is directly convertable
    vidx = (double*) mxGetData((mxArray*)prhs[1]);
  } else {
    mxArray* conv;                              // otherwise callback convert to uint32
    mexCallMATLAB(1,&conv,1,(mxArray**) &(prhs[1]),"double");
    vidx = (double*) mxGetData(conv);
  }
  size_t maxIdx = 0;
  if (F1.vars().nvar() > 0) maxIdx = F1.vars().last().label()+1;
  mex::vector<size_t> sub(maxIdx);        // allocate storage for subscripts
  for (size_t i=0;i<m;++i) { 
    for (size_t v=0;v<F1.nvar();++v) sub[ F1.vars()[v] ] = vidx[v*m+i]-1;  // TODO: better 0/1 indexing solution?
    idx[i] = F1[sub2ind(F1.vars(), sub)];
  }
  //size_t idx = sub2ind(F1.vars(),sub)+1; 
  //plhs[0]=mxCreateDoubleScalar( idx ); 


#endif 

#ifdef CONDITION
  if(nrhs > 3 || nrhs < 3)
    mexErrMsgTxt("Takes 3 input arguments");
  if(nlhs > 1)
    mexErrMsgTxt("Outputs only one result");
  mex::Factor F1,F2; F1.mxSet((mxArray*)prhs[0]); F2.mxGet();
  mex::VarSet vs; 
  if (mxGetClassID(prhs[1])==mxUINT32_CLASS) {    // check that our variable list is directly convertable
    vs.mxSet((mxArray*)prhs[1]);
  } else {
    mxArray* conv;                  // otherwise callback convert to uint32
    mexCallMATLAB(1,&conv,1,(mxArray**) &(prhs[1]),"uint32");
    vs.mxSet(conv);
  }
mexErrMsgTxt("something");
  double* vidx; size_t m,n;
  m = mxGetM((mxArray*)prhs[2]); n = mxGetN((mxArray*)prhs[2]);
  if (n != vs.size()) mexErrMsgTxt("The number of variables should match the tuple length");
  if (m != 1) mexErrMsgTxt("The number of tuples should be 1");
  if (mxGetClassID(prhs[2])==mxDOUBLE_CLASS) {  // check that our variable list is directly convertable
    vidx = (double*) mxGetData((mxArray*)prhs[2]);
  } else {
    mxArray* conv;                              // otherwise callback convert to uint32
    mexCallMATLAB(1,&conv,1,(mxArray**) &(prhs[2]),"double");
    vidx = (double*) mxGetData(conv);
  }
  size_t maxIdx = 0;
  if (F1.vars().nvar() > 0) maxIdx = F1.vars().last().label()+1;
  mex::vector<size_t> sub(maxIdx);        // allocate storage for subscripts
  for (size_t v=0;v<vs.nvar();++v) sub[ vs[v] ] = vidx[v]-1;  // TODO: better 0/1 indexing solution?
  F2 = F1.condition(vs, sub2ind(vs,sub));
  plhs[0] = F2.mxGet(); 
#endif 


#ifdef SAMPLE
  // Matlab callback to get radn(1,nSamples)
  // then do sampling with those random numbers?
#endif



/*
  // check for the right number of arguments
  int minrhs,maxrhs;
#ifdef CONDITION
  if(nrhs !=3)
    mexErrMsgTxt("Takes 3 input arguments");
#else
  if(nrhs > 2 || nrhs < 1)				// note -- only checking for one!
    mexErrMsgTxt("Takes 2 input arguments");
#endif
  if(nlhs > 1)
    mexErrMsgTxt("Outputs only one result"); 

  mex::Factor F1; F1.mxSet((mxArray*)prhs[0]); 

  // If output is a factor?
  mex::Factor F2; F2.mxGet();
  //F2.mxCopy(prhs[0]);		// maybe faster than creating a new one? depends on size...

// if 2nd argument is a variable set
  mex::VarSet v;
  if (mxGetClassID(prhs[1])==mxUINT32_CLASS) {	// check that our variable list is directly convertable
    v.mxSet((mxArray*)prhs[1]);
  } else {
    mxArray* conv;				// otherwise callback convert to uint32
    mexCallMATLAB(1,&conv,1,(mxArray**) &(prhs[1]),"uint32");
    v.mxSet(conv);
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
  if (nrhs < 2) mexErrMsgTxt("Takes 2 input arguments");
  F2 = F1.marginal( F1.variables() & v );	// need to fix up dimensions (required?) & compute marginal
#endif
#ifdef MAXMARGINAL
  if (nrhs < 2) mexErrMsgTxt("Takes 2 input arguments");
  F2 = F1.maxmarginal( F1.variables() & v );//    or max, or min
#endif
#ifdef MINMARGINAL
  if (nrhs < 2) mexErrMsgTxt("Takes 2 input arguments");
  F2 = F1.minmarginal( F1.variables() & v );
#endif

  if (nrhs > 1) 
    plhs[0]=F2.mxGet();						// return the pointer to the new copy
  else
    plhs[0]=mxCreateDoubleScalar( F2[0] );  // or, return scalar if no vars specified
*/

}
