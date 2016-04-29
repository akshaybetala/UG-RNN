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
#include "../../cpp/src/Factor.cpp"
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

// Create MEX versions of WMB operators:
//   wmbFwd, wmbBwd, wmbImpSample

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  // check for the right number of arguments
  int minrhs,maxrhs;
  if(nrhs > 2 || nrhs < 1)				// note -- only checking for one!
    mexErrMsgTxt("Takes 2 input arguments");

  if(nlhs > 2)
    mexErrMsgTxt("Outputs 1-2 results"); 

  size_t nBatch = 1;
  if (nrhs == 2) nBatch = mxGetScalar(prhs[1]);


  const mxArray* m_wmb = prhs[0];
  const mxArray* m_alg = mxGetField(m_wmb,0,"Alg");

  vector< midx<uint32_t> > _order;             _order.mxSet( (mxArray*) mxGetField(m_alg,0,"order") );
  vector< vector<midx<uint32_t> > > _clique;   _clique.mxSet( (mxArray*) mxGetField(m_alg,0,"clique") );
  vector< Factor > _theta;                     _theta.mxSet( (mxArray*) mxGetField(m_alg,0,"theta") );
  //vector< Factor > _belief;
  vector<double> _parent;                      _parent.mxSet( (mxArray*) mxGetField(m_alg,0,"parent") );
  vector< vector<midx<uint32_t> > > _children; _children.mxSet( (mxArray*) mxGetField(m_alg,0,"children") );
  vector<double> _wt;                          _wt.mxSet( (mxArray*) mxGetField(m_alg,0,"wt") );
  vector< Factor > _msgFwd;                    _msgFwd.mxSet( (mxArray*) mxGetField(m_alg,0,"msgFwd") );
  vector< Factor > _msgBwd;                    _msgBwd.mxSet( (mxArray*) mxGetField(m_alg,0,"msgBwd") );

  // TODO: to uint32?
  vector< vector<midx<double> > > _nodes;           _nodes.mxSet( (mxArray*) mxGetField(m_alg,0,"nodes") ); 
  //vector< vector< vector<midx<uint32_t> > > > _match; _match.mxSet( (mxArray*) mxGetField(m_alg,0,"match") );
  
  //_factors.mxSet( mxGetField(m_alg,0,"factors") );

  size_t maxvar = 0;
  for (size_t i=0;i<_order.size();++i) if (maxvar < _order[i]) maxvar=_order[i]; 

  vector< size_t > idx; idx.resize(nBatch);
  vector< double > sample_wts;                    sample_wts.mxGet();  sample_wts.resize(nBatch);  
  vector< vector<midx<uint32_t> > > sample_vals;  sample_vals.mxGet(); sample_vals.resize(nBatch); 
  for (size_t s=0;s<nBatch;++s) sample_vals[s].resize(maxvar);

  for (size_t i=_order.size()-1;i<_order.size();--i) {
    Var X = Var( _order[i] , 0 );
    if (_nodes[i].size() == 0) continue;
    size_t nNodes = _nodes[i].size(); 

    // Sample which mini-bucket each sample is drawn from:
    vector<double> nodeWts; nodeWts.reserve(nNodes);
    for (size_t j=0;j<nNodes;++j) nodeWts.push_back(_wt[ _nodes[i][j] ]);
    Factor qMix(Var(0,nNodes),&nodeWts[0]); 
    for (size_t s=0;s<nBatch;++s) idx[s] = qMix.sample(); 

    // Compute the beliefs at each mini-bucket:
    vector<Factor> bel(nNodes);
    for (size_t j=0;j<nNodes;++j) {
      size_t n = _nodes[i][j];
      bel[j] = _theta[n];
      for (size_t cc=0;cc<_children[n].size();++cc) {
        size_t c = _children[n][cc];
        bel[j] += _msgFwd[c];
      }
      bel[j] -= _msgFwd[n];
    } // done computing beliefs

    // Now draw the samples from their respective mini-buckets:
    for (size_t s=0;s<nBatch;++s) {
      VarSet vCond = bel[idx[s]].vars()-X;
      Factor pr1 = bel[idx[s]].condition( vCond, sub2ind(vCond,sample_vals[s]) );
      pr1 *= (1.0/_wt[_nodes[i][idx[s]]]);
      pr1.exp(); 
      ind2sub(pr1.vars(),pr1.sample(),sample_vals[s]);  // write sample value into configuration tuple

      double qs=0.0;   // and compute (log) f(x) and q(x) contributions from this minibucket
      for (size_t j=0;j<nNodes;++j) {
        double tmp = bel[j][ sub2ind(bel[j].vars(), sample_vals[s]) ];
        sample_wts[s] += tmp;
        qs += std::exp(tmp*(1.0/_wt[_nodes[i][j]]))*_wt[_nodes[i][j]];
      }
      sample_wts[s] -= std::log( qs + 1e-300 );
    }
  }

  if (nlhs >  1) plhs[1] = sample_vals.mxGet();
  if (nlhs > -1) {
    //plhs[0] = sample_wts.mxGet();   
    mxArray *weights, *weightsT;        // Should transpose the vector for compatibility
    weights = sample_wts.mxGet();
    mexCallMATLAB(1,&weightsT,1,&weights,"transpose");
    plhs[0] = weightsT;
  }

}
