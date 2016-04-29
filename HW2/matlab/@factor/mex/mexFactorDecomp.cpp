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
#include <string.h>

using namespace std;

// Create MEX versions of mean-like operators ("mean", "geomean")

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  // check for the right number of arguments
  if(nrhs != 3)
    mexErrMsgTxt("Takes 3 input arguments");
  if(nlhs > 1)
    mexErrMsgTxt("Outputs only one result"); 

	mex::Factor F1; F1.mxCopy(prhs[0]);                      // Get the original factor

	mxArray const* vlistCell=prhs[1];
	int nV=mxGetNumberOfElements(vlistCell);            // Get the number of factors to decompose into
	std::vector<mex::VarSet> vlist(nV);
	for (int i=0;i<nV;i++) {                            // Load in a list of variables from a cell array
	  vlist[i].mxSet(mxGetCell(vlistCell,i));           //   linking to each var list
		vlist[i].setDimsFromGlobal(F1.variables());       //   Also need to assign dimensions
	}

	char methodS[20];
	mex::Factor::DecompType method;
	bool found=false;
	mxGetString( prhs[2] , methodS, 20);
	for (int j=0;j<20;j++) methodS[j]=tolower(methodS[j]);
	if (!strcmp(methodS,"l2"))     { method = mex::Factor::DECOMP_L2;     found=true;}
	if (!strcmp(methodS,"l2+hpm")) { method = mex::Factor::DECOMP_L2HPM;  found=true;}
	if (!strcmp(methodS,"l2+mas")) { method = mex::Factor::DECOMP_L2MAS;  found=true;}
	if (!found) mexErrMsgTxt("Decomposition method not implemented. (mex)");

#ifdef SUM
	std::vector<mex::Factor> Flist = F1.decompSum( vlist , method );
#endif
#ifdef PROD
	std::vector<mex::Factor> Flist = F1.decompProd( vlist , method );
#endif
  plhs[0]=mxCreateCellMatrix( 1, Flist.size() );
	for (int j=0;j<Flist.size();j++) mxSetCell(plhs[0],j,Flist[j].mxGet());

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

