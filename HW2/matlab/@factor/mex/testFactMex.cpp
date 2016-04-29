// Basic testing function for mex code / compilation issues
//
// Written by Alex Ihler 
//

#define MEX
#include "mex.h"
#include "Variables.h"
#include "Factor.h"
#include <cstdio>
#include <iostream>
#include <cfloat>
#include <cmath>
#include <limits>
#include <stdint.h>

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  printf("Running %s...\n",mexFunctionName());

  // check for the right number of arguments
  //
  //
  //
  //
  vindex Dims[] = {4,2,3,2,2};
	mex::Variables v1(2);
	mex::Variables v2(3);
  v1[0]=0; v1[1]=3;
  v2[0]=0; v2[1]=1; v2[2]=2;
  v1.setDimsFromGlobal(Dims); v2.setDimsFromGlobal(Dims);

	mex::Variables v3;
  v3=v2;

  const vsize *d = v3.dims();
  for (int i=0;i<v3.nvar();i++) printf("%d ",(int)d[i]); printf("\n");
  //for (int i=0;i<v3.nvar();i++) cout << d[i] << " "; printf("\n");

	mex::Factor F1(v1,0.0);
	mex::Factor F2(v2,1.0);
  F2.abs();
	mex::Factor F=F1+F2;
  printf("%f\n",F.partition());

	//std::cout << (double) mxGetScalar(prhs[0]) << "\n";
	//std::cout << (double) mxGetScalar(prhs[1]) << "\n";
	//std::cout << (double) mxGetScalar(prhs[2]) << "\n";

}
