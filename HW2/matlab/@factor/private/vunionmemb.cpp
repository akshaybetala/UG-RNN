//
// Matlab MEX union for sorted uint32 arrays
//

#include "mex.h"
#include <stdint.h>

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
  if (nrhs != 2) mexErrMsgTxt("Takes 2 input arguments");
  int N1=mxGetN(prhs[0]); int N2=mxGetN(prhs[1]);
  if (!mxIsUint32(prhs[0]) && N1>0) mexErrMsgTxt("Arguments must be sorted uint32s");
  if (!mxIsUint32(prhs[1]) && N2>0) mexErrMsgTxt("Arguments must be sorted uint32s");

  uint32_t *src, *cmp, *dest;
  uint8_t  *mem1, *mem2;
  src = (uint32_t *) mxGetData(prhs[0]);
  cmp = (uint32_t *) mxGetData(prhs[1]);
  
  plhs[0] = mxCreateNumericMatrix(1,N1+N2,mxUINT32_CLASS,mxREAL); 
  dest=(uint32_t *) mxGetData(plhs[0]);
  plhs[1] = mxCreateLogicalMatrix(1,N1+N2); 
  mem1=(uint8_t *) mxGetData(plhs[1]);
  plhs[2] = mxCreateLogicalMatrix(1,N1+N2); 
  mem2=(uint8_t *) mxGetData(plhs[2]);

  unsigned int i,j,k;
  for (i=0,j=0,k=0; i<N1 && j<N2; ) {
    if (src[i] < cmp[j] && i<N1) { mem1[k]=1; dest[k++]=src[i++]; }
    else if (src[i]==cmp[j]) { mem1[k]=mem2[k]=1; dest[k++]=src[i++]; j++; }
    else if (src[i] > cmp[j] && j<N2) { mem2[k]=1; dest[k++]=cmp[j++]; }
  }
  while (i<N1) { mem1[k]=1; dest[k++]=src[i++]; }
  while (j<N2) { mem2[k]=1; dest[k++]=cmp[j++]; }

  mxSetN(plhs[0],k);
  mxSetN(plhs[1],k);
  mxSetN(plhs[2],k);
  
}
