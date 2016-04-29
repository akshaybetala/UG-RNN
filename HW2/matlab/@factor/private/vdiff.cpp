//
// Matlab MEX setdiff for sorted uint32 arrays
//

#include <stdint.h>
#include <algorithm>

#include "mex.h"
void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
  if (nrhs != 2) mexErrMsgTxt("Takes 2 input arguments");
  int N1=mxGetN(prhs[0]); int N2=mxGetN(prhs[1]);
  if (!mxIsUint32(prhs[0]) && N1>0) mexErrMsgTxt("Arguments must be sorted uint32s");
  if (!mxIsUint32(prhs[1]) && N2>0) mexErrMsgTxt("Arguments must be sorted uint32s");

  uint32_t *src, *cmp, *dest;
  src = (uint32_t *) mxGetData(prhs[0]); 
  cmp = (uint32_t *) mxGetData(prhs[1]);
  
  plhs[0] = mxCreateNumericMatrix(1,N1,mxUINT32_CLASS,mxREAL); 
  dest=(uint32_t *) mxGetData(plhs[0]);

  uint32_t* end = std::set_difference( src, src+N1, cmp, cmp+N2, dest );
  mxSetN(plhs[0],end-dest);
/*	
  unsigned int i,j,k;
  for (i=0,j=0,k=0; i<N1 && j<N2; ) {
    if (src[i] < cmp[j]) dest[k++]=src[i++]; 
    else if (src[i]==cmp[j]) { i++; }
    else if (src[i] > cmp[j] && j<N2) j++;
  }
  while (i<N1) dest[k++]=src[i++];
  mxSetN(plhs[0],k);
*/
}

// Alternative octave file; does not play well with ".m" file in same directory...
#ifdef OCTAVE
#include <octave/oct.h>
DEFUN_DLD (vdiff, args, nargout, "Set difference for sorted uint32 vectors") {
  int nargin = args.length();
  if (nargin != 2) print_usage();
  else {
    uint32NDArray srcM = args(0).uint32_array_value();  
    octave_int<uint32_t>* src=srcM.fortran_vec();
    uint32NDArray cmpM = args(1).uint32_array_value();  
    octave_int<uint32_t>* cmp=cmpM.fortran_vec();
    size_t N1=srcM.nelem(), N2=cmpM.nelem();
    if (!error_state) {
      uint32NDArray dstM = srcM;
      octave_int<uint32_t>* dst=dstM.fortran_vec();
      octave_int<uint32_t>* end = std::set_difference( src, src+N1, cmp, cmp+N2, dst );
      dstM.resize(1,end-dst);
      return octave_value( dstM );
    }
  }
  return octave_value_list();
}
#endif

