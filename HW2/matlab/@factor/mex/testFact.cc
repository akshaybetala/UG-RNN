#include <stdio.h>
#include "Variables.h"
#include "Factor.h"
#include <iostream>

using namespace std;

int main(void) {

 vindex Dims[] = {4,2,3,2,2};
 Variables v1(2);
 Variables v2(3);
 v1[0]=0; v1[1]=3;
 v2[0]=0; v2[1]=1; v2[2]=2;
 v1.setDimsFromGlobal(Dims); v2.setDimsFromGlobal(Dims); 

 Variables v3;
 v3=v2;

 const vsize *d = v3.dims();
 for (int i=0;i<v3.nvar();i++) cout << d[i] << " "; printf("\n");

 Factor F1(v1,0.1);
 //Factor F2(1.0);
 Factor F2(v2,1.1);
 F2^=0.25;
 printf("+: %f, %f\n",(F1+F2).sum(),(F1+3).sum());
 printf("-: %f, %f\n",(F1-F2).sum(),(F1-3).sum());
 printf("*: %f, %f\n",(F1*F2).sum(),(F1*3).sum());
 printf("/: %f, %f\n",(F1/F2).sum(),(F1/3).sum());
 printf("nvar: %d, numel: %d\n",F1.nvar(),F1.numel());
 printf("Is empty/nan/finite/scalar: %d %d %d %d\n",F1.isempty(),F1.isnan(),F1.isfinite(),F1.isscalar());

 printf("logsumexp: %f : %f\n",F1.logsumexp(F1.variables())[0], std::log((exp(F1)).sum()));
 printf("Norms: %f , %f \n",F1.norm(), F1.distance(Factor(v1,0.0)));

 // in-place arithmetic

 // Decomposition functions decompSum, decompProd

 // entropy, logpartition, normalize, normalized

 // minmarginal, maxmarginal, marginal, sum(v), max(v), sumPower


 return 0;
}
