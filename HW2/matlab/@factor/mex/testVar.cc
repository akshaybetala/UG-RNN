#include <stdio.h>
#include "Variables.h"
#include "subindex.h"

int main(void) {

 Variables v1(2);
 Variables v2(3);
 v1[0]=0; v1[1]=3;
 v2[0]=0; v2[1]=1; v2[2]=2;
 //printf("V1: %d entries\n",v.nvar());
 //for (int i=0;i<3;i++) printf("%d ",v2[i]); printf("\n");
 printf("V1: "); for (vindex i=0;i<v1.nvar();i++) printf("%d ",v1[i]); printf("\n");
 printf("V2: "); for (vindex i=0;i<v2.nvar();i++) printf("%d ",v2[i]); printf("\n");
 Variables U=v1+v2;
 printf("Union: "); for (vindex i=0;i<U.nvar();i++) printf("%d ",U[i]); printf("\n");
 Variables I=v1&v2;
 printf("Inter: "); for (vindex i=0;i<I.nvar();i++) printf("%d ",I[i]); printf("\n");
 Variables D=v1-v2;
 printf("D1: "); for (vindex i=0;i<D.nvar();i++) printf("%d ",D[i]); printf("\n");
 D=v2-v1;
 printf("D2: "); for (vindex i=0;i<D.nvar();i++) printf("%d ",D[i]); printf("\n");

 return 0;
}
