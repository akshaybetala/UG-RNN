#include <stdio.h>
#include "Variables.h"
#include "subindex.h"

int main(void) {

 vindex full[] = {4,2,3,2};
 vindex sub[]  = {4,1,3,1};

 subindex idx(4,full,sub);

 for (int i=0;i<4;i++) printf("%d ",idx.skipped[i]); printf("\n");
 for (int i=0;i<4;i++) printf("%d ",idx.add[i]);     printf("\n");
 for (int i=0;i<4;i++) printf("%d ",idx.subtract[i]);printf("\n");

 for (int i=0;i<48;i++,++idx) printf("%d ",(int)((size_t)idx)); printf("\n");

 return 0;
}
