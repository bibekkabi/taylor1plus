#include "tilings.h"  /* Include the header (not strictly necessary here) */

// circular shift
void circshift1 (int* ind, int* indnew, int count1, int shift)
{


int j,cnt=0;

for (j=shift;j<count1;j++){
   indnew[cnt]=ind[j];
   cnt=cnt+1;
}

int itr=count1-cnt-1;

for(j=0;j<=itr;j++){
  indnew[cnt]=ind[j];
  cnt=cnt+1;
}

}
