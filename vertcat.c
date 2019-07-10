#include "tilings.h"  /* Include the header (not strictly necessary here) */


// vertical concatenation
void vertcat (int** sgns, int** negsgns, struct zonotope *zon2){

int i,j;

for(i=0;i<2*zon2->column;i++){
    for(j=0;j<zon2->column;j++){
       if (i>=zon2->column)
          sgns[i][j]=negsgns[i-zon2->column][j];
        }
}

}
