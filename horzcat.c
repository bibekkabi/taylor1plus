#include "tilings.h"  /* Include the header (not strictly necessary here) */

// horizontal concatenation
void horzcat (struct zonotope *zon1, double** negmat){

int i,j;

for(i=0;i<zon1->row;i++){
    for(j=0;j<zon1->column;j++){
       negmat[i][j]=-1*zon1->mat[i][j];
        }
}

for(i=0;i<zon1->row;i++){
    for(j=0;j<2*zon1->column;j++){
       if (j>=zon1->column)
          zon1->mat[i][j]=negmat[i][j-zon1->column];
        }
}

}
