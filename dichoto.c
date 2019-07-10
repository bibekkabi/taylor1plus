#include "tilings.h"  /* Include the header (not strictly necessary here) */


void dichoto(struct zonotope *zon1, struct match *matchindexcol)
{

int i;
int test;

//for(i=0;i<zon1->row;i++)
  //  zon1->mat[i][matchindexcol->matchcolumn]=matchindexcol->chck*zon1->mat[i][matchindexcol->matchcolumn];


for(i=0;i<zon1->row;i++)
    zon1->center[i]=zon1->center[i]+(matchindexcol->chck*zon1->mat[i][matchindexcol->matchcolumn]);

zon1->column--;

for(i=0;i<zon1->row;i++){
   test=matchindexcol->matchcolumn;
    while(test<zon1->column)
{
    zon1->mat[i][test]=zon1->mat[i][test+1];
    test++;
}

}


}


