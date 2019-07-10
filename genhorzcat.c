#include<stdio.h>
#include "tilings.h"  /* Include the header (not strictly necessary here) */

// horizontal concatenation
void genhorzcat (struct zonotope *p1, struct match *matchindexcol, int sgns[][100]){

int i,j;

int test=p1->column;


if(matchindexcol->chck==1){
       while(test>0){
          for(j=0;j<2*p1->column;j++){
               sgns[j][test]=sgns[j][test-1];
            }
          test--;
        }
       for(j=0;j<2*p1->column;j++)
            sgns[j][0]=1;
}


if(matchindexcol->chck==-1){
       while(test>0){
          for(j=0;j<2*p1->column;j++){
               sgns[j][test]=sgns[j][test-1];
            }
          test--;
        }
       for(j=0;j<2*p1->column;j++)
            sgns[j][0]=-1;
}


}
