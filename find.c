#include<stdio.h>
#include "tilings.h"  /* Include the header (not strictly necessary here) */


// find the row which is equal to 1 or -1
void find (int sgns[][100], struct match *matchindexcol){

int i,j,cnt;

for(i=0;i<matchindexcol->matchcolumn;i++){
cnt=0;
    for(j=0;j<2*matchindexcol->matchcolumn;j++){
        if(sgns[j][i]==matchindexcol->chck)
          cnt=cnt+1;
       }
    if(cnt>=matchindexcol->matchcolumn){
       break;}
    else {
        cnt=0;
        for(j=0;j<2*matchindexcol->matchcolumn;j++){
          if(sgns[j][i]==-1*matchindexcol->chck)
             cnt=cnt+1;
         }
        if(cnt>=matchindexcol->matchcolumn){
           matchindexcol->chck=-1*matchindexcol->chck;
           break;}     
       } 
  }

matchindexcol->matchcolumn=i;

//printf("%d \t", matchindexcol->matchcolumn);


//return stud;

}
