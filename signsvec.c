#include <stdio.h>
#include "tilings.h"  /* Include the header (not strictly necessary here) */
#include <math.h>
#include <stdlib.h>

#define PI 3.14159265


void signsvec(struct zonotope *zon1, int** sgns){

//double angle1[100],angle2[100];
//int negsgns[100][100];
//int ind[100],indnew[100];
int c,d,j,i,k;
double t;

struct zonotope zon2;
zon2=*zon1;

double* angle1= (double* )calloc(2*zon2.column,sizeof(double));
double* angle2= (double* )calloc(2*zon2.column,sizeof(double));
int **negsgns = (int **) calloc (2*zon2.column, sizeof (int *));
    for (i = 0; i < 2*zon2.column; i++)
       negsgns[i] = (int *) calloc (zon2.column, sizeof (int));
int* ind= (int* )calloc(2*zon2.column,sizeof(int));
int* indnew= (int* )calloc(2*zon2.column,sizeof(int));


/*for(i=0;i<zon2.row;i++){
    for(j=0;j<zon2.column;j++){
        printf("%f \t",zon2.mat[i][j]);}

printf("\n");}*/

for (j=1;j<2*zon2.column;j++){
    angle1[j]=atan2((zon2.mat[0][0]*zon2.mat[1][j])-(zon2.mat[1][0]*zon2.mat[0][j]),(zon2.mat[0][0]*zon2.mat[0][j])+(zon2.mat[1][0]*zon2.mat[1][j])) * (180/PI);
    if (angle1[j]<0)
        angle1[j]=180+180+angle1[j];
//printf("%f \t", angle1[j]);
}

for(j=0;j<2*zon2.column;j++){
   angle2[j]=angle1[j];
}

int count=2*zon2.column;

for(c=1;c<count;++c){
  for(d=count-1; d>=c; --d){
     if (angle2[d-1] > angle2[d]){
        t=angle2[d-1];
        angle2[d-1]=angle2[d];
        angle2[d]=t;
       }
     }
  }

int count1=0;

for(j=0;j<2*zon2.column;j++){
   for(i=0;i<2*zon2.column;i++){
       if (angle2[j]==angle1[i]){
          ind[j]=i+1;
          count1=count1+1;
       }
      }
   }

for(j=0;j<2*zon2.column;j++){
   indnew[j]=0;
   //printf("%u \t", ind[j]);
}

int cnt=0;

for(j=0;j<2*zon2.column;j++){
  if (ind[j]<=zon2.column){ 
     circshift1(ind, indnew, count1, j);
     for(k=0;k<zon2.column;k++){
        if(indnew[k]<=zon2.column)
          sgns[cnt][indnew[k]-1]=1;
        else
          sgns[cnt][indnew[k]-zon2.column-1]=-1;
        }
     cnt=cnt+1;
     }
  }

for(i=0;i<zon2.column;i++){
    for(j=0;j<zon2.column;j++){
       negsgns[i][j]=-1*sgns[i][j]; 
        }
}

vertcat(sgns, negsgns, &zon2);


//zon1=&zon2;

/*for(i=0;i<2*zon2.column;i++){
    for(j=0;j<zon2.column;j++){
        printf("%d \t",sgns[i][j]);}

printf("\n");}*/

}
