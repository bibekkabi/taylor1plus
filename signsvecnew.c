#include <stdio.h>
#include "tilings.h"  
#include <math.h>
#include <stdlib.h>

#define PI 3.14159265


void signsvecnew(struct zonotope *zon1, int** sgns){

//double angle1[100],angle2[100];
//int negsgns[100][100];
//int ind[100],indnew[100];
int i,j,k;

struct zonotope zon2;
zon2=*zon1;

double* angle1= (double* )calloc(2*zon2.column,sizeof(double));
int **negsgns = (int **) calloc (2*zon2.column, sizeof (int *));
    for (i = 0; i < 2*zon2.column; i++)
       negsgns[i] = (int *) calloc (zon2.column, sizeof (int));

int* id= (int* )calloc(zon2.column,sizeof(int));

int* idneg= (int* )calloc(zon2.column,sizeof(int));

int* index= (int* )calloc(2*zon2.column,sizeof(int));

int* ind= (int* )calloc(2*zon2.column,sizeof(int));

int* index1= (int* )calloc(zon2.column,sizeof(int));

int* s= (int* )calloc(zon2.column,sizeof(int));

double temp1,temp2;


/*for(i=0;i<zon2.row;i++){
    for(j=0;j<zon2.column;j++){
        printf("%f \t",zon2.mat[i][j]);}

printf("\n");}*/

/*id=zeros(1,size(G,2));
id(1)=1;
for i=2:size(G,2)
    id(i)=id(i-1)+1;
end
idneg=-1*id;
index=horzcat(id,idneg);*/

id[0]=1;

for(i=1;i<zon2.column;i++)
   id[i]=id[i-1]+1;

for(i=0;i<zon2.column;i++){
   idneg[i]=-1*id[i];} // -(id[i]) 

for(i=0;i<2*zon2.column;i++){
   if(i>=zon2.column){
     index[i]=idneg[i-zon2.column];}
   else{
     index[i]=id[i];}}

for (j=1;j<2*zon2.column;j++){
    angle1[j]=atan2((zon2.mat[0][0]*zon2.mat[1][j])-(zon2.mat[1][0]*zon2.mat[0][j]),(zon2.mat[0][0]*zon2.mat[0][j])+(zon2.mat[1][0]*zon2.mat[1][j])) * (180/PI);
    if (angle1[j]<0)
        angle1[j]=180+180+angle1[j];
//printf("%f \t", angle1[j]);
}

/*for i=1:length(angle2)
for j=1:length(angle2)-i
if angle2(j)>=angle2(j+1)
temp1=angle2(j);
temp2=id(j);
angle2(j)=angle2(j+1);
id(j)=id(j+1);
angle2(j+1)=temp1;
id(j+1)=temp2;
end
end
end*/

ind[0]=1;

for(i=1;i<2*zon2.column;i++)
   ind[i]=ind[i-1]+1;

for(i=0;i<2*zon2.column;++i){
   for(j=0;j<(2*zon2.column)-i-1;++j){
      if (angle1[j] >= angle1[j+1]){
         temp1 =  angle1[j];
         temp2 = ind[j];
         angle1[j] = angle1[j+1];
         ind[j]=ind[j+1];
         angle1[j+1] = temp1;
         ind[j+1]=temp2;}}}

/*for i=1:length(ind)/2
    index1(i)=index(ind(i));
end*/

for(i=0;i<zon2.column;++i){
   index1[i]=index[ind[i]-1];}

/*index2=-1*index1;
index3=horzcat(index1,index2);*/

for(i=0;i<zon2.column;i++)
   idneg[i]=-1*index1[i];

for(i=0;i<2*zon2.column;i++){
   if(i>=zon2.column){
     index[i]=idneg[i-zon2.column];}
   else{
     index[i]=index1[i];}}

/*for i=1:size(G,2)
    s=zeros(1,size(G,2));
    s(abs(index3(i)))=sign(index3(i));
    for j=i+1:size(G,2)+i-1
        s(abs(index3(j)))=sign(index3(j));
    end
    s1(i,:)=s;
    s2(i,:)=-1*s;
end
sgns=vertcat(s1,s2);*/

for(i=0;i<zon2.column;++i){
   for(j=0;j<zon2.column;++j){
      s[j]=0;}
   s[abs(index[i])-1]=(index[i] > 0) ? 1 : -1;
   for(j=i+1;j<zon2.column+i;++j){
      s[abs(index[j])-1]=(index[j] > 0) ? 1 : -1;}
   for(j=0;j<zon2.column;++j){
      sgns[i][j]=s[j];
      negsgns[i][j]=-1*s[j];}}

vertcat(sgns, negsgns, &zon2);            

//zon1=&zon2;

/*for(i=0;i<2*zon2.column;i++){
    for(j=0;j<zon2.column;j++){
        printf("%d \t",sgns[i][j]);}

printf("\n");}*/

}
