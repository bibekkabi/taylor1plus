#include <stdio.h>
#include "tilings.h"  /* Include the header (not strictly necessary here) */
#include <math.h>
#include "stack.h"

#define PI 3.14159265

void signs(struct zonotope *zon1, int sgns[][100]);

// generating signs 
void signsvec (struct zonotope *zon){

int i,j,k,l,m,n,p,q,r,s;
double negmat[100][100];
int sgns[100][100]; 
int sgnsold[100][100];
int ind=0;
int sgnsarr[100];
int flag,flag1;
int sn,ss,sg,sh,sy;
int cnt=0;
int cntfori;
int sgnsarr1[100][100];

struct zonotope zon1;
zon1=*zon;

if (zon1.column==2){
printf("o");
//tilings[ind][ind]=zon1;
//tilings[ind][ind].nbtiles=1;

for (i = 0; i < zon1.row; i++) {
    for (j = 0; j < zon1.column; j++) {
      printf ("%f", zon1.mat[i][j]);
    }

    printf ("\n");
  }

struct zonotope p2[2];

p2[0].row=2;
p2[0].column=2;

// center

double norm[2];
int id;

norm[0]=zon1.mat[0][0]*zon1.mat[0][0]+zon1.mat[1][0]*zon1.mat[1][0];
norm[1]=zon1.mat[0][1]*zon1.mat[0][1]+zon1.mat[1][1]*zon1.mat[1][1];

norm[0]=sqrt(norm[0]);
norm[1]=sqrt(norm[1]);

if (norm[0]>norm[1])
   id=0;
else
   id=1;

    p2[0].center[0]=zon1.center[0]-(zon1.mat[0][id])/2;
    p2[0].center[1]=zon1.center[1]-(zon1.mat[1][id])/2;
    //printf("%f \t", p2[n].center[k]);}

// generator
if (id==0)
p2[0].mat[0][0]=zon1.mat[0][0]/2;
else
p2[0].mat[0][0]=zon1.mat[0][0];
if (id==1)
p2[0].mat[0][1]=zon1.mat[0][1]/2;
else
p2[0].mat[0][1]=zon1.mat[0][1];
if (id==0)
p2[0].mat[1][0]=zon1.mat[1][0]/2;
else
p2[0].mat[1][0]=zon1.mat[1][0];
if (id==1)
p2[0].mat[1][1]=zon1.mat[1][1]/2;
else
p2[0].mat[1][1]=zon1.mat[1][1];

push(&p2[0],1,1);  // there will be two tilings 

tilings[ind][ind]=stack[top[1]][1];

p2[1].row=2;
p2[1].column=2;

// center

    p2[1].center[0]=zon1.center[0]+(zon1.mat[0][id])/2;
    p2[1].center[1]=zon1.center[1]+(zon1.mat[1][id])/2;
    //printf("%f \t", p2[n].center[k]);}

// generator

p2[1].mat[0][0]=p2[0].mat[0][0];
p2[1].mat[0][1]=p2[0].mat[0][1];
p2[1].mat[1][0]=p2[0].mat[1][0];
p2[1].mat[1][1]=p2[0].mat[1][1];

push(&p2[1],1,1);  // there will be two tilings 

ind++;
tilings[ind][ind]=stack[top[1]][1];

tilings[0][0].nbtiles=2;

return 0;
}

horzcat(&zon1, negmat);

signs(&zon1, sgns);

struct match matchindexcol;
matchindexcol.chck=1;
matchindexcol.matchcolumn=zon->column;

find (sgns, &matchindexcol);

if(zon1.column==zon1.totnocols)
push(&zon1,0,0);

dichoto(&zon1, &matchindexcol);

//printf("%f", zon1.center[0]);
//printf("%f", zon1.center[1]);

/*for(i=0;i<zon1.row;i++){
    for(j=0;j<zon1.column;j++){
        printf("%f \t",zon1.mat[i][j]);}

printf("\n");}*/

push(&zon1,0,0);


if(zon1.column!=2)
signsvec(&zon1);
else{

//printf("%d", stack[top[0]][0].column);


struct zonotope p1; //=peek(0);

int kint=1;

while(stack[top[0]][0].column!=stack[top[0]][0].totnocols){

p1=peek(0);

push(&p1,kint,kint-1);

pop(0);

horzcat(&stack[top[0]][0], negmat);

signs(&stack[top[0]][0],sgns);

for(i=0;i<2*stack[top[0]][0].column;i++){
   for(j=0;j<stack[top[0]][0].column;j++){
       sgnsold[i][j]=sgns[i][j];
  }
}

matchindexcol.matchcolumn=stack[top[0]][0].column;

find (sgns, &matchindexcol);

horzcat(&p1, negmat);
signs(&p1, sgns);

// general horzcat  sgns=horzcat(ones(size(sgns,1),1),sgns); || sgns=horzcat(ones(size(sgns,1),1),sgns);  

genhorzcat (&p1, &matchindexcol, sgns);

/*for(i=0;i<2*p1.column;i++){
    for(j=0;j<=p1.column;j++){
        printf("%d \t",sgns[i][j]);}

printf("\n");}

for(i=0;i<2*stack[top[0]][0].column;i++){
   for(j=0;j<stack[top[0]][0].column;j++){
      printf("%d \t",sgnsold[i][j]);}

printf("\n");}*/


if(p1.column==2){

tilings[ind][ind]=p1;

/*for(i=0;i<p1.row;i++){
    for(j=0;j<p1.column;j++){
        printf("%f \t",p1.mat[i][j]);}

printf("\n");}*/

//setxor [~,r]=setxor(sgns,signsvec(s{1}.top()),'rows');
cntfori=0;
for(i=0;i<2*p1.column;i++){
for(j=0;j<=p1.column;j++){
   sgnsarr[j]=sgns[i][j];}
for(k=0;k<2*stack[top[0]][0].column;k++){
   flag=0;
   for(j=0;j<stack[top[0]][0].column;j++){
      if(sgnsarr[j]==sgnsold[k][j])
         flag++;
      }
   if(flag==stack[top[0]][0].column)
     break;
}
if(k==2*stack[top[0]][0].column){
//break;
//vertcat 
cntfori++;
if (cntfori==1){
for(j=0;j<stack[top[0]][0].column;j++){
    sgnsold[2*stack[top[0]][0].column][j]=sgns[i][j];
}
}
}
}


for(i=0;i<=2*stack[top[0]][0].column;i++){
   for(j=0;j<stack[top[0]][0].column;j++){
      sgns[i][j]=sgnsold[i][j];}}

/*for(i=0;i<=2*stack[top[0]][0].column;i++){
   for(j=0;j<stack[top[0]][0].column;j++){
      printf("%d \t",sgns[i][j]);}

printf("\n");}*/


// subset temp=find(ismember(s{1}.top().Z(:,2:end)',p1.Z(:,2:end)','rows')) in the matlab code this part is before the if condition 

int temp[2]={0};
for(j=0;j<p1.column;j++){
   flag=0;
   for(k=0;k<stack[top[0]][0].column;k++){
      for(i=0;i<stack[top[0]][0].row;i++){   
         if(p1.mat[i][j]==stack[top[0]][0].mat[i][k])
            flag++;
         }
      if(flag==p1.row)
        break;
      }
      temp[j]=k;
      //printf("%d \t", flag);
}


int combs[2]={-1,1};


// if length(find(ismember(sgns(:,temp(j)),combs(m,:),'rows')))==4

int len_temp= sizeof(temp)/sizeof(temp[0]);
int len_combs= sizeof(combs)/sizeof(combs[0]);
int dim_combs= sizeof(combs)/(len_combs*sizeof(combs[0]));

//printf("%d \t", len_combs);

struct zonotope p2[len_temp];
for(n=0;n<len_temp;n++){
   for(m=0;m<len_combs;m++){
   flag=0;
   for(k=0;k<=2*stack[top[0]][0].column;k++){   
         if(combs[m]==sgns[k][temp[n]]){
            flag++;}
         if(flag==4){
           sn=combs[m];
           break;}
         }

//printf("%d \t", sn);

if(flag==4){ // this is just a cross check

p2[cnt].row=2;
p2[cnt].column=stack[top[0]][0].column;

// c=center(s{1}.top())+(sn(1))*G(:,temp(j))

for(k=0;k<stack[top[0]][0].row;k++){
    p2[cnt].center[k]=stack[top[0]][0].center[k]+sn*stack[top[0]][0].mat[k][temp[n]];}
    //printf("%f \t", p2[n].center[k]);}

// G(:,temp(j))=[];

for(i=0;i<stack[top[0]][0].row;i++){
   for(j=0;j<stack[top[0]][0].column;j++){
      p2[cnt].mat[i][j]=stack[top[0]][0].mat[i][j];}}

p2[cnt].column--;

for(i=0;i<p2[cnt].row;i++){
    flag=temp[n];
    while(flag<stack[top[0]][0].column)
{
    p2[cnt].mat[i][flag]=p2[cnt].mat[i][flag+1];
    flag++;
}

}

push(&p2[n],1,1);  // there will be two tilings 

cnt=cnt+1;
ind++;
tilings[ind][ind]=stack[top[1]][1];
tilings[0][0].nbtiles=tilings[0][0].nbtiles+cnt;

/*for(i=0;i<stack[top[1]][1].row;i++){
    for(j=0;j<stack[top[1]][1].column;j++){
        printf("%f \t",stack[top[1]][1].mat[i][j]);}

printf("\n");}*/


}  // end of check for the flag if it is equal to 4 if(flag==4)

} // for(m=0;m<len_combs;m++){

} // for(n=0;n<len_temp;n++){

} // end of if(p1.column==2){ 

// for temp==3 


/*for(i=0;i<2*p1.column;i++){
    for(j=0;j<=p1.column;j++){
        printf("%d \t",sgns[i][j]);}

printf("\n");}*/

//printf("\n");

/*for(i=0;i<2*stack[top[0]][0].column;i++){
   for(j=0;j<stack[top[0]][0].column;j++){
      printf("%d \t",sgnsold[i][j]);}

printf("\n");}*/


if(p1.column==3){
tilings[0][0].nbtiles=0;
cnt=0;
//setxor [~,r]=setxor(sgns,signsvec(s{1}.top()),'rows');
cntfori=0;
for(i=0;i<2*p1.column;i++){
for(j=0;j<=p1.column;j++){
   sgnsarr[j]=sgns[i][j];}
for(k=0;k<2*stack[top[0]][0].column;k++){
   flag=0;
   for(j=0;j<stack[top[0]][0].column;j++){
      if(sgnsarr[j]==sgnsold[k][j])
         flag++;
      }
   if(flag==stack[top[0]][0].column)
     break;
}
if(k==2*stack[top[0]][0].column){
//break;
//vertcat 
cntfori++;
if (cntfori==1){
for(j=0;j<stack[top[0]][0].column;j++){
    sgnsold[2*stack[top[0]][0].column][j]=sgns[i][j];
}
}
if (cntfori==2){
for(j=0;j<stack[top[0]][0].column;j++){
    sgnsold[2*(stack[top[0]][0].column)+1][j]=sgns[i][j];
}
}
}
}

for(i=0;i<=2*(stack[top[0]][0].column)+1;i++){
   for(j=0;j<stack[top[0]][0].column;j++){
      sgns[i][j]=sgnsold[i][j];}}

/*for(i=0;i<=2*(stack[top[0]][0].column)+1;i++){
   for(j=0;j<stack[top[0]][0].column;j++){
      printf("%d \t",sgns[i][j]);}

printf("\n");}*/

//cnt=0;

// subset temp=find(ismember(s{1}.top().Z(:,2:end)',p1.Z(:,2:end)','rows')) in the matlab code this part is before the if condition 

int temp[3]={0};
for(j=0;j<p1.column;j++){
   flag=0;
   for(k=0;k<stack[top[0]][0].column;k++){
      for(i=0;i<stack[top[0]][0].row;i++){   
         if(p1.mat[i][j]==stack[top[0]][0].mat[i][k])
            flag++;
         }
      if(flag==p1.row)
        break;
      }
      temp[j]=k;
      //printf("%d \t", temp[j]);
}


int combs[4][2]={-1,-1,-1,1,1,-1,1,1};

// if length(find(ismember(sgns(:,[temp(j) temp(p)]),combs(m,:),'rows')))==4

int len_temp= sizeof(temp)/sizeof(temp[0]);
int len_combs= sizeof(combs)/sizeof(combs[0]);
int dim_combs= sizeof(combs)/(len_combs*sizeof(combs[0][0]));
struct zonotope p2[len_temp];

//printf("%d \t", dim_combs);

for(n=0;n<len_temp;n++){
   for(p=0;p<len_temp;p++){
      if(n<p){
         for(k=0;k<=2*(stack[top[0]][0].column)+1;k++){
         sgnsarr1[k][0]=sgns[k][temp[n]];
         sgnsarr1[k][1]=sgns[k][temp[p]];}
         for(m=0;m<len_combs;m++){
         flag1=0;
         for(k=0;k<=2*(stack[top[0]][0].column)+1;k++){ 
            flag=0;
            for(j=0;j<dim_combs;j++){ 
               if(combs[m][j]==sgnsarr1[k][j])
                  flag++;}
               if(flag==2){
                  flag1++;}
                  //printf("%d \t", k);}
               if (flag1==4){
                  sn=combs[m][0];
                  ss=combs[m][1];
                  //printf("%d \t", sn);
                  //printf("%d \t", ss);
                  break;}
         }


if(flag1==4){ // this is just a cross check 

p2[cnt].row=2;
p2[cnt].column=stack[top[0]][0].column;

// c=center(s{1}.top())+(sn(1))*G(:,temp(j))+(ss(1))*G(:,temp(p));

for(k=0;k<stack[top[0]][0].row;k++){
    p2[cnt].center[k]=stack[top[0]][0].center[k]+sn*stack[top[0]][0].mat[k][temp[n]]+ss*stack[top[0]][0].mat[k][temp[p]];}
    //printf("%f \t", p2[cnt].center[k]);}


// G(:,temp(j))=[];

for(i=0;i<stack[top[0]][0].row;i++){
   for(j=0;j<stack[top[0]][0].column;j++){
      p2[cnt].mat[i][j]=stack[top[0]][0].mat[i][j];}}

p2[cnt].column=p2[cnt].column-2;

for(i=0;i<p2[cnt].row;i++){
    flag=temp[p];
    while(flag<stack[top[0]][0].column)
{
    p2[cnt].mat[i][flag]=p2[cnt].mat[i][flag+1];
    flag++;
}
}
for(i=0;i<p2[cnt].row;i++){
    flag=temp[n];
    while(flag<stack[top[0]][0].column-1)
{
    p2[cnt].mat[i][flag]=p2[cnt].mat[i][flag+1];
    flag++;
}
}


push(&p2[cnt],2,2);  // there will be three tilings

cnt=cnt+1;
ind++;
tilings[ind][ind]=stack[top[2]][2];
tilings[0][0].nbtiles=tilings[0][0].nbtiles+cnt;

/*for(i=0;i<stack[top[2]][2].row;i++){
    for(j=0;j<stack[top[2]][2].column;j++){
        printf("%f \t",stack[top[2]][2].mat[i][j]);}

printf("\n");}*/


        }
       }  
      }
    }
}


} // end of if(p1.column==3){


// for temp==4 


/*for(i=0;i<2*p1.column;i++){
    for(j=0;j<=p1.column;j++){
        printf("%d \t",sgns[i][j]);}

printf("\n");}*/

//printf("\n");

/*for(i=0;i<2*stack[top[0]][0].column;i++){
   for(j=0;j<stack[top[0]][0].column;j++){
      printf("%d \t",sgnsold[i][j]);}

printf("\n");}*/

if(p1.column==4){
tilings[0][0].nbtiles=0;
cnt=0;
//setxor [~,r]=setxor(sgns,signsvec(s{1}.top()),'rows');
cntfori=0;
for(i=0;i<2*p1.column;i++){
for(j=0;j<=p1.column;j++){
   sgnsarr[j]=sgns[i][j];}
for(k=0;k<2*stack[top[0]][0].column;k++){
   flag=0;
   for(j=0;j<stack[top[0]][0].column;j++){
      if(sgnsarr[j]==sgnsold[k][j])
         flag++;
      }
   if(flag==stack[top[0]][0].column)
     break;
}
if(k==2*stack[top[0]][0].column){
//break;
//vertcat 
cntfori++;
if (cntfori==1){
for(j=0;j<stack[top[0]][0].column;j++){
    sgnsold[2*stack[top[0]][0].column][j]=sgns[i][j];
}
}
if (cntfori==2){
for(j=0;j<stack[top[0]][0].column;j++){
    sgnsold[2*(stack[top[0]][0].column)+1][j]=sgns[i][j];
}
}
if (cntfori==3){
for(j=0;j<stack[top[0]][0].column;j++){
    sgnsold[2*(stack[top[0]][0].column)+2][j]=sgns[i][j];
}
}
}
}

for(i=0;i<=2*(stack[top[0]][0].column)+2;i++){
   for(j=0;j<stack[top[0]][0].column;j++){
      sgns[i][j]=sgnsold[i][j];}}

/*for(i=0;i<=2*(stack[top[0]][0].column)+2;i++){
   for(j=0;j<stack[top[0]][0].column;j++){
      printf("%d \t",sgns[i][j]);}

printf("\n");}*/

//cnt=0;

// subset temp=find(ismember(s{1}.top().Z(:,2:end)',p1.Z(:,2:end)','rows')) in the matlab code this part is before the if condition 

int temp[4]={0};
for(j=0;j<p1.column;j++){
   flag=0;
   for(k=0;k<stack[top[0]][0].column;k++){
      for(i=0;i<stack[top[0]][0].row;i++){   
         if(p1.mat[i][j]==stack[top[0]][0].mat[i][k])
            flag++;
         }
      if(flag==p1.row)
        break;
      }
      temp[j]=k;
      //printf("%d \t", temp[j]);
}


int combs[8][3]={-1,-1,-1,-1,-1,1,-1,1,-1,-1,1,1,1,-1,-1,1,-1,1,1,1,-1,1,1,1};

// if length(find(ismember(sgns(:,[temp(j) temp(p)]),combs(m,:),'rows')))==4

int len_temp= sizeof(temp)/sizeof(temp[0]);
int len_combs= sizeof(combs)/sizeof(combs[0]);
int dim_combs= sizeof(combs)/(len_combs*sizeof(combs[0][0]));
struct zonotope p2[len_temp];

//printf("%d \t", dim_combs);

for(n=0;n<len_temp;n++){
   for(p=0;p<len_temp;p++){
      for(q=0;q<len_temp;q++){
      if(n<p && p<q && n<q){
         for(k=0;k<=2*(stack[top[0]][0].column)+2;k++){
         sgnsarr1[k][0]=sgns[k][temp[n]];
         sgnsarr1[k][1]=sgns[k][temp[p]];
         sgnsarr1[k][2]=sgns[k][temp[q]];}
         for(m=0;m<len_combs;m++){
         flag1=0;
         for(k=0;k<=2*(stack[top[0]][0].column)+2;k++){ 
            flag=0;
            for(j=0;j<dim_combs;j++){ 
               if(combs[m][j]==sgnsarr1[k][j])
                  flag++;}
               if(flag==3){
                  flag1++;}
                  //printf("%d \t", k);}
               if (flag1==4){
                  sn=combs[m][0];
                  ss=combs[m][1];
                  sg=combs[m][2];
                  //printf("%d \t", sn);
                  //printf("%d \t", ss);
                  //printf("%d \t", sg);
                  break;}
         }


if(flag1==4){ // this is just a cross check 

p2[cnt].row=2;
p2[cnt].column=stack[top[0]][0].column;

// c=center(s{1}.top())+(sn(1))*G(:,temp(j))+(ss(1))*G(:,temp(p));

for(k=0;k<stack[top[0]][0].row;k++){
    p2[cnt].center[k]=stack[top[0]][0].center[k]+sn*stack[top[0]][0].mat[k][temp[n]]+ss*stack[top[0]][0].mat[k][temp[p]]+sg*stack[top[0]][0].mat[k][temp[q]];}
    //printf("%f \t", p2[cnt].center[k]);}

// G(:,temp(j))=[];

for(i=0;i<stack[top[0]][0].row;i++){
   for(j=0;j<stack[top[0]][0].column;j++){
      p2[cnt].mat[i][j]=stack[top[0]][0].mat[i][j];}}

p2[cnt].column=p2[cnt].column-3;

for(i=0;i<p2[cnt].row;i++){
    flag=temp[q];
    while(flag<stack[top[0]][0].column)
{
    p2[cnt].mat[i][flag]=p2[cnt].mat[i][flag+1];
    flag++;
}
}
for(i=0;i<p2[cnt].row;i++){
    flag=temp[p];
    while(flag<stack[top[0]][0].column-1)
{
    p2[cnt].mat[i][flag]=p2[cnt].mat[i][flag+1];
    flag++;
}
}
for(i=0;i<p2[cnt].row;i++){
    flag=temp[n];
    while(flag<stack[top[0]][0].column-2)
{
    p2[cnt].mat[i][flag]=p2[cnt].mat[i][flag+1];
    flag++;
}
}


push(&p2[cnt],3,3);  // there will be three tilings

cnt=cnt+1;
ind++;
tilings[ind][ind]=stack[top[3]][3];
tilings[0][0].nbtiles=tilings[0][0].nbtiles+cnt;

/*for(i=0;i<stack[top[3]][3].row;i++){
    for(j=0;j<stack[top[3]][3].column;j++){
        printf("%f \t",stack[top[3]][3].mat[i][j]);}

printf("\n");}*/


        }
       }  
      }
     }
    }
}


} // end of if(p1.column==4){

// for temp==5 

/*for(i=0;i<2*p1.column;i++){
    for(j=0;j<=p1.column;j++){
        printf("%d \t",sgns[i][j]);}

printf("\n");}*/

//printf("\n");

/*for(i=0;i<2*stack[top[0]][0].column;i++){
   for(j=0;j<stack[top[0]][0].column;j++){
      printf("%d \t",sgnsold[i][j]);}

printf("\n");}*/
if(p1.column==5){
tilings[0][0].nbtiles=0;
cnt=0;
//setxor [~,r]=setxor(sgns,signsvec(s{1}.top()),'rows');
cntfori=0;
for(i=0;i<2*p1.column;i++){
for(j=0;j<=p1.column;j++){
   sgnsarr[j]=sgns[i][j];}
for(k=0;k<2*stack[top[0]][0].column;k++){
   flag=0;
   for(j=0;j<stack[top[0]][0].column;j++){
      if(sgnsarr[j]==sgnsold[k][j])
         flag++;
      }
   if(flag==stack[top[0]][0].column)
     break;
}
if(k==2*stack[top[0]][0].column){
//break;
//vertcat 
cntfori++;
if (cntfori==1){
for(j=0;j<stack[top[0]][0].column;j++){
    sgnsold[2*stack[top[0]][0].column][j]=sgns[i][j];
}
}
if (cntfori==2){
for(j=0;j<stack[top[0]][0].column;j++){
    sgnsold[2*(stack[top[0]][0].column)+1][j]=sgns[i][j];
}
}
if (cntfori==3){
for(j=0;j<stack[top[0]][0].column;j++){
    sgnsold[2*(stack[top[0]][0].column)+2][j]=sgns[i][j];
}
}
if (cntfori==4){
for(j=0;j<stack[top[0]][0].column;j++){
    sgnsold[2*(stack[top[0]][0].column)+3][j]=sgns[i][j];
}
}
}
}


for(i=0;i<=2*(stack[top[0]][0].column)+3;i++){
   for(j=0;j<stack[top[0]][0].column;j++){
      sgns[i][j]=sgnsold[i][j];}}

/*for(i=0;i<=2*(stack[top[0]][0].column)+3;i++){
   for(j=0;j<stack[top[0]][0].column;j++){
      printf("%d \t",sgns[i][j]);}

printf("\n");}*/

//cnt=0

// subset temp=find(ismember(s{1}.top().Z(:,2:end)',p1.Z(:,2:end)','rows')) in the matlab code this part is before the if condition 

int temp[5]={0};
for(j=0;j<p1.column;j++){
   flag=0;
   for(k=0;k<stack[top[0]][0].column;k++){
      for(i=0;i<stack[top[0]][0].row;i++){   
         if(p1.mat[i][j]==stack[top[0]][0].mat[i][k])
            flag++;
         }
      if(flag==p1.row)
        break;
      }
      temp[j]=k;
      //printf("%d \t", temp[j]);
}


int combs[16][4]={-1,-1,-1,-1,-1,-1,-1,1,-1,-1,1,-1,-1,-1,1,1,-1,1,-1,-1,-1,1,-1,1,-1,1,1,-1,-1,1,1,1,1,-1,-1,-1,1,-1,-1,1,1,-1,1,-1,1,-1,1,1,1,1,-1,-1,1,1,-1,1,1,1,1,-1,1,1,1,1};

// if length(find(ismember(sgns(:,[temp(j) temp(p)]),combs(m,:),'rows')))==4

int len_temp= sizeof(temp)/sizeof(temp[0]);
int len_combs= sizeof(combs)/sizeof(combs[0]);
int dim_combs= sizeof(combs)/(len_combs*sizeof(combs[0][0]));
struct zonotope p2[len_temp];

//printf("%d \t", dim_combs);

for(n=0;n<len_temp;n++){
   for(p=0;p<len_temp;p++){
      for(q=0;q<len_temp;q++){
        for(r=0;r<len_temp;r++){
      if(n<p && p<q && q<r && p<r && n<r && n<q){
         for(k=0;k<=2*(stack[top[0]][0].column)+3;k++){
         sgnsarr1[k][0]=sgns[k][temp[n]];
         sgnsarr1[k][1]=sgns[k][temp[p]];
         sgnsarr1[k][2]=sgns[k][temp[q]];
         sgnsarr1[k][3]=sgns[k][temp[r]];}
         for(m=0;m<len_combs;m++){
         flag1=0;
         for(k=0;k<=2*(stack[top[0]][0].column)+3;k++){ 
            flag=0;
            for(j=0;j<dim_combs;j++){ 
               if(combs[m][j]==sgnsarr1[k][j])
                  flag++;}
               if(flag==4){
                  flag1++;}
                  //printf("%d \t", k);}
               if (flag1==4){
                  sn=combs[m][0];
                  ss=combs[m][1];
                  sg=combs[m][2];
                  sh=combs[m][3];
                  //printf("%d \t", sn);
                  //printf("%d \t", ss);
                  //printf("%d \t", sg);
                  break;}
         }


if(flag1==4){ // this is just a cross check 

p2[cnt].row=2;
p2[cnt].column=stack[top[0]][0].column;

// c=center(s{1}.top())+(sn(1))*G(:,temp(j))+(ss(1))*G(:,temp(p));

for(k=0;k<stack[top[0]][0].row;k++){
    p2[cnt].center[k]=stack[top[0]][0].center[k]+sn*stack[top[0]][0].mat[k][temp[n]]+ss*stack[top[0]][0].mat[k][temp[p]]+sg*stack[top[0]][0].mat[k][temp[q]]+sh*stack[top[0]][0].mat[k][temp[r]];}
    //printf("%f \t", p2[cnt].center[k]);}

// G(:,temp(j))=[];

for(i=0;i<stack[top[0]][0].row;i++){
   for(j=0;j<stack[top[0]][0].column;j++){
      p2[cnt].mat[i][j]=stack[top[0]][0].mat[i][j];}}

p2[cnt].column=p2[cnt].column-4;

for(i=0;i<p2[cnt].row;i++){
    flag=temp[r];
    while(flag<stack[top[0]][0].column)
{
    p2[cnt].mat[i][flag]=p2[cnt].mat[i][flag+1];
    flag++;
}
}
for(i=0;i<p2[cnt].row;i++){
    flag=temp[q];
    while(flag<stack[top[0]][0].column-1)
{
    p2[cnt].mat[i][flag]=p2[cnt].mat[i][flag+1];
    flag++;
}
}
for(i=0;i<p2[cnt].row;i++){
    flag=temp[p];
    while(flag<stack[top[0]][0].column-2)
{
    p2[cnt].mat[i][flag]=p2[cnt].mat[i][flag+1];
    flag++;
}
}
for(i=0;i<p2[cnt].row;i++){
    flag=temp[n];
    while(flag<stack[top[0]][0].column-3)
{
    p2[cnt].mat[i][flag]=p2[cnt].mat[i][flag+1];
    flag++;
}
}


push(&p2[cnt],4,4);  // there will be three tilings

cnt=cnt+1;
ind++;
tilings[ind][ind]=stack[top[4]][4];

tilings[0][0].nbtiles=tilings[0][0].nbtiles+cnt;

/*for(i=0;i<stack[top[4]][4].row;i++){
    for(j=0;j<stack[top[4]][4].column;j++){
        printf("%f \t",stack[top[4]][4].mat[i][j]);}

printf("\n");}*/


        }
       }  
       }
      }
     }
    }
}


} // end of if(p1.column==5){

// for temp==6 


/*for(i=0;i<2*p1.column;i++){
    for(j=0;j<=p1.column;j++){
        printf("%d \t",sgns[i][j]);}

printf("\n");}*/

//printf("\n");

/*for(i=0;i<2*stack[top[0]][0].column;i++){
   for(j=0;j<stack[top[0]][0].column;j++){
      printf("%d \t",sgnsold[i][j]);}

printf("\n");}*/
if(p1.column==6){
tilings[0][0].nbtiles=0;
cnt=0;
//setxor [~,r]=setxor(sgns,signsvec(s{1}.top()),'rows');
cntfori=0;
for(i=0;i<2*p1.column;i++){
for(j=0;j<=p1.column;j++){
   sgnsarr[j]=sgns[i][j];}
for(k=0;k<2*stack[top[0]][0].column;k++){
   flag=0;
   for(j=0;j<stack[top[0]][0].column;j++){
      if(sgnsarr[j]==sgnsold[k][j])
         flag++;
      }
   if(flag==stack[top[0]][0].column)
     break;
}
if(k==2*stack[top[0]][0].column){
//break;
//vertcat 
cntfori++;
if (cntfori==1){
for(j=0;j<stack[top[0]][0].column;j++){
    sgnsold[2*stack[top[0]][0].column][j]=sgns[i][j];
}
}
if (cntfori==2){
for(j=0;j<stack[top[0]][0].column;j++){
    sgnsold[2*(stack[top[0]][0].column)+1][j]=sgns[i][j];
}
}
if (cntfori==3){
for(j=0;j<stack[top[0]][0].column;j++){
    sgnsold[2*(stack[top[0]][0].column)+2][j]=sgns[i][j];
}
}
if (cntfori==4){
for(j=0;j<stack[top[0]][0].column;j++){
    sgnsold[2*(stack[top[0]][0].column)+3][j]=sgns[i][j];
}
}
if (cntfori==5){
for(j=0;j<stack[top[0]][0].column;j++){
    sgnsold[2*(stack[top[0]][0].column)+4][j]=sgns[i][j];
}
}
}
}


for(i=0;i<=2*(stack[top[0]][0].column)+4;i++){
   for(j=0;j<stack[top[0]][0].column;j++){
      sgns[i][j]=sgnsold[i][j];}}

/*for(i=0;i<=2*(stack[top[0]][0].column)+4;i++){
   for(j=0;j<stack[top[0]][0].column;j++){
      printf("%d \t",sgns[i][j]);}

printf("\n");}*/

//cnt=0

// subset temp=find(ismember(s{1}.top().Z(:,2:end)',p1.Z(:,2:end)','rows')) in the matlab code this part is before the if condition 

int temp[6]={0};
for(j=0;j<p1.column;j++){
   flag=0;
   for(k=0;k<stack[top[0]][0].column;k++){
      for(i=0;i<stack[top[0]][0].row;i++){   
         if(p1.mat[i][j]==stack[top[0]][0].mat[i][k])
            flag++;
         }
      if(flag==p1.row)
        break;
      }
      temp[j]=k;
      //printf("%d \t", temp[j]);
}


int combs[32][5]={-1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,-1,-1,1,-1,-1,-1,-1,1,1,-1,-1,1,-1,-1,-1,-1,1,-1,1,-1,-1,1,1,-1,-1,-1,1,1,1,-1,1,-1,-1,-1,-1,1,-1,-1,1,-1,1,-1,1,-1,-1,1,-1,1,1,-1,1,1,-1,-1,-1,1,1,-1,1,-1,1,1,1,-1,-1,1,1,1,1,1,-1,-1,-1,-1,1,-1,-1,-1,1,1,-1,-1,1,-1,1,-1,-1,1,1,1,-1,1,-1,-1,1,-1,1,-1,1,1,-1,1,1,-1,1,-1,1,1,1,1,1,-1,-1,-1,1,1,-1,-1,1,1,1,-1,1,-1,1,1,-1,1,1,1,1,1,-1,-1,1,1,1,-1,1,1,1,1,1,-1,1,1,1,1,1};

// if length(find(ismember(sgns(:,[temp(j) temp(p)]),combs(m,:),'rows')))==4

int len_temp= sizeof(temp)/sizeof(temp[0]);
int len_combs= sizeof(combs)/sizeof(combs[0]);
int dim_combs= sizeof(combs)/(len_combs*sizeof(combs[0][0]));
struct zonotope p2[len_temp];

//printf("%d \t", dim_combs);

for(n=0;n<len_temp;n++){
   for(p=0;p<len_temp;p++){
      for(q=0;q<len_temp;q++){
        for(r=0;r<len_temp;r++){
          for(s=0;s<len_temp;s++){
      if(n<p && p<q && q<r && r<s && q<s && p<s && n<s && p<r && n<r && n<q){
         for(k=0;k<=2*(stack[top[0]][0].column)+4;k++){
         sgnsarr1[k][0]=sgns[k][temp[n]];
         sgnsarr1[k][1]=sgns[k][temp[p]];
         sgnsarr1[k][2]=sgns[k][temp[q]];
         sgnsarr1[k][3]=sgns[k][temp[r]];
         sgnsarr1[k][4]=sgns[k][temp[s]];}
         for(m=0;m<len_combs;m++){
         flag1=0;
         for(k=0;k<=2*(stack[top[0]][0].column)+4;k++){ 
            flag=0;
            for(j=0;j<dim_combs;j++){ 
               if(combs[m][j]==sgnsarr1[k][j])
                  flag++;}
               if(flag==5){
                  flag1++;}
                  //printf("%d \t", k);}
               if (flag1==4){
                  sn=combs[m][0];
                  ss=combs[m][1];
                  sg=combs[m][2];
                  sh=combs[m][3];
                  sy=combs[m][4];
                  //printf("%d \t", sn);
                  //printf("%d \t", ss);
                  //printf("%d \t", sg);
                  break;}
         }


if(flag1==4){ // this is just a cross check 

p2[cnt].row=2;
p2[cnt].column=stack[top[0]][0].column;

// c=center(s{1}.top())+(sn(1))*G(:,temp(j))+(ss(1))*G(:,temp(p));

for(k=0;k<stack[top[0]][0].row;k++){
    p2[cnt].center[k]=stack[top[0]][0].center[k]+sn*stack[top[0]][0].mat[k][temp[n]]+ss*stack[top[0]][0].mat[k][temp[p]]+sg*stack[top[0]][0].mat[k][temp[q]]+sh*stack[top[0]][0].mat[k][temp[r]]+sy*stack[top[0]][0].mat[k][temp[s]];}
    //printf("%f \t", p2[cnt].center[k]);}

// G(:,temp(j))=[];

for(i=0;i<stack[top[0]][0].row;i++){
   for(j=0;j<stack[top[0]][0].column;j++){
      p2[cnt].mat[i][j]=stack[top[0]][0].mat[i][j];}}

p2[cnt].column=p2[cnt].column-5;

for(i=0;i<p2[cnt].row;i++){
    flag=temp[s];
    while(flag<stack[top[0]][0].column)
{
    p2[cnt].mat[i][flag]=p2[cnt].mat[i][flag+1];
    flag++;
}
}
for(i=0;i<p2[cnt].row;i++){
    flag=temp[r];
    while(flag<stack[top[0]][0].column-1)
{
    p2[cnt].mat[i][flag]=p2[cnt].mat[i][flag+1];
    flag++;
}
}
for(i=0;i<p2[cnt].row;i++){
    flag=temp[q];
    while(flag<stack[top[0]][0].column-2)
{
    p2[cnt].mat[i][flag]=p2[cnt].mat[i][flag+1];
    flag++;
}
}
for(i=0;i<p2[cnt].row;i++){
    flag=temp[p];
    while(flag<stack[top[0]][0].column-3)
{
    p2[cnt].mat[i][flag]=p2[cnt].mat[i][flag+1];
    flag++;
}
}
for(i=0;i<p2[cnt].row;i++){
    flag=temp[n];
    while(flag<stack[top[0]][0].column-4)
{
    p2[cnt].mat[i][flag]=p2[cnt].mat[i][flag+1];
    flag++;
}
}


push(&p2[cnt],5,5);  // there will be three tilings

cnt=cnt+1;
ind++;
tilings[ind][ind]=stack[top[5]][5];

tilings[0][0].nbtiles=tilings[0][0].nbtiles+cnt;

/*for(i=0;i<stack[top[5]][5].row;i++){
    for(j=0;j<stack[top[5]][5].column;j++){
        printf("%f \t",stack[top[5]][5].mat[i][j]);}

printf("\n");}*/



        }
       }  
       }
       }
      }
     }
    }
}


} // end of if(p1.column==5){

kint++;

} // while part

} // end of else of the if part i.e. if(zon1.column!=2)



} // emd of the function signsvec


void signs(struct zonotope *zon1, int sgns[][100]){

double angle1[100],angle2[100];
int negsgns[100][100];
int ind[100],indnew[100];
int c,d,j,i,k;
double t;

struct zonotope zon2;
zon2=*zon1;



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
