//#ifndef STACK_H_   /* Include guard */
//#define STACK_H_
#include <stdio.h>
#include "tilings.h"
#include "stack.h"

int MAXSIZE = 100;       
struct zonotope stack[100][100];    

int top[5]={-1,-1,-1,-1,-1};            

int isempty(int k) {

   if(top[k] == -1)
      return 1;
   else
      return 0;
}
   
int isfull(int k) {

   if(top[k] == MAXSIZE)
      return 1;
   else
      return 0;
}

struct zonotope peek(int k) {
   return stack[top[k]][k];
}

void pop(int k) {	
   if(!isempty(k)) {
      top[k] = top[k] - 1;   
      //return stack[top][k];
   } else {
      printf("Could not retrieve data, Stack is empty.\n");
   }
}

void push(struct zonotope *zon1, int k, int dernierk) {
if(top[k]!=-1 && k!=0 && top[k]!=k && k!=dernierk){
top[k]=-1;
if(!isfull(k)) {
      top[k] = top[k] + 1;   
      stack[top[k]][k] = *zon1;
      //printf("%d \t",stack[top][k].column);
   } else {
      printf("Could not insert data, Stack is full.\n");
   }
}

if(top[k]==-1 || k==0 || k==dernierk) {
   if(!isfull(k)) {
      top[k] = top[k] + 1;   
      stack[top[k]][k] = *zon1;
      //printf("%d \t",stack[top][k].column);
   } else {
      printf("Could not insert data, Stack is full.\n");
   }
}

}

//#endif

/*int main() {
   // push items on to the stack 
   push(3);
   push(5);
   push(9);
   push(1);
   push(12);
   push(15);

   printf("Element at top of the stack: %d\n" ,peek());
   printf("Elements: \n");

   // print stack data 
   while(!isempty()) {
      int data = pop();
      printf("%d\n",data);
   }

   printf("Stack full: %s\n" , isfull()?"true":"false");
   printf("Stack empty: %s\n" , isempty()?"true":"false");
   
   return 0;
}*/


