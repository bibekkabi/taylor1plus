#ifndef STACK_H_   /* Include guard */
#define STACK_H_
#include <stdio.h>
#include "tilings.h"

int MAXSIZE;

int top[5];         

struct zonotope stack[100][100];

struct zonotope tilings[100][100];

int isempty(int k); 
   
int isfull(int k);

struct zonotope peek(int k) ;

void pop(int k);

void push(struct zonotope *zon1, int k, int dernierk); 


#endif

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


