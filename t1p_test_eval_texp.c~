#include <time.h>
#include <stdlib.h>
#include <limits.h>
#include "ap_global1.h"
#include "ap_interval.h"

#include "num.h"
#include "itv.h"
#include "t1p_internal.h"
#include "t1p_fun.h"
#include "t1p_constructor.h"
#include "t1p_representation.h"
#include "t1p_meetjoin.h"
#include "t1p_tilings.h"
#include "tilings.h"
#include "stack.h"


void
test_display (t1p_internal_t * pr1, t1p_t * env)
{
  printf ("### zonotope ###\n");
  ap_manager_t *man_tiles = t1p_manager_alloc ();
  t1p_internal_t *pr = t1p_init_from_manager (man_tiles, AP_FUNID_UNKNOWN);
  t1p_aff_check_free (pr, env->paf[0]);
  itv_t coeff;
  itv_init (coeff);
  num_t num;
  num_init (num);
  env->paf[0] = t1p_aff_alloc_init (pr);
  num_set_double (num, 20);
  itv_set_num (coeff, num);
  itv_set (env->paf[0]->c, coeff);
  num_set_double (num, -3);
  itv_set_num (coeff, num);
  t1p_nsym_t *symb1 = t1p_nsym_add (pr, IN);
  t1p_aff_nsym_add (pr, env->paf[0], coeff, symb1);
  //t1p_aff_nsym_create(pr, env->paf[0], coeff, symb1->type);
  //t1p_aff_build(pr, env->paf[0], coeff, 0);
  num_set_double (num, 5);
  itv_set_num (coeff, num);
  t1p_nsym_t *symb2 = t1p_nsym_add (pr, IN);
  t1p_aff_nsym_add (pr, env->paf[0], coeff, symb2);
  //t1p_aff_nsym_create(pr, env->paf[0], coeff, symb2->type);
  //t1p_aff_build(pr, env->paf[0], coeff, 1);
  num_set_double (num, 2);
  itv_set_num (coeff, num);
  t1p_nsym_t *symb3 = t1p_nsym_add (pr, IN);
  t1p_aff_nsym_add (pr, env->paf[0], coeff, symb3);
  //t1p_aff_nsym_create(pr, env->paf[0], coeff, symb3->type);
  //t1p_aff_build(pr, env->paf[0], coeff, 2);
  num_set_double (num, 1);
  itv_set_num (coeff, num);
  t1p_nsym_t *symb4 = t1p_nsym_add (pr, IN);
  t1p_aff_nsym_add (pr, env->paf[0], coeff, symb4);
  //t1p_aff_nsym_create(pr, env->paf[0], coeff, symb4->type);
  //t1p_aff_build(pr, env->paf[0], coeff, 3);
  num_set_double (num, 3);
  itv_set_num (coeff, num);
  t1p_nsym_t *symb5 = t1p_nsym_add (pr, IN);
  t1p_aff_nsym_add (pr, env->paf[0], coeff, symb5);
  //t1p_aff_nsym_create(pr, env->paf[0], coeff, symb5->type);
  //t1p_aff_build(pr, env->paf[0], coeff, 4);

  /*num_set_double(num, 0);
     itv_set_num(coeff, num);
     t1p_nsym_t* symb6=t1p_nsym_add(pr, IN);
     t1p_aff_nsym_add(pr, env->paf[0], coeff, symb6);
     /*num_set_double(num, 0);
     itv_set_num(coeff, num);
     t1p_nsym_t* symb7=t1p_nsym_add(pr, IN);
     t1p_aff_nsym_add(pr, env->paf[0], coeff, symb7);
     num_set_double(num, 0);
     itv_set_num(coeff, num);
     t1p_nsym_t* symb8=t1p_nsym_add(pr, IN);
     t1p_aff_nsym_add(pr, env->paf[0], coeff, symb8);
     num_set_double(num, 0);
     itv_set_num(coeff, num);
     t1p_nsym_t* symb9=t1p_nsym_add(pr, IN);
     t1p_aff_nsym_add(pr, env->paf[0], coeff, symb9);
     num_set_double(num, 0);
     itv_set_num(coeff, num);
     t1p_nsym_t* symb10=t1p_nsym_add(pr, IN);
     t1p_aff_nsym_add(pr, env->paf[0], coeff, symb10); */


  //t1p_aff_reduce(pr, env->paf[0]);


  //if (t1p_aff_is_perfectly_eq(pr, env->paf[2], env->paf[0])) printf("o");
  //else printf("n");


  t1p_aff_check_free (pr, env->paf[1]);
  env->paf[1] = t1p_aff_alloc_init (pr);
  num_set_double (num, 10);
  itv_set_num (coeff, num);
  itv_set (env->paf[1]->c, coeff);
  num_set_double (num, -4);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, env->paf[1], coeff, 0);
  num_set_double (num, 2);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, env->paf[1], coeff, 1);
  num_set_double (num, 1);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, env->paf[1], coeff, 2);
  num_set_double (num, 0);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, env->paf[1], coeff, 3);
  num_set_double (num, 5);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, env->paf[1], coeff, 4);
/*num_set_double(num, 6);
    itv_set_num(coeff, num);
    t1p_aff_build(pr, env->paf[1], coeff, 5);*/
/*    num_set_double(num, -4);
    itv_set_num(coeff, num);
    t1p_aff_build(pr, env->paf[1], coeff, 6);
num_set_double(num, 0);
    itv_set_num(coeff, num);
    t1p_aff_build(pr, env->paf[1], coeff, 7);
    num_set_double(num, 2);
    itv_set_num(coeff, num);
    t1p_aff_build(pr, env->paf[1], coeff, 8);
    //num_set_double(num, 0);
    //itv_set_num(coeff, num);
    //t1p_aff_build(pr, env->paf[1], coeff, 2);
    num_set_double(num, 1);
    itv_set_num(coeff, num);
    t1p_aff_build(pr, env->paf[1], coeff, 9);*/

  t1p_t *env2 = t1p_top (man_tiles, 0, 2);

  env2->paf[0] = t1p_aff_alloc_init (pr);
  num_set_double (num, 20);
  itv_set_num (coeff, num);
  itv_set (env2->paf[0]->c, coeff);
  num_set_double (num, 8);
  itv_set_num (coeff, num);
  t1p_nsym_t *symb6 = t1p_nsym_add (pr, IN);
  t1p_aff_nsym_add (pr, env2->paf[0], coeff, symb3);
  //t1p_aff_nsym_create(pr, env->paf[0], coeff, symb1->type);
  //t1p_aff_build(pr, env->paf[0], coeff, 0);
  num_set_double (num, 0);
  itv_set_num (coeff, num);
  t1p_nsym_t *symb7 = t1p_nsym_add (pr, IN);
  t1p_aff_nsym_add (pr, env2->paf[0], coeff, symb4);

  env2->paf[1] = t1p_aff_alloc_init (pr);
  num_set_double (num, 10);
  itv_set_num (coeff, num);
  itv_set (env2->paf[1]->c, coeff);
  num_set_double (num, 0);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, env2->paf[1], coeff, 2);
  num_set_double (num, 6);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, env2->paf[1], coeff, 3);

  printf ("%u", pr->dim);

  printf ("%u", pr->epsilon[8]->index);

/*t1p_t* tjoin;

tjoin= t1p_join(man_tiles, 'true',env,env2);

itv_print(tjoin->paf[0]->q->coeff);
itv_print(tjoin->paf[1]->q->coeff);*/


  bool res;

  res = t1p_is_leq (man_tiles, env, env2);

  printf ("%d\n", res);

  //itv_t** ptitv[env->dims][pr->dim];

  int i, j;

  itv_t **ptitv = (itv_t **) calloc (env->dims, sizeof (itv_t *));
  for (i = 0; i < env->dims; i++)
    ptitv[i] = (itv_t *) calloc (pr->dim, sizeof (itv_t));

  itv_t coeff1;
  itv_init (coeff1);



  for (i = 0; i < env->dims; i++) {
    for (j = 0; j < pr->dim; j++) {
      itv_init (ptitv[i][j]);	//=coeff1;
      //itv_print(ptitv[i][j]);
    }
  }

//itv_print(ptitv[0][0]);

  t1p_aaterm_t *p = NULL;
  t1p_aaterm_t *r = NULL;

  int *vec = (int *) calloc (pr->dim, sizeof (int));

  for (i = 0; i < pr->dim; i++) {
    vec[i] = -1;
  }

  int k = 0;

  for (i = 0; i < env->dims; i++) {
    for (p = env->paf[i]->q; p; p = p->n) {
      itv_init (ptitv[i][p->pnsym->index]);
      itv_set (ptitv[i][p->pnsym->index], p->coeff);
      itv_print (ptitv[i][p->pnsym->index]);
      if (vec[p->pnsym->index] == -1) {
	vec[p->pnsym->index] = k;
	k++;
      }
    }
  }

  printf ("%u", k);

  itv_t **ptitv1 = (itv_t **) calloc (env->dims, sizeof (itv_t *));
  for (i = 0; i < env->dims; i++)
    ptitv1[i] = (itv_t *) calloc (k, sizeof (itv_t));

  for (i = 0; i < env->dims; i++) {
    for (j = 0; j < k; j++) {
      itv_init (ptitv1[i][j]);	//=coeff1;
      //itv_print(ptitv[i][j]);
    }
  }

  t1p_aaterm_t *p1 = NULL;

  for (i = 0; i < env->dims; i++) {
    for (p1 = env->paf[i]->q; p1; p1 = p1->n) {
      itv_init (ptitv1[i][vec[p1->pnsym->index]]);
      itv_set (ptitv1[i][vec[p1->pnsym->index]], p1->coeff);
      itv_print (ptitv1[i][p1->pnsym->index]);

    }
  }

  itv_print (ptitv1[1][0]);
  itv_print (ptitv1[1][1]);
  itv_print (ptitv1[1][2]);

/// test ////

  itv_t **ptitv2 = (itv_t **) calloc (env->dims, sizeof (itv_t *));
  for (i = 0; i < env->dims; i++)
    ptitv2[i] = (itv_t *) calloc (pr->dim, sizeof (itv_t));

  int h;

  h = t1p_aff_get_generator_matrix_size (pr, env, ptitv2);

  printf ("%u", h);

  itv_t **ptitv3 = (itv_t **) calloc (env->dims, sizeof (itv_t *));
  for (i = 0; i < env->dims; i++)
    ptitv3[i] = (itv_t *) calloc (h, sizeof (itv_t));

  num_t num2;
  num_init (num2);
  num_set_double (num2, 0.4);

  t1p_aff_get_generator_matrix_num (pr, env, ptitv2, ptitv3, num2);

  itv_print (ptitv3[0][0]);
  itv_print (ptitv3[0][1]);
  itv_print (ptitv3[1][0]);
  itv_print (ptitv3[1][1]);



//printf("%u",sizeof(ptitv[0][1]));

//itv_print(ptitv[0][4]);


  /*printf("%u",env->paf[1]->q->n->n->pnsym->index);

     itv_print(env->paf[0]->q->n->n->coeff);

     itv_t** ptitv1[env->dims][pr->dim];

     t1p_aff_get_gen(pr, env, ptitv1);

     itv_print(ptitv1[0][0]);
     itv_print(ptitv1[0][1]);
     itv_print(ptitv1[0][2]);
     itv_print(ptitv1[0][3]);
     itv_print(ptitv1[0][4]);
     itv_print(ptitv1[1][0]);
     itv_print(ptitv1[1][1]);
     itv_print(ptitv1[1][2]);
     itv_print(ptitv1[1][3]);
     itv_print(ptitv1[1][4]); */

  //t1p_aff_reduce(pr, env->paf[1]);

  //printf("%u", pr->dim);

  //printf("%u",env->paf[0]->l);

  //printf("%u",env->paf[0]->end->pnsym->index);

  //ap_interval_t** intarr = ap_interval_array_alloc(pr->dim);

  //itv_internal_t* intern = itv_internal_alloc();
  //itv_t* ptitv= itv_array_alloc(2*pr->dim);
  //itv_array_set_ap_interval_array(intern, ptitv, intarr, pr->dim);
  //itv_t* res=t1p_aff_get_coeff(pr, env->paf[0], env->paf[0]->q->pnsym->index);
  //itv_print(res[0]);



  /*struct zonotope zon1;
     itv_magnitude(num, env->paf[0]->c);
     double_set_num(&zon1.center[0], num);
     itv_magnitude(num, env->paf[1]->c);
     double_set_num(&zon1.center[1], num);



     itv_magnitude(num, env->paf[0]->q->coeff);
     double_set_num(&zon1.mat[0][0], num);
     itv_magnitude(num, env->paf[1]->q->coeff);
     double_set_num(&zon1.mat[1][0], num);

     if (itv_is_neg(env->paf[0]->q->coeff)==1)
     zon1.mat[0][0]=-1*zon1.mat[0][0];
     if (itv_is_neg(env->paf[1]->q->coeff)==1)
     zon1.mat[1][0]=-1*zon1.mat[1][0];

     itv_magnitude(num, env->paf[0]->q->n->coeff);
     double_set_num(&zon1.mat[0][1], num);
     itv_magnitude(num, env->paf[1]->q->n->coeff);
     double_set_num(&zon1.mat[1][1], num);

     if (itv_is_neg(env->paf[0]->q->n->coeff)==1)
     zon1.mat[0][1]=-1*zon1.mat[0][1];
     if (itv_is_neg(env->paf[1]->q->n->coeff)==1)
     zon1.mat[1][1]=-1*zon1.mat[1][1];


     itv_magnitude(num, env->paf[0]->q->n->n->coeff);
     double_set_num(&zon1.mat[0][2], num);
     itv_magnitude(num, env->paf[1]->q->n->n->coeff);
     double_set_num(&zon1.mat[1][2], num);

     if (itv_is_neg(env->paf[0]->q->n->n->coeff)==1)
     zon1.mat[0][2]=-1*zon1.mat[0][2];
     if (itv_is_neg(env->paf[1]->q->n->n->coeff)==1)
     zon1.mat[1][2]=-1*zon1.mat[1][2];

     itv_magnitude(num, env->paf[0]->q->n->n->n->coeff);
     double_set_num(&zon1.mat[0][3], num);
     itv_magnitude(num, env->paf[1]->q->n->n->n->coeff);
     double_set_num(&zon1.mat[1][3], num);

     if (itv_is_neg(env->paf[0]->q->n->n->n->coeff)==1)
     zon1.mat[0][3]=-1*zon1.mat[0][3];
     if (itv_is_neg(env->paf[1]->q->n->n->n->coeff)==1)
     zon1.mat[1][3]=-1*zon1.mat[1][3];

     itv_magnitude(num, env->paf[0]->q->n->n->n->n->coeff);
     double_set_num(&zon1.mat[0][4], num);
     itv_magnitude(num, env->paf[1]->q->n->n->n->n->coeff);
     double_set_num(&zon1.mat[1][4], num);

     if (itv_is_neg(env->paf[0]->q->n->coeff)==1)
     zon1.mat[0][4]=-1*zon1.mat[0][4];
     if (itv_is_neg(env->paf[1]->q->n->coeff)==1)
     zon1.mat[1][4]=-1*zon1.mat[1][4];


     ap_abstract0_array_t s;

     ap_abstract0_t* env1 = ap_abstract0_cons(man_tiles,env);

     s = t1p_tilings(man_tiles, env1);

     t1p_t* tiles1 = (t1p_t*) (s.abs[0]->value);

     t1p_t* tiles2 = (t1p_t*) (s.abs[7]->value); */

/*itv_print(tiles1->paf[0]->c);
itv_print(tiles1->paf[1]->c);
itv_print(tiles1->paf[0]->q->coeff);
printf("%u",env->paf[0]->q->pnsym->index);
itv_print(tiles1->paf[1]->q->coeff);
//printf("%u",tiles1->paf[1]->q->pnsym->index);
itv_print(tiles1->paf[0]->q->n->coeff);
//printf("%u",tiles1->paf[0]->q->n->pnsym->index);
itv_print(tiles1->paf[1]->q->n->coeff);
printf("%u",tiles1->paf[1]->end->pnsym->index);*/


/*ap_interval_t** b = ap_interval_array_alloc(2);

ap_interval_set_double(b[0],12.0,28.0);
ap_interval_set_double(b[1],4.0,16.0);

ap_manager_t* manbx = t1p_manager_alloc();
t1p_internal_t * prbx = t1p_init_from_manager(manbx, AP_FUNID_OF_BOX);

//t1p_t* tb;
t1p_t* tb=t1p_of_box(manbx, 0, 2, b);

num_t num2;
num_init(num2);
num_set_double(num2, 0.4);

itv_t** ptitv2[tb->dims][prbx->dim];
t1p_aff_get_gen_num(prbx, tb, ptitv2,num2);

itv_print(ptitv2[0][0]);
itv_print(ptitv2[0][1]);*/

/*itv_t coeff1;
itv_init(coeff1);

int i,j;

for (i=0;i<tb->dims;i++){
   for (j=0;j<prbx->dim;j++){
       ptitv[i][j]=coeff1;
       //itv_print(ptitv[i][j]);
     }
}

//printf("%u",tb->paf[0]->q->pnsym->index);
     
for (i=0;i<tb->dims;i++){
   t1p_aaterm_t* cell = tb->paf[i]->q;
   for (j=0;j<prbx->dim;j++){
       if (cell != NULL){
          ptitv[i][cell->pnsym->index]=cell->coeff;
          cell = cell->n;
      }
     }
}

itv_print(ptitv[0][1]);*/


/*int k=0;
     int l;
  for(l=0;l<tb->dims;l++){
     t1p_aaterm_t* cell = tb->paf[l]->q;
     while (cell != NULL){
     itv_set(ptitv[k],cell->coeff);
     itv_print(ptitv[k]);
     k++;
     cell = cell->n;
 }
}*/

//printf("%u",tb->paf[0]->l);

/*ap_abstract0_array_t s1;

ap_abstract0_t* env2 = ap_abstract0_cons(manbx,tb);

s1 = t1p_tilings(manbx, env2);

t1p_t* tiles2 = (t1p_t*) (s1.abs[1]->value);*/

/////////////// meet operation //////////////////////


/*itv_print(tb->paf[0]->c);
itv_print(tb->paf[1]->c);
itv_print(tb->paf[0]->q->coeff);
//printf("%u",tb->paf[0]->q->pnsym->index);
itv_print(tb->paf[1]->q->coeff);
//printf("%u",tb->paf[1]->q->pnsym->index);

printf("%u",tb->paf[0]->end->pnsym->index);
printf("%u",tb->paf[1]->end->pnsym->index);*/

/*ap_interval_t** c = ap_interval_array_alloc(2);

ap_interval_set_double(c[0],-10.0,10.0);
ap_interval_set_double(c[1],-3.0,3.0);

t1p_t* tc;
tc=t1p_of_box(man_tiles, 0, 2, c);

t1p_t* tmeet;

tmeet= t1p_join(man_tiles, 'true',tb,tc);

itv_print(tmeet->paf[0]->q->coeff);
itv_print(tmeet->paf[1]->q->coeff);

printf("%u", prbx->dim);*/


//itv_print(tb->paf[0]->q->n->coeff);
//itv_print(tiles2->paf[1]->q->n->coeff);*/

/*itv_print(tiles2->paf[0]->c);
itv_print(tiles2->paf[1]->c);
itv_print(tiles2->paf[0]->q->coeff);
itv_print(tiles2->paf[1]->q->coeff);
itv_print(tiles2->paf[0]->q->n->coeff);
itv_print(tiles2->paf[1]->q->n->coeff);*/

//printf("%u",tiles2->paf[0]->l);


/*t1p_t* tiles0 = t1p_top(man_tiles, 0, 2);

tiles0->paf[0] = t1p_aff_alloc_init(pr);
    num_set_double(num, 28);
    itv_set_num(coeff, num);
    itv_set(tiles0->paf[0]->c, coeff);
    num_set_double(num, -3);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles0->paf[0], coeff, 0);
    num_set_double(num, 3);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles0->paf[0], coeff, 1);
tiles0->paf[1] = t1p_aff_alloc_init(pr);
    num_set_double(num, 13);
    itv_set_num(coeff, num);
    itv_set(tiles0->paf[1]->c, coeff);
    num_set_double(num, -4);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles0->paf[1], coeff, 0);
    num_set_double(num, 5);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles0->paf[1], coeff, 1);*/

//// ****** testing of inclusion **********

/* itv_t coeff1;
    itv_init(coeff1);
    num_t num1;
    num_init(num1);

num_set_double(num, -1);
itv_set_num(coeff, num);

num_set_double(num1, 1);
itv_set_num(coeff1, num1);

size_t i,j,k;

t1p_aaterm_t* p;
t1p_aaterm_t* ptr = t1p_aff_alloc_init(pr1);
t1p_aaterm_t* diffc = t1p_aff_alloc_init(pr1);
t1p_aaterm_t* exp1 = t1p_aff_alloc_init(pr1);
t1p_aaterm_t* exp2 = t1p_aff_alloc_init(pr1);

t1p_aaterm_t* ptmp1[100]; //= t1p_aff_alloc_init(pr1);
t1p_aaterm_t* ptmp2[100];
t1p_aaterm_t* lefteq= t1p_aff_alloc_init(pr1);
t1p_aaterm_t* righteq= t1p_aff_alloc_init(pr1);
t1p_aaterm_t* righteq1= t1p_aff_alloc_init(pr1);
t1p_aaterm_t* righteq2= t1p_aff_alloc_init(pr1);
t1p_aaterm_t* tmp= t1p_aff_alloc_init(pr1);

//// inequality's left part ///// for first generator

p=tiles2->paf[0]->q;
//for (p=tiles2->paf[0]->q;p;) {
for (k=0;k<tiles2->paf[0]->l;k++) {
    for (i=0;i<tiles2->dims;i++){
           if (i==0){
              p=tiles2->paf[0]->q;
              if (k!=0){
              p = p->n;}
              itv_mul(pr1->itv, ptr->coeff, coeff1, p->coeff);
              //itv_print(ptr->coeff);
              itv_sub(diffc->coeff, tiles2->paf[1]->c, tiles1->paf[1]->c);
              //itv_print(diffc->coeff);
              itv_mul(pr1->itv, exp1->coeff, diffc->coeff, ptr->coeff);
              //ptmp1[j]=t1p_aff_alloc_init(pr1);
              ptmp1[i]=exp1;
              printf("%zu",i);
              itv_print(ptmp1[i]->coeff);
              printf("%zu",i);
           }
           else {
              p=tiles2->paf[1]->q;
              if (k!=0){
              p = p->n;}
              //itv_print(ptmp1[0]->coeff);
              itv_mul(pr1->itv, ptr->coeff, coeff, p->coeff);
              //itv_print(ptr->coeff);
              itv_sub(diffc->coeff, tiles2->paf[0]->c, tiles1->paf[0]->c);
              //itv_print(diffc->coeff);
              itv_mul(pr1->itv, exp2->coeff, diffc->coeff, ptr->coeff);
              //ptmp2[j]=t1p_aff_alloc_init(pr1);
              ptmp1[i]=exp2;
              printf("%zu",i);
              itv_print(ptmp1[i]->coeff);
              printf("%zu",i);
           } 
         //itv_print(ptr->coeff);
               
}
//p = p->n;
itv_add(lefteq->coeff, ptmp1[0]->coeff, ptmp1[1]->coeff);
if (itv_is_neg(lefteq->coeff)==1)
         itv_mul(pr1->itv, lefteq->coeff, coeff, lefteq->coeff);
itv_print(lefteq->coeff);

p=tiles2->paf[0]->q;
    for (i=0;i<tiles2->dims;i++){
           if (i==0){
              p=tiles2->paf[1]->q;
              if (k!=0){
              p = p->n;}
              itv_mul(pr1->itv, ptr->coeff, coeff1, p->coeff);
              p=tiles2->paf[0]->q;
              itv_mul(pr1->itv, exp1->coeff, p->coeff, ptr->coeff);
              ptmp1[i]=exp1;
              itv_print(ptmp1[i]->coeff);
           }
           else {
              p=tiles2->paf[0]->q;
              if (k!=0){
              p = p->n;}
              itv_mul(pr1->itv, ptr->coeff, coeff, p->coeff);
              p=tiles2->paf[1]->q;
              itv_mul(pr1->itv, exp2->coeff, p->coeff, ptr->coeff);
              ptmp1[i]=exp2;
              itv_print(ptmp1[i]->coeff);
           } 
         //itv_print(ptr->coeff);
               
}

itv_add(righteq1->coeff, ptmp1[0]->coeff, ptmp1[1]->coeff);
if (itv_is_neg(righteq1->coeff)==1)
         itv_mul(pr1->itv, righteq1->coeff, coeff, righteq1->coeff);
itv_print(righteq1->coeff);

p=tiles2->paf[0]->q;
    for (i=0;i<tiles2->dims;i++){
           if (i==0){
              p=tiles2->paf[1]->q;
              if (k!=0){
              p = p->n;}
              itv_mul(pr1->itv, ptr->coeff, coeff, p->coeff);
              p=tiles2->paf[0]->q;
              p = p->n;
              itv_mul(pr1->itv, exp1->coeff, p->coeff, ptr->coeff);
              ptmp1[i]=exp1;
              itv_print(ptmp1[i]->coeff);
           }
           else {
              p=tiles2->paf[0]->q;
              if (k!=0){
              p = p->n;}
              itv_mul(pr1->itv, ptr->coeff, coeff1, p->coeff);
              p=tiles2->paf[1]->q;
              p = p->n;
              itv_mul(pr1->itv, exp2->coeff, p->coeff, ptr->coeff);
              ptmp1[i]=exp2;
              itv_print(ptmp1[i]->coeff);
           } 
         //itv_print(ptr->coeff);
               
}

itv_add(tmp->coeff, ptmp1[0]->coeff, ptmp1[1]->coeff);
if (itv_is_neg(tmp->coeff)==1)
         itv_mul(pr1->itv, tmp->coeff, coeff, tmp->coeff);
itv_add(righteq1->coeff, tmp->coeff, righteq1->coeff);
itv_print(righteq1->coeff);

//////////////////////////////////////////////////////////////////////

p=tiles1->paf[0]->q;
    for (i=0;i<tiles2->dims;i++){
           if (i==0){
              p=tiles2->paf[1]->q;
              if (k!=0){
              p = p->n;}
              itv_mul(pr1->itv, ptr->coeff, coeff1, p->coeff);
              p=tiles1->paf[0]->q;
              itv_mul(pr1->itv, exp1->coeff, p->coeff, ptr->coeff);
              ptmp1[i]=exp1;
              itv_print(ptmp1[i]->coeff);
           }
           else {
              p=tiles2->paf[0]->q;
              if (k!=0){
              p = p->n;}
              itv_mul(pr1->itv, ptr->coeff, coeff, p->coeff);
              p=tiles1->paf[1]->q;
              itv_mul(pr1->itv, exp2->coeff, p->coeff, ptr->coeff);
              ptmp1[i]=exp2;
              itv_print(ptmp1[i]->coeff);
           } 
         //itv_print(ptr->coeff);
               
}

itv_add(righteq2->coeff, ptmp1[0]->coeff, ptmp1[1]->coeff);
if (itv_is_neg(righteq2->coeff)==1)
         itv_mul(pr1->itv, righteq2->coeff, coeff, righteq2->coeff);
itv_print(righteq2->coeff);

p=tiles1->paf[0]->q;
    for (i=0;i<tiles2->dims;i++){
           if (i==0){
              p=tiles2->paf[1]->q;
              if (k!=0){
              p = p->n;}
              itv_mul(pr1->itv, ptr->coeff, coeff, p->coeff);
              p=tiles1->paf[0]->q;
              p = p->n;
              itv_mul(pr1->itv, exp1->coeff, p->coeff, ptr->coeff);
              ptmp1[i]=exp1;
              itv_print(ptmp1[i]->coeff);
           }
           else {
              p=tiles2->paf[0]->q;
              if (k!=0){
              p = p->n;}
              itv_mul(pr1->itv, ptr->coeff, coeff1, p->coeff);
              p=tiles1->paf[1]->q;
              p = p->n;
              itv_mul(pr1->itv, exp2->coeff, p->coeff, ptr->coeff);
              ptmp1[i]=exp2;
              itv_print(ptmp1[i]->coeff);
           } 
         //itv_print(ptr->coeff);
               
}

itv_add(tmp->coeff, ptmp1[0]->coeff, ptmp1[1]->coeff);
if (itv_is_neg(tmp->coeff)==1)
         itv_mul(pr1->itv, tmp->coeff, coeff, tmp->coeff);

itv_add(righteq2->coeff, tmp->coeff, righteq2->coeff);
itv_print(righteq2->coeff);

itv_sub(righteq->coeff, righteq1->coeff, righteq2->coeff);
itv_print(righteq->coeff);

//itv_print(lefteq->coeff);

if (itv_is_leq(righteq->coeff,lefteq->coeff))
   printf("o");

}*/



//t1p_t* tiles1 = (t1p_t*) (s.abs[5]->value);


/*itv_add(lefteq->coeff, ptmp1[0]->coeff, ptmp1[1]->coeff);
if (itv_is_neg(lefteq->coeff)==1)
         itv_mul(pr1->itv, lefteq->coeff, coeff, lefteq->coeff);
itv_print(lefteq->coeff);

//// inequality's right part ///// for first generator

p=tiles2->paf[0]->q;
    for (i=0;i<tiles2->dims;i++){
           if (i==0){
              p=tiles2->paf[1]->q;
              itv_mul(pr1->itv, ptr->coeff, coeff1, p->coeff);
              p=tiles2->paf[0]->q;
              itv_mul(pr1->itv, exp1->coeff, p->coeff, ptr->coeff);
              ptmp1[i]=exp1;
              itv_print(ptmp1[i]->coeff);
           }
           else {
              p=tiles2->paf[0]->q;
              itv_mul(pr1->itv, ptr->coeff, coeff, p->coeff);
              p=tiles2->paf[1]->q;
              itv_mul(pr1->itv, exp2->coeff, p->coeff, ptr->coeff);
              ptmp1[i]=exp2;
              itv_print(ptmp1[i]->coeff);
           } 
         //itv_print(ptr->coeff);
               
}

itv_add(righteq1->coeff, ptmp1[0]->coeff, ptmp1[1]->coeff);
if (itv_is_neg(righteq1->coeff)==1)
         itv_mul(pr1->itv, righteq1->coeff, coeff, righteq1->coeff);
itv_print(righteq1->coeff);

p=tiles2->paf[0]->q;
    for (i=0;i<tiles2->dims;i++){
           if (i==0){
              p=tiles2->paf[1]->q;
              itv_mul(pr1->itv, ptr->coeff, coeff, p->coeff);
              p=tiles2->paf[0]->q;
              p = p->n;
              itv_mul(pr1->itv, exp1->coeff, p->coeff, ptr->coeff);
              ptmp1[i]=exp1;
              itv_print(ptmp1[i]->coeff);
           }
           else {
              p=tiles2->paf[0]->q;
              itv_mul(pr1->itv, ptr->coeff, coeff1, p->coeff);
              p=tiles2->paf[1]->q;
              p = p->n;
              itv_mul(pr1->itv, exp2->coeff, p->coeff, ptr->coeff);
              ptmp1[i]=exp2;
              itv_print(ptmp1[i]->coeff);
           } 
         //itv_print(ptr->coeff);
               
}

itv_add(tmp->coeff, ptmp1[0]->coeff, ptmp1[1]->coeff);
if (itv_is_neg(tmp->coeff)==1)
         itv_mul(pr1->itv, tmp->coeff, coeff, tmp->coeff);
itv_add(righteq1->coeff, tmp->coeff, righteq1->coeff);
itv_print(righteq1->coeff);

//////////////////////////////////////////////////////////////////////

p=tiles1->paf[0]->q;
    for (i=0;i<tiles2->dims;i++){
           if (i==0){
              p=tiles2->paf[1]->q;
              itv_mul(pr1->itv, ptr->coeff, coeff1, p->coeff);
              p=tiles1->paf[0]->q;
              itv_mul(pr1->itv, exp1->coeff, p->coeff, ptr->coeff);
              ptmp1[i]=exp1;
              itv_print(ptmp1[i]->coeff);
           }
           else {
              p=tiles2->paf[0]->q;
              itv_mul(pr1->itv, ptr->coeff, coeff, p->coeff);
              p=tiles1->paf[1]->q;
              itv_mul(pr1->itv, exp2->coeff, p->coeff, ptr->coeff);
              ptmp1[i]=exp2;
              itv_print(ptmp1[i]->coeff);
           } 
         //itv_print(ptr->coeff);
               
}

itv_add(righteq2->coeff, ptmp1[0]->coeff, ptmp1[1]->coeff);
if (itv_is_neg(righteq2->coeff)==1)
         itv_mul(pr1->itv, righteq2->coeff, coeff, righteq2->coeff);
itv_print(righteq2->coeff);

p=tiles1->paf[0]->q;
    for (i=0;i<tiles2->dims;i++){
           if (i==0){
              p=tiles2->paf[1]->q;
              itv_mul(pr1->itv, ptr->coeff, coeff, p->coeff);
              p=tiles1->paf[0]->q;
              p = p->n;
              itv_mul(pr1->itv, exp1->coeff, p->coeff, ptr->coeff);
              ptmp1[i]=exp1;
              itv_print(ptmp1[i]->coeff);
           }
           else {
              p=tiles2->paf[0]->q;
              itv_mul(pr1->itv, ptr->coeff, coeff1, p->coeff);
              p=tiles1->paf[1]->q;
              p = p->n;
              itv_mul(pr1->itv, exp2->coeff, p->coeff, ptr->coeff);
              ptmp1[i]=exp2;
              itv_print(ptmp1[i]->coeff);
           } 
         //itv_print(ptr->coeff);
               
}

itv_add(tmp->coeff, ptmp1[0]->coeff, ptmp1[1]->coeff);
if (itv_is_neg(tmp->coeff)==1)
         itv_mul(pr1->itv, tmp->coeff, coeff, tmp->coeff);

itv_add(righteq2->coeff, tmp->coeff, righteq2->coeff);
itv_print(righteq2->coeff);

itv_sub(righteq->coeff, righteq1->coeff, righteq2->coeff);
itv_print(righteq->coeff);




//// inequality's left part ///// for second generator

p=tiles2->paf[0]->q;
p = p->n;
for (i=0;i<tiles2->dims;i++){
           if (i==0){
              itv_mul(pr1->itv, ptr->coeff, coeff1, p->coeff);
              //itv_print(ptr->coeff);
              itv_sub(diffc->coeff, tiles2->paf[1]->c, tiles1->paf[1]->c);
              //itv_print(diffc->coeff);
              itv_mul(pr1->itv, exp1->coeff, diffc->coeff, ptr->coeff);
              //ptmp1[j]=t1p_aff_alloc_init(pr1);
              ptmp2[i]=exp1;
              printf("%zu",i);
              itv_print(ptmp2[i]->coeff);
              printf("%zu",i);
           }
           else {
              p=tiles2->paf[1]->q;
              p = p->n;
              //itv_print(ptmp1[0]->coeff);
              itv_mul(pr1->itv, ptr->coeff, coeff, p->coeff);
              //itv_print(ptr->coeff);
              itv_sub(diffc->coeff, tiles2->paf[0]->c, tiles1->paf[0]->c);
              //itv_print(diffc->coeff);
              itv_mul(pr1->itv, exp2->coeff, diffc->coeff, ptr->coeff);
              //ptmp2[j]=t1p_aff_alloc_init(pr1);
              ptmp2[i]=exp2;
              printf("%zu",i);
              itv_print(ptmp2[i]->coeff);
              printf("%zu",i);
           } 
         //itv_print(ptr->coeff);
               
}*/




//////// ******* end of test *********



/*t1p_aff_t *exp1 = tiles[6]->paf[0];
t1p_aff_t *exp2 = tiles[6]->paf[1];
void * ptmp = NULL;

ptmp = (void*)exp1;
exp1 = exp2;
exp2 = (t1p_aff_t*)ptmp;

tiles[6]->paf[0]=exp1;
tiles[6]->paf[1]=exp2;

itv_print(tiles[6]->paf[0]->q->coeff);
itv_print(tiles[6]->paf[1]->q->coeff);
itv_print(tiles[6]->paf[0]->q->n->coeff);
itv_print(tiles[6]->paf[1]->q->n->coeff);

t1p_aff_t* res = t1p_aff_alloc_init(pr1);

itv_print(tiles[6]->paf[0]->c);
itv_print(tiles[6]->paf[1]->c);

itv_print(tiles[5]->paf[0]->c);
itv_print(tiles[5]->paf[1]->c);

//for (i=0; i<tiles[6]->dims; i++){
   itv_sub(res->c, tiles[6]->paf[1]->c, tiles[5]->paf[1]->c);
   itv_print(res->c);*/
//}

//itv_print(tiles[6]->paf[0]->c);
//itv_print(tiles[6]->paf[1]->c);


/*sum1=0;
for j=1:length(Gnormal)
    %u(:,j)=Gnormal(:,j);%./norm(Gnormal(:,j),1);
    for i=1:length(c1)
        left(i)=(c2(i)-c1(i))*Gnormal(i,j);
    end
    lefteq(j)=abs(sum(left));*/





  //if (t1p_aff_is_perfectly_eq(pr, tiles[6]->paf[1], tiles0->paf[1])) printf("o");
  //else printf("n");


  //printf("%f",zon1.mat[0][3]);
  //printf("%f",zon1.mat[1][3]);

/*     zon1.row=2;
    zon1.column=5;
    zon1.totnocols=5;
    
     signsvec(&zon1);

     int i,j,k;
for(k=0;k<10;k++){
  for(i=0;i<tilings[k][k].row;i++){
    for(j=0;j<tilings[k][k].column;j++){
        printf("%f \t",tilings[k][k].mat[i][j]);}

printf("\n");}
}

for(k=0;k<tilings[0][0].row;k++){
printf("%f \t", tilings[0][0].center[k]);}

for(k=0;k<tilings[3][3].row;k++){
printf("%f \t", tilings[3][3].center[k]);}

ap_manager_t* man_tiles = t1p_manager_alloc();
    t1p_internal_t * pr1 = t1p_init_from_manager(man_tiles, AP_FUNID_UNKNOWN);



t1p_t* tiles0 = t1p_top(man_tiles, 0, 2);
t1p_t* tiles1 = t1p_top(man_tiles, 0, 2);
t1p_t* tiles2 = t1p_top(man_tiles, 0, 2);
t1p_t* tiles3 = t1p_top(man_tiles, 0, 2);
t1p_t* tiles4 = t1p_top(man_tiles, 0, 2);
t1p_t* tiles5 = t1p_top(man_tiles, 0, 2);
t1p_t* tiles6 = t1p_top(man_tiles, 0, 2);
t1p_t* tiles7 = t1p_top(man_tiles, 0, 2);
t1p_t* tiles8 = t1p_top(man_tiles, 0, 2);
t1p_t* tiles9 = t1p_top(man_tiles, 0, 2);


tiles0->paf[0] = t1p_aff_alloc_init(pr1);
num_set_double(num, tilings[0][0].center[0]);
    itv_set_num(coeff, num);
    itv_set(tiles0->paf[0]->c, coeff);
    num_set_double(num, tilings[0][0].mat[0][0]);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles0->paf[0], coeff, 0);
    num_set_double(num, tilings[0][0].mat[0][1]);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles0->paf[0], coeff, 1);
t1p_aff_reduce(pr1, tiles0->paf[0]);
tiles0->paf[1] = t1p_aff_alloc_init(pr1);
num_set_double(num, tilings[0][0].center[1]);
    itv_set_num(coeff, num);
    itv_set(tiles0->paf[1]->c, coeff);
    num_set_double(num, tilings[0][0].mat[1][0]);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles0->paf[1], coeff, 0);
    num_set_double(num, tilings[0][0].mat[1][1]);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles0->paf[1], coeff, 1);
t1p_aff_reduce(pr1, tiles0->paf[1]);

tiles1->paf[0] = t1p_aff_alloc_init(pr1);
num_set_double(num, tilings[1][1].center[0]);
    itv_set_num(coeff, num);
    itv_set(tiles1->paf[0]->c, coeff);
    num_set_double(num, tilings[1][1].mat[0][0]);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles1->paf[0], coeff, 0);
    num_set_double(num, tilings[1][1].mat[0][1]);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles1->paf[0], coeff, 1);
t1p_aff_reduce(pr1, tiles1->paf[0]);
tiles1->paf[1] = t1p_aff_alloc_init(pr1);
num_set_double(num, tilings[1][1].center[1]);
    itv_set_num(coeff, num);
    itv_set(tiles1->paf[1]->c, coeff);
    num_set_double(num, tilings[1][1].mat[1][0]);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles1->paf[1], coeff, 0);
    num_set_double(num, tilings[1][1].mat[1][1]);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles1->paf[1], coeff, 1);
t1p_aff_reduce(pr1, tiles1->paf[1]);

tiles2->paf[0] = t1p_aff_alloc_init(pr1);
num_set_double(num, tilings[2][2].center[0]);
    itv_set_num(coeff, num);
    itv_set(tiles2->paf[0]->c, coeff);
    num_set_double(num, tilings[2][2].mat[0][0]);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles2->paf[0], coeff, 0);
    num_set_double(num, tilings[2][2].mat[0][1]);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles2->paf[0], coeff, 1);
t1p_aff_reduce(pr1, tiles2->paf[0]);
tiles2->paf[1] = t1p_aff_alloc_init(pr1);
num_set_double(num, tilings[2][2].center[1]);
    itv_set_num(coeff, num);
    itv_set(tiles2->paf[1]->c, coeff);
    num_set_double(num, tilings[2][2].mat[1][0]);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles2->paf[1], coeff, 0);
    num_set_double(num, tilings[2][2].mat[1][1]);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles2->paf[1], coeff, 1);
t1p_aff_reduce(pr1, tiles2->paf[1]);

tiles3->paf[0] = t1p_aff_alloc_init(pr1);
num_set_double(num, tilings[3][3].center[0]);
    itv_set_num(coeff, num);
    itv_set(tiles3->paf[0]->c, coeff);
    num_set_double(num, tilings[3][3].mat[0][0]);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles3->paf[0], coeff, 0);
    num_set_double(num, tilings[3][3].mat[0][1]);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles3->paf[0], coeff, 1);
t1p_aff_reduce(pr1, tiles3->paf[0]);
tiles3->paf[1] = t1p_aff_alloc_init(pr1);
num_set_double(num, tilings[3][3].center[1]);
    itv_set_num(coeff, num);
    itv_set(tiles3->paf[1]->c, coeff);
    num_set_double(num, tilings[3][3].mat[1][0]);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles3->paf[1], coeff, 0);
    num_set_double(num, tilings[3][3].mat[1][1]);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles3->paf[1], coeff, 1);
t1p_aff_reduce(pr1, tiles3->paf[1]);

tiles4->paf[0] = t1p_aff_alloc_init(pr1);
num_set_double(num, tilings[4][4].center[0]);
    itv_set_num(coeff, num);
    itv_set(tiles4->paf[0]->c, coeff);
    num_set_double(num, tilings[4][4].mat[0][0]);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles4->paf[0], coeff, 0);
    num_set_double(num, tilings[4][4].mat[0][1]);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles4->paf[0], coeff, 1);
t1p_aff_reduce(pr1, tiles4->paf[0]);
tiles4->paf[1] = t1p_aff_alloc_init(pr1);
num_set_double(num, tilings[4][4].center[1]);
    itv_set_num(coeff, num);
    itv_set(tiles4->paf[1]->c, coeff);
    num_set_double(num, tilings[4][4].mat[1][0]);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles4->paf[1], coeff, 0);
    num_set_double(num, tilings[4][4].mat[1][1]);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles4->paf[1], coeff, 1);
t1p_aff_reduce(pr1, tiles4->paf[1]);

tiles5->paf[0] = t1p_aff_alloc_init(pr1);
num_set_double(num, tilings[5][5].center[0]);
    itv_set_num(coeff, num);
    itv_set(tiles5->paf[0]->c, coeff);
    num_set_double(num, tilings[5][5].mat[0][0]);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles5->paf[0], coeff, 0);
    num_set_double(num, tilings[5][5].mat[0][1]);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles5->paf[0], coeff, 1);
t1p_aff_reduce(pr1, tiles5->paf[0]);
tiles5->paf[1] = t1p_aff_alloc_init(pr1);
num_set_double(num, tilings[5][5].center[1]);
    itv_set_num(coeff, num);
    itv_set(tiles5->paf[1]->c, coeff);
    num_set_double(num, tilings[5][5].mat[1][0]);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles5->paf[1], coeff, 0);
    num_set_double(num, tilings[5][5].mat[1][1]);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles5->paf[1], coeff, 1);
t1p_aff_reduce(pr1, tiles5->paf[1]);


tiles6->paf[0] = t1p_aff_alloc_init(pr1);
num_set_double(num, tilings[6][6].center[0]);
    itv_set_num(coeff, num);
    itv_set(tiles6->paf[0]->c, coeff);
    num_set_double(num, tilings[6][6].mat[0][0]);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles6->paf[0], coeff, 0);
    num_set_double(num, tilings[6][6].mat[0][1]);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles6->paf[0], coeff, 1);
t1p_aff_reduce(pr1, tiles6->paf[0]);
tiles6->paf[1] = t1p_aff_alloc_init(pr1);
num_set_double(num, tilings[6][6].center[1]);
    itv_set_num(coeff, num);
    itv_set(tiles6->paf[1]->c, coeff);
    num_set_double(num, tilings[6][6].mat[1][0]);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles6->paf[1], coeff, 0);
    num_set_double(num, tilings[6][6].mat[1][1]);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles6->paf[1], coeff, 1);
t1p_aff_reduce(pr1, tiles6->paf[1]);


tiles7->paf[0] = t1p_aff_alloc_init(pr1);
num_set_double(num, tilings[7][7].center[0]);
    itv_set_num(coeff, num);
    itv_set(tiles7->paf[0]->c, coeff);
    num_set_double(num, tilings[7][7].mat[0][0]);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles7->paf[0], coeff, 0);
    num_set_double(num, tilings[7][7].mat[0][1]);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles7->paf[0], coeff, 1);
t1p_aff_reduce(pr1, tiles7->paf[0]);
tiles7->paf[1] = t1p_aff_alloc_init(pr1);
num_set_double(num, tilings[7][7].center[1]);
    itv_set_num(coeff, num);
    itv_set(tiles7->paf[1]->c, coeff);
    num_set_double(num, tilings[7][7].mat[1][0]);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles7->paf[1], coeff, 0);
    num_set_double(num, tilings[7][7].mat[1][1]);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles7->paf[1], coeff, 1);
t1p_aff_reduce(pr1, tiles7->paf[1]);


tiles8->paf[0] = t1p_aff_alloc_init(pr1);
num_set_double(num, tilings[8][8].center[0]);
    itv_set_num(coeff, num);
    itv_set(tiles8->paf[0]->c, coeff);
    num_set_double(num, tilings[8][8].mat[0][0]);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles8->paf[0], coeff, 0);
    num_set_double(num, tilings[8][8].mat[0][1]);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles8->paf[0], coeff, 1);
t1p_aff_reduce(pr1, tiles8->paf[0]);
tiles8->paf[1] = t1p_aff_alloc_init(pr1);
num_set_double(num, tilings[8][8].center[1]);
    itv_set_num(coeff, num);
    itv_set(tiles8->paf[1]->c, coeff);
    num_set_double(num, tilings[0][0].mat[1][0]);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles8->paf[1], coeff, 0);
    num_set_double(num, tilings[0][0].mat[1][1]);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles8->paf[1], coeff, 1);
t1p_aff_reduce(pr1, tiles8->paf[1]);


tiles9->paf[0] = t1p_aff_alloc_init(pr1);
num_set_double(num, tilings[9][9].center[0]);
    itv_set_num(coeff, num);
    itv_set(tiles9->paf[0]->c, coeff);
    num_set_double(num, tilings[9][9].mat[0][0]);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles9->paf[0], coeff, 0);
    num_set_double(num, tilings[9][9].mat[0][1]);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles9->paf[0], coeff, 1);
t1p_aff_reduce(pr1, tiles9->paf[0]);
tiles9->paf[1] = t1p_aff_alloc_init(pr1);
num_set_double(num, tilings[9][9].center[1]);
    itv_set_num(coeff, num);
    itv_set(tiles9->paf[1]->c, coeff);
    num_set_double(num, tilings[9][9].mat[1][0]);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles9->paf[1], coeff, 0);
    num_set_double(num, tilings[9][9].mat[1][1]);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles9->paf[1], coeff, 1);
t1p_aff_reduce(pr1, tiles9->paf[1]);


/*tiles1->paf[0] = t1p_aff_alloc_init(pr1);
num_set_double(num, 24);
    itv_set_num(coeff, num);
    itv_set(tiles1->paf[0]->c, coeff);
    num_set_double(num, 1);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles1->paf[0], coeff, 0);
    num_set_double(num, 3);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles1->paf[0], coeff, 1);
t1p_aff_reduce(pr1, tiles1->paf[0]);
tiles1->paf[1] = t1p_aff_alloc_init(pr1);
num_set_double(num, 8);
    itv_set_num(coeff, num);
    itv_set(tiles1->paf[1]->c, coeff);
    num_set_double(num, 1);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles1->paf[1], coeff, 0);
    num_set_double(num, 5);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles1->paf[1], coeff, 1);
t1p_aff_reduce(pr1, tiles1->paf[1]);

tiles2->paf[0] = t1p_aff_alloc_init(pr1);
num_set_double(num, 14);
    itv_set_num(coeff, num);
    itv_set(tiles2->paf[0]->c, coeff);
    num_set_double(num, 5);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles2->paf[0], coeff, 0);
    num_set_double(num, 3);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles2->paf[0], coeff, 1);
t1p_aff_reduce(pr1, tiles2->paf[0]);
tiles2->paf[1] = t1p_aff_alloc_init(pr1);
num_set_double(num, 5);
    itv_set_num(coeff, num);
    itv_set(tiles2->paf[1]->c, coeff);
    num_set_double(num, 2);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles2->paf[1], coeff, 0);
    num_set_double(num, 5);
    itv_set_num(coeff, num);
    t1p_aff_build(pr1, tiles2->paf[1], coeff, 1);
t1p_aff_reduce(pr1, tiles2->paf[1]);


bool res = true;

for (i=0; i<tiles1->dims; i++) {
	    if (tiles1->paf[i] != tiles2->paf[i] && itv_is_eq(tiles1->box[i], tiles2->box[i])) {
		res = t1p_aff_is_eq(pr1, tiles1->paf[i], tiles2->paf[i]);
		if (!res) break;
	    } else {
		res = false;
		break;
	    }
	}
    
printf("%d\n", res);

if (t1p_is_leq(man_tiles,tiles1,tiles2)) printf("o");
else
printf("n");*/


/*bool res = true;

for (i=0; i<tiles8->dims; i++) {
	    if (tiles8->paf[i] != tiles9->paf[i] && itv_is_eq(tiles8->box[i], tiles9->box[i])) {
		res = t1p_aff_is_eq(pr1, tiles8->paf[i], tiles9->paf[i]);
		if (!res) break;
	    } else {
		res = false;
		break;
	    }
	}
    
printf("%d\n", res);*/


}

t1p_aff_t *
afexpr_random_cst (t1p_internal_t * pr, unsigned int seed)
{
  itv_t itv;
  num_t num_dbl;
  double dbl;
  t1p_aff_t *afexpr = t1p_aff_alloc_init (pr);
  srand (seed);
  dbl = (double) ((rand () - rand ()) / (1.0 * rand ()));
  num_set_double (num_dbl, dbl);
  itv_set_num (itv, num_dbl);
  itv_set (afexpr->c, itv);
  itv_set (afexpr->itv, itv);
  return afexpr;
}

t1p_aff_t *
afexpr_random_nb (t1p_internal_t * pr, unsigned int seed, size_t Nb)
{
  itv_t itv;
  itv_init (itv);
  num_t num_dbl;
  num_init (num_dbl);
  double dbl;
  size_t tmp;
  t1p_aff_t *afexpr = t1p_aff_alloc_init (pr);
  dbl = (double) ((rand () - rand ()) / (1.0 * rand ()));
  num_set_double (num_dbl, dbl);
  itv_set_num (itv, num_dbl);
  itv_set (afexpr->c, itv);
  itv_set (afexpr->itv, itv);
  srand (seed);
  while (afexpr->l < Nb) {
    dbl = (double) ((rand () - rand ()) / (1.0 * rand ()));
    num_set_double (num_dbl, dbl);
    itv_set_num (itv, num_dbl);
    if (afexpr->l < Nb - 3)
      t1p_aff_nsym_create (pr, afexpr, itv, IN);
    else
      t1p_aff_nsym_create (pr, afexpr, itv, UN);
  }
  itv_clear (itv);
  num_clear (num_dbl);
  return afexpr;
}

/* multiplication */
/* 
   - 0 x random = 0
   - bottom x random = bottom
   - top x random = top
   - exp1 x exp2 = exp3 (compared with Fluctuat)
 */
void
test_mul_non_constrained (t1p_internal_t * pr, t1p_t * env)
{
  printf ("### aff x aff ###\n");
  t1p_aff_check_free (pr, env->paf[0]);
  t1p_aff_check_free (pr, env->paf[1]);
  t1p_aff_check_free (pr, env->paf[2]);
  itv_t coeff;
  itv_init (coeff);
  num_t num;
  num_init (num);
  env->paf[0] = t1p_aff_alloc_init (pr);
  long int seed48 = time (NULL);
  env->paf[1] = afexpr_random_nb (pr, (unsigned int) (seed48 % UINT_MAX), 10);
  /* 0 x random = 0 */
  env->paf[2] = t1p_aff_mul (pr, env->paf[0], env->paf[1], env);
  if (t1p_aff_is_zero (pr, env->paf[2]))
    printf ("o");
  else
    printf ("n");
  t1p_aff_check_free (pr, env->paf[2]);
  /* bottom x random = bottom */
  env->paf[2] = t1p_aff_mul (pr, pr->bot, env->paf[1], env);
  if (t1p_aff_is_bottom (pr, env->paf[2]))
    printf ("o");
  else
    printf ("n");
  t1p_aff_check_free (pr, env->paf[2]);
  /* top x random = top */
  env->paf[2] = t1p_aff_mul (pr, pr->top, env->paf[1], env);
  if (t1p_aff_is_top (pr, env->paf[2]))
    printf ("o");
  else
    printf ("n");
  t1p_aff_check_free (pr, env->paf[2]);
  /* exp1 = -2.5 + .5\eps1 + \eps2 ; [-4,-1]
     exp2 = exp1 x exp1
     (fluctuat) exp2 = 6.875 - 2.5\eps1 -5\eps2 + 1.625\eps_f ; [1,16]
     (t1+) exp2 = 6.25 - 2.5\eps1 -5\eps2 + [0,2.25]  ; [1,16]
     = 7.375  - 2.5\eps1 -5\eps2 + 1.125\eps_f        ; [1,16]
   */
  t1p_aff_check_free (pr, env->paf[1]);
  env->paf[1] = t1p_aff_alloc_init (pr);
  num_set_double (num, -2.5);
  itv_set_num (coeff, num);
  itv_set (env->paf[1]->c, coeff);
  num_set_double (num, 0.5);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, env->paf[1], coeff, 0);
  num_set_double (num, 1);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, env->paf[1], coeff, 1);
  itv_set_int2 (env->paf[1]->itv, -4, -1);
  env->paf[2] = t1p_aff_mul (pr, env->paf[1], env->paf[1], env);
  t1p_aff_reduce (pr, env->paf[2]);

  t1p_aff_t *aff_exact = t1p_aff_alloc_init (pr);
  num_set_double (num, 7.375);
  itv_set_num (coeff, num);
  itv_set (aff_exact->c, coeff);
  num_set_double (num, -2.5);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, aff_exact, coeff, 0);
  num_set_double (num, -5);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, aff_exact, coeff, 1);
  num_set_double (num, 1.125);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, aff_exact, coeff, (-1 + pr->dim));
  itv_set_int2 (aff_exact->itv, 1, 16);
  if (t1p_aff_is_perfectly_eq (pr, env->paf[2], aff_exact))
    printf ("o");
  else
    printf ("n");
  t1p_aff_check_free (pr, env->paf[2]);

  /* exp1 = -2.5 + .5\eps1 + \eps2 ; [-4,-1]
     exp0 = 1.5 + eps1 - 0.5\eps3 ; [0,3]
     exp2 = exp1 x exp0
     (fluctuat, t1+) exp2 = -3.5 - 1.75\eps1 + 1.5\eps2 + 1.25\eps3 + 2\eps_f ; [-10,0]
   */
  t1p_aff_check_free (pr, env->paf[0]);
  env->paf[0] = t1p_aff_alloc_init (pr);
  num_set_double (num, 1.5);
  itv_set_num (coeff, num);
  itv_set (env->paf[0]->c, coeff);
  num_set_double (num, 1.0);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, env->paf[0], coeff, 0);
  num_set_double (num, -0.5);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, env->paf[0], coeff, 2);
  itv_set_int2 (env->paf[0]->itv, 0, 3);
  env->paf[2] = t1p_aff_mul (pr, env->paf[0], env->paf[1], env);
  t1p_aff_reduce (pr, env->paf[2]);

  t1p_aff_t *aff_exact1 = t1p_aff_alloc_init (pr);
  num_set_double (num, -3.5);
  itv_set_num (coeff, num);
  itv_set (aff_exact1->c, coeff);
  num_set_double (num, -1.75);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, aff_exact1, coeff, 0);
  num_set_double (num, 1.5);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, aff_exact1, coeff, 1);
  num_set_double (num, 1.25);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, aff_exact1, coeff, 2);
  num_set_double (num, 2);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, aff_exact1, coeff, (-1 + pr->dim));
  itv_set_int2 (aff_exact1->itv, -10, 0);

  if (t1p_aff_is_perfectly_eq (pr, env->paf[2], aff_exact1))
    printf ("o");
  else
    printf ("n");
  t1p_aff_check_free (pr, env->paf[2]);

  /* exp1 = -2.5 + .5\eps1 + \eps2 ; [-4,-1]
     lambda = [1,1]
     exp2 = lambda x exp1
     exp2 = exp1
   */
  itv_t lambda;
  itv_init (lambda);
  itv_set_int (lambda, 1);
  env->paf[2] = t1p_aff_mul_itv (pr, env->paf[1], lambda);
  if (t1p_aff_is_perfectly_eq (pr, env->paf[2], env->paf[1]))
    printf ("o");
  else
    printf ("n");

  /* exp1 = 2 + eps0 - eps1 + 1.5eps3 ; [-1.5,5.5]
     exp2 = -1 + eps0 + 0.5eps2 + 0.75eps3 + eps4 ; [1.25,5.25]
     exp3 = exp1*exp2 
     = (fluctuat) -0.9375 + eps0 + eps1 + eps2 + 2eps4 + 10.3125epsf; [-16.25,12.375]
   */
  t1p_aff_check_free (pr, env->paf[0]);
  t1p_aff_check_free (pr, env->paf[1]);
  t1p_aff_check_free (pr, env->paf[2]);

  env->paf[0] = t1p_aff_alloc_init (pr);
  num_set_double (num, 2);
  itv_set_num (coeff, num);
  itv_set (env->paf[0]->c, coeff);
  num_set_double (num, 1.0);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, env->paf[0], coeff, 0);
  num_set_double (num, -1.0);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, env->paf[0], coeff, 1);
  num_set_double (num, 1.5);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, env->paf[0], coeff, 3);
  num_set_double (num, 1.5);
  bound_set_num (env->paf[0]->itv->inf, num);
  num_set_double (num, 5.5);
  bound_set_num (env->paf[0]->itv->sup, num);

  env->paf[1] = t1p_aff_alloc_init (pr);
  num_set_double (num, -1);
  itv_set_num (coeff, num);
  itv_set (env->paf[1]->c, coeff);
  num_set_double (num, 1.0);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, env->paf[1], coeff, 0);
  num_set_double (num, 0.5);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, env->paf[1], coeff, 2);
  num_set_double (num, 0.75);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, env->paf[1], coeff, 3);
  num_set_double (num, 1);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, env->paf[1], coeff, 4);
  num_set_double (num, 4.25);
  bound_set_num (env->paf[1]->itv->inf, num);
  num_set_double (num, 2.25);
  bound_set_num (env->paf[1]->itv->sup, num);
  env->paf[2] = t1p_aff_mul (pr, env->paf[0], env->paf[1], env);
  t1p_aff_reduce (pr, env->paf[2]);

  t1p_aff_check_free (pr, aff_exact);
  aff_exact = t1p_aff_alloc_init (pr);
  num_set_double (num, -0.9375);
  itv_set_num (coeff, num);
  itv_set (aff_exact->c, coeff);
  num_set_double (num, 1.);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, aff_exact, coeff, 0);
  num_set_double (num, 1.);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, aff_exact, coeff, 1);
  num_set_double (num, 1.);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, aff_exact, coeff, 2);
  num_set_double (num, 2);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, aff_exact, coeff, 4);
  num_set_double (num, 10.3125);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, aff_exact, coeff, (-1 + pr->dim));
  num_set_double (num, 16.25);
  bound_set_num (aff_exact->itv->inf, num);
  num_set_double (num, 12.375);
  bound_set_num (aff_exact->itv->sup, num);

  if (t1p_aff_is_perfectly_eq (pr, env->paf[2], aff_exact))
    printf ("o");
  else
    printf ("n");

  num_clear (num);
  itv_clear (coeff);
  itv_clear (lambda);
  t1p_aff_free (pr, aff_exact);
  t1p_aff_free (pr, aff_exact1);
  printf ("\n");
}

void
test_mul_constrained (t1p_internal_t * pr, t1p_t * env)
{
  printf ("### aff x aff (constrained) ###\n");
  t1p_aff_check_free (pr, env->paf[0]);
  t1p_aff_check_free (pr, env->paf[1]);
  t1p_aff_check_free (pr, env->paf[2]);
  itv_t coeff;
  itv_init (coeff);
  num_t num;
  num_init (num);

  env->paf[0] = t1p_aff_alloc_init (pr);
  num_set_double (num, 1.094658404232323301);
  itv_set_num (coeff, num);
  itv_set (env->paf[0]->c, coeff);
  num_set_double (num, -0.48202531945707516314);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, env->paf[0], coeff, 0);
  num_set_double (num, 0.93829687951457496631);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, env->paf[0], coeff, 1);
  num_set_double (num, -0.33333333333333331483);
  bound_set_num (env->paf[0]->itv->inf, num);
  num_set_double (num, 2.5149806032039729864);
  bound_set_num (env->paf[0]->itv->sup, num);
  /*
     num_set_double(num, 1);
     itv_set_num(coeff, num);
     itv_set(env->paf[0]->c, coeff);
     num_set_double(num, -0.5);
     itv_set_num(coeff, num);
     t1p_aff_build(pr, env->paf[0], coeff, 0);
     num_set_double(num, 1);
     itv_set_num(coeff, num);
     t1p_aff_build(pr, env->paf[0], coeff, 1);
     num_set_double(num, 0.0);
     bound_set_num(env->paf[0]->itv->inf,num);
     num_set_double(num, 2.5);
     bound_set_num(env->paf[0]->itv->sup,num);
   */

  t1p_aff_t *aff = t1p_aff_alloc_init (pr);
  itv_set_int (coeff, -1);
  t1p_aff_build (pr, aff, coeff, 0);
  ap_linexpr0_t *linexpr0 = t1p_ap_linexpr0_set_aff (pr, aff, env);
  ap_lincons0_t lincons0;
  lincons0.constyp = AP_CONS_SUPEQ;
  lincons0.linexpr0 = linexpr0;
  lincons0.scalar = NULL;
  ap_lincons0_array_t eps_lincons_array;
  eps_lincons_array.size = 1;
  eps_lincons_array.p = &lincons0;
  ap_abstract0_meet_lincons_array (pr->manNS, true, env->abs, &eps_lincons_array);
  t1p_update_nsymcons_gamma (pr, env);

  env->paf[2] = t1p_aff_mul (pr, env->paf[0], env->paf[0], env);
  t1p_aff_reduce (pr, env->paf[2]);

  t1p_aff_t *aff_exact1 = t1p_aff_alloc_init (pr);
  num_set_double (num, 1.8355754144819342866);
  itv_set_num (coeff, num);
  itv_set (aff_exact1->c, coeff);
  num_set_double (num, -1.2876545425906109621);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, aff_exact1, coeff, 0);
  num_set_double (num, 2.5065119827447759349);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, aff_exact1, coeff, 1);
  num_set_double (num, 0.69538549467490184952);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, aff_exact1, coeff, (-1 + pr->dim));
  num_set_double (num, -0.11111111111111109107);
  bound_set_num (aff_exact1->itv->inf, num);
  num_set_double (num, 6.3251274344922201465);
  bound_set_num (aff_exact1->itv->sup, num);
  if (t1p_aff_is_perfectly_eq (pr, env->paf[2], aff_exact1))
    printf ("o");
  if (t1p_aff_is_tild_eq (pr, env->paf[2], aff_exact1))
    printf ("O");
  else
    printf ("n");

  t1p_aff_check_free (pr, env->paf[0]);
  t1p_aff_check_free (pr, env->paf[2]);

  env->paf[0] = t1p_aff_alloc_init (pr);
  num_set_double (num, 2.5);
  itv_set_num (coeff, num);
  itv_set (env->paf[0]->c, coeff);
  num_set_double (num, 20);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, env->paf[0], coeff, 0);
  num_set_double (num, 2.5);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, env->paf[0], coeff, 1);
  num_set_double (num, 20);
  bound_set_num (env->paf[0]->itv->inf, num);
  num_set_double (num, 5);
  bound_set_num (env->paf[0]->itv->sup, num);

  env->paf[1] = t1p_aff_alloc_init (pr);
  num_set_double (num, -1.25);
  itv_set_num (coeff, num);
  itv_set (env->paf[1]->c, coeff);
  num_set_double (num, 10);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, env->paf[1], coeff, 0);
  num_set_double (num, 1.25);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, env->paf[1], coeff, 2);
  num_set_double (num, 12.5);
  bound_set_num (env->paf[1]->itv->inf, num);
  num_set_double (num, 0);
  bound_set_num (env->paf[1]->itv->sup, num);

  env->paf[2] = t1p_aff_mul (pr, env->paf[0], env->paf[1], env);
  t1p_aff_reduce (pr, env->paf[2]);

  t1p_aff_t *aff_exact = t1p_aff_alloc_init (pr);
  num_set_double (num, -28.125);
  itv_set_num (coeff, num);
  itv_set (aff_exact->c, coeff);
  num_set_double (num, -200);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, aff_exact, coeff, 0);
  num_set_double (num, -15.625);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, aff_exact, coeff, 1);
  num_set_double (num, -9.375);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, aff_exact, coeff, 2);
  num_set_double (num, 53.125);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, aff_exact, coeff, (-1 + pr->dim));
  num_set_double (num, 62.5);
  bound_set_num (aff_exact->itv->inf, num);
  num_set_double (num, 250);
  bound_set_num (aff_exact->itv->sup, num);
  if (t1p_aff_is_perfectly_eq (pr, env->paf[2], aff_exact))
    printf ("o");
  else
    printf ("n");

  num_clear (num);
  itv_clear (coeff);
  //itv_clear(lambda);
  t1p_aff_free (pr, aff_exact);
  t1p_aff_free (pr, aff);
  ap_lincons0_clear (&lincons0);
  t1p_aff_free (pr, aff_exact1);
  printf ("\n");
}


/* division */
/* 
   - random / 0 = top
   - random / 1 = random
   - exp1 / exp2 = exp3 (compared with Fluctuat)
   - exp1 / [1, +oo] = exp4 (compared with Fluctuat)
   - exp1 / [-oo, -1] = exp5 (compared with Fluctuat)
 */

void
test_div_non_constrained (t1p_internal_t * pr, t1p_t * env)
{
  printf ("### aff / aff ###\n");
  t1p_aff_check_free (pr, env->paf[0]);
  t1p_aff_check_free (pr, env->paf[1]);
  t1p_aff_check_free (pr, env->paf[2]);
  itv_t coeff, lambda;
  itv_init (coeff);
  itv_init (lambda);
  num_t num;
  num_init (num);
  env->paf[0] = t1p_aff_alloc_init (pr);
  long int seed48 = time (NULL);
  env->paf[1] = afexpr_random_nb (pr, (unsigned int) (seed48 % UINT_MAX), 10);
  /* random / top = top */
  env->paf[2] = t1p_aff_div (pr, env->paf[1], pr->top, env);
  if (t1p_aff_is_top (pr, env->paf[2]))
    printf ("o");
  else
    printf ("n");
  t1p_aff_check_free (pr, env->paf[2]);

  /* random / 0 = top */
  env->paf[2] = t1p_aff_div (pr, env->paf[1], env->paf[0], env);
  if (t1p_aff_is_top (pr, env->paf[2]))
    printf ("o");
  else
    printf ("n");
  t1p_aff_check_free (pr, env->paf[2]);

  /* randim / 1 = random */
  t1p_aff_t *dummy = t1p_aff_alloc_init (pr);
  itv_set_int (dummy->c, 1);
  itv_set_int (dummy->itv, 1);
  env->paf[2] = t1p_aff_div (pr, env->paf[1], dummy, env);
  if (t1p_aff_is_perfectly_eq (pr, env->paf[1], env->paf[2]))
    printf ("o");
  else
    printf ("n");
  t1p_aff_check_free (pr, env->paf[2]);

  /* exp1 = -2.5 + .5\eps1 + \eps2 ; [-4,-1]
     exp2 = exp1 / [1, +oo]
     (fluctuat) exp2 = -1.25 + 2.75\eps_U ; [-4,0]
     (t1+) exp2 = -1.25 + 0.25\eps1 + 0.5\eps2 + 2\eps_f ; [-4,0]
   */
  itv_set_top (dummy->c);
  bound_set_int (dummy->c->inf, -1);
  itv_set (dummy->itv, dummy->c);
  t1p_aff_check_free (pr, env->paf[1]);
  env->paf[1] = t1p_aff_alloc_init (pr);
  num_set_double (num, -2.5);
  itv_set_num (coeff, num);
  itv_set (env->paf[1]->c, coeff);
  num_set_double (num, 0.5);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, env->paf[1], coeff, 0);
  num_set_double (num, 1);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, env->paf[1], coeff, 1);
  itv_set_int2 (env->paf[1]->itv, -4, -1);
  env->paf[2] = t1p_aff_div (pr, env->paf[1], dummy, env);
  t1p_aff_reduce (pr, env->paf[2]);

  t1p_aff_t *aff_exact = t1p_aff_alloc_init (pr);
  num_set_double (num, -1.25);
  itv_set_num (coeff, num);
  itv_set (aff_exact->c, coeff);
  num_set_double (num, 0.25);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, aff_exact, coeff, 0);
  num_set_double (num, 0.5);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, aff_exact, coeff, 1);
  num_set_double (num, 2);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, aff_exact, coeff, (-1 + pr->dim));
  itv_set_int2 (aff_exact->itv, -4, 0);
  if (t1p_aff_is_perfectly_eq (pr, env->paf[2], aff_exact))
    printf ("o");
  else
    printf ("n");
  t1p_aff_check_free (pr, env->paf[2]);

  /* exp1 = -2.5 + .5\eps1 + \eps2 ; [-4,-1]
     exp0 = 1.5 + eps1 - 0.5\eps3 ; [0,3]
     exp2 = exp0 / exp1
     (fluctuat, t1+) exp2 = -0.91 - 0.7\eps1 - 0.24\eps2 + 0.29\eps3 + 0.86\eps_f ; [-3,0]
   */
  num_set_double (num, 1.5);
  itv_set_num (coeff, num);
  itv_set (env->paf[0]->c, coeff);
  num_set_double (num, 1.0);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, env->paf[0], coeff, 0);
  num_set_double (num, -0.5);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, env->paf[0], coeff, 2);
  itv_set_int2 (env->paf[0]->itv, 0, 3);
  env->paf[2] = t1p_aff_div (pr, env->paf[0], env->paf[1], env);
  t1p_aff_reduce (pr, env->paf[2]);

  t1p_aff_t *aff_exact1 = t1p_aff_alloc_init (pr);
  num_set_double (num, -0.91);
  itv_set_num (coeff, num);
  itv_set (aff_exact1->c, coeff);
  num_set_double (num, -0.7);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, aff_exact1, coeff, 0);
  num_set_double (num, -0.24);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, aff_exact1, coeff, 1);
  num_set_double (num, 0.29);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, aff_exact1, coeff, 2);
  num_set_double (num, 0.86);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, aff_exact1, coeff, (-1 + pr->dim));
  itv_set_int2 (aff_exact1->itv, -3, 0);
  if (t1p_aff_is_perfectly_eq (pr, env->paf[2], aff_exact1))
    printf ("o");
  else if (t1p_aff_is_tild_eq (pr, env->paf[2], aff_exact1))
    printf ("O");
  else
    printf ("n");

  num_clear (num);
  itv_clear (coeff);
  itv_clear (lambda);
  t1p_aff_free (pr, dummy);
  t1p_aff_free (pr, aff_exact);
  t1p_aff_free (pr, aff_exact1);
  printf ("\n");
}

/* sqrt */
/*
 - sqrt([-1,0]) = bottom
 - sqrt(bot) = bot
 - sqrt(top) = top
 - exp2 = sqrt(exp0) (compared with fluctuat)
 */
void
test_sqrt_non_constrained (t1p_internal_t * pr, t1p_t * env)
{
  printf ("### sqrt(aff) ###\n");
  t1p_aff_check_free (pr, env->paf[0]);
  t1p_aff_check_free (pr, env->paf[1]);
  t1p_aff_check_free (pr, env->paf[2]);
  itv_t coeff, lambda;
  itv_init (coeff);
  itv_init (lambda);
  num_t num;
  num_init (num);
  env->paf[0] = t1p_aff_alloc_init (pr);
  /* sqrt(0) = 0 */
  env->paf[2] = t1p_aff_sqrt (pr, env->paf[0], env);
  if (t1p_aff_is_zero (pr, env->paf[2]))
    printf ("o");
  else
    printf ("n");
  t1p_aff_check_free (pr, env->paf[2]);
  /* exp1 < 0; sqrt(exp1) = bottom */
  env->paf[1] = t1p_aff_alloc_init (pr);
  num_set_double (num, -2.5);
  itv_set_num (coeff, num);
  itv_set (env->paf[1]->c, coeff);
  num_set_double (num, 0.5);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, env->paf[1], coeff, 0);
  num_set_double (num, 1);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, env->paf[1], coeff, 1);
  itv_set_int2 (env->paf[1]->itv, -4, -1);
  env->paf[2] = t1p_aff_sqrt (pr, env->paf[1], env);
  t1p_aff_reduce (pr, env->paf[2]);
  if (t1p_aff_is_bottom (pr, env->paf[2]))
    printf ("o");
  else
    printf ("n");
  t1p_aff_check_free (pr, env->paf[2]);
  /* exp0 = -exp1 = 2.5 - .5\eps1 - \eps2 ; [1,4]
     exp2 = sqrt(exp0)
     (fluctuat, t1+) exp2 = 1.527 -0.158\eps1 -0.316\eps2 + 0.053\eps_f ; [1,2]
   */
  t1p_aff_check_free (pr, env->paf[0]);
  env->paf[0] = t1p_aff_neg (pr, env->paf[1]);
  env->paf[2] = t1p_aff_sqrt (pr, env->paf[0], env);
  t1p_aff_reduce (pr, env->paf[2]);

  t1p_aff_t *aff_exact = t1p_aff_alloc_init (pr);
  num_set_double (num, 1.5277402395547232);
  itv_set_num (coeff, num);
  itv_set (aff_exact->c, coeff);
  num_set_double (num, -0.1581138830084189);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, aff_exact, coeff, 0);
  num_set_double (num, -0.316227766016837);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, aff_exact, coeff, 1);
  num_set_double (num, 0.053398590529466938);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, aff_exact, coeff, (-1 + pr->dim));
  itv_set_int2 (aff_exact->itv, 1, 2);
  if (t1p_aff_is_perfectly_eq (pr, env->paf[2], aff_exact))
    printf ("o");
  else if (t1p_aff_is_tild_eq (pr, env->paf[2], aff_exact))
    printf ("O");
  else
    printf ("n");
  t1p_aff_check_free (pr, env->paf[2]);
  /* exp0 = 1.5 + 1\eps1 - 0.5\eps2 ; [0,3]
     exp2 = sqrt(exp0)
     (fluctuat, t1+) exp2 = 0.918 + 0.40\eps1 -0.20\eps2 + 0.30\eps_f ; [0,1.73...]
   */
  t1p_aff_check_free (pr, env->paf[0]);
  env->paf[0] = t1p_aff_alloc_init (pr);
  num_set_double (num, 1.5);
  itv_set_num (coeff, num);
  itv_set (env->paf[0]->c, coeff);
  num_set_double (num, 1.0);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, env->paf[0], coeff, 0);
  num_set_double (num, -0.5);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, env->paf[0], coeff, 2);
  itv_set_int2 (env->paf[0]->itv, 0, 3);
  env->paf[2] = t1p_aff_sqrt (pr, env->paf[0], env);
  t1p_aff_reduce (pr, env->paf[2]);

  t1p_aff_t *aff_exact1 = t1p_aff_alloc_init (pr);
  num_set_double (num, 0.91855865354369);
  itv_set_num (coeff, num);
  itv_set (aff_exact1->c, coeff);
  num_set_double (num, 0.408248290463863);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, aff_exact1, coeff, 0);
  num_set_double (num, -0.2041241452319315);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, aff_exact1, coeff, 1);
  num_set_double (num, 0.306186217847897);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, aff_exact1, coeff, (-1 + pr->dim));
  num_set_double (num, 1.7320508075688);
  itv_set_num (aff_exact1->itv, num);
  bound_set_int (aff_exact1->itv->inf, 0);
  if (t1p_aff_is_perfectly_eq (pr, env->paf[2], aff_exact1))
    printf ("o");
  else if (t1p_aff_is_tild_eq (pr, env->paf[2], aff_exact1))
    printf ("O");
  else
    printf ("n");
  /* sqrt(bottom) = bottom */
  t1p_aff_check_free (pr, env->paf[0]);
  t1p_aff_check_free (pr, env->paf[2]);
  env->paf[0] = pr->bot;
  env->paf[2] = t1p_aff_sqrt (pr, env->paf[0], env);
  if (t1p_aff_is_bottom (pr, env->paf[2]))
    printf ("o");
  else
    printf ("n");
  /* sqrt(top) = [0,+oo] + caveat */
  t1p_aff_check_free (pr, env->paf[0]);
  t1p_aff_check_free (pr, env->paf[2]);
  env->paf[0] = pr->top;
  env->paf[2] = t1p_aff_sqrt (pr, env->paf[0], env);
  itv_t o_oo;
  itv_init (o_oo);
  itv_set_top (o_oo);
  bound_set_int (o_oo->inf, 0);
  if (itv_is_leq (o_oo, env->paf[2]->itv))
    printf ("o");
  else
    printf ("n");

  num_clear (num);
  itv_clear (coeff);
  itv_clear (lambda);
  itv_clear (o_oo);
  t1p_aff_free (pr, aff_exact);
  t1p_aff_free (pr, aff_exact1);
  printf ("\n");
}

/* addition / subtraction */
// TODO tester que rand1 + rand2 = rand2 + rand1
void
test_pm_non_constrained (t1p_internal_t * pr, t1p_t * env)
{
  printf ("### aff +/- aff ###\n");
  t1p_aff_check_free (pr, env->paf[0]);
  t1p_aff_check_free (pr, env->paf[1]);
  t1p_aff_check_free (pr, env->paf[2]);
  itv_t coeff;
  itv_init (coeff);
  num_t num;
  num_init (num);
  long int seed48 = time (NULL);
  env->paf[0] = afexpr_random_nb (pr, (unsigned int) (seed48 % UINT_MAX), 10);
  /* random - random  = 0 */
  env->paf[2] = t1p_aff_sub (pr, env->paf[0], env->paf[0], env);
  if (t1p_aff_is_zero (pr, env->paf[2]))
    printf ("o");
  else
    printf ("n");
  t1p_aff_check_free (pr, env->paf[2]);

  env->paf[1] = t1p_aff_neg (pr, env->paf[0]);
  env->paf[2] = t1p_aff_add (pr, env->paf[0], env->paf[1], env);
  if (t1p_aff_is_zero (pr, env->paf[2]))
    printf ("o");
  else
    printf ("n");

  num_clear (num);
  itv_clear (coeff);
  //itv_clear(lambda);
  //t1p_aff_free(pr, aff_exact);
  //t1p_aff_free(pr, aff_exact1);
  printf ("\n");
}

/* some texpr */
void
test_eval_texpr_non_constrained (t1p_internal_t * pr, t1p_t * abs)
{
  printf ("### eval texp ###\n");
  t1p_aff_check_free (pr, abs->paf[0]);
//    t1p_aff_check_free(pr, abs->paf[1]);
  t1p_aff_check_free (pr, abs->paf[2]);
  itv_t coeff, lambda;
  itv_init (coeff);
  itv_init (lambda);
  num_t num;
  num_init (num);
  /* x = 0
     y = sqrt(x)*sqrt(x) 
     y = 0
   */
  abs->paf[0] = t1p_aff_alloc_init (pr);
  abs->paf[0]->pby++;
  ap_texpr0_t *texp = ap_texpr0_binop (AP_TEXPR_MUL,
				       ap_texpr0_unop (AP_TEXPR_SQRT,
						       ap_texpr0_dim (0),
						       AP_RTYPE_DOUBLE, AP_RDIR_UP),
				       ap_texpr0_unop (AP_TEXPR_SQRT,
						       ap_texpr0_dim (0),
						       AP_RTYPE_DOUBLE, AP_RDIR_UP),
				       AP_RTYPE_DOUBLE, AP_RDIR_UP);
  abs->paf[2] = t1p_aff_eval_ap_texpr0 (pr, texp, abs);
  if (t1p_aff_is_zero (pr, abs->paf[2]))
    printf ("o");
  else
    printf ("n");
  t1p_aff_check_free (pr, abs->paf[2]);
  t1p_aff_check_free (pr, abs->paf[0]);
  /* x = bottom
     y = sqrt(x)*sqrt(x) 
     y = bottom
   */
  abs->paf[0] = pr->bot;
  abs->paf[0]->pby++;
  abs->paf[2] = t1p_aff_eval_ap_texpr0 (pr, texp, abs);
  //t1p_aff_fprint(pr, stdout, abs->paf[2]);
  if (t1p_aff_is_bottom (pr, abs->paf[2]))
    printf ("o");
  else
    printf ("n");

  /* rump polynome
     x = 77617
     y = 33096
     z = 1335y^6/4 + x*x(11x^2y^2 - y^6 - 121y^4 - 2) + 11y^8/2 + x/2y
     z = a + b + c + d
   */
  t1p_aff_check_free (pr, abs->paf[0]);
  t1p_aff_check_free (pr, abs->paf[1]);
  t1p_aff_check_free (pr, abs->paf[2]);

  abs->paf[0] = t1p_aff_alloc_init (pr);
  abs->paf[0]->pby++;
  abs->paf[1] = t1p_aff_alloc_init (pr);
  abs->paf[1]->pby++;

  itv_set_int (abs->paf[0]->c, 77617);
  itv_set_int (abs->paf[0]->itv, 77617);
  itv_set_int (abs->paf[1]->c, 33096);
  itv_set_int (abs->paf[1]->itv, 33096);

  ap_texpr0_t *a = ap_texpr0_binop (AP_TEXPR_DIV,
				    ap_texpr0_binop (AP_TEXPR_MUL,
						     ap_texpr0_cst_scalar_int (1335),
						     ap_texpr0_binop (AP_TEXPR_MUL,
								      ap_texpr0_binop
								      (AP_TEXPR_MUL,
								       ap_texpr0_binop
								       (AP_TEXPR_MUL,
									ap_texpr0_binop
									(AP_TEXPR_MUL,
									 ap_texpr0_binop
									 (AP_TEXPR_MUL,
									  ap_texpr0_dim
									  (1),
									  ap_texpr0_dim
									  (1),
									  AP_RTYPE_DOUBLE,
									  AP_RDIR_UP),
									 ap_texpr0_dim
									 (1),
									 AP_RTYPE_DOUBLE,
									 AP_RDIR_UP),
									ap_texpr0_dim
									(1),
									AP_RTYPE_DOUBLE,
									AP_RDIR_UP),
								       ap_texpr0_dim
								       (1),
								       AP_RTYPE_DOUBLE,
								       AP_RDIR_UP),
								      ap_texpr0_dim
								      (1),
								      AP_RTYPE_DOUBLE,
								      AP_RDIR_UP),
						     AP_RTYPE_DOUBLE, AP_RDIR_UP),
				    ap_texpr0_cst_scalar_int (4),
				    AP_RTYPE_DOUBLE, AP_RDIR_UP);

  ap_texpr0_t *bb = ap_texpr0_binop (AP_TEXPR_SUB,
				     ap_texpr0_binop (AP_TEXPR_SUB,
						      ap_texpr0_binop (AP_TEXPR_MUL,
								       ap_texpr0_cst_scalar_int
								       (11),
								       ap_texpr0_binop
								       (AP_TEXPR_MUL,
									ap_texpr0_dim
									(0),
									ap_texpr0_binop
									(AP_TEXPR_MUL,
									 ap_texpr0_dim
									 (0),
									 ap_texpr0_binop
									 (AP_TEXPR_MUL,
									  ap_texpr0_dim
									  (1),
									  ap_texpr0_dim
									  (1),
									  AP_RTYPE_DOUBLE,
									  AP_RDIR_UP),
									 AP_RTYPE_DOUBLE,
									 AP_RDIR_UP),
									AP_RTYPE_DOUBLE,
									AP_RDIR_UP),
								       AP_RTYPE_DOUBLE,
								       AP_RDIR_UP),
						      ap_texpr0_binop (AP_TEXPR_MUL,
								       ap_texpr0_dim
								       (1),
								       ap_texpr0_binop
								       (AP_TEXPR_MUL,
									ap_texpr0_dim
									(1),
									ap_texpr0_binop
									(AP_TEXPR_MUL,
									 ap_texpr0_dim
									 (1),
									 ap_texpr0_binop
									 (AP_TEXPR_MUL,
									  ap_texpr0_dim
									  (1),
									  ap_texpr0_binop
									  (AP_TEXPR_MUL,
									   ap_texpr0_dim
									   (1),
									   ap_texpr0_dim
									   (1),
									   AP_RTYPE_DOUBLE,
									   AP_RDIR_UP),
									  AP_RTYPE_DOUBLE,
									  AP_RDIR_UP),
									 AP_RTYPE_DOUBLE,
									 AP_RDIR_UP),
									AP_RTYPE_DOUBLE,
									AP_RDIR_UP),
								       AP_RTYPE_DOUBLE,
								       AP_RDIR_UP),
						      AP_RTYPE_DOUBLE, AP_RDIR_UP),
				     ap_texpr0_binop (AP_TEXPR_ADD,
						      ap_texpr0_binop (AP_TEXPR_MUL,
								       ap_texpr0_cst_scalar_int
								       (121),
								       ap_texpr0_binop
								       (AP_TEXPR_MUL,
									ap_texpr0_dim
									(1),
									ap_texpr0_binop
									(AP_TEXPR_MUL,
									 ap_texpr0_dim
									 (1),
									 ap_texpr0_binop
									 (AP_TEXPR_MUL,
									  ap_texpr0_dim
									  (1),
									  ap_texpr0_dim
									  (1),
									  AP_RTYPE_DOUBLE,
									  AP_RDIR_UP),
									 AP_RTYPE_DOUBLE,
									 AP_RDIR_UP),
									AP_RTYPE_DOUBLE,
									AP_RDIR_UP),
								       AP_RTYPE_DOUBLE,
								       AP_RDIR_UP),
						      ap_texpr0_cst_scalar_int (2),
						      AP_RTYPE_DOUBLE, AP_RDIR_UP),
				     AP_RTYPE_DOUBLE, AP_RDIR_UP);

  ap_texpr0_t *b = ap_texpr0_binop (AP_TEXPR_MUL,
				    ap_texpr0_dim (0),
				    ap_texpr0_binop (AP_TEXPR_MUL,
						     ap_texpr0_dim (0),
						     bb,
						     AP_RTYPE_DOUBLE, AP_RDIR_UP),
				    AP_RTYPE_DOUBLE, AP_RDIR_UP);

  ap_texpr0_t *c = ap_texpr0_binop (AP_TEXPR_DIV,
				    ap_texpr0_binop (AP_TEXPR_MUL,
						     ap_texpr0_cst_scalar_int (11),
						     ap_texpr0_binop (AP_TEXPR_MUL,
								      ap_texpr0_binop
								      (AP_TEXPR_MUL,
								       ap_texpr0_binop
								       (AP_TEXPR_MUL,
									ap_texpr0_binop
									(AP_TEXPR_MUL,
									 ap_texpr0_binop
									 (AP_TEXPR_MUL,
									  ap_texpr0_binop
									  (AP_TEXPR_MUL,
									   ap_texpr0_binop
									   (AP_TEXPR_MUL,
									    ap_texpr0_dim
									    (1),
									    ap_texpr0_dim
									    (1),
									    AP_RTYPE_DOUBLE,
									    AP_RDIR_UP),
									   ap_texpr0_dim
									   (1),
									   AP_RTYPE_DOUBLE,
									   AP_RDIR_UP),
									  ap_texpr0_dim
									  (1),
									  AP_RTYPE_DOUBLE,
									  AP_RDIR_UP),
									 ap_texpr0_dim
									 (1),
									 AP_RTYPE_DOUBLE,
									 AP_RDIR_UP),
									ap_texpr0_dim
									(1),
									AP_RTYPE_DOUBLE,
									AP_RDIR_UP),
								       ap_texpr0_dim
								       (1),
								       AP_RTYPE_DOUBLE,
								       AP_RDIR_UP),
								      ap_texpr0_dim
								      (1),
								      AP_RTYPE_DOUBLE,
								      AP_RDIR_UP),
						     AP_RTYPE_DOUBLE, AP_RDIR_UP),
				    ap_texpr0_cst_scalar_int (2),
				    AP_RTYPE_DOUBLE, AP_RDIR_UP);

  ap_texpr0_t *d = ap_texpr0_binop (AP_TEXPR_DIV,
				    ap_texpr0_dim (0),
				    ap_texpr0_binop (AP_TEXPR_MUL,
						     ap_texpr0_cst_scalar_int (2),
						     ap_texpr0_dim (1),
						     AP_RTYPE_DOUBLE, AP_RDIR_UP),
				    AP_RTYPE_DOUBLE, AP_RDIR_UP);

  ap_texpr0_t *rump = ap_texpr0_binop (AP_TEXPR_ADD,
				       ap_texpr0_binop (AP_TEXPR_ADD,
							a, b,
							AP_RTYPE_DOUBLE, AP_RDIR_UP),
				       ap_texpr0_binop (AP_TEXPR_ADD,
							c, d,
							AP_RTYPE_DOUBLE, AP_RDIR_UP),
				       AP_RTYPE_DOUBLE, AP_RDIR_UP);

  abs->paf[2] = t1p_aff_eval_ap_texpr0 (pr, rump, abs);
  /* x = -2.5 + 0.5*eps1 + eps2;
     z = 2x - x - x;
   */
  t1p_aff_check_free (pr, abs->paf[1]);
  t1p_aff_check_free (pr, abs->paf[2]);

  abs->paf[1] = t1p_aff_alloc_init (pr);
  abs->paf[1]->pby++;
  num_set_double (num, -2.5);
  itv_set_num (coeff, num);
  itv_set (abs->paf[1]->c, coeff);
  num_set_double (num, 0.5);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, abs->paf[1], coeff, 0);
  num_set_double (num, 1);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, abs->paf[1], coeff, 1);
  itv_set_int2 (abs->paf[1]->itv, -4, -1);

  ap_texpr0_free (texp);
  texp = ap_texpr0_binop (AP_TEXPR_SUB,
			  ap_texpr0_binop (AP_TEXPR_SUB,
					   ap_texpr0_binop (AP_TEXPR_MUL,
							    ap_texpr0_cst_scalar_int
							    (2), ap_texpr0_dim (1),
							    AP_RTYPE_DOUBLE,
							    AP_RDIR_UP),
					   ap_texpr0_dim (1), AP_RTYPE_DOUBLE,
					   AP_RDIR_UP), ap_texpr0_dim (1),
			  AP_RTYPE_DOUBLE, AP_RDIR_UP);
  abs->paf[2] = t1p_aff_eval_ap_texpr0 (pr, texp, abs);
  if (t1p_aff_is_zero (pr, abs->paf[2]))
    printf ("o");
  else
    printf ("n");

  ap_texpr0_free (texp);
  texp = ap_texpr0_binop (AP_TEXPR_MUL,
			  ap_texpr0_dim (0),
			  ap_texpr0_binop (AP_TEXPR_SUB,
					   ap_texpr0_binop (AP_TEXPR_MUL,
							    ap_texpr0_cst_scalar_int
							    (2), ap_texpr0_dim (1),
							    AP_RTYPE_DOUBLE,
							    AP_RDIR_UP),
					   ap_texpr0_dim (1), AP_RTYPE_DOUBLE,
					   AP_RDIR_UP), AP_RTYPE_DOUBLE, AP_RDIR_UP);

  t1p_aff_check_free (pr, abs->paf[0]);
  abs->paf[0] = t1p_aff_alloc_init (pr);
  abs->paf[0]->pby++;
  num_set_double (num, -2.5);
  itv_set_num (coeff, num);
  itv_set (abs->paf[0]->c, coeff);
  num_set_double (num, 1.5);
  itv_set_num (coeff, num);
  t1p_aff_nsym_create (pr, abs->paf[0], coeff, UN);
  itv_set_int2 (abs->paf[0]->itv, -4, -1);

  t1p_aff_check_free (pr, abs->paf[2]);
  abs->paf[2] = t1p_aff_eval_ap_texpr0 (pr, texp, abs);
  t1p_aff_reduce (pr, abs->paf[2]);

  t1p_aff_t *aff_exact = t1p_aff_alloc_init (pr);
  num_set_double (num, 6.25);
  itv_set_num (coeff, num);
  itv_set (aff_exact->c, coeff);
  num_set_double (num, -1.25);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, aff_exact, coeff, 0);
  num_set_double (num, -2.5);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, aff_exact, coeff, 1);
  num_set_double (num, -3.75);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, aff_exact, coeff, (-2 + pr->dim));
  num_set_double (num, 2.25);
  itv_set_num (coeff, num);
  t1p_aff_build (pr, aff_exact, coeff, (-1 + pr->dim));
  itv_set_int2 (aff_exact->itv, 1, 16);
  if (t1p_aff_is_perfectly_eq (pr, abs->paf[2], aff_exact))
    printf ("o");
  else
    printf ("n");

  num_clear (num);
  itv_clear (coeff);
  itv_clear (lambda);
  ap_texpr0_free (texp);
  ap_texpr0_free (rump);
  //t1p_aff_free(pr, dummy);
  t1p_aff_free (pr, aff_exact);
  //t1p_aff_free(pr, aff_exact1);
  printf ("\n");
}

int
main (void)
{
  //  mtrace();
  //long int seed48 = time(NULL);
  //unsigned int seed = (unsigned int) (seed48 % UINT_MAX);
  //srand(seed);
  ap_manager_t *man = t1p_manager_alloc ();
  t1p_internal_t *pr = t1p_init_from_manager (man, AP_FUNID_UNKNOWN);
  t1p_t *abs = t1p_top (man, 0, 2);
  //t1p_t* abs1=t1p_copy(man,abs);

  test_display (pr, abs);
  //test_pm_non_constrained(pr, abs);
  /*test_mul_non_constrained(pr, abs);
     test_div_non_constrained(pr, abs);
     test_sqrt_non_constrained(pr, abs);
     test_eval_texpr_non_constrained(pr, abs);

     test_mul_constrained(pr, abs); */

  /* free managers */
  t1p_free (man, abs);
  ap_manager_free (man);
  return 0;
}
