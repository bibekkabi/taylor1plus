/*
   APRON Library / Taylor1+ Domain (beta version)
   Copyright (C) 2009-2011 Khalil Ghorbal

*/


#include <stdlib.h>

#include "num.h"
#include "itv.h"
#include "itv_linexpr.h"
#include "itv_linearize.h"

#include "ap_manager.h"
#include "ap_tcons0.h"
#include "ap_lincons0.h"
#include "ap_generator0.h"
#include "ap_generic.h"

#include "t1p_internal.h"
#include "t1p_representation.h"
#include "t1p_constructor.h"

#include "t1p_itv_utils.h"
#include "t1p_fun.h"

#include "t1p_meetjoin.h"
#include <glpk.h>


////////// old meet /////////////////////

t1p_t* t1p_meet(ap_manager_t* man, bool destructive, t1p_t* a1, t1p_t* a2)
/* TODO: destructive not used */
{
    CALL();
    t1p_internal_t* pr = t1p_init_from_manager(man,AP_FUNID_MEET);
    arg_assert(a1->dims==a2->dims && a1->intdim==a2->intdim,abort(););
#ifdef _T1P_DEBUG
    fprintf(stdout, "### MEET OPERANDS (destructive %d)###\n",destructive);
    t1p_fprint(stdout, man, a1, 0x0);
    t1p_fprint(stdout, man, a2, 0x0);
    fprintf(stdout, "### ### ###\n");
#endif
    t1p_t* res = NULL;
    size_t i = 0;
    size_t realdim = a1->dims - a1->intdim;
    size_t intdim = a1->intdim;
    man->result.flag_best = tbool_true;
    man->result.flag_exact = tbool_true;
    if (t1p_is_eq(man, a1, a2)) res = t1p_copy(man, a1);
    else if (tbool_or(t1p_is_bottom(man, a1), t1p_is_bottom(man, a2)) == tbool_true) {
        res = t1p_bottom(man, intdim, realdim);
    } else if (t1p_is_top(man, a1) == tbool_true) {
        res = t1p_copy(man, a2);
    } else if (t1p_is_top(man, a2) == tbool_true) {
        res = t1p_copy(man, a1);
    } else {
        res = t1p_alloc(man, intdim, realdim);
        //ap_abstract0_free(pr->manNS, res->abs);
        size_t k = 0;
        ap_dim_t j = 0;
        size_t dims1 = t1p_nsymcons_get_dimension(pr, a1);
        size_t dims2 = t1p_nsymcons_get_dimension(pr, a2);
        if (dims1 && dims2) {
            size_t dim2 = 0;
            /* au pire il y a toute la liste de l'un a rajouter dans l'autre */
            ap_dimchange_t* dimchange2 = ap_dimchange_alloc(0, dims2);
            if (dims1 > res->size) {
                res->nsymcons = (ap_dim_t*)realloc(res->nsymcons, (dims1)*sizeof(ap_dim_t));
                res->gamma = (ap_interval_t**)realloc(res->gamma, (dims1)*sizeof(ap_interval_t*));
                for (k=res->size;k<dims1;k++) res->gamma[k] = NULL;
                res->size = dims1;
            }
            res->nsymcons = memcpy((void *)res->nsymcons, (const void *)a1->nsymcons, dims1*sizeof(ap_dim_t));
            for (i=0; i<dims1; i++) res->gamma[i] = ap_interval_alloc_set(a1->gamma[i]);
            ap_abstract0_free(pr->manNS, res->abs);
            res->abs = ap_abstract0_copy(pr->manNS, a1->abs);
            for (k=0; k<dims1; k++) {
                if (!t1p_nsymcons_get_dimpos(pr, &j, a1->nsymcons[k], a2)) {
                    dimchange2->dim[dim2] = j;
                    dim2++;
                }
            }
            dimchange2->realdim = dim2;
            size_t dim1 = 0;
            for (k=0; k<dims2; k++) t1p_insert_constrained_nsym(pr, &j, a2->nsymcons[k], res);

            /* destructive, without projection (new dimension set to top) */
            ap_abstract0_add_dimensions(pr->manNS, true, a2->abs, dimchange2, false);

            ap_abstract0_meet(pr->manNS, true, res->abs, a2->abs);

            /* update res->gamma */
            t1p_update_nsymcons_gamma(pr, res);

            ap_dimchange_add_invert(dimchange2);
            ap_abstract0_remove_dimensions(pr->manNS, true, a2->abs, dimchange2);

            dimchange2->realdim = dims2; ap_dimchange_free(dimchange2);
        } else if (dims1) {
            if (dims1 > res->size) {
                res->nsymcons = (ap_dim_t*)realloc(res->nsymcons, (dims1)*sizeof(ap_dim_t));
                res->gamma = (ap_interval_t**)realloc(res->gamma, (dims1)*sizeof(ap_interval_t*));
                for (k=res->size;k<dims1;k++) res->gamma[k] = NULL;
                res->size = dims1;
            }
            res->nsymcons = memcpy((void *)res->nsymcons, (const void *)a1->nsymcons, dims1*sizeof(ap_dim_t));
            for (i=0; i<dims1; i++) res->gamma[i] = ap_interval_alloc_set(a1->gamma[i]);
            res->abs = ap_abstract0_copy(pr->manNS, a1->abs);
        } else if (dims2) {
            if (dims2 > res->size) {
                res->nsymcons = (ap_dim_t*)realloc(res->nsymcons, (dims2)*sizeof(ap_dim_t));
                res->gamma = (ap_interval_t**)realloc(res->gamma, (dims2)*sizeof(ap_interval_t*));
                for (k=res->size;k<dims2;k++) res->gamma[k] = NULL;
                res->size = dims2;
            }
            res->nsymcons = memcpy((void *)res->nsymcons, (const void *)a2->nsymcons, dims2*sizeof(ap_dim_t));
            for (i=0; i<dims2; i++) res->gamma[i] = ap_interval_alloc_set(a2->gamma[i]);
            res->abs = ap_abstract0_copy(pr->manNS, a2->abs);
        } else {
            /* non constrained abstract objects */
        }

        //res->abs = ap_abstract0_meet(pr->manNS, false, a1->abs, a2->abs);
        /* update res->gamma */
        //t1p_update_nsymcons_gamma(pr, res);
        itv_t tmp; itv_init(tmp);
        for (i=0; i<intdim+realdim; i++) {
            /* update res->box */
            itv_meet(pr->itv, res->box[i], a1->box[i], a2->box[i]);
            if ((a1->paf[i]->q == NULL) && (a2->paf[i]->q == NULL)) {
                itv_meet(pr->itv, tmp, a1->paf[i]->c, a2->paf[i]->c);
                if (itv_is_bottom(pr->itv, tmp)) {
                    t1p_free(man, res);
                    res = t1p_bottom(man, intdim, realdim);
                    break;
                } else if (itv_is_top(tmp)) res->paf[i] = pr->top;
                else if (itv_has_finite_bound(tmp)) {
                    res->paf[i] = t1p_aff_alloc_init(pr);
                    t1p_aff_add_itv(pr, res->paf[i], tmp, IN);
                } else {
                    res->paf[i] = t1p_aff_alloc_init(pr);
                    itv_set(res->paf[i]->c, tmp);
                }
                res->paf[i]->pby++;
            } else if (t1p_aff_is_eq(pr, a1->paf[i], a2->paf[i])) {
                res->paf[i] = a1->paf[i];
                res->paf[i]->pby++;
            } else if (itv_is_point(pr->itv, res->box[i])) {
                res->paf[i] = t1p_aff_alloc_init(pr);
                itv_set(res->paf[i]->c, res->box[i]);
                // The following line was missing which caused a double free bug. Detected and fixed by Alexandra-Olimpia Bugariu on May 31st 2018 
                res->paf[i]->pby++;
            } else {
                fprintf(stdout, "destructive ? %d\n",destructive);
                t1p_fprint(stdout, man, a1, 0x0);
                t1p_fprint(stdout, man, a2, 0x0);
                //not_implemented();
                /* return a top instead */
                res->paf[i] = pr->top;
                res->paf[i]->pby++;
                itv_set_top(res->box[i]);
            }
        }
        itv_clear(tmp);
    }
#ifdef _T1P_DEBUG
    fprintf(stdout, "### RESULT of MEET ###\n");
    t1p_fprint(stdout, man, res, 0x0);
    fprintf(stdout, "### ### ###\n");
#endif
    return res;
}


/* Meet and Join */
/*****************/

/************************************************/
/* 1.Meet (bounding box meet)*/
/************************************************/

#if 0
t1p_t *
t1p_meet (ap_manager_t * man, bool destructive, t1p_t * a1, t1p_t * a2)
/* TODO: destructive not used */
{
 CALL();
 int i;
 t1p_internal_t* prbx1 = t1p_init_from_manager(man,AP_FUNID_MEET);
 #ifdef _T1P_DEBUG
    fprintf(stdout, "### MEET OPERANDS (destructive %d)###\n",destructive);
    t1p_fprint(stdout, man, a1, 0x0);
    t1p_fprint(stdout, man, a2, 0x0);
    fprintf(stdout, "### ### ###\n");
 #endif
 ap_interval_t** int1=t1p_to_box(man, a1);
 ap_interval_t** int2=t1p_to_box(man, a2);

 ap_manager_t* man_box = box_manager_alloc();

 box_t* box1=box_of_box(man_box, 0, 2, int1);
 box_t* box2=box_of_box(man_box, 0, 2, int2);

 box_t* bxmeet=box_meet(man_box, true , box1, box2);
 ap_interval_t** intmeet=box_to_box(man_box, bxmeet);
 t1p_t* tpmeet=t1p_of_box(man, 0, 2, intmeet);
 /* update tpmeet->box */  
 for (i=0; i<tpmeet->dims; i++){ 
      itv_meet(prbx1->itv, tpmeet->box[i], a1->box[i], a2->box[i]);
      itv_set(tpmeet->paf[i]->itv, tpmeet->box[i]);
      tpmeet->paf[i]->pby++;}
 #ifdef _T1P_DEBUG
    fprintf(stdout, "### RESULT of MEET ###\n");
    t1p_fprint(stdout, man, tpmeet, 0x0);
    fprintf(stdout, "### ### ###\n");
 #endif
 return tpmeet;
}
#endif

#if 0

/************************************************/
/* 1.Meet as in weighted minkowski sum (without GLPK)*/   //// corner (motivational) example/////////
/************************************************/

t1p_t *
t1p_meet (ap_manager_t * man, bool destructive, t1p_t * a1, t1p_t * a2)
/* TODO: destructive not used */
{
  //int b;
  //scanf("%f",&b);
  CALL ();
  t1p_internal_t *prbx1 = t1p_init_from_manager (man, AP_FUNID_MEET);
  #ifdef _T1P_DEBUG
    fprintf(stdout, "### MEET OPERANDS (destructive %d)###\n",destructive);
    t1p_fprint(stdout, man, a1, 0x0);
    t1p_fprint(stdout, man, a2, 0x0);
    fprintf(stdout, "### ### ###\n");
  #endif
  scanf("%f",&b);
  itv_t coeff;
  itv_init (coeff);
  num_t num;
  num_init (num);
  num_set_double (num, 0);
  itv_set_num (coeff, num);
  num_t num2;
  num_init (num2);
  num_set_double (num2, 1.05);       //// 0.1925 ///////
  itv_t num2coeff;
  itv_init (num2coeff);
  itv_set_num (num2coeff,num2);
  num_t num4;
  num_init (num4);
  num_set_double (num4, 1.15);
  itv_t num4coeff;
  itv_init (num4coeff);
  itv_set_num (num4coeff,num4);
  num_t num6;
  num_init (num6);
  num_set_double (num6, 0.15);    ////////// -0.2122 //////////
  itv_t num6coeff;
  itv_init (num6coeff);
  itv_set_num (num6coeff,num6);
  num_t num7;
  num_init (num7);
  num_set_double (num7, 0.025);
  itv_t num7coeff;
  itv_init (num7coeff);
  itv_set_num (num7coeff,num7);
  num_t num8;
  num_init (num8);
  num_set_double (num8, 1.075);
  itv_t num8coeff;
  itv_init (num8coeff);
  itv_set_num (num8coeff,num8);
  num_t num9;
  num_init (num9);
  num_set_double (num9, 0.15);
  itv_t num9coeff;
  itv_init (num9coeff);
  itv_set_num (num9coeff,num9);
  

  int i, j, size; //, size1, size2;

///////// first generator matrix ///////////////////

  /*itv_t **ptitv1 = (itv_t **) calloc (a1->dims, sizeof (itv_t *));
  for (i = 0; i < a1->dims; i++)
    ptitv1[i] = (itv_t *) calloc (prbx1->dim, sizeof (itv_t));

  size1 = t1p_aff_get_generator_matrix_size (prbx1, a1, ptitv1);

//printf("%u", size1);

  itv_t **ptitv2 = (itv_t **) calloc (a1->dims, sizeof (itv_t *));
  for (i = 0; i < a1->dims; i++)
    ptitv2[i] = (itv_t *) calloc (size1, sizeof (itv_t));

  t1p_aff_get_generator_matrix_num (prbx1, a1, ptitv1, ptitv2, num2);

//////////////// second generator matrix //////////////////

  itv_t **ptitv3 = (itv_t **) calloc (a2->dims, sizeof (itv_t *));
  for (i = 0; i < a2->dims; i++)
    ptitv3[i] = (itv_t *) calloc (prbx1->dim, sizeof (itv_t));

  size2 = t1p_aff_get_generator_matrix_size (prbx1, a2, ptitv3);

//printf("%u", size2);

  itv_t **ptitv4 = (itv_t **) calloc (a2->dims, sizeof (itv_t *));
  for (i = 0; i < a2->dims; i++)
    ptitv4[i] = (itv_t *) calloc (size2, sizeof (itv_t));

  t1p_aff_get_generator_matrix_num (prbx1, a2, ptitv3, ptitv4, num4);*/


  size = 4; //size1 + size2;

  itv_t **ptitv_meet = (itv_t **) calloc (a2->dims, sizeof (itv_t *));
  for (i = 0; i < a2->dims; i++)
    ptitv_meet[i] = (itv_t *) calloc (size, sizeof (itv_t));

/////////// fixed size buffer /////////////////

  //itv_t **ptitv_meet[a1->dims][size];

//printf("%u", size1);

//printf("%u", size2);

  //scanf("%f",&b);


 /* for (i = 0; i < a1->dims; i++) {
    for (j = 0; j < size; j++) {
      if (j >= (2)) {
	itv_set (ptitv_meet[i][j], ptitv4[i][j - size1]);
	//itv_print(ptitv_meet[i][j]);
      }
      else {
	itv_set (ptitv_meet[i][j], ptitv2[i][j]);
	//itv_print(ptitv_meet[i][j]);
      }
      //itv_print(ptitv[i][j]);
    }
  }*/

itv_set (ptitv_meet[0][0], num2coeff);
itv_set (ptitv_meet[0][1], coeff);
itv_set (ptitv_meet[0][2], num4coeff);
itv_set (ptitv_meet[0][3], num6coeff);
//itv_set (ptitv_meet[0][4], num7coeff);
itv_set (ptitv_meet[1][0], num9coeff);
itv_set (ptitv_meet[1][1], num2coeff);
itv_set (ptitv_meet[1][2], num7coeff);
itv_set (ptitv_meet[1][3], num8coeff);
//itv_set (ptitv_meet[1][4], coeff);

  //scanf("%f",&b);

//itv_print(ptitv3[0][2]);

//int k=2*prbx1->dim;
//prbx1->dim=0;

//printf("%u",prbx1->dim);

  t1p_t *tpmeet = t1p_alloc (man, 0, 2);
  tpmeet->paf[0] = t1p_aff_alloc_init (prbx1);
  tpmeet->paf[1] = t1p_aff_alloc_init (prbx1);

  t1p_nsym_t *nsymeet1;		//[size];
  //for (i=0;i<tb->dims;i++){
  t1p_aaterm_t *p1 = NULL;
  p1 = a1->paf[0]->q;
  /*t1p_aaterm_t *p2 = NULL;
     p2=a2->paf[0]->q; */
  t1p_aaterm_t *p3 = NULL;
  p3 = a1->paf[1]->q;

for (j = 0; j < size; j++) {
      nsymeet1 = t1p_nsym_add (prbx1, UN);  //// IN is for eps and UN is for eta
      t1p_aff_nsym_add (prbx1, tpmeet->paf[0], ptitv_meet[0][j], nsymeet1);	//////// creating new noise symbols ///////////
      //scanf("%f",&b);
      //printf("%u", prbx1->dim);
      //scanf("%f",&b);
      t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], prbx1->dim - 1);	//// using the index of the new noise symbols created above /////////////////
}

  /*for (j = 0; j < size; j++) {
    if (j >= size1) {
      //if (p2 != NULL){
      //t1p_aff_build (prbx1, tpmeet->paf[0], ptitv_meet[0][j], p2->pnsym->index);
      nsymeet1 = t1p_nsym_add (prbx1, UN);  //// IN is for eps and UN is for eta
      t1p_aff_nsym_add (prbx1, tpmeet->paf[0], ptitv_meet[0][j], nsymeet1);	//////// creating new noise symbols ///////////
      //scanf("%f",&b);
      //printf("%u", prbx1->dim);
      //scanf("%f",&b);
      t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], prbx1->dim - 1);	//// using the index of the new noise symbols created above /////////////////
      //printf("%u", p2->pnsym->index);
      //p2=p2->n;}
    }
    else {
      if (p1 != NULL && itv_is_eq (ptitv_meet[0][j], coeff) == 0) {
	t1p_aff_build (prbx1, tpmeet->paf[0], ptitv_meet[0][j], p1->pnsym->index);
	//printf("%u", p1->pnsym->index);
	p1 = p1->n;
      }
      if (p3 != NULL && itv_is_eq (ptitv_meet[1][j], coeff) == 0) {
	t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], p3->pnsym->index);
	//printf("%u", p3->pnsym->index);
	p3 = p3->n;
      }
    }

    //nsymeet1 = t1p_nsym_add(prbx1, IN);
    //t1p_aff_nsym_add(prbx1, tpmeet->paf[0], ptitv_meet[0][j], nsymeet1);    //////// creating new noise symbols ///////////
    //scanf("%f",&b);
    //printf("%u", prbx1->dim);
    //scanf("%f",&b);
    //t1p_aff_build(prbx1, tpmeet->paf[1], ptitv_meet[1][j], prbx1->dim-1);  //// using the index of the new noise symbols created above /////////////////

    //t1p_aff_build (prbx1, tpmeet->paf[0], ptitv_meet[0][j], j);  //////////// creating new noise symbols but indexed from j= to size ///////////////// 
  }*/

  /*t1p_aaterm_t *p3 = NULL;
     p3=a1->paf[1]->q;
     /*t1p_aaterm_t *p4 = NULL;
     p4=a2->paf[1]->q; */

  //for (j = 0; j < size; j++) {
  //t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], j);
  /*if (j >= size1){
     if (p4 !=NULL){
     t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], p4->pnsym->index);
     printf("%u", p4->pnsym->index);
     p4=p4->n;}}
     else {
     if (p3 !=NULL){
     t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], p3->pnsym->index);
     printf("%u", p3->pnsym->index);
     p3=p3->n;}} */
  //}

  //printf ("%u", prbx1->dim);

  /*itv_t a1_centre1_paf0;
  itv_init (a1_centre1_paf0);
  itv_t a1_centre2_paf1;
  itv_init (a1_centre2_paf1);
  itv_t a2_centre1_paf0;
  itv_init (a2_centre1_paf0);
  itv_t a2_centre2_paf1;
  itv_init (a2_centre2_paf1);*/

  /*itv_mul_num (a1_centre1_paf0, a1->paf[0]->c, num2);
  itv_mul_num (a1_centre2_paf1, a1->paf[1]->c, num2);
  itv_mul_num (a2_centre1_paf0, a2->paf[0]->c, num4);
  itv_mul_num (a2_centre2_paf1, a2->paf[1]->c, num4);
  itv_add (tpmeet->paf[0]->c, a1_centre1_paf0, a2_centre1_paf0);
  itv_add (tpmeet->paf[1]->c, a1_centre2_paf1, a2_centre2_paf1);*/

  
  itv_set_num (tpmeet->paf[0]->c,num);
  itv_set_num (tpmeet->paf[1]->c,num);

  

  /* update tpmeet->box */  
  for (i=0; i<tpmeet->dims; i++){ 
      //itv_meet(prbx1->itv, tpmeet->box[i], a1->box[i], a2->box[i]);
      t1p_aff_boxize (prbx1, tpmeet->box[i], tpmeet->paf[i], tpmeet);
      itv_set(tpmeet->paf[i]->itv, tpmeet->box[i]);
      tpmeet->paf[i]->pby++;}

  //scanf("%f",&b);
  //itv_print(tpmeet->paf[1]->end->coeff);
  //printf("%u", tpmeet->paf[1]->end->pnsym->index);
  //scanf("%f",&b);
//itv_print(tpmeet->paf[0]->q->n->coeff);
//itv_print(tpmeet->paf[0]->q->n->n->coeff);
//itv_print(tpmeet->paf[0]->end->coeff);
//printf("%u", tpmeet->paf[0]->end->pnsym->index);
//itv_print(tpmeet->paf[1]->q->coeff);
//printf("%u", tpmeet->paf[1]->q->pnsym->index);
//itv_print(tpmeet->paf[1]->q->n->coeff);
//itv_print(tpmeet->paf[1]->q->n->n->coeff);
//itv_print(tpmeet->paf[1]->end->coeff);
//printf("%u", tpmeet->paf[1]->end->pnsym->index);

  //itv_print (a1->paf[0]->q->coeff);

  /*free (ptitv1);
  free (ptitv2);
  free (ptitv3);
  free (ptitv4);*/
  num_clear (num);
  num_clear (num2);
  num_clear (num4);
  num_clear (num6);
  num_clear (num7);
  num_clear (num8);
  num_clear (num9);
  itv_clear (coeff);
  itv_clear (num2coeff);
  itv_clear (num4coeff);
  itv_clear (num6coeff);
  itv_clear (num7coeff);
  itv_clear (num8coeff);
  itv_clear (num9coeff);
  for (i = 0; i < a2->dims; i++) {
    for (j = 0; j < size; j++) {
      itv_clear (ptitv_meet[i][j]);}}
  for (i = 0; i < a2->dims; i++) {
     free(ptitv_meet[i]);}
  free (ptitv_meet);
  
  #ifdef _T1P_DEBUG
    fprintf(stdout, "### RESULT of MEET ###\n");
    t1p_fprint(stdout, man, tpmeet, 0x0);
    fprintf(stdout, "### ### ###\n");
  #endif
  scanf("%f",&b);

  return tpmeet;

}

#endif

/************************************************/
/* 1.Meet as in weighted minkowski sum (without GLPK)*/   //// filter /////////
/************************************************/
#if 0

t1p_t *
t1p_meet (ap_manager_t * man, bool destructive, t1p_t * a1, t1p_t * a2)
/* TODO: destructive not used */
{
  int b;
  //scanf("%f",&b);
  CALL ();
  t1p_internal_t *prbx1 = t1p_init_from_manager (man, AP_FUNID_MEET);
  #ifdef _T1P_DEBUG
    fprintf(stdout, "### MEET OPERANDS (destructive %d)###\n",destructive);
    t1p_fprint(stdout, man, a1, 0x0);
    t1p_fprint(stdout, man, a2, 0x0);
    fprintf(stdout, "### ### ###\n");
  #endif
  scanf("%f",&b);
  itv_t coeff;
  itv_init (coeff);
  num_t num;
  num_init (num);
  num_set_double (num, 0);
  itv_set_num (coeff, num);
  num_t num2;
  num_init (num2);
  num_set_double (num2, 0.1925);       //// 0.1925 ///////
  itv_t num2coeff;
  itv_init (num2coeff);
  itv_set_num (num2coeff,num2);
  num_t num4;
  num_init (num4);
  num_set_double (num4, 3.595);
  itv_t num4coeff;
  itv_init (num4coeff);
  itv_set_num (num4coeff,num4);
  num_t num6;
  num_init (num6);
  num_set_double (num6, -0.2125);    ////////// -0.2122 //////////
  itv_t num6coeff;
  itv_init (num6coeff);
  itv_set_num (num6coeff,num6);
  /*num_t num7;
  num_init (num7);
  num_set_double (num7, -0.0003);
  itv_t num7coeff;
  itv_init (num7coeff);
  itv_set_num (num7coeff,num7);*/
  num_t num8;
  num_init (num8);
  num_set_double (num8, 1.6033);
  itv_t num8coeff;
  itv_init (num8coeff);
  itv_set_num (num8coeff,num8);
  num_t num9;
  num_init (num9);
  num_set_double (num9, 2.3967);
  itv_t num9coeff;
  itv_init (num9coeff);
  itv_set_num (num9coeff,num9);
  

  int i, j, size, size1, size2;

///////// first generator matrix ///////////////////

  /*itv_t **ptitv1 = (itv_t **) calloc (a1->dims, sizeof (itv_t *));
  for (i = 0; i < a1->dims; i++)
    ptitv1[i] = (itv_t *) calloc (prbx1->dim, sizeof (itv_t));

  size1 = t1p_aff_get_generator_matrix_size (prbx1, a1, ptitv1);

//printf("%u", size1);

  itv_t **ptitv2 = (itv_t **) calloc (a1->dims, sizeof (itv_t *));
  for (i = 0; i < a1->dims; i++)
    ptitv2[i] = (itv_t *) calloc (size1, sizeof (itv_t));

  t1p_aff_get_generator_matrix_num (prbx1, a1, ptitv1, ptitv2, num2);

//////////////// second generator matrix //////////////////

  itv_t **ptitv3 = (itv_t **) calloc (a2->dims, sizeof (itv_t *));
  for (i = 0; i < a2->dims; i++)
    ptitv3[i] = (itv_t *) calloc (prbx1->dim, sizeof (itv_t));

  size2 = t1p_aff_get_generator_matrix_size (prbx1, a2, ptitv3);

//printf("%u", size2);

  itv_t **ptitv4 = (itv_t **) calloc (a2->dims, sizeof (itv_t *));
  for (i = 0; i < a2->dims; i++)
    ptitv4[i] = (itv_t *) calloc (size2, sizeof (itv_t));

  t1p_aff_get_generator_matrix_num (prbx1, a2, ptitv3, ptitv4, num4);*/


  size = 4; //size1 + size2;

  itv_t **ptitv_meet = (itv_t **) calloc (a2->dims, sizeof (itv_t *));
  for (i = 0; i < a2->dims; i++)
    ptitv_meet[i] = (itv_t *) calloc (size, sizeof (itv_t));

/////////// fixed size buffer /////////////////

  //itv_t **ptitv_meet[a1->dims][size];

//printf("%u", size1);

//printf("%u", size2);

  //scanf("%f",&b);


 /* for (i = 0; i < a1->dims; i++) {
    for (j = 0; j < size; j++) {
      if (j >= (2)) {
	itv_set (ptitv_meet[i][j], ptitv4[i][j - size1]);
	//itv_print(ptitv_meet[i][j]);
      }
      else {
	itv_set (ptitv_meet[i][j], ptitv2[i][j]);
	//itv_print(ptitv_meet[i][j]);
      }
      //itv_print(ptitv[i][j]);
    }
  }*/

itv_set (ptitv_meet[0][0], num2coeff);
itv_set (ptitv_meet[0][1], coeff);
itv_set (ptitv_meet[0][2], num4coeff);
itv_set (ptitv_meet[0][3], num6coeff);
//itv_set (ptitv_meet[0][4], num7coeff);
itv_set (ptitv_meet[1][0], coeff);
itv_set (ptitv_meet[1][1], num8coeff);
itv_set (ptitv_meet[1][2], num9coeff);
itv_set (ptitv_meet[1][3], coeff);
//itv_set (ptitv_meet[1][4], coeff);

  //scanf("%f",&b);

//itv_print(ptitv3[0][2]);

//int k=2*prbx1->dim;
//prbx1->dim=0;

//printf("%u",prbx1->dim);

  t1p_t *tpmeet = t1p_alloc (man, 0, 2);
  tpmeet->paf[0] = t1p_aff_alloc_init (prbx1);
  tpmeet->paf[1] = t1p_aff_alloc_init (prbx1);

  t1p_nsym_t *nsymeet1;		//[size];
  //for (i=0;i<tb->dims;i++){
  t1p_aaterm_t *p1 = NULL;
  p1 = a1->paf[0]->q;
  /*t1p_aaterm_t *p2 = NULL;
     p2=a2->paf[0]->q; */
  t1p_aaterm_t *p3 = NULL;
  p3 = a1->paf[1]->q;

  for (j = 0; j < size; j++) {
      nsymeet1 = t1p_nsym_add (prbx1, UN);  //// IN is for eps and UN is for eta
      t1p_aff_nsym_add (prbx1, tpmeet->paf[0], ptitv_meet[0][j], nsymeet1);	//////// creating new noise symbols ///////////
      //scanf("%f",&b);
      //printf("%u", prbx1->dim);
      //scanf("%f",&b);
      t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], prbx1->dim - 1);	//// using the index of the new noise symbols created above /////////////////
}

  /*for (j = 0; j < size; j++) {
    if (j >= size1) {
      //if (p2 != NULL){
      //t1p_aff_build (prbx1, tpmeet->paf[0], ptitv_meet[0][j], p2->pnsym->index);
      nsymeet1 = t1p_nsym_add (prbx1, UN);  //// IN is for eps and UN is for eta
      t1p_aff_nsym_add (prbx1, tpmeet->paf[0], ptitv_meet[0][j], nsymeet1);	//////// creating new noise symbols ///////////
      //scanf("%f",&b);
      //printf("%u", prbx1->dim);
      //scanf("%f",&b);
      t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], prbx1->dim - 1);	//// using the index of the new noise symbols created above /////////////////
      //printf("%u", p2->pnsym->index);
      //p2=p2->n;}
    }
    else {
      if (p1 != NULL && itv_is_eq (ptitv_meet[0][j], coeff) == 0) {
	t1p_aff_build (prbx1, tpmeet->paf[0], ptitv_meet[0][j], p1->pnsym->index);
	//printf("%u", p1->pnsym->index);
	p1 = p1->n;
      }
      if (p3 != NULL && itv_is_eq (ptitv_meet[1][j], coeff) == 0) {
	t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], p3->pnsym->index);
	//printf("%u", p3->pnsym->index);
	p3 = p3->n;
      }
    }

    //nsymeet1 = t1p_nsym_add(prbx1, IN);
    //t1p_aff_nsym_add(prbx1, tpmeet->paf[0], ptitv_meet[0][j], nsymeet1);    //////// creating new noise symbols ///////////
    //scanf("%f",&b);
    //printf("%u", prbx1->dim);
    //scanf("%f",&b);
    //t1p_aff_build(prbx1, tpmeet->paf[1], ptitv_meet[1][j], prbx1->dim-1);  //// using the index of the new noise symbols created above /////////////////

    //t1p_aff_build (prbx1, tpmeet->paf[0], ptitv_meet[0][j], j);  //////////// creating new noise symbols but indexed from j= to size ///////////////// 
  }*/

  /*t1p_aaterm_t *p3 = NULL;
     p3=a1->paf[1]->q;
     /*t1p_aaterm_t *p4 = NULL;
     p4=a2->paf[1]->q; */

  //for (j = 0; j < size; j++) {
  //t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], j);
  /*if (j >= size1){
     if (p4 !=NULL){
     t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], p4->pnsym->index);
     printf("%u", p4->pnsym->index);
     p4=p4->n;}}
     else {
     if (p3 !=NULL){
     t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], p3->pnsym->index);
     printf("%u", p3->pnsym->index);
     p3=p3->n;}} */
  //}

  //printf ("%u", prbx1->dim);

  /*itv_t a1_centre1_paf0;
  itv_init (a1_centre1_paf0);
  itv_t a1_centre2_paf1;
  itv_init (a1_centre2_paf1);
  itv_t a2_centre1_paf0;
  itv_init (a2_centre1_paf0);
  itv_t a2_centre2_paf1;
  itv_init (a2_centre2_paf1);*/

  /*itv_mul_num (a1_centre1_paf0, a1->paf[0]->c, num2);
  itv_mul_num (a1_centre2_paf1, a1->paf[1]->c, num2);
  itv_mul_num (a2_centre1_paf0, a2->paf[0]->c, num4);
  itv_mul_num (a2_centre2_paf1, a2->paf[1]->c, num4);
  itv_add (tpmeet->paf[0]->c, a1_centre1_paf0, a2_centre1_paf0);
  itv_add (tpmeet->paf[1]->c, a1_centre2_paf1, a2_centre2_paf1);*/

  
  itv_set_num (tpmeet->paf[0]->c,num);
  itv_set_num (tpmeet->paf[1]->c,num);

  

  /* update tpmeet->box */  
  for (i=0; i<tpmeet->dims; i++){ 
      //itv_meet(prbx1->itv, tpmeet->box[i], a1->box[i], a2->box[i]);
      t1p_aff_boxize (prbx1, tpmeet->box[i], tpmeet->paf[i], tpmeet);
      itv_set(tpmeet->paf[i]->itv, tpmeet->box[i]);
      tpmeet->paf[i]->pby++;}

  //scanf("%f",&b);
  //itv_print(tpmeet->paf[1]->end->coeff);
  //printf("%u", tpmeet->paf[1]->end->pnsym->index);
  //scanf("%f",&b);
//itv_print(tpmeet->paf[0]->q->n->coeff);
//itv_print(tpmeet->paf[0]->q->n->n->coeff);
//itv_print(tpmeet->paf[0]->end->coeff);
//printf("%u", tpmeet->paf[0]->end->pnsym->index);
//itv_print(tpmeet->paf[1]->q->coeff);
//printf("%u", tpmeet->paf[1]->q->pnsym->index);
//itv_print(tpmeet->paf[1]->q->n->coeff);
//itv_print(tpmeet->paf[1]->q->n->n->coeff);
//itv_print(tpmeet->paf[1]->end->coeff);
//printf("%u", tpmeet->paf[1]->end->pnsym->index);

  //itv_print (a1->paf[0]->q->coeff);

  /*free (ptitv1);
  free (ptitv2);
  free (ptitv3);
  free (ptitv4);*/
  num_clear (num);
  num_clear (num2);
  num_clear (num4);
  num_clear (num6);
  //num_clear (num7);
  num_clear (num8);
  num_clear (num9);
  itv_clear (coeff);
  itv_clear (num2coeff);
  itv_clear (num4coeff);
  itv_clear (num6coeff);
  //itv_clear (num7coeff);
  itv_clear (num8coeff);
  itv_clear (num9coeff);
  for (i = 0; i < a2->dims; i++) {
    for (j = 0; j < size; j++) {
      itv_clear (ptitv_meet[i][j]);}}
  for (i = 0; i < a2->dims; i++) {
     free(ptitv_meet[i]);}
  free (ptitv_meet);
  
 #ifdef _T1P_DEBUG
    fprintf(stdout, "### RESULT of MEET ###\n");
    t1p_fprint(stdout, man, tpmeet, 0x0);
    fprintf(stdout, "### ### ###\n");
  #endif
  scanf("%f",&b);

  return tpmeet;

}
#endif


/************************************************/
/* 1.Meet as in weighted minkowski sum (without GLPK)*/  ///// sine /////////
/************************************************/

#if 0

t1p_t *
t1p_meet (ap_manager_t * man, bool destructive, t1p_t * a1, t1p_t * a2)
/* TODO: destructive not used */
{
  int b;
  //scanf("%f",&b);
  CALL ();
  t1p_internal_t *prbx1 = t1p_init_from_manager (man, AP_FUNID_MEET);
  #ifdef _T1P_DEBUG
    fprintf(stdout, "### MEET OPERANDS (destructive %d)###\n",destructive);
    t1p_fprint(stdout, man, a1, 0x0);
    t1p_fprint(stdout, man, a2, 0x0);
    fprintf(stdout, "### ### ###\n");
  #endif
  scanf("%f",&b);

 //////////////////// sine [-1.05,1.05] bounds ///////////////
  itv_t coeff;
  itv_init (coeff);
  num_t num;
  num_init (num);
  num_set_double (num, 0);
  itv_set_num (coeff, num);
  num_t num2;
  num_init (num2);
  num_set_double (num2, 0.4844);  //// 0.4824 for bounds [1.1,1.1] and 0.4844 for bounds [1.05,1.05] 
  itv_t num2coeff;
  itv_init (num2coeff);
  itv_set_num (num2coeff,num2);
  num_t num4;
  num_init (num4);
  num_set_double (num4, 0.5413);
  itv_t num4coeff;
  itv_init (num4coeff);
  itv_set_num (num4coeff,num4);
  num_t num6;
  num_init (num6);
  num_set_double (num6, 0.0573); ///// 0.0763 for bounds [1.1,1.1] and (0.0573 or 0.0569) for bounds [1.05,1.05] 
  itv_t num6coeff;
  itv_init (num6coeff);
  itv_set_num (num6coeff,num6);
  num_t num7;
  num_init (num7);
  num_set_double (num7, 1.2555);
  itv_t num7coeff;
  itv_init (num7coeff);
  itv_set_num (num7coeff,num7);
  num_t num8;
  num_init (num8);
  num_set_double (num8, 0.7445);
  itv_t num8coeff;
  itv_init (num8coeff);
  itv_set_num (num8coeff,num8);
  /*num_t num9;
  num_init (num9);
  num_set_double (num9, 2.3967);
  itv_t num9coeff;
  itv_init (num9coeff);
  itv_set_num (num9coeff,num9);


  //////////////////// sine [-1.00921,1.00921] bounds ///////////////
  /*itv_t coeff;
  itv_init (coeff);
  num_t num;
  num_init (num);
  num_set_double (num, 0);
  itv_set_num (coeff, num);
  num_t num2;
  num_init (num2);
  num_set_double (num2, 0.5144);  //// 0.4824 for bounds [1.1,1.1] and 0.4844 for bounds [1.05,1.05] 
  itv_t num2coeff;
  itv_init (num2coeff);
  itv_set_num (num2coeff,num2);
  num_t num4;
  num_init (num4);
  num_set_double (num4, 0.5113);
  itv_t num4coeff;
  itv_init (num4coeff);
  itv_set_num (num4coeff,num4);
  num_t num6;
  num_init (num6);
  num_set_double (num6, 0.0002); ///// 0.0763 for bounds [1.1,1.1] and (0.0573 or 0.0569) for bounds [1.05,1.05] 
  itv_t num6coeff;
  itv_init (num6coeff);
  itv_set_num (num6coeff,num6);
  num_t num7;
  num_init (num7);
  num_set_double (num7, 1.2555);
  itv_t num7coeff;
  itv_init (num7coeff);
  itv_set_num (num7coeff,num7);
  num_t num8;
  num_init (num8);
  num_set_double (num8, 0.7445);
  itv_t num8coeff;
  itv_init (num8coeff);
  itv_set_num (num8coeff,num8);*/
  

  int i, j, size, size1, size2;

///////// first generator matrix ///////////////////

  /*itv_t **ptitv1 = (itv_t **) calloc (a1->dims, sizeof (itv_t *));
  for (i = 0; i < a1->dims; i++)
    ptitv1[i] = (itv_t *) calloc (prbx1->dim, sizeof (itv_t));

  size1 = t1p_aff_get_generator_matrix_size (prbx1, a1, ptitv1);

//printf("%u", size1);

  itv_t **ptitv2 = (itv_t **) calloc (a1->dims, sizeof (itv_t *));
  for (i = 0; i < a1->dims; i++)
    ptitv2[i] = (itv_t *) calloc (size1, sizeof (itv_t));

  t1p_aff_get_generator_matrix_num (prbx1, a1, ptitv1, ptitv2, num2);

//////////////// second generator matrix //////////////////

  itv_t **ptitv3 = (itv_t **) calloc (a2->dims, sizeof (itv_t *));
  for (i = 0; i < a2->dims; i++)
    ptitv3[i] = (itv_t *) calloc (prbx1->dim, sizeof (itv_t));

  size2 = t1p_aff_get_generator_matrix_size (prbx1, a2, ptitv3);

//printf("%u", size2);

  itv_t **ptitv4 = (itv_t **) calloc (a2->dims, sizeof (itv_t *));
  for (i = 0; i < a2->dims; i++)
    ptitv4[i] = (itv_t *) calloc (size2, sizeof (itv_t));

  t1p_aff_get_generator_matrix_num (prbx1, a2, ptitv3, ptitv4, num4);*/


  size = 4; //size1 + size2;

  itv_t **ptitv_meet = (itv_t **) calloc (a2->dims, sizeof (itv_t *));
  for (i = 0; i < a2->dims; i++)
    ptitv_meet[i] = (itv_t *) calloc (size, sizeof (itv_t));

/////////// fixed size buffer /////////////////

  //itv_t **ptitv_meet[a1->dims][size];

//printf("%u", size1);

//printf("%u", size2);

  //scanf("%f",&b);


 /* for (i = 0; i < a1->dims; i++) {
    for (j = 0; j < size; j++) {
      if (j >= (2)) {
	itv_set (ptitv_meet[i][j], ptitv4[i][j - size1]);
	//itv_print(ptitv_meet[i][j]);
      }
      else {
	itv_set (ptitv_meet[i][j], ptitv2[i][j]);
	//itv_print(ptitv_meet[i][j]);
      }
      //itv_print(ptitv[i][j]);
    }
  }*/

itv_set (ptitv_meet[0][0], num6coeff);
itv_set (ptitv_meet[0][1], num2coeff);
itv_set (ptitv_meet[0][2], coeff);
itv_set (ptitv_meet[0][3], num4coeff);
itv_set (ptitv_meet[1][0], coeff);
itv_set (ptitv_meet[1][1], coeff);
itv_set (ptitv_meet[1][2], num7coeff);
itv_set (ptitv_meet[1][3], num8coeff);

  //scanf("%f",&b);

//itv_print(ptitv3[0][2]);

//int k=2*prbx1->dim;
//prbx1->dim=0;

//printf("%u",prbx1->dim);

  t1p_t *tpmeet = t1p_alloc (man, 0, 2);
  tpmeet->paf[0] = t1p_aff_alloc_init (prbx1);
  tpmeet->paf[1] = t1p_aff_alloc_init (prbx1);

  t1p_nsym_t *nsymeet1;		//[size];
  //for (i=0;i<tb->dims;i++){
  t1p_aaterm_t *p1 = NULL;
  p1 = a1->paf[0]->q;
  /*t1p_aaterm_t *p2 = NULL;
     p2=a2->paf[0]->q; */
  t1p_aaterm_t *p3 = NULL;
  p3 = a1->paf[1]->q;

  for (j = 0; j < size; j++) {
    if (j >= size1) {
      //if (p2 != NULL){
      //t1p_aff_build (prbx1, tpmeet->paf[0], ptitv_meet[0][j], p2->pnsym->index);
      nsymeet1 = t1p_nsym_add (prbx1, UN);  //// IN is for eps and UN is for eta
      t1p_aff_nsym_add (prbx1, tpmeet->paf[0], ptitv_meet[0][j], nsymeet1);	//////// creating new noise symbols ///////////
      //scanf("%f",&b);
      //printf("%u", prbx1->dim);
      //scanf("%f",&b);
      t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], prbx1->dim - 1);	//// using the index of the new noise symbols created above /////////////////
      //printf("%u", p2->pnsym->index);
      //p2=p2->n;}
    }
    else {
      if (p1 != NULL && itv_is_eq (ptitv_meet[0][j], coeff) == 0) {
	t1p_aff_build (prbx1, tpmeet->paf[0], ptitv_meet[0][j], p1->pnsym->index);
	//printf("%u", p1->pnsym->index);
	p1 = p1->n;
      }
      if (p3 != NULL && itv_is_eq (ptitv_meet[1][j], coeff) == 0) {
	t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], p3->pnsym->index);
	//printf("%u", p3->pnsym->index);
	p3 = p3->n;
      }
    }

    //nsymeet1 = t1p_nsym_add(prbx1, IN);
    //t1p_aff_nsym_add(prbx1, tpmeet->paf[0], ptitv_meet[0][j], nsymeet1);    //////// creating new noise symbols ///////////
    //scanf("%f",&b);
    //printf("%u", prbx1->dim);
    //scanf("%f",&b);
    //t1p_aff_build(prbx1, tpmeet->paf[1], ptitv_meet[1][j], prbx1->dim-1);  //// using the index of the new noise symbols created above /////////////////

    //t1p_aff_build (prbx1, tpmeet->paf[0], ptitv_meet[0][j], j);  //////////// creating new noise symbols but indexed from j= to size ///////////////// 
  }

  /*t1p_aaterm_t *p3 = NULL;
     p3=a1->paf[1]->q;
     /*t1p_aaterm_t *p4 = NULL;
     p4=a2->paf[1]->q; */

  //for (j = 0; j < size; j++) {
  //t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], j);
  /*if (j >= size1){
     if (p4 !=NULL){
     t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], p4->pnsym->index);
     printf("%u", p4->pnsym->index);
     p4=p4->n;}}
     else {
     if (p3 !=NULL){
     t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], p3->pnsym->index);
     printf("%u", p3->pnsym->index);
     p3=p3->n;}} */
  //}

  //printf ("%u", prbx1->dim);

  /*itv_t a1_centre1_paf0;
  itv_init (a1_centre1_paf0);
  itv_t a1_centre2_paf1;
  itv_init (a1_centre2_paf1);
  itv_t a2_centre1_paf0;
  itv_init (a2_centre1_paf0);
  itv_t a2_centre2_paf1;
  itv_init (a2_centre2_paf1);

  itv_mul_num (a1_centre1_paf0, a1->paf[0]->c, num);
  itv_mul_num (a1_centre2_paf1, a1->paf[1]->c, num);
  itv_mul_num (a2_centre1_paf0, a2->paf[0]->c, num);
  itv_mul_num (a2_centre2_paf1, a2->paf[1]->c, num);
  itv_add (tpmeet->paf[0]->c, a1_centre1_paf0, a2_centre1_paf0);
  itv_add (tpmeet->paf[1]->c, a1_centre2_paf1, a2_centre2_paf1);*/


  itv_set_num (tpmeet->paf[0]->c,num);
  itv_set_num (tpmeet->paf[1]->c,num);

  /* update tpmeet->box */  
  for (i=0; i<tpmeet->dims; i++){ 
      //itv_meet(prbx1->itv, tpmeet->box[i], a1->box[i], a2->box[i]);
      t1p_aff_boxize (prbx1, tpmeet->box[i], tpmeet->paf[i], tpmeet);
      itv_set(tpmeet->paf[i]->itv, tpmeet->box[i]);
      tpmeet->paf[i]->pby++;}

  //scanf("%f",&b);
  //itv_print(tpmeet->paf[1]->end->coeff);
  //printf("%u", tpmeet->paf[1]->end->pnsym->index);
  //scanf("%f",&b);
//itv_print(tpmeet->paf[0]->q->n->coeff);
//itv_print(tpmeet->paf[0]->q->n->n->coeff);
//itv_print(tpmeet->paf[0]->end->coeff);
//printf("%u", tpmeet->paf[0]->end->pnsym->index);
//itv_print(tpmeet->paf[1]->q->coeff);
//printf("%u", tpmeet->paf[1]->q->pnsym->index);
//itv_print(tpmeet->paf[1]->q->n->coeff);
//itv_print(tpmeet->paf[1]->q->n->n->coeff);
//itv_print(tpmeet->paf[1]->end->coeff);
//printf("%u", tpmeet->paf[1]->end->pnsym->index);

  //itv_print (a1->paf[0]->q->coeff);

  /*free (ptitv1);
  free (ptitv2);
  free (ptitv3);
  free (ptitv4);*/
  num_clear (num);
  num_clear (num2);
  num_clear (num4);
  num_clear (num6);
  num_clear (num7);
  num_clear (num8);
  //num_clear (num9);
  itv_clear (coeff);
  itv_clear (num2coeff);
  itv_clear (num4coeff);
  itv_clear (num6coeff);
  itv_clear (num7coeff);
  itv_clear (num8coeff);
  //itv_clear (num9coeff);
  for (i = 0; i < a2->dims; i++) {
    for (j = 0; j < size; j++) {
      itv_clear (ptitv_meet[i][j]);}}
  for (i = 0; i < a2->dims; i++) {
     free(ptitv_meet[i]);}
  free (ptitv_meet);
  
  #ifdef _T1P_DEBUG
    fprintf(stdout, "### RESULT of MEET ###\n");
    t1p_fprint(stdout, man, tpmeet, 0x0);
    fprintf(stdout, "### ### ###\n");
  #endif

  return tpmeet;

}
#endif

/************************************************/
/* 1.Meet as in weighted minkowski sum (without GLPK)*/  ///// newton /////////
/************************************************/

#if 0

t1p_t *
t1p_meet (ap_manager_t * man, bool destructive, t1p_t * a1, t1p_t * a2)
/* TODO: destructive not used */
{
  int b;
  //scanf("%f",&b);
  CALL ();
  t1p_internal_t *prbx1 = t1p_init_from_manager (man, AP_FUNID_MEET);
  #ifdef _T1P_DEBUG
    fprintf(stdout, "### MEET OPERANDS (destructive %d)###\n",destructive);
    t1p_fprint(stdout, man, a1, 0x0);
    t1p_fprint(stdout, man, a2, 0x0);
    fprintf(stdout, "### ### ###\n");
  #endif
  scanf("%f",&b);

  /////////////// 12 iterations, 12 reachable //////////////

  /*itv_t coeff;
  itv_init (coeff);
  num_t num;
  num_init (num);
  num_set_double (num, 0);
  itv_set_num (coeff, num);
  num_t num2;
  num_init (num2);
  num_set_double (num2, 0.6947);  //// 0.4824 ///////
  itv_t num2coeff;
  itv_init (num2coeff);
  itv_set_num (num2coeff,num2);
  num_t num4;
  num_init (num4);
  num_set_double (num4, -0.3060);
  itv_t num4coeff;
  itv_init (num4coeff);
  itv_set_num (num4coeff,num4);
  num_t num6;
  num_init (num6);
  num_set_double (num6, -0.0641);
  itv_t num6coeff;
  itv_init (num6coeff);
  itv_set_num (num6coeff,num6);
  num_t num7;
  num_init (num7);
  num_set_double (num7, -0.0105);
  itv_t num7coeff;
  itv_init (num7coeff);
  itv_set_num (num7coeff,num7);
  num_t num8;
  num_init (num8);
  num_set_double (num8, 0.2315);
  itv_t num8coeff;
  itv_init (num8coeff);
  itv_set_num (num8coeff,num8);*/
  /*num_t num9;
  num_init (num9);
  num_set_double (num9, 2.3967);
  itv_t num9coeff;
  itv_init (num9coeff);
  itv_set_num (num9coeff,num9);*/


  //////////////////// 11 iterations, 9 reachable //////////

 
  itv_t coeff;
  itv_init (coeff);
  num_t num;
  num_init (num);
  num_set_double (num, 0);
  itv_set_num (coeff, num);
  num_t num2;
  num_init (num2);
  num_set_double (num2, 0.7256);  //// 0.6956 for [0.6,0.6] and 0.7256 for [0.56,0.56]
  itv_t num2coeff;
  itv_init (num2coeff);
  itv_set_num (num2coeff,num2);
  num_t num4;
  num_init (num4);
  num_set_double (num4, -0.2869); /// -0.3069 for [0.6,0.6] and -0.2869 for [0.56,0.56]
  itv_t num4coeff;
  itv_init (num4coeff);
  itv_set_num (num4coeff,num4);
  num_t num6;
  num_init (num6);
  num_set_double (num6, -0.0642);
  itv_t num6coeff;
  itv_init (num6coeff);
  itv_set_num (num6coeff,num6);
  num_t num7;
  num_init (num7);
  num_set_double (num7, -0.0105);
  itv_t num7coeff;
  itv_init (num7coeff);
  itv_set_num (num7coeff,num7);
  num_t num8;
  num_init (num8);
  num_set_double (num8, 0.2); /// 0.2322 for [0.6,0.6] and 0.2 for [0.56,0.56]
  itv_t num8coeff;
  itv_init (num8coeff);
  itv_set_num (num8coeff,num8);
  /*num_t num9;
  num_init (num9);
  num_set_double (num9, 2.3967);
  itv_t num9coeff;
  itv_init (num9coeff);
  itv_set_num (num9coeff,num9);*/
  

  int i, j, size, size1, size2;

///////// first generator matrix ///////////////////

  /*itv_t **ptitv1 = (itv_t **) calloc (a1->dims, sizeof (itv_t *));
  for (i = 0; i < a1->dims; i++)
    ptitv1[i] = (itv_t *) calloc (prbx1->dim, sizeof (itv_t));

  size1 = t1p_aff_get_generator_matrix_size (prbx1, a1, ptitv1);

//printf("%u", size1);

  itv_t **ptitv2 = (itv_t **) calloc (a1->dims, sizeof (itv_t *));
  for (i = 0; i < a1->dims; i++)
    ptitv2[i] = (itv_t *) calloc (size1, sizeof (itv_t));

  t1p_aff_get_generator_matrix_num (prbx1, a1, ptitv1, ptitv2, num2);

//////////////// second generator matrix //////////////////

  itv_t **ptitv3 = (itv_t **) calloc (a2->dims, sizeof (itv_t *));
  for (i = 0; i < a2->dims; i++)
    ptitv3[i] = (itv_t *) calloc (prbx1->dim, sizeof (itv_t));

  size2 = t1p_aff_get_generator_matrix_size (prbx1, a2, ptitv3);

//printf("%u", size2);

  itv_t **ptitv4 = (itv_t **) calloc (a2->dims, sizeof (itv_t *));
  for (i = 0; i < a2->dims; i++)
    ptitv4[i] = (itv_t *) calloc (size2, sizeof (itv_t));

  t1p_aff_get_generator_matrix_num (prbx1, a2, ptitv3, ptitv4, num4);*/


  size = 5; //size1 + size2;

  itv_t **ptitv_meet = (itv_t **) calloc (a2->dims, sizeof (itv_t *));
  for (i = 0; i < a2->dims; i++)
    ptitv_meet[i] = (itv_t *) calloc (size, sizeof (itv_t));

/////////// fixed size buffer /////////////////

  //itv_t **ptitv_meet[a1->dims][size];

//printf("%u", size1);

//printf("%u", size2);

  //scanf("%f",&b);


 /* for (i = 0; i < a1->dims; i++) {
    for (j = 0; j < size; j++) {
      if (j >= (2)) {
	itv_set (ptitv_meet[i][j], ptitv4[i][j - size1]);
	//itv_print(ptitv_meet[i][j]);
      }
      else {
	itv_set (ptitv_meet[i][j], ptitv2[i][j]);
	//itv_print(ptitv_meet[i][j]);
      }
      //itv_print(ptitv[i][j]);
    }
  }*/

itv_set (ptitv_meet[0][0], coeff);
itv_set (ptitv_meet[0][1], num4coeff);
itv_set (ptitv_meet[0][2], num6coeff);
itv_set (ptitv_meet[0][3], num7coeff);
itv_set (ptitv_meet[0][4], num8coeff);
itv_set (ptitv_meet[1][0], num2coeff);
itv_set (ptitv_meet[1][1], coeff);
itv_set (ptitv_meet[1][2], num6coeff);
itv_set (ptitv_meet[1][3], num7coeff);
itv_set (ptitv_meet[1][4], num8coeff);

  //scanf("%f",&b);

//itv_print(ptitv3[0][2]);

//int k=2*prbx1->dim;
//prbx1->dim=0;

//printf("%u",prbx1->dim);

  t1p_t *tpmeet = t1p_alloc (man, 0, 2);
  tpmeet->paf[0] = t1p_aff_alloc_init (prbx1);
  tpmeet->paf[1] = t1p_aff_alloc_init (prbx1);

  t1p_nsym_t *nsymeet1;		//[size];
  //for (i=0;i<tb->dims;i++){
  t1p_aaterm_t *p1 = NULL;
  p1 = a1->paf[0]->q;
  /*t1p_aaterm_t *p2 = NULL;
     p2=a2->paf[0]->q; */
  t1p_aaterm_t *p3 = NULL;
  p3 = a1->paf[1]->q;

  for (j = 0; j < size; j++) {
    if (j >= size1) {
      //if (p2 != NULL){
      //t1p_aff_build (prbx1, tpmeet->paf[0], ptitv_meet[0][j], p2->pnsym->index);
      nsymeet1 = t1p_nsym_add (prbx1, UN);  //// IN is for eps and UN is for eta
      t1p_aff_nsym_add (prbx1, tpmeet->paf[0], ptitv_meet[0][j], nsymeet1);	//////// creating new noise symbols ///////////
      //scanf("%f",&b);
      //printf("%u", prbx1->dim);
      //scanf("%f",&b);
      t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], prbx1->dim - 1);	//// using the index of the new noise symbols created above /////////////////
      //printf("%u", p2->pnsym->index);
      //p2=p2->n;}
    }
    else {
      if (p1 != NULL && itv_is_eq (ptitv_meet[0][j], coeff) == 0) {
	t1p_aff_build (prbx1, tpmeet->paf[0], ptitv_meet[0][j], p1->pnsym->index);
	//printf("%u", p1->pnsym->index);
	p1 = p1->n;
      }
      if (p3 != NULL && itv_is_eq (ptitv_meet[1][j], coeff) == 0) {
	t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], p3->pnsym->index);
	//printf("%u", p3->pnsym->index);
	p3 = p3->n;
      }
    }

    //nsymeet1 = t1p_nsym_add(prbx1, IN);
    //t1p_aff_nsym_add(prbx1, tpmeet->paf[0], ptitv_meet[0][j], nsymeet1);    //////// creating new noise symbols ///////////
    //scanf("%f",&b);
    //printf("%u", prbx1->dim);
    //scanf("%f",&b);
    //t1p_aff_build(prbx1, tpmeet->paf[1], ptitv_meet[1][j], prbx1->dim-1);  //// using the index of the new noise symbols created above /////////////////

    //t1p_aff_build (prbx1, tpmeet->paf[0], ptitv_meet[0][j], j);  //////////// creating new noise symbols but indexed from j= to size ///////////////// 
  }

  /*t1p_aaterm_t *p3 = NULL;
     p3=a1->paf[1]->q;
     /*t1p_aaterm_t *p4 = NULL;
     p4=a2->paf[1]->q; */

  //for (j = 0; j < size; j++) {
  //t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], j);
  /*if (j >= size1){
     if (p4 !=NULL){
     t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], p4->pnsym->index);
     printf("%u", p4->pnsym->index);
     p4=p4->n;}}
     else {
     if (p3 !=NULL){
     t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], p3->pnsym->index);
     printf("%u", p3->pnsym->index);
     p3=p3->n;}} */
  //}

  //printf ("%u", prbx1->dim);

  /*itv_t a1_centre1_paf0;
  itv_init (a1_centre1_paf0);
  itv_t a1_centre2_paf1;
  itv_init (a1_centre2_paf1);
  itv_t a2_centre1_paf0;
  itv_init (a2_centre1_paf0);
  itv_t a2_centre2_paf1;
  itv_init (a2_centre2_paf1);

  itv_mul_num (a1_centre1_paf0, a1->paf[0]->c, num);
  itv_mul_num (a1_centre2_paf1, a1->paf[1]->c, num);
  itv_mul_num (a2_centre1_paf0, a2->paf[0]->c, num);
  itv_mul_num (a2_centre2_paf1, a2->paf[1]->c, num);
  itv_add (tpmeet->paf[0]->c, a1_centre1_paf0, a2_centre1_paf0);
  itv_add (tpmeet->paf[1]->c, a1_centre2_paf1, a2_centre2_paf1);*/

   itv_set_num (tpmeet->paf[0]->c,num);
  itv_set_num (tpmeet->paf[1]->c,num);

  /* update tpmeet->box */  
  for (i=0; i<tpmeet->dims; i++){ 
      //itv_meet(prbx1->itv, tpmeet->box[i], a1->box[i], a2->box[i]);
      t1p_aff_boxize (prbx1, tpmeet->box[i], tpmeet->paf[i], tpmeet);
      itv_set(tpmeet->paf[i]->itv, tpmeet->box[i]);
      tpmeet->paf[i]->pby++;}

  //scanf("%f",&b);
  //itv_print(tpmeet->paf[1]->end->coeff);
  //printf("%u", tpmeet->paf[1]->end->pnsym->index);
  //scanf("%f",&b);
//itv_print(tpmeet->paf[0]->q->n->coeff);
//itv_print(tpmeet->paf[0]->q->n->n->coeff);
//itv_print(tpmeet->paf[0]->end->coeff);
//printf("%u", tpmeet->paf[0]->end->pnsym->index);
//itv_print(tpmeet->paf[1]->q->coeff);
//printf("%u", tpmeet->paf[1]->q->pnsym->index);
//itv_print(tpmeet->paf[1]->q->n->coeff);
//itv_print(tpmeet->paf[1]->q->n->n->coeff);
//itv_print(tpmeet->paf[1]->end->coeff);
//printf("%u", tpmeet->paf[1]->end->pnsym->index);

  //itv_print (a1->paf[0]->q->coeff);

  /*free (ptitv1);
  free (ptitv2);
  free (ptitv3);
  free (ptitv4);*/
  num_clear (num);
  num_clear (num2);
  num_clear (num4);
  num_clear (num6);
  num_clear (num7);
  num_clear (num8);
  //num_clear (num9);
  itv_clear (coeff);
  itv_clear (num2coeff);
  itv_clear (num4coeff);
  itv_clear (num6coeff);
  itv_clear (num7coeff);
  itv_clear (num8coeff);
  //itv_clear (num9coeff);
  for (i = 0; i < a2->dims; i++) {
    for (j = 0; j < size; j++) {
      itv_clear (ptitv_meet[i][j]);}}
  for (i = 0; i < a2->dims; i++) {
     free(ptitv_meet[i]);}
  free (ptitv_meet);
  
  #ifdef _T1P_DEBUG
    fprintf(stdout, "### RESULT of MEET ###\n");
    t1p_fprint(stdout, man, tpmeet, 0x0);
    fprintf(stdout, "### ### ###\n");
  #endif

    scanf("%f",&b);

  return tpmeet;

}
#endif

/************************************************/
/* 1.Meet as in weighted minkowski sum (without GLPK)*/  ///// newton2 /////////
/************************************************/

#if 0

t1p_t *
t1p_meet (ap_manager_t * man, bool destructive, t1p_t * a1, t1p_t * a2)
/* TODO: destructive not used */
{
  int b;
  //scanf("%f",&b);
  CALL ();
  t1p_internal_t *prbx1 = t1p_init_from_manager (man, AP_FUNID_MEET);
  #ifdef _T1P_DEBUG
    fprintf(stdout, "### MEET OPERANDS (destructive %d)###\n",destructive);
    t1p_fprint(stdout, man, a1, 0x0);
    t1p_fprint(stdout, man, a2, 0x0);
    fprintf(stdout, "### ### ###\n");
  #endif
  scanf("%f",&b);

  /////////////// 12 iterations, 12 reachable //////////////

  /*itv_t coeff;
  itv_init (coeff);
  num_t num;
  num_init (num);
  num_set_double (num, 0);
  itv_set_num (coeff, num);
  num_t num2;
  num_init (num2);
  num_set_double (num2, 0.6947);  //// 0.4824 ///////
  itv_t num2coeff;
  itv_init (num2coeff);
  itv_set_num (num2coeff,num2);
  num_t num4;
  num_init (num4);
  num_set_double (num4, -0.3060);
  itv_t num4coeff;
  itv_init (num4coeff);
  itv_set_num (num4coeff,num4);
  num_t num6;
  num_init (num6);
  num_set_double (num6, -0.0641);
  itv_t num6coeff;
  itv_init (num6coeff);
  itv_set_num (num6coeff,num6);
  num_t num7;
  num_init (num7);
  num_set_double (num7, -0.0105);
  itv_t num7coeff;
  itv_init (num7coeff);
  itv_set_num (num7coeff,num7);
  num_t num8;
  num_init (num8);
  num_set_double (num8, 0.2315);
  itv_t num8coeff;
  itv_init (num8coeff);
  itv_set_num (num8coeff,num8);*/
  /*num_t num9;
  num_init (num9);
  num_set_double (num9, 2.3967);
  itv_t num9coeff;
  itv_init (num9coeff);
  itv_set_num (num9coeff,num9);*/


  //////////////////// 11 iterations, 9 reachable //////////

 
  itv_t coeff;
  itv_init (coeff);
  num_t num;
  num_init (num);
  num_set_double (num, 0);
  itv_set_num (coeff, num);
  num_t num2;
  num_init (num2);
  num_set_double (num2, 0.9256);  //// 0.6956 for [0.6,0.6] and 0.7256 for [0.56,0.56]
  itv_t num2coeff;
  itv_init (num2coeff);
  itv_set_num (num2coeff,num2);
  num_t num4;
  num_init (num4);
  num_set_double (num4, -0.0869); /// -0.3069 for [0.6,0.6] and -0.2869 for [0.56,0.56]
  itv_t num4coeff;
  itv_init (num4coeff);
  itv_set_num (num4coeff,num4);
  num_t num6;
  num_init (num6);
  num_set_double (num6, -0.0542);
  itv_t num6coeff;
  itv_init (num6coeff);
  itv_set_num (num6coeff,num6);
  num_t num7;
  num_init (num7);
  num_set_double (num7, -0.01);
  itv_t num7coeff;
  itv_init (num7coeff);
  itv_set_num (num7coeff,num7);
  num_t num8;
  num_init (num8);
  num_set_double (num8, 0.02); /// 0.2322 for [0.6,0.6] and 0.2 for [0.56,0.56]
  itv_t num8coeff;
  itv_init (num8coeff);
  itv_set_num (num8coeff,num8);
  /*num_t num9;
  num_init (num9);
  num_set_double (num9, 2.3967);
  itv_t num9coeff;
  itv_init (num9coeff);
  itv_set_num (num9coeff,num9);*/
  

  int i, j, size, size1, size2;

///////// first generator matrix ///////////////////

  /*itv_t **ptitv1 = (itv_t **) calloc (a1->dims, sizeof (itv_t *));
  for (i = 0; i < a1->dims; i++)
    ptitv1[i] = (itv_t *) calloc (prbx1->dim, sizeof (itv_t));

  size1 = t1p_aff_get_generator_matrix_size (prbx1, a1, ptitv1);

//printf("%u", size1);

  itv_t **ptitv2 = (itv_t **) calloc (a1->dims, sizeof (itv_t *));
  for (i = 0; i < a1->dims; i++)
    ptitv2[i] = (itv_t *) calloc (size1, sizeof (itv_t));

  t1p_aff_get_generator_matrix_num (prbx1, a1, ptitv1, ptitv2, num2);

//////////////// second generator matrix //////////////////

  itv_t **ptitv3 = (itv_t **) calloc (a2->dims, sizeof (itv_t *));
  for (i = 0; i < a2->dims; i++)
    ptitv3[i] = (itv_t *) calloc (prbx1->dim, sizeof (itv_t));

  size2 = t1p_aff_get_generator_matrix_size (prbx1, a2, ptitv3);

//printf("%u", size2);

  itv_t **ptitv4 = (itv_t **) calloc (a2->dims, sizeof (itv_t *));
  for (i = 0; i < a2->dims; i++)
    ptitv4[i] = (itv_t *) calloc (size2, sizeof (itv_t));

  t1p_aff_get_generator_matrix_num (prbx1, a2, ptitv3, ptitv4, num4);*/


  size = 5; //size1 + size2;

  itv_t **ptitv_meet = (itv_t **) calloc (a2->dims, sizeof (itv_t *));
  for (i = 0; i < a2->dims; i++)
    ptitv_meet[i] = (itv_t *) calloc (size, sizeof (itv_t));

/////////// fixed size buffer /////////////////

  //itv_t **ptitv_meet[a1->dims][size];

//printf("%u", size1);

//printf("%u", size2);

  //scanf("%f",&b);


 /* for (i = 0; i < a1->dims; i++) {
    for (j = 0; j < size; j++) {
      if (j >= (2)) {
	itv_set (ptitv_meet[i][j], ptitv4[i][j - size1]);
	//itv_print(ptitv_meet[i][j]);
      }
      else {
	itv_set (ptitv_meet[i][j], ptitv2[i][j]);
	//itv_print(ptitv_meet[i][j]);
      }
      //itv_print(ptitv[i][j]);
    }
  }*/

itv_set (ptitv_meet[0][0], coeff);
itv_set (ptitv_meet[0][1], num4coeff);
itv_set (ptitv_meet[0][2], num6coeff);
itv_set (ptitv_meet[0][3], num7coeff);
itv_set (ptitv_meet[0][4], num8coeff);
itv_set (ptitv_meet[1][0], num2coeff);
itv_set (ptitv_meet[1][1], coeff);
itv_set (ptitv_meet[1][2], num6coeff);
itv_set (ptitv_meet[1][3], num7coeff);
itv_set (ptitv_meet[1][4], num8coeff);

  //scanf("%f",&b);

//itv_print(ptitv3[0][2]);

//int k=2*prbx1->dim;
//prbx1->dim=0;

//printf("%u",prbx1->dim);

  t1p_t *tpmeet = t1p_alloc (man, 0, 2);
  tpmeet->paf[0] = t1p_aff_alloc_init (prbx1);
  tpmeet->paf[1] = t1p_aff_alloc_init (prbx1);

  t1p_nsym_t *nsymeet1;		//[size];
  //for (i=0;i<tb->dims;i++){
  t1p_aaterm_t *p1 = NULL;
  p1 = a1->paf[0]->q;
  /*t1p_aaterm_t *p2 = NULL;
     p2=a2->paf[0]->q; */
  t1p_aaterm_t *p3 = NULL;
  p3 = a1->paf[1]->q;

  for (j = 0; j < size; j++) {
    if (j >= size1) {
      //if (p2 != NULL){
      //t1p_aff_build (prbx1, tpmeet->paf[0], ptitv_meet[0][j], p2->pnsym->index);
      nsymeet1 = t1p_nsym_add (prbx1, UN);  //// IN is for eps and UN is for eta
      t1p_aff_nsym_add (prbx1, tpmeet->paf[0], ptitv_meet[0][j], nsymeet1);	//////// creating new noise symbols ///////////
      //scanf("%f",&b);
      //printf("%u", prbx1->dim);
      //scanf("%f",&b);
      t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], prbx1->dim - 1);	//// using the index of the new noise symbols created above /////////////////
      //printf("%u", p2->pnsym->index);
      //p2=p2->n;}
    }
    else {
      if (p1 != NULL && itv_is_eq (ptitv_meet[0][j], coeff) == 0) {
	t1p_aff_build (prbx1, tpmeet->paf[0], ptitv_meet[0][j], p1->pnsym->index);
	//printf("%u", p1->pnsym->index);
	p1 = p1->n;
      }
      if (p3 != NULL && itv_is_eq (ptitv_meet[1][j], coeff) == 0) {
	t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], p3->pnsym->index);
	//printf("%u", p3->pnsym->index);
	p3 = p3->n;
      }
    }

    //nsymeet1 = t1p_nsym_add(prbx1, IN);
    //t1p_aff_nsym_add(prbx1, tpmeet->paf[0], ptitv_meet[0][j], nsymeet1);    //////// creating new noise symbols ///////////
    //scanf("%f",&b);
    //printf("%u", prbx1->dim);
    //scanf("%f",&b);
    //t1p_aff_build(prbx1, tpmeet->paf[1], ptitv_meet[1][j], prbx1->dim-1);  //// using the index of the new noise symbols created above /////////////////

    //t1p_aff_build (prbx1, tpmeet->paf[0], ptitv_meet[0][j], j);  //////////// creating new noise symbols but indexed from j= to size ///////////////// 
  }

  /*t1p_aaterm_t *p3 = NULL;
     p3=a1->paf[1]->q;
     /*t1p_aaterm_t *p4 = NULL;
     p4=a2->paf[1]->q; */

  //for (j = 0; j < size; j++) {
  //t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], j);
  /*if (j >= size1){
     if (p4 !=NULL){
     t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], p4->pnsym->index);
     printf("%u", p4->pnsym->index);
     p4=p4->n;}}
     else {
     if (p3 !=NULL){
     t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], p3->pnsym->index);
     printf("%u", p3->pnsym->index);
     p3=p3->n;}} */
  //}

  //printf ("%u", prbx1->dim);

  /*itv_t a1_centre1_paf0;
  itv_init (a1_centre1_paf0);
  itv_t a1_centre2_paf1;
  itv_init (a1_centre2_paf1);
  itv_t a2_centre1_paf0;
  itv_init (a2_centre1_paf0);
  itv_t a2_centre2_paf1;
  itv_init (a2_centre2_paf1);

  itv_mul_num (a1_centre1_paf0, a1->paf[0]->c, num);
  itv_mul_num (a1_centre2_paf1, a1->paf[1]->c, num);
  itv_mul_num (a2_centre1_paf0, a2->paf[0]->c, num);
  itv_mul_num (a2_centre2_paf1, a2->paf[1]->c, num);
  itv_add (tpmeet->paf[0]->c, a1_centre1_paf0, a2_centre1_paf0);
  itv_add (tpmeet->paf[1]->c, a1_centre2_paf1, a2_centre2_paf1);*/

   itv_set_num (tpmeet->paf[0]->c,num);
  itv_set_num (tpmeet->paf[1]->c,num);

  /* update tpmeet->box */  
  for (i=0; i<tpmeet->dims; i++){ 
      //itv_meet(prbx1->itv, tpmeet->box[i], a1->box[i], a2->box[i]);
      t1p_aff_boxize (prbx1, tpmeet->box[i], tpmeet->paf[i], tpmeet);
      itv_set(tpmeet->paf[i]->itv, tpmeet->box[i]);
      tpmeet->paf[i]->pby++;}

  //scanf("%f",&b);
  //itv_print(tpmeet->paf[1]->end->coeff);
  //printf("%u", tpmeet->paf[1]->end->pnsym->index);
  //scanf("%f",&b);
//itv_print(tpmeet->paf[0]->q->n->coeff);
//itv_print(tpmeet->paf[0]->q->n->n->coeff);
//itv_print(tpmeet->paf[0]->end->coeff);
//printf("%u", tpmeet->paf[0]->end->pnsym->index);
//itv_print(tpmeet->paf[1]->q->coeff);
//printf("%u", tpmeet->paf[1]->q->pnsym->index);
//itv_print(tpmeet->paf[1]->q->n->coeff);
//itv_print(tpmeet->paf[1]->q->n->n->coeff);
//itv_print(tpmeet->paf[1]->end->coeff);
//printf("%u", tpmeet->paf[1]->end->pnsym->index);

  //itv_print (a1->paf[0]->q->coeff);

  /*free (ptitv1);
  free (ptitv2);
  free (ptitv3);
  free (ptitv4);*/
  num_clear (num);
  num_clear (num2);
  num_clear (num4);
  num_clear (num6);
  num_clear (num7);
  num_clear (num8);
  //num_clear (num9);
  itv_clear (coeff);
  itv_clear (num2coeff);
  itv_clear (num4coeff);
  itv_clear (num6coeff);
  itv_clear (num7coeff);
  itv_clear (num8coeff);
  //itv_clear (num9coeff);
  for (i = 0; i < a2->dims; i++) {
    for (j = 0; j < size; j++) {
      itv_clear (ptitv_meet[i][j]);}}
  for (i = 0; i < a2->dims; i++) {
     free(ptitv_meet[i]);}
  free (ptitv_meet);
  
  #ifdef _T1P_DEBUG
    fprintf(stdout, "### RESULT of MEET ###\n");
    t1p_fprint(stdout, man, tpmeet, 0x0);
    fprintf(stdout, "### ### ###\n");
  #endif

    scanf("%f",&b);

  return tpmeet;

}
#endif


#if 0

/************************************************/
/* 1.Meet as in weighted minkowski sum (without GLPK)*/ ////////// Harmonic oscillator ////////// Roux_FMSD15_Ex8.txt
/************************************************/

t1p_t *
t1p_meet (ap_manager_t * man, bool destructive, t1p_t * a1, t1p_t * a2)
/* TODO: destructive not used */
{
  //int b;
  //scanf("%f",&b);
  CALL ();
  t1p_internal_t *prbx1 = t1p_init_from_manager (man, AP_FUNID_MEET);
  #ifdef _T1P_DEBUG
    fprintf(stdout, "### MEET OPERANDS (destructive %d)###\n",destructive);
    t1p_fprint(stdout, man, a1, 0x0);
    t1p_fprint(stdout, man, a2, 0x0);
    fprintf(stdout, "### ### ###\n");
  #endif
  itv_t coeff;
  itv_init (coeff);
  num_t num;
  num_init (num);
  num_set_double (num, 0);
  itv_set_num (coeff, num);
  num_t num2;
  num_init (num2);
  num_set_double (num2, 0.3048);
  itv_t num2coeff;
  itv_init (num2coeff);
  itv_set_num (num2coeff,num2);
  num_t num4;
  num_init (num4);
  num_set_double (num4, 0.8734);
  itv_t num4coeff;
  itv_init (num4coeff);
  itv_set_num (num4coeff,num4);
  num_t num6;
  num_init (num6);
  num_set_double (num6, 0.0918);
  itv_t num6coeff;
  itv_init (num6coeff);
  itv_set_num (num6coeff,num6);
  num_t num7;
  num_init (num7);
  num_set_double (num7, 0.3040);
  itv_t num7coeff;
  itv_init (num7coeff);
  itv_set_num (num7coeff,num7);
  num_t num8;
  num_init (num8);
  num_set_double (num8, -0.0919);
  itv_t num8coeff;
  itv_init (num8coeff);
  itv_set_num (num8coeff,num8);
  num_t num9;
  num_init (num9);
  num_set_double (num9, 0.8740);
  itv_t num9coeff;
  itv_init (num9coeff);
  itv_set_num (num9coeff,num9);
  

  int i, j, size, size1, size2;

///////// first generator matrix ///////////////////

  /*itv_t **ptitv1 = (itv_t **) calloc (a1->dims, sizeof (itv_t *));
  for (i = 0; i < a1->dims; i++)
    ptitv1[i] = (itv_t *) calloc (prbx1->dim, sizeof (itv_t));

  size1 = t1p_aff_get_generator_matrix_size (prbx1, a1, ptitv1);

//printf("%u", size1);

  itv_t **ptitv2 = (itv_t **) calloc (a1->dims, sizeof (itv_t *));
  for (i = 0; i < a1->dims; i++)
    ptitv2[i] = (itv_t *) calloc (size1, sizeof (itv_t));

  t1p_aff_get_generator_matrix_num (prbx1, a1, ptitv1, ptitv2, num2);

//////////////// second generator matrix //////////////////

  itv_t **ptitv3 = (itv_t **) calloc (a2->dims, sizeof (itv_t *));
  for (i = 0; i < a2->dims; i++)
    ptitv3[i] = (itv_t *) calloc (prbx1->dim, sizeof (itv_t));

  size2 = t1p_aff_get_generator_matrix_size (prbx1, a2, ptitv3);

//printf("%u", size2);

  itv_t **ptitv4 = (itv_t **) calloc (a2->dims, sizeof (itv_t *));
  for (i = 0; i < a2->dims; i++)
    ptitv4[i] = (itv_t *) calloc (size2, sizeof (itv_t));

  t1p_aff_get_generator_matrix_num (prbx1, a2, ptitv3, ptitv4, num4);*/


  size = 4; //size1 + size2;

  itv_t **ptitv_meet = (itv_t **) calloc (a2->dims, sizeof (itv_t *));
  for (i = 0; i < a2->dims; i++)
    ptitv_meet[i] = (itv_t *) calloc (size, sizeof (itv_t));

/////////// fixed size buffer /////////////////

  //itv_t **ptitv_meet[a1->dims][size];

//printf("%u", size1);

//printf("%u", size2);

  //scanf("%f",&b);


 /* for (i = 0; i < a1->dims; i++) {
    for (j = 0; j < size; j++) {
      if (j >= (2)) {
	itv_set (ptitv_meet[i][j], ptitv4[i][j - size1]);
	//itv_print(ptitv_meet[i][j]);
      }
      else {
	itv_set (ptitv_meet[i][j], ptitv2[i][j]);
	//itv_print(ptitv_meet[i][j]);
      }
      //itv_print(ptitv[i][j]);
    }
  }*/

itv_set (ptitv_meet[0][0], num2coeff);
itv_set (ptitv_meet[0][1], coeff);
itv_set (ptitv_meet[0][2], num4coeff);
itv_set (ptitv_meet[0][3], num6coeff);
itv_set (ptitv_meet[1][0], coeff);
itv_set (ptitv_meet[1][1], num7coeff);
itv_set (ptitv_meet[1][2], num8coeff);
itv_set (ptitv_meet[1][3], num9coeff);
  //scanf("%f",&b);

//itv_print(ptitv3[0][2]);

//int k=2*prbx1->dim;
//prbx1->dim=0;

//printf("%u",prbx1->dim);

  t1p_t *tpmeet = t1p_alloc (man, 0, 2);
  tpmeet->paf[0] = t1p_aff_alloc_init (prbx1);
  tpmeet->paf[1] = t1p_aff_alloc_init (prbx1);

  t1p_nsym_t *nsymeet1;		//[size];
  //for (i=0;i<tb->dims;i++){
  t1p_aaterm_t *p1 = NULL;
  p1 = a1->paf[0]->q;
  /*t1p_aaterm_t *p2 = NULL;
     p2=a2->paf[0]->q; */
  t1p_aaterm_t *p3 = NULL;
  p3 = a1->paf[1]->q;

  for (j = 0; j < size; j++) {
    if (j >= size1) {
      //if (p2 != NULL){
      //t1p_aff_build (prbx1, tpmeet->paf[0], ptitv_meet[0][j], p2->pnsym->index);
      nsymeet1 = t1p_nsym_add (prbx1, UN);  //// IN is for eps and UN is for eta
      t1p_aff_nsym_add (prbx1, tpmeet->paf[0], ptitv_meet[0][j], nsymeet1);	//////// creating new noise symbols ///////////
      //scanf("%f",&b);
      //printf("%u", prbx1->dim);
      //scanf("%f",&b);
      t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], prbx1->dim - 1);	//// using the index of the new noise symbols created above /////////////////
      //printf("%u", p2->pnsym->index);
      //p2=p2->n;}
    }
    else {
      if (p1 != NULL && itv_is_eq (ptitv_meet[0][j], coeff) == 0) {
	t1p_aff_build (prbx1, tpmeet->paf[0], ptitv_meet[0][j], p1->pnsym->index);
	//printf("%u", p1->pnsym->index);
	p1 = p1->n;
      }
      if (p3 != NULL && itv_is_eq (ptitv_meet[1][j], coeff) == 0) {
	t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], p3->pnsym->index);
	//printf("%u", p3->pnsym->index);
	p3 = p3->n;
      }
    }

    //nsymeet1 = t1p_nsym_add(prbx1, IN);
    //t1p_aff_nsym_add(prbx1, tpmeet->paf[0], ptitv_meet[0][j], nsymeet1);    //////// creating new noise symbols ///////////
    //scanf("%f",&b);
    //printf("%u", prbx1->dim);
    //scanf("%f",&b);
    //t1p_aff_build(prbx1, tpmeet->paf[1], ptitv_meet[1][j], prbx1->dim-1);  //// using the index of the new noise symbols created above /////////////////

    //t1p_aff_build (prbx1, tpmeet->paf[0], ptitv_meet[0][j], j);  //////////// creating new noise symbols but indexed from j= to size ///////////////// 
  }

  /*t1p_aaterm_t *p3 = NULL;
     p3=a1->paf[1]->q;
     /*t1p_aaterm_t *p4 = NULL;
     p4=a2->paf[1]->q; */

  //for (j = 0; j < size; j++) {
  //t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], j);
  /*if (j >= size1){
     if (p4 !=NULL){
     t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], p4->pnsym->index);
     printf("%u", p4->pnsym->index);
     p4=p4->n;}}
     else {
     if (p3 !=NULL){
     t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], p3->pnsym->index);
     printf("%u", p3->pnsym->index);
     p3=p3->n;}} */
  //}

  //printf ("%u", prbx1->dim);

  itv_t a1_centre1_paf0;
  itv_init (a1_centre1_paf0);
  itv_t a1_centre2_paf1;
  itv_init (a1_centre2_paf1);
  itv_t a2_centre1_paf0;
  itv_init (a2_centre1_paf0);
  itv_t a2_centre2_paf1;
  itv_init (a2_centre2_paf1);

  itv_mul_num (a1_centre1_paf0, a1->paf[0]->c, num2);
  itv_mul_num (a1_centre2_paf1, a1->paf[1]->c, num2);
  itv_mul_num (a2_centre1_paf0, a2->paf[0]->c, num4);
  itv_mul_num (a2_centre2_paf1, a2->paf[1]->c, num4);
  itv_add (tpmeet->paf[0]->c, a1_centre1_paf0, a2_centre1_paf0);
  itv_add (tpmeet->paf[1]->c, a1_centre2_paf1, a2_centre2_paf1);

  /* update tpmeet->box */  
  for (i=0; i<tpmeet->dims; i++){ 
      //itv_meet(prbx1->itv, tpmeet->box[i], a1->box[i], a2->box[i]);
      t1p_aff_boxize (prbx1, tpmeet->box[i], tpmeet->paf[i], tpmeet);
      itv_set(tpmeet->paf[i]->itv, tpmeet->box[i]);
      tpmeet->paf[i]->pby++;}

  //scanf("%f",&b);
  //itv_print(tpmeet->paf[1]->end->coeff);
  //printf("%u", tpmeet->paf[1]->end->pnsym->index);
  //scanf("%f",&b);
//itv_print(tpmeet->paf[0]->q->n->coeff);
//itv_print(tpmeet->paf[0]->q->n->n->coeff);
//itv_print(tpmeet->paf[0]->end->coeff);
//printf("%u", tpmeet->paf[0]->end->pnsym->index);
//itv_print(tpmeet->paf[1]->q->coeff);
//printf("%u", tpmeet->paf[1]->q->pnsym->index);
//itv_print(tpmeet->paf[1]->q->n->coeff);
//itv_print(tpmeet->paf[1]->q->n->n->coeff);
//itv_print(tpmeet->paf[1]->end->coeff);
//printf("%u", tpmeet->paf[1]->end->pnsym->index);

  //itv_print (a1->paf[0]->q->coeff);

  /*free (ptitv1);
  free (ptitv2);
  free (ptitv3);
  free (ptitv4);*/
  free (ptitv_meet);
  
  #ifdef _T1P_DEBUG
    fprintf(stdout, "### RESULT of MEET ###\n");
    t1p_fprint(stdout, man, tpmeet, 0x0);
    fprintf(stdout, "### ### ###\n");
  #endif

  return tpmeet;

}

#endif

#if 0
/************************************************/
/* 1.Meet as in weighted minkowski sum (without GLPK)*/ ////////// Dampened oscillator ////////// Roux_FMSD15_Ex7.txt
/************************************************/

t1p_t *
t1p_meet (ap_manager_t * man, bool destructive, t1p_t * a1, t1p_t * a2)
/* TODO: destructive not used */
{
  //int b;
  //scanf("%f",&b);
  CALL ();
  t1p_internal_t *prbx1 = t1p_init_from_manager (man, AP_FUNID_MEET);
  #ifdef _T1P_DEBUG
    fprintf(stdout, "### MEET OPERANDS (destructive %d)###\n",destructive);
    t1p_fprint(stdout, man, a1, 0x0);
    t1p_fprint(stdout, man, a2, 0x0);
    fprintf(stdout, "### ### ###\n");
  #endif
  itv_t coeff;
  itv_init (coeff);
  num_t num;
  num_init (num);
  num_set_double (num, 0);
  itv_set_num (coeff, num);
  num_t num2;
  num_init (num2);
  num_set_double (num2, 0.3746);
  itv_t num2coeff;
  itv_init (num2coeff);
  itv_set_num (num2coeff,num2);
  num_t num4;
  num_init (num4);
  num_set_double (num4, 0.7082);
  itv_t num4coeff;
  itv_init (num4coeff);
  itv_set_num (num4coeff,num4);
  num_t num6;
  num_init (num6);
  num_set_double (num6, 0.0172);
  itv_t num6coeff;
  itv_init (num6coeff);
  itv_set_num (num6coeff,num6);
  num_t num7;
  num_init (num7);
  num_set_double (num7, 1.7385);
  itv_t num7coeff;
  itv_init (num7coeff);
  itv_set_num (num7coeff,num7);
  num_t num8;
  num_init (num8);
  num_set_double (num8, -0.0708);
  itv_t num8coeff;
  itv_init (num8coeff);
  itv_set_num (num8coeff,num8);
  num_t num9;
  num_init (num9);
  num_set_double (num9, 1.7007);
  itv_t num9coeff;
  itv_init (num9coeff);
  itv_set_num (num9coeff,num9);
  

  int i, j, size, size1, size2;

///////// first generator matrix ///////////////////

  /*itv_t **ptitv1 = (itv_t **) calloc (a1->dims, sizeof (itv_t *));
  for (i = 0; i < a1->dims; i++)
    ptitv1[i] = (itv_t *) calloc (prbx1->dim, sizeof (itv_t));

  size1 = t1p_aff_get_generator_matrix_size (prbx1, a1, ptitv1);

//printf("%u", size1);

  itv_t **ptitv2 = (itv_t **) calloc (a1->dims, sizeof (itv_t *));
  for (i = 0; i < a1->dims; i++)
    ptitv2[i] = (itv_t *) calloc (size1, sizeof (itv_t));

  t1p_aff_get_generator_matrix_num (prbx1, a1, ptitv1, ptitv2, num2);

//////////////// second generator matrix //////////////////

  itv_t **ptitv3 = (itv_t **) calloc (a2->dims, sizeof (itv_t *));
  for (i = 0; i < a2->dims; i++)
    ptitv3[i] = (itv_t *) calloc (prbx1->dim, sizeof (itv_t));

  size2 = t1p_aff_get_generator_matrix_size (prbx1, a2, ptitv3);

//printf("%u", size2);

  itv_t **ptitv4 = (itv_t **) calloc (a2->dims, sizeof (itv_t *));
  for (i = 0; i < a2->dims; i++)
    ptitv4[i] = (itv_t *) calloc (size2, sizeof (itv_t));

  t1p_aff_get_generator_matrix_num (prbx1, a2, ptitv3, ptitv4, num4);*/


  size = 4; //size1 + size2;

  itv_t **ptitv_meet = (itv_t **) calloc (a2->dims, sizeof (itv_t *));
  for (i = 0; i < a2->dims; i++)
    ptitv_meet[i] = (itv_t *) calloc (size, sizeof (itv_t));

/////////// fixed size buffer /////////////////

  //itv_t **ptitv_meet[a1->dims][size];

//printf("%u", size1);

//printf("%u", size2);

  //scanf("%f",&b);


 /* for (i = 0; i < a1->dims; i++) {
    for (j = 0; j < size; j++) {
      if (j >= (2)) {
	itv_set (ptitv_meet[i][j], ptitv4[i][j - size1]);
	//itv_print(ptitv_meet[i][j]);
      }
      else {
	itv_set (ptitv_meet[i][j], ptitv2[i][j]);
	//itv_print(ptitv_meet[i][j]);
      }
      //itv_print(ptitv[i][j]);
    }
  }*/

itv_set (ptitv_meet[0][0], num2coeff);
itv_set (ptitv_meet[0][1], coeff);
itv_set (ptitv_meet[0][2], num4coeff);
itv_set (ptitv_meet[0][3], num6coeff);
itv_set (ptitv_meet[1][0], coeff);
itv_set (ptitv_meet[1][1], num7coeff);
itv_set (ptitv_meet[1][2], num8coeff);
itv_set (ptitv_meet[1][3], num9coeff);
  //scanf("%f",&b);

//itv_print(ptitv3[0][2]);

//int k=2*prbx1->dim;
//prbx1->dim=0;

//printf("%u",prbx1->dim);

  t1p_t *tpmeet = t1p_alloc (man, 0, 2);
  tpmeet->paf[0] = t1p_aff_alloc_init (prbx1);
  tpmeet->paf[1] = t1p_aff_alloc_init (prbx1);

  t1p_nsym_t *nsymeet1;		//[size];
  //for (i=0;i<tb->dims;i++){
  t1p_aaterm_t *p1 = NULL;
  p1 = a1->paf[0]->q;
  /*t1p_aaterm_t *p2 = NULL;
     p2=a2->paf[0]->q; */
  t1p_aaterm_t *p3 = NULL;
  p3 = a1->paf[1]->q;

  for (j = 0; j < size; j++) {
    if (j >= size1) {
      //if (p2 != NULL){
      //t1p_aff_build (prbx1, tpmeet->paf[0], ptitv_meet[0][j], p2->pnsym->index);
      nsymeet1 = t1p_nsym_add (prbx1, UN);  //// IN is for eps and UN is for eta
      t1p_aff_nsym_add (prbx1, tpmeet->paf[0], ptitv_meet[0][j], nsymeet1);	//////// creating new noise symbols ///////////
      //scanf("%f",&b);
      //printf("%u", prbx1->dim);
      //scanf("%f",&b);
      t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], prbx1->dim - 1);	//// using the index of the new noise symbols created above /////////////////
      //printf("%u", p2->pnsym->index);
      //p2=p2->n;}
    }
    else {
      if (p1 != NULL && itv_is_eq (ptitv_meet[0][j], coeff) == 0) {
	t1p_aff_build (prbx1, tpmeet->paf[0], ptitv_meet[0][j], p1->pnsym->index);
	//printf("%u", p1->pnsym->index);
	p1 = p1->n;
      }
      if (p3 != NULL && itv_is_eq (ptitv_meet[1][j], coeff) == 0) {
	t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], p3->pnsym->index);
	//printf("%u", p3->pnsym->index);
	p3 = p3->n;
      }
    }

    //nsymeet1 = t1p_nsym_add(prbx1, IN);
    //t1p_aff_nsym_add(prbx1, tpmeet->paf[0], ptitv_meet[0][j], nsymeet1);    //////// creating new noise symbols ///////////
    //scanf("%f",&b);
    //printf("%u", prbx1->dim);
    //scanf("%f",&b);
    //t1p_aff_build(prbx1, tpmeet->paf[1], ptitv_meet[1][j], prbx1->dim-1);  //// using the index of the new noise symbols created above /////////////////

    //t1p_aff_build (prbx1, tpmeet->paf[0], ptitv_meet[0][j], j);  //////////// creating new noise symbols but indexed from j= to size ///////////////// 
  }

  /*t1p_aaterm_t *p3 = NULL;
     p3=a1->paf[1]->q;
     /*t1p_aaterm_t *p4 = NULL;
     p4=a2->paf[1]->q; */

  //for (j = 0; j < size; j++) {
  //t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], j);
  /*if (j >= size1){
     if (p4 !=NULL){
     t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], p4->pnsym->index);
     printf("%u", p4->pnsym->index);
     p4=p4->n;}}
     else {
     if (p3 !=NULL){
     t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], p3->pnsym->index);
     printf("%u", p3->pnsym->index);
     p3=p3->n;}} */
  //}

  //printf ("%u", prbx1->dim);

  itv_t a1_centre1_paf0;
  itv_init (a1_centre1_paf0);
  itv_t a1_centre2_paf1;
  itv_init (a1_centre2_paf1);
  itv_t a2_centre1_paf0;
  itv_init (a2_centre1_paf0);
  itv_t a2_centre2_paf1;
  itv_init (a2_centre2_paf1);

  itv_mul_num (a1_centre1_paf0, a1->paf[0]->c, num2);
  itv_mul_num (a1_centre2_paf1, a1->paf[1]->c, num2);
  itv_mul_num (a2_centre1_paf0, a2->paf[0]->c, num4);
  itv_mul_num (a2_centre2_paf1, a2->paf[1]->c, num4);
  itv_add (tpmeet->paf[0]->c, a1_centre1_paf0, a2_centre1_paf0);
  itv_add (tpmeet->paf[1]->c, a1_centre2_paf1, a2_centre2_paf1);

  /* update tpmeet->box */  
  for (i=0; i<tpmeet->dims; i++){ 
      //itv_meet(prbx1->itv, tpmeet->box[i], a1->box[i], a2->box[i]);
      t1p_aff_boxize (prbx1, tpmeet->box[i], tpmeet->paf[i], tpmeet);
      itv_set(tpmeet->paf[i]->itv, tpmeet->box[i]);
      tpmeet->paf[i]->pby++;}

  //scanf("%f",&b);
  //itv_print(tpmeet->paf[1]->end->coeff);
  //printf("%u", tpmeet->paf[1]->end->pnsym->index);
  //scanf("%f",&b);
//itv_print(tpmeet->paf[0]->q->n->coeff);
//itv_print(tpmeet->paf[0]->q->n->n->coeff);
//itv_print(tpmeet->paf[0]->end->coeff);
//printf("%u", tpmeet->paf[0]->end->pnsym->index);
//itv_print(tpmeet->paf[1]->q->coeff);
//printf("%u", tpmeet->paf[1]->q->pnsym->index);
//itv_print(tpmeet->paf[1]->q->n->coeff);
//itv_print(tpmeet->paf[1]->q->n->n->coeff);
//itv_print(tpmeet->paf[1]->end->coeff);
//printf("%u", tpmeet->paf[1]->end->pnsym->index);

  //itv_print (a1->paf[0]->q->coeff);

  /*free (ptitv1);
  free (ptitv2);
  free (ptitv3);
  free (ptitv4);*/
  free (ptitv_meet);
  
  #ifdef _T1P_DEBUG
    fprintf(stdout, "### RESULT of MEET ###\n");
    t1p_fprint(stdout, man, tpmeet, 0x0);
    fprintf(stdout, "### ### ###\n");
  #endif

  return tpmeet;

}
#endif


/************************************************/
/* 1.Meet as in weighted minkowski sum (without GLPK)*/
/************************************************/
#if 0

t1p_t *
t1p_meet (ap_manager_t * man, bool destructive, t1p_t * a1, t1p_t * a2)
/* TODO: destructive not used */
{

  if(t1p_is_top (man, a1)){
     return a2;}
  if(t1p_is_top (man, a2)){
    return a1;}
  int b;
  //scanf("%f",&b);*/
  CALL ();
  t1p_internal_t *prbx1 = t1p_init_from_manager (man, AP_FUNID_MEET);
  #ifdef _T1P_DEBUG
    fprintf(stdout, "### MEET OPERANDS (destructive %d)###\n",destructive);
    t1p_fprint(stdout, man, a1, 0x0);
    t1p_fprint(stdout, man, a2, 0x0);
    fprintf(stdout, "### ### ###\n");
  #endif
  scanf("%f",&b);
  itv_t coeff;
  itv_init (coeff);
  num_t num;
  num_init (num);
  num_set_double (num, 0);
  itv_set_num (coeff, num);
  num_t num2;
  num_init (num2);
  num_set_double (num2, 0.5);
  num_t num4;
  num_init (num4);
  num_set_double (num4, 0.5);

  int i, j, size, size1, size2;

///////// first generator matrix ///////////////////

  itv_t **ptitv1 = (itv_t **) calloc (a1->dims, sizeof (itv_t *));
  for (i = 0; i < a1->dims; i++)
    ptitv1[i] = (itv_t *) calloc (prbx1->dim, sizeof (itv_t));

  //size1 = t1p_aff_get_generator_matrix_size (prbx1, a1, ptitv1);

//printf("%u", size1);

/*for (i=0;i<a1->dims;i++){
   for (j=0;j<prbx1->dim;j++){
       itv_init(ptitv1[i][j]);//=coeff1;
       //itv_print(ptitv[i][j]);
     }
}*/

//itv_print(ptitv[0][0]);

t1p_aaterm_t *p = NULL;

int* vec= (int* )calloc(prbx1->dim,sizeof(int));

for (i=0;i<prbx1->dim;i++){
vec[i]=-1;
}

int k=0;

for(i=0;i<a1->dims;i++) {
 for(p=a1->paf[i]->q;p;p=p->n) {
		//itv_init(ptitv1[i][p->pnsym->index]);
		itv_set(ptitv1[i][p->pnsym->index],p->coeff);
                //itv_print(ptitv[i][p->pnsym->index]);
                if (vec[p->pnsym->index]==-1){
                   vec[p->pnsym->index]=k;
                   k++;}
   }
}

size1=k;

  itv_t **ptitv2 = (itv_t **) calloc (a1->dims, sizeof (itv_t *));
  for (i = 0; i < a1->dims; i++)
    ptitv2[i] = (itv_t *) calloc (size1, sizeof (itv_t));

  t1p_aff_get_generator_matrix_num (prbx1, a1, ptitv1, ptitv2, num2);

/*for (i=0;i<a1->dims;i++){
   for (j=0;j<k;j++){
       itv_init(ptitv2[i][j]);//=coeff1;
       //itv_print(ptitv[i][j]);
     }
}*/

t1p_aaterm_t *p1 = NULL;

itv_t tmp;

for(i=0;i<a1->dims;i++) {
 for(p1=a1->paf[i]->q;p1;p1=p1->n) {
		//itv_init(ptitv2[i][vec[p1->pnsym->index]]);
                itv_init(tmp);
                itv_mul_num(tmp, p1->coeff, num2);
		itv_set(ptitv2[i][vec[p1->pnsym->index]],tmp);
                //itv_print(ptitv1[i][p1->pnsym->index]);
                  
	   }
}

//////////////// second generator matrix //////////////////

  itv_t **ptitv3 = (itv_t **) calloc (a2->dims, sizeof (itv_t *));
  for (i = 0; i < a2->dims; i++)
    ptitv3[i] = (itv_t *) calloc (prbx1->dim, sizeof (itv_t));

  //size2 = t1p_aff_get_generator_matrix_size (prbx1, a2, ptitv3);

//printf("%u", size2);

/*for (i=0;i<a2->dims;i++){
   for (j=0;j<prbx1->dim;j++){
       itv_init(ptitv3[i][j]);//=coeff1;
       //itv_print(ptitv[i][j]);
     }
}*/

//itv_print(ptitv[0][0]);

//t1p_aaterm_t *p = NULL;

//int* vec= (int* )calloc(prbx1->dim,sizeof(int));

for (i=0;i<prbx1->dim;i++){
vec[i]=-1;
}

k=0;

for(i=0;i<a2->dims;i++) {
 for(p=a2->paf[i]->q;p;p=p->n) {
		//itv_init(ptitv3[i][p->pnsym->index]);
		itv_set(ptitv3[i][p->pnsym->index],p->coeff);
                //itv_print(ptitv[i][p->pnsym->index]);
                if (vec[p->pnsym->index]==-1){
                   vec[p->pnsym->index]=k;
                   k++;}
   }
}

size2=k;

  itv_t **ptitv4 = (itv_t **) calloc (a2->dims, sizeof (itv_t *));
  for (i = 0; i < a2->dims; i++)
    ptitv4[i] = (itv_t *) calloc (size2, sizeof (itv_t));

  t1p_aff_get_generator_matrix_num (prbx1, a2, ptitv3, ptitv4, num4);

/*for (i=0;i<a2->dims;i++){
   for (j=0;j<k;j++){
       itv_init(ptitv4[i][j]);//=coeff1;
       //itv_print(ptitv[i][j]);
     }
}*/

//t1p_aaterm_t *p1 = NULL;

//itv_t tmp;

for(i=0;i<a2->dims;i++) {
 for(p1=a2->paf[i]->q;p1;p1=p1->n) {
		//itv_init(ptitv4[i][vec[p1->pnsym->index]]);
                itv_init(tmp);
                itv_mul_num(tmp, p1->coeff, num4);
		itv_set(ptitv4[i][vec[p1->pnsym->index]],tmp);
                //itv_print(ptitv1[i][p1->pnsym->index]);
                  
	   }
}


  size = size1 + size2;

  itv_t **ptitv_meet = (itv_t **) calloc (a2->dims, sizeof (itv_t *));
  for (i = 0; i < a2->dims; i++)
    ptitv_meet[i] = (itv_t *) calloc (size, sizeof (itv_t));

/////////// fixed size buffer /////////////////

  //itv_t **ptitv_meet[a1->dims][size];

//printf("%u", size1);

//printf("%u", size2);

  //scanf("%f",&b);


  for (i = 0; i < a1->dims; i++) {
    for (j = 0; j < size; j++) {
      if (j >= (size1)) {
	itv_set (ptitv_meet[i][j], ptitv4[i][j - size1]);
	//itv_print(ptitv_meet[i][j]);
      }
      else {
	itv_set (ptitv_meet[i][j], ptitv2[i][j]);
	//itv_print(ptitv_meet[i][j]);
      }
      //itv_print(ptitv[i][j]);
    }
  }

  //scanf("%f",&b);

//itv_print(ptitv3[0][2]);

//int k=2*prbx1->dim;
//prbx1->dim=0;

//printf("%u",prbx1->dim);

  t1p_t *tpmeet = t1p_alloc (man, 0, 2);
  tpmeet->paf[0] = t1p_aff_alloc_init (prbx1);
  tpmeet->paf[1] = t1p_aff_alloc_init (prbx1);

  t1p_nsym_t *nsymeet1;		//[size];
  //for (i=0;i<tb->dims;i++){
  t1p_aaterm_t *p4 = NULL;
  p4 = a1->paf[0]->q;
  /*t1p_aaterm_t *p2 = NULL;
     p2=a2->paf[0]->q; */
  t1p_aaterm_t *p3 = NULL;
  p3 = a1->paf[1]->q;

  for (j = 0; j < size; j++) {
    if (j >= size1) {
      //if (p2 != NULL){
      //t1p_aff_build (prbx1, tpmeet->paf[0], ptitv_meet[0][j], p2->pnsym->index);
      nsymeet1 = t1p_nsym_add (prbx1, UN);  //// IN is for eps and UN is for eta
      t1p_aff_nsym_add (prbx1, tpmeet->paf[0], ptitv_meet[0][j], nsymeet1);	//////// creating new noise symbols ///////////
      //scanf("%f",&b);
      //printf("%u", prbx1->dim);
      //scanf("%f",&b);
      t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], prbx1->dim - 1);	//// using the index of the new noise symbols created above /////////////////
      //printf("%u", p2->pnsym->index);
      //p2=p2->n;}
    }
    else {
      if (p4 != NULL && itv_is_eq (ptitv_meet[0][j], coeff) == 0) {
	t1p_aff_build (prbx1, tpmeet->paf[0], ptitv_meet[0][j], p4->pnsym->index);
	//printf("%u", p1->pnsym->index);
	p4 = p4->n;
      }
      if (p3 != NULL && itv_is_eq (ptitv_meet[1][j], coeff) == 0) {
	t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], p3->pnsym->index);
	//printf("%u", p3->pnsym->index);
	p3 = p3->n;
      }
    }

    //nsymeet1 = t1p_nsym_add(prbx1, IN);
    //t1p_aff_nsym_add(prbx1, tpmeet->paf[0], ptitv_meet[0][j], nsymeet1);    //////// creating new noise symbols ///////////
    //scanf("%f",&b);
    //printf("%u", prbx1->dim);
    //scanf("%f",&b);
    //t1p_aff_build(prbx1, tpmeet->paf[1], ptitv_meet[1][j], prbx1->dim-1);  //// using the index of the new noise symbols created above /////////////////

    //t1p_aff_build (prbx1, tpmeet->paf[0], ptitv_meet[0][j], j);  //////////// creating new noise symbols but indexed from j= to size ///////////////// 
  }

  /*t1p_aaterm_t *p3 = NULL;
     p3=a1->paf[1]->q;
     /*t1p_aaterm_t *p4 = NULL;
     p4=a2->paf[1]->q; */

  //for (j = 0; j < size; j++) {
  //t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], j);
  /*if (j >= size1){
     if (p4 !=NULL){
     t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], p4->pnsym->index);
     printf("%u", p4->pnsym->index);
     p4=p4->n;}}
     else {
     if (p3 !=NULL){
     t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], p3->pnsym->index);
     printf("%u", p3->pnsym->index);
     p3=p3->n;}} */
  //}

  //printf ("%u", prbx1->dim);

  itv_t a1_centre1_paf0;
  itv_init (a1_centre1_paf0);
  itv_t a1_centre2_paf1;
  itv_init (a1_centre2_paf1);
  itv_t a2_centre1_paf0;
  itv_init (a2_centre1_paf0);
  itv_t a2_centre2_paf1;
  itv_init (a2_centre2_paf1);

  itv_mul_num (a1_centre1_paf0, a1->paf[0]->c, num2);
  itv_mul_num (a1_centre2_paf1, a1->paf[1]->c, num2);
  itv_mul_num (a2_centre1_paf0, a2->paf[0]->c, num4);
  itv_mul_num (a2_centre2_paf1, a2->paf[1]->c, num4);
  itv_add (tpmeet->paf[0]->c, a1_centre1_paf0, a2_centre1_paf0);
  itv_add (tpmeet->paf[1]->c, a1_centre2_paf1, a2_centre2_paf1);

  /* update tpmeet->box */  
  for (i=0; i<tpmeet->dims; i++){ 
      //itv_meet(prbx1->itv, tpmeet->box[i], a1->box[i], a2->box[i]);
      t1p_aff_boxize (prbx1, tpmeet->box[i], tpmeet->paf[i], tpmeet);
      itv_set(tpmeet->paf[i]->itv, tpmeet->box[i]);
      tpmeet->paf[i]->pby++;}

  //scanf("%f",&b);
  //itv_print(tpmeet->paf[1]->end->coeff);
  //printf("%u", tpmeet->paf[1]->end->pnsym->index);
  //scanf("%f",&b);
//itv_print(tpmeet->paf[0]->q->n->coeff);
//itv_print(tpmeet->paf[0]->q->n->n->coeff);
//itv_print(tpmeet->paf[0]->end->coeff);
//printf("%u", tpmeet->paf[0]->end->pnsym->index);
//itv_print(tpmeet->paf[1]->q->coeff);
//printf("%u", tpmeet->paf[1]->q->pnsym->index);
//itv_print(tpmeet->paf[1]->q->n->coeff);
//itv_print(tpmeet->paf[1]->q->n->n->coeff);
//itv_print(tpmeet->paf[1]->end->coeff);
//printf("%u", tpmeet->paf[1]->end->pnsym->index);

  //itv_print (a1->paf[0]->q->coeff);
   
  itv_clear (coeff);
  num_clear (num);
  num_clear (num2);
  num_clear (num4);
  itv_clear (tmp);
  for (i = 0; i < a1->dims; i++) {
    for (j = 0; j < prbx1->dim; j++) {
      itv_clear (ptitv1[i][j]);}}
  for (i = 0; i < a1->dims; i++) {
     free(ptitv1[i]);}
  free (ptitv1);
  for (i = 0; i < a1->dims; i++) {
    for (j = 0; j < size1; j++) {
      itv_clear (ptitv2[i][j]);}}
  for (i = 0; i < a1->dims; i++) {
     free(ptitv2[i]);}
  free (ptitv2);
  for (i = 0; i < a2->dims; i++) {
    for (j = 0; j < prbx1->dim; j++) {
      itv_clear (ptitv3[i][j]);}}
  for (i = 0; i < a2->dims; i++) {
     free(ptitv3[i]);}
  free (ptitv3);
  for (i = 0; i < a2->dims; i++) {
    for (j = 0; j < size2; j++) {
      itv_clear (ptitv4[i][j]);}}
  for (i = 0; i < a2->dims; i++) {
     free(ptitv4[i]);}
  free (ptitv4);
  for (i = 0; i < a2->dims; i++) {
    for (j = 0; j < size; j++) {
      itv_clear (ptitv_meet[i][j]);}}
  for (i = 0; i < a2->dims; i++) {
     free(ptitv_meet[i]);}
  free (ptitv_meet);

  
  #ifdef _T1P_DEBUG
    fprintf(stdout, "### RESULT of MEET ###\n");
    t1p_fprint(stdout, man, tpmeet, 0x0);
    fprintf(stdout, "### ### ###\n");
  #endif

  //scanf("%f",&b);

  return tpmeet;

}
#endif 

#if 0

/************************************************/
/* 1.Meet as in weighted minkowski sum (with GLPK)					*/
/************************************************/

t1p_t *
t1p_meet (ap_manager_t * man, bool destructive, t1p_t * a1, t1p_t * a2)
/* TODO: destructive not used */
{
  //int b;
  //scanf("%f",&b);
  #ifdef _T1P_DEBUG
    fprintf(stdout, "### MEET OPERANDS (destructive %d)###\n",destructive);
    t1p_fprint(stdout, man, a1, 0x0);
    t1p_fprint(stdout, man, a2, 0x0);
    fprintf(stdout, "### ### ###\n");
  #endif
  int i, j,status;  
  double numfx1_glp,numfx2_glp;
  itv_t coeff;
  itv_init (coeff);
  num_t num;
  num_init (num);
  num_set_double (num, 0);
  itv_set_num (coeff, num);
  num_t num1;
  num_init (num1);
  itv_t *diffc = (itv_t *) calloc (2, sizeof (itv_t));
  for (i = 0; i < a1->dims; i++) {
      itv_init (diffc[i]);
      itv_sub (diffc[i], a2->paf[i]->c, a1->paf[i]->c);}
  
  itv_magnitude (num1, diffc[0]);
  double_set_num (&numfx1_glp, num1);
  if (itv_is_neg (diffc[0]) == 1)
     numfx1_glp = -1 * numfx1_glp;
  itv_magnitude (num1, diffc[1]);
  double_set_num (&numfx2_glp, num1);
  if (itv_is_neg (diffc[1]) == 1)
     numfx2_glp = -1 * numfx2_glp;
  
  glp_prob *lp;
  s1:   lp = glp_create_prob();
  s2:   glp_set_obj_dir(lp, GLP_MIN);
  s3:   glp_add_rows(lp, 2);
  s4:   glp_set_row_bnds(lp, 1, GLP_FX, numfx1_glp, numfx1_glp);
  s5:   glp_set_row_bnds(lp, 2, GLP_FX, numfx2_glp, numfx2_glp);

  //printf("%f", numfx1_glp);
  //printf("%f", numfx2_glp);

  CALL ();
  t1p_internal_t *prbx1 = t1p_init_from_manager (man, AP_FUNID_MEET);
  
  //num_t num2;
  //num_init (num2);
  //num_set_double (num2, 0.4);
  num_t num4;
  num_init (num4);
  num_set_double (num4, 1.0);

  int size, size1, size2;

  ///////// first generator matrix ///////////////////

  itv_t **ptitv1 = (itv_t **) calloc (a1->dims, sizeof (itv_t *));
  for (i = 0; i < a1->dims; i++)
    ptitv1[i] = (itv_t *) calloc (prbx1->dim, sizeof (itv_t));

  size1 = t1p_aff_get_generator_matrix_size (prbx1, a1, ptitv1);

  //printf("%u", size1);

  itv_t **ptitv2 = (itv_t **) calloc (a1->dims, sizeof (itv_t *));
  for (i = 0; i < a1->dims; i++)
    ptitv2[i] = (itv_t *) calloc (size1, sizeof (itv_t));

  t1p_aff_get_generator_matrix(prbx1, a1, ptitv1, ptitv2);

  //////////////// second generator matrix //////////////////

  itv_t **ptitv3 = (itv_t **) calloc (a2->dims, sizeof (itv_t *));
  for (i = 0; i < a2->dims; i++)
    ptitv3[i] = (itv_t *) calloc (prbx1->dim, sizeof (itv_t));

  size2 = t1p_aff_get_generator_matrix_size (prbx1, a2, ptitv3);

  //printf("%u", size2);

  itv_t **ptitv4 = (itv_t **) calloc (a2->dims, sizeof (itv_t *));
  for (i = 0; i < a2->dims; i++)
    ptitv4[i] = (itv_t *) calloc (size2, sizeof (itv_t));

  t1p_aff_get_generator_matrix_num (prbx1, a2, ptitv3, ptitv4, num4);


  size = size1 + size2;

  itv_t **ptitv_meet = (itv_t **) calloc (a2->dims, sizeof (itv_t *));
  for (i = 0; i < a2->dims; i++)
    ptitv_meet[i] = (itv_t *) calloc (size, sizeof (itv_t));

  for (i = 0; i < a1->dims; i++) {
    for (j = 0; j < size; j++) {
      if (j >= (size1)) {
	itv_set (ptitv_meet[i][j], ptitv4[i][j - size1]);
	//itv_print(ptitv_meet[i][j]);
      }
      else {
	itv_set (ptitv_meet[i][j], ptitv2[i][j]);
	//itv_print(ptitv_meet[i][j]);
      }
      //itv_print(ptitv[i][j]);
    }
  }

  s6:  glp_add_cols(lp, size);
  for (j = 0; j < size; j++) {
  s7:  glp_set_col_bnds(lp, j+1, GLP_DB, -1.0, 1.0);}
  
  int ia[1+1000], ja[1+1000];
  double *ar;
  ar = calloc(1+1000, sizeof(double));
  for (i = 0; i < a1->dims; i++) {
    for (j = 0; j < size; j++) {
        itv_magnitude (num1, ptitv_meet[i][j]);
        double_set_num (&numfx1_glp, num1);
        if (itv_is_neg (ptitv_meet[i][j]) == 1){
           numfx1_glp = -1 * numfx1_glp;}
        if (i==0){
           s8:  ia[1+j] = i+1, ja[1+j] = j+1, ar[1+j] =  numfx1_glp;}
           //printf("%f", ar[1+j]);}
        if (i==1){
           s9:  ia[1+j+size] = i+1, ja[1+j+size] = j+1, ar[1+j+size] =  numfx1_glp;}
           //printf("%f", ar[1+j+size]);}
    }
  }
  
  ////// 2*size because the dimension is 2 ///////////////
  s10:  glp_load_matrix(lp, 2*size, ia, ja, ar);
  s11:  glp_simplex(lp, NULL);
  s12:  status=glp_get_status(lp);

  //printf("%u", status);
  
  if (status!=2 && status!=5){
     t1p_t* bot=t1p_bottom(man, 0, 2);
     s13: glp_delete_prob(lp);
     return bot;}
  else{
      num_t num2;
      num_init (num2);
      num_set_double (num2, 0.4);
      num_t num4;
      num_init (num4);
      num_set_double (num4, 0.6);
      t1p_aff_get_generator_matrix_num(prbx1, a1, ptitv1, ptitv2, num2);
      t1p_aff_get_generator_matrix_num (prbx1, a2, ptitv3, ptitv4, num4);
      for (i = 0; i < a1->dims; i++) {
         for (j = 0; j < size; j++) {
            if (j >= (size1)) {
	       itv_set (ptitv_meet[i][j], ptitv4[i][j - size1]);
	       //itv_print(ptitv_meet[i][j]);
             }
            else {
	       itv_set (ptitv_meet[i][j], ptitv2[i][j]);
	       //itv_print(ptitv_meet[i][j]);
             }
            //itv_print(ptitv[i][j]);
          }
       }
  
  //scanf("%f",&b);

//itv_print(ptitv3[0][2]);

//int k=2*prbx1->dim;
//prbx1->dim=0;

//printf("%u",prbx1->dim);

  t1p_t *tpmeet = t1p_alloc (man, 0, 2);
  tpmeet->paf[0] = t1p_aff_alloc_init (prbx1);
  tpmeet->paf[1] = t1p_aff_alloc_init (prbx1);


  t1p_nsym_t *nsymeet1;		//[size];
  //for (i=0;i<tb->dims;i++){
  t1p_aaterm_t *p1 = NULL;
  p1 = a1->paf[0]->q;
  /*t1p_aaterm_t *p2 = NULL;
     p2=a2->paf[0]->q; */
  t1p_aaterm_t *p3 = NULL;
  p3 = a1->paf[1]->q;

  for (j = 0; j < size; j++) {
    if (j >= size1) {
      //if (p2 != NULL){
      //t1p_aff_build (prbx1, tpmeet->paf[0], ptitv_meet[0][j], p2->pnsym->index);
      nsymeet1 = t1p_nsym_add (prbx1, UN);  //// IN is for eps and UN is for eta
      t1p_aff_nsym_add (prbx1, tpmeet->paf[0], ptitv_meet[0][j], nsymeet1);	//////// creating new noise symbols ///////////
      //scanf("%f",&b);
      //printf("%u", prbx1->dim);
      //scanf("%f",&b);
      t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], prbx1->dim - 1);	//// using the index of the new noise symbols created above /////////////////
      //printf("%u", p2->pnsym->index);
      //p2=p2->n;}
    }
    else {
      if (p1 != NULL && itv_is_eq (ptitv_meet[0][j], coeff) == 0) {
	t1p_aff_build (prbx1, tpmeet->paf[0], ptitv_meet[0][j], p1->pnsym->index);
	//printf("%u", p1->pnsym->index);
	p1 = p1->n;
      }
      if (p3 != NULL && itv_is_eq (ptitv_meet[1][j], coeff) == 0) {
	t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], p3->pnsym->index);
	//printf("%u", p3->pnsym->index);
	p3 = p3->n;
      }
    }

    //nsymeet1 = t1p_nsym_add(prbx1, IN);
    //t1p_aff_nsym_add(prbx1, tpmeet->paf[0], ptitv_meet[0][j], nsymeet1);    //////// creating new noise symbols ///////////
    //scanf("%f",&b);
    //printf("%u", prbx1->dim);
    //scanf("%f",&b);
    //t1p_aff_build(prbx1, tpmeet->paf[1], ptitv_meet[1][j], prbx1->dim-1);  //// using the index of the new noise symbols created above /////////////////

    //t1p_aff_build (prbx1, tpmeet->paf[0], ptitv_meet[0][j], j);  //////////// creating new noise symbols but indexed from j= to size ///////////////// 
  }

  /*t1p_aaterm_t *p3 = NULL;
     p3=a1->paf[1]->q;
     /*t1p_aaterm_t *p4 = NULL;
     p4=a2->paf[1]->q; */

  //for (j = 0; j < size; j++) {
  //t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], j);
  /*if (j >= size1){
     if (p4 !=NULL){
     t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], p4->pnsym->index);
     printf("%u", p4->pnsym->index);
     p4=p4->n;}}
     else {
     if (p3 !=NULL){
     t1p_aff_build (prbx1, tpmeet->paf[1], ptitv_meet[1][j], p3->pnsym->index);
     printf("%u", p3->pnsym->index);
     p3=p3->n;}} */
  //}

  //printf ("%u", prbx1->dim);

  itv_t a1_centre1_paf0;
  itv_init (a1_centre1_paf0);
  itv_t a1_centre2_paf1;
  itv_init (a1_centre2_paf1);
  itv_t a2_centre1_paf0;
  itv_init (a2_centre1_paf0);
  itv_t a2_centre2_paf1;
  itv_init (a2_centre2_paf1);

  itv_mul_num (a1_centre1_paf0, a1->paf[0]->c, num2);
  itv_mul_num (a1_centre2_paf1, a1->paf[1]->c, num2);
  itv_mul_num (a2_centre1_paf0, a2->paf[0]->c, num4);
  itv_mul_num (a2_centre2_paf1, a2->paf[1]->c, num4);
  itv_add (tpmeet->paf[0]->c, a1_centre1_paf0, a2_centre1_paf0);
  itv_add (tpmeet->paf[1]->c, a1_centre2_paf1, a2_centre2_paf1);

  /* update tpmeet->box */  
  for (i=0; i<tpmeet->dims; i++){ 
      //itv_meet(prbx1->itv, tpmeet->box[i], a1->box[i], a2->box[i]);
      t1p_aff_boxize (prbx1, tpmeet->box[i], tpmeet->paf[i], tpmeet);
      itv_set(tpmeet->paf[i]->itv, tpmeet->box[i]);
      tpmeet->paf[i]->pby++;}

  //scanf("%f",&b);
  //itv_print(tpmeet->paf[1]->end->coeff);
  //printf("%u", tpmeet->paf[1]->end->pnsym->index);
  //scanf("%f",&b);
//itv_print(tpmeet->paf[0]->q->n->coeff);
//itv_print(tpmeet->paf[0]->q->n->n->coeff);
//itv_print(tpmeet->paf[0]->end->coeff);
//printf("%u", tpmeet->paf[0]->end->pnsym->index);
//itv_print(tpmeet->paf[1]->q->coeff);
//printf("%u", tpmeet->paf[1]->q->pnsym->index);
//itv_print(tpmeet->paf[1]->q->n->coeff);
//itv_print(tpmeet->paf[1]->q->n->n->coeff);
//itv_print(tpmeet->paf[1]->end->coeff);
//printf("%u", tpmeet->paf[1]->end->pnsym->index);

  //itv_print (a1->paf[0]->q->coeff);

  free (ptitv1);
  free (ptitv2);
  free (ptitv3);
  free (ptitv4);
  free (ptitv_meet);
  
  #ifdef _T1P_DEBUG
    fprintf(stdout, "### RESULT of MEET ###\n");
    t1p_fprint(stdout, man, tpmeet, 0x0);
    fprintf(stdout, "### ### ###\n");
  #endif

  return tpmeet;

} // else part 

}

#endif

#if 0

t1p_t* t1p_meet(ap_manager_t* man, bool destructive, t1p_t* a1, t1p_t* a2)
/* TODO: destructive not used */
{
    CALL();
    t1p_internal_t* pr = t1p_init_from_manager(man,AP_FUNID_MEET);
    arg_assert(a1->dims==a2->dims && a1->intdim==a2->intdim,abort(););
#ifdef _T1P_DEBUG
    fprintf(stdout, "### MEET OPERANDS (destructive %d)###\n",destructive);
    t1p_fprint(stdout, man, a1, 0x0);
    t1p_fprint(stdout, man, a2, 0x0);
    fprintf(stdout, "### ### ###\n");
#endif
    t1p_t* res = NULL;
    size_t i = 0;
    size_t realdim = a1->dims - a1->intdim;
    size_t intdim = a1->intdim;
    man->result.flag_best = tbool_true;
    man->result.flag_exact = tbool_true;
    if (t1p_is_eq(man, a1, a2)) res = t1p_copy(man, a1);
    else if (tbool_or(t1p_is_bottom(man, a1), t1p_is_bottom(man, a2)) == tbool_true) {
	res = t1p_bottom(man, intdim, realdim);
    } else if (t1p_is_top(man, a1) == tbool_true) {
	res = t1p_copy(man, a2);
    } else if (t1p_is_top(man, a2) == tbool_true) {
	res = t1p_copy(man, a1);
    } else {
	res = t1p_alloc(man, intdim, realdim);
	//ap_abstract0_free(pr->manNS, res->abs);
	size_t k = 0;
	ap_dim_t j = 0;
	size_t dims1 = t1p_nsymcons_get_dimension(pr, a1);
	size_t dims2 = t1p_nsymcons_get_dimension(pr, a2);
	if (dims1 && dims2) {
	    size_t dim2 = 0;
	    /* au pire il y a toute la liste de l'un a rajouter dans l'autre */
	    ap_dimchange_t* dimchange2 = ap_dimchange_alloc(0, dims2);
		if (dims1 > res->size) {
		    res->nsymcons = (ap_dim_t*)realloc(res->nsymcons, (dims1)*sizeof(ap_dim_t));
		    res->gamma = (ap_interval_t**)realloc(res->gamma, (dims1)*sizeof(ap_interval_t*));
		    for (k=res->size;k<dims1;k++) res->gamma[k] = NULL;
		    res->size = dims1;
		}
	    res->nsymcons = memcpy((void *)res->nsymcons, (const void *)a1->nsymcons, dims1*sizeof(ap_dim_t));
	    for (i=0; i<dims1; i++) res->gamma[i] = ap_interval_alloc_set(a1->gamma[i]);
	    ap_abstract0_free(pr->manNS, res->abs);
	    res->abs = ap_abstract0_copy(pr->manNS, a1->abs);
	    for (k=0; k<dims1; k++) {
		if (!t1p_nsymcons_get_dimpos(pr, &j, a1->nsymcons[k], a2)) {
		    dimchange2->dim[dim2] = j;
		    dim2++;
		}
	    }
	    dimchange2->realdim = dim2;
	    size_t dim1 = 0;
	    for (k=0; k<dims2; k++) t1p_insert_constrained_nsym(pr, &j, a2->nsymcons[k], res);

	    /* destructive, without projection (new dimension set to top) */
	    ap_abstract0_add_dimensions(pr->manNS, true, a2->abs, dimchange2, false);

	    ap_abstract0_meet(pr->manNS, true, res->abs, a2->abs);

	    /* update res->gamma */
	    t1p_update_nsymcons_gamma(pr, res);

	    ap_dimchange_add_invert(dimchange2);
	    ap_abstract0_remove_dimensions(pr->manNS, true, a2->abs, dimchange2);

	    dimchange2->realdim = dims2; ap_dimchange_free(dimchange2);
	} else if (dims1) {
		if (dims1 > res->size) {
		    res->nsymcons = (ap_dim_t*)realloc(res->nsymcons, (dims1)*sizeof(ap_dim_t));
		    res->gamma = (ap_interval_t**)realloc(res->gamma, (dims1)*sizeof(ap_interval_t*));
		    for (k=res->size;k<dims1;k++) res->gamma[k] = NULL;
		    res->size = dims1;
		}
	    res->nsymcons = memcpy((void *)res->nsymcons, (const void *)a1->nsymcons, dims1*sizeof(ap_dim_t));
	    for (i=0; i<dims1; i++) res->gamma[i] = ap_interval_alloc_set(a1->gamma[i]);
	    res->abs = ap_abstract0_copy(pr->manNS, a1->abs);
	} else if (dims2) {
		if (dims2 > res->size) {
		    res->nsymcons = (ap_dim_t*)realloc(res->nsymcons, (dims2)*sizeof(ap_dim_t));
		    res->gamma = (ap_interval_t**)realloc(res->gamma, (dims2)*sizeof(ap_interval_t*));
		    for (k=res->size;k<dims2;k++) res->gamma[k] = NULL;
		    res->size = dims2;
		}
	    res->nsymcons = memcpy((void *)res->nsymcons, (const void *)a2->nsymcons, dims2*sizeof(ap_dim_t));
	    for (i=0; i<dims2; i++) res->gamma[i] = ap_interval_alloc_set(a2->gamma[i]);
	    res->abs = ap_abstract0_copy(pr->manNS, a2->abs);
	} else {
	    /* non constrained abstract objects */
	}

	//res->abs = ap_abstract0_meet(pr->manNS, false, a1->abs, a2->abs);
	/* update res->gamma */
	//t1p_update_nsymcons_gamma(pr, res);
	itv_t tmp; itv_init(tmp);
	for (i=0; i<intdim+realdim; i++) {
	    /* update res->box */
	    itv_meet(pr->itv, res->box[i], a1->box[i], a2->box[i]);
	    if ((a1->paf[i]->q == NULL) && (a2->paf[i]->q == NULL)) {
		itv_meet(pr->itv, tmp, a1->paf[i]->c, a2->paf[i]->c);
		if (itv_is_bottom(pr->itv, tmp)) {
		    t1p_free(man, res);
		    res = t1p_bottom(man, intdim, realdim);
		    break;
		} else if (itv_is_top(tmp)) res->paf[i] = pr->top;
		else if (itv_has_finite_bound(tmp)) {
		    res->paf[i] = t1p_aff_alloc_init(pr);
		    t1p_aff_add_itv(pr, res->paf[i], tmp, IN);
		} else {
		    res->paf[i] = t1p_aff_alloc_init(pr);
		    itv_set(res->paf[i]->c, tmp);
		}
		res->paf[i]->pby++;
	    } else if (t1p_aff_is_eq(pr, a1->paf[i], a2->paf[i])) {
		res->paf[i] = a1->paf[i];
		res->paf[i]->pby++;
	    } else if (itv_is_point(pr->itv, res->box[i])) {
		res->paf[i] = t1p_aff_alloc_init(pr);
		itv_set(res->paf[i]->c, res->box[i]);
	    } else {
		fprintf(stdout, "destructive ? %d\n",destructive);
		t1p_fprint(stdout, man, a1, 0x0);
		t1p_fprint(stdout, man, a2, 0x0);
		//not_implemented();
		/* return a top instead */
		res->paf[i] = pr->top;
		res->paf[i]->pby++;
		itv_set_top(res->box[i]);
	    }
	}
	itv_clear(tmp);
    }
#ifdef _T1P_DEBUG
    fprintf(stdout, "### RESULT of MEET ###\n");
    t1p_fprint(stdout, man, res, 0x0);
    fprintf(stdout, "### ### ###\n");
#endif
    return res;
}

#endif

/**************/
/* meet array */
/**************/
t1p_t *
t1p_meet_array (ap_manager_t * man, t1p_t ** tab, size_t size)
{
  CALL ();
  return (t1p_t *) ap_generic_meet_array (man, (void **) tab, size);
}

static inline bool
t1p_eval_lincons_array (t1p_internal_t * pr, size_t * hash,
			t1p_aff_t ** eps_linexp_array,
			ap_lincons0_array_t * eps_array,
			ap_lincons0_array_t * var_lincons_array,
			ap_lincons0_array_t * array, t1p_t * a)
{
  CALL ();
  /* Requires: an array of varibale constraints */
  /* Ensures: an array of epsilons constraints */
  /* Ensures: an array of var constraints */
  bool res = true;
  size_t i = 0;
  size_t j = 0;
  size_t eps_dim = 0;		/* count the number of constraints that involves nsym */
  size_t var_dim = 0;		/* count the number of constraints that doesn't involve nsym, either simplified when evaluating the constraint or doesn't exist from the begining */
  ap_linexpr0_t *linexpr0 = NULL;
  t1p_aff_t *tmp1, *tmp2, *tmp3, *sum;
  for (i = 0; i < array->size; i++) {
    if (ap_lincons0_is_unsat (&(array->p[i]))) {
      res = false;
      break;			/* break the for loop */
    }
    else {
      /* recuperer l'expression lineaire et l'interpreter en eps */
      linexpr0 = array->p[i].linexpr0;
      sum = t1p_aff_alloc_init (pr);
      itv_set_ap_coeff (pr->itv, sum->c, &(linexpr0->cst));
      tmp1 = t1p_aff_alloc_init (pr);
      switch (linexpr0->discr) {
      case AP_LINEXPR_DENSE:
	for (j = 0; j < linexpr0->size; j++) {
	  /* j corresponds to the dimension of the variable */
	  if (ap_coeff_zero (&(linexpr0->p.coeff[j]))) {;
	  }
	  else {
	    itv_set_ap_coeff (pr->itv, tmp1->c, &(linexpr0->p.coeff[j]));
	    tmp2 = t1p_aff_mul_constrained (pr, tmp1, a->paf[j], a);
	    tmp3 = t1p_aff_add (pr, tmp2, sum, a);
	    t1p_aff_free (pr, sum);
	    sum = tmp3;
	    t1p_aff_free (pr, tmp1);
	    t1p_aff_free (pr, tmp2);
	  }
	}
	break;
      case AP_LINEXPR_SPARSE:
	for (j = 0; j < linexpr0->size; j++) {
	  /* j corresponds to the index of the linterm */
	  if (ap_coeff_zero (&(linexpr0->p.linterm[j].coeff))) {;
	  }
	  else {
	    itv_set_ap_coeff (pr->itv, tmp1->c, &(linexpr0->p.linterm[j].coeff));
	    tmp2 =
	      t1p_aff_mul_constrained (pr, tmp1, a->paf[linexpr0->p.linterm[j].dim],
				       a);
	    tmp3 = t1p_aff_add (pr, tmp2, sum, a);
	    t1p_aff_free (pr, sum);
	    sum = tmp3;
	    t1p_aff_free (pr, tmp1);
	    t1p_aff_free (pr, tmp2);
	  }
	}
	break;
      default:
	fatal ("unknown linexpr0 Discriminant");
	break;
      }
      eps_linexp_array[i] = sum;
      if (sum->l == 0) {
	/* this is a variable constraint */
	/* TODO: now duplicate, not optimal  */
	if (ap_lincons0_is_unsat (&array->p[i])) {
	  res = false;
	  break;
	}
	else {
	  var_lincons_array->p[var_dim] = ap_lincons0_copy (&(array->p[i]));
	  var_dim++;
	}
      }
      else {
	/* transformer cette forme en une contrainte lineaire (DENSE ou SPARSE a voir) */
	/* "sum" contains an affine form to transform into a constraint */
	ap_linexpr0_t *eps_linexpr0 = NULL;
	ap_coeff_t *coeff = ap_coeff_alloc (AP_COEFF_INTERVAL);
	ap_coeff_set_itv (pr->itv, coeff, sum->c);
	eps_linexpr0 = ap_linexpr0_alloc (AP_LINEXPR_SPARSE, 0);
	ap_linexpr0_set_cst (eps_linexpr0, coeff);
	t1p_aaterm_t *p;
	/* sum->l donne la taille de la forme affine (ne comptabilise pas la constante) */
	size_t epscons_size = t1p_nsymcons_get_dimension (pr, a);
	/* au pire tous les symboles de la contrainte seront \E0 rajouter */
	size_t new = 0;
	eps_linexpr0->p.linterm =
	  (ap_linterm_t *) malloc (sum->l * sizeof (ap_linterm_t));
	eps_linexpr0->size = sum->l;
	size_t k = 0;
	ap_dim_t dim = 0;
	for (p = sum->q; p; p = p->n) {
	  ap_coeff_init (&eps_linexpr0->p.linterm[k].coeff, AP_COEFF_INTERVAL);
	  ap_coeff_set_itv (pr->itv, &eps_linexpr0->p.linterm[k].coeff, p->coeff);
	  t1p_insert_constrained_nsym (pr, &dim, p->pnsym->index, a);
	  eps_linexpr0->p.linterm[k].dim = dim;
	  k++;
	}
	ap_scalar_t *scalar = NULL;
	if (array->p[i].scalar)
	  scalar = ap_scalar_alloc_set (array->p[i].scalar);
	eps_array->p[eps_dim] =
	  ap_lincons0_make (array->p[i].constyp, eps_linexpr0, scalar);
	hash[eps_dim] = i;
	eps_dim++;
	ap_coeff_free (coeff);
      }
    }
  }
  ap_lincons0_array_resize (eps_array, eps_dim);
  ap_lincons0_array_resize (var_lincons_array, var_dim);
  return res;
}

/**********************/
/* meet lincons array */
/**********************/
t1p_t *
t1p_meet_lincons_array (ap_manager_t * man, bool destructive, t1p_t * a,
			ap_lincons0_array_t * array)
{
  CALL ();
  t1p_internal_t *pr = t1p_init_from_manager (man, AP_FUNID_MEET_LINCONS_ARRAY);
  arg_assert (a && array, abort ();
    );
#ifdef _T1P_DEBUG
  fprintf (stdout, "### MEET LINCONS ARRAY (des %d) ###\n", destructive);
  t1p_fprint (stdout, man, a, NULL);
  ap_lincons0_array_fprint (stdout, array, NULL);
  fprintf (stdout, "### ### ###\n");
#endif

  size_t i = 0;
  ap_tcons0_array_t tcons0_array;
  tcons0_array.size = array->size;
  tcons0_array.p = (ap_tcons0_t *) malloc (array->size * sizeof (ap_tcons0_t));
  for (i = 0; i < array->size; i++) {
    tcons0_array.p[i] = ap_tcons0_from_lincons0 (&array->p[i]);
  }

  t1p_t *res = t1p_meet_tcons_array (man, destructive, a, &tcons0_array);
#ifdef _T1P_DEBUG
  fprintf (stdout, "### RESULT OF MEET LINCONS ARRAY (des %d) ###\n", destructive);
  t1p_fprint (stdout, man, res, NULL);
  fprintf (stdout, "### ### ###\n");
#endif
  return res;
}

/********************/
/* meet tcons array */
/********************/
t1p_t *
t1p_meet_tcons_array (ap_manager_t * man, bool destructive, t1p_t * a,
		      ap_tcons0_array_t * array)
{
  CALL ();
  t1p_internal_t *pr = t1p_init_from_manager (man, AP_FUNID_MEET_TCONS_ARRAY);
/*#ifdef _T1P_DEBUG
  fprintf (stdout, "### MEET TCONS ARRAY (des %d) ###\n", destructive);
  t1p_fprint (stdout, man, a, NULL);
  ap_tcons0_array_fprint (stdout, array, NULL);
  fprintf (stdout, "### ### ###\n");
#endif*/
  size_t i = 0;

  /*
     size_t size = array1->size;
     ap_tcons0_array_t arraybis = ap_tcons0_array_make(array1->size);
     for (i=0; i<array1->size; i++) {
     arraybis.p[i] = ap_tcons0_copy(&array1->p[i]);
     }
     ap_tcons0_array_t* array = &arraybis;
   */

  t1p_t *res = destructive ? a : t1p_copy (man, a);
  bool is_bottom = false;
  size_t kmax = 2;		/* specifies the maximum number of iterations */

  bool *tchange = (bool *) calloc (2 * a->dims, sizeof (bool));
  itv_t box;
  itv_init (box);
  itv_lincons_array_t itv_lincons_array;
  itv_lincons_array.size = 0;
  itv_lincons_array.p = NULL;
  if (!itv_intlinearize_ap_tcons0_array
      (pr->itv, &itv_lincons_array, array, res->box, res->intdim)) {
    size_t kmax = 2;		/* specifies the maximum number of iterations */
    /* intervalonly is set to false which means try to improve all dimensions, not only the ones with an interval coefficient */
    if (itv_boxize_lincons_array
	(pr->itv, res->box, tchange, &itv_lincons_array, res->box, res->intdim, kmax,
	 false)) {
      /* there is some inferred bounds */
      for (i = 0; i < res->dims; i++) {
	if (tchange[2 * i] || tchange[2 * i + 1]) {
	  if (itv_is_bottom (pr->itv, res->box[i])) {
	    is_bottom = true;
	    break;
	  }
	  else if (res->paf[i]->q == NULL) {
	    t1p_aff_check_free (pr, res->paf[i]);
	    res->paf[i] = t1p_aff_alloc_init (pr);

	    if (itv_has_finite_bound (res->box[i])) {
	      t1p_aff_add_itv (pr, res->paf[i], res->box[i], IN);
	    }
	    else {
	      itv_set (res->paf[i]->c, res->box[i]);
	    }

	    res->paf[i]->pby++;
	  }			/*else if (itv_is_point(pr->itv, res->box[i])) {
				   ap_tcons0_array_resize(array, 1+(array->size));
				   ap_interval_t* ap_itv = ap_interval_alloc();
				   ap_interval_set_itv(pr->itv, ap_itv, res->box[i]);
				   array->p[-1+(array->size)] = ap_tcons0_make(AP_CONS_EQ,ap_texpr0_binop(AP_TEXPR_SUB,ap_texpr0_dim(i),ap_texpr0_cst_interval(ap_itv),AP_RTYPE_REAL,AP_RDIR_UP),NULL);
				   ap_interval_free(ap_itv);
				   } */
	}
      }
    }
    else {
      /* nothing change */
    }
    itv_lincons_array_clear (&itv_lincons_array);
  }
  else {
    /* bottom */
    is_bottom = true;
  }
  if (!is_bottom) {
    /* texpr -> aff forme */
    ap_texpr0_t *texpr0 = NULL;
    t1p_aff_t **aff = malloc (array->size * sizeof (t1p_aff_t *));
    ap_linexpr0_t *linexpr0 = NULL;
    ap_lincons0_t lincons0;
    ap_interval_t cst;
    for (i = 0; i < a->dims; i++)
      itv_set (a->paf[i]->itv, a->box[i]);
    for (i = 0; i < array->size; i++) {
      texpr0 = array->p[i].texpr0;
      aff[i] = t1p_aff_eval_ap_texpr0 (pr, texpr0, a);
      if (aff[i]->q == NULL) {
	/* only the centers are involved in this constraint, already treated while updating res->box */
      }
      else {
	/* infer constraints on noise symbols */
	linexpr0 = t1p_ap_linexpr0_set_aff (pr, aff[i], res);
	lincons0.constyp = array->p[i].constyp;
	lincons0.linexpr0 = linexpr0;
	lincons0.scalar = array->p[i].scalar;
	ap_lincons0_array_t eps_lincons_array;
	eps_lincons_array.size = 1;
	eps_lincons_array.p = &lincons0;
	ap_abstract0_meet_lincons_array (pr->manNS, true, res->abs,
					 &eps_lincons_array);
	if (ap_abstract0_is_bottom (pr->manNS, res->abs)) {
	  is_bottom = true;
	  break;
	}
      }
      ap_linexpr0_free (linexpr0);
    }
    /* update res->gamma */
    t1p_update_nsymcons_gamma (pr, res);

    for (i = 0; i < array->size; i++) {
      /* update the abstract object with the new affine forms */
      if (array->p[i].constyp == AP_CONS_EQ) {
	itv_t dummy;
	itv_init (dummy);
	t1p_aff_t *tmp, *tmp1;
	size_t j = 0;
	for (j = 0; j < res->dims; j++) {
	  if (res->paf[j]->q) {
	    t1p_aff_cons_eq_lambda (pr, &dummy, res->paf[j], aff[i], res);
	    tmp = t1p_aff_mul_itv (pr, aff[i], dummy);
	    itv_set (res->paf[j]->itv, res->box[j]);
	    tmp1 = t1p_aff_add (pr, tmp, res->paf[j], res);
	    t1p_aff_check_free (pr, res->paf[j]);
	    res->paf[j] = tmp1;
	    /* update res->box */
	    itv_meet (pr->itv, res->box[j], res->box[j], res->paf[j]->itv);
	    res->paf[j]->pby++;
	    t1p_aff_free (pr, tmp);
	  }
	}
	itv_clear (dummy);
      }
      else {
	/* do nothing, just update res->box */
	size_t j = 0;
	for (j = 0; j < res->dims; j++) {
	  t1p_aff_bound (pr, box, res->paf[j], res);
	  itv_meet (pr->itv, res->box[j], res->box[j], box);
	}

      }
      t1p_aff_check_free (pr, aff[i]);
    }
    free (aff);
  }
  if (is_bottom) {
    size_t intdim = res->intdim;
    size_t realdim = res->dims - res->intdim;
    t1p_free (man, res);
    res = t1p_bottom (man, intdim, realdim);
  }
  //ap_tcons0_array_clear(array);
  free (tchange);
  itv_clear (box);
/*#ifdef _T1P_DEBUG
  fprintf (stdout, "### RESULT OF MEET TCONS ARRAY (des %d) ###\n", destructive);
  t1p_fprint (stdout, man, res, NULL);
  fprintf (stdout, "### ### ###\n");
#endif*/
  return res;
}

/************************************************/
/* 2.Join					*/
/************************************************/
/* local join */
t1p_t *
t1p_join (ap_manager_t * man, bool destructive, t1p_t * a1, t1p_t * a2)
    /* TODO destructive not used  */
{
  CALL ();
  size_t i = 0;
  t1p_internal_t *pr = (t1p_internal_t *) t1p_init_from_manager (man, AP_FUNID_JOIN);
  arg_assert (a1->dims == a2->dims && a1->intdim == a2->intdim, abort ();
    );
#ifdef _T1P_DEBUG
  fprintf (stdout, "### JOIN OPERANDS (des %d) (%tx U %tx) ###\n", destructive,
	   (intptr_t) a1, (intptr_t) a2);
  t1p_fprint (stdout, man, a1, NULL);
  t1p_fprint (stdout, man, a2, NULL);
  fprintf (stdout, "### ### ###\n");
#endif
  /* FILE* stream = fopen("/home/donquijote/taylor1p/taylor1plus/taylor1plus/demo/log", "a+"); */
  /* fprintf(stream,"** %zu **\n", pr->dim); */
  t1p_t *res;
  size_t intdim = a1->intdim;
  size_t realdim = a1->dims - a1->intdim;
  if (t1p_is_eq (man, a1, a2)) {
    if (destructive)
      res = a1;
    else
      res = t1p_copy (man, a1);
  }
  else if (tbool_or (t1p_is_top (man, a1), t1p_is_top (man, a2)) == tbool_true) {
    if (destructive) {
      t1p_free (man, a1);
      res = a1 = t1p_top (man, intdim, realdim);
    }
    else {
      res = t1p_top (man, intdim, realdim);
    }
  }
  else if (t1p_is_bottom (man, a1)) {
    if (destructive) {
      t1p_free (man, a1);
      res = a1 = t1p_copy (man, a2);
    }
    else {
      res = t1p_copy (man, a2);
    }
  }
  else if (t1p_is_bottom (man, a2)) {
    if (destructive)
      res = a1;
    else
      res = t1p_copy (man, a1);
  }
  else {
    /* TODO: destructive not yet supported */
    itv_t tmp;
    itv_init (tmp);
    res = t1p_alloc (man, intdim, realdim);
    /* update res->box */
    for (i = 0; i < (intdim + realdim); i++)
      itv_join (res->box[i], a1->box[i], a2->box[i]);

    if (a1->hypercube && a2->hypercube) {
      for (i = 0; i < (intdim + realdim); i++) {
	//printf("%d: ",i);
	if (t1p_aff_is_bottom (pr, a1->paf[i]))
	  res->paf[i] = a2->paf[i];
	else if (t1p_aff_is_bottom (pr, a2->paf[i]))
	  res->paf[i] = a1->paf[i];
	else if (t1p_aff_is_top (pr, a1->paf[i]) || t1p_aff_is_top (pr, a2->paf[i]))
	  res->paf[i] = pr->top;
	/*
	   else if (t1p_aff_is_eq(pr, a1->paf[i], a2->paf[i])) {
	   fprintf(stream,"** res = x **\n");
	   res->paf[i] = a1->paf[i];
	   } else if (t1p_aff_is_leq_constrained(pr, a2->paf[i], a1->paf[i], a2,a1)) {
	   fprintf(stream,"** res = x **\n");
	   res->paf[i] = a1->paf[i];
	   } else if (t1p_aff_is_leq_constrained(pr, a1->paf[i], a2->paf[i], a1, a2)) {
	   fprintf(stream,"** res = y **\n");
	   res->paf[i] = a2->paf[i];
	   *
	   res->paf[i] = t1p_aff_alloc_init(pr);
	   itv_set(res->paf[i]->c, a2->paf[i]->c);
	   itv_set(res->paf[i]->itv, a2->paf[i]->itv);
	   t1p_aaterm_t* p = NULL;
	   for(p = a2->paf[i]->q; p; p=p->n) {
	   t1p_aff_nsym_add(pr, res->paf[i], p->coeff, p->pnsym);
	   if (p->pnsym->type == UN) {
	   res->paf[i]->end->pnsym = pr->mubGlobal.p[p->pnsym->index].y;
	   }
	   }
	   *
	   } */
	else {
	  if (itv_has_infty_bound (a1->box[i]) || itv_has_infty_bound (a2->box[i])) {
	    /* Do nothing, the join of concretisations is already done and stored in res->box */
	    res->paf[i] = t1p_aff_alloc_init (pr);
	    itv_set (res->paf[i]->c, res->box[i]);
	  }
	  else {
	    /* join two affine form expressions */
	    itv_set (a1->paf[i]->itv, a1->box[i]);
	    itv_set (a2->paf[i]->itv, a2->box[i]);
	    res->paf[i] =
	      t1p_aff_join_constrained6 (pr, a1->paf[i], a2->paf[i], a1, a2, res);
	  }
	}
	//printf("%d",i);itv_print(a1->box[i]); printf("\t");itv_print(a2->box[i]);printf("\n");
	res->paf[i]->pby++;
      }

    }
    else {
      size_t k = 0;
      ap_dim_t j = 0;
      size_t dims1 = t1p_nsymcons_get_dimension (pr, a1);
      size_t dims2 = t1p_nsymcons_get_dimension (pr, a2);
      if (dims1 && dims2) {
	size_t dim2 = 0;
	/* au pire il y a toute la liste de l'un a rajouter dans l'autre */
	ap_dimchange_t *dimchange2 = ap_dimchange_alloc (0, dims1);
	/* par defaut, res->size = 128 lors de sa creation. verifier que la taille est ok avant de faire la copie */
	if (dims1 > res->size) {
	  res->nsymcons =
	    (ap_dim_t *) realloc (res->nsymcons, (dims1) * sizeof (ap_dim_t));
	  res->gamma =
	    (ap_interval_t **) realloc (res->gamma,
					(dims1) * sizeof (ap_interval_t *));
	  for (k = res->size; k < dims1; k++)
	    res->gamma[k] = NULL;
	  res->size = dims1;
	}
	res->nsymcons =
	  memcpy ((void *) res->nsymcons, (const void *) a1->nsymcons,
		  dims1 * sizeof (ap_dim_t));
	//for (i=0; i<dims1; i++) res->gamma[i] = ap_interval_alloc_set(a1->gamma[i]);
	//for (i=0; i<dims1; i++) res->gamma[i] = NULL;
	ap_abstract0_free (pr->manNS, res->abs);

	res->abs = ap_abstract0_copy (pr->manNS, a1->abs);
	for (k = 0; k < dims1; k++) {
	  if (!t1p_nsymcons_get_dimpos (pr, &j, a1->nsymcons[k], a2)) {
	    dimchange2->dim[dim2] = j;
	    dim2++;
	  }
	}
	dimchange2->realdim = dim2;
	size_t dim1 = 0;
	for (k = 0; k < dims2; k++)
	  t1p_insert_constrained_nsym (pr, &j, a2->nsymcons[k], res);

	/* destructive, without projection (new dimension set to top) */
	ap_abstract0_add_dimensions (pr->manNS, true, a2->abs, dimchange2, false);
	/* ajouter les contraintes <= 1, >= -1 aux nouvelles dimensions ajoutees a a2 */
	ap_dimchange_add_invert (dimchange2);
	for (k = 0; k < dim2; k++) {
	  t1p_set_lincons_dim (pr, dimchange2->dim[k]);
	  ap_abstract0_meet_lincons_array (pr->manNS, true, a2->abs, &pr->moo);
	}

	ap_abstract0_join (pr->manNS, true, res->abs, a2->abs);
	/*if (ap_abstract0_is_leq(pr->manNS, res->abs, a2->abs)) {
	   for (k=0; k<dim2; k++) t1p_delete_constrained_nsym(pr, dimchange2->dim[k], res) ...;
	   } */

	/* update res->gamma */
	t1p_update_nsymcons_gamma (pr, res);

	ap_abstract0_remove_dimensions (pr->manNS, true, a2->abs, dimchange2);

	dimchange2->realdim = dims2;
	ap_dimchange_free (dimchange2);

	size_t nsymcons_size = t1p_nsymcons_get_dimension (pr, res);
	pr->dimtoremove =
	  (ap_dim_t *) realloc (pr->dimtoremove,
				(nsymcons_size) * sizeof (ap_dim_t));
	memset ((void *) pr->dimtoremove, (int) 0, nsymcons_size * sizeof (int));
      }
      else {
	/* res->abs is a hypercube */
      }
      for (i = 0; i < (intdim + realdim); i++) {
	if (t1p_aff_is_bottom (pr, a1->paf[i]))
	  res->paf[i] = a2->paf[i];
	else if (t1p_aff_is_bottom (pr, a2->paf[i]))
	  res->paf[i] = a1->paf[i];
	else if (t1p_aff_is_top (pr, a1->paf[i]) || t1p_aff_is_top (pr, a2->paf[i]))
	  res->paf[i] = pr->top;
	else if (t1p_aff_is_eq (pr, a1->paf[i], a2->paf[i]))
	  res->paf[i] = a1->paf[i];
	else {
	  if (itv_has_infty_bound (a1->box[i]) || itv_has_infty_bound (a2->box[i])) {
	    /* Do nothing, the join of concretisations is already done and stored in res->box */
	    res->paf[i] = t1p_aff_alloc_init (pr);
	    itv_set (res->paf[i]->c, res->box[i]);
	  }
	  else {
	    /* join two affine form expressions */
	    itv_set (a1->paf[i]->itv, a1->box[i]);
	    itv_set (a2->paf[i]->itv, a2->box[i]);
	    res->paf[i] =
	      t1p_aff_join_constrained6 (pr, a1->paf[i], a2->paf[i], a1, a2, res);
	  }
	}
	res->paf[i]->pby++;
      }

      man->result.flag_best = tbool_top;
      man->result.flag_exact = tbool_top;
    }
    pr->mubGlobal.cx = NULL;
    pr->mubGlobal.cy = NULL;
    free (pr->mubGlobal.p);
    man->result.flag_best = tbool_true;
    man->result.flag_exact = tbool_top;
    itv_clear (tmp);
  }
  man->result.flag_best = tbool_true;
  man->result.flag_exact = tbool_true;


#ifdef _T1P_DEBUG
  fprintf (stdout, "### RESULT of JOIN (des %d) [%tx] ###\n", destructive,
	   (intptr_t) res);
  t1p_fprint (stdout, man, res, NULL);
  /* for (i=0; i<intdim+realdim;i++) */
  /*   {fprintf(stream,"%zu :",i);itv_fprint(stream,res->box[i]); fprintf(stream,"\n");} */
  /* fprintf(stream,"** %zu **\n", pr->dim); */
  /* fprintf(stream, "### ### ###\n"); */
  /* fclose(stream); */
#endif

  return res;
}

/* pour essayer avec constrained8/8bis  */
t1p_t *
t1p_join_faux (ap_manager_t * man, bool destructive, t1p_t * a1, t1p_t * a2)
    /* TODO destructive not used  */
{
  CALL ();
  size_t i = 0;
  t1p_internal_t *pr = (t1p_internal_t *) t1p_init_from_manager (man, AP_FUNID_JOIN);
  arg_assert (a1->dims == a2->dims && a1->intdim == a2->intdim, abort ();
    );
#ifdef _T1P_DEBUG
  fprintf (stdout, "### JOIN OPERANDS (des %d) (%tx U %tx) ###\n", destructive,
	   (intptr_t) a1, (intptr_t) a2);
  //t1p_fprint(stdout, man, a1, NULL);
  //t1p_fprint(stdout, man, a2, NULL);
  fprintf (stdout, "### ### ###\n");
#endif
  /* FILE* stream = fopen("/home/donquijote/taylor1p/taylor1plus/taylor1plus/demo/log", "a+"); */
  /* fprintf(stream,"** %zu **\n", pr->dim); */
  t1p_t *res;
  size_t intdim = a1->intdim;
  size_t realdim = a1->dims - a1->intdim;
  if (t1p_is_eq (man, a1, a2)) {
    if (destructive)
      res = a1;
    else
      res = t1p_copy (man, a1);
  }
  else if (tbool_or (t1p_is_top (man, a1), t1p_is_top (man, a2)) == tbool_true) {
    if (destructive) {
      t1p_free (man, a1);
      res = a1 = t1p_top (man, intdim, realdim);
    }
    else {
      res = t1p_top (man, intdim, realdim);
    }
  }
  else if (t1p_is_bottom (man, a1)) {
    if (destructive) {
      t1p_free (man, a1);
      res = a1 = t1p_copy (man, a2);
    }
    else {
      res = t1p_copy (man, a2);
    }
  }
  else if (t1p_is_bottom (man, a2)) {
    if (destructive)
      res = a1;
    else
      res = t1p_copy (man, a1);
  }
  else {
    /* TODO: destructive not yet supported */
    itv_t tmp;
    itv_init (tmp);
    res = t1p_alloc (man, intdim, realdim);
    /* update res->box */
    for (i = 0; i < (intdim + realdim); i++)
      itv_join (res->box[i], a1->box[i], a2->box[i]);

    size_t old = pr->dim;
    for (i = 0; i < pr->epssize; i++) {
      if (i != pr->inputns[i]) {
	t1p_nsym_t *tmp;
	tmp = pr->epsilon[i];
	pr->epsilon[i] = pr->epsilon[pr->inputns[i]];
	pr->epsilon[pr->inputns[i]] = tmp;
	pr->epsilon[i]->index = i;
	pr->epsilon[pr->inputns[i]]->index = pr->inputns[i];
	pr->inputns[i] = i;
      }
    }
    pr->mubGlobal.p = (Tobj1 *) calloc (old, sizeof (Tobj1));
    for (i = pr->epssize; i < old; i++) {
      pr->mubGlobal.p[i].x = pr->epsilon[i];
      pr->mubGlobal.p[i].y = t1p_nsym_add (pr, UN);
    }
    pr->mubGlobal.cx = t1p_nsym_add (pr, UN);
    for (i = 0; i < pr->epssize; i++) {
      pr->mubGlobal.p[i].x = t1p_nsym_add (pr, UN);
    }
    pr->mubGlobal.cy = t1p_nsym_add (pr, UN);
    for (i = 0; i < pr->epssize; i++) {
      pr->mubGlobal.p[i].y = t1p_nsym_add (pr, UN);
    }
    if (a1->hypercube && a2->hypercube) {
      for (i = 0; i < (intdim + realdim); i++) {
	t1p_aff_canonical (pr, a1->paf[i]);
	t1p_aff_canonical (pr, a2->paf[i]);
	//printf("%d: ",i);
	if (t1p_aff_is_bottom (pr, a1->paf[i]))
	  res->paf[i] = a2->paf[i];
	else if (t1p_aff_is_bottom (pr, a2->paf[i]))
	  res->paf[i] = a1->paf[i];
	else if (t1p_aff_is_top (pr, a1->paf[i]) || t1p_aff_is_top (pr, a2->paf[i]))
	  res->paf[i] = pr->top;
	/*
	   else if (t1p_aff_is_eq(pr, a1->paf[i], a2->paf[i])) {
	   fprintf(stream,"** res = x **\n");
	   res->paf[i] = a1->paf[i];
	   } else if (t1p_aff_is_leq_constrained(pr, a2->paf[i], a1->paf[i], a2,a1)) {
	   fprintf(stream,"** res = x **\n");
	   res->paf[i] = a1->paf[i];
	   } else if (t1p_aff_is_leq_constrained(pr, a1->paf[i], a2->paf[i], a1, a2)) {
	   fprintf(stream,"** res = y **\n");
	   res->paf[i] = a2->paf[i];
	   *
	   res->paf[i] = t1p_aff_alloc_init(pr);
	   itv_set(res->paf[i]->c, a2->paf[i]->c);
	   itv_set(res->paf[i]->itv, a2->paf[i]->itv);
	   t1p_aaterm_t* p = NULL;
	   for(p = a2->paf[i]->q; p; p=p->n) {
	   t1p_aff_nsym_add(pr, res->paf[i], p->coeff, p->pnsym);
	   if (p->pnsym->type == UN) {
	   res->paf[i]->end->pnsym = pr->mubGlobal.p[p->pnsym->index].y;
	   }
	   }
	   *
	   } */
	else {
	  if (itv_has_infty_bound (a1->box[i]) || itv_has_infty_bound (a2->box[i])) {
	    /* Do nothing, the join of concretisations is already done and stored in res->box */
	    res->paf[i] = t1p_aff_alloc_init (pr);
	    itv_set (res->paf[i]->c, res->box[i]);
	  }
	  else {
	    /* join two affine form expressions */
	    itv_set (a1->paf[i]->itv, a1->box[i]);
	    itv_set (a2->paf[i]->itv, a2->box[i]);
	    res->paf[i] =
	      t1p_aff_join_constrained8 (pr, a1->paf[i], a2->paf[i], a1, a2, res);
	  }
	}
	//printf("%d",i);itv_print(a1->box[i]); printf("\t");itv_print(a2->box[i]);printf("\n");
	res->paf[i]->pby++;
      }
      size_t k = 0;
      res->g = (itv_t **) calloc (1 + pr->dim, sizeof (itv_t *));	/* centre + au max chaque symbole lui est associ\E9 un g\E9n\E9rateur */
      res->gn = pr->dim;

      t1p_aaterm_t *p = NULL;
      itv_t tmp1;
      itv_init (tmp1);

      /* init res->g (structure creuse) */
      for (i = 0; i < intdim + realdim; i++) {
	if (!res->g[0])
	  res->g[0] = (itv_t *) calloc (intdim + realdim, sizeof (itv_t));
	itv_init (res->g[0][i]);
	itv_set (res->g[0][i], res->paf[i]->c);
	for (p = res->paf[i]->q; p; p = p->n) {
	  if (!res->g[1 + p->pnsym->index])
	    res->g[1 + p->pnsym->index] =
	      (itv_t *) calloc (intdim + realdim, sizeof (itv_t));
	  itv_init (res->g[1 + p->pnsym->index][i]);
	  itv_set (res->g[1 + p->pnsym->index][i], p->coeff);
	}
      }
      /* virer les vecteurs colineaires */
      size_t j = 0;
      itv_t un;
      itv_init (un);
      itv_set_int (un, 1);
      for (i = 1 + pr->epssize; i < 1 + (-1 + pr->dim); i++) {
	for (j = i + 1; j < 1 + pr->dim; j++) {
	  if (res->g[i] && res->g[j]) {
	    if (res->g[j][0] && res->g[i][0])
	      itv_div (pr->itv, tmp, res->g[j][0], res->g[i][0]);
	    else
	      break;
	    for (k = 1; k < intdim + realdim;) {
	      if (res->g[j][k] && res->g[i][k])
		itv_div (pr->itv, tmp1, res->g[j][k], res->g[i][k]);
	      else
		break;
	      if (itv_is_eq (tmp, tmp1))
		k++;
	      else
		break;
	    }
	    if (k == intdim + realdim) {
	      itv_abs (tmp1, tmp1);
	      itv_add (tmp1, tmp1, un);
	      for (k = 0; k < intdim + realdim; k++) {
		itv_mul (pr->itv, res->g[i][k], tmp1, res->g[i][k]);
		itv_clear (res->g[j][k]);
	      }
	      free (res->g[j]);
	      res->g[j] = NULL;
	    }
	  }
	}
      }
      itv_clear (un);
      itv_clear (tmp1);

      /* FILE* poly = fopen("/home/donquijote/taylor1p/taylor1plus/taylor1plus/demo/draw.poly", "a"); */
      /*  fprintf(poly,"$z%zu=zonotope(new Matrix<Rational>([",pr->it); */
      /* for (i=0;i<pr->dim;i++) { */
      /*     if (res->g[i+1]) { */
      /*     fprintf(poly,"[0,"); */
      /*      for (k=0;k<3;k++) { */
      /*          if (res->g[i+1][k] == NULL) {fprintf(poly,"0");} */
      /*          else if (itv_is_zero(res->g[i+1][k])) {fprintf(poly,"0");} */
      /*          else {bound_fprint(poly,res->g[i+1][k]->sup);} */
      /*          if (k<2) fprintf(poly,","); */
      /*      } */
      /*     fprintf(poly,"],"); */
      /*     } */
      /* } */
      /* fprintf(poly,"]));\n"); */
      /*  fprintf(poly,"$p%zu=new Polytope<Rational>(POINTS=>$z%zu);\n",pr->it,pr->it); */
      /*  fprintf(poly,"$pp%zu=transform($p%zu,new Matrix<Rational>([[1,",pr->it,pr->it); */
      /*      for (k=0;k<3;k++) { */
      /*          if (res->g[0][k] == NULL) {fprintf(poly,"0");} */
      /*          else if (itv_is_zero(res->g[0][k])) {fprintf(poly,"0");} */
      /*          else {bound_fprint(poly,res->g[0][k]->sup);} */
      /*          if (k<2) fprintf(poly,","); */
      /*      } */
      /*  fprintf(poly,"],[0,1,0,0],[0,0,1,0],[0,0,0,1]]));\n"); */
      /*  fprintf(poly,"javaview($pp%zu->VISUAL,File=>\"pp%zu\");\n",pr->it,pr->it); */

      /* fclose(poly); */
      pr->it++;

    }
    else {
      size_t k = 0;
      ap_dim_t j = 0;
      size_t dims1 = t1p_nsymcons_get_dimension (pr, a1);
      size_t dims2 = t1p_nsymcons_get_dimension (pr, a2);
      if (dims1 && dims2) {
	size_t dim2 = 0;
	/* au pire il y a toute la liste de l'un a rajouter dans l'autre */
	ap_dimchange_t *dimchange2 = ap_dimchange_alloc (0, dims2);
	if (dims1 > res->size) {
	  res->nsymcons =
	    (ap_dim_t *) realloc (res->nsymcons, (dims1) * sizeof (ap_dim_t));
	  res->gamma =
	    (ap_interval_t **) realloc (res->gamma,
					(dims1) * sizeof (ap_interval_t *));
	  for (k = res->size; k < dims1; k++)
	    res->gamma[k] = NULL;
	  res->size = dims1;
	}
	res->nsymcons =
	  memcpy ((void *) res->nsymcons, (const void *) a1->nsymcons,
		  dims1 * sizeof (ap_dim_t));
	for (i = 0; i < dims1; i++)
	  res->gamma[i] = ap_interval_alloc_set (a1->gamma[i]);
	ap_abstract0_free (pr->manNS, res->abs);
	res->abs = ap_abstract0_copy (pr->manNS, a1->abs);
	for (k = 0; k < dims1; k++) {
	  if (!t1p_nsymcons_get_dimpos (pr, &j, a1->nsymcons[k], a2)) {
	    dimchange2->dim[dim2] = j;
	    dim2++;
	  }
	}
	dimchange2->realdim = dim2;
	size_t dim1 = 0;
	for (k = 0; k < dims2; k++)
	  t1p_insert_constrained_nsym (pr, &j, a2->nsymcons[k], res);

	/* destructive, without projection (new dimension set to top) */
	ap_abstract0_add_dimensions (pr->manNS, true, a2->abs, dimchange2, false);

	ap_abstract0_join (pr->manNS, true, res->abs, a2->abs);

	/* update res->gamma */
	t1p_update_nsymcons_gamma (pr, res);

	ap_dimchange_add_invert (dimchange2);
	ap_abstract0_remove_dimensions (pr->manNS, true, a2->abs, dimchange2);

	dimchange2->realdim = dims2;
	ap_dimchange_free (dimchange2);

	size_t nsymcons_size = t1p_nsymcons_get_dimension (pr, res);
	pr->dimtoremove =
	  (ap_dim_t *) realloc (pr->dimtoremove,
				(nsymcons_size) * sizeof (ap_dim_t));
	memset ((void *) pr->dimtoremove, (int) 0, nsymcons_size * sizeof (int));
      }
      else {
	/* res->abs is a hypercube */
      }
      for (i = 0; i < (intdim + realdim); i++) {
	if (t1p_aff_is_bottom (pr, a1->paf[i]))
	  res->paf[i] = a2->paf[i];
	else if (t1p_aff_is_bottom (pr, a2->paf[i]))
	  res->paf[i] = a1->paf[i];
	else if (t1p_aff_is_top (pr, a1->paf[i]) || t1p_aff_is_top (pr, a2->paf[i]))
	  res->paf[i] = pr->top;
	else if (t1p_aff_is_eq (pr, a1->paf[i], a2->paf[i]))
	  res->paf[i] = a1->paf[i];
	else {
	  if (itv_has_infty_bound (a1->box[i]) || itv_has_infty_bound (a2->box[i])) {
	    /* Do nothing, the join of concretisations is already done and stored in res->box */
	    res->paf[i] = t1p_aff_alloc_init (pr);
	    itv_set (res->paf[i]->c, res->box[i]);
	  }
	  else {
	    /* join two affine form expressions */
	    itv_set (a1->paf[i]->itv, a1->box[i]);
	    itv_set (a2->paf[i]->itv, a2->box[i]);
	    res->paf[i] =
	      t1p_aff_join_constrained8 (pr, a1->paf[i], a2->paf[i], a1, a2, res);
	  }
	}
	res->paf[i]->pby++;
      }

      man->result.flag_best = tbool_top;
      man->result.flag_exact = tbool_top;
    }
    pr->mubGlobal.cx = NULL;
    pr->mubGlobal.cy = NULL;
    free (pr->mubGlobal.p);
    man->result.flag_best = tbool_true;
    man->result.flag_exact = tbool_top;
    itv_clear (tmp);
  }
  man->result.flag_best = tbool_true;
  man->result.flag_exact = tbool_true;


#ifdef _T1P_DEBUG
  fprintf (stdout, "### RESULT of JOIN (des %d) [%tx] ###\n", destructive,
	   (intptr_t) res);
  //t1p_fprint(stdout, man, res, NULL);
  /* for (i=0; i<intdim+realdim;i++) */
  /*   {fprintf(stream,"%zu :",i);itv_fprint(stream,res->box[i]); fprintf(stream,"\n");} */
  /* fprintf(stream,"** %zu **\n", pr->dim); */
  /* fprintf(stream, "### ### ###\n"); */
  /* fclose(stream); */
#endif

  return res;
}

/* ub global (fausse) */
t1p_t *
t1p_join_bub (ap_manager_t * man, bool destructive, t1p_t * a1, t1p_t * a2)
    /* TODO destructive not used  */
{
  CALL ();
  size_t i = 0;
  t1p_internal_t *pr = (t1p_internal_t *) t1p_init_from_manager (man, AP_FUNID_JOIN);
  arg_assert (a1->dims == a2->dims && a1->intdim == a2->intdim, abort ();
    );
#ifdef _T1P_DEBUG
  fprintf (stdout, "### JOIN OPERANDS (des %d) (%tx U %tx) ###\n", destructive,
	   (intptr_t) a1, (intptr_t) a2);
  //t1p_fprint(stdout, man, a1, NULL);
  //t1p_fprint(stdout, man, a2, NULL);
  fprintf (stdout, "### ### ###\n");
  //fprintf(stream,"** %d **\n", pr->dim);
#endif
  /* FILE* stream = fopen("/home/donquijote/taylor1p/taylor1plus/taylor1plus/demo/log", "a+"); */
  t1p_t *res;
  size_t intdim = a1->intdim;
  size_t realdim = a1->dims - a1->intdim;
  if (t1p_is_eq (man, a1, a2)) {
    if (destructive)
      res = a1;
    else
      res = t1p_copy (man, a1);
  }
  else if (tbool_or (t1p_is_top (man, a1), t1p_is_top (man, a2)) == tbool_true) {
    if (destructive) {
      t1p_free (man, a1);
      res = a1 = t1p_top (man, intdim, realdim);
    }
    else {
      res = t1p_top (man, intdim, realdim);
    }
  }
  else if (t1p_is_bottom (man, a1)) {
    if (destructive) {
      t1p_free (man, a1);
      res = a1 = t1p_copy (man, a2);
    }
    else {
      res = t1p_copy (man, a2);
    }
  }
  else if (t1p_is_bottom (man, a2)) {
    if (destructive)
      res = a1;
    else
      res = t1p_copy (man, a1);
  }
  else {
    /* TODO: destructive not yet supported */
    itv_t tmp;
    itv_init (tmp);
    res = t1p_alloc (man, intdim, realdim);
    /* update res->box */
    for (i = 0; i < (intdim + realdim); i++)
      itv_join (res->box[i], a1->box[i], a2->box[i]);

    size_t old = pr->dim;
    for (i = 0; i < pr->epssize; i++) {
      if (i != pr->inputns[i]) {
	t1p_nsym_t *tmp;
	tmp = pr->epsilon[i];
	pr->epsilon[i] = pr->epsilon[pr->inputns[i]];
	pr->epsilon[pr->inputns[i]] = tmp;
	pr->epsilon[i]->index = i;
	pr->epsilon[pr->inputns[i]]->index = pr->inputns[i];
	pr->inputns[i] = i;
      }
    }
    pr->mubGlobal.p = (Tobj1 *) calloc (old, sizeof (Tobj1));
    /*
       for (i=pr->epssize;i<old;i++) {
       pr->mubGlobal.p[i].x = pr->epsilon[i];
       pr->mubGlobal.p[i].y = t1p_nsym_add(pr, UN);
       }
     */
    pr->mubGlobal.cx = t1p_nsym_add (pr, UN);
    for (i = 0; i < pr->epssize; i++) {
      pr->mubGlobal.p[i].x = t1p_nsym_add (pr, UN);
    }
    /*
       pr->mubGlobal.cy = t1p_nsym_add(pr, UN);
       for (i=0;i<pr->epssize;i++) {
       pr->mubGlobal.p[i].y = t1p_nsym_add(pr, UN);
       }
     */
    if (a1->hypercube && a2->hypercube) {
      for (i = 0; i < (intdim + realdim); i++) {
	t1p_aff_canonical (pr, a1->paf[i]);
	t1p_aff_canonical (pr, a2->paf[i]);
	//printf("%d: ",i);
	if (t1p_aff_is_bottom (pr, a1->paf[i]))
	  res->paf[i] = a2->paf[i];
	else if (t1p_aff_is_bottom (pr, a2->paf[i]))
	  res->paf[i] = a1->paf[i];
	else if (t1p_aff_is_top (pr, a1->paf[i]) || t1p_aff_is_top (pr, a2->paf[i]))
	  res->paf[i] = pr->top;
	/*
	   else if (t1p_aff_is_eq(pr, a1->paf[i], a2->paf[i])) {
	   res->paf[i] = a1->paf[i];
	   } else if (t1p_aff_is_leq_constrained(pr, a2->paf[i], a1->paf[i], a2,a1)) {
	   res->paf[i] = a1->paf[i];
	   } else if (t1p_aff_is_leq_constrained(pr, a1->paf[i], a2->paf[i], a1, a2)) {
	   res->paf[i] = a2->paf[i];
	   } */
	else {
	  if (itv_has_infty_bound (a1->box[i]) || itv_has_infty_bound (a2->box[i])) {
	    /* Do nothing, the join of concretisations is already done and stored in res->box */
	    res->paf[i] = t1p_aff_alloc_init (pr);
	    itv_set (res->paf[i]->c, res->box[i]);
	  }
	  else {
	    /* join two affine form expressions */
	    itv_set (a1->paf[i]->itv, a1->box[i]);
	    itv_set (a2->paf[i]->itv, a2->box[i]);
	    //res->paf[i] = t1p_aff_join_arXiv2(pr, a1->paf[i], a2->paf[i], a1, a2, res);
	    //res->paf[i] = t1p_aff_join_arXiv2bis(pr, a1->paf[i], a2->paf[i], a1, a2, res);
	    res->paf[i] =
	      t1p_aff_join_arXiv2ter (pr, a1->paf[i], a2->paf[i], a1, a2, res);
	  }
	}
	//fprintf(stream,"********************************************\n");
	//itv_fprint(stream,res->paf[i]->end->coeff);fprintf(stream,"\n");
	//fprintf(stream,"********************************************\n");
	res->paf[i]->pby++;
      }
    }
    else {
      size_t k = 0;
      ap_dim_t j = 0;
      size_t dims1 = t1p_nsymcons_get_dimension (pr, a1);
      size_t dims2 = t1p_nsymcons_get_dimension (pr, a2);
      if (dims1 && dims2) {
	size_t dim2 = 0;
	/* au pire il y a toute la liste de l'un a rajouter dans l'autre */
	ap_dimchange_t *dimchange2 = ap_dimchange_alloc (0, dims2);
	if (dims1 > res->size) {
	  res->nsymcons =
	    (ap_dim_t *) realloc (res->nsymcons, (dims1) * sizeof (ap_dim_t));
	  res->gamma =
	    (ap_interval_t **) realloc (res->gamma,
					(dims1) * sizeof (ap_interval_t *));
	  for (k = res->size; k < dims1; k++)
	    res->gamma[k] = NULL;
	  res->size = dims1;
	}
	res->nsymcons =
	  memcpy ((void *) res->nsymcons, (const void *) a1->nsymcons,
		  dims1 * sizeof (ap_dim_t));
	for (i = 0; i < dims1; i++)
	  res->gamma[i] = ap_interval_alloc_set (a1->gamma[i]);
	ap_abstract0_free (pr->manNS, res->abs);
	res->abs = ap_abstract0_copy (pr->manNS, a1->abs);
	for (k = 0; k < dims1; k++) {
	  if (!t1p_nsymcons_get_dimpos (pr, &j, a1->nsymcons[k], a2)) {
	    dimchange2->dim[dim2] = j;
	    dim2++;
	  }
	}
	dimchange2->realdim = dim2;
	size_t dim1 = 0;
	for (k = 0; k < dims2; k++)
	  t1p_insert_constrained_nsym (pr, &j, a2->nsymcons[k], res);

	/* destructive, without projection (new dimension set to top) */
	ap_abstract0_add_dimensions (pr->manNS, true, a2->abs, dimchange2, false);

	ap_abstract0_join (pr->manNS, true, res->abs, a2->abs);

	/* update res->gamma */
	t1p_update_nsymcons_gamma (pr, res);

	ap_dimchange_add_invert (dimchange2);
	ap_abstract0_remove_dimensions (pr->manNS, true, a2->abs, dimchange2);

	dimchange2->realdim = dims2;
	ap_dimchange_free (dimchange2);

	size_t nsymcons_size = t1p_nsymcons_get_dimension (pr, res);
	pr->dimtoremove =
	  (ap_dim_t *) realloc (pr->dimtoremove,
				(nsymcons_size) * sizeof (ap_dim_t));
	memset ((void *) pr->dimtoremove, (int) 0, nsymcons_size * sizeof (int));
      }
      else {
	/* res->abs is a hypercube */
      }
      for (i = 0; i < (intdim + realdim); i++) {
	if (t1p_aff_is_bottom (pr, a1->paf[i]))
	  res->paf[i] = a2->paf[i];
	else if (t1p_aff_is_bottom (pr, a2->paf[i]))
	  res->paf[i] = a1->paf[i];
	else if (t1p_aff_is_top (pr, a1->paf[i]) || t1p_aff_is_top (pr, a2->paf[i]))
	  res->paf[i] = pr->top;
	else if (t1p_aff_is_eq (pr, a1->paf[i], a2->paf[i]))
	  res->paf[i] = a1->paf[i];
	else {
	  if (itv_has_infty_bound (a1->box[i]) || itv_has_infty_bound (a2->box[i])) {
	    /* Do nothing, the join of concretisations is already done and stored in res->box */
	    res->paf[i] = t1p_aff_alloc_init (pr);
	    itv_set (res->paf[i]->c, res->box[i]);
	  }
	  else {
	    /* join two affine form expressions */
	    itv_set (a1->paf[i]->itv, a1->box[i]);
	    itv_set (a2->paf[i]->itv, a2->box[i]);
	    res->paf[i] =
	      t1p_aff_join_constrained8 (pr, a1->paf[i], a2->paf[i], a1, a2, res);
	  }
	}
	res->paf[i]->pby++;
      }

      man->result.flag_best = tbool_top;
      man->result.flag_exact = tbool_top;
    }
    man->result.flag_best = tbool_true;
    man->result.flag_exact = tbool_top;
    itv_clear (tmp);
  }
  man->result.flag_best = tbool_true;
  man->result.flag_exact = tbool_true;


#ifdef _T1P_DEBUG
  fprintf (stdout, "### RESULT of JOIN (des %d) [%tx] ###\n", destructive,
	   (intptr_t) res);
  //t1p_fprint(stdout, man, res, NULL);
  /*
     for (i=0; i<intdim+realdim;i++)
     {fprintf(stream,"%d :",i);itv_fprint(stream,res->box[i]); fprintf(stream,"\n");}
     fprintf(stream,"** %d **\n", pr->dim);
     fprintf(stream, "### ### ###\n");
   */
  //fprintf(stream,"-");num_fprint(stream,bound_numref(res->box[0]->inf));fprintf(stream,"\t");num_fprint(stream,bound_numref(res->box[0]->sup));fprintf(stream,"\t%zu\n",pr->dim);
  //t1p_fprint(stream,man,res,0);
  //fclose(stream);
#endif

  return res;
}

/* ub avec formule de l'identit\E9 du //  */
/* pour tester le global + drawzonotopes */
t1p_t *
t1p_join_global (ap_manager_t * man, bool destructive, t1p_t * a1, t1p_t * a2)
    /* TODO destructive not used  */
{
  CALL ();
  size_t i = 0;
  t1p_internal_t *pr = (t1p_internal_t *) t1p_init_from_manager (man, AP_FUNID_JOIN);
  arg_assert (a1->dims == a2->dims && a1->intdim == a2->intdim, abort ();
    );
#ifdef _T1P_DEBUG
  fprintf (stdout, "### JOIN OPERANDS (des %d) (%tx U %tx) ###\n", destructive,
	   (intptr_t) a1, (intptr_t) a2);
  //t1p_fprint(stdout, man, a1, NULL);
  //t1p_fprint(stdout, man, a2, NULL);
  fprintf (stdout, "### ### ###\n");
  //fprintf(stream,"** %d **\n", pr->dim);
#endif
  /* FILE* stream = fopen("/home/donquijote/taylor1p/taylor1plus/taylor1plus/demo/log", "a+"); */
  t1p_t *res;
  size_t intdim = a1->intdim;
  size_t realdim = a1->dims - a1->intdim;
  if (t1p_is_eq (man, a1, a2)) {
    if (destructive)
      res = a1;
    else
      res = t1p_copy (man, a1);
  }
  else if (tbool_or (t1p_is_top (man, a1), t1p_is_top (man, a2)) == tbool_true) {
    if (destructive) {
      t1p_free (man, a1);
      res = a1 = t1p_top (man, intdim, realdim);
    }
    else {
      res = t1p_top (man, intdim, realdim);
    }
  }
  else if (t1p_is_bottom (man, a1)) {
    if (destructive) {
      t1p_free (man, a1);
      res = a1 = t1p_copy (man, a2);
    }
    else {
      res = t1p_copy (man, a2);
    }
  }
  else if (t1p_is_bottom (man, a2)) {
    if (destructive)
      res = a1;
    else
      res = t1p_copy (man, a1);
  }
  else {
    /* TODO: destructive not yet supported */
    itv_t tmp;
    itv_init (tmp);
    itv_t tmp1;
    itv_init (tmp1);
    res = t1p_alloc (man, intdim, realdim);
    /* update res->box */
    for (i = 0; i < (intdim + realdim); i++)
      itv_join (res->box[i], a1->box[i], a2->box[i]);

    /*
       a1->g = (itv_t**)calloc(1+pr->dim,sizeof(itv_t*));
       a2->g = (itv_t**)calloc(1+pr->dim,sizeof(itv_t*));
       res->g = (itv_t**)calloc(1+pr->dim,sizeof(itv_t*));
       a1->gn = a2->gn = res->gn = pr->dim;

       t1p_aaterm_t *p = NULL;
       for(i=0;i<intdim+realdim;i++) {
       if (!a1->g[0]) a1->g[0] = (itv_t*)calloc(intdim+realdim,sizeof(itv_t));
       itv_init(a1->g[0][i]);
       itv_set(a1->g[0][i],a1->paf[i]->c);
       for(p=a1->paf[i]->q;p;p=p->n) {
       if (!a1->g[p->pnsym->index]) a1->g[p->pnsym->index] = (itv_t*)calloc(intdim+realdim,sizeof(itv_t));
       itv_init(a1->g[p->pnsym->index][i]);
       itv_set(a1->g[p->pnsym->index][i],p->coeff);
       }

       if (!a2->g[0]) a2->g[0] = (itv_t*)calloc(intdim+realdim,sizeof(itv_t));
       itv_init(a2->g[0][i]);
       itv_set(a2->g[0][i],a2->paf[i]->c);
       for(p=a2->paf[i]->q;p;p=p->n) {
       if (!a2->g[p->pnsym->index]) a2->g[p->pnsym->index] = (itv_t*)calloc(intdim+realdim,sizeof(itv_t));
       itv_init(a2->g[p->pnsym->index][i]);
       itv_set(a2->g[p->pnsym->index][i],p->coeff);
       }
       }
     */

    size_t old = pr->dim;
    for (i = 0; i < pr->epssize; i++) {
      if (i != pr->inputns[i]) {
	t1p_nsym_t *tmp;
	tmp = pr->epsilon[i];
	pr->epsilon[i] = pr->epsilon[pr->inputns[i]];
	pr->epsilon[pr->inputns[i]] = tmp;
	pr->epsilon[i]->index = i;
	pr->epsilon[pr->inputns[i]]->index = pr->inputns[i];
	pr->inputns[i] = i;
      }
    }
    pr->mubGlobal.p = (Tobj1 *) calloc (old, sizeof (Tobj1));
    for (i = pr->epssize; i < old; i++) {
      pr->mubGlobal.p[i].x = pr->epsilon[i];	/* P^X + P^Y */
      pr->mubGlobal.p[i].y = t1p_nsym_add (pr, UN);	/* P^X - P^Y */
    }
    pr->mubGlobal.cx = t1p_nsym_add (pr, UN);
    for (i = 0; i < pr->epssize; i++) {
      pr->mubGlobal.p[i].x = t1p_nsym_add (pr, UN);	/* C^X - C^Y */
    }
    if (a1->hypercube && a2->hypercube) {
      for (i = 0; i < (intdim + realdim); i++) {
	t1p_aff_canonical (pr, a1->paf[i]);
	t1p_aff_canonical (pr, a2->paf[i]);
	//printf("%d: ",i);
	if (t1p_aff_is_bottom (pr, a1->paf[i]))
	  res->paf[i] = a2->paf[i];
	else if (t1p_aff_is_bottom (pr, a2->paf[i]))
	  res->paf[i] = a1->paf[i];
	else if (t1p_aff_is_top (pr, a1->paf[i]) || t1p_aff_is_top (pr, a2->paf[i]))
	  res->paf[i] = pr->top;
	/*
	   else if (t1p_aff_is_eq(pr, a1->paf[i], a2->paf[i])) {
	   res->paf[i] = a1->paf[i];
	   } else if (t1p_aff_is_leq_constrained(pr, a2->paf[i], a1->paf[i], a2,a1)) {
	   res->paf[i] = a1->paf[i];
	   } else if (t1p_aff_is_leq_constrained(pr, a1->paf[i], a2->paf[i], a1, a2)) {
	   res->paf[i] = a2->paf[i];
	   } */
	else {
	  if (itv_has_infty_bound (a1->box[i]) || itv_has_infty_bound (a2->box[i])) {
	    /* Do nothing, the join of concretisations is already done and stored in res->box */
	    res->paf[i] = t1p_aff_alloc_init (pr);
	    itv_set (res->paf[i]->c, res->box[i]);
	  }
	  else {
	    /* join two affine form expressions */
	    itv_set (a1->paf[i]->itv, a1->box[i]);
	    itv_set (a2->paf[i]->itv, a2->box[i]);
	    //res->paf[i] = t1p_aff_join_arXiv2(pr, a1->paf[i], a2->paf[i], a1, a2, res);
	    //res->paf[i] = t1p_aff_join_arXiv2bis(pr, a1->paf[i], a2->paf[i], a1, a2, res);
	    //res->paf[i] = t1p_aff_join_arXiv2ter(pr, a1->paf[i], a2->paf[i], a1, a2, res);
	    res->paf[i] = t1p_aff_join_bub (pr, a1->paf[i], a2->paf[i], a1, a2, res);
	    //res->paf[i] = t1p_aff_join_bub_argmin(pr, a1->paf[i], a2->paf[i], a1, a2, res);
	  }
	}
	//fprintf(stream,"********************************************\n");
	//itv_fprint(stream,res->paf[i]->end->coeff);fprintf(stream,"\n");
	//fprintf(stream,"********************************************\n");
	res->paf[i]->pby++;
      }
      size_t k = 0;
      res->g = (itv_t **) calloc (1 + pr->dim, sizeof (itv_t *));	/* centre + au max chaque symbole lui est associ\E9 un g\E9n\E9rateur */
      res->gn = pr->dim;

      t1p_aaterm_t *p = NULL;

      /* init res->g (structure creuse) */
      for (i = 0; i < intdim + realdim; i++) {
	if (!res->g[0])
	  res->g[0] = (itv_t *) calloc (intdim + realdim, sizeof (itv_t));
	itv_init (res->g[0][i]);
	itv_set (res->g[0][i], res->paf[i]->c);
	for (p = res->paf[i]->q; p; p = p->n) {
	  if (!res->g[1 + p->pnsym->index])
	    res->g[1 + p->pnsym->index] =
	      (itv_t *) calloc (intdim + realdim, sizeof (itv_t));
	  itv_init (res->g[1 + p->pnsym->index][i]);
	  itv_set (res->g[1 + p->pnsym->index][i], p->coeff);
	}
      }
      /* virer les vecteurs colineaires */
      size_t j = 0;
      itv_t un;
      itv_init (un);
      itv_set_int (un, 1);
      for (i = 1 + pr->epssize; i < 1 + (-1 + pr->dim); i++) {
	for (j = i + 1; j < 1 + pr->dim; j++) {
	  if (res->g[i] && res->g[j]) {
	    if (res->g[j][0] && res->g[i][0])
	      itv_div (pr->itv, tmp, res->g[j][0], res->g[i][0]);
	    else
	      break;
	    for (k = 1; k < intdim + realdim;) {
	      if (res->g[j][k] && res->g[i][k])
		itv_div (pr->itv, tmp1, res->g[j][k], res->g[i][k]);
	      else
		break;
	      if (itv_is_eq (tmp, tmp1))
		k++;
	      else
		break;
	    }
	    if (k == intdim + realdim) {
	      itv_abs (tmp1, tmp1);
	      itv_add (tmp1, tmp1, un);
	      for (k = 0; k < intdim + realdim; k++) {
		itv_mul (pr->itv, res->g[i][k], tmp1, res->g[i][k]);
		itv_clear (res->g[j][k]);
	      }
	      free (res->g[j]);
	      res->g[j] = NULL;
	    }
	  }
	}
      }
      itv_clear (un);

      /* FILE* poly = fopen("/home/donquijote/taylor1p/taylor1plus/taylor1plus/demo/draw.poly", "a"); */
      /*  fprintf(poly,"$z%zu=zonotope(new Matrix<Rational>([",pr->it); */
      /*  for (i=0;i<pr->dim;i++) { */
      /*      if (res->g[i+1]) { */
      /*      fprintf(poly,"[0,"); */
      /*          for (k=0;k<3;k++) { */
      /*              if (res->g[i+1][k] == NULL) {fprintf(poly,"0");} */
      /*              else if (itv_is_zero(res->g[i+1][k])) {fprintf(poly,"0");} */
      /*              else {bound_fprint(poly,res->g[i+1][k]->sup);} */
      /*              if (k<2) fprintf(poly,","); */
      /*          } */
      /*      fprintf(poly,"],"); */
      /*      } */
      /*  } */
      /*  fprintf(poly,"]));\n"); */
      /*  fprintf(poly,"$p%zu=new Polytope<Rational>(POINTS=>$z%zu);\n",pr->it,pr->it); */
      /*  fprintf(poly,"$pp%zu=transform($p%zu,new Matrix<Rational>([[1,",pr->it,pr->it); */
      /*      for (k=0;k<3;k++) { */
      /*          if (res->g[0][k] == NULL) {fprintf(poly,"0");} */
      /*          else if (itv_is_zero(res->g[0][k])) {fprintf(poly,"0");} */
      /*          else {bound_fprint(poly,res->g[0][k]->sup);} */
      /*          if (k<2) fprintf(poly,","); */
      /*      } */
      /*  fprintf(poly,"],[0,1,0,0],[0,0,1,0],[0,0,0,1]]));\n"); */
      /*  fprintf(poly,"javaview($pp%zu->VISUAL,File=>\"pp%zu\");\n",pr->it,pr->it); */

      /* fclose(poly); */

    }
    else {
      size_t k = 0;
      ap_dim_t j = 0;
      size_t dims1 = t1p_nsymcons_get_dimension (pr, a1);
      size_t dims2 = t1p_nsymcons_get_dimension (pr, a2);
      if (dims1 && dims2) {
	size_t dim2 = 0;
	/* au pire il y a toute la liste de l'un a rajouter dans l'autre */
	ap_dimchange_t *dimchange2 = ap_dimchange_alloc (0, dims2);
	if (dims1 > res->size) {
	  res->nsymcons =
	    (ap_dim_t *) realloc (res->nsymcons, (dims1) * sizeof (ap_dim_t));
	  res->gamma =
	    (ap_interval_t **) realloc (res->gamma,
					(dims1) * sizeof (ap_interval_t *));
	  for (k = res->size; k < dims1; k++)
	    res->gamma[k] = NULL;
	  res->size = dims1;
	}
	res->nsymcons =
	  memcpy ((void *) res->nsymcons, (const void *) a1->nsymcons,
		  dims1 * sizeof (ap_dim_t));
	for (i = 0; i < dims1; i++)
	  res->gamma[i] = ap_interval_alloc_set (a1->gamma[i]);
	ap_abstract0_free (pr->manNS, res->abs);
	res->abs = ap_abstract0_copy (pr->manNS, a1->abs);
	for (k = 0; k < dims1; k++) {
	  if (!t1p_nsymcons_get_dimpos (pr, &j, a1->nsymcons[k], a2)) {
	    dimchange2->dim[dim2] = j;
	    dim2++;
	  }
	}
	dimchange2->realdim = dim2;
	size_t dim1 = 0;
	for (k = 0; k < dims2; k++)
	  t1p_insert_constrained_nsym (pr, &j, a2->nsymcons[k], res);

	/* destructive, without projection (new dimension set to top) */
	ap_abstract0_add_dimensions (pr->manNS, true, a2->abs, dimchange2, false);

	ap_abstract0_join (pr->manNS, true, res->abs, a2->abs);

	/* update res->gamma */
	t1p_update_nsymcons_gamma (pr, res);

	ap_dimchange_add_invert (dimchange2);
	ap_abstract0_remove_dimensions (pr->manNS, true, a2->abs, dimchange2);

	dimchange2->realdim = dims2;
	ap_dimchange_free (dimchange2);

	size_t nsymcons_size = t1p_nsymcons_get_dimension (pr, res);
	pr->dimtoremove =
	  (ap_dim_t *) realloc (pr->dimtoremove,
				(nsymcons_size) * sizeof (ap_dim_t));
	memset ((void *) pr->dimtoremove, (int) 0, nsymcons_size * sizeof (int));
      }
      else {
	/* res->abs is a hypercube */
      }
      for (i = 0; i < (intdim + realdim); i++) {
	if (t1p_aff_is_bottom (pr, a1->paf[i]))
	  res->paf[i] = a2->paf[i];
	else if (t1p_aff_is_bottom (pr, a2->paf[i]))
	  res->paf[i] = a1->paf[i];
	else if (t1p_aff_is_top (pr, a1->paf[i]) || t1p_aff_is_top (pr, a2->paf[i]))
	  res->paf[i] = pr->top;
	else if (t1p_aff_is_eq (pr, a1->paf[i], a2->paf[i]))
	  res->paf[i] = a1->paf[i];
	else {
	  if (itv_has_infty_bound (a1->box[i]) || itv_has_infty_bound (a2->box[i])) {
	    /* Do nothing, the join of concretisations is already done and stored in res->box */
	    res->paf[i] = t1p_aff_alloc_init (pr);
	    itv_set (res->paf[i]->c, res->box[i]);
	  }
	  else {
	    /* join two affine form expressions */
	    itv_set (a1->paf[i]->itv, a1->box[i]);
	    itv_set (a2->paf[i]->itv, a2->box[i]);
	    res->paf[i] =
	      t1p_aff_join_constrained8 (pr, a1->paf[i], a2->paf[i], a1, a2, res);
	  }
	}
	res->paf[i]->pby++;
      }

      man->result.flag_best = tbool_top;
      man->result.flag_exact = tbool_top;

    }
    man->result.flag_best = tbool_true;
    man->result.flag_exact = tbool_top;
    itv_clear (tmp);
    itv_clear (tmp1);
  }
  man->result.flag_best = tbool_true;
  man->result.flag_exact = tbool_true;


#ifdef _T1P_DEBUG
  fprintf (stdout, "### RESULT of JOIN (des %d) [%tx] ###\n", destructive,
	   (intptr_t) res);
  //t1p_fprint(stdout, man, res, NULL);
  /*
     for (i=0; i<intdim+realdim;i++)
     {fprintf(stream,"%d :",i);itv_fprint(stream,res->box[i]); fprintf(stream,"\n");}
     fprintf(stream,"** %d **\n", pr->dim);
     fprintf(stream, "### ### ###\n");
   */
  double inf, sup;
  double_set_num (&inf, bound_numref (res->box[0]->inf));
  double_set_num (&sup, bound_numref (res->box[0]->sup));
  //fprintf(stream,"-%f \t %f \t %zu\n",inf, sup, pr->dim);
  //num_fprint(stream,bound_numref(res->box[0]->inf));fprintf(stream,"\t");num_fprint(stream,bound_numref(res->box[0]->sup));fprintf(stream,"\t%d\n",pr->dim);
  //t1p_fprint(stream,man,res,0);
  //fclose(stream);

#endif

  pr->it++;
  return res;
}

/**************/
/* join array */
/**************/
t1p_t *
t1p_join_array (ap_manager_t * man, t1p_t ** tab, size_t size)
{
  CALL ();
  t1p_internal_t *pr = t1p_init_from_manager (man, AP_FUNID_JOIN_ARRAY);
  t1p_t *res = (t1p_t *) ap_generic_join_array (man, (void **) tab, size);
  return res;
}

/*****************/
/* add ray array */
/*****************/
t1p_t *
t1p_add_ray_array (ap_manager_t * man, bool destructive, t1p_t * a,
		   ap_generator0_array_t * array)
{
  CALL ();
  t1p_internal_t *pr = t1p_init_from_manager (man, AP_FUNID_ADD_RAY_ARRAY);
  not_implemented ();
}
