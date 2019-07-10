/*
   APRON Library / Taylor1+ Domain (beta version)
   Copyright (C) 2009-2011 Khalil Ghorbal

*/


#include "t1p_internal.h"
#include "t1p_representation.h"
//#include "t1p_constructor.h"

#include "ap_dimension.h"
#include "ap_lincons0.h"
#include "ap_manager.h"
#include "ap_interval.h"
#include "ap_tcons0.h"
#include "ap_generator0.h"

#include "tilings.h"
#include <glpk.h>

#define M 3
#define N 2
#define LDA M
#define LDU M
#define LDVT N

/****************/
/* Constructors */
/****************/
/* 1.Basic constructors */
t1p_t *
t1p_bottom (ap_manager_t * man, size_t intdim, size_t realdim)
{
  CALL ();
  size_t i;
  size_t dims = intdim + realdim;
  t1p_internal_t *pr = t1p_init_from_manager (man, AP_FUNID_BOTTOM);
  t1p_t *res = t1p_alloc (man, intdim, realdim);
  for (i = 0; i < dims; i++) {
    res->paf[i] = pr->bot;
    res->paf[i]->pby++;
    itv_set_bottom (res->box[i]);
  }
  man->result.flag_best = tbool_true;
  man->result.flag_exact = tbool_true;
  return res;
}

t1p_t *
t1p_top (ap_manager_t * man, size_t intdim, size_t realdim)
{
  CALL ();
  size_t i;
  size_t dims = intdim + realdim;
  t1p_internal_t *pr = t1p_init_from_manager (man, AP_FUNID_TOP);
  t1p_t *res = t1p_alloc (man, intdim, realdim);
  for (i = 0; i < dims; i++) {
    res->paf[i] = pr->top;
    res->paf[i]->pby++;
    itv_set_top (res->box[i]);
  }
  man->result.flag_best = tbool_true;
  man->result.flag_exact = tbool_true;
  return res;
}

/* Abstract an hypercube defined by the array of intervals of size intdim+realdim */
t1p_t *
t1p_of_box (ap_manager_t * man, size_t intdim, size_t realdim,
	    ap_interval_t ** tinterval)
{
  CALL ();
  t1p_internal_t *pr = t1p_init_from_manager (man, AP_FUNID_OF_BOX);
  itv_t *itv_array = itv_array_alloc (intdim + realdim);
  itv_array_set_ap_interval_array (pr->itv, &itv_array, tinterval, intdim + realdim);
  t1p_t *res = t1p_alloc (man, intdim, realdim);
  size_t i = 0;
  for (i = 0; i < intdim + realdim; i++) {
    itv_set (res->box[i], itv_array[i]);
    res->paf[i] = t1p_aff_alloc_init (pr);
    if (itv_is_bottom (pr->itv, itv_array[i]))
      res->paf[i] = pr->bot;
    else if (itv_is_top (itv_array[i]))
      res->paf[i] = pr->top;
    else if (itv_is_point (pr->itv, itv_array[i]))
      itv_set (res->paf[i]->c, itv_array[i]);
    else if (itv_has_finite_bound (itv_array[i]))
      t1p_aff_add_itv (pr, res->paf[i], itv_array[i], UN);	//// modified to UN, earlier it was IN //// (modified by Bibek)
    else
      itv_set (res->paf[i]->c, itv_array[i]);
    res->paf[i]->pby++;
  }
  itv_array_free (itv_array, intdim + realdim);
  man->result.flag_best = tbool_true;
  man->result.flag_exact = tbool_true;
  return res;
}

/* 2.Accessors */
/* return dimension (type ap_dimension_t) */
ap_dimension_t
t1p_dimension (ap_manager_t * man, t1p_t * a)
{
  CALL ();
  ap_dimension_t res;
  res.intdim = a->intdim;
  res.realdim = a->dims - a->intdim;
  man->result.flag_best = tbool_true;
  man->result.flag_exact = tbool_true;
  return res;
}

/* 3.Tests */
tbool_t
t1p_is_bottom (ap_manager_t * man, t1p_t * a)
{
  CALL ();
  size_t i;
  t1p_internal_t *pr = t1p_init_from_manager (man, AP_FUNID_IS_BOTTOM);
  bool res = itv_is_bottom (pr->itv, a->box[0]);
  for (i = 1; i < a->dims; i++) {
    res &= itv_is_bottom (pr->itv, a->box[i]);
  }
  man->result.flag_best = tbool_true;
  man->result.flag_exact = tbool_true;
  return tbool_of_bool (res);
}

tbool_t
t1p_is_top (ap_manager_t * man, t1p_t * a)
{
  CALL ();
  t1p_internal_t *pr = t1p_init_from_manager (man, AP_FUNID_IS_TOP);
  size_t i;
  bool res = itv_is_top (a->box[0]);
  for (i = 1; i < a->dims; i++) {
    res &= itv_is_top (a->box[i]);
  }
  man->result.flag_best = tbool_true;
  man->result.flag_exact = tbool_true;
  return tbool_of_bool (res);
}

tbool_t
t1p_is_eq (ap_manager_t * man, t1p_t * a, t1p_t * b)
{
  CALL ();
  t1p_internal_t *pr = t1p_init_from_manager (man, AP_FUNID_IS_EQ);
  arg_assert (a && b && (a->dims == b->dims), abort ();
    );
  man->result.flag_best = tbool_true;
  man->result.flag_exact = tbool_true;
  ap_dimension_t dim1 = ap_abstract0_dimension (pr->manNS, a->abs);
  ap_dimension_t dim2 = ap_abstract0_dimension (pr->manNS, b->abs);
  if ((dim1.intdim != dim2.intdim) || (dim1.realdim != dim2.realdim))
    return tbool_of_bool (false);
  else if (!ap_abstract0_is_eq (pr->manNS, a->abs, b->abs))
    return tbool_of_bool (false);
  else if (a == b)
    return tbool_of_bool (true);
  else {
    size_t i = 0;
    bool res = true;
    for (i = 0; i < a->dims; i++) {
      if (a->paf[i] != b->paf[i] && itv_is_eq (a->box[i], b->box[i])) {
	res = t1p_aff_is_eq (pr, a->paf[i], b->paf[i]);
	if (!res)
	  break;
      }
      else {
	res = false;
	break;
      }
    }
    return tbool_of_bool (res);
  }
}

void comb(int m, int n, unsigned char *c, int arr1[100])
{
	int i;
	for (i = 0; i < n; i++) c[i] = n - i;
        int p=0;
 
	while (1) {
		for (i = n; i--;){
			//printf("%d%c", c[i], i ? ' ': '\n');
                        arr1[p]=c[i];
                        //printf("%d", c[i]);
                        //printf("%d",arr[p]);
                        p++;}

           
                 
		/* this check is not strictly necessary, but if m is not close to n,
		   it makes the whole thing quite a bit faster */
		i = 0;
		if (c[i]++ < m) continue;
 
		for (; c[i] >= m - i;) if (++i >= n) return;
		for (c[i]++; i; i--) c[i-1] = c[i] + 1;
	}

}

int binomialCoeff(int n, int k)
{
  // Base Cases
  if (k==0 || k==n)
    return 1;
 
  // Recur
  return  binomialCoeff(n-1, k-1) + binomialCoeff(n-1, k);
}

tbool_t
t1p_is_leq (ap_manager_t * man, t1p_t * a, t1p_t * b)
{
  //int b1;
  //scanf ("%f", &b1);
  CALL ();
  t1p_internal_t *pr = t1p_init_from_manager (man, AP_FUNID_IS_LEQ);
  int c;
  #ifdef _T1P_DEBUG
    t1p_fprint(stdout, man, a, 0x0);
    t1p_fprint(stdout, man, b, 0x0);
    fprintf(stdout, "### ### ###\n");
  #endif
  scanf("%f",&c);
  bool res;

  if(t1p_is_top (man, b)){
     res=true;
     return res;}
  if(t1p_is_top (man, a)){
    res=false;
    return res;}

  //t1p_aff_update(pr, a);
  //t1p_aff_update(pr, b);

  int size, size1, size2, i, j, k, l;
  int sum = 0;
  double num2, num3;

  //itv_print (a->paf[0]->q->coeff);

  ///////// first generator matrix ///////////////////
  itv_t **ptitv1 = (itv_t **) calloc (a->dims, sizeof (itv_t *));
  for (i = 0; i < a->dims; i++)
    ptitv1[i] = (itv_t *) calloc (pr->dim, sizeof (itv_t));
  //size1 = t1p_aff_get_generator_matrix_size (pr, a, ptitv1);

  /*for (i=0;i<a->dims;i++){
   for (j=0;j<pr->dim;j++){
       itv_init(ptitv1[i][j]);//=coeff1;
       //itv_print(ptitv[i][j]);
     }
}*/

//itv_print(ptitv[0][0]);

t1p_aaterm_t *p = NULL;

int* vec= (int* )calloc(pr->dim,sizeof(int));

for (i=0;i<pr->dim;i++){
vec[i]=-1;
}

k=0;

for(i=0;i<a->dims;i++) {
 for(p=a->paf[i]->q;p;p=p->n) {
		//itv_init(ptitv1[i][p->pnsym->index]);
		itv_set(ptitv1[i][p->pnsym->index],p->coeff);
                //itv_print(ptitv[i][p->pnsym->index]);
                if (vec[p->pnsym->index]==-1){
                   vec[p->pnsym->index]=k;
                   k++;}
   }
}

size1=k;

  //printf ("%u", size1);
  itv_t **ptitv2 = (itv_t **) calloc (a->dims, sizeof (itv_t *));
  for (i = 0; i < a->dims; i++)
    ptitv2[i] = (itv_t *) calloc (size1, sizeof (itv_t));

  t1p_aff_get_generator_matrix (pr, a, ptitv1, ptitv2);

/*for (i=0;i<a->dims;i++){
   for (j=0;j<k;j++){
       itv_init(ptitv2[i][j]);//=coeff1;
       //itv_print(ptitv[i][j]);
     }
}*/

t1p_aaterm_t *p1 = NULL;

for(i=0;i<a->dims;i++) {
 for(p1=a->paf[i]->q;p1;p1=p1->n) {
		//itv_init(ptitv2[i][vec[p1->pnsym->index]]);
		itv_set(ptitv2[i][vec[p1->pnsym->index]],p1->coeff);
                //itv_print(ptitv1[i][p1->pnsym->index]);
                  
	   }
}

  //// printing the generator matrices //////////
  /*scanf ("%f", &b1);
     itv_print(a->paf[0]->c);
     itv_print(a->paf[1]->c);
     scanf ("%f", &b1);
     for (i = 0; i < a->dims; i++) {
     for (j = 0; j < size1; j++) {
     itv_print (ptitv2[i][j]);
     }
     } 
     scanf ("%f", &b1); */

  //////////////// second generator matrix //////////////////
  itv_t **ptitv3 = (itv_t **) calloc (b->dims, sizeof (itv_t *));
  for (i = 0; i < b->dims; i++)
    ptitv3[i] = (itv_t *) calloc (pr->dim, sizeof (itv_t));
  //size2 = t1p_aff_get_generator_matrix_size (pr, b, ptitv3);

  /*for (i=0;i<b->dims;i++){
   for (j=0;j<pr->dim;j++){
       itv_init(ptitv3[i][j]);//=coeff1;
       //itv_print(ptitv[i][j]);
     }
}*/

//itv_print(ptitv[0][0]);

//t1p_aaterm_t *p = NULL;

//int* vec= (int* )calloc(pr->dim,sizeof(int));

for (i=0;i<pr->dim;i++){
vec[i]=-1;
}

k=0;

for(i=0;i<b->dims;i++) {
 for(p=b->paf[i]->q;p;p=p->n) {
		//itv_init(ptitv3[i][p->pnsym->index]);
		itv_set(ptitv3[i][p->pnsym->index],p->coeff);
                //itv_print(ptitv[i][p->pnsym->index]);
                if (vec[p->pnsym->index]==-1){
                   vec[p->pnsym->index]=k;
                   k++;}
   }
}

size2=k;

  //printf ("%u", size2);
  itv_t **ptitv4 = (itv_t **) calloc (b->dims, sizeof (itv_t *));
  for (i = 0; i < b->dims; i++)
    ptitv4[i] = (itv_t *) calloc (size2, sizeof (itv_t));

  t1p_aff_get_generator_matrix (pr, b, ptitv3, ptitv4);

  /*for (i=0;i<b->dims;i++){
   for (j=0;j<k;j++){
       itv_init(ptitv4[i][j]);//=coeff1;
       //itv_print(ptitv[i][j]);
     }
}*/

//t1p_aaterm_t *p1 = NULL;

for(i=0;i<b->dims;i++) {
 for(p1=b->paf[i]->q;p1;p1=p1->n) {
		//itv_init(ptitv4[i][vec[p1->pnsym->index]]);
		itv_set(ptitv4[i][vec[p1->pnsym->index]],p1->coeff);
                //itv_print(ptitv1[i][p1->pnsym->index]);
                  
	   }
}

  //// printing the generator matrices //////////
  /*scanf ("%f", &b1);
     itv_print(b->paf[0]->c);
     itv_print(b->paf[1]->c);
     scanf ("%f", &b1);
     for (i = 0; i < b->dims; i++) {
     for (j = 0; j < size2; j++) {
     itv_print (ptitv4[i][j]);
     }
     }
     scanf ("%f", &b1); */

//////////////// tranpose of the matrices /////////////////

  itv_t **ptitv2trnps = (itv_t **) calloc (size1, sizeof (itv_t *));
  for (i = 0; i < size1; i++)
    ptitv2trnps[i] = (itv_t *) calloc (a->dims, sizeof (itv_t));

  for (i = 0; i < a->dims; i++) {
    for (j = 0; j < size1; j++) {
      //itv_init (ptitv2trnps[j][i]);
      itv_set (ptitv2trnps[j][i], ptitv2[i][j]);
    }
  }

  itv_t **ptitv4trnps = (itv_t **) calloc (size2, sizeof (itv_t *));
  for (i = 0; i < size2; i++)
    ptitv4trnps[i] = (itv_t *) calloc (b->dims, sizeof (itv_t));

  for (i = 0; i < b->dims; i++) {
    for (j = 0; j < size2; j++) {
      //itv_init (ptitv4trnps[j][i]);
      itv_set (ptitv4trnps[j][i], ptitv4[i][j]);
    }
  }

  //scanf ("%f", &b1);

//// printing the generator matrices //////////
  /*for (i = 0; i < a->dims; i++) {
     for (j = 0; j < size1; j++) {
     itv_print (ptitv2[i][j]);
     }
     }

     for (i = 0; i < b->dims; i++) {
     for (j = 0; j < size2; j++) {
     itv_print (ptitv4[i][j]);
     }
     }

     scanf ("%f", &b1);

     scanf ("%f", &b1); */

//// printing the transpose of generator matrices //////////
  /*for (i = 0; i < size1; i++) {
     for (j = 0; j < a->dims; j++) {
     itv_print (ptitv2trnps[i][j]);
     }
     }

     for (i = 0; i < size2; i++) {
     for (j = 0; j < b->dims; j++) {
     itv_print (ptitv4trnps[i][j]);
     }
     }

     scanf ("%f", &b1); */

/////////////  G2(2,:)=-1*G2(2,:);  /////////////////////

  itv_t coeff;
  itv_init (coeff);
  num_t num;
  num_init (num);
  num_set_double (num, -1);
  itv_set_num (coeff, num);
  num_t num1;
  num_init (num1);
  itv_t tmp;
  itv_init (tmp);
  itv_t tmp1;
  itv_init (tmp1);
  itv_t tmp2;
  itv_init (tmp2);
  itv_t tmp3;
  itv_init (tmp3);
  itv_t tmp4;
  itv_init (tmp4);
  itv_t tmp5;
  itv_init (tmp5);
  itv_t tmp6;
  itv_init (tmp6);

        int nocoeffs;
        nocoeffs=binomialCoeff(size2, 2);
        //printf("%d", nocoeffs);
        int arr1[100];
	unsigned char buf[100];
	comb(size2, 2, buf, arr1);
        int combsmat[nocoeffs][2];
      for(j=0;j<nocoeffs;j++){
        for (i=0;i<2;i++){
              combsmat[j][i]=arr1[j*2+i];}}

   itv_t **gnormal = (itv_t **) calloc (b->dims, sizeof (itv_t *));
   for (i = 0; i < b->dims; i++)
     gnormal[i] = (itv_t *) calloc (nocoeffs, sizeof (itv_t));

       //// cross product for finding normal vectors /////////////
       for (i=0; i<nocoeffs; i++) {                  
          itv_mul(pr->itv, tmp1, ptitv4[1][combsmat[i][0]-1], ptitv4[2][combsmat[i][1]-1]);
          itv_mul(pr->itv, tmp2, ptitv4[2][combsmat[i][0]-1], ptitv4[1][combsmat[i][1]-1]);
          itv_sub(gnormal[0][i], tmp1, tmp2);

          itv_mul(pr->itv, tmp3, ptitv4[2][combsmat[i][0]-1], ptitv4[0][combsmat[i][1]-1]);
          itv_mul(pr->itv, tmp4, ptitv4[0][combsmat[i][0]-1], ptitv4[2][combsmat[i][1]-1]);
          itv_sub(gnormal[1][i], tmp3, tmp4);

          itv_mul(pr->itv, tmp5, ptitv4[0][combsmat[i][0]-1], ptitv4[1][combsmat[i][1]-1]);
          itv_mul(pr->itv, tmp6, ptitv4[1][combsmat[i][0]-1], ptitv4[0][combsmat[i][1]-1]);
          itv_sub(gnormal[2][i], tmp5, tmp6);}

          /*for(i=0;i<3;i++){
             for(j=0;j<nocoeffs;j++){
             itv_print(gnormal[i][j]);}
             printf("\n");}

       scanf("%f",&c);*/

  itv_t *diffc = (itv_t *) calloc (3, sizeof (itv_t));
  itv_t *left = (itv_t *) calloc (3, sizeof (itv_t));
  itv_t *lefteq = (itv_t *) calloc (size2, sizeof (itv_t));
  itv_t *righteq = (itv_t *) calloc (size2, sizeof (itv_t));
  itv_t *righteq1 = (itv_t *) calloc (size1, sizeof (itv_t));
  itv_t *righteq2 = (itv_t *) calloc (size2, sizeof (itv_t));

  //itv_t *diffc[2];            //// 2, here is the no. of variables or two affine forms (dimension)
  //itv_t *left[2];
  //itv_t *lefteq[size2];
  //itv_t *righteq[size2];
  //itv_t *righteq1[size1];
  /*for (i = 0; i < size1; i++) {
    itv_init (righteq1[i]);
  }*/
  itv_t normrighteq1;
  itv_init (normrighteq1);
  //itv_t *righteq2[size2];
  /*for (i = 0; i < size2; i++) {
    itv_init (righteq2[i]);
  }*/
  itv_t normrighteq2;
  itv_init (normrighteq2);

  // \Bigg\lvert\langle u_{i},c_{x}-c_{y} \rangle\Bigg\rvert \leq {||Y_{+}u_{i}||}_{1} - {||X_{+}u_{i}||}_{1}, \forall i=1,\cdots,k (FMSD 2015) ////

  for (j = 0; j < nocoeffs; j++) {
    for (i = 0; i < a->dims; i++) {
      itv_init (diffc[i]);
      itv_sub (diffc[i], b->paf[i]->c, a->paf[i]->c);
      itv_init (left[i]);
      itv_mul (pr->itv, left[i], diffc[i], gnormal[i][j]);
      //scanf ("%f", &b1);
      //itv_print (left[i]);
      //scanf ("%f", &b1);
    }
    ////////// lefteq(j)=abs(sum(left)) /////////////
    itv_init (lefteq[j]);
    itv_add (lefteq[j], left[0], left[1]);
    if (itv_is_neg (lefteq[j]) == 1)
      itv_mul (pr->itv, lefteq[j], coeff, lefteq[j]);
    //printf ("lefteq");
    //scanf ("%f", &b1);
    //itv_print (lefteq[j]);
    //scanf ("%f", &b1);
    itv_init (normrighteq1);
    for (k = 0; k < size1; k++) {    /// size1? doubtful
      itv_init (righteq1[k]);
      for (l = 0; l < a->dims; l++) {
	itv_mul (pr->itv, tmp, ptitv2trnps[k][l], gnormal[l][j]);
	//itv_print (tmp);
	itv_add (righteq1[k], righteq1[k], tmp);
	//printf ("righteq1");
	//scanf ("%f", &b1);
	//itv_print (righteq1[k]);
	//scanf ("%f", &b1);
      }
      //righteq1[k] += ptitv2trnps[l][k] * ptitv4[l][j];}
      if (itv_is_neg (righteq1[k]) == 1) {
	itv_mul (pr->itv, righteq1[k], coeff, righteq1[k]);
      }
      itv_add (normrighteq1, normrighteq1, righteq1[k]);
      //printf ("normrighteq1");
      //scanf ("%f", &b1);
      //itv_print (normrighteq1);
      //scanf ("%f", &b1);
    }
    itv_init (normrighteq2);
    for (k = 0; k < size2; k++) {    /// size2? doubtful
      itv_init (righteq2[k]);
      for (l = 0; l < b->dims; l++) {
	itv_mul (pr->itv, tmp, ptitv4trnps[k][l], gnormal[l][j]);
	//itv_print (tmp);
	itv_add (righteq2[k], righteq2[k], tmp);
	//printf ("righteq2");
	//scanf ("%f", &b1);
	//itv_print (righteq2[k]);
	//scanf ("%f", &b1);
      }
      //righteq2[k] += ptitv4trnps[l][k] * ptitv4[l][j];}
      if (itv_is_neg (righteq2[k]) == 1) {
	itv_mul (pr->itv, righteq2[k], coeff, righteq2[k]);
      }
      itv_add (normrighteq2, normrighteq2, righteq2[k]);
      //printf ("normrighteq2");
      //scanf ("%f", &b1);
      //itv_print (normrighteq2);
      //scanf ("%f", &b1);
    }
    //printf ("righteq");
    itv_init (righteq[j]);
    itv_sub (righteq[j], normrighteq2, normrighteq1);
    itv_magnitude (num, lefteq[j]);
    double_set_num (&num2, num);
    if (itv_is_neg (lefteq[j]) == 1)
      num2 = -1 * num2;

    //itv_print (righteq[j]);
    itv_magnitude (num1, righteq[j]);
    double_set_num (&num3, num1);
    if (itv_is_neg (righteq[j]) == 1)
      num3 = -1 * num3;
    if (num2 <= num3) {
      sum++;
    }
    else {
      res = false;
      //printf ("%d", res);
      //scanf ("%f", &c);
      itv_clear (coeff);
      itv_clear (tmp);
      itv_clear (tmp1);
      itv_clear (tmp2);
      itv_clear (tmp3);
      itv_clear (tmp4);
      itv_clear (tmp5);
      itv_clear (tmp6);
      itv_clear (normrighteq1);
      itv_clear (normrighteq2);      
      for (i = 0; i < a->dims; i++) {
    for (j = 0; j < pr->dim; j++) {
      itv_clear (ptitv1[i][j]);}}
      for (i = 0; i < a->dims; i++){
         free(ptitv1[i]);}
      free (ptitv1);      
      for (i = 0; i < a->dims; i++) {
    for (j = 0; j < size1; j++) {
      itv_clear (ptitv2[i][j]);
      itv_clear (ptitv2trnps[j][i]);}}
      for (i = 0; i < a->dims; i++){
         free(ptitv2[i]);
         free(ptitv2trnps[i]);}
      free (ptitv2);
      free (ptitv2trnps);      
      for (i = 0; i < b->dims; i++) {
    for (j = 0; j < pr->dim; j++) {
      itv_clear (ptitv3[i][j]);}}
      for (i = 0; i < b->dims; i++){
         free(ptitv3[i]);}
      free (ptitv3);
      for (i = 0; i < b->dims; i++) {
      for (j = 0; j < size2; j++) {
      itv_clear (ptitv4[i][j]);
      itv_clear (ptitv4trnps[j][i]);}}
      for (i = 0; i < b->dims; i++){
         free(ptitv4[i]);
         free(ptitv4trnps[i]);}
      free (ptitv4);
            free (ptitv4trnps);
      for (i = 0; i < b->dims; i++) {
      for (j = 0; j < nocoeffs; j++) {
      itv_clear (gnormal[i][j]);}}
      for (i = 0; i < b->dims; i++){
         free(gnormal[i]);}
      free(gnormal);
      /*for (i = 0; i < a->dims; i++) {
    for (j = 0; j < size1; j++) {
      itv_clear (ptitv2trnps[j][i]);}}
      for (i = 0; i < size1; i++){
         free(ptitv2trnps[i]);}          
      free (ptitv2trnps);
      for (i = 0; i < b->dims; i++) {
    for (j = 0; j < size2; j++) {
      itv_clear (ptitv4trnps[j][i]);}}
      for (i = 0; i < size2; i++){
         free(ptitv4trnps[i]);}
      free (ptitv4trnps);*/
      for (i = 0; i < size1; i++) {
	itv_clear (righteq1[i]);
      }
      free (righteq1);
      for (i = 0; i < size2; i++) {
	itv_clear (righteq2[i]);
	itv_clear (lefteq[i]);
	itv_clear (righteq[i]);
      }
      free (lefteq);
      free (righteq);
      free (righteq2);
      for (i = 0; i < a->dims; i++) {
	itv_clear (diffc[i]);
	itv_clear (left[i]);
      }
      free (diffc);
      free (left);
      free (vec);
      return res;
    }
    if (sum == 3)  //// nocoeffs? Because, that's the number of new normal vectors obtained after doing cross product or SVD
      res = true;


  }

  //printf ("%d", res);
  //scanf ("%f", &c);
  itv_clear (coeff);
      itv_clear (tmp);
      itv_clear (tmp1);
      itv_clear (tmp2);
      itv_clear (tmp3);
      itv_clear (tmp4);
      itv_clear (tmp5);
      itv_clear (tmp6);
      itv_clear (normrighteq1);
      itv_clear (normrighteq2);      
      for (i = 0; i < a->dims; i++) {
    for (j = 0; j < pr->dim; j++) {
      itv_clear (ptitv1[i][j]);}}
      for (i = 0; i < a->dims; i++){
         free(ptitv1[i]);}
      free (ptitv1);      
      for (i = 0; i < a->dims; i++) {
    for (j = 0; j < size1; j++) {
      itv_clear (ptitv2[i][j]);
      itv_clear (ptitv2trnps[j][i]);}}
      for (i = 0; i < a->dims; i++){
         free(ptitv2[i]);
         free(ptitv2trnps[i]);}
      free (ptitv2); 
            free (ptitv2trnps);     
      for (i = 0; i < b->dims; i++) {
    for (j = 0; j < pr->dim; j++) {
      itv_clear (ptitv3[i][j]);}}
      for (i = 0; i < b->dims; i++){
         free(ptitv3[i]);}
      free (ptitv3);
      for (i = 0; i < b->dims; i++) {
    for (j = 0; j < size2; j++) {
      itv_clear (ptitv4[i][j]);
      itv_clear (ptitv4trnps[j][i]);}}
      for (i = 0; i < b->dims; i++){
         free(ptitv4[i]);
          free(ptitv4trnps[i]);}
      free (ptitv4);
            free (ptitv4trnps);
      for (i = 0; i < b->dims; i++) {
      for (j = 0; j < nocoeffs; j++) {
      itv_clear (gnormal[i][j]);}}
      for (i = 0; i < b->dims; i++){
         free(gnormal[i]);}
      free(gnormal);
      /*for (i = 0; i < a->dims; i++) {
    for (j = 0; j < size1; j++) {
      itv_clear (ptitv2trnps[j][i]);}}
      for (i = 0; i < size1; i++){
         free(ptitv2trnps[i]);}          
      free (ptitv2trnps);
      for (i = 0; i < b->dims; i++) {
    for (j = 0; j < size2; j++) {
      itv_clear (ptitv4trnps[j][i]);}}
      for (i = 0; i < size2; i++){
         free(ptitv4trnps[i]);}
      free (ptitv4trnps);*/
      for (i = 0; i < size1; i++) {
	itv_clear (righteq1[i]);
      }
      free (righteq1);
      for (i = 0; i < size2; i++) {
	itv_clear (righteq2[i]);
	itv_clear (lefteq[i]);
	itv_clear (righteq[i]);
      }
      free (lefteq);
      free (righteq);
      free (righteq2);
      for (i = 0; i < a->dims; i++) {
	itv_clear (diffc[i]);
	itv_clear (left[i]);
      }
      free (diffc);
      free (left);
      free (vec);
  return res;

}

tbool_t
t1p_is_dimension_unconstrained (ap_manager_t * man, t1p_t * a, ap_dim_t dim)
{
  CALL ();
  t1p_internal_t *pr =
    t1p_init_from_manager (man, AP_FUNID_IS_DIMENSION_UNCONSTRAINED);
  not_implemented ();
}

tbool_t
t1p_sat_tcons (ap_manager_t * man, t1p_t * a, ap_tcons0_t * tcons)
{
  CALL ();
  t1p_internal_t *pr = t1p_init_from_manager (man, AP_FUNID_SAT_TCONS);
  not_implemented ();
}

tbool_t
t1p_sat_interval (ap_manager_t * man, t1p_t * a, ap_interval_t * interval)
{
  CALL ();
  t1p_internal_t *pr = t1p_init_from_manager (man, AP_FUNID_SAT_INTERVAL);
  not_implemented ();
}

tbool_t
t1p_sat_lincons (ap_manager_t * man, t1p_t * a, ap_lincons0_t * lincons)
{
  CALL ();
  t1p_internal_t *pr = t1p_init_from_manager (man, AP_FUNID_SAT_LINCONS);
  not_implemented ();
}

/* 4.Extraction of properties */
ap_interval_t *
t1p_bound_texpr (ap_manager_t * man, t1p_t * a, ap_texpr0_t * expr)
{
  CALL ();
  t1p_internal_t *pr = t1p_init_from_manager (man, AP_FUNID_BOUND_TEXPR);
  not_implemented ();
}

ap_interval_t *
t1p_bound_dimension (ap_manager_t * man, t1p_t * a, ap_dim_t dim)
{
  CALL ();
  t1p_internal_t *pr = t1p_init_from_manager (man, AP_FUNID_BOUND_DIMENSION);
  //itv_t tmp; itv_init(tmp);
  ap_interval_t *ap_itv = ap_interval_alloc ();
  //t1p_aff_boxize(pr, tmp, a->paf[dim], a);
  //ap_interval_set_itv(pr->itv, ap_itv, tmp);
  ap_interval_set_itv (pr->itv, ap_itv, a->box[dim]);
  //itv_clear(tmp);
  //ap_interval_free(ap_itv);
  return ap_itv;
}

ap_interval_t *
t1p_bound_linexpr (ap_manager_t * man, t1p_t * a, ap_linexpr0_t * expr)
{
  t1p_internal_t *pr = t1p_init_from_manager (man, AP_FUNID_BOUND_LINEXPR);
  not_implemented ();
}

ap_interval_t **
t1p_to_box (ap_manager_t * man, t1p_t * a)
{
  /* operation couteuse */
  /* la transformation ap_interval <-> itv est couteuse */
  /* TODO:savoir si la transformation est exact et mettre a jour la valeur de flag_exact */
  CALL ();
  t1p_internal_t *pr = t1p_init_from_manager (man, AP_FUNID_TO_BOX);
  ap_interval_t **res =
    (ap_interval_t **) malloc (a->dims * sizeof (ap_interval_t *));
  size_t i = 0;
  itv_t tmp;
  itv_init (tmp);
  for (i = 0; i < a->dims; i++) {
    t1p_aff_boxize (pr, tmp, a->paf[i], a);
    res[i] = ap_interval_alloc ();
    //ap_interval_set_itv(pr->itv, res[i], a->box[i]);
    ap_interval_set_itv (pr->itv, res[i], tmp);
  }
  itv_clear (tmp);
  man->result.flag_best = tbool_true;
  man->result.flag_exact = tbool_true;
  return res;
}

ap_tcons0_array_t
t1p_to_tcons_array (ap_manager_t * man, t1p_t * a)
{
  CALL ();
  t1p_internal_t *pr = t1p_init_from_manager (man, AP_FUNID_TO_TCONS_ARRAY);
  not_implemented ();
}

ap_lincons0_array_t
t1p_to_lincons_array (ap_manager_t * man, t1p_t * a)
    /* Taken from box_to_lincons_array */
    /* TODO: use constraints in eps domain to deduce constraints on variable ? */
{
  CALL ();
  size_t i;
  t1p_internal_t *pr = t1p_init_from_manager (man, AP_FUNID_TO_LINCONS_ARRAY);
  ap_lincons0_array_t array;

  size_t nbdims = a->dims;

  man->result.flag_best = tbool_true;
  man->result.flag_exact = tbool_true;
  if (nbdims == 0) {
    array = ap_lincons0_array_make (0);
  }
  else if (t1p_is_bottom (man, a) == tbool_true) {
    array = ap_lincons0_array_make (1);
    array.p[0] = ap_lincons0_make_unsat ();
  }
  else {
    size_t size;
    ap_linexpr0_t *expr;
    ap_scalar_t *scalar;
    bool point;

    size = 0;
    //itv_t* tmp = itv_array_alloc(nbdims);
    for (i = 0; i < nbdims; i++) {
      //t1p_aff_boxize(pr, tmp[i], a->paf[i], a);
      //if (!bound_infty(tmp[i]->inf)) size++;
      if (!bound_infty (a->box[i]->inf))
	size++;
      point = itv_is_point (pr->itv, a->box[i]);
      if (!point && !bound_infty (a->box[i]->sup))
	size++;
      //point = itv_is_point(pr->itv,tmp[i]);
      //if (!point && !bound_infty(tmp[i]->sup)) size++;
    }
    array = ap_lincons0_array_make (size);
    size = 0;
    for (i = 0; i < nbdims; i++) {
      point = false;
      //if (!bound_infty(tmp[i]->inf)){
      if (!bound_infty (a->box[i]->inf)) {
	expr = ap_linexpr0_alloc (AP_LINEXPR_SPARSE, 1);
	ap_coeff_set_scalar_int (&expr->p.linterm[0].coeff, 1);
	expr->p.linterm[0].dim = i;

	ap_coeff_reinit (&expr->cst, AP_COEFF_SCALAR, AP_SCALAR_DOUBLE);
	scalar = expr->cst.val.scalar;
	//ap_scalar_set_bound(scalar,tmp[i]->inf);
	ap_scalar_set_bound (scalar, a->box[i]->inf);
	point = itv_is_point (pr->itv, a->box[i]);
	//point = itv_is_point(pr->itv,tmp[i]);
	array.p[size].constyp = point ? AP_CONS_EQ : AP_CONS_SUPEQ;
	array.p[size].linexpr0 = expr;
	size++;
      }
      //if (!point && !bound_infty(tmp[i]->sup)){
      if (!point && !bound_infty (a->box[i]->sup)) {
	expr = ap_linexpr0_alloc (AP_LINEXPR_SPARSE, 1);
	ap_coeff_set_scalar_int (&expr->p.linterm[0].coeff, -1);
	expr->p.linterm[0].dim = i;

	ap_coeff_reinit (&expr->cst, AP_COEFF_SCALAR, AP_SCALAR_DOUBLE);
	ap_scalar_set_bound (expr->cst.val.scalar, a->box[i]->sup);
	//ap_scalar_set_bound(expr->cst.val.scalar,tmp[i]->sup);

	array.p[size].constyp = AP_CONS_SUPEQ;
	array.p[size].linexpr0 = expr;
	size++;
      }
    }
    //itv_array_free(tmp, nbdims);
  }
  //fprintf(stdout, "end lincons\n");
  //ap_lincons0_array_fprint (stdout , array.p[size] , 0x0);

  return array;
}

ap_generator0_array_t
t1p_to_generator_array (ap_manager_t * man, t1p_t * a)
{
  CALL ();
  t1p_internal_t *pr = t1p_init_from_manager (man, AP_FUNID_TO_GENERATOR_ARRAY);
  int c, b;
  /*#ifdef _T1P_DEBUG
    t1p_fprint(stdout, man, a, 0x0);
    fprintf(stdout, "### ### ###\n");
  #endif*/
  
  //t1p_aff_update(pr, a);

  /*#ifdef _T1P_DEBUG
    t1p_fprint(stdout, man, a, 0x0);
    fprintf(stdout, "### ### ###\n");
  #endif*/
    
  //not_implemented ();
  size_t i, j, k, size;
  size_t nbcoeffs, nbvertices, v, sizegen;
  ap_generator0_array_t array;
  ap_linexpr0_t *vertex;
  ap_scalar_t scalar;

  //scanf("%f",&b);

  size = a->dims;
  if (size == 0) {
    array = ap_generator0_array_make (0);
    return array;
  }

  ///////// generator matrix ///////////////////

  itv_t **ptitv1 = (itv_t **) calloc (a->dims, sizeof (itv_t *));
  for (i = 0; i < a->dims; i++)
    ptitv1[i] = (itv_t *) calloc (pr->dim, sizeof (itv_t));

  //sizegen = t1p_aff_get_generator_matrix_size (pr, a, ptitv1);

  for (i=0;i<a->dims;i++){
   for (j=0;j<pr->dim;j++){
       itv_init(ptitv1[i][j]);//=coeff1;
       //itv_print(ptitv[i][j]);
     }
}

//itv_print(ptitv[0][0]);

t1p_aaterm_t *p = NULL;

int* vec= (int* )calloc(pr->dim,sizeof(int));

for (i=0;i<pr->dim;i++){
vec[i]=-1;
}

k=0;

for(i=0;i<a->dims;i++) {
 for(p=a->paf[i]->q;p;p=p->n) {
		itv_init(ptitv1[i][p->pnsym->index]);
		itv_set(ptitv1[i][p->pnsym->index],p->coeff);
                //itv_print(ptitv[i][p->pnsym->index]);
                if (vec[p->pnsym->index]==-1){
                   vec[p->pnsym->index]=k;
                   k++;}
   }
}

sizegen=k;

  itv_t **ptitv2 = (itv_t **) calloc (a->dims, sizeof (itv_t *));
  for (i = 0; i < a->dims; i++)
    ptitv2[i] = (itv_t *) calloc (sizegen, sizeof (itv_t));

  //t1p_aff_get_generator_matrix (pr, a, ptitv1, ptitv2);

  for (i=0;i<a->dims;i++){
   for (j=0;j<k;j++){
       itv_init(ptitv2[i][j]);//=coeff1;
       //itv_print(ptitv[i][j]);
     }
}

t1p_aaterm_t *p1 = NULL;

for(i=0;i<a->dims;i++) {
 for(p1=a->paf[i]->q;p1;p1=p1->n) {
		itv_init(ptitv2[i][vec[p1->pnsym->index]]);
		itv_set(ptitv2[i][vec[p1->pnsym->index]],p1->coeff);
                //itv_print(ptitv1[i][p1->pnsym->index]);
                  
	   }
}

  num_t num;
  num_init (num);

  struct zonotope zon;
  zon.row = 2;
  zon.column = sizegen;
  zon.totnocols = sizegen;
  zon.nbtiles = 0;
  zon.mat = (double **) calloc (zon.row, sizeof (double *));
  for (i = 0; i < zon.row; i++)
    zon.mat[i] = (double *) calloc (2 * zon.column, sizeof (double));
  zon.center = (double *) calloc (zon.row, sizeof (double));

  itv_magnitude (num, a->paf[0]->c);
  double_set_num (&zon.center[0], num);
  itv_magnitude (num, a->paf[1]->c);
  double_set_num (&zon.center[1], num);
  if (itv_is_neg (a->paf[0]->c) == 1)
    zon.center[0] = -1 * zon.center[0];
  if (itv_is_neg (a->paf[1]->c) == 1)
    zon.center[1] = -1 * zon.center[1];

   //printf("%f", zon.center[0]);
   //printf("%f", zon.center[1]);


  //ap_abstract0_fprint(stdout, pr->manNS, a->abs, name_of_ns);
  //scanf("%f",&c);

  for (j = 0; j < sizegen; j++) {
    itv_magnitude (num, ptitv2[0][j]);
    double_set_num (&zon.mat[0][j], num);
    itv_magnitude (num, ptitv2[1][j]);
    double_set_num (&zon.mat[1][j], num);
    if (itv_is_neg (ptitv2[0][j]) == 1)
      zon.mat[0][j] = -1 * zon.mat[0][j];
    if (itv_is_neg (ptitv2[1][j]) == 1)
      zon.mat[1][j] = -1 * zon.mat[1][j];
  }

  int **sgns = (int **) calloc (2 * zon.column, sizeof (int *));
  for (i = 0; i < 2 * zon.column; i++)
    sgns[i] = (int *) calloc (zon.column, sizeof (int));
  double **negmat = (double **) calloc (2 * zon.column, sizeof (double *));
  for (i = 0; i < 2 * zon.column; i++)
    negmat[i] = (double *) calloc (zon.column, sizeof (double));

  horzcat (&zon, negmat);
  signsvecnew (&zon, sgns);

  double tmp;

  double **cordinates = (double **) calloc (zon.row, sizeof (double *));
  for (i = 0; i < zon.row; i++)
    cordinates[i] = (double *) calloc (2 * zon.column, sizeof (double));

  /*for (i=0;i<size;i++){
     printf("%f", zon.center[i]);} */

  /*for (i=0;i<size;i++){
     for(k=0;k<zon.column;k++){
     printf("%f",zon.mat[i][k]);}
     printf("\n");}
     scanf("%f",&b); */

  for (i = 0; i < size; i++) {
    for (j = 0; j < 2 * zon.column; j++) {
      tmp = 0.0;
      for (k = 0; k < zon.column; k++) {
	tmp = tmp + sgns[j][k] * zon.mat[i][k];
      }
      cordinates[i][j] = zon.center[i] + tmp;
    }
  }
  //printf("%f", cordinates[i][j]);}
  //printf("\n");}

  /*scanf("%f",&b);

     for (i=0;i<size;i++){
     for (j=0;j<2*zon.column;j++){
     printf("%f", cordinates[i][j]);}
     printf("\n");}

     scanf("%f",&b); */

  nbcoeffs = 2;
  nbvertices = 2 * zon.column;
  /* Preparation */
  array = ap_generator0_array_make (nbvertices);
  ap_scalar_init (&scalar, AP_SCALAR_DOUBLE);
  ap_scalar_set_double (&scalar, 0.0);
  /* Let's go now ! */
  v = 2 * zon.column;
  /* Creates the vertices */
  vertex = ap_linexpr0_alloc (AP_LINEXPR_DENSE, nbcoeffs);	// nbcoeffs=2
  for (i = 0; i < size; i++) {
    ap_linexpr0_set_coeff_scalar (vertex, i, &scalar);
  }

  for (j = 0; j < v; j++) {
    if (j == 0)
      array.p[j] = ap_generator0_make (AP_GEN_VERTEX, vertex);
    else
      array.p[j] = ap_generator0_copy (&array.p[j - 1]);
    for (i = 0; i < size; i++) {
      ap_scalar_set_double (&scalar, cordinates[i][j]);
      ap_linexpr0_set_coeff_scalar (array.p[j].linexpr0, i, &scalar);

    }
//ap_linexpr0_print(array.p[j].linexpr0, 0x0);  

  }
/*
  for (i=0; i<size; i++){
      for (j=0; j<v; j++){
        if (j>0){
           array.p[j] =ap_generator0_copy(&array.p[j-1]);}
	ap_scalar_set_double(&scalar,cordinates[i][j]);
	ap_linexpr0_set_coeff_scalar(array.p[j].linexpr0,
				     i,&scalar);
        scanf("%f",&b);
        ap_linexpr0_print(array.p[j].linexpr0, 0x0);
        scanf("%f",&b);
      }
    }*/
  /*for(j=0; j < v; j++)
     {
     ap_linexpr0_print(array.p[j].linexpr0, 0x0);
     printf("\n");
     }
     scanf("%f",&b); */
  /* Clear things */
  ap_scalar_clear (&scalar);
  num_clear (num);
  for (i = 0; i < 2 * zon.totnocols; i++){
      free(sgns[i]);}
  free (sgns);
  for (i = 0; i < 2 * zon.totnocols; i++){
      free(negmat[i]);}
  free (negmat);
  for (i = 0; i < zon.row; i++)
    free(cordinates[i]);
  free (cordinates);
  free (zon.center);
      for (i = 0; i < a->dims; i++) {
    for (j = 0; j < pr->dim; j++) {
      itv_clear (ptitv1[i][j]);}}
      for (i = 0; i < a->dims; i++) {
          free(ptitv1[i]);}
      free(ptitv1);
      for (i = 0; i < a->dims; i++) {
    for (j = 0; j < sizegen; j++) {
      itv_clear (ptitv2[i][j]);}}
      for (i = 0; i < a->dims; i++) {
          free(ptitv2[i]);}
      free(ptitv2);
    for (i = 0; i < zon.row; i++){
        free(zon.mat[i]);}
    free(zon.mat);
  return array;


}

tbool_t
t1p_is_intersect (ap_manager_t * man, t1p_t * a1, t1p_t * a2)
{

  bool res;
  if(t1p_is_top (man, a1)){
     res=true;
     return res;}
  if(t1p_is_top (man, a2)){
    res=true;
    return res;}

  CALL ();
  t1p_internal_t *prbx1 = t1p_init_from_manager (man, AP_FUNID_IS_INTERSECT);

  /*int c;
  #ifdef _T1P_DEBUG
    t1p_fprint(stdout, man, a1, 0x0);
    t1p_fprint(stdout, man, a2, 0x0);
    fprintf(stdout, "### ### ###\n");
  #endif
  scanf("%f",&c);*/

  //t1p_aff_update(prbx1, a1);
  //t1p_aff_update(prbx1, a2);

  int i, j, status;
  double numfx1_glp, numfx2_glp, numfx3_glp;
  itv_t coeff;
  itv_init (coeff);
  num_t num;
  num_init (num);
  num_set_double (num, 0);
  itv_set_num (coeff, num);
  num_t num1;
  num_init (num1);
  itv_t *diffc = (itv_t *) calloc (3, sizeof (itv_t));   /// 3? Because, the dimension is 3
  for (i = 0; i < a1->dims; i++) {
    itv_init (diffc[i]);
    itv_sub (diffc[i], a2->paf[i]->c, a1->paf[i]->c);
  }

  itv_magnitude (num1, diffc[0]);
  double_set_num (&numfx1_glp, num1);
  if (itv_is_neg (diffc[0]) == 1)
    numfx1_glp = -1 * numfx1_glp;
  itv_magnitude (num1, diffc[1]);
  double_set_num (&numfx2_glp, num1);
  if (itv_is_neg (diffc[1]) == 1)
    numfx2_glp = -1 * numfx2_glp;
  itv_magnitude (num1, diffc[2]);
  double_set_num (&numfx3_glp, num1);
  if (itv_is_neg (diffc[2]) == 1)
    numfx3_glp = -1 * numfx3_glp;

  glp_prob *lp;
s1:lp = glp_create_prob ();
s2:glp_set_obj_dir (lp, GLP_MIN);
s3:glp_add_rows (lp, 3);
s4:glp_set_row_bnds (lp, 1, GLP_FX, numfx1_glp, numfx1_glp);
s5:glp_set_row_bnds (lp, 2, GLP_FX, numfx2_glp, numfx2_glp);
s6:glp_set_row_bnds (lp, 3, GLP_FX, numfx3_glp, numfx3_glp);

  //printf("%f", numfx1_glp);
  //printf("%f", numfx2_glp);

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

  t1p_aff_get_generator_matrix (prbx1, a1, ptitv1, ptitv2);

  /*for (i=0;i<a1->dims;i++){
   for (j=0;j<k;j++){
       itv_init(ptitv2[i][j]);//=coeff1;
       //itv_print(ptitv[i][j]);
     }
}*/

t1p_aaterm_t *p1 = NULL;

for(i=0;i<a1->dims;i++) {
 for(p1=a1->paf[i]->q;p1;p1=p1->n) {
		//itv_init(ptitv2[i][vec[p1->pnsym->index]]);
                itv_set(ptitv2[i][vec[p1->pnsym->index]],p1->coeff);
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

itv_t tmp;

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

  //printf("%d", size);
  //scanf("%f",&c);

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

  //scanf("%f",&c);

s7:glp_add_cols (lp, size);
  for (j = 0; j < size; j++) {
  s8:glp_set_col_bnds (lp, j + 1, GLP_DB, -1.0, 1.0);
  }

  int ia[1 + 1000], ja[1 + 1000];
  double *ar;
  ar = calloc (1 + 1000, sizeof (double));
  for (i = 0; i < a1->dims; i++) {
    for (j = 0; j < size; j++) {
      itv_magnitude (num1, ptitv_meet[i][j]);
      double_set_num (&numfx1_glp, num1);
      if (itv_is_neg (ptitv_meet[i][j]) == 1) {
	numfx1_glp = -1 * numfx1_glp;
      }
      if (i == 0) {
      s9:ia[1 + j] = i + 1, ja[1 + j] = j + 1, ar[1 + j] = numfx1_glp;
      }
      //printf("%f", ar[1+j]);}
      if (i == 1) {
      s10:ia[1 + j + size] = i + 1, ja[1 + j + size] = j + 1, ar[1 + j + size] =
	  numfx1_glp;
      }
      //printf("%f", ar[1+j+size]);}
    }
  }

  ////// 2*size because the dimension is 2? doubtful, will look into it ///////////////
s11:glp_load_matrix (lp, 2 * size, ia, ja, ar);
s12:glp_simplex (lp, NULL);
s13:status = glp_get_status (lp);

  //printf("%u", status);
  //scanf("%f",&c);

  if (status == 2 || status == 5) {
    res = true;
    itv_clear (coeff);
    num_clear (num);
    num_clear (num1);
    num_clear (num4);
    for (i = 0; i < a1->dims; i++) {
	itv_clear (diffc[i]);}
    free (diffc);
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
  free(vec);
  free(ar);
    return res;
  }
  else {
    res = false;
    itv_clear (coeff);
    num_clear (num);
    num_clear (num1);
    num_clear (num4);
    for (i = 0; i < a1->dims; i++) {
	itv_clear (diffc[i]);}
    free (diffc);
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
  free(vec);
  free(ar);
    return res;
  }

}
