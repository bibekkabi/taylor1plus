#ifndef TILINGS_H_   /* Include guard */
#define TILINGS_H_


struct zonotope
{
    int row;
    int column;
    int totnocols;
    double **mat;
    double *center;
    int nbtiles;
    
};

struct match
{
    int chck;
    int matchcolumn;
};

void signsvec(struct zonotope *zon1, int** sgns);

void signsvecnew(struct zonotope *zon1, int** sgns);

void horzcat (struct zonotope *zon1, double** negmat); 

void circshift1 (int* ind, int* indnew, int count1, int shift);

void vertcat (int** sgns, int** negsgns, struct zonotope *zon);

void find (int sgns[][100], struct match *matchindexcol);  /// not being currently used

void dichoto(struct zonotope *zon1, struct match *matchindexcol);  /// not being currently used 

//int noofcolumns(struct zonotope *zon);

void genhorzcat (struct zonotope *p1, struct match *matchindexcol, int sgns[][100]); /// not being currently used

//void combsoneminusone (int len_temp, int combs[]);



#endif // TILINGS_H_
