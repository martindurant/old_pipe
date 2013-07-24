/*
Compile with:  > gcc image_register2.c util.c analyze_io.c splines.c -lm -o output_executable
=====================================================================

  File:  image_register.c
  Author:        D C Barber Date:        15/01/04

	Program: image_register

	  Usage: image_register  -L LAMDA

		=====================================================================
		Public Routines:  medlab_interface  matrix_interface
		=====================================================================

		  History:  This revision  15/01/04


This version is now a fixed revision and should be archived as such.

*/

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "util.h"

/* ===================  global variables  ======================== */

int MSLICES;	/* number of slices in the images */
int MROWS;	/* number of rows in the image */
int MCOLS;	/* number of columns in the images*/
int CMIN, CMAX, RMIN, RMAX, SMIN, SMAX; /* limits of region of interest */
int XP, YP, ZP;	/* dimensions of node grid */
int NODES_PER_ELEMENT;	/* either 4 (2D) or 8 (3D)*/
int ELEMENTS_PER_NODE;	/* either 4	(2D) or 8 (£D) */
int **ELEMENT_NODES;	/* defines which nodes define each element */
int **NODE_ELEMENTS;	/* defines which elements are connected to each node */
float *XC, *YC, *ZC;	/* co-ordinate arrays for the grid */
/* XC, YC, ZC store the pixel locations of the nodes. i.e. node i corresponds
 * to the location (XC[i],YC[i],ZC[i]) */
int NUMBER_OF_NODES;	/* number of nodes in full grid*/
int NUMBER_OF_ELEMENTS;	/* number of elements in full grid */
int NN;		/* number of nodes in active grid */
int **N;	/* indices of non-zero elements of the P array */
int *NP;	/* number of non-zero elements on each row of P = T'T (and N) */;
double **P;	/* sparse array for storing the P matrix */
double *D;	/* diagonal values of P + lamda*Q */
double *QD;	/* diagonal elements of Q = L'L */
float **Q;	/* sparse array for Q */
int SIDE;	/* spacing between nodes in units of pixel dimensions: chosen by g flag */
int FSM;	/* number of smooths for fixed image: chosem by s parameter */
int MSM;	/* number of smooths for moved image: chosen by m parameter */
float *UT, *VT, *WT, *AT;	/* arrays for storing total displacement values */
float *MAP;
float *RES;
int *PMAP;
float CLIMIT;
int REDUCE;
int ITERM;
int SETLAM;
int PACK;
int S;
float Z;
int FAST;
int m;
int NN;
float ***MOVED;
float ***FIXED;
float ***ROI;
float ***ITEMP;
float LAMDA;
float BI;
float *BL;
float ALAM;
int ND;
int QC;
int QO;
int *NO;
int QS;
int QL;
int QR;
int MR;
float K;
float B;
float BK;
double LINTOL;
float **TDH;
float MIMAX;
int MUTUAL;

double *RANDOM;
double *TEMPA;
double *TEMPB;
double *TEMPC;
double *TEMPD;
double *LAPB;
double *LAPC;

float ORIGSEP;
float SLICESEP;


float ***new_image(int dimslice,int dimrow,int dimcol)
{
	float ***ppp;
	int row, slice;
	float *p;

	p=(float *)calloc((size_t)(dimslice*dimrow*dimcol),(size_t)sizeof(float));
	if(!p)
		return(NULL);

/* allocate pointer to whole buffer */

	ppp=(float ***)calloc((size_t)dimslice,(size_t)sizeof(float **));
	if(!ppp)
	{
		return(NULL);
	}
/* allocate pointers to slices */

	for (slice = 0; slice < dimslice; ++slice)
	{
		ppp[slice] = (float **)calloc((size_t)dimrow,(size_t)sizeof(float *));
		if(!ppp[slice])
			return(NULL);

/* allocate pointers to rows */

		for (row = 0; row < dimrow; ++row)
		{
			ppp[slice][row] = p+row*dimcol+slice*dimrow*dimcol;
			if(!ppp[slice][row])
				return(NULL);
		}
	}
	return(ppp);
}


image *fixed,*moving,*out,*mask;

void make_fixed(float *data, int rows,int cols,int slices, float sepx,float sepy,float sepz)
{
    extern image *fixed;
    fixed = (image *)malloc(sizeof(image));
    fixed->rows = rows;
    fixed->cols = cols;
    fixed->slices = slices;
    fixed->sepx = sepx;
    fixed->sepy = sepy;
    fixed->sepz = sepz;
    fixed->data = data;
 }

void make_moving(float *data, int rows,int cols,int slices, float sepx,float sepy,float sepz)
{
    extern image *moving;
    moving = (image *)malloc(sizeof(image));
    moving->rows = rows;
    moving->cols = cols;
    moving->slices = slices;
    moving->sepx = sepx;
    moving->sepy = sepy;
    moving->sepz = sepz;
    moving->data = data;
}

void make_mask(float *data, int rows,int cols,int slices, float sepx,float sepy,float sepz)
{
    extern image *mask;
    mask = (image *)malloc(sizeof(image));
    mask->rows = rows;
    mask->cols = cols;
    mask->slices = slices;
    mask->sepx = sepx;
    mask->sepy = sepy;
    mask->sepz = sepz;
    mask->data = data;
}

float ***load_image_alt(int *cols, int *rows, int *slices, int choice)
{
	extern float ORIGSEP;
	extern float SLICESEP;
    extern image *fixed,*moving, *mask;
	float ***p;
	int i,j,k;
	int interp=0;
	image *interp_image,*t_image;
printf("Choice %d\n",choice);
	switch (choice){
        case 1: t_image = fixed; break;
        case 2: t_image = moving; break;
        case 3: t_image = mask;
	}

	ORIGSEP = t_image->sepz;

//	if (t_image->sepz != t_image->sepx && t_image->slices > 1) {
//			interp_image = interpolate_slices(t_image, t_image->sepx);
//			t_image = interp_image;
//			interp=1;
//	}
	printf("interped\n");
	SLICESEP = t_image->sepz;

	*cols = t_image->cols;
	*rows = t_image->rows;
	*slices = t_image->slices;

	p = new_image(*slices,*rows,*cols);
	printf("newed\n");
 	for (i = 0; i<*slices; i++)
			for (j=0; j<*rows; j++)
					for (k=0; k<*cols; k++)
							p[i][j][k] = t_image->data[i*t_image->rows*t_image->cols + j*t_image->cols + k];
	if(interp) free_image(t_image);
	return(p);
}

/* These are some useful initialisation routines derived from similar routines in Numerical Recipies.
They are essentially routines for setting up matrices and multi-dimensional arrays.
*/

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define posneg(a) ((a >= 0.0 ? 1: 0))

float *vector(int nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
	float *v;

	v=(float *)calloc((size_t)(nh+1), (size_t)sizeof(float));
	return v;
}

int *ivector(int nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)calloc((size_t)(nh+1), (size_t)sizeof(int));
	return v;
}


double *dvector(long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)calloc((size_t)(nh+1), (size_t)sizeof(double));
	return v;
}

float **matrix(int nrh, int nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	int i, nrow=nrh,ncol=nch;
	float **m;

	/* allocate pointers to rows */
	m=(float **) calloc((size_t)(nrow+1), (size_t)sizeof(float*));

	/* allocate rows and set pointers to them */
	m[1]=(float *) calloc((size_t)(nrow*ncol+1), (size_t)sizeof(float));

	for(i=2;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

double **dmatrix(long nrh, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	int i, nrow=nrh,ncol=nch;
	double **m;

	/* allocate pointers to rows */
	m=(double **) calloc((size_t)(nrow+1), (size_t)sizeof(double*));

	/* allocate rows and set pointers to them */
	m[1]=(double *) calloc((size_t)(nrow*ncol+1), (size_t)sizeof(double));

	for(i=2;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

int **imatrix(long nrh, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	int i, nrow=nrh,ncol=nch;
	int **m;

	/* allocate pointers to rows */
	m=(int **) calloc((size_t)(nrow+1), (size_t)sizeof(int*));

	/* allocate rows and set pointers to them */
	m[1]=(int *) calloc((size_t)(nrow*ncol+1), (size_t)sizeof(int));

	for(i=2;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

void free_matrix(float **m)
/* free a float matrix allocated by matrix() */
{
	free(m[1]);
	free(m);
}

void free_dmatrix(double **m)
/* free a double matrix allocated by dmatrix() */
{
	free(m[1]);
	free(m);
}

void free_imatrix(int **m)
/* free an int matrix allocated by imatrix() */
{
	free(m[1]);
	free(m);
}

int **iarray2(int dimcol,int dimrow)
{
	int **pp;
	int row;
	int *p;

	p=(int *)malloc(dimrow*dimcol*sizeof(int));

	if(!p)
		return(0);

/* allocate pointer to whole buffer */

	pp=(int **)malloc(dimrow*sizeof(int *));
	if(!pp)
		return(NULL);

/* allocate pointers to rows */

	for (row = 0; row < dimrow; ++row)
	{
		pp[row] = p+row*dimcol;
			if(!pp[row])
				return(NULL);

	}
	return(pp);
}

float **array2(dimcol,dimrow)

int dimcol, dimrow;

{
	float **pp;
	int row;
	float *p;

	p=(float *)malloc(dimrow*dimcol*sizeof(float));

	if(!p)
		return(0);

/* allocate pointer to whole buffer */

	pp=(float **)malloc(dimrow*sizeof(float *));
	if(!pp)
	{
		return(NULL);
	}
/* allocate pointers to rows */

	for (row = 0; row < dimrow; ++row)
	{
		pp[row] = p+row*dimcol;
		if(!pp[row])
			return(NULL);

	}
	return(pp);
}

double **darray2(dimcol,dimrow)

int dimcol, dimrow;

{
	double **pp;
	int row;
	double *p;

	p=(double *)malloc(dimrow*dimcol*sizeof(double));

	if(!p)
		return(0);

/* allocate pointer to whole buffer */

	pp=(double **)malloc(dimrow*sizeof(double *));
	if(!pp)
	{
		return(NULL);
	}
/* allocate pointers to rows */

	for (row = 0; row < dimrow; ++row)
	{
		pp[row] = p+row*dimcol;
		if(!pp[row])
			return(NULL);

	}
	return(pp);
}


void free_array_2d(float **pp)
{
	float *p;
	p = pp[0];
	free(p);
	free(pp);
}


void free_double_array_2d(double **pp)
{
	double *p;
	p = pp[0];
	free(p);
	free(pp);
}


double *make_double_array_1d(dimcol)
int dimcol;
{
    double *p;

	p=(double *)malloc(dimcol*sizeof(double));
	if(!p)
		return(0);
	return(p);
}

void free_double_array_1d(double *p)
{
	free(p);
}



char *copy(char *a)
{
	char *c;
	char *p;

	c = (char *)calloc((strlen(a)+1), sizeof(char));
	p = c;

	while((*p++ = *a++))
		;

	return(c);
}



double snrm(int n, double sx[])
/* Returns L_2 norm of array sx */
{
	int i;
	double ans;

	ans = 0.0;
	for (i=1;i<=n;i++) ans += sx[i]*sx[i];
	return sqrt(ans);
}
/* (C) Copr. 1986-92 Numerical Recipes Software $#">A2. */

void delete_image(float ***ppp)
{
	float *p;

	p = ppp[0][0];
	free(p);
	free(ppp);
}

float *array(int dimcol)
{
 	float *p;

	p=(float *)calloc((size_t)dimcol,(size_t)sizeof(float));
	if(!p)
		return(0);

	return(p);
}

void free_array(float *p)
{
	free(p);
}

int *iarray(int dimcol)
{
 	int *p;

	p=(int *)calloc((size_t)dimcol, (size_t)sizeof(int));
	if(!p)
		return(0);

	return(p);
}

void free_iarray(int *p)
{
	free(p);
}


/* ====================  functions ================================*/
#include <math.h>

/* ***************************************************************
 *
 * int pnt(int i, int j, int k)
 *
 * Input: int i,j,k corresponding to node coordinates.
 *
 * Output: Number of node corresponding to grid coordinate (i,j,k).
 *
 * This function appears to be the inverse of get_cords()
 *
 * **************************************************************/
int pnt(int i, int j, int k)
{
	if (i >= 0 && i < XP && j >= 0 && j < YP && k >= 0 && k < ZP)
		return(i + j*XP + k*XP*YP);
	else
		return(-1);
}

void lapsqm(double *a, double *b)
{
	int i,j,k;
	int np;
	int *PN;
	float *PQ;

	extern int NN;
	extern float **Q;
	extern int **N;
	extern int *NP;
	extern double *LAPB;
	extern double *LAPC;
	extern int ND;
	extern float *BL;


	for (k = 0; k < ND; k++)
	{
		np = NN*k;

		for (j = 1; j <= NN; j++)
			LAPC[j] = a[j+np]*BL[j+np];

		for (j = 1; j <= NN; j++)
		{
			LAPB[j] = 0.0;
			PN = &N[j][1];
			PQ = &Q[j][1];
			for (i = 1; i <= NP[j]/ND; i++)
			{
				LAPB[j] += LAPC[*PN]*(*PQ);
				PN++;
				PQ++;
			}

		}

		for (j = 1; j <= NN; j++)
		{
			LAPC[j] = 0.0;
			PN = &N[j][1];
			PQ = &Q[j][1];
			for (i = 1; i <= NP[j]/ND; i++)
			{
				LAPC[j] += LAPB[*PN]*(*PQ);
				PN++;
				PQ++;
			}

		}
		for (j = 1; j <= NN; j++)
			b[j+np] = LAPC[j]*BL[j+np];

	}

}

void lapsq(float *a, float *b, int k)
{
	int i,j;
	int p;
	int np;

	extern int NN;
	extern float **Q;
	extern int **N;
	extern double *TEMPB;
	extern double *TEMPC;
	extern int ND;
	extern float *BL;
	extern int *NP;

	np = k*NN;
	for (j = 1; j <= NN; j++)
		TEMPC[j] = a[j-1]*BL[j+np];

	for (j = 1; j <= NN; j++)
	{
		TEMPB[j] = 0.0;
		for (i = 1; i <= NP[j]/ND; i++)
		{
				TEMPB[j] += TEMPC[N[j][i]]*Q[j][i];
		}

	}

	for (j = 1; j <= NN; j++)
	{
		TEMPC[j] = 0.0;
		for (i = 1; i <= NP[j]/ND; i++)
		{
			p = N[j][i];
			if (p > 0)
				TEMPC[j] += TEMPB[N[j][i]]*Q[j][i];
		}

	}
	for (j = 1; j <= NN; j++)
		b[j-1] = TEMPC[j]*BL[j+np];
}


/***********************************************************
 *
 * void lap_line(int ip, int imax, int jmax, int kmax, int *p)
 *
 * Input: int ip
 * 		int imax,jmax,kmax
 *
 * 	Output: Fills in array p (must be allocated elsewhere)
 *
 * 	Uses: ND -- Augmented dimension of registration (i.e. 4 if 3D, or 3 if
 * 		2D.)
 *
 * 	Notes: In all calling instances, imax, jmax, and kmax are the numbers
 * 	of columns, rows, and slices respectively in the node grid. Using this
 * 	assumption, it would appear that ip corresponds to a "node number".
 * 	It looks like p contains the node numbers of the 3x3 (or 3x3x3) node square
 * 	(or cube) centred at the node given by ip.
 *
 * 	********************************************************/

void lap_line(int ip, int imax, int jmax, int kmax, int *p)
{

	int x,y,z;
	int i,j,k;
	int ipp;

	extern int ND;

	if (ND == 4)
	{
		ipp = ip - 1;
		k = ipp/(imax*jmax);
		ipp = ipp - (k*imax*jmax);
		j = ipp/(imax);
		i = ipp - j*imax;
		ipp = 0;
		for (z = k-1; z <= k+1; z++)
			for (y = j-1; y <= j+1; y++)
				for (x = i-1; x <= i+1; x++)
				{
					if (x >= 0 && x < imax && y >= 0 && y < jmax && z >= 0 && z < kmax)
						p[ipp] = x + (y + z*jmax)*imax +1;
					else
						p[ipp] = 0;
					ipp++;
				}
	}
	else
	{
		ipp = ip - 1;
		j = ipp/(imax);
		i = ipp - j*imax;
		ipp = 0;
		for (y = j-1; y <= j+1; y++)
			for (x = i-1; x <= i+1; x++)
			{
				if (x >= 0 && x < imax && y >= 0 && y < jmax )
					p[ipp] = x + y*imax +1;
				else
					p[ipp] = 0;
				ipp++;
			}
	}
}

void make_pointers()
{
	int i,k,j;
	int pp[27];
	int ppp[108];
	int n, m;

	extern int **N;
	extern int *NP;
	extern int XP, YP, ZP;
	extern int NUMBER_OF_NODES;
	extern int ND;
	extern int QC;
	extern int QL;
	extern int QO;
	extern int *NO;

	for (i = 1; i <= NUMBER_OF_NODES; i++)
	{
		lap_line(i, XP, YP, ZP, pp);
		for (k = 0; k < QC; k++)
		{
			ppp[k] = pp[k];
			if (pp[k] != 0)
				for (j = 1; j < ND; j++)
					ppp[k+j*QC] = pp[k] + j*NUMBER_OF_NODES;
			else
				for (j = 1; j < ND; j++)
					ppp[k+j*QC] = 0;

		}
		ppp[QO-1] = -ppp[QO-1];

		n = 0;
		for (k = 0; k < QL; k++)
			if (ppp[k] != 0)
				ppp[n++] = ppp[k];

		for (k = n; k < QL; k++)
			ppp[k] = 0;

		NP[i] = n;


		m = 0;
		for (j = 0; j <  NP[i]; j++)
			if (ppp[j] < 0)
				m = j+1;

		NO[i] = m;
		ppp[m-1] = -ppp[m-1];


		for (j = 1; j <= NP[i]; j++)
			for (k = 0; k < ND; k++)
				N[k*NUMBER_OF_NODES+i][j] = ppp[j-1];



	}
	for (k = 0; k < ND; k++)
		for (j = 1; j <= NUMBER_OF_NODES; j++)
		{
			NP[j + k*NUMBER_OF_NODES] = NP[j];
			NO[j + k*NUMBER_OF_NODES] = NO[j] + k*(NP[j]/ND);
		}
}

void asolve(int n, double b[], double x[])
{
	extern double *D;
	int i;

	for(i=1;i<=n;i++) x[i]=b[i]/D[i];

}
/* (C) Copr. 1986-92 Numerical Recipes Software $#">A2. */

/**
 * Outputs z = \lambda^2 Q^TQx, and
 * r = (P^TP + \lambda^2 Q^TQ + BK*I) x
 */
void atimes(int n, double x[], double r[], double z[])
{
	int k, i;
	double *pp;
	int *pn;

	extern double **P;
	extern int **N;
	extern int *NP;
	extern float BK;

	for (k = 1;k <= n; k++)
	{
		r[k] = 0.0;
		pn = &N[k][1];
		pp = &P[k][1];
		for (i = 1; i <= NP[k]; i++)
		{
			r[k] += (*pp)*x[*pn];
			pp++;
			pn++;
		}
	}

	lapsqm(x, z);

	for (k = 1; k <= n; k++)
	{
		r[k] += (z[k]+BK*x[k]);
	}


}
/* (C) Copr. 1986-92 Numerical Recipes Software $#">A2. */

int linbcg(int n, double b[], double x[], double tol, int itmax, int *iter, double *err)
{
	void asolve(int n, double b[], double x[]);
	void atimes(int n, double x[], double r[], double temp[]);
	double snrm(int n, double sx[]);
	int j;
	double ak,akden,bk,bkden,bknum,bnrm;
	double *p,*r,*z, *temp;

	p=dvector(n);
	r=dvector(n);
	z=dvector(n);
	temp=dvector(n);

	*iter=0;
	atimes(n,x,r,temp);
	for (j=1;j<=n;j++) {
		r[j]=b[j]-r[j];
	}
	bnrm=snrm(n,b);


	asolve(n,r,z);

	while (*iter <= itmax)
	{
		++(*iter);
		for (bknum=0.0,j=1;j<=n;j++) bknum += z[j]*r[j];

		if (*iter == 1) {
			for (j=1;j<=n;j++) {
				p[j]=z[j];
			}
		}
		else {
			bk=bknum/bkden;
			for (j=1;j<=n;j++) {
				p[j]=bk*p[j]+z[j];
			}
		}
		bkden=bknum;
		atimes(n,p,z,temp);

		for (akden=0.0,j=1;j<=n;j++) akden += z[j]*p[j];
		ak=bknum/akden;
		for (j=1;j<=n;j++) {
			x[j] += ak*p[j];
			r[j] -= ak*z[j];
		}
		asolve(n,r,z);
		*err=snrm(n,r)/bnrm;
		if (*err <= tol) break;
	}

	free(p);
	free(r);
	free(z);
	free(temp);

	return(*iter);
}

void first_comp()
{
	int k, j, i;
	double sa;

	extern int NN;
	extern double *TEMPA;
	extern float LAMDA;
	extern double **P;
	extern double *D;
	extern float ALAM;
	extern double *QD;
	extern double *TEMPB;
	extern double *TEMPC;
	extern int ND;
	extern int QO;
	extern int QC;
	extern int QL;
	extern float *BL;
	extern float K;
	extern float B;
	extern float BK;

	for (k = 1; k <= NN*(ND-1); k++)
		BL[k] = ALAM*sqrt(LAMDA);

	for (k = NN*(ND-1)+1; k <= NN*ND; k++)
		BL[k] = ALAM*sqrt(K*LAMDA);

	for (j = 0; j < ND; j++)
		for (k = 1; k <= NN; k++)
			D[k+j*NN] = P[k+j*NN][QO+j*QC] + BL[k+j*NN]*BL[k+j*NN]*QD[k+j*NN]+BK;

	for (j = 1; j <= NN*ND; j++)
	{
		TEMPC[j] = 0.0;
		for (i = 1; i <= QL; i++)
			TEMPC[j] += P[j][i];
	}

	for (k = 0; k < 20;k++)
	{
		atimes(NN*ND,TEMPC,TEMPB,TEMPA);
		sa = snrm(NN*ND,TEMPB);
		for (j = 1;j <= NN*ND;j++)
			TEMPC[j] = TEMPB[j]/sa;
	}
	sa = snrm(NN*ND,TEMPB);
	BK = sa*B;
}

float cond(float x)
{
	int k, j, i;
	double sa, sb;
	double err;
	float t;
	int iter;

	extern int NN;
	extern int NN;
	extern double *TEMPA;
	extern double *RANDOM;
	extern float LAMDA;
	extern double **P;
	extern double *D;
	extern float ALAM;
	extern double *QD;
	extern double *TEMPB;
	extern double *TEMPC;
	extern float *BL;
	extern int ND;
	extern int QL;
	extern float *BL;
	extern float K;
	extern float B;
	extern float BK;
	extern double LINTOL;
	extern int *NO;

	LAMDA = x;
	BK = 0.0;

	for (k = 1; k <= NN*(ND-1); k++)
		BL[k] = ALAM*sqrt(LAMDA);

	for (k = NN*(ND-1)+1; k <= NN*ND; k++)
		BL[k] = ALAM*sqrt(K*LAMDA);

	for (j = 0; j < ND; j++)
		for (k = 1; k <= NN; k++)
		{
			D[k+j*NN] = P[k+j*NN][NO[k+j*NN]] + BL[k+j*NN]*BL[k+j*NN]*QD[k+j*NN] + BK;
		}

	for (j = 1; j <= NN*ND; j++)
	{
		TEMPC[j] = 0.0;
		for (i = 1; i <= QL; i++)
			TEMPC[j] += P[j][i];
	}

	for (k = 0; k < 20;k++)
	{
		atimes(NN*ND,TEMPC,TEMPB,TEMPA);
		sa = snrm(NN*ND,TEMPB);
		for (j = 1;j <= NN*ND;j++)
			TEMPC[j] = TEMPB[j]/sa;
	}
	sa = snrm(NN*ND,TEMPB);

		BK = sa*B;

	for (j = 1; j <= NN*ND; j++)
	{
		TEMPC[j] = 0.0;
		for (i = 1; i <= QL; i++)
			TEMPC[j] += P[j][i];
	}

	for (k = 0; k < 20;k++)
	{
		atimes(NN*ND,TEMPC,TEMPB,TEMPA);
		sa = snrm(NN*ND,TEMPB);
		for (j = 1;j <= NN*ND;j++)
			TEMPC[j] = TEMPB[j]/sa;
	}
	sa = snrm(NN*ND,TEMPB);

	for (k = 1; k <= NN*ND; k++)
	{
		TEMPB[k] = 0;
		TEMPC[k] = RANDOM[k]/0.5774;
	}


	for (j = 0; j < 3; j++)
	{
		linbcg(NN*ND, TEMPC, TEMPB, LINTOL, 500, &iter, &err);
		for (k = 1; k <= NN*ND; k++)
		{
			TEMPC[k] = TEMPB[k];
			TEMPB[k] = 0.0;
		}
	}
		sb = 1.0/snrm(NN*ND,TEMPB);


	linbcg(NN*ND, TEMPC, TEMPB, LINTOL, 500, &iter, &err);

	sb = 1.0/snrm(NN*ND,TEMPB);
	sb = sqrt(sqrt(sb));
	t = sa/sb;
	printf("cond = %f %f %f LAMDA = %f  iter = %d\n",t, sa, sb, LAMDA, iter);

	return(t);
}


#define NRANSI
#define ITMAX 15
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

float brent(float ax, float bx, float cx, float tol)
{
	int iter;
	float a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	float e=0.0;

	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=w=v=bx;
	fw=fv=fx=cond(x);
	for (iter=1;iter<=ITMAX;iter++) {
		xm=0.5*(a+b);
		tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
			return(x);
		}
		if (fabs(e) > tol1) {
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if (q > 0.0) p = -p;
			q=fabs(q);
			etemp=e;
			e=d;
			if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
				d=CGOLD*(e=(x >= xm ? a-x : b-x));
			else {
				d=p/q;
				u=x+d;
				if (u-a < tol2 || b-u < tol2)
					d=SIGN(tol1,xm-x);
			}
		} else {
			d=CGOLD*(e=(x >= xm ? a-x : b-x));
		}
		u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
		fu=cond(u);
		if (fu <= fx) {
			if (u >= x) a=x; else b=x;
			SHFT(v,w,x,u)
			SHFT(fv,fw,fx,fu)
		} else {
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			} else if (fu <= fv || v == x || v == w) {
				v=u;
				fv=fu;
			}
		}
	}
	return(x);
}
#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software $#">A2. */

/* *********************************************************
 * void invert(double **a, double **b, int size)
 *
 * input: size*size matrix a
 * ouput: size*size matrix b corresponding to a^{-1}
 *
 * Performs Gaussian elimination to invert input matrix a
 * ********************************************************/

void invert(double **a,double **b, int size)
{
	double **h;
	double *em;
	double am,ainv,firv;

	int i,j,k;

	h = darray2(size*2, size);
	em = make_double_array_1d(size*2);
	for (i = 0; i < size; ++i)
		for(j = 0; j < size; ++j)
			h[i][j] = a[i][j];

	for (i = 0; i < size; ++i)
		for (j = 0; j < size; ++j)
			h[i][j+size] = 0.0;

	for (k = 0; k < size; ++k)
		h[k][k+size] = 1.0;


	for (k = 0; k < size; ++k) {
		am = h[k][k];
		ainv = 1.0/am;
		for (j = 0; j < size*2; ++j)
			h[k][j] = ainv*h[k][j];

		for (i = 0; i < size; ++i)
			if (i != k) {
				firv = (-h[i][k]);
				for (j = 0; j < size*2; ++j)
					em[j] = firv*h[k][j];
				for (j = 0; j < size*2; ++j)
					h[i][j] = h[i][j]+em[j];
			}

	}

	for (j=0; j < size; ++j)
		for (k = 0; k < size; ++k)
			b[j][k] = h[j][k+size];

	free_double_array_2d(h);
	free_double_array_1d(em);
	return;
}

void get_cords(int i, int X, int Y, int Z, int *xp, int *yp, int *zp)
{
	*zp = i/(X*Y);
	*yp = (i - (*zp)*X*Y)/X;
	*xp = i - (*yp)*X - (*zp)*X*Y;
}

/* ***********************************************************
 *
 * void make_mesh_data()
 *
 * Input: int ND -- dimension of image (plus one)
 * 		int MCOLS, MROWS, MSLICES -- size of image
 * 		int SIDE -- blocks in grid are SIDE x SIDE (x SIDE)
 *
 * 	Output:
 * 		int XP, YP, ZP -- dimensions of node grid
 * 		int NUMBER_OF_NODES -- XP*YP*ZP
 * 		int NUMBER_OF_ELEMENTS -- (XP-1)*(YP-1)*(ZP-1)
 * 		ELEMENT_NODES -- 2-D array. For each element, gives the index
 * 			of the 4 (or 8) surrounding nodes.
 * 		NODE_ELEMENTS -- 2-D array. For each node, gives the index
 * 			of the 4 (or 8) adjacent elements.
 *
 * This function is responsible for the initialization of many of the
 * global arrays used by the program, including BL, P, NP, N, Q, D, QD,
 * NO, RANDOM, TEMPA, TEMPB, TEMPC, TEMPD, LAPB, LAPC, XC, YC, ZC, MAP,
 * RES, UT, VT, WT, and AT. Most of these are simply allocated, to be
 * populated elsewhere.
 *
 * **********************************************************/
void make_mesh_data()
{
	float x, y, z;
	int col, row, slice;
	float xmin, xmax;
	float ymin, ymax;
	float zmin, zmax;
	int p;
	int kp, jp, ip;
	int k, j, i, n;
	float step;
	float xp, yp, zp;

	extern int MSLICES;
	extern int MROWS;
	extern int MCOLS;
	extern int XP, YP, ZP;
	extern int NODES_PER_ELEMENT;
	extern int ELEMENTS_PER_NODE;
	extern int **ELEMENT_NODES;
	extern int **NODE_ELEMENTS;
	extern float *XC, *YC, *ZC;
	extern int NUMBER_OF_NODES;
	extern int NUMBER_OF_ELEMENTS;
	extern double **P;
	extern double *D;
	extern double *QD;
	extern int **N;
	extern int *NP;
	extern float **Q;
	extern int SIDE;
	extern float *UT, *VT, *WT, *AT;
	extern float *MAP;
	extern float *RES;

	extern double *TEMPA;
	extern double *TEMPB;
	extern double *TEMPC;
	extern double *TEMPD;
	extern double *RANDOM;
	extern double *LAPB;
	extern double *LAPC;

	extern int ND;
	extern int QC;
	extern int QL;
	extern int QO;
	extern int QR;
	extern int *NO;
	extern float *BL;

	if (ND == 4)
	{
		QL = 108;
		QC = 27;
		QO = 14;
		QR = 32;
		NODES_PER_ELEMENT = 8;
		ELEMENTS_PER_NODE = 8;
	}
	else
	{
		QL = 27;
		QC = 9;
		QO = 5;
		QR = 12;
		NODES_PER_ELEMENT = 4;
		ELEMENTS_PER_NODE = 4;
	}


	xp = ceil((MCOLS/((float)SIDE)))/2.0;
	xmax = xp*SIDE;
	xmin = -xmax;
	XP = 2.0*xp + 1;

	yp = ceil((MROWS/((float)SIDE)))/2.0;
	ymax = yp*SIDE;
	ymin = -ymax;
	YP = 2.0*yp + 1;

	if (ND == 4)
	{
		zp = ceil((MSLICES/((float)SIDE)))/2.0;
		zmax = zp*SIDE;
		zmin = -zmax;
		ZP = 2.0*zp + 1;
		NUMBER_OF_ELEMENTS = (XP-1)*(YP-1)*(ZP-1);
		NUMBER_OF_NODES = XP*YP*ZP;
	}
	else
	{
		ZP = 1;
		NUMBER_OF_ELEMENTS = (XP-1)*(YP-1);
		NUMBER_OF_NODES = XP*YP;
	}

	BL = vector(NUMBER_OF_NODES*ND);
	P = dmatrix(NUMBER_OF_NODES*ND,QL);

	NP = ivector(NUMBER_OF_NODES*ND);
	N = imatrix(NUMBER_OF_NODES*ND,QL);
	Q = matrix(NUMBER_OF_NODES,QC);
	D = dvector(NUMBER_OF_NODES*ND);
	QD =dvector(NUMBER_OF_NODES*ND);
	NO =ivector(NUMBER_OF_NODES*ND);

	TEMPA = dvector(NUMBER_OF_NODES*ND);
	RANDOM = dvector(NUMBER_OF_NODES*ND);
	TEMPB = dvector(NUMBER_OF_NODES*ND);
	TEMPC = dvector(NUMBER_OF_NODES*ND);
	TEMPD = dvector(NUMBER_OF_NODES*ND);
	LAPB = dvector(NUMBER_OF_NODES*ND);
	LAPC = dvector(NUMBER_OF_NODES*ND);

	for (j = 1; j <= NUMBER_OF_NODES*ND; j++)
		for (i = 1; i <= QL; i++)
			N[j][i] = 0;


	for (k = 1; k <= NUMBER_OF_NODES*ND; k++)
		RANDOM[k] = (2.0*((float)rand())/RAND_MAX)-1.0;

	make_pointers();


	XC = array(NUMBER_OF_NODES);
	YC = array(NUMBER_OF_NODES);
	ZC = array(NUMBER_OF_NODES);
	MAP = array(NUMBER_OF_NODES);
	RES = array(NUMBER_OF_NODES);

	UT = array(NUMBER_OF_NODES);
	VT = array(NUMBER_OF_NODES);
	WT = array(NUMBER_OF_NODES);
	AT = array(NUMBER_OF_NODES);

	ELEMENT_NODES = iarray2(NODES_PER_ELEMENT, NUMBER_OF_ELEMENTS);
	NODE_ELEMENTS = iarray2(ELEMENTS_PER_NODE, NUMBER_OF_NODES);


	if (ND == 4)
	{
		p = 0;
		step = SIDE;
		for (z = zmin; z <= zmax+step/2; z = z + step)
			for (y = ymin; y <= ymax+step/2; y = y + step)
				for (x = xmin; x <= xmax+step/2; x = x + step)
				{
					XC[p] = x + MCOLS/2 + 0.5;
					YC[p] = y + MROWS/2 + 0.5;
					ZC[p] = z + MSLICES/2 + 0.5;
					p++;
				}
		p = 0;
		for (kp = 0; kp < ZP-1; kp++)
			for (jp = 0; jp < YP-1; jp++)
				for (ip = 0; ip < XP-1; ip++)
				{
					ELEMENT_NODES[p][0] = pnt(ip, jp, kp);
					ELEMENT_NODES[p][1] = pnt(ip+1, jp, kp);
					ELEMENT_NODES[p][2] = pnt(ip, jp+1, kp);
					ELEMENT_NODES[p][3] = pnt(ip+1, jp+1, kp);
					ELEMENT_NODES[p][4] = pnt(ip, jp, kp+1);
					ELEMENT_NODES[p][5] = pnt(ip+1, jp, kp+1);
					ELEMENT_NODES[p][6] = pnt(ip, jp+1, kp+1);
					ELEMENT_NODES[p][7] = pnt(ip+1, jp+1, kp+1);
					p++;
				}

		for (n = 0; n < NUMBER_OF_NODES; n++)
		{
			get_cords(n, XP, YP, ZP, &ip, &jp, &kp);
			p = 0;
			for (k = -1; k <= 0; k++)
				for (j = -1; j <= 0; j++)
					for (i = -1; i <= 0; i++)
					{
						col = ip + i;
						row = jp + j;
						slice = kp + k;
						if (col < 0)col = 0;
						if (col > XP - 2)col = XP - 2;
						if (row < 0)row = 0;
						if (row > YP - 2)row = YP - 2;
						if (slice < 0)slice = 0;
						if (slice > ZP - 2)slice = ZP - 2;
						NODE_ELEMENTS[n][p] = col + row*(XP-1) + slice*(XP-1)*(YP-1);
						p++;
					}

		}
	}
	else
	{
		p = 0;
		step = SIDE;
		for (y = ymin; y <= ymax+step/2; y = y + step)
			for (x = xmin; x <= xmax+step/2; x = x + step)
			{
				XC[p] = x + MCOLS/2 + 0.5;
				YC[p] = y + MROWS/2 + 0.5;
				ZC[p] = 0.0;
				p++;
			}

		p = 0;
		for (jp = 0; jp < YP-1; jp++)
			for (ip = 0; ip < XP-1; ip++)
			{
				ELEMENT_NODES[p][0] = pnt(ip, jp, 0);
				ELEMENT_NODES[p][1] = pnt(ip+1, jp, 0);
				ELEMENT_NODES[p][2] = pnt(ip, jp+1, 0);
				ELEMENT_NODES[p][3] = pnt(ip+1, jp+1, 0);
				p++;
			}

		for (n = 0; n < NUMBER_OF_NODES; n++)
		{
			get_cords(n, XP, YP, ZP, &ip, &jp, &kp);
			p = 0;
			for (j = -1; j <= 0; j++)
				for (i = -1; i <= 0; i++)
				{
					col = ip + i;
					row = jp + j;
					if (col < 0)col = 0;
					if (col > XP - 2)col = XP - 2;
					if (row < 0)row = 0;
					if (row > YP - 2)row = YP - 2;
					NODE_ELEMENTS[n][p] = col + row*(XP-1);
					p++;
				}
		}
	}
}


void make_q()
{
	int s[27];
	float QT[27];
	int j,i,n;
	float t;
	float *w;
	float wt_3d[] = {1,2,1,2,4,2,1,2,1,2,4,2,4,-56,4,2,4,2,1,2,1,2,4,2,1,2,1};
	float wt_2d[] = {0.5, 1.0, 0.5, 1.0, -6.0, 1.0, 0.5, 1.0, 0.5};


	extern float **Q;
	extern double *QD;
	extern int REDUCE;
	extern int NUMBER_OF_NODES;
	extern int XP, YP, ZP;
	extern int ND;

	extern int QC;
	extern int QO;


	for (j = 1; j <= NUMBER_OF_NODES; j++)
		for (i = 1; i <= QC; i++)
		{
			Q[j][i] = 0.0;
		}

	if (ND == 4)
		w = wt_3d;
	else
		w = wt_2d;

	for (j = 1; j <= NUMBER_OF_NODES; j++)
	{
		lap_line(j, XP, YP, ZP, s);

		for (i = 0; i < QC; i++)
			QT[i] = 0.0;


		t = 0;
		for (i = 0; i < QC; i++)
			if (s[i] != 0)
			{
				QT[i] = w[i];
				t += w[i];
			}

		if (REDUCE == 0)
			QT[QO-1] -= t;

		n = 0;
		for (i = 0; i < QC; i++)
			if (QT[i] != 0)
				QT[n++] = QT[i];

		QD[j] = 0.0;
		for (i = 1; i <= n; i++)
			{
				Q[j][i] = QT[i-1];
				QD[j] += QT[i-1]*QT[i-1];
			}

	}

	for (i = 0; i < ND; i++)
		for (j = 1; j <= NUMBER_OF_NODES; j++)
			QD[j+i*NUMBER_OF_NODES] = QD[j];

}

void update_3d(float *U, float *V, float *W, float *A)
{
	double **X;
	double **XI;
	float UVW[4][16];
	float AA[4][16];
	float **ut;
	int k,n,j,m,i;
	float x, y, z, a;
	float xc, yc, zc, ac;
	float step;
	float Y[16];

	extern int MCOLS, MROWS, MSLICES;
	extern int **ELEMENT_NODES;
	extern int **NODE_ELEMENTS;
	extern int NUMBER_OF_NODES;
	extern float *XC, *YC, *ZC;
	extern int SIDE;
	extern float *UT, *VT, *WT, *AT;
	extern float *RES;


	for (k = 0; k < NUMBER_OF_NODES; k++)
	{
		XC[k] = XC[k] - MCOLS/2 - 0.5;
		YC[k] = YC[k] - MROWS/2 - 0.5;
		ZC[k] = ZC[k] - MSLICES/2 - 0.5;
	}
	step = (float)SIDE;

	X = darray2(16,16);
	XI = darray2(16,16);
	ut = array2(4,NUMBER_OF_NODES);

	for (k = 0; k < NUMBER_OF_NODES;k++)
	{


		n = ((U[k] >= 0.0)?1:0) + 2*((V[k] >= 0.0)?1:0) + 4*((W[k] >= 0.0)?1:0);
		n = NODE_ELEMENTS[k][n];

		for (j = 0; j < 8; j++)
		{
			m = ELEMENT_NODES[n][j];
			xc = XC[m];
			yc = YC[m];
			zc = ZC[m];

			if (A[k] >= 0.0)
				ac = step;
			else
				ac = -step;


			X[0][j] = 1.0;
			X[1][j] = xc;
			X[2][j] = yc;
			X[3][j] = xc*yc;
			X[4][j] = zc;
			X[5][j] = zc*xc;
			X[6][j] = zc*yc;
			X[7][j] = xc*yc*zc;
			X[8][j] = 0.0;
			X[9][j] = 0.0;
			X[10][j] = 0.0;
			X[11][j] = 0.0;
			X[12][j] = 0.0;
			X[13][j] = 0.0;
			X[14][j] = 0.0;
			X[15][j] = 0.0;

			X[0][j+8] = 1.0;
			X[1][j+8] = xc;
			X[2][j+8] = yc;
			X[3][j+8] = xc*yc;
			X[4][j+8] = zc;
			X[5][j+8] = zc*xc;
			X[6][j+8] = zc*yc;
			X[7][j+8] = xc*yc*zc;
			X[8][j+8] = ac;
			X[9][j+8] = ac*xc;
			X[10][j+8] = ac*yc;
			X[11][j+8] = ac*xc*yc;
			X[12][j+8] = ac*zc;
			X[13][j+8] = ac*zc*xc;
			X[14][j+8] = ac*zc*yc;
			X[15][j+8] = ac*xc*yc*zc;
		}

		invert(X,XI,16);

		x = U[k] + XC[k];
		y = V[k] + YC[k];
		z = W[k] + ZC[k];
		a = A[k];

		Y[0] = 1.0;
		Y[1] = x;
		Y[2] = y;
		Y[3] = x*y;
		Y[4] = z;
		Y[5] = z*x;
		Y[6] = z*y;
		Y[7] = x*y*z;
		Y[8] = a;
		Y[9] = a*x;
		Y[10] = a*y;
		Y[11] = a*x*y;
		Y[12] = a*z;
		Y[13] = a*z*x;
		Y[14] = a*z*y;
		Y[15] = a*x*y*z;



		for (j = 0; j < 8; j++)
		{
			m = ELEMENT_NODES[n][j];
			xc = XC[m];
			yc = YC[m];
			zc = ZC[m];

			if (A[k] >= 0.0)
				ac = step;
			else
				ac = -step;

			UVW[0][j] = UT[m]+xc;
			UVW[1][j] = VT[m]+yc;
			UVW[2][j] = WT[m]+zc;
			UVW[3][j] = AT[m];

			UVW[0][j+8] = UT[m]+xc;
			UVW[1][j+8] = VT[m]+yc;
			UVW[2][j+8] = WT[m]+zc;
			UVW[3][j+8] = AT[m]+ac;
		}



		for (j = 0; j < 4; j++)
		{
			for (i = 0; i < 16; i++)
			{
				AA[j][i] = 0.0;
				for (n = 0; n < 16; n++)
					AA[j][i] += UVW[j][n]*XI[n][i];
			}
		}


		for (j = 0; j <4; j++)
		{
			ut[k][j] = 0.0;
			for (i = 0; i < 16; i++)
				ut[k][j] += AA[j][i]*Y[i];
		}



	}

	for (k = 0; k < NUMBER_OF_NODES; k++)
	{

		ut[k][0] -= XC[k];
		ut[k][1] -= YC[k];
		ut[k][2] -= ZC[k];

		RES[k] = sqrt((ut[k][0]-UT[k])*(ut[k][0]-UT[k]) + (ut[k][1]-VT[k])*(ut[k][1]-VT[k])
			+ (ut[k][2]-WT[k])*(ut[k][2]-WT[k]));


		UT[k] = ut[k][0];
		VT[k] = ut[k][1];
		WT[k] = ut[k][2];
		AT[k] = ut[k][3];


		XC[k] = XC[k] + MCOLS/2 + 0.5;
		YC[k] = YC[k] + MROWS/2 + 0.5;
		ZC[k] = ZC[k] + MSLICES/2 + 0.5;
	}

	free_double_array_2d(X);
	free_double_array_2d(XI);


}

float normalise(float ***p, float ***roi, float val)
{
	int col, row, slice;
	float tot;
	long num;

	extern int MCOLS, MROWS, MSLICES;
	extern int CMIN, CMAX, RMIN, RMAX, SMIN, SMAX;

	tot = 0.0;
	num = 0;
	for (slice = SMIN; slice <= SMAX; slice++)
		for (row = RMIN; row <= RMAX; row++)
			for (col = CMIN; col <= CMAX; col++)
				if (roi[slice][row][col] != 0)
				{
					tot += fabs(p[slice][row][col]);
					num++;
				}
	tot = tot/num;
	tot = val/tot;

	for (slice = 0; slice < MSLICES; slice++)
		for (row = 0; row < MROWS; row++)
			for (col = 0; col < MCOLS; col++)
				p[slice][row][col] *= tot;
	return(tot);
}

float trilint(float r, float s, float t, float ***p)

{
    int nr,ns,nt;
    int nru,nsu,ntu;
    float wr,ws,wt;
    float va,vau;
    float vb,vbu;
    float vc,vcu;
    float res;

    extern int MCOLS, MROWS, MSLICES;

    nr = (int)r;
	ns = (int)s;
    nt = (int)t;

	if (r < 0)
	{
		nr = 0;
		r = nr;
	}
	if (nr >= (MCOLS - 2))
	{
		nr = MCOLS -2;
		r = nr+1;
	}

	if (s < 0)
	{
		ns = 0;
		s = ns;
	}
	if (ns >= (MROWS - 2))
	{
		ns = MROWS -2;
		s = ns+1;
	}

	if (t < 0)
	{
		nt = 0;
		t = nt;
	}
	if (nt >= (MSLICES - 2))
	{
		nt = MSLICES -2;
		t = nt+1;
	}



	nru=nr+1;
	nsu=ns+1;
	ntu=nt+1;
	wr=r-nr;
	ws=s-ns;
	wt=t-nt;
	va=p[nt][ns][nr]+(p[nt][ns][nru]-p[nt][ns][nr])*wr;
	vau=p[nt+1][ns][nr]+(p[ntu][ns][nru]-p[ntu][ns][nr])*wr;
	vc=va+(vau-va)*wt;

	vb=p[nt][nsu][nr]+(p[nt][nsu][nru]-p[nt][nsu][nr])*wr;
	vbu=p[ntu][nsu][nr]+(p[ntu][nsu][nru]-p[ntu][nsu][nr])*wr;
	vcu=vb+(vbu-vb)*wt;

	res = vc+(vcu-vc)*ws;

	return (res);


}

void interpolate_image_3d(float ***pic)
{
    int row, col, slice;
    float sa, ra, ta;
	float va, vb, wr, ws, wt, vca, vcb;
    float x,y,z;
	int ip, jp, kp;
	int na, nb, nc, nd, ne, nf, ng, nh;

	extern float *UT, *VT, *WT;
	extern int SIDE;
	extern float *XC, *YC, *ZC;
	extern float ***MOVED;
	extern int XP, YP, ZP;
	extern int MCOLS, MROWS, MSLICES;


	for (slice = 0; slice < MSLICES; slice++)
		for (row = 0; row < MROWS; row++)
			for (col = 0; col < MCOLS; col++)
				pic[slice][row][col] = 0.0;

	for (kp = 0; kp < ZP-1; kp++)
	for (jp = 0; jp < YP-1; jp++)
	for (ip = 0; ip < XP-1; ip++)
	{
		na = pnt(ip, jp, kp);
		nb = pnt(ip+1, jp, kp);
		nc = pnt(ip, jp+1, kp);
		nd = pnt(ip+1, jp+1, kp);
		ne = pnt(ip, jp, kp+1);
		nf = pnt(ip+1, jp, kp+1);
		ng = pnt(ip, jp+1, kp+1);
		nh = pnt(ip+1, jp+1, kp+1);
		x = XC[na] - 0.5;
		y = YC[na] - 0.5;
		z = ZC[na] - 0.5;

		for (slice = z; slice < z+SIDE; slice++)
		for (row = y; row < y+SIDE; row++)
			for (col = x; col < x+SIDE; col++)
				if (col >= 0 && col < MCOLS && row >= 0 && row < MROWS && slice >= 0 && slice < MSLICES)
				{
					wt=(slice-z+0.5)/SIDE;
					ws=(row-y+0.5)/SIDE;
					wr=(col-x+0.5)/SIDE;

					va=UT[na] + (UT[nb] - UT[na])*wr;
					vb=UT[nc] + (UT[nd] - UT[nc])*wr;

					vca = va+(vb-va)*ws;

					va=UT[ne] + (UT[nf] - UT[ne])*wr;
					vb=UT[ng] + (UT[nh] - UT[ng])*wr;

					vcb = va+(vb-va)*ws;

					ra = vca + (vcb - vca)*wt;

					va=VT[na] + (VT[nb] - VT[na])*wr;
					vb=VT[nc] + (VT[nd] - VT[nc])*wr;

					vca = va+(vb-va)*ws;

					va=VT[ne] + (VT[nf] - VT[ne])*wr;
					vb=VT[ng] + (VT[nh] - VT[ng])*wr;

					vcb = va+(vb-va)*ws;

					sa = vca + (vcb - vca)*wt;

					va=WT[na] + (WT[nb] - WT[na])*wr;
					vb=WT[nc] + (WT[nd] - WT[nc])*wr;

					vca = va+(vb-va)*ws;

					va=WT[ne] + (WT[nf] - WT[ne])*wr;
					vb=WT[ng] + (WT[nh] - WT[ng])*wr;

					vcb = va+(vb-va)*ws;

					ta = vca + (vcb - vca)*wt;


					pic[slice][row][col] = trilint(ra + col, sa + row, ta + slice, MOVED);

				}
		}
}


void interpolate_image_amp_3d(float ***pic)
{
    int row, col, slice;
    float aa;
	float va, vb, wr, ws, wt, vca, vcb;
    float x,y,z;
	int ip, jp, kp;
	int na, nb, nc, nd, ne, nf, ng, nh;

	extern float ***MOVED;
	extern int XP, YP, ZP;
	extern float *AT;
	extern int SIDE;
	extern float *XC, *YC, *ZC;
	extern float BI;

	for (kp = 0; kp < ZP-1; kp++)
	for (jp = 0; jp < YP-1; jp++)
	for (ip = 0; ip < XP-1; ip++)
	{
		na = pnt(ip, jp, kp);
		nb = pnt(ip+1, jp, kp);
		nc = pnt(ip, jp+1, kp);
		nd = pnt(ip+1, jp+1, kp);
		ne = pnt(ip, jp, kp+1);
		nf = pnt(ip+1, jp, kp+1);
		ng = pnt(ip, jp+1, kp+1);
		nh = pnt(ip+1, jp+1, kp+1);
		x = XC[na] - 0.5;
		y = YC[na] - 0.5;
		z = ZC[na] - 0.5;

		for (slice = z; slice < z+SIDE; slice++)
		for (row = y; row < y+SIDE; row++)
			for (col = x; col < x+SIDE; col++)
				if (col >= 0 && col < MCOLS && row >= 0 && row < MROWS && slice >= 0 && slice < MSLICES)
				{
					wt=(slice-z+0.5)/SIDE;
					ws=(row-y+0.5)/SIDE;
					wr=(col-x+0.5)/SIDE;


					va=AT[na] + (AT[nb] - AT[na])*wr;
					vb=AT[nc] + (AT[nd] - AT[nc])*wr;

					vca = va+(vb-va)*ws;

					va=AT[ne] + (AT[nf] - AT[ne])*wr;
					vb=AT[ng] + (AT[nh] - AT[ng])*wr;

					vcb = va+(vb-va)*ws;

					aa = vca + (vcb - vca)*wt;

					pic[slice][row][col] -= aa*BI*(FIXED[slice][row][col]+MOVED[slice][row][col])/2.0;

				}
		}


}


void where_in_list(int *N, int *qp, int *pt)
{
	int t, p, k;

	extern int QL;
	extern int QR;

	p = 0;
	t = qp[p];
	for (k = 1; k <= QL; k++)
	{
		if (t == N[k])
		{
			pt[p] = k;
			p++;
			if (p >= QR)return;
			t = qp[p];
		}
	}
	return;
}

float MI(float ***fixed, float***roi, float ***registered)

{
	int row,col,slice;
	float pa[64], pb[64];
 	int i,j,k,m,N;
	float PA, PB, tdh, h;
	float mi;
	float maxf, maxr;
	int max;
	float minf, minr;

	extern int MSLICES, MCOLS, MROWS;
	extern int CMIN, CMAX, RMIN, RMAX, SMIN, SMAX;
 	extern float **TDH;

	max = 64;

	for (j = 0; j < max; j++)
	{
		pa[j] = 0.0;
		pb[j] = 0.0;
		for (i = 0; i < max; i++)
			TDH[i][j] = 0;
	}

	maxf = maxr = minf = minr = 0.0;
	for (slice = SMIN; slice <= SMAX; slice++)
		for (row = RMIN; row <= RMAX; row++)
			for (col = CMIN; col <= CMAX; col++)
				if (roi[slice][row][col] != 0)
				{
					if (fixed[slice][row][col] > maxf)
						maxf = fixed[slice][row][col];
					if (fixed[slice][row][col] < minf)
						minf = fixed[slice][row][col];
					if (registered[slice][row][col] > maxr)
						maxr = registered[slice][row][col];
					if (registered[slice][row][col] < minr)
						minr = registered[slice][row][col];
				}


	for (slice = SMIN; slice <= SMAX; slice++)
		for (row = RMIN; row <= RMAX; row++)
			for (col = CMIN; col <= CMAX; col++)
				if (roi[slice][row][col] != 0)
				{
					k = (int)((registered[slice][row][col]-minr)*62/(maxr-minr));
					m = (int)((fixed[slice][row][col]-minf)*62/(maxf-minf));
					TDH[k][m] += 1;
				}

	N = MROWS*MCOLS*MSLICES;

	for (j = 0; j < max; j++)
		for (i = 0; i < max; i++)
			TDH[j][i] /= N;


	for (j = 0; j < max; j++)
	{
		h = 0.0;
		for (i = 0; i < max; i++)
			h += TDH[j][i];
		pa[j] = h;
	}

	for (i = 0; i < max; i++)
	{
		h = 0.0;
		for (j = 0; j < max; j++)
			h += TDH[j][i];
		pb[i] = h;
	}

	mi = 0.0;

	for (j = 0; j < max; j++)
		for (i = 0; i < max; i++)
		{
			tdh = TDH[j][i];
			PA = pa[j];
			PB = pb[i];
			if (tdh*PA*PB != 0.0)
				mi += tdh*log(tdh/(PA*PB));
		}

	return(mi);
}

void pack_matrices()
{
	int a,b,k,j, kp, jp;
	extern double **P;
	extern double *QD;
	extern int NN;
	extern int **N;
	extern int ND;
	extern int QL;
	extern int *PMAP;


	a = 1;
 	for (k = 1; k <= NUMBER_OF_NODES; k++)
	{
		kp = PMAP[k];
		if (kp > 0)
		{
			b = 1;
			for (j = 1; j <= NP[k]/ND; j++)
			{
				jp = PMAP[N[k][j]];
				if (jp > 0)
				{
					Q[a][b] = Q[k][j];
					b++;
				}
			}
			a++;
		}
	}

	a = 1;
	for (k = 1; k < NUMBER_OF_NODES*ND; k++)
	{
		kp = PMAP[k];
		if (kp > 0)
		{
			QD[a] = QD[kp];
			a++;
		}
	}

	a = 1;
	for (k = 1; k <= NUMBER_OF_NODES*ND; k++)
	{
		kp = PMAP[k];
		if (kp > 0)
		{
			b = 1;
			for (j = 1; j <= NP[k]; j++)
			{
				jp = PMAP[N[k][j]];
				if (jp > 0)
				{
					P[a][b] = P[k][j];
					b++;
				}
			}
			a++;
		}
	}

	a = 1;
	for (k = 1; k <= NUMBER_OF_NODES*ND; k++)
	{
		kp = PMAP[k];
		if (kp > 0)
		{
			b = 1;
			for (j = 1; j <= NP[k]; j++)
			{
				jp = PMAP[N[k][j]];
				if (jp > 0)
				{
					if (k == N[k][j])
						NO[a] = b;
					N[a][b] = jp;
					b++;
				}
			}
				for (j = b; j <= 108; j++)
					N[a][j] = 0;
				a++;
		}
	}
	NN = (a-1)/ND;

	for (k = 1; k <= NN*ND; k++)
	{
		a = 0;
		for (j = 1; j <= QL; j++)
			if (N[k][j] > 0)
				a++;
		NP[k] = a;
	}
}


int register_image_full_3d(float ***fixed, float ***roi, float ***registered)
{
    int row, col, slice;
    extern int m;
    int i, j, k;
    double *TD;
	double *TDX;
    float q[32];
    int qj[32];
	int qi[32][32];
	float PTP[32][32];
	double TDT[32];
    float grad[4];
    float *U, *V, *W, *A;   /* matrices containing the incremental x,y,and z displacements at each node*/
    float um, umm ,ua, nm, uat;
    float *w;
    int *pp;
    int p;
    int a;
	float fmd;
	int na, nb, nc, nd, ne, nf, ng, nh;
	float wr, ws, wt;
	float xoff, yoff, zoff;
	int ip, jp, kp;
	int iter;
	double err;
	float qqj;
	float mi, mit;
	int mc;
	float *UO, *VO, *WO, *AO;
	float *UL, *VL, *WL, *AL;
	float *US, *VS, *WS, *AS;
	int fast;

	extern double **P;
	extern double *D;
	extern double *QD;
	extern int **N;
	extern int ND;
	extern int QL;
	extern int QR;
	extern int MR;
	extern int *PMAP;

    extern int XP, YP, ZP;
    extern int NUMBER_OF_NODES;
	extern int NN;
	extern float *UT, *VT, *WT, *AT;
	extern float *RES;
	extern float LAMDA;
	extern int REDUCE;
	extern float *BL;
	extern float BI;
	extern double LINTOL;
	extern int FAST;
	extern float MIMAX;
	extern int MUTUAL;

    TD = dvector(NUMBER_OF_NODES*ND);
    TDX = dvector(NUMBER_OF_NODES*ND);

    U = array(NUMBER_OF_NODES);
    V = array(NUMBER_OF_NODES);
    W = array(NUMBER_OF_NODES);
    A = array(NUMBER_OF_NODES);

	US = array(NUMBER_OF_NODES);
    VS = array(NUMBER_OF_NODES);
    WS = array(NUMBER_OF_NODES);
    AS = array(NUMBER_OF_NODES);

    UL = array(NN);
    VL = array(NN);
    WL = array(NN);
    AL = array(NN);

    UO = array(NN);
    VO = array(NN);
    WO = array(NN);
    AO = array(NN);

	w = array(NODES_PER_ELEMENT);
    pp = iarray(NODES_PER_ELEMENT);

	m = 0;
	mc = 0;

	for (slice = 0; slice < MSLICES; slice++)
		for (row = 0; row < MROWS; row++)
			for (col = 0; col < MCOLS; col++)
				registered[slice][row][col] = MOVED[slice][row][col];

	fast = 0;
	do
	{
		if (REDUCE == 1)
		{
			make_pointers();
			make_q();
		}

		for (j = 1; j <= NUMBER_OF_NODES*ND; j++)
			TD[j] = 0.0;

		if (fast == 0)
			for (j = 1; j <= NUMBER_OF_NODES*ND; j++)
				for (i = 1; i <= QL; i++)
					P[j][i] = 0.0;

		for (kp = 0; kp < ZP-1; kp++)
			for (jp = 0; jp < YP-1; jp++)
				for (ip = 0; ip < XP-1; ip++)
				{
					na = pnt(ip, jp, kp);
					nb = pnt(ip+1, jp, kp);
					nc = pnt(ip, jp+1, kp);
					nd = pnt(ip+1, jp+1, kp);
					ne = pnt(ip, jp, kp+1);
					nf = pnt(ip+1, jp, kp+1);
					ng = pnt(ip, jp+1, kp+1);
					nh = pnt(ip+1, jp+1, kp+1);
					xoff = XC[na] - 0.5;
					yoff = YC[na] - 0.5;
					zoff = ZC[na] - 0.5;
					pp[0] = na;
					pp[1] = nb;
					pp[2] = nc;
					pp[3] = nd;
					pp[4] = ne;
					pp[5] = nf;
					pp[6] = ng;
					pp[7] = nh;

					a = 0;
					for (i = 0; i < 4; i++)
						for (j = 0; j < 8; j++)
							qj[a++] = pp[j] + i*NUMBER_OF_NODES +1;


					for (j = 0; j < QR; j++)
						where_in_list(N[qj[j]],qj,qi[j]);

					for (j=0; j < 32; j++)
						TDT[j]= 0;

					if (fast == 0)
						for (j = 0; j < 32; j++)
							for (i = 0; i < 32; i++)
								PTP[j][i] = 0.0;

					for (slice = zoff; slice < zoff+SIDE; slice++)
						for (row = yoff; row < yoff+SIDE; row++)
							for (col = xoff; col < xoff+SIDE; col++)
								if (col >= 0 && col < MCOLS && row >= 0 && row < MROWS && slice >= 0 && slice < MSLICES)
								{
									wt=(slice-zoff+0.5)/SIDE;
									ws=(row-yoff+0.5)/SIDE;
									wr=(col-xoff+0.5)/SIDE;
									w[0] = (1-ws)*(1-wr)*(1-wt);
									w[1] = wr*(1-ws)*(1-wt);
									w[2] = ws*(1-wr)*(1-wt);
									w[3] = ws*wr*(1-wt);
									w[4] = (1-ws)*(1-wr)*wt;
									w[5] = wr*(1-ws)*wt;
									w[6] = ws*(1-wr)*wt;
									w[7] = ws*wr*wt;
									if (roi[slice][row][col] >= 1.0)
									{
										grad[0] = (fixed[slice][row][col+1]-fixed[slice][row][col-1]);
										grad[0] += (registered[slice][row][col+1]-registered[slice][row][col-1]);
										grad[0] /= 4.0;
										grad[1] = (fixed[slice][row+1][col]-fixed[slice][row-1][col]);
										grad[1] += (registered[slice][row+1][col]-registered[slice][row-1][col]);
										grad[1] /= 4.0;
										grad[2] = (fixed[slice+1][row][col]-fixed[slice-1][row][col]);
										grad[2] += (registered[slice+1][row][col]-registered[slice-1][row][col]);
										grad[2] /= 4.0;
										grad[3] = -BI*(fixed[slice][row][col] + registered[slice][row][col])/2.0;

										fmd = fixed[slice][row][col] -registered[slice][row][col];

										a = 0;
										for (i = 0; i < 4; i++)
											for (j = 0; j < 8; j++)
												q[a++] = w[j]*grad[i];

										for (j = 0; j < QR; j++)
											TDT[j] += q[j]*fmd;

										if (fast == 0)
											for (j = 0; j < QR; j++)
											{
												qqj = q[j];
												for (i = j; i < QR; i++)
													PTP[j][i] += qqj*q[i];
											}
									}

								}
								for (j = 0; j < QR; j++)
									TD[qj[j]] += TDT[j];

								if (fast == 0)
								{
									for (j = 0; j < QR; j++)
										for (i = j; i < QR; i++)
											PTP[i][j] = PTP[j][i];

									for (j = 0; j < QR; j++)
										for (i = 0; i < QR; i++)
											P[qj[j]][qi[j][i]] += (double)PTP[j][i];
								}
				}

				if (fast == 0)
					NN = NUMBER_OF_NODES;

				if (REDUCE == 1 && fast == 0)
				{
					pack_matrices();
					for (k = 1; k <= NUMBER_OF_NODES*ND; k++)
					{
						kp = PMAP[k];
						if (kp > 0)
							TD[kp] = TD[k];
					}
				}

				for (j = 0; j < ND; j++)
					for (k = 1; k <= NN; k++)
						D[k+j*NN] = P[k+j*NN][NO[k+j*NN]] + BL[k+j*NN]*BL[k+j*NN]*QD[k+j*NN] +BK;


				lapsq(UO,UL,0);
				lapsq(VO,VL,1);
				lapsq(WO,WL,2);
				lapsq(AO,AL,3);
				p = 1;
				for (k = 0; k < NN; k++)
					TD[p++] -= UL[k];
				for (k = 0; k < NN; k++)
					TD[p++] -= VL[k];
				for (k = 0; k < NN; k++)
					TD[p++] -= WL[k];
				for (k = 0; k < NN; k++)
					TD[p++] -= AL[k];

				for (k = 1; k <= NN*ND; k++)
					TDX[k] = 0.0;

				linbcg(NN*ND, TD, TDX, LINTOL, 300, &iter, &err);

				p = 1;

				for (k = 0; k < NN; k++)
					UL[k] = TDX[p++];
				for (k = 0; k < NN; k++)
					VL[k] = TDX[p++];
				for (k = 0; k < NN; k++)
					WL[k] = TDX[p++];
				for (k = 0; k < NN; k++)
					AL[k] = TDX[p++];

				for (k = 1; k <= NUMBER_OF_NODES; k++)
				{
					kp = PMAP[k];
					if (kp > 0)
						U[k-1] = UL[kp-1],V[k-1] = VL[kp-1],W[k-1] = WL[kp-1],A[k-1] = AL[kp-1];
				}
				for (k = 0; k < NUMBER_OF_NODES; k++)
					RES[k] = 0.0,US[k] = UT[k],VS[k] = VT[k],WS[k] = WT[k],AS[k] = AT[k];

				update_3d(U,V,W,A);

				um = ua = nm = 0.0;
				for (k = 0; k < NUMBER_OF_NODES; k++)
				{
					if (MAP[k] > 0.0)
					{
						umm = RES[k];
						if (umm > um)
						{
							um = umm;
							j = k;
						}
						ua += umm*MAP[k];
						nm += MAP[k];
					}
				}

				ua /= nm;
				uat = ua;

				for (k = 1; k <= NUMBER_OF_NODES; k++)
				{
					kp = PMAP[k];
					if (kp > 0)
						UO[kp-1] = UT[k-1],VO[kp-1] = VT[k-1],WO[kp-1] = WT[k-1],AO[kp-1] = AT[k-1];
				}


/* interpolation is always from the original MOVED image using the current total displavcement*/


				interpolate_image_3d(registered);

				mi = MI(fixed,roi,registered);
				mit = mi;

				if (FAST == 1) {
					if (ua < 1.0)
						fast = 1;
					else
						fast = 0;
				}

				if (mit < MIMAX && MUTUAL == 1)
				{
					for (k = 0; k < NUMBER_OF_NODES; k++)
						UT[k] = US[k],VT[k] = VS[k],WT[k] = WS[k],AT[k] = AS[k];

					uat = 0.0;

					interpolate_image_3d(registered);
					interpolate_image_amp_3d(registered);

					if (mc == 0) {
						if (m == 0)
							return(1);
						else
							return(0);
					}
				}
				else
				{
					MIMAX = mi;
					interpolate_image_amp_3d(registered);

					//printf("m = %d, mc = %d um = %f ua = %f  MI = %f %f %d\n", m, mc, um, ua, mi, LAMDA, iter);
				}


				for (k = 1; k <= NUMBER_OF_NODES; k++)
				{
					kp = PMAP[k];
					if (kp > 0)
						UO[kp-1] = UT[k-1],VO[kp-1] = VT[k-1],WO[kp-1] = WT[k-1],AO[kp-1] = AT[k-1];
				}

				if (um > SIDE && m > 0)
				{
					printf("registration not converging\n");
					return(1);
				}


				if ((uat <= CLIMIT))
				{
					if (m == 0)
						return(1);

					mc = 0;
					if (ITERM >= 1)
					{
						LAMDA /= 2.0;
						for (k = 1; k <= NN*ND; k++)
							BL[k] /= sqrt(2.0);
						ITERM--;
						fast = 0;
						MUTUAL = 1;
						m++;
					}
				}
				else
				{
					if (ITERM == 0)
						m = MR;
					mc++;
					m++;
				}

	} while (m < MR && ITERM > 0);

	return(0);
}


void get_map_3d(float ***fixed, float ***roi)
{
    int row, col, slice;
    int k;
    int *pp;
	int p;
 	int na, nb, nc, nd, ne, nf, ng, nh;
	float xoff, yoff, zoff;
	int ip, jp, kp;

	extern int NUMBER_OF_NODES;
	extern int NN;
	extern int MSLICES;
 	extern int MROWS;
    extern int MCOLS;
    extern int XP, YP, ZP;
	extern float *XC, *YC, *ZC;
	extern float *MAP;
	extern int *PMAP;
	extern int SIDE;
	extern int REDUCE;

    pp = iarray(NODES_PER_ELEMENT);
	PMAP = ivector(NUMBER_OF_NODES*ND);

	for (k = 1; k <= NUMBER_OF_NODES*ND; k++)
		PMAP[k] = 0;

	for (k = 0; k < NUMBER_OF_NODES; k++)
		MAP[k] = 0.0;

	for (kp = 0; kp < ZP-1; kp++)
		for (jp = 0; jp < YP-1; jp++)
			for (ip = 0; ip < XP-1; ip++)
			{
				na = pnt(ip, jp, kp);
				nb = pnt(ip+1, jp, kp);
				nc = pnt(ip, jp+1, kp);
				nd = pnt(ip+1, jp+1, kp);
				ne = pnt(ip, jp, kp+1);
				nf = pnt(ip+1, jp, kp+1);
				ng = pnt(ip, jp+1, kp+1);
				nh = pnt(ip+1, jp+1, kp+1);
				xoff = XC[na] - 0.5;
				yoff = YC[na] - 0.5;
				zoff = ZC[na] - 0.5;

				for (slice = zoff; slice < zoff+SIDE; slice++)
					for (row = yoff; row < yoff+SIDE; row++)
						for (col = xoff; col < xoff+SIDE; col++)
							if (col >= 0 && col < MCOLS && row >= 0 && row < MROWS && slice >= 0 && slice < MSLICES)
							{
								if (roi[slice][row][col] > 0.0)
								{
									MAP[na]++;
									MAP[nb]++;
									MAP[nc]++;
									MAP[nd]++;
									MAP[ne]++;
									MAP[nf]++;
									MAP[ng]++;
									MAP[nh]++;
								}
							}
			}

	if (REDUCE == 0)
	{
		for (k = 1; k <= NUMBER_OF_NODES*ND; k++)
			PMAP[k] = k;
		NN = NUMBER_OF_NODES;
		return;
	}

	for (k = 1; k <= NUMBER_OF_NODES; k++)
		if (MAP[k-1] > 0)
			PMAP[k] = PMAP[k+NUMBER_OF_NODES] = PMAP[k+NUMBER_OF_NODES*2] = PMAP[k+NUMBER_OF_NODES*3] = 1;

	p = 1;
	for (k = 1; k <= NUMBER_OF_NODES*ND; k++)
		if (PMAP[k] > 0)
		{
			PMAP[k] = p;
			p++;
		}
	NN = p-1;
}

void get_B_and_lamda_3d(float ***fixed, float ***roi)
{
    int row, col, slice;
    int i, j, k;
    float q[32];
    int qj[32];
	int qi[32][32];
	double PTP[32][32];
    float grad[4];
    float tota, totc;
    float *w;
    int *pp;
    int a;
	float totlam, ttot;
	int na, nb, nc, nd, ne, nf, ng, nh;
	float wr, ws, wt;
	float xoff, yoff, zoff;
	int ip, jp, kp;
	float bs, bss;
	double qqj;

	extern double **P;
	extern double *QD;
	extern int **N;
	extern int ND;
	extern int *NP;
	extern int QL;
	extern int QR;
	extern float BI;

	extern int MSLICES;
 	extern int MROWS;
    extern int MCOLS;
    extern int XP, YP, ZP;
    extern int NUMBER_OF_NODES;
	extern int NN;
	extern int NN;
	extern float LAMDA;
	extern int SETLAM;
	extern int REDUCE;
	extern float ALAM;
	extern float K;
	extern float *BL;
	extern int *NO;

	get_map_3d(fixed, roi);
	make_q();

    w = array(NODES_PER_ELEMENT);
    pp = iarray(NODES_PER_ELEMENT);

	for (j = 1; j <= NUMBER_OF_NODES*ND; j++)
		for (i = 1; i <= QL; i++)
			P[j][i] = 0.0;

	for (kp = 0; kp < ZP-1; kp++)
	  for (jp = 0; jp < YP-1; jp++)
		for (ip = 0; ip < XP-1; ip++)
		{
			na = pnt(ip, jp, kp);
			nb = pnt(ip+1, jp, kp);
			nc = pnt(ip, jp+1, kp);
			nd = pnt(ip+1, jp+1, kp);
			ne = pnt(ip, jp, kp+1);
			nf = pnt(ip+1, jp, kp+1);
			ng = pnt(ip, jp+1, kp+1);
			nh = pnt(ip+1, jp+1, kp+1);
			xoff = XC[na] - 0.5;
			yoff = YC[na] - 0.5;
			zoff = ZC[na] - 0.5;

			pp[0] = na;
			pp[1] = nb;
			pp[2] = nc;
			pp[3] = nd;
			pp[4] = ne;
			pp[5] = nf;
			pp[6] = ng;
			pp[7] = nh;

			a = 0;
			for (i = 0; i < 4; i++)
				for (j = 0; j < 8; j++)
					qj[a++] = pp[j] + i*NUMBER_OF_NODES +1;

			for (j = 0; j < 32; j++)
				for (i = 0; i < 32; i++)
					PTP[j][i] = 0.0;

			for (j = 0; j < QR; j++)
				where_in_list(N[qj[j]],qj,qi[j]);

			for (slice = zoff; slice < zoff+SIDE; slice++)
				for (row = yoff; row < yoff+SIDE; row++)
					for (col = xoff; col < xoff+SIDE; col++)
						if (col >= 0 && col < MCOLS && row >= 0 && row < MROWS && slice >= 0 && slice < MSLICES)
						{
							wt=(slice-zoff+0.5)/SIDE;
							ws=(row-yoff+0.5)/SIDE;
							wr=(col-xoff+0.5)/SIDE;
							w[0] = (1-ws)*(1-wr)*(1-wt);
							w[1] = wr*(1-ws)*(1-wt);
							w[2] = ws*(1-wr)*(1-wt);
							w[3] = ws*wr*(1-wt);
							w[4] = (1-ws)*(1-wr)*wt;
							w[5] = wr*(1-ws)*wt;
							w[6] = ws*(1-wr)*wt;
							w[7] = ws*wr*wt;
							if (roi[slice][row][col] >= 1.0)
							{
								grad[0] = (fixed[slice][row][col+1]-fixed[slice][row][col-1])/(float)2.0;
								grad[1] = (fixed[slice][row+1][col]-fixed[slice][row-1][col])/(float)2.0;
								grad[2] = (fixed[slice+1][row][col]-fixed[slice-1][row][col])/(float)2.0;
								grad[3] = -fixed[slice][row][col];

								a = 0;
								for (i = 0; i < 4; i++)
									for (j = 0; j < 8; j++)
										q[a++] = w[j]*grad[i];


								for (j = 0; j < QR; j++)
								{
									qqj = q[j];
									for (i = j; i < QR; i++)
										PTP[j][i] += qqj*q[i];
								}
							}
						}

						for (j = 0; j < QR; j++)
							for (i = j; i < QR; i++)
								PTP[i][j] = PTP[j][i];

						for (j = 0; j < QR; j++)
							for (i = 0; i < QR; i++)
								P[qj[j]][qi[j][i]] += PTP[j][i];
		}
		NN = NUMBER_OF_NODES;

		if (REDUCE == 1)
			pack_matrices();

		tota = 0.0;
		for (j = 0; j < ND-1; j++)
			for (k = 1; k <= NN; k++)
				tota += P[k+j*NN][NO[k+j*NN]];

		totc = 0.0;
		for (k = 1; k <= NN; k++)
			totc += P[k+NN*(ND-1)][NO[k+(ND-1)*NN]];

		bs = tota/(totc*(ND-1));
		bss = sqrt(bs);
		BI = bss;

		for (k = 1; k <= NN*(ND-1); k++)
			for (j = 1; j <= NP[k]; j++)
				if (N[k][j] > NN*(ND-1))
					P[k][j] *= bss;

		for (k = NN*(ND-1) + 1; k <= NN*ND; k++)
			for (j = 1; j <= NP[k]; j++)
				if (N[k][j] <= NN*(ND-1))
					P[k][j] *= bss;
				else
					P[k][j] *= bs;

		ttot = 0.0;
		for (j = 0; j < ND; j++)
			for (k = 1; k <= NN; k++)
				ttot += P[k+j*NN][NO[k+j*NN]];

		totlam = 0.0;
		for (k = 1; k <= NN*ND; k++)
			totlam += QD[k];

		ALAM = ttot/totlam;
		ALAM = sqrt(ALAM);

		if (SETLAM == 1)
			LAMDA = brent((float)1, (float)10, (float)100, (float)0.01);

		SETLAM = 0;

		for (k = 1; k <= NN*(ND-1); k++)
			BL[k] = ALAM*sqrt(LAMDA);

		for (k = NN*(ND-1)+1; k <= NN*ND; k++)
			BL[k] = ALAM*sqrt(K*LAMDA);

}

void crop_roi(float ***roi)
{
	int col, row, slice;

	extern int ND;
	extern int MCOLS, MROWS, MSLICES;

		for (row = 0; row < MROWS; row++)
			for (col = 0; col < MCOLS; col++)
				roi[0][row][col] = roi[MSLICES-1][row][col] = roi[1][row][col] = roi[MSLICES-2][row][col] = 0;

	for (slice = 0; slice < MSLICES; slice++)
		for (col = 0; col < MCOLS; col++)
			roi[slice][0][col] = roi[slice][MROWS-1][col] = roi[slice][1][col] = roi[slice][MROWS-2][col] = 0;

	for (row = 0; row < MROWS; row++)
		for (slice = 0; slice < MSLICES; slice++)
			roi[slice][row][0] = roi[slice][row][MCOLS-1] = roi[slice][row][1] = roi[slice][row][MCOLS-2] = 0;

	for (slice = 0; slice < MSLICES; slice++)
		for (row = 0; row < MROWS; row++)
			for (col = 0; col < MCOLS; col++)
				if (roi[slice][row][col] > 0.5)
					roi[slice][row][col] = 1.0;
				else
					roi[slice][row][col] = 0.0;
}

void get_limits_roi(float ***roi)
{
    int col, row, slice;
    int count;

	extern int CMIN, CMAX, RMIN, RMAX, SMIN, SMAX;
    extern int MSLICES;
	extern int MROWS;
	extern int MCOLS;
	extern int ND;


		SMIN = -1;
		SMAX = -1;
		for (slice = 0; slice < MSLICES; ++slice)
		{
			count = 0;
			for (row = 0; row < MROWS; ++row)
				for (col = 0; col < MCOLS; ++col)
					if (roi[slice][row][col] != 0)
						count++;

			if (count != 0 && SMIN < 0)
				SMIN = slice;
			else if (count == 0 && SMIN >= 0)
			{
				SMAX = slice -1;
				break;
			}
		}
		if (SMAX == -1)
				SMAX = MSLICES-1;


    RMIN = -1;
    RMAX = -1;
    for (row = 0; row < MROWS; ++row)
    {
        count = 0;
		for (slice = 0; slice < MSLICES; ++slice)
			for (col = 0; col < MCOLS; ++col)
				if (roi[slice][row][col] != 0)
					count++;

		if (count != 0 && RMIN < 0)
			RMIN = row;
		else if (count == 0 && RMIN >= 0)
		{
			RMAX = row -1;
			break;
		}
    }
	if (RMAX == -1)
			RMAX = MROWS-1;

    CMIN = -1;
    CMAX = -1;
    for (col = 0; col < MCOLS; ++col)
    {
        count = 0;
		for (row = 0; row < MROWS; ++row)
			for (slice = 0; slice < MSLICES; ++slice)
				if (roi[slice][row][col] != 0)
					count++;

		if (count != 0 && CMIN < 0)
			CMIN = col;
		else if (count == 0 && CMIN >= 0)
		{
			CMAX = col -1;
			break;
		}
    }
	if (CMAX == -1)
			CMAX = MCOLS-1;

}

void init()
{
    extern int MCOLS, MROWS, MSLICES;
	extern int SIDE;
	extern float CLIMIT;
	extern int NUMBER_OF_NODES;
	extern float *UT, *VT, *WT, *AT;
	extern float LAMDA;
	extern int REDUCE, SETLAM;
	extern int PACK;
	extern int FAST;
	extern float ***MOVED;
	extern float B;
	extern int ND;
	extern float K;
	extern double LINTOL;
	extern float *XC;
	extern int FSM;
	extern int S;
	extern float **TDH;
	extern float *BL;
	extern float MIMAX;
	extern int MUTUAL;
	REDUCE = 0;
	ITERM = 0;
	SETLAM = 0;
	MUTUAL = 0;
	LINTOL = 0.01;

	SIDE = 16;     /* -g opt */
	LAMDA = 64.0; /* -L opt */
	PACK = 0;     /* -p opt */
	FSM = 0;      /* -s opt */
	TDH = array2(64,64);
	MR = 100;                                /* -X opt */
	CLIMIT = 0.1;                            /* -t opt */
	K = 1.0;
 	ND = 4;
	S = 1;

}

void pack(float ***pic, float ***ppic)
{
	int col, row, slice;
	int pcols, prows, pslices;
	float mx;

	extern int MCOLS, MROWS, MSLICES;

	pslices = MSLICES/2;
	prows = MROWS/2;
	pcols = MCOLS/2;

	if (pslices ==0)
	{
		pslices = 1;
		mx = 4;
	}
	else
		mx = 8;


	if (pslices > 1)
		MSLICES = pslices*2;

	MROWS = prows*2;
	MCOLS = pcols*2;

	for (slice = 0; slice < pslices; slice++)
		for (row = 0; row < prows; row++)
			for (col = 0; col < pcols; col++)
				ppic[slice][row][col] = 0;

	for (slice = 0; slice < MSLICES; slice++)
		for (row = 0; row < MROWS; row++)
			for (col = 0; col < MCOLS; col++)
				ppic[slice/2][row/2][col/2] += pic[slice][row][col];

	for (slice = 0; slice < pslices; slice++)
		for (row = 0; row < prows; row++)
			for (col = 0; col < pcols; col++)
				ppic[slice][row][col] /= mx;
}


void resample()
{
    float ***registered;
    image *t_image,*interp_image;
	int i,j,k;
	int nslices;
	extern float ***ROI,***FIXED,***MOVED;

	ROI = load_image_alt(&MCOLS, &MROWS, &MSLICES, 3);
    FIXED =  load_image_alt(&MCOLS, &MROWS, &MSLICES, 1);
	MOVED = load_image_alt(&MCOLS, &MROWS, &MSLICES, 2);
    registered = new_image(MSLICES, MROWS, MCOLS);
	interpolate_image_3d(registered);
	t_image = alloc_image(MCOLS, MROWS, MSLICES);
	for (i = 0; i<MSLICES; i++)
			for (j=0; j<MROWS; j++)
					for (k=0; k<MCOLS; k++)
						t_image->data[i*t_image->rows*t_image->cols + j*t_image->cols + k] = registered[i][j][k];

	t_image->sepx = t_image->sepy = t_image->sepz = SLICESEP;

//	if (ORIGSEP != SLICESEP && MSLICES > 1) {
//			interp_image = interpolate_slices(t_image,ORIGSEP);
//			free_image(t_image);
//			t_image = interp_image;
//	}
	nslices = t_image->slices;
	for (i = 0; i<nslices; i++)
			for (j=0; j<MROWS; j++)
					for (k=0; k<MCOLS; k++)
						mask->data[i*t_image->rows*t_image->cols + j*t_image->cols + k] = t_image->data[i*t_image->rows*t_image->cols + j*t_image->cols + k] ;
    //free_image(t_image);
}

void go(int argc, char *argv[])
{
    float ***registered;
	float ***ROI;
	float off;
	float mi;

	extern float ***FIXED;
    extern int MCOLS, MROWS, MSLICES;
	extern int SIDE;
	extern float CLIMIT;
	extern int NUMBER_OF_NODES;
	extern float *UT, *VT, *WT, *AT;
	extern float LAMDA;
	extern int REDUCE, SETLAM;
	extern int PACK;
	extern int FAST;
	extern float ***MOVED;
	extern float B;
	extern int ND;
	extern float K;
	extern double LINTOL;
	extern float *XC;
	extern int FSM;
	extern int S;
	extern float **TDH;
	extern float *BL;
	extern float MIMAX;
	extern int MUTUAL;
	extern int m;
	extern image *mask;
	int k;
	int pcols, prows, pslices;

    FIXED =  load_image_alt(&MCOLS, &MROWS, &MSLICES, 1);
	MOVED = load_image_alt(&MCOLS, &MROWS, &MSLICES, 2);
	ROI = load_image_alt(&MCOLS, &MROWS, &MSLICES, 3);
	printf("loaded\n");
	S = 1;
	if (PACK > 0)
	{
		for (k = 0; k < PACK; k++)
		{
        printf("pack\n");
			pslices = MSLICES/2;
			prows = MROWS/2;
			pcols = MCOLS/2;

			if (pslices ==0)
				pslices = 1;

			ITEMP = new_image(pslices, prows, pcols);
			pack(ROI, ITEMP);
			delete_image(ROI);
			ROI = ITEMP;

			ITEMP = new_image(pslices, prows, pcols);
			pack(FIXED, ITEMP);
			delete_image(FIXED);
			FIXED = ITEMP;

			ITEMP = new_image(pslices, prows, pcols);
			pack(MOVED, ITEMP);
			delete_image(MOVED);
			MOVED = ITEMP;

			MSLICES = pslices;
			MROWS = prows;
			MCOLS = pcols;

			S = S*2;
		}
	}

	crop_roi(ROI);
	printf("crop\n");
	get_limits_roi(ROI);
	printf("lims\n");
	SIDE = (int)(SIDE/S);
	CLIMIT /= S;

	make_mesh_data();
	printf("mesh\n");

	normalise(FIXED, ROI, 1.0);
	printf("norm1\n");
	get_B_and_lamda_3d(FIXED,ROI);
	printf("Blam\n");


	normalise(MOVED, ROI, 1.0);

	registered = new_image(MSLICES, MROWS, MCOLS);
	first_comp();
	printf("comp\n");

	MIMAX = MI(FIXED, ROI, MOVED);
	mi = MIMAX;
	register_image_full_3d(FIXED, ROI, registered);
	//MIMAX = MI(FIXED, ROI, MOVED);
	printf("MI = %f -> %f  in %d iterations\n",mi,MI(FIXED, ROI, registered),m);

	delete_image(FIXED);
	delete_image(ROI);
	delete_image(MOVED);

/* now do final interpolation */

	MOVED = load_image_alt(&MCOLS, &MROWS, &MSLICES, 2);

	off = 0.5;
	for (k = 0; k < NUMBER_OF_NODES; k++)
	{
		AT[k] *= BI;
		XC[k] = (float)(XC[k] - off)*S + off;
		YC[k] = (float)(YC[k] - off)*S + off;
		ZC[k] = (float)(ZC[k] - off)*S + off;
		UT[k] *= S;
		VT[k] *= S;
		WT[k] *= S;
	}

	SIDE = SIDE*S;
}

