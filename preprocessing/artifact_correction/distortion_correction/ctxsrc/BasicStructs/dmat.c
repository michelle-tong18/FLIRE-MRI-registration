///////////////////////////////////////////////////////////////////////
// This source code is Copyright ©2007-2010 CorTechs Labs
// Inc. Permission is granted by CorTechs Labs to copy, distribute and
// modify this source code, provided that it is used for research
// purposes only, and only by not-for-profit organizations, and
// further provided that:
// 
//    * this copyright notice appears in all copies.
//    * you indemnify and hold harmless CorTechs Labs and its
//      successors and associates from any and all liability from any
//      use of the information.
//    * you ensure that these conditions apply to any further
//      redistribution of the source code or derived software
// 
// For commercial entities and for commercial applications, the right
// to copy, distribute or modify this source code is retained by
// CorTechs Labs, and must be obtained through completion of a written
// license agreement with CorTechs Labs Inc.
///////////////////////////////////////////////////////////////////////

#ifdef _WIN32
#define DLLTYPE  __declspec(dllexport)
#endif
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
//#include <mathimf.h>
#include <math.h>
#include "hvecmat.h"
#include "dmat.h"
#include <mkl.h>
//#include "../f2c.h"
#define EPS 1.0e-10

dmat* dMatNew(int rows, int cols, double *data)
{
	dmat *mat = (dmat *) malloc(sizeof(dmat));
	mat->rows = rows;
	mat->cols = cols;
	if (data)
	{
		mat->data = data;
		mat->ballocated = 0;
	}
	else
	{
		mat->data = (double *) calloc(rows*cols, sizeof(double));
		mat->ballocated = 1;
	}
	return mat;
}

void dMatDel(dmat* mat)
{
	if (mat->ballocated == 1)
		free(mat->data);
	mat->data = NULL;
	free(mat);
	
}

dmat* dMatCopy(dmat* mat, int bcopydata)
{
	dmat *matc = dMatNew(mat->rows, mat->cols, NULL);
	if (bcopydata)
	{
		int tsize= mat->rows * mat->cols;
		for (int i=0;i<tsize;i++)
			matc->data[i] = mat->data[i];
	}
	return matc;
}

dmat* dMatEye(int n)
{
	dmat *mat = dMatNew(n, n, NULL);
	for (int i=0;i<n;i++)
		dMatSetVal(mat, i, i, 1.);
	return mat;
}

void dMatPrint(dmat *mat, char *fn)
{
	int i, j;
	if (fn == NULL)
	{
		for (i=0;i<mat->rows;i++)
		{
			for (j=0;j<mat->cols;j++)
				printf("%f ", mat->data[i*mat->cols+j]);
			 printf("\n");
		}	
	}
	else
	{
		FILE *fp = fopen(fn, "w+");
    	if (fp)
    	{
			for (i=0;i<mat->rows;i++)
			{
				for (j=0;j<mat->cols;j++)
					fprintf(fp, "%f ", mat->data[i*mat->cols+j]);
			}
			fclose(fp);
    	}	
	}
}

dmat* dMatfromFile(int rows, int cols, char *fn)
{
	FILE *fp = fopen(fn, "r");
	dmat *m = NULL;
	if (fp)
	{
		m = dMatNew(rows, cols, NULL);
		char tmp[256];
		for (int i=0;i<rows;i++)
		{
			for (int j=0;j<cols;j++)
			{
				fscanf(fp, "%s ", &tmp);
				//printf("\n (%d %d) = %s", i, j, tmp);
				dMatSetVal(m, i, j, atof(tmp));
			}
		}
		fclose(fp);
		
	}
	return m;
}

void dMatSetVal(dmat * mat, int i, int j, double val)
{
  	mat->data[i*mat->cols+j] = val;
}

double dMatGetVal(dmat * mat, int i, int j)
{
	return mat->data[i*mat->cols+j];
}


dmat* dMatGetRow(dmat * mat, int i)
{
	dmat *r = dMatNew(1, mat->cols, NULL);
	for (int j=0;j<mat->rows;j++)
		r->data[j] = dMatGetVal(mat, i, j);
	return r;
}

dmat* dMatGetCol(dmat * mat, int j)
{
	dmat *c = dMatNew(mat->rows, 1,  NULL);
	for (int i=0;i<mat->cols;i++)
		c->data[i] = dMatGetVal(mat, i, j);
	return c;
}


int dMatSameSize(dmat * m1, dmat *m2)
{
	if (m1->rows != m2->rows)
		return 0;
	if (m1->cols != m2->cols)
		return 0;
	return 1;
}

dmat* dMatMultiply(dmat * matL, dmat *matR, double alpha, double beta)
{
	if (matL->cols != matR->rows)
		return NULL;
	dmat *m = dMatNew(matL->rows, matR->cols, NULL);
	char T='N';
	dgemm( &T, &T, &(matL->rows), &(matR->cols), &(matL->cols), &alpha, 
		 matL->data, &(matL->rows), matR->data, &(matR->rows), &beta, m->data, &(m->rows));
	return m;
}


dmat* dMatInv(dmat * mat)
{
	if (mat->rows != mat->cols)
		return NULL;
	dmat *m = dMatCopy(mat, 1);
	int n = mat->cols;
	int info;
	int *ipiv = (int*)calloc(n, sizeof(int));
	double *work = (double*)calloc(2*n, sizeof(double));
	int lwork =2*n;
	dgetrf(&n, &n, m->data, &n, ipiv, &info);
	dgetri(&n, m->data, &n, ipiv, work, &lwork, &info);
	free(work);
	free(ipiv);
	return m;
}

dmat * dMatAdd(dmat *m1, dmat *m2, double a)
{
	//ma =a*m1+m2;
	dmat *ma = NULL;
	if (dMatSameSize(m1,m2))
	{
		ma = dMatCopy(m2, 1);
		int n = m2->rows*m1->cols;
		int inc= 1;
		daxpy(&n, &a, m1->data, &inc, ma->data, &inc);
	}
	return ma;
}

double dMatNorm2(dmat *m, int n)
{
	if (n < 0)
		n = m->rows*m->cols;
	int inc= 1;
	return dnrm2(&n, m->data, &inc);
}

double dMatDot(dmat *m1, dmat *m2, int n)
{
	double rt=0.;
	if (dMatSameSize(m1,m2))
	{
		if (n < 0)
			n = m1->rows*m1->cols;
		int inc=1;
		rt = ddot(&n, m1->data, &inc, m2->data, &inc);
	}
	return rt;
}

dmat * dMatScal(dmat *m,  double a)
{
	dmat *ms = dMatCopy(m, 1);
	int inc = 1;
	int n = m->rows*m->cols;
	dscal(&n, &a, ms->data, &inc);
	return ms;
}

dmat* dMatCross3(dmat *m1, dmat *m2)
{
	dmat *mc=NULL;
	if (dMatSameSize(m1, m2))
	{
		mc = dMatCopy(m1, 0);
		mc->data[0]=m1->data[1]*m2->data[2]-m1->data[2]*m2->data[1];
		mc->data[1]=m1->data[2]*m2->data[0]-m1->data[0]*m2->data[2];
	    mc->data[2]=m1->data[0]*m2->data[1]-m1->data[1]*m2->data[0];
	    mc->data[3] =1.;
	}
	
	return mc;
}

double dMatComputeAngle(dmat *m1, dmat *m2, int n)
{
	
	if (dMatSameSize(m1, m2))
	{
		if (n<0)
			n = m1->rows*m1->cols;
		double n1 = dMatNorm2(m1, n);
		double n2 = dMatNorm2(m1, n);
		double dot12 = dMatDot(m1, m2, n);
		double stuff = dot12/(n1*n2);
		/* "stuff" is the argument to acos: make sure it is in [-1,1] */
	  
	  	if	( stuff < -1 )
	    {
	      if( stuff > -(1+EPS) )
		  	stuff = -1;
	    }
	
	  	if( stuff > +1 )
	    {
	      if( stuff < +(1+EPS) )
		  	stuff = +1;
	    }
	  
	  	/* Calc theta */ 
	  	return acos(stuff);
	}
	else
		return 0.;
}

dmat* dMatRotI(double radian)
{
	dmat *mat = dMatEye(4);
	double sa = sin(radian);
	double ca = cos(radian);
	dMatSetVal(mat,1,1,ca);
	dMatSetVal(mat,1,2,sa);
	dMatSetVal(mat,2,1,-sa);
	dMatSetVal(mat,2,2,ca);
	return mat;
}

dmat* dMatRotJ(double radian)
{
	dmat *mat = dMatEye(4);
	double sa = sin(radian);
	double ca = cos(radian);
	dMatSetVal(mat,0,0,ca);
	dMatSetVal(mat,0,2,sa);
	dMatSetVal(mat,2,0,-sa);
	dMatSetVal(mat,2,2,ca);
	return mat;
}

dmat* dMatRotK(double radian)
{
	dmat *mat = dMatEye(4);
	double sa = sin(radian);
	double ca = cos(radian);
	dMatSetVal(mat,0,0,ca);
	dMatSetVal(mat,0,1,sa);
	dMatSetVal(mat,1,0,-sa);
	dMatSetVal(mat,1,1,ca);
	return mat;	
}

dmat* dMatTranslate(double ti, double tj, double tk)
{
	dmat *mat = dMatEye(4);
	dMatSetVal(mat,0,3, ti);
	dMatSetVal(mat,1,3, tj);
	dMatSetVal(mat,2,3, tk);
	return mat;		
}


dmat* dMatInvhMat(hmat *m, hmat *mi)
{
	dmat *dm = dMatNew(4, 4, NULL);
	int i,j;
	for (i=0;i<4;i++)
	for (j=0;j<4;j++)
		dMatSetVal(dm, i, j, hMatGetVal(m, i, j));
	
	dmat *dmi = dMatInv(dm);
	
	for (i=0;i<4;i++)
	for (j=0;j<4;j++)
		hMatSetVal(mi, i, j, dMatGetVal(dmi, i, j));
		
	dMatDel(dmi);
	dMatDel(dm);
}
