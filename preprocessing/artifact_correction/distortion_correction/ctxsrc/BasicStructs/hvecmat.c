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

/****************************************************
 *  hvecmat.c                                         
 *  Created on: 13-Dec-2006 08:32:09                      
 *  Implementation of the Class hvecmat       
 *  Original author: Gennan Chen                     
 ****************************************************/

#ifdef _WIN32
#define DLLTYPE  __declspec(dllexport)
#endif
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
//#include <mathimf.h>
#include <math.h>
#include "hvecmat.h"
#define EPS 1.0e-10


void hVecInit(hvec *v1,float val1, float val2, float val3, float val4)
{
	float *fp = (float*)&(v1->fVect4);
	*(fp) = val1;
	*(fp+1) = val2;
	*(fp+2) = val3;
	*(fp+3) = val4;
} 

void hVec2fPtr(hvec *v1,float * fptr)
{
	float *fp = (float*)&(v1->fVect4);
	int i;
	for (i=0;i<4;i++)
		fptr[i]=*(fp+i);
}


float hVecDot(hvec *v1, hvec *v2)
{
	__m128 fvm=_mm_mul_ps(v1->fVect4, v2->fVect4);
	__m128 fa1 = _mm_add_ps(fvm, _mm_movehl_ps(fvm, fvm));
	__m128 ftemp = _mm_add_ss(fa1, _mm_shuffle_ps(fa1, fa1, 1));
	float *fp = (float*)&ftemp;
	return 	*fp;

}

float hVecDot3(hvec *v1, hvec *v2)
{
	float *f1 = (float*)&(v1->fVect4);
	float *f2 = (float*)&(v2->fVect4);
	return 	(hVecDot(v1,v2) -(f1[3]*f2[3]));

}

void hVecCrossProduct3(hvec *v1, hvec *v2, hvec *vout)
{
	float *pu =  (float*)&(v1->fVect4);
	float *pv =  (float*)&(v2->fVect4);
	float *pc =  (float*)&(vout->fVect4);
	pc[0]=pu[1]*pv[2]-pu[2]*pv[1];
	pc[1]=pu[2]*pv[0]-pu[0]*pv[2];
	pc[2]=pu[0]*pv[1]-pu[1]*pv[0];
	pc[3] =1.;
}

double hVecComputeAngle(hvec *v1, hvec *v2)
{
	double m1 = sqrt(hVecDot3(v1,v1));
	double m2 = sqrt(hVecDot3(v2,v2));
	double dot12 = hVecDot3(v1, v2);
	double stuff = dot12/(m1*m2);
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

void hVecFromfPtr(float* fptr, hvec* v1)
{
	float *fp = (float*)&(v1->fVect4);
	for (int i=0;i<4;i++)
		*(fp+i)=fptr[i];
	
}


void hVecPrint(hvec *v)
{
	float *fp = (float*)&(v->fVect4);
	printf("\n %f %f %f %f \n", fp[0], fp[1], fp[2], fp[3]);
}

void hVecSetVal(hvec *v1, int i, float val)
{
	float *fp = (float*)&(v1->fVect4);
	*(fp+i) = val;
}

float hVecGetVal(hvec *v1, int i)
{
	  float *fp = (float*)&(v1->fVect4);
	  return *(fp+i);
}

float hVecNorm(hvec *v)
{
	return sqrtf(hVecDot(v,v));
}


float* hVecGetfPtr(hvec *v)
{
	return  (float*)&(v->fVect4);
}

 
void hVecMultiplyfVal(hvec *v, float val, hvec* vo)
{
	__m128 fv = _mm_set_ps(val, val, val, val);
	vo->fVect4 = _mm_mul_ps(v->fVect4, fv);
}

void hVecMultiplyhVec(hvec *v1, hvec *v2, hvec* vo)
{
	vo->fVect4 = _mm_mul_ps(v1->fVect4, v2->fVect4);
}

void hVecAdd(hvec *v1, hvec *v2, hvec* vo)
{
	vo->fVect4 = _mm_add_ps(v1->fVect4, v2->fVect4);
}

void hVecSubtracts(hvec *v1, hvec *v2, hvec* vo)
{
    vo->fVect4 = _mm_sub_ps(v1->fVect4, v2->fVect4);
}

float hVecSum(hvec *v)
{
	__m128 fa1 = _mm_add_ps(v->fVect4, _mm_movehl_ps(v->fVect4, v->fVect4));
	__m128 ftemp = _mm_add_ss(fa1, _mm_shuffle_ps(fa1, fa1, 1));

	float *fp = (float*)&ftemp;
	return 	*fp;
}

void hMat2fPtr(hmat *mat,float * fptr)
{
	int indx=0,i,j;
	for (i=0;i<4;i++)
	{
		float *fp = (float*)&(mat->rows[i]);
		for (j=0;j<4;j++)
		{
			 *(fptr+indx) = *(fp+j);
			indx++;
		}
		
	}
}

void hMatFromfPtr(float * fptr, hmat *mat)
{
	int indx=0,i,j;
	for (i=0;i<4;i++)
	{
		float *fp = (float*)&(mat->rows[i]);
		for (j=0;j<4;j++)
		{
			*(fp+j) = *(fptr+indx);
			indx++;
		}
		
	}
}

void hMatEye(hmat *mat)
{
	for (int i=0;i<4;i++)
	{
		float *fp = (float*)&(mat->rows[i]);
		for (int j=0;j<4;j++)
		{
			if (i==j)
			 	*(fp+j) = 1.;
			else
				*(fp+j) = 0.; 
		}
		
	}
}
void hMatSetVal(hmat *mat, int i, int j, float val)
{
	float *fp = (float*)&(mat->rows[i]);
	*(fp+j) = val;
}

float hMatGetVal(hmat *mat, int i, int j)
{
	float *fp = (float*)&(mat->rows[i]);
	return *(fp+j);
}


void hMatMultiplyhVec(hmat *mat, hvec *vin, hvec *vout)
{
	float *fpout = (float*)&(vout->fVect4);
    for (int i=0;i<4;i++)
    {
        hvec htemp;
        htemp.fVect4 = mat->rows[i];
        fpout[i] = hVecDot(&htemp, vin);
    }
	
}

void hMatMultiplyhMat(hmat *matL, hmat *matR, hmat *matout)
{
	hmat matRT;
	hMatTranspose(matR, &matRT);

	for (int i=0;i<4;i++)
	{
		float *fp = (float*)&(matout->rows[i]);
		for (int j=0;j<4;j++)
        {
            hvec htemp;
            htemp.fVect4 = matL->rows[i];
            hvec htemp1;
            htemp1.fVect4 = matRT.rows[j];
            fp[j] = hVecDot(&htemp, &htemp1);
        }
	}
	
}

void hMatInverse(hmat *matin, hmat *matout)
{
	__m128 tmp,tmp1;
	float r0[4];
	int i;
	float *fp = (float*)&(matin->rows[0]);
	for (i=0;i<4;i++)
	    r0[i]=fp[i];
	    
	float r1[4];
	fp = (float*)&(matin->rows[1]);
	for (i=0;i<4;i++)
	    r1[i]=fp[i];
	        
	float r2[4];
	fp = (float*)&(matin->rows[2]);
	for (i=0;i<4;i++)
	    r2[i]=fp[i];
	    
	float r3[4];
	fp = (float*)&(matin->rows[3]);
	for (i=0;i<4;i++)
	    r3[i]=fp[i];    
	
	__m128 tr0,tr1,tr2,tr3;

	//Not the typical Transfpose
	tmp=_mm_loadl_pi(tmp,(__m64*)r0);
	tmp=_mm_loadh_pi(tmp,(__m64*)r1);
	tmp1=_mm_loadl_pi(tmp1,(__m64*)r2);
	tmp1=_mm_loadh_pi(tmp1,(__m64*)r3);
	
	tr0=_mm_shuffle_ps(tmp, tmp1, 0x88);
	tr1=_mm_shuffle_ps(tmp1, tmp, 0xDD);

	tmp=_mm_loadl_pi(tmp,(__m64*)(r0+2));
	tmp=_mm_loadh_pi(tmp,(__m64*)(r1+2));
	tmp1=_mm_loadl_pi(tmp1,(__m64*)(r2+2));
	tmp1=_mm_loadh_pi(tmp1,(__m64*)(r3+2));

	tr2=_mm_shuffle_ps(tmp, tmp1, 0x88);
	tr3=_mm_shuffle_ps(tmp1, tmp, 0xDD);

	__m128 t1,t2;
	__m128 m0,m1,m2,m3,det;

	//
	t1=_mm_mul_ps(tr2,tr3);
	t1=_mm_shuffle_ps(t1,t1,0xB1);

	m0=_mm_mul_ps(tr1,t1);
	m1=_mm_mul_ps(tr0,t1);

	t1=_mm_shuffle_ps(t1,t1,0x4E);
	t2=_mm_mul_ps(tr1,t1);
	m0=_mm_sub_ps(t2,m0);
	t2=_mm_mul_ps(tr0,t1);
	m1=_mm_sub_ps(t2,m1);
	m1=_mm_shuffle_ps(m1,m1,0x4E);

	//
	t1=_mm_mul_ps(tr1,tr2);
	t1=_mm_shuffle_ps(t1,t1,0xB1);

	m0=_mm_add_ps(_mm_mul_ps(tr3,t1),m0);
	m3=_mm_mul_ps(tr0,t1);

	t1=_mm_shuffle_ps(t1,t1,0x4E);
	t2=_mm_mul_ps(tr3,t1);
	m0=_mm_sub_ps(m0,t2);
	t2=_mm_mul_ps(tr0,t1);
	m3=_mm_sub_ps(t2,m3);
	m3=_mm_shuffle_ps(m3,m3,0x4E);

	//
	t1=_mm_mul_ps(_mm_shuffle_ps(tr1,tr1,0x4E),tr3);
	t1=_mm_shuffle_ps(t1,t1,0xB1);
	tr2=_mm_shuffle_ps(tr2,tr2,0x4E);

	m0=_mm_add_ps(_mm_mul_ps(tr2,t1),m0);
	m2=_mm_mul_ps(tr0,t1);

	t1=_mm_shuffle_ps(t1,t1,0x4E);
	t2=_mm_mul_ps(tr2,t1);
	m0=_mm_sub_ps(m0,t2);
	t2=_mm_mul_ps(tr0,t1);
	m2=_mm_sub_ps(t2,m2);
	m2=_mm_shuffle_ps(m2,m2,0x4E);

	//
	t1=_mm_mul_ps(tr0,tr1);
	t1=_mm_shuffle_ps(t1,t1,0xB1);

	m2=_mm_add_ps(_mm_mul_ps(tr3,t1),m2);
	m3=_mm_sub_ps(_mm_mul_ps(tr2,t1),m3);
	t1=_mm_shuffle_ps(t1,t1,0x4E);
	m2=_mm_sub_ps(_mm_mul_ps(tr3,t1),m2);
	m3=_mm_sub_ps(m3,_mm_mul_ps(tr2,t1));

	//
	t1=_mm_mul_ps(tr0,tr3);
	t1=_mm_shuffle_ps(t1,t1,0xB1);

	m1=_mm_sub_ps(m1,_mm_mul_ps(tr2,t1));
	m2=_mm_add_ps(_mm_mul_ps(tr1,t1),m2);
	t1=_mm_shuffle_ps(t1,t1,0x4E);
	m1=_mm_add_ps(_mm_mul_ps(tr2,t1),m1);
	m2=_mm_sub_ps(m2,_mm_mul_ps(tr1,t1));

	//
	t1=_mm_mul_ps(tr0,tr2);
	t1=_mm_shuffle_ps(t1,t1,0xB1);

	m1=_mm_add_ps(m1,_mm_mul_ps(tr3,t1));
	m3=_mm_sub_ps(m3,_mm_mul_ps(tr1,t1));
	t1=_mm_shuffle_ps(t1,t1,0x4E);
	m1=_mm_sub_ps(m1,_mm_mul_ps(tr3,t1));
	m3=_mm_add_ps(_mm_mul_ps(tr1,t1),m3);

	//
	det=_mm_mul_ps(tr0,m0);
	det=_mm_add_ps(_mm_shuffle_ps(det,det,0x4E),det);
	det=_mm_add_ps(_mm_shuffle_ps(det,det,0xB1),det);
	t1=_mm_rcp_ss(det);
	det= _mm_sub_ss(_mm_add_ss(t1,t1),_mm_mul_ss(det,_mm_mul_ss(t1,t1)));
	det=_mm_shuffle_ps(det,det,0x00);

	//
	matout->rows[0]=_mm_mul_ps(det,m0);
	matout->rows[1]=_mm_mul_ps(det,m1);
	matout->rows[2]=_mm_mul_ps(det,m2);
	matout->rows[3]=_mm_mul_ps(det,m3);

	
}

void hMatTranspose(hmat *matin, hmat *matout)
{
	__m128 tmp,tmp1;
	float r0[4];
	int i;
	float *fp = (float*)&(matin->rows[0]);
	for (i=0;i<4;i++)
	    r0[i]=fp[i];
	    
	float r1[4];
	fp = (float*)&(matin->rows[1]);
	for (i=0;i<4;i++)
	    r1[i]=fp[i];
	        
	float r2[4];
	fp = (float*)&(matin->rows[2]);
	for (i=0;i<4;i++)
	    r2[i]=fp[i];
	    
	float r3[4];
	fp = (float*)&(matin->rows[3]);
	for (i=0;i<4;i++)
	    r3[i]=fp[i];    
	    

	tmp=_mm_loadl_pi(tmp,(__m64*)r0);
	tmp=_mm_loadh_pi(tmp,(__m64*)r1);
	tmp1=_mm_loadl_pi(tmp1,(__m64*)r2);
	tmp1=_mm_loadh_pi(tmp1,(__m64*)r3);
	
	matout->rows[0]= _mm_shuffle_ps(tmp, tmp1, _MM_SHUFFLE(2,0,2,0));
	matout->rows[1]= _mm_shuffle_ps(tmp, tmp1, _MM_SHUFFLE(3,1,3,1));
	

	tmp=_mm_loadl_pi(tmp,(__m64*)(r0+2));
	tmp=_mm_loadh_pi(tmp,(__m64*)(r1+2));
	tmp1=_mm_loadl_pi(tmp1,(__m64*)(r2+2));
	tmp1=_mm_loadh_pi(tmp1,(__m64*)(r3+2));
	matout->rows[2]= _mm_shuffle_ps(tmp, tmp1, _MM_SHUFFLE(2,0,2,0));
	matout->rows[3]= _mm_shuffle_ps(tmp, tmp1, _MM_SHUFFLE(3,1,3,1));

}

void hMatRotI(hmat *mat, double radian)
{
	double sa = sin(radian);
	double ca = cos(radian);
	hMatEye(mat);
	hMatSetVal(mat,1,1,ca);
	hMatSetVal(mat,1,2,sa);
	hMatSetVal(mat,2,1,-sa);
	hMatSetVal(mat,2,2,ca);
}

void hMatRotJ(hmat *mat, double radian)
{
	double sa = sin(radian);
	double ca = cos(radian);
	hMatEye(mat);
	hMatSetVal(mat,0,0,ca);
	hMatSetVal(mat,0,2,sa);
	hMatSetVal(mat,2,0,-sa);
	hMatSetVal(mat,2,2,ca);
}

void hMatRotK(hmat *mat, double radian)
{
	double sa = sin(radian);
	double ca = cos(radian);
	hMatEye(mat);
	hMatSetVal(mat,0,0,ca);
	hMatSetVal(mat,0,1,sa);
	hMatSetVal(mat,1,0,-sa);
	hMatSetVal(mat,1,1,ca);
}

void hMatTranslate(hmat *mat, double ti, double tj, double tk)
{
	hMatEye(mat);
	hMatSetVal(mat,0,3, ti);
	hMatSetVal(mat,1,3, tj);
	hMatSetVal(mat,2,3, tk);
}


void hMatPrint(hmat *mat)
{
	int i,j;
	printf("\n");
	for (i=0;i<4;i++)
	{
		float *fp = (float*)&(mat->rows[i]);
		for (j=0;j<4;j++)
			printf(" % 3.6f", *(fp+j));
		 printf("\n");
	}	
}

void hMatSetRow(hmat *mat, hvec *v, int i)
{
	float *fp = (float*)&(mat->rows[i]);
	float *fi = (float*)&(v->fVect4);
	for (int j=0;j<4;j++)
		*(fp+j) = *(fi+j);
}

void hMatSetCol(hmat *mat, hvec *v, int j)
{
	float *fi = (float*)&(v->fVect4);
	for (int i=0;i<4;i++)
	{
		float *fp = (float*)&(mat->rows[i]);
		*(fp+j) = *(fi+i);
	}
}

void hMatGetRow(hmat *mat, hvec *v, int i)
{
	float *fp = (float*)&(mat->rows[i]);
	float *fi = (float*)&(v->fVect4);
	for (int j=0;j<4;j++)
		*(fi+j) = *(fp+j);
}

void hMatGetCol(hmat *mat, hvec *v, int j)
{
	float *fi = (float*)&(v->fVect4);
	for (int i=0;i<4;i++)
	{
		float *fp = (float*)&(mat->rows[i]);
		*(fi+i) = *(fp+j);
	}
}

float * hMatGetRowfPtr(hmat *mat, int i)
{
	return (float*)&(mat->rows[i]);
}

int hMat2File(hmat *mat, char * fn)
{
    FILE *fp = fopen(fn, "w+");
    if (fp)
    {
        for (int i=0;i<4;i++)
        {
            float *prow = hMatGetRowfPtr(mat, i);
            fprintf(fp, "\n %f %f %f %f", prow[0], prow[1],prow[2], prow[3]);
            prow = NULL;
        }
        fclose(fp);
        return 1;
    }
    else
        return 0;
}

int hMatFromFile(char *fn, hmat *mat)
{
    FILE *fp = fopen(fn, "r");
    if (fp)
    {
        for (int i=0;i<4;i++)
        {
            float *prow = hMatGetRowfPtr(mat, i);
            fscanf(fp, "\n %f %f %f %f", &prow[0], &prow[1], &prow[2], &prow[3]);
        }
        fclose(fp);
        return 1;
    }
    else
        return 0;
}

void hMatCopy(hmat *mi, hmat *mc)
{
    int i,j;
    for (i=0;i<4;i++)
    {
        for (j=0;j<4;j++)
            hMatSetVal(mc, i,  j, hMatGetVal(mi, i, j));
    }
}

void hMatGetNormCol(hmat *mat, hvec *v, int j)
{
    hvec tmp;
    float *fi = (float*)&(tmp.fVect4);
    for (int i=0;i<4;i++)
    {
        float *fp = (float*)&(mat->rows[i]);
        *(fi+i) = *(fp+j);
    }
    double normi = 1./hVecNorm(&tmp);
    hVecMultiplyfVal(&tmp, normi, v);
}
    
