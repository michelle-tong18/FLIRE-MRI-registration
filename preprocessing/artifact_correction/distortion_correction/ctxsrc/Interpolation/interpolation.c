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
 *  interpolation.c                                         
 *  Created on: 14-Dec-2006 14:50:44                      
 *  Implementation of the Class interpolation       
 ****************************************************/
#ifdef _WIN32
#define DLLTYPE  __declspec(dllexport)
#endif

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
//#include <mathimf.h>
#include <math.h>

#include "interpolation.h"
#define MINJD 0
#define MAXJD 3
#define EPS 2e-16

int round4(double val)
{
    double val4 = val*10000;
    return  round(0.0001*round(val4));
}
void proj(UNWARPTYPE uwtype, hvec *pvec, hvec * prow, hvec * pcol, hvec *pdep, hvec *ppjvec)
{
	
    int i;
    if (uwtype==THROUGHPLANE) // through plane unwarpping
    {
       
       double pjval= hVecDot(pvec,pdep);
       hVecMultiplyfVal(pdep, pjval, ppjvec);
       //hVecSetVal(ppjvec, 3, 0.);
       /*for (i=0;i<3;i++)
           pjval=pjval+pvec[i]*pdep[i];
       for (i=0;i<3;i++)
           ppjvec[i]=pjval*pdep[i];       
       ppjvec[3]=0.;*/
    } 
    else if (uwtype==INPLANE) // inplane unwarping
    {
        //printf("\n inplan unwarping");
        double pjvalr=hVecDot(pvec,prow);
        double pjvalc=hVecDot(pvec,pcol);
        hvec t1,t2;
        hVecMultiplyfVal(prow, pjvalr, &t1);
        hVecMultiplyfVal(pcol, pjvalc, &t2);
        hVecAdd(&t1, &t2, ppjvec);
        //hVecSetVal(ppjvec, 3, 0.);
        /*for (i=0;i<3;i++)
        {
           pjvalr=pjvalr+pvec[i]*prow[i];
           pjvalc=pjvalc+pvec[i]*pcol[i];
         }
        for (i=0;i<3;i++)
            ppjvec[i]=pjvalr*prow[i]+pjvalc*pcol[i];     
        ppjvec[3]=0.;*/
    }
    else
    {
        hVecFromfPtr(hVecGetfPtr(pvec), ppjvec);
        /*for (i=0;i<4;i++)
            ppjvec[i]=pvec[i];*/
    }
	
}

void getuwlph(hvec *plph, fvols *dfield,  hmat* Mlph2vxl_g, hvec *puwlph)
{
    
    hvec dv, vxl;
    hVecInit(&dv, 0, 0, 0, 0);
    hMatMultiplyhVec(Mlph2vxl_g, plph, &vxl);
    GetVxlVals(&vxl, hVecGetfPtr(&dv), dfield, LINEAR, THROWAWAY, 0); 
    hVecAdd(plph, &dv, puwlph);
	/*int i;
    double dx=0.;
    double dv[3], vxl[4];
    mapvxl(Mlph2vxl_g, plph, vxl);
    getvxlval(vxl, dL, 1, 1, &dv[0]);
    getvxlval(vxl, dP, 1, 1, &dv[1]);
    getvxlval(vxl, dH, 1, 1, &dv[2]);
    for (i=0;i<3;i++)
        puwlph[i]=plph[i]+dv[i];
    puwlph[3]=1.;*/
}

double getJac(hvec *plph, UNWARPTYPE uwtype, fvols* dfield,  hmat* Mlph2vxl_g, double d,
			  hvec *prow, hvec *pcol , hvec *pdep)
{
    double jac =1.;
    hvec puw0, puw1, puw2, puw3;
    hvec ptemp, ptemp0;
    hvec pdiff, pdiff1, pdiff2, pdiff3;
    switch (uwtype)
    {
        case THROUGHPLANE:
            getuwlph(plph, dfield, Mlph2vxl_g, &puw0);
            hVecMultiplyfVal(pdep, d, &ptemp0);
            hVecAdd(plph, &ptemp0, &ptemp);
            hVecSetVal(&ptemp, 3, 1);
            getuwlph(&ptemp, dfield, Mlph2vxl_g, &puw1);
            hVecSubtracts(&puw1, &puw0, &pdiff);
            jac = fabs(hVecDot3(&pdiff, pdep)/d);
            break;
        case INPLANE:
          
            getuwlph(plph, dfield, Mlph2vxl_g, &puw0);
            
            hVecMultiplyfVal(prow, d, &ptemp0);
            hVecAdd(plph, &ptemp0, &ptemp);
            hVecSetVal(&ptemp, 3, 1);
            getuwlph(&ptemp, dfield, Mlph2vxl_g, &puw1);
            hVecSubtracts(&puw1, &puw0, &pdiff1);
            
            hVecMultiplyfVal(pcol, d, &ptemp0);
            hVecAdd(plph, &ptemp0, &ptemp);
            hVecSetVal(&ptemp, 3, 1);
            getuwlph(&ptemp, dfield, Mlph2vxl_g, &puw2);
            hVecSubtracts(&puw2, &puw0, &pdiff2);
            
            hVecCrossProduct3(&pdiff1, &pdiff2, &ptemp);
            jac = fabs(hVecDot3(&ptemp, pdep)/(d*d));
            break;
            
        default:
           
            
            //Obatin the morph point
            getuwlph(plph, dfield, Mlph2vxl_g, &puw0);
            
            // add a d times in row direction and map to morph space;
            hVecMultiplyfVal(prow, d, &ptemp0);
            hVecAdd(plph, &ptemp0, &ptemp);
            hVecSetVal(&ptemp, 3, 1);
            getuwlph(&ptemp, dfield, Mlph2vxl_g, &puw1);
            
            // Calculate the difference between them
            hVecSubtracts(&puw1, &puw0, &pdiff1);
            
            // add a d times in col direction and map to morph space;
            hVecMultiplyfVal(pcol, d, &ptemp0);
            hVecAdd(plph, &ptemp0, &ptemp);
            hVecSetVal(&ptemp, 3, 1);
            getuwlph(&ptemp, dfield, Mlph2vxl_g, &puw2);
            
            // Calculate the difference between them 
            hVecSubtracts(&puw2, &puw0, &pdiff2);
            
            // add a d times in dep direction and map to morph space;
            hVecMultiplyfVal(pdep, d, &ptemp0);
            hVecAdd(plph, &ptemp0, &ptemp);
            hVecSetVal(&ptemp, 3, 1);
            getuwlph(&ptemp, dfield, Mlph2vxl_g, &puw3);
             
            // Calculate the difference between them
            hVecSubtracts(&puw3, &puw0, &pdiff3);
            
            // Calcuate the volume ratio.
            hVecCrossProduct3(&pdiff1, &pdiff2, &ptemp);
            jac = fabs(hVecDot3(&ptemp, &pdiff3)/(d*d*d));
    }
   
    //clamping
    if (jac < MINJD)
        jac=MINJD;
    else if (jac > MAXJD)
        jac=MAXJD;

    return jac;
  
   
}

float fDataGetValPAD(volinfo *info, float* data, int i, int j, int k, INTERPMETHOD method, INTERPEDGE edge)
{
	int ind;
	float val = 0;
	switch (method)
	{
		case NEAREST:
			ind = sub2ind3(info, i,j,k);
			val = data[ind];
			break;
		case LINEAR:
		case CUBIC:
			switch (edge)
			{
                case THROWAWAY:
                    val = data[sub2ind3(info, i,j,k)];
                    break;
				case ZEROPAD:
					if (i > (info->dim[0]-1))
						{val=0;break;}
					if (j > (info->dim[1]-1))
						{val=0;break;}
					if (k > (info->dim[2]-1))
						{val=0;break;}
					if (i  < 0)
						{val=0;break;}
					if (j < 0)
						{val=0;break;}
					if (k < 0)
						{val=0;break;}
					val = data[sub2ind3(info, i,j,k)];
					break;
				default: //MIRROR
					if (i > (info->dim[0]-1))
						i = info->dim[0]-1-(i-info->dim[0]);
					if (j > (info->dim[1]-1))
						j = info->dim[1]-1-(j-info->dim[1]);
					if (k > (info->dim[2])-1)
						k = info->dim[2]-1-(k-info->dim[2]);
					if (i < 0)
						i = 0-i;
					if (j < 0)
						j = 0-j;
					if (k < 0)
						k = 0-k;
					val = data[sub2ind3(info, i,j,k)];
			}
			break;
	}
	return val;
	
}
//Nearest Neighbor
int getVxlVal_NN(float *pvxl, volinfo *info, float* data,  float* val)
{
	int i=round4(pvxl[0]);
    int j=round4(pvxl[1]);
    int k=round4(pvxl[2]);
	*val = 0;
   	if (i>(info->dim[0]-1))
		return 0;
	if (j>(info->dim[1]-1))
		return 0;
	if (k>(info->dim[2]-1))
		return 0;
	if (i<0)
		return 0;
	if (j<0)
		return 0;
	if (k<0)
		return 0;
    *val = data[sub2ind3(info, i,j,k)];
	return 1;
}

int getVxlVals_NN(float *pvxl, volinfo *info, float** data,  float* val)
{
    int i=round4(pvxl[0]);
    int j=round4(pvxl[1]);
    int k=round4(pvxl[2]);
    for (int l=0;l<info->dim[3];l++)
        val[l]=0;
    if (i>(info->dim[0]-1))
        return 0;
    if (j>(info->dim[1]-1))
        return 0;
    if (k>(info->dim[2]-1))
        return 0;
    if (i<0)
        return 0;
    if (j<0)
        return 0;
    if (k<0)
        return 0;
    int ind = sub2ind3(info, i,j,k);
    for (int l=0;l<info->dim[3];l++)
    {
        float *fp = data[l];
        val[l] = fp[ind];
    }
    return 1;
}

//Linear Interpolation
int getVxlVals_Linear(float *pvxl, volinfo *info, float** data, INTERPEDGE edge, float* val)
{
 	int l, im=0, ip=0, jm=0, jp=0, km=0, kp=0;
	im=(int)(pvxl[0]);
    ip=im+1;
    jm=(int)(pvxl[1]);
    jp=jm+1;
    km=(int)(pvxl[2]);
    kp=km+1;
	// check if out of bound
    if (edge == THROWAWAY)
    {
        if (ip>(info->dim[0]-1))
            return 0;
        if (jp>(info->dim[1]-1))
            return 0;
        if (kp>(info->dim[2]-1))
            return 0;
        if (im<0)
            return 0;
        if (jm<0)
            return 0;
        if (km<0)
            return 0;
    }
    else
    {
    	if (ip>(info->dim[0]))
    		return 0;
    	if (jp>(info->dim[1]))
    		return 0;
    	if (kp>(info->dim[2]))
    		return 0;
    	if (im<-1)
    		return 0;
    	if (jm<-1)
    		return 0;
    	if (km<-1)
    		return 0;
    }	
	float tim=pvxl[0]-im;
    float tip=1-tim;
    
    float tjm=pvxl[1]-jm;
    float tjp=1-tjm;
    
    float tkm=pvxl[2]-km;
    float tkp=1-tkm;
  
    double t3[4], t4[4];
    t3[0]=tip;t3[1]=tim;t3[2]=tim;t3[3]=tip;
    t4[0]=tjp;t4[1]=tjp;t4[2]=tjm;t4[3]=tjm;
    for (int l=0;l<info->dim[3];l++)
    {
        float *fp = data[l];   
    	double tf1[4],tf2[4];
    	tf1[0]=fDataGetValPAD(info, fp, im, jm, km, LINEAR, edge);
    	tf1[1]=fDataGetValPAD(info, fp, ip, jm, km, LINEAR, edge);
    	tf1[2]=fDataGetValPAD(info, fp, ip, jp, km, LINEAR, edge);
    	tf1[3]=fDataGetValPAD(info, fp, im, jp, km, LINEAR, edge);
    	tf2[0]=fDataGetValPAD(info, fp, im, jm, kp, LINEAR, edge);
    	tf2[1]=fDataGetValPAD(info, fp, ip, jm, kp, LINEAR, edge);
    	tf2[2]=fDataGetValPAD(info, fp, ip, jp, kp, LINEAR, edge);
    	tf2[3]=fDataGetValPAD(info, fp, im, jp, kp, LINEAR, edge);
    
        double tf3[4];
        double tval=0;
    	int i;
        for (i=0;i<4;i++)
            tf3[i]=tkp*tf1[i]+tkm*tf2[i];
    
        for (i=0;i<4;i++)
            tval+=tf3[i]*t3[i]*t4[i];
        val[l] = tval;
    }
	return 1;
}

int getVxlVal_Linear(float *pvxl, volinfo *info, float* data, INTERPEDGE edge, float* val)
{
    int im=0, ip=0, jm=0, jp=0, km=0, kp=0;
    *val =0.;
    im=(int)(pvxl[0]);
    ip=im+1;
    jm=(int)(pvxl[1]);
    jp=jm+1;
    km=(int)(pvxl[2]);
    kp=km+1;
    // check if out of bound
    if (edge == THROWAWAY)
    {
        if (ip>(info->dim[0]-1))
            return 0;
        if (jp>(info->dim[1]-1))
            return 0;
        if (kp>(info->dim[2]-1))
            return 0;
        if (im<0)
            return 0;
        if (jm<0)
            return 0;
        if (km<0)
            return 0;
    }
    else
    {
        if (ip>(info->dim[0]))
            return 0;
        if (jp>(info->dim[1]))
            return 0;
        if (kp>(info->dim[2]))
            return 0;
        if (im<-1)
            return 0;
        if (jm<-1)
            return 0;
        if (km<-1)
            return 0;
    }   
        
    float tim=pvxl[0]-im;
    float tip=1-tim;
    
    float tjm=pvxl[1]-jm;
    float tjp=1-tjm;
    
    float tkm=pvxl[2]-km;
    float tkp=1-tkm;
  

    double tf1[4],tf2[4];
    tf1[0]=fDataGetValPAD(info, data, im, jm, km, LINEAR, edge);
    tf1[1]=fDataGetValPAD(info, data, ip, jm, km, LINEAR, edge);
    tf1[2]=fDataGetValPAD(info, data, ip, jp, km, LINEAR, edge);
    tf1[3]=fDataGetValPAD(info, data, im, jp, km, LINEAR, edge);
    tf2[0]=fDataGetValPAD(info, data, im, jm, kp, LINEAR, edge);
    tf2[1]=fDataGetValPAD(info, data, ip, jm, kp, LINEAR, edge);
    tf2[2]=fDataGetValPAD(info, data, ip, jp, kp, LINEAR, edge);
    tf2[3]=fDataGetValPAD(info, data, im, jp, kp, LINEAR, edge);

    double tf3[4],t3[4],t4[4];
    int i;
    double tval=0;
    for (i=0;i<4;i++)
        tf3[i]=tkp*tf1[i]+tkm*tf2[i];
    t3[0]=tip;t3[1]=tim;t3[2]=tim;t3[3]=tip;
    t4[0]=tjp;t4[1]=tjp;t4[2]=tjm;t4[3]=tjm;

    for (i=0;i<4;i++)
        tval+=tf3[i]*t3[i]*t4[i];
    (*val) = tval;
    return 1;
}

// Cubic Interpolation
double cubicInterp(double x, const double *y)
{
    double w0=y[1];
    double w1=(-y[0]+y[2])/2.;
    double w2=(2*y[0]-5*y[1]+4*y[2]-y[3])/2.;
    double w3=(-y[0]+3*y[1]-3*y[2]+y[3])/2.;
    return  w0+w1*x+w2*x*x+w3*x*x*x;
}

int getVxlVal_Cubic(float *pvxl, volinfo *info, float* data, INTERPEDGE edge, float* val)
{
	int iL=0, iH=0, jL=0, jH=0, kL=0, kH=0;
	iL=(int)(pvxl[0])-1;
	iH=iL+3;
	jL=(int)(pvxl[1])-1;
	jH=jL+3;
	kL=(int)(pvxl[2])-1;
	kH=kL+3;
	*val = 0 ;
    // check if out of bound
    if (edge == THROWAWAY)
    {
        if (iH>(info->dim[0]-1))
            return 0;
        if (jH>(info->dim[1]-1))
            return 0;
        if (kH>(info->dim[2]-1))
            return 0;
        if (iL<0)
            return 0;
        if (jL<0)
            return 0;
        if (kL<0)
            return 0;
    }
    else
    {
            // check if inbound
        if (iH>(info->dim[0]+1))
            return 0;
        if (jH>(info->dim[1]+1))  // Dominic changed "if (jH>(info->dim[0]+1))"  to  "if (jH>(info->dim[1]+1))"
            return 0;
        if (kH>(info->dim[2]+1))  // Dominic changed "if (kH>(info->dim[0]+1))"  to  "if (kH>(info->dim[2]+1))"
            return 0;
        if (iL<-2)
            return 0;
        if (jL<-2)
            return 0;
        if (kL<-2)
            return 0;   
    }   
	
	
	//Interp
	int i,j,k;
    int im=(int)pvxl[0];
    int ip=im+1;
    float tim=pvxl[0]-im;
    float tip=1-tim;
    
    int jm=(int)pvxl[1];
    int jp=jm+1;
    float tjm=pvxl[1]-jm;
    float tjp=1-tjm;
    
    int km=(int)pvxl[2];
    int kp=km+1;
    float tkm=pvxl[2]-km;
    float tkp=1-tkm;
    
    double fj[4], fi[4], fk[4];
    for (k=km-1;k<=km+2;k++)
    {
        for (i=im-1;i<=im+2;i++)
        {
            for (j=jm-1;j<=jm+2;j++)
                fj[j-jm+1]=fDataGetValPAD(info, data, i, j, k, CUBIC, edge);
            fi[i-im+1]=cubicInterp(tjm, fj);
        }
        fk[k-km+1]=cubicInterp(tim, fi);
    }
    *val = cubicInterp(tkm, fk);
	return 1;
}

int getVxlVals_Cubic(float *pvxl, volinfo *info, float** data, INTERPEDGE edge, float* val)
{
    int iL=0, iH=0, jL=0, jH=0, kL=0, kH=0;
    iL=(int)(pvxl[0])-1;
    iH=iL+3;
    jL=(int)(pvxl[1])-1;
    jH=jL+3;
    kL=(int)(pvxl[2])-1;
    kH=kL+3;
    
    // check if out of bound
    if (edge == THROWAWAY)
    {
        if (iH>(info->dim[0]-1))
            return 0;
        if (jH>(info->dim[1]-1))
            return 0;
        if (kH>(info->dim[2]-1))
            return 0;
        if (iL<0)
            return 0;
        if (jL<0)
            return 0;
        if (kL<0)
            return 0;
    }
    else
    {
            // check if inbound
        if (iH>(info->dim[0]+1))
            return 0;
        if (jH>(info->dim[1]+1))  // Dominic changed "if (jH>(info->dim[0]+1))"  to  "if (jH>(info->dim[1]+1))"
            return 0;
        if (kH>(info->dim[2]+1))  // Dominic changed "if (kH>(info->dim[0]+1))"  to  "if (kH>(info->dim[2]+1))"
            return 0;
        if (iL<-2)
            return 0;
        if (jL<-2)
            return 0;
        if (kL<-2)
            return 0;   
    }   
    
    //Interp
    int i,j,k;
    int im=(int)pvxl[0];
    int ip=im+1;
    float tim=pvxl[0]-im;
    float tip=1-tim;
    
    int jm=(int)pvxl[1];
    int jp=jm+1;
    float tjm=pvxl[1]-jm;
    float tjp=1-tjm;
    
    int km=(int)pvxl[2];
    int kp=km+1;
    float tkm=pvxl[2]-km;
    float tkp=1-tkm;
    
    double fj[4], fi[4], fk[4];
    for (int l=0;l<info->dim[3];l++)
    {
        float *fp = data[l];
        for (k=km-1;k<=km+2;k++)
        {
            for (i=im-1;i<=im+2;i++)
            {
                for (j=jm-1;j<=jm+2;j++)
                    fj[j-jm+1]=fDataGetValPAD(info, fp, i, j, k, CUBIC, edge);
                fi[i-im+1]=cubicInterp(tjm, fj);
            }
            fk[k-km+1]=cubicInterp(tim, fi);
        }
        val[l] = cubicInterp(tkm, fk);
    }
    return 1;
}

/**
 * Get Vxl Value from Morph Space
 */
void fVolGetMorphVal(int n, float* pos, float* val, int* inbound, fvol* vol, fvols* dfield, INTERPMETHOD method, INTERPEDGE edge, int bPositive, int bJacobian)
{
	hvec lphr, vxlr, dlph, lpho, vxlo ;
    // get direcion info
    hvec dirrow,dircol, dirdep;
    hMatGetNormCol(&(vol->info.M_vxl2lph), &dirrow, 0);
    hMatGetNormCol(&(vol->info.M_vxl2lph), &dircol, 1);
    hMatGetNormCol(&(vol->info.M_vxl2lph), &dirdep, 2);
   
    hmat Mlph2vxl_g;
    //hMatInverse(&(dfield->info.M_vxl2lph), &Mlph2vxl_g);
    dMatInvhMat(&(dfield->info.M_vxl2lph), &Mlph2vxl_g);
	hmat mi, mi1;
    dMatInvhMat(&(dfield->info.M_vxl2lph), &mi);
    dMatInvhMat(&(vol->info.M_vxl2lph), &mi1);
	float  dl, dp, dh, fv;
    double jac;
	int i;
	#pragma omp parallel  private(i,fv, dl, dp, dh, vxlr, lphr, dlph, lpho, vxlo,jac)
	#pragma omp for schedule(dynamic)
	for (i=0;i<n;i++)
	{
		hVecInit(&lphr, pos[i*3], pos[i*3+1], pos[i*3+2], 1.);
		hMatMultiplyhVec(&mi, &lphr, &vxlr);
        GetVxlVals(&vxlr, hVecGetfPtr(&dlph), dfield, LINEAR, ZEROPAD, 0);
		hVecInit(&dlph, dl, dp, dh, 0.);
		hVecAdd(&lphr, &dlph, &lpho);
		hMatMultiplyhVec(&mi1, &lpho, &vxlo);
		inbound[i] = GetVxlVal(&vxlo, &fv, vol, method, edge, bPositive);
        // Obtain Jacobian
        if (bJacobian)
        {
            jac = getJac(&lpho, ALL, dfield,  &Mlph2vxl_g, 1., &dirrow, &dircol, &dirdep);
            fv*=jac;
        }
		val[i] = fv;
	}
}

/**
 * Morph fvol into another fvol
 */
void fVolMorph(fvol* vol, fvol* volr, int bsame, fvols *dfield, ucvol *mask, INTERPMETHOD method, INTERPEDGE edge, int bPositive, int bJacobian)
{
	int i,j,k, l;
	hvec vxlr, lphr, dlph, lpho, vxlo ;
	hmat mi;
    // get direcion info
    hvec dirrow,dircol, dirdep;
    hMatGetNormCol(&(vol->info.M_vxl2lph), &dirrow, 0);
    hMatGetNormCol(&(vol->info.M_vxl2lph), &dircol, 1);
    hMatGetNormCol(&(vol->info.M_vxl2lph), &dirdep, 2);
   
    hmat Mlph2vxl_g;
    dMatInvhMat(&(dfield->info.M_vxl2lph), &Mlph2vxl_g);
    
	dMatInvhMat(&(vol->info.M_vxl2lph), &mi);
	float val, dl, dp, dh;
    float *pfdlph= NULL;
    double jac;
	#pragma omp parallel  private(i,j, k, l, val, vxlr, lphr, dlph, lpho, vxlo, pfdlph, jac)
	#pragma omp for schedule(dynamic)
	for (k=0;k<volr->info.dim[2];k++)
	for (i=0;i<volr->info.dim[0];i++)
	for (j=0;j<volr->info.dim[1];j++)
	{
        if (ucVolGetVal(mask, i, j, k))
        {
    		hVecInit(&vxlr, i, j, k, 1.);
    		hMatMultiplyhVec(&(volr->info.M_vxl2lph), &vxlr, &lphr);
    		if (bsame) //volr and dL* are the same vol struct
                fVolsGetVals(dfield, i, j, k, hVecGetfPtr(&dlph));
    		else
    			GetVxlVals(&vxlr, hVecGetfPtr(&dlph), dfield, LINEAR, ZEROPAD, 0);
 
    		hVecSetVal(&dlph, 3, 0);
    		hVecAdd(&lphr, &dlph, &lpho);
    		hMatMultiplyhVec(&mi, &lpho, &vxlo);
    		GetVxlVal(&vxlo, &val, vol, method, edge, bPositive);
            // Obtain Jacobian
            if (bJacobian)
            {
                jac = getJac(&lpho, ALL, dfield,  &Mlph2vxl_g, 1., &dirrow, &dircol, &dirdep);
                val*=jac;
            }
            
    		fVolSetVal(volr, i, j, k, val);
        }
        else
            fVolSetVal(volr, i, j, k, 0);
	}
}

/**
 * Get Vxl Value for a list of points
 */
void fVolGetVxlVal(int n, float * pos, float* val, int *inbound, fvol* vol, hmat* M_r2o, INTERPMETHOD method, INTERPEDGE edge, int bPositive)
{
	hvec lphr, vxlo;
	hmat m1, m2;
    dMatInvhMat(&(vol->info.M_vxl2lph), &m2);
	hMatMultiplyhMat(&m2, M_r2o, &m1);
	float fv;
	int i;
	#pragma omp parallel private(i, lphr, vxlo, fv)
	#pragma omp for schedule(dynamic)
	for (i=0;i<n;i++)
	{
		hVecInit(&lphr, pos[i*3], pos[i*3+1], pos[i*3+2], 1.);
		hMatMultiplyhVec(&m1, &lphr, &vxlo);
		inbound[i] = GetVxlVal(&vxlo, &fv, vol, method, edge, bPositive);
		val[i] = fv;
	}
}

/**
 * Get Vxl Value for a single point from Volumns
 */
int  GetVxlVals(hvec * vxl, float* val, fvols *vols, INTERPMETHOD method, INTERPEDGE edge, int bPositive)
{  
    int l;
    int inbound = 0;
    for (l=0;l<vols->info.dim[3];l++)
          val[l]=0.;
    switch (method)
    {
        case NEAREST:
            inbound = getVxlVals_NN(hVecGetfPtr(vxl),  &(vols->info),  vols->data, val);
            break;
        case LINEAR:
            inbound = getVxlVals_Linear(hVecGetfPtr(vxl), &(vols->info),  vols->data, edge, val);
            break;
        case CUBIC:
            inbound = getVxlVals_Cubic(hVecGetfPtr(vxl), &(vols->info),  vols->data, edge, val);
            break;
    }
    
    //clamping
    for (l=0;l<vols->info.dim[3];l++)
    {
        if ((bPositive) && (val[l]<0))
            val[l]=0.;
    }
    
    return inbound;
}


/**
 * Get Vxl Value for a single point from A volumn
 */
int GetVxlVal(hvec * vxl, float* val, fvol *vol, INTERPMETHOD method, INTERPEDGE edge, int bPositive)
{
	int inbound=0;	
    *val =0;
	switch (method)
	{
		case NEAREST:
			inbound = getVxlVal_NN(hVecGetfPtr(vxl),  &(vol->info),  vol->data, val);
			break;
		case LINEAR:
			inbound = getVxlVal_Linear(hVecGetfPtr(vxl), &(vol->info),  vol->data, edge, val);
			break;
		case CUBIC:
			inbound = getVxlVal_Cubic(hVecGetfPtr(vxl), &(vol->info),  vol->data, edge, val);
			break;
		default:
			*val = 0;
		
	}
	
	//clamping
	if (bPositive)
	{
		if (*val <0)
			*val = 0;
	}
	
	return inbound;
}

/**
 * Resample fvol into another fvol
 */
void fVolResample(fvol* vol, fvol* volr, hmat* M_r2o, INTERPMETHOD method, INTERPEDGE edge, int bPositive)
{
	int i,j,k;
	hvec vxlr, vxlo;
	hmat m1, m2, m3;
	hMatMultiplyhMat(M_r2o, &(volr->info.M_vxl2lph), &m2);
    //hMatPrint(&m2);
    dMatInvhMat(&(vol->info.M_vxl2lph), &m1);
    //hMatPrint(&m1);
	hMatMultiplyhMat(&m1, &m2, &m3);
    //hMatPrint(&m3);
	float val;
	
	#pragma omp parallel  private(i,j, k, val, vxlr, vxlo)
	#pragma omp for schedule(dynamic)
    //for (k=21;k<22;k++)
    //for (i=0;i<1;i++)
    //for (j=14;j<15;j++)
	for (k=0;k<volr->info.dim[2];k++)
	for (i=0;i<volr->info.dim[0];i++)
	for (j=0;j<volr->info.dim[1];j++)
	{
		hVecInit(&vxlr, i, j, k, 1.);
        //hVecPrint(&vxlr);
		hMatMultiplyhVec(&m3, &vxlr, &vxlo);
        //hVecPrint(&vxlo);
		GetVxlVal(&vxlo, &val, vol, method, edge, bPositive);
        //printf("\n val=%f", val);
		fVolSetVal(volr, i, j, k, val);
	}
	
	
} 

void fVolGradUnwarp(fvol* vol, fvol* volr, fvols *dfield, UNWARPTYPE uwtype, 
                   int bJacobian, INTERPMETHOD method, INTERPEDGE edge, int bPositive)
{
    // calcuate mapping from gradint to vxl and vice vesa
    hmat Mvxl2grad, Mgrad2vxl;
    
    hMatMultiplyhMat(&(vol->info.M_pat2grad), &(vol->info.M_vxl2lph), &(Mvxl2grad));
    dMatInvhMat(&(Mvxl2grad), &(Mgrad2vxl));
   
    // get direcion info
    hvec dirrow,dircol, dirdep;
    hMatGetNormCol(&(Mvxl2grad), &dirrow, 0);
    hMatGetNormCol(&(Mvxl2grad), &dircol, 1);
    hMatGetNormCol(&(Mvxl2grad), &dirdep, 2);
    int i,j,k, l;
    hvec vxluw, lphuw, dlph, lpho, vxlo ,dpjvec, vxlg;
    hmat Mlph2vxl_g;
    dMatInvhMat(&(dfield->info.M_vxl2lph), &Mlph2vxl_g);
    double jac=1;
    float val;
    #pragma omp parallel  private(i,j, k, val, vxluw, lphuw, dlph, lpho, vxlo, dpjvec, jac,vxlg)
    #pragma omp for schedule(dynamic)
    //for (k=40;k<41;k++)
    //for (i=100;i<101;i++)
    //for (j=130;j<131;j++)
    for (k=0;k<vol->info.dim[2];k++)
    for (i=0;i<vol->info.dim[0];i++)
    for (j=0;j<vol->info.dim[1];j++)
    {
            hVecInit(&vxlo, i, j, k, 1.);
            hMatMultiplyhVec(&(Mvxl2grad), &vxlo, &lpho);
            hMatMultiplyhVec(&(Mlph2vxl_g), &lpho, &vxlg);
            GetVxlVals(&vxlg, hVecGetfPtr(&dlph), dfield, LINEAR, THROWAWAY, 0);
            // Projection and Get voxel Value
            proj(uwtype, &dlph, &dirrow, &dircol, &dirdep, &dpjvec);
            hVecSetVal(&dpjvec, 3, 0);       
            hVecAdd(&lpho, &dpjvec, &lphuw); 
            hMatMultiplyhVec(&(Mgrad2vxl), &lphuw, &vxluw);
            GetVxlVal(&vxluw, &val, vol, method, edge, bPositive);
            // Obtain Jacobian
            if (bJacobian)
            {
                jac = getJac(&lpho, uwtype, dfield,  &Mlph2vxl_g, 1., &dirrow, &dircol, &dirdep);
                val*=jac;
            }
            fVolSetVal(volr, i, j, k, val);
    }
    
}

void fVolGetValsOnRay(int nv, double *vertice, double *dvec,int np, double *steps, double *vals, int numSmooth, 
  				   	  fvol* vol, hmat* M_r2o, INTERPMETHOD method, INTERPEDGE edge, int bPositive)
{
	float *pos = (float *)malloc(nv*3*sizeof(float));
	float *val = (float *)malloc(nv*sizeof(float));
	int *inb = (int *)malloc(nv*sizeof(int));
	int i,j, k;
	
	// Get Vals
	for (k=0;k<np;k++)
	{
		for (i=0;i<nv;i++)
		for (j=0;j<3;j++)
			pos[3*i+j] = vertice[3*i+j]+steps[k]*dvec[3*i+j];
		
		fVolGetVxlVal(nv, pos, val, inb, vol, M_r2o, method, edge, bPositive);
		
		for (i=0;i<nv;i++)
			vals[np*i+k] = val[i];
	}
	
	//smooth
	float *valt = (float *)malloc(np*sizeof(float));
	for (k=0;k<numSmooth;k++)
	{
		for (i=0;i<nv;i++)
		{
			for (j=0;j<np;j++)
				valt[j] = vals[np*i+j];
			
			vals[np*i] = (2.*valt[0]+valt[1]) /3.;
			
			for (j=1;j<np-1;j++)
				vals[np*i+j] = (valt[j-1]+valt[j]+valt[j+1])/3.;
				
			vals[np*i+np-1] = (2.*valt[np-1]+valt[np-2]) /3.;
			
		}	
	}
	free(pos);
	free(val);
	free(inb);
	free(valt);
	
}

void fVolMorphWithVxlMap(fvol* vol, fvol* volr, fvol *rVi, fvol *rVj,  fvol *rVk, 
						 INTERPMETHOD method, INTERPEDGE edge, int bPositive)
{
	int i, subr[3];
	hvec vxl;
	float val;
	#pragma omp parallel  private(i,subr, vxl, val)
    #pragma omp for schedule(dynamic)
	for (i=0;i<volr->info.vxlcount;i++)
	{
		if (fabs(rVi->data[i])>EPS)
		{
			ind2sub(&(volr->info), i, subr);
			hVecInit(&vxl, rVi->data[i], rVj->data[i], rVk->data[i], 1.);
			if (GetVxlVal(&vxl, &val, vol, method, edge, bPositive))
				volr->data[i] = val;
		}
	}
}
