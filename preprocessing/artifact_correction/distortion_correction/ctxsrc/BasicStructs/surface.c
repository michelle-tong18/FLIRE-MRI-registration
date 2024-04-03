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

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
//#include <mathimf.h>
#include <math.h>
#include <stdbool.h>

#ifdef _WIN32
#define DLLTYPE  __declspec(dllexport)
#endif
#include "hvecmat.h"
#include "imgvol.h"
#include "surface.h"


int mesh2CXS(char * filename, mesh *pmesh)
{
    FILE *fid = fopen(filename, "wb");
    if (fid)
    {
        char h[6]={'C', 'T', 'X', 'S', 'U', 'F'};
        fwrite((char *)h, sizeof(char), 6, fid);
        fwrite((char *)(&(pmesh->numV)), sizeof(unsigned int), 1, fid);
        fwrite((char *)(&(pmesh->numT)), sizeof(unsigned int), 1, fid);
        fwrite((char *)(pmesh->vertices), sizeof(double), 3*pmesh->numV, fid);
        fwrite((char *)(pmesh->triangles), sizeof(unsigned int), 3*pmesh->numT, fid);
    }
    else
        printf("\n cannot open file: %s", filename);
    return fclose(fid);
}


mesh* meshFromCXS(char * filename)
{
    FILE *fid = fopen(filename, "rb");
    mesh *pmesh = NULL;
    if (fid)
    {
        char H[6]={'C', 'T', 'X', 'S', 'U', 'F'};
        char h[6];
        fread((char *)h, sizeof(char), 6, fid);
        if (strncmp(h, H, 6)!=0)
        {
            printf("\n %s is not a cxs file", filename);
            fclose(fid);
            return pmesh;
        }
        
        int numv, numt;
        fread((char *)(&numv), sizeof(unsigned int), 1, fid);
        fread((char *)(&numt), sizeof(unsigned int), 1, fid);
        // Allocate memory and read data
        pmesh = meshNew(numv, numt, NULL, NULL);
        if (pmesh)
        {
            fread((char *)(pmesh->vertices), sizeof(double), 3*pmesh->numV, fid);
            fread((char *)(pmesh->triangles), sizeof(unsigned int), 3*pmesh->numT, fid);
        }
        
    }
    else
        printf("\n cannot open file: %s", filename);
    fclose(fid);
    return pmesh;
}

mesh* meshNew(unsigned int numv, unsigned int numt, double * pvert, unsigned int * ptri)
{
    mesh *pmesh = (mesh *)malloc(sizeof(mesh));
    pmesh->numV = numv;
    if (pvert)
    {
        pmesh->vertices = pvert;
        pmesh->ballocV = 0;
    }
    else
    {
        pmesh->ballocV = 1;
        pmesh->vertices = (double *)calloc(3*pmesh->numV, sizeof(double));
    }
        
    pmesh->numT = numt;
    if (ptri)
    {
        pmesh->triangles = ptri;
        pmesh->ballocT = 0;
    }
    else
    {
        pmesh->ballocT = 1;
        pmesh->triangles = (unsigned int *)calloc(3*pmesh->numT, sizeof(unsigned int));
    }
    
    return pmesh;
    
}

void meshDelete(mesh* pmesh)
{
    if (pmesh->ballocV)
        free(pmesh->vertices);
     pmesh->vertices = NULL;
    if (pmesh->ballocT)
        free(pmesh->triangles);
     pmesh->triangles = NULL; 
    free(pmesh);
    pmesh = NULL;    
}

maskindx* maskindxNew(unsigned int numi,  unsigned int * pind)
{
    maskindx *pmi = (maskindx *) malloc (sizeof(maskindx));
    pmi->numi = numi;
    if (pind)
    {
        pmi->balloc = 0;
        pmi->indx = pind;
    }
    else
    {
        pmi->balloc = 1;
        pmi->indx = (unsigned int *)calloc(pmi->numi, sizeof(unsigned int));
    }
    
    return pmi;
}

void maskindDelete(maskindx* pmi)
{
    if (pmi->balloc)
        free(pmi->indx);
    free(pmi);
    pmi = NULL;   
}

//index from volume threshold
maskindx* maskindxFromfVolTh(fvol *vol, float Il, float Ih)
{
	unsigned int *aindx = (unsigned int *)malloc(sizeof(unsigned int)*(vol->info.vxlcount));
	unsigned int icount=0;
	int i;
	for (i=0;i<vol->info.vxlcount;i++)
	{
		if ((vol->data[i] >Il) & (vol->data[i] <Ih))
		{
			aindx[icount]=i;
			icount++;
		}
		
	}
	
	maskindx *mindx =NULL;
    mindx = maskindxNew(icount, NULL);
        for (i=0;i<icount;i++)
            mindx->indx[i] = aindx[i];
	free(aindx);
	return mindx;
}


maskindx* maskindxFromucVolTh(ucvol *vol, int Il, int Ih)
{
	unsigned int *aindx = (unsigned int *)malloc(sizeof(unsigned int)*(vol->info.vxlcount));
	unsigned int icount=0;
	int i;
	for (i=0;i<vol->info.vxlcount;i++)
	{
		if ((vol->data[i] >Il) & (vol->data[i] <Ih))
		{
			aindx[icount]=i;
			icount++;
		}
		
	}
	
	maskindx *mindx =NULL;
    mindx = maskindxNew(icount, NULL);
    for (i=0;i<icount;i++)
        mindx->indx[i] = aindx[i];
	free(aindx);
	return mindx;
}

void meshGetVert(mesh *pmesh, int i, hvec *v)
{
    hVecInit(v, pmesh->vertices[i*3], pmesh->vertices[i*3+1], pmesh->vertices[i*3+2], 1.);
}

void meshSetVert(mesh *pmesh, int i, hvec *v)
{
    pmesh->vertices[i*3] = hVecGetVal(v, 0);
    pmesh->vertices[i*3+1] = hVecGetVal(v, 1);
    pmesh->vertices[i*3+2] = hVecGetVal(v, 2);
}

//Mesh to vol


// this kind of simple functions are not worthy for parallel....
double getmin(int isize,  float *pval)
{
    double minv=1000000.;
    int i;
    for (int i=0;i<isize;i++)
    {
        if (pval[i]<minv)
            minv=pval[i];
    }
    return minv;

}

double getmax(int isize,  float *pval)
{
    double maxv=-1000000.;
    int i;
    for (i=0;i<isize;i++)
    {
        if (pval[i]>maxv)
            maxv=pval[i];
    }
    return maxv;

}

void cross (const double *pa, const double *pb, double *pc)
{
        pc[0]=pa[1]*pb[2]-pa[2]*pb[1];
        pc[1]=pa[2]*pb[0]-pa[0]*pb[2];
        pc[2]=pa[0]*pb[1]-pa[1]*pb[0];
}

double dot(const double *pa, const double *pb)
{
    double sum=0.;
    for (int i=0;i<3;i++)
        sum+=pa[i]*pb[i];
    return sum;
}

bool isSameSide(const double *p1, const double *p2, const double *a, const double *b)
{
    double cp1[3],cp2[3];
    double dba[3], dp1a[3], dp2a[3];

    for (int i=0;i<3;i++)
    {
        dba[i]=b[i]-a[i];
        dp1a[i]=p1[i]-a[i];
        dp2a[i]=p2[i]-a[i];
    }

    cross(dba, dp1a, cp1);
    cross(dba, dp2a, cp2);
    if (dot(cp1,cp2)>=0)
        return true;
    else
        return false;
}

bool isPntInTriangle(const double *p,const double *a, const double *b,const double *c) 
{
    bool t1=isSameSide(p,a,b,c);
    bool t2=isSameSide(p,b,a,c);
    bool t3=isSameSide(p,c,a,b);

    if (t1 && t2 && t3)
        return true;
    else
        return false;
}

void sort(int n,   double *pval)
{
    for (int i=0; i<n-1; i++) 
    for (int j=0; j<n-1-i; j++)
    {
        if (pval[j+1] < pval[j])
        {  /* compare the two neighbors */
            double tmp = pval[j];         /* swap a[j] and a[j+1]      */
            pval[j] = pval[j+1];
            pval[j+1] = tmp;
        }
    }
}

maskindx * meshGetMaskVol(fvol *vol, mesh *pmesh, hmat *Mreg, int brind, fvol* maskvol, ucvol *mask)
{
    
    hmat Mlph2vxl, mi;
    hMatInverse(&(vol->info.M_vxl2lph), &mi);
    hMatMultiplyhMat(&mi, Mreg, &Mlph2vxl);
    
    //Create bounding boxes for each Triangles and the whole mesh
    double *ptribox=(double *)malloc(sizeof(double)*(pmesh->numT*6));
    double *pmeshbox=(double *)malloc(sizeof(double)*(6));
    pmeshbox[0]=-100000.;
    pmeshbox[1]=100000;
    pmeshbox[2]=-100000.;
    pmeshbox[3]=100000;
    pmeshbox[4]=-100000.;
    pmeshbox[5]=100000.;
    int i,j,k,l;
    hvec lph, vxl;
    float vxlsi[3], vxlsj[3], vxlsk[3];
 
    for (i=0;i<pmesh->numT;i++)
    {
        for (j=0;j<3;j++)
        {
            meshGetVert(pmesh, pmesh->triangles[3*i+j], &lph);
            hMatMultiplyhVec(&Mlph2vxl, &lph, &vxl);
            vxlsi[j]= hVecGetVal(&vxl, 0);
            vxlsj[j]= hVecGetVal(&vxl, 1);
            vxlsk[j]= hVecGetVal(&vxl, 2);
        }
        ptribox[6*i]=getmax(3, vxlsi);
        ptribox[6*i+1]=getmin(3, vxlsi);
        if (pmeshbox[0]<ptribox[6*i])
            pmeshbox[0]=ptribox[6*i];
        if (pmeshbox[1]>ptribox[6*i+1])
            pmeshbox[1]=ptribox[6*i+1];

        ptribox[6*i+2]=getmax(3, vxlsj);
        ptribox[6*i+3]=getmin(3, vxlsj);
        if (pmeshbox[2]<ptribox[6*i+2])
            pmeshbox[2]=ptribox[6*i+2];
        if (pmeshbox[3]>ptribox[6*i+3])
            pmeshbox[3]=ptribox[6*i+3];

        ptribox[6*i+4]=getmax(3, vxlsk);
        ptribox[6*i+5]=getmin(3, vxlsk);
        if (pmeshbox[4]<ptribox[6*i+4])
            pmeshbox[4]=ptribox[6*i+4];
        if (pmeshbox[5]>ptribox[6*i+5])
            pmeshbox[5]=ptribox[6*i+5];
    }
    
    int maxi= ceil(pmeshbox[0]);
    int mini= floor(pmeshbox[1]);
    int maxj= ceil(pmeshbox[2]);
    int minj= floor(pmeshbox[3]);
    int maxk= ceil(pmeshbox[4]);
    int mink= floor(pmeshbox[5]);
    
    if (maxi>vol->info.dim[0]-1)
        maxi = vol->info.dim[0]-1;
    if (mini<0)
        mini=0;
    if (maxj>vol->info.dim[1]-1)
        maxj = vol->info.dim[1]-1;
    if (minj<0)
        minj=0;
    if (maxk>vol->info.dim[2]-1)
        maxk = vol->info.dim[2]-1;
    if (mink<0)
        mink=0;
    //printf("\n Meshbox mini=%i maxi=%i", mini, maxi);
    //printf("\n Meshbox minj=%i maxj=%i", minj, maxj);
    //printf("\n Meshbox mink=%i maxk=%i", mink, maxk);
    
    //scan line
    double jp[4096];
    for (k=0;k<4096;k++)
        jp[k]=0;
    
    int icount =0;
    unsigned int *aindx = (unsigned int *)malloc(sizeof(unsigned int)*(vol->info.vxlcount));
    for (k=mink;k<=maxk;k++)
    for (i=mini;i<=maxi;i++)
    {
        int ic=0;
        int sets=0;
        //mexPrintf("\n k=%d i=%d ",k, i);
        
        for (int t=0;t<pmesh->numT;t++)
        {
            if ( (i>=ptribox[6*t+1]) & (i<=ptribox[6*t]) & (k>=ptribox[6*t+5]) & (k<=ptribox[6*t+4]))
            {
                
                //get vertice indexing
                hvec lpha, vxla, lphb, vxlb,lphc, vxlc;
                meshGetVert(pmesh, pmesh->triangles[3*t], &lpha);
                hMatMultiplyhVec(&Mlph2vxl, &lpha, &vxla);
                double paj = hVecGetVal(&vxla, 1);
                hVecSetVal(&vxla, 1, 0.);
                
                meshGetVert(pmesh, pmesh->triangles[3*t+1], &lphb);
                hMatMultiplyhVec(&Mlph2vxl, &lphb, &vxlb);
                double pbj = hVecGetVal(&vxlb, 1);
                hVecSetVal(&vxlb, 1, 0.);
                
                meshGetVert(pmesh, pmesh->triangles[3*t+2], &lphc);
                hMatMultiplyhVec(&Mlph2vxl, &lphc, &vxlc);
                double pcj = hVecGetVal(&vxlc, 1);
                hVecSetVal(&vxlc, 1, 0.);
                
                double pavxl[3], pbvxl[3], pcvxl[3];
                for (l=0;l<3;l++)
                {
                    pavxl[l]=hVecGetVal(&vxla,l);
                    pbvxl[l]=hVecGetVal(&vxlb,l);
                    pcvxl[l]=hVecGetVal(&vxlc,l);
                }

                double pv[3];
                pv[0]=i;pv[1]=0.;pv[2]=k;
                if (isPntInTriangle(pv, pavxl, pbvxl, pcvxl))
                {
                    pcvxl[1]=pcj;
                    pbvxl[1]=pbj;
                    pavxl[1]=paj;
                    double vba[3];
                    double vca[3];
                    int ii;
                    for (ii=0;ii<3;ii++)
                    {
                        vba[ii]=pbvxl[ii]-pavxl[ii];
                        vca[ii]=pcvxl[ii]-pavxl[ii];
                    }

                    double N[3];
                    cross(vba,vca,N);
                    double D[3];
                    D[0]=0.;D[1]=1.;D[2]=0.;
                    double d=-dot(pavxl,N);
                   // mexPrintf("\n k=%d i=%d t=%d dot(N,D)=%f",k, i, t, dot(N,D));
                    double td=-(d+dot(N,pv))/(dot(N,D));
                    double ip[3];
                    for (ii=0;ii<3;ii++)
                    {
                        ip[ii]=pv[ii]+D[ii]*td;
                    }
                    ic=ic+1;
                    int mic=(int)fmod(ic,2.);
                    if(mic==1)
                        jp[ic-1]=(ip[1]);
                    else
                    {
                        jp[ic-1]=(ip[1]);
                        sets=sets+1;
                    }
                    
                }

            }
        }
        if (sets>0)
        {
            double *jps=(double *)malloc(sizeof(double)*(ic));
            for (l=0;l<ic;l++)
                jps[l]=jp[l];
            sort(ic, jps);

            
            for (int s=0;s<sets;s++)
            {
                int ic_s=2*s;
                int ic_e=2*s+1;
                int jss=(int)ceil(jps[ic_s]);
                int jse=(int)floor(jps[ic_e]);
                for (j=jss;j<=jse;j++)
                {
                    aindx[icount] = ind3(vol->info.dim, i, j, k);                    
                    if (maskvol)
                        maskvol->data[aindx[icount]] = vol->data[aindx[icount]];
                    if (mask)
                        mask->data[aindx[icount]] = 1;
                    icount++;
                }
            }
            free(jps);
        }
        
        for (int ii=0;ii<ic-1;ii++)
            jp[ii]=0;
    }
    
    maskindx *mindx =NULL;
    if (brind)
    {
        mindx = maskindxNew(icount, NULL);
        for (i=0;i<icount;i++)
            mindx->indx[i] = aindx[i];
    }
    
    free(ptribox);
    free(pmeshbox);
    free(aindx);
    return mindx;
}

double meshGetVolSize(mesh *pmesh)
{
    double volsize =0;
    int i;
    int nA, nB, nC;
    hvec crossBC, vecA, vecB, vecC;
    #pragma omp parallel  private(i,vecA, vecB, vecC, crossBC) reduction(+:volsize)
    #pragma omp for schedule(dynamic)
    for (i=0;i<pmesh->numT;i++)
    {
        // get vertice indexing
        meshGetVert(pmesh, pmesh->triangles[3*i], &vecA);
        meshGetVert(pmesh, pmesh->triangles[3*i+1], &vecB);
        meshGetVert(pmesh, pmesh->triangles[3*i+2], &vecC);
        
        hVecCrossProduct3(&vecB, &vecC, &crossBC);
        volsize+=hVecDot3(&vecA, &crossBC);
    }
    return fabs(volsize)/6.;
}


void meshComputeNormals(mesh *pmesh, double *pnormal)
{	
	int i,j,k, ind[3];
	hvec v0, v1, v2, v01, v02, nF, nFn;
	double mag, theta[3];
	for (i=0;i<pmesh->numT;i++)
	{
		// get vertice indexing
		ind[0] = pmesh->triangles[3*i];
		ind[1] = pmesh->triangles[3*i+1];
		ind[2] = pmesh->triangles[3*i+2];
        meshGetVert(pmesh, ind[0], &v0);
        meshGetVert(pmesh, ind[1], &v1);
        meshGetVert(pmesh, ind[2], &v2);
        
        //Calculate normal for this face
		hVecSubtracts(&v1, &v0, &v01);
		hVecSubtracts(&v2, &v0, &v02);
		hVecCrossProduct3(&v01, &v02, &nF);
		mag = sqrt(hVecDot3(&nF, &nF));
		hVecMultiplyfVal(&nF, 1./mag, &nFn);
		
		// Calculate theta
		theta[0] = hVecComputeAngle(&v01, &v02);
		
		hVecSubtracts(&v2, &v1, &v01);
		hVecSubtracts(&v0, &v1, &v02);
		theta[1] = hVecComputeAngle(&v01, &v02);
		
		hVecSubtracts(&v0, &v2, &v01);
		hVecSubtracts(&v1, &v2, &v02);
		theta[2] = hVecComputeAngle(&v01, &v02);
		
		for (j=0;j<3;j++)
		for (k=0;k<3;k++)
			pnormal[ind[j]*3+k] += theta[j]*hVecGetVal(&nFn, k);
	}
	
	for (i=0;i<pmesh->numV;i++)
	{
		mag =0;
		for (j=0;j<3;j++)
			mag  += pnormal[i*3+j]*pnormal[i*3+j];
		mag =sqrt(mag);
		for (j=0;j<3;j++)
			pnormal[i*3+j] = pnormal[i*3+j]/ mag;
	}
}

double meshTriArea(mesh *pmesh, double *parea)
{
	int i;
	hvec v0, v1, v2, v01, v02, nF;
	double mag;
	double tsurface = 0.;
	#pragma omp parallel  private(i,v0, v1, v2,v01, v02, nF, mag) reduction(+:tsurface)
    #pragma omp for schedule(dynamic)
	for (i=0;i<pmesh->numT;i++)
	{
		// get vertice indexing

        meshGetVert(pmesh, pmesh->triangles[3*i], &v0);
        meshGetVert(pmesh, pmesh->triangles[3*i+1], &v1);
        meshGetVert(pmesh, pmesh->triangles[3*i+2], &v2);
        
        //Calculate normal for this face
		hVecSubtracts(&v1, &v0, &v01);
		hVecSubtracts(&v2, &v0, &v02);
		hVecCrossProduct3(&v01, &v02, &nF);
		mag = sqrt(hVecDot3(&nF, &nF));
		
		if (parea)
			parea[i] = mag /2.;
		tsurface += mag /2.;
	}
	
	return tsurface;

}

void meshGetCOM(mesh *pmesh, double* pcom)
{
	int i,j;
	pcom[0]=pcom[1]=pcom[2]=0.;
	for (i=0;i<pmesh->numV;i++)
	{
		for (j=0;j<3;j++)
			pcom[j] += pmesh->vertices[3*i+j];
	}
	for (j=0;j<3;j++)
		pcom[j] = pcom[j] / pmesh->numV;
}


void meshGetMeanEdgeVec(mesh *pmesh, double *pmev)
{
	int i,j, k, ind[3];
	int *vtc = (int *)calloc(pmesh->numV, sizeof (int));
	
	for (i=0;i<pmesh->numV;i++)
	{
		for (j=0;j<3;j++)
			pmev[i*3+j] = 0.;
	}
	
	for (i=0;i<pmesh->numT;i++)
	{
		// get vertice indexing
		ind[0] = pmesh->triangles[3*i];
		ind[1] = pmesh->triangles[3*i+1];
		ind[2] = pmesh->triangles[3*i+2];
        for (j=0;j<3;j++)
        	vtc[ind[j]] = vtc[ind[j]]+2;
		
		for (j=0;j<3;j++)
		{
			pmev[ind[0]*3+j] += pmesh->vertices[ind[1]*3+j]+pmesh->vertices[ind[2]*3+j]-2*pmesh->vertices[ind[0]*3+j];
			pmev[ind[1]*3+j] += pmesh->vertices[ind[0]*3+j]+pmesh->vertices[ind[2]*3+j]-2*pmesh->vertices[ind[1]*3+j];
			pmev[ind[2]*3+j] += pmesh->vertices[ind[0]*3+j]+pmesh->vertices[ind[1]*3+j]-2*pmesh->vertices[ind[2]*3+j];
		}
	}
	
	for (i=0;i<pmesh->numV;i++)
	{
		for (j=0;j<3;j++)
			pmev[i*3+j] = pmev[i*3+j]/ vtc[i];// -  pmesh->vertices[3*i+j];
	}
	
	free(vtc);
}

void meshInflateQ(mesh* pmesh, mesh *pmeshif, int numIts, double stepSize, int Normalize)
{
	int i, j, k;	
	double *pmev = (double *)calloc(3*pmesh->numV, sizeof(double));
	double surfAreaOrig = meshTriArea(pmesh, NULL);
	
	// initialize
	for (i=0;i<pmesh->numV;i++)
	{
		for (j=0;j<3;j++)
			pmeshif->vertices[i*3+j] = pmesh->vertices[i*3+j];
	} 	
		
	for (k=0;k<numIts;k++)
	{
		meshGetMeanEdgeVec(pmeshif, pmev);
		
		for (i=0;i<pmesh->numV;i++)
		{
			for (j=0;j<3;j++)
				pmeshif->vertices[i*3+j] += stepSize * pmev[i*3+j];
		} 	
		
		
	}
	
	if (Normalize)
	{
		double pcom[3];
		double surfAreaNow =  meshTriArea(pmeshif, NULL);
		meshGetCOM(pmeshif, pcom);
		double sf = sqrt(surfAreaOrig/surfAreaNow);
		for (i=0;i<pmesh->numV;i++)
		{
			for (j=0;j<3;j++)
			   pmeshif->vertices[i*3+j] = sf*(pmeshif->vertices[i*3+j]- pcom[j])+pcom[j];
		} 	
	}
	free(pmev);
}

void meshSmoothVertexValue(mesh *pmesh, int numIts, double * pval, double * pvals)
{
	int i,j, k, ind[3];
	int *vtc = (int *)malloc(pmesh->numV*sizeof (int));
	double *pnsum = (double *)malloc(pmesh->numV*sizeof (double));
	for (i=0;i<pmesh->numV;i++)
		pvals[i] = pval[i];
	
	for (k=0;k<numIts;k++)
	{
		// Reset
		for (i=0;i<pmesh->numV;i++)
		{
			vtc[i]=0;
			pnsum[i]=0.;
		}
		
		// Summing neighbors
		for (i=0;i<pmesh->numT;i++)
		{
			// get vertice indexing
			ind[0] = pmesh->triangles[3*i];
			ind[1] = pmesh->triangles[3*i+1];
			ind[2] = pmesh->triangles[3*i+2];
	        for (j=0;j<3;j++)
	        	vtc[ind[j]] = vtc[ind[j]]+1;
			pnsum[ind[0]] += pvals[ind[1]]+pvals[ind[2]];
			pnsum[ind[1]] += pvals[ind[0]]+pvals[ind[2]];
			pnsum[ind[2]] += pvals[ind[0]]+pvals[ind[1]];
		}
		
		// add center 
		for (i=0;i<pmesh->numV;i++)
			pvals[i] = (0.5*pnsum[i]+pvals[i]) / (vtc[i]+1);
	}

	
	free(vtc);
	free(pnsum);
}
