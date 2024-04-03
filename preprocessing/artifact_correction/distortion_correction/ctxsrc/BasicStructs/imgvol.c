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
 *  imgvol.c                                         
 *  Created on: 13-Dec-2006 18:22:57                      
 *  Implementation of the Class imgvol       
 ****************************************************/
#ifdef _WIN32
#define DLLTYPE  __declspec(dllexport)
#endif
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
//#include <mathimf.h>
#include <math.h>
#include <unistd.h> // For swab()

#include <errno.h>
extern int errno;

#include "imgvol.h"

int ind4(const int *ps, int i,int j, int k, int l)
{return (l)*ps[2]*ps[1]*ps[0]+(k)*ps[1]*ps[0]+j*ps[1]+i;}

int ind2(const int *ps, int i,int j)
{return i*ps[1]+j;}

int ind3(const int *ps, int i,int j, int k)
{return (k)*ps[1]*ps[0]+(i)*ps[1]+j;}


int GetNumberOfBytes(DATATYPE dtype)
{
    int nb;
    switch (dtype)
    {
        case UINT8:
            nb = sizeof(unsigned char);
            break;
        case INT:
            nb = sizeof(int);
            break;
        case LONG:
            nb = sizeof(long);
            break;
        case FLOAT:
            nb = sizeof(float);
            break;
        case UINT16:
            nb = sizeof(unsigned short);
            break;
        default:
            nb = 0;
            break;
    }
    return nb;
}

void ByteSwap(char * pin, char* pout, int n, DATATYPE dtype)
{
    int nb = GetNumberOfBytes(dtype);
    int i,j;
    //#pragma omp parallel  private(i,j)
    //#pragma omp for schedule(dynamic)
    #pragma omp parallel for
    for (i=0;i<n;i++)
    {
        for (j=0;j<nb;j++)
            pout[i*nb+(nb-1-j)] = pin[i*nb+j];
    }
}

/*
 * Casting dloat to differnet data type back/forth
 */
void CastfData(int tsize, float *fdata, char * data, DATATYPE dtype, int bmgh2vol)
{
  int    i;
  int   *dpc;
  int   *dpi;
  long  *dpl;
  float *dpf;
  unsigned short *dps;
  
  if (bmgh2vol) {
    switch (dtype) {
    case UINT8:
#pragma omp parallel for
      for (i=0;i<tsize;i++)
	fdata[i] = (float)data[i];
      break;
    case INT:
      //int *dpi = (int *)data;
      dpi = (int *)data;
#pragma omp parallel for
      for (i=0;i<tsize;i++)
	fdata[i] = (float)dpi[i];
      break;
    case LONG:
      //long *dpl = (long *)data;
      dpl = (long *)data;
#pragma omp parallel for
      for (i=0;i<tsize;i++)
	fdata[i] = (float)dpl[i];
      break;
    case FLOAT:
      //float *dpf = (float *)data;
      dpf = (float *)data;
#pragma omp parallel for
      for (i=0;i<tsize;i++)
	fdata[i] = dpf[i];
      break;
    case UINT16:
      //unsigned short *dps = (unsigned short *)data;
      dps = (unsigned short *)data;
#pragma omp parallel for
      for (i=0;i<tsize;i++)
	fdata[i] = (float)dps[i];
      break;
    }
  }
  else {
    switch (dtype) {
    case UINT8:
      //int *dpc = (int *)data;
      dpc = (int *)data;
#pragma omp parallel for
      for (i=0;i<tsize;i++)
	dpc[i] = (unsigned char)(fdata[i]);
      break;
    case INT:
      //int *dpi = (int *)data;
      dpi = (int *)data;
#pragma omp parallel for
      for (i=0;i<tsize;i++)
	dpi[i] = (int)(fdata[i]);
      break;
    case LONG:
      //long *dpl = (long *)data;
      dpl = (long *)data;
#pragma omp parallel for
      for (i=0;i<tsize;i++)
	dpl[i] = (long)(fdata[i]);
      break;
    case FLOAT:
      //float *dpf = (float *)data;
      dpf = (float *)data;
#pragma omp parallel for
      for (i=0;i<tsize;i++)
	dpf[i] = (fdata[i]);
      break;
    case UINT16:
      //unsigned short *dps = (unsigned short *)data;
      dps = (unsigned short *)data;
#pragma omp parallel for
      for (i=0;i<tsize;i++)
	dps[i] = (unsigned short)(fdata[i]);
      break;
    }
  }
}

/**
 * Fill in Volume inof header
 */
 
void VolInfoFill(volinfo *info, int ndim, int* dim, hmat * mvxl2lph, hmat * mpat2grad, DATATYPE dtype, int ballocated)
{
    int i;
    info->dtype = dtype;
    info->nbytes = GetNumberOfBytes(dtype);
    info->ndim = ndim;
    info->vxlcount =1;
    for (i=0;i<4;i++)
        info->dim[i] = 1;

    for (i=0;i<info->ndim;i++)
    {
        info->dim[i] = dim[i];
        info->vxlcount = info->vxlcount* dim[i];
    }
    info->volvxlcount = info->vxlcount / info->dim[3];
    info->bytecount = info->vxlcount*info->nbytes;
    info->volbytecount = info->volvxlcount*info->nbytes;
    info->ballocated = ballocated;
    hMatCopy(mvxl2lph, &(info->M_vxl2lph));
    hMatCopy(mpat2grad, &(info->M_pat2grad));
    
    //Get Voxel size
    for (i=0;i<3;i++)
    {
        hvec dirvec;
        hMatGetCol(&(info->M_vxl2lph), &dirvec, i);
        info->vxlsize[i] = hVecNorm(&dirvec);
    }
    
}

void VolInfoCopy(volinfo *in, volinfo *inc)
{
    int i;
    inc->dtype = in->dtype;
    inc->nbytes = in->nbytes;
    inc->ndim = in->ndim;
    inc->vxlcount =in->vxlcount;
    inc->volvxlcount = in->volvxlcount;
    for (i=0;i<in->ndim;i++)
        inc->dim[i] = in->dim[i];
    inc->bytecount = in->bytecount;
    inc->volbytecount = in->volbytecount;
    inc->ballocated = in->ballocated;
    hMatCopy(&(in->M_vxl2lph), &(inc->M_vxl2lph));
    hMatCopy(&(in->M_pat2grad), &(inc->M_pat2grad));
    
    //Get Voxel size
    for (i=0;i<3;i++)
        inc->vxlsize[i] = in->vxlsize[i];
}

void VolInfoPrint(volinfo *info)
{
    printf("\n Image Volume info:");
    printf("\n Number of Rows: %d ", info->dim[0]);
    printf("\n Number of Columns: %d ", info->dim[1]);
    printf("\n Number of Slices: %d ", info->dim[2]);
    printf("\n Number of Volumes: %d ", info->dim[3]);
    printf("\n DataType: %d  Bytes for each data: %d ", info->dtype, info->nbytes);
    printf("\n Total Voxel count: %d Total Bytes Count: %d", info->vxlcount, info->bytecount);
    printf("\n Total Voxel count per Volume: %d Total Bytes Count per Volume: %d", info->volvxlcount, info->volbytecount);
    printf("\n Voxel size in mm :  %6.3f %6.3f %6.3f \n", info->vxlsize[0],info->vxlsize[1],info->vxlsize[2]);
    printf("\n Mapping From Voxel to Patient LPH (M_vxl2lph): ");
    hMatPrint(&(info->M_vxl2lph));
    printf("\n Mapping From Patient LPH to Gradient LPH   (M_pat2grad): ");
    hMatPrint(&(info->M_pat2grad));
    printf("\n Memory allocated: %d ",info->ballocated);
}

/*Fill in MGH Header
 */
void VolInfo2MGH(volinfo *info, FILE *fid, DATATYPE dtype)
{
    int tmp_b[7], tmp_l[7];
    tmp_l[0] = 1;
    tmp_l[1] = info->dim[1];
    tmp_l[2] = info->dim[0];
    tmp_l[3] = info->dim[2];
    tmp_l[4] = 1;
    tmp_l[5] = dtype;
    tmp_l[6] = 0;
    ByteSwap((char*)tmp_l, (char*)tmp_b, 7, INT);
    fwrite((char *)tmp_b, sizeof(int), 7, fid);
    
    short ras_good_l = 1;
    short ras_good_b;
    swab((char *) &ras_good_l, (char *) &ras_good_b, 2); // swap adjacent bytes?
    fwrite((char *)&ras_good_b, sizeof(short), 1, fid);
    
    float ftmp_b[15], ftmp_l[15];
    ftmp_l[0] = info->vxlsize[1];
    ftmp_l[1] = info->vxlsize[0];
    ftmp_l[2] = info->vxlsize[2];
    
    ftmp_l[3] = - hMatGetVal(&(info->M_vxl2lph), 0, 1)/info->vxlsize[1];
    ftmp_l[4] = - hMatGetVal(&(info->M_vxl2lph), 1, 1)/info->vxlsize[1];
    ftmp_l[5] = hMatGetVal(&(info->M_vxl2lph), 2, 1)/info->vxlsize[1];
    ftmp_l[6] = - hMatGetVal(&(info->M_vxl2lph), 0, 0)/info->vxlsize[0];
    ftmp_l[7] = - hMatGetVal(&(info->M_vxl2lph), 1, 0)/info->vxlsize[0];
    ftmp_l[8] = hMatGetVal(&(info->M_vxl2lph), 2, 0)/info->vxlsize[0];
    ftmp_l[9] = - hMatGetVal(&(info->M_vxl2lph), 0, 2)/info->vxlsize[2];
    ftmp_l[10] = - hMatGetVal(&(info->M_vxl2lph), 1, 2)/info->vxlsize[2];
    ftmp_l[11] = hMatGetVal(&(info->M_vxl2lph), 2, 2)/info->vxlsize[2];
    hvec cent_vxl, cent_lph;
    hVecInit(&cent_vxl, (info->dim[0]-1.)/2., (info->dim[1]-1.)/2., (info->dim[2]-1.)/2., 1.);
    hMatMultiplyhVec(&(info->M_vxl2lph), &cent_vxl, &cent_lph);
    ftmp_l[12] = -hVecGetVal(&cent_lph, 0);
    ftmp_l[13] = -hVecGetVal(&cent_lph, 1);
    ftmp_l[14] = hVecGetVal(&cent_lph, 2);
    ByteSwap((char*)ftmp_l, (char*)ftmp_b, 15, FLOAT);
    fwrite((char *)ftmp_b, sizeof(float), 15, fid);
    char unused[194];
    fwrite((char *)unused, sizeof(char), 194, fid); 
}

void VolInfoFromMGH(volinfo *info, FILE *fid)
{
    int tmp_b[7], tmp_l[7];
    fread((char *)tmp_b, sizeof(int), 7, fid);
    ByteSwap((char*)tmp_b, (char*)tmp_l, 7, INT);
    int version = tmp_l[0];
    int width = tmp_l[1];
    int height = tmp_l[2];
    int depth = tmp_l[3];
    int dim[3];
    dim[0]=height;dim[1]=width;dim[2]=depth;
    int nframes = tmp_l[4];
    int dtype  = tmp_l[5];
    int dof = tmp_l[6];
    short ras_good;
    fread((char *)(&ras_good), sizeof(short), 1, fid);
    float ftmp_b[15], ftmp_l[15];
    fread((char *)ftmp_b, sizeof(float), 15, fid);
    ByteSwap((char*)ftmp_b, (char*)ftmp_l, 15, FLOAT);
    float xsize=ftmp_l[0];
    float ysize=ftmp_l[1];
    float zsize=ftmp_l[2];
    float x_r=ftmp_l[3];
    float x_a=ftmp_l[4];
    float x_s=ftmp_l[5];
    float y_r=ftmp_l[6];
    float y_a=ftmp_l[7];
    float y_s=ftmp_l[8];
    float z_r=ftmp_l[9];
    float z_a=ftmp_l[10];
    float z_s=ftmp_l[11];
    float c_r=ftmp_l[12];
    float c_a=ftmp_l[13];
    float c_s=ftmp_l[14];
    char unused[194];
    fread((char *)unused, sizeof(char), 194, fid);
    // Calculate M_vxl2lph;
    hmat Mvxl2lph, Mpat2grad;
    hMatEye(&Mpat2grad);
    hMatEye(&Mvxl2lph);
    hvec rvec, cvec, dvec;
    hVecInit(&rvec, -ysize*y_r, -ysize*y_a, ysize*y_s, 0);
    hMatSetCol(&Mvxl2lph, &rvec, 0);
    hVecInit(&cvec, -xsize*x_r, -xsize*x_a, xsize*x_s, 0);
    hMatSetCol(&Mvxl2lph, &cvec, 1);
    hVecInit(&dvec, -zsize*z_r, -zsize*z_a, zsize*z_s, 0);
    hMatSetCol(&Mvxl2lph, &dvec, 2);
    hvec cent_vxl, cent_lph;
    hVecInit(&cent_vxl, -(height-1.)/2., -(width-1.)/2., -(depth-1.)/2., 1.);
    hVecInit(&cent_lph, -c_r,-c_a, c_s, 0);
    hvec T1, T;
    hMatMultiplyhVec(&Mvxl2lph, &cent_vxl, &T1);
    hVecAdd(&T1, &cent_lph, &T);
    hMatSetCol(&Mvxl2lph, &T, 3);
    
    DATATYPE datatype;
    switch(dtype) {
    case 0: datatype = UINT8;
      break;
    case 1: datatype = INT;
      break;
    case 2: datatype = LONG;
      break;
    case 3: datatype = FLOAT;
      break;
    case 4: datatype = UINT16;
    }
    
    //VolInfoFill(info, 3, dim, &Mvxl2lph, &Mpat2grad, dtype, 0);
    VolInfoFill(info, 3, dim, &Mvxl2lph, &Mpat2grad, datatype, 0);
}

void VolInfo2CXV(volinfo *info, FILE *fid)
{
    unsigned short ver=1;
    fwrite((char *)(&ver), sizeof(unsigned short), 1, fid);
    fwrite((char *)(&(info->dtype)), sizeof(int), 1, fid);
    fwrite((char *)(&(info->ndim)), sizeof(int), 1, fid);
    fwrite((char *)(info->dim), sizeof(int), 4, fid);
    float m[16];
    hMat2fPtr(&(info->M_vxl2lph), m);
    fwrite((char *)(m), sizeof(float), 16, fid);
    hMat2fPtr(&(info->M_pat2grad), m);
    fwrite((char *)(m), sizeof(float), 16, fid);
    char unused[256];
    fwrite((char *)unused, sizeof(char), 256, fid);
}


void VolInfoFromCXV(volinfo *info, FILE *fid)
{
    unsigned short ver;
    fread((char *)(&ver), sizeof(unsigned short), 1, fid);
    int dtype;
    fread((char *)(&dtype), sizeof(int), 1, fid);
    int ndim;
    fread((char *)(&ndim), sizeof(int), 1, fid);
    int dim[4];
    fread((char *)(dim), sizeof(int), 4, fid);
    float m[16];
    fread((char *)(m), sizeof(float), 16, fid);
    hmat Mvxl2lph, Mpat2grad;
    hMatFromfPtr(m, &(Mvxl2lph));
    fread((char *)(m), sizeof(float), 16, fid);
    hMatFromfPtr(m, &(Mpat2grad));
    char unused[256];
    fread((char *)unused, sizeof(char), 256, fid);

    DATATYPE datatype;
    switch(dtype) {
    case 0: datatype = UINT8;
      break;
    case 1: datatype = INT;
      break;
    case 2: datatype = LONG;
      break;
    case 3: datatype = FLOAT;
      break;
    case 4: datatype = UINT16;
    }
    
    //VolInfoFill(info, ndim, dim, &Mvxl2lph, &Mpat2grad, dtype, 0);
    VolInfoFill(info, ndim, dim, &Mvxl2lph, &Mpat2grad, datatype, 0);
}


void ind2sub(volinfo *info, int ind, int *sub)
{
	int ndx;
    ndx = ind;
    if (info->ndim >3)
    {
        sub[3] = (int)floor(ndx / (info->dim[2]*info->dim[1]*info->dim[0])) ;
        ndx = (int)fmod(ndx, info->dim[2]*info->dim[1]*info->dim[0]);
    }
    
    sub[2] = (int)floor(ndx / (info->dim[1]*info->dim[0])) ;
    ndx = (int)fmod(ndx, info->dim[1]*info->dim[0]);
    
    sub[0]= (int)floor(ndx / info->dim[1]);
    ndx = (int)fmod(ndx, info->dim[1]);
     
    sub[1]=ndx;
}

int sub2ind3(volinfo *info, int i, int j, int k)
{
    return ind3(info->dim, i, j, k);
        
}



/**
 * Create a fvol struct
 */
fvol * fVolNew(int* dim, hmat * mvxl2lph, hmat * mpat2grad, float* img)
{
	fvol *vol = (fvol *) malloc(sizeof(fvol)) ;
	VolInfoFill(&(vol->info), 3, dim, mvxl2lph, mpat2grad, FLOAT, 0);
	
	if (img == NULL)
	{
		vol->info.ballocated = 1;
		vol->data = (float *)malloc(vol->info.bytecount);
		if (vol->data == NULL)
			return NULL;
	}
	else
	{
        vol->data = img;
		vol->info.ballocated = 0;
	}
	
	return vol;
}

/**
 * Free fVol and its allocated memory if necessary
 */
void fVolDelete(fvol* vol)
{
	if (vol != NULL)
	{
		if (vol->info.ballocated == 1)
			free(vol->data);
		free(vol);
		vol = NULL;
	}
		
}

/**
 * Create another copy of fVol struct. If bcopydata==1, it  will copy voxel value
 * to the new volume.
 */
fvol* fVolCopy(fvol* vol, int bcopydata)
{
	fvol *volc = (fvol *) malloc(sizeof(fvol)); ;
	int i,j;
	VolInfoCopy(&(vol->info), &(volc->info));
	// allocate memory
	volc->info.ballocated = 1;
	volc->data = (float *)malloc(vol->info.bytecount);
	if (volc->data == NULL)
			return NULL;
			
	if (bcopydata == 1)
	{
		#pragma omp parallel for
		for (i=0;i<vol->info.vxlcount;i++)
			volc->data[i] = vol->data[i];
	}
	
	return volc;
}


/**
 * Copy an fVol into another already existing and allocated fVol of
 * the same voxel dimensions. ADDED BY DOMINIC
 */
char
fVol2volCopy(fvol* vol, fvol* volc) {
  
  if(vol->info.dim[1] != volc->info.dim[1] ||
     vol->info.dim[0] != volc->info.dim[0] ||
     vol->info.dim[2] != volc->info.dim[2] ) {
    printf("\n %s: volumes are different sizes!\n",__FILE__);
    return 1; // ==> error. Caller must check error code.
  }
  
  VolInfoCopy(&(vol->info), &(volc->info));
  int i,j;
#pragma omp parallel for
  for (i=0;i<vol->info.vxlcount;i++)
    volc->data[i] = vol->data[i];
  return 0;   // no error.
}


void fVolPrint(fvol* vol)
{
	VolInfoPrint(&(vol->info));
	printf("\n Data address: %u \n",vol->data);
}

void fVolFill(fvol *vol, float val)
{
    int i;
    int tsize = vol->info.vxlcount;	
	#pragma omp parallel for
	for (i=0;i<tsize;i++)
		vol->data[i] = val;
}


void fVolRescaleI(fvol *vol, float val) { // ADDED BY DOMINIC
  int i;
  int tsize = vol->info.vxlcount;	
  
  if(val != 1.0 ) {
#pragma omp parallel for
    for (i=0;i<tsize;i++)
      vol->data[i] *= val;
  }
}


void fVolSetVal(fvol *vol, int i, int j, int k, float val)
{
	vol->data[sub2ind3(&(vol->info), i, j, k)] = val;
}

float fVolGetVal(fvol *vol, int i, int j, int k)
{
	return vol->data[sub2ind3(&(vol->info), i, j, k)];
}



static int (*myclose)(FILE *);  // declare function pointer

int fVol2MGH(fvol *vol, char *fn, DATATYPE dtype) {

  // static int (*myclose)(FILE *);  // declare function pointer  
  FILE *fid = 0;
  char *ext;
  int  gzipped = 0;
  
  ext = strrchr(fn, '.') ;
  if(ext) {
    char command[512];
    ++ext;
    if( !strcmp(ext, "mgz") || strstr(fn, "mgh.gz") ) { // if mgz, then it is compressed;
                                                        // route stdout to a file
      gzipped = 1;
      myclose = pclose; // assign function pointer for closing
      // pipe writeto "gzip" open
      // pipe can executed under shell and thus understands >
      strcpy(command, "gzip -f -c > ");
      strcat(command, fn);
      errno = 0;
      fid = popen(command, "w");
      if(!fid) {
	errno = 0;
	printf("\n %s: cannot open file %s for writing\n", __FILE__, fn);
	fflush(stdout);
	return EOF;
      }
      if(errno) {
	pclose(fid);
	errno = 0;
	printf("\n %s: gzip had error writing file %s\n", __FILE__, fn);
	fflush(stdout);
	return EOF;
      }
    }
    else if( !strcmp(ext, "mgh") ) {
      fid = fopen(fn, "wb") ;
      myclose = fclose; // assign function pointer for closing
      if(!fid) {
	errno = 0;
	printf("\n %s: cannot open file %s for writing\n", __FILE__, fn);
	fflush(stdout);
	return EOF;
      }
    }
  }
  
  
    size_t sizeT_Written;
    //FILE *fid = fopen(fn, "wb");
    if (fid)
    {
        VolInfo2MGH(&(vol->info), fid, dtype);      
        // Write  data
        int nb = GetNumberOfBytes(dtype);
        int tsize = vol->info.vxlcount;
        char *data_b, *data_l;
        data_b = (char *) malloc(tsize*nb);
        data_l = (char *) malloc(tsize*nb);
        CastfData(tsize, vol->data, data_l, dtype, 0);
        ByteSwap(data_l, data_b, tsize, dtype);
        sizeT_Written = fwrite((char *)data_b, nb, tsize, fid);
	//printf("%s:%d  sizeT_Written = %d\n",__FILE__, __LINE__, sizeT_Written);
	//fflush(stdout);
        free(data_b);
        free(data_l);
    }
    else {
      printf("\n cannot open file: %s\n", fn);
      fflush(stdout);
    }
    //return fclose(fid);
    return myclose(fid);
}



fvol* fVolFromMGH(char *fn, int bLoadData) {
  //static int (*myclose)(FILE *);  // declare function pointer
  FILE *fid = 0;
  char *ext;
  int  gzipped = 0;
  char command[512];
  
  ext = strrchr(fn, '.');
  if (ext) {
    ++ext;
    
    if( !strcmp(ext, "mgz") || strstr(fn, "mgh.gz") ) {
      gzipped = 1;
      myclose = pclose;          // assign function pointer for closing
      strcpy(command, "zcat ");  // OR strcpy(command, "gunzip -c ");
      strcat(command, fn);
      
      errno = 0;
      fid = popen(command, "r"); // popen(): read and write to a unix pipe. E.g., fid = popen("ls -l", "r");
      if(!fid) {
	errno = 0;
        printf("\n%s command: %s\n", __FILE__, command);
      }
      if(errno) {
	pclose(fid);
	errno = 0;
        printf("\n%s command: %s\n", __FILE__, command);
      }
    }
    else if(!strcmp(ext, "mgh")) {
      myclose = fclose;          // assign function pointer for closing
      fid = fopen(fn, "rb");
      if(!fid) {
	errno = 0;
      }
    }
  }
  
  
  //FILE *fid = fopen(fn, "rb");
    fvol *vol = NULL;
    if (fid)
    {
        volinfo info;
        VolInfoFromMGH(&info, fid);
        // Allocate memory and read data
        vol = fVolNew(info.dim, &(info.M_vxl2lph), &(info.M_pat2grad), NULL);
        if (vol && bLoadData)
        {
            int nb = GetNumberOfBytes(info.dtype);
            int tsize = vol->info.vxlcount;
            char *data_b, *data_l;
            data_b = (char *) malloc(tsize*nb);
            data_l = (char *) malloc(tsize*nb);
            fread((char *)data_b, nb, tsize, fid);
            ByteSwap(data_b, data_l, tsize, info.dtype);
            CastfData(tsize, vol->data, data_l, info.dtype, 1);
            free(data_b);
            free(data_l);
        }
        //fclose(fid);
	myclose(fid);
    }
    else
        printf("\n cannot open file: %s", fn);
    return vol;
}


/*
 * New a fvol in atlas space: 256^3 1mm^3 and coronal view
 */
fvol* fVolNewAtlVol()
{
    int dim[3]={256, 256, 256};
    hvec dv;
    hmat mvxl2lph, mpat2grad;
    hMatEye(&mpat2grad);
    hVecInit(&dv, 0., 0., -1., 0.);
    hMatSetCol(&mvxl2lph, &dv, 0);
    hVecInit(&dv, 1., 0., 0., 0.);
    hMatSetCol(&mvxl2lph, &dv, 1);
    hVecInit(&dv, 0., -1., 0., 0.);
    hMatSetCol(&mvxl2lph, &dv, 2);
    hVecInit(&dv, -127.5, 127.5, 127.5, 1.);
    hMatSetCol(&mvxl2lph, &dv, 3);
    return fVolNew(dim, &mvxl2lph, &mpat2grad, 0);
}

int fVol2CXV(fvol *vol, char *fn)
{
    FILE *fid = fopen(fn, "wb");
    if (fid)
    {
        char h[6]={'C', 'T', 'X', 'V', 'O', 'L'};
        fwrite((char *)h, sizeof(char), 6, fid);
        VolInfo2CXV(&(vol->info), fid);
        fwrite((char *)vol->data, vol->info.nbytes, vol->info.vxlcount, fid);
    }
    else
        printf("\n cannot open file: %s", fn);
    return fclose(fid);
}

fvol* fVolFromCXV(char *fn, int bLoadData)
{
    FILE *fid = fopen(fn, "rb");
    fvol *vol = NULL;
    if (fid)
    {
        char H[6]={'C', 'T', 'X', 'V', 'O', 'L'};
        char h[6];
        fread((char *)h, sizeof(char), 6, fid);
        if (strncmp(h, H, 6)!=0)
        {
            printf("\n %s is not a cxv file", fn);
            fclose(fid);
            return vol;
        }
        
        volinfo info;
        VolInfoFromCXV(&info, fid);
        if (info.ndim !=3) 
        {
            printf("\n %s is not a 3-D volume. Number of dimension is:%d", fn, info.ndim);
            fclose(fid);
            return vol;
        }
        
        if (info.dtype != FLOAT) 
        {
            printf("\n %s's datatype is not float. it is type:%d", fn, info.dtype);
            fclose(fid);
            return vol;
        }
        
        // Allocate memory and read data
        vol = fVolNew(info.dim, &(info.M_vxl2lph), &(info.M_pat2grad), NULL);
        if (vol && bLoadData)
            fread((char *)(vol->data), vol->info.nbytes, vol->info.vxlcount, fid);
        
    }
    else
        printf("\n cannot open file: %s", fn);
    fclose(fid);
    return vol;
}

int fVolGetVxlCountTH(fvol *vol, int th, int bgreater, unsigned int *ind)
{
    int i, count=0;
    
    if (bgreater)
    {
        for (i=0;i<vol->info.vxlcount;i++)
        {
            if (vol->data[i]>th)
            {
                ind[count]=i;
                count++;
            }
        }
    }
    else
    {
        for (i=0;i<vol->info.vxlcount;i++)
        {
            if (vol->data[i]<th)
            {
                ind[count]=i;
                count++;
            }
        }
    }
    
    return count;
}


ucvol * ucVolNew(int* dim, hmat * mvxl2lph, hmat * mpat2grad, unsigned char* img)
{
    ucvol *vol = (ucvol *) malloc(sizeof(ucvol)) ;
    VolInfoFill(&(vol->info), 3, dim, mvxl2lph, mpat2grad, UINT8, 0);
    int i,j;
    
    
    if (img == NULL)
    {
        vol->info.ballocated = 1;
        vol->data = (unsigned char *)malloc(vol->info.bytecount);
        if (vol->data == NULL)
            return NULL;
    }
    else
    {
        vol->data = img;
        vol->info.ballocated = 0;
    }
    
    return vol;
}

void ucVolDelete(ucvol* vol)
{
    if (vol != NULL)
    {
        if (vol->info.ballocated == 1)
            free(vol->data);
        free(vol);
        vol = NULL;
    }
}

ucvol* ucVolCopy(ucvol* vol, int bcopydata)
{
    ucvol *volc = (ucvol *) malloc(sizeof(ucvol)); ;
    int i,j;
    VolInfoCopy(&(vol->info), &(volc->info));
    // allocate memory
    volc->info.ballocated = 1;
    volc->data = (unsigned char *)malloc(volc->info.bytecount);
    if (volc->data == NULL)
            return NULL;
            
    if (bcopydata == 1)
    {
        #pragma omp parallel for
        for (i=0;i<vol->info.vxlcount;i++)
            volc->data[i] = vol->data[i];
    }
    
    return volc;
}

void ucVolPrint(ucvol* vol)
{
    VolInfoPrint(&(vol->info));
    printf("\n Data address: %u \n",vol->data);
}

void ucVolFill(ucvol *vol, unsigned char val)
{
    int tsize = vol->info.vxlcount;
    int i;
    
    #pragma omp parallel for
    for (i=0;i<tsize;i++)
        vol->data[i] = val;
}

void ucVolSetVal(ucvol *vol, int i, int j, int k, unsigned char val)
{
   vol->data[sub2ind3(&(vol->info), i, j, k)] = val;
}

unsigned char ucVolGetVal(ucvol *vol, int i, int j, int k)
{
    return vol->data[sub2ind3(&(vol->info), i, j, k)];
}

int ucVol2MGH(ucvol *vol, char *fn)
{
    FILE *fid = fopen(fn, "wb");
    if (fid)
    {
        VolInfo2MGH(&(vol->info), fid, UINT8);
        int nb = GetNumberOfBytes(UINT8);
        int tsize = vol->info.vxlcount; 
        fwrite((char *)vol->data, nb, tsize, fid);     
    }
    else
        printf("\n cannot open file: %s", fn);
    return fclose(fid);
}

ucvol* ucVolFromMGH(char *fn, int bLoadData)
{
    FILE *fid = fopen(fn, "rb");
    ucvol *vol = NULL;
    if (fid)
    {
        volinfo info;
        VolInfoFromMGH(&info, fid);
        if (info.dtype == UINT8)
        {
            // Allocate memory and read data
            vol = ucVolNew(info.dim, &(info.M_vxl2lph), &(info.M_pat2grad), NULL);
            if (vol && bLoadData)
            {
                int nb = GetNumberOfBytes(UINT8);
                int tsize = vol->info.vxlcount;
                fread((char *)vol->data, nb, tsize, fid);
            }
        }
        else
          printf("\n %s is not a UINT8 type file", fn);  
        
        fclose(fid);
    }
    else
        printf("\n cannot open file: %s", fn);
    
    return vol;
}

ucvol* ucVolNewAtlVol()
{
    int dim[3]={256, 256, 256};
    hvec dv;
    hmat mvxl2lph, mpat2grad;
    hMatEye(&mpat2grad);
    hVecInit(&dv, 0., 0., -1., 0.);
    hMatSetCol(&mvxl2lph, &dv, 0);
    hVecInit(&dv, 1., 0., 0., 0.);
    hMatSetCol(&mvxl2lph, &dv, 1);
    hVecInit(&dv, 0., -1., 0., 0.);
    hMatSetCol(&mvxl2lph, &dv, 2);
    hVecInit(&dv, -127.5, 127.5, 127.5, 1.);
    hMatSetCol(&mvxl2lph, &dv, 3);
    return ucVolNew(dim, &mvxl2lph, &mpat2grad, 0);
}



int ucVol2CXV(ucvol *vol, char *fn)
{
    FILE *fid = fopen(fn, "wb");
    if (fid)
    {
        char h[6]={'C', 'T', 'X', 'V', 'O', 'L'};
        fwrite((char *)h, sizeof(char), 6, fid);
        VolInfo2CXV(&(vol->info), fid);
        fwrite((char *)vol->data, vol->info.nbytes, vol->info.vxlcount, fid);
    }
    else
        printf("\n cannot open file: %s", fn);
    return fclose(fid);
}

ucvol* ucVolFromCXV(char *fn, int bLoadData)
{
    FILE *fid = fopen(fn, "rb");
    ucvol *vol = NULL;
    if (fid)
    {
        char H[6]={'C', 'T', 'X', 'V', 'O', 'L'};
        char h[6];
        fread((char *)h, sizeof(char), 6, fid);
        if (strncmp(h, H, 6)!=0)
        {
            printf("\n %s is not a cxv file", fn);
            fclose(fid);
            return vol;
        }
        
        volinfo info;
        VolInfoFromCXV(&info, fid);
        if (info.ndim !=3) 
        {
            printf("\n %s is not a 3-D volume. Number of dimension is:%d", fn, info.ndim);
            fclose(fid);
            return vol;
        }
        
        if (info.dtype != UINT8) 
        {
            printf("\n %s's datatype is not UINT8. it is type:%d", fn, info.dtype);
            fclose(fid);
            return vol;
        }
        
        // Allocate memory and read data
        vol = ucVolNew(info.dim, &(info.M_vxl2lph), &(info.M_pat2grad), NULL);
        if (vol && bLoadData)
            fread((char *)(vol->data), vol->info.nbytes, vol->info.vxlcount, fid);
    }
    else
        printf("\n cannot open file: %s", fn);
    fclose(fid);
    
    return vol;
}

int ucVolGetVxlCountTH(ucvol *vol, int th, int bgreater, unsigned int *ind)
{
    int i, count=0;
    
    if (bgreater)
    {
        for (i=0;i<vol->info.vxlcount;i++)
        {
            if (vol->data[i]>th)
            {
                ind[count]=i;
                count++;
            }
        }
    }
    else
    {
        for (i=0;i<vol->info.vxlcount;i++)
        {
            if (vol->data[i]<th)
            {
                ind[count]=i;
                count++;
            }
        }
    }
    
    return count;
}


//////////////////////////////////////////////////////////////////////////


/*Fill in MGH Header
 */
void VolsInfo2MGH(volinfo *info, FILE *fid, DATATYPE dtype) { // Added by Dominic
    int tmp_b[7], tmp_l[7];
    tmp_l[0] = 1;
    tmp_l[1] = info->dim[1];
    tmp_l[2] = info->dim[0];
    tmp_l[3] = info->dim[2];
    tmp_l[4] = info->dim[3];
    tmp_l[5] = dtype;
    tmp_l[6] = 0;
    ByteSwap((char*)tmp_l, (char*)tmp_b, 7, INT);
    fwrite((char *)tmp_b, sizeof(int), 7, fid);
    
    short ras_good_l = 1;
    short ras_good_b;
    swab((char *) &ras_good_l, (char *) &ras_good_b, 2); // swap adjacent bytes?
    fwrite((char *)&ras_good_b, sizeof(short), 1, fid);
    
    float ftmp_b[15], ftmp_l[15];
    ftmp_l[0] = info->vxlsize[1];
    ftmp_l[1] = info->vxlsize[0];
    ftmp_l[2] = info->vxlsize[2];
    
    ftmp_l[3] = - hMatGetVal(&(info->M_vxl2lph), 0, 1)/info->vxlsize[1];
    ftmp_l[4] = - hMatGetVal(&(info->M_vxl2lph), 1, 1)/info->vxlsize[1];
    ftmp_l[5] = hMatGetVal(&(info->M_vxl2lph), 2, 1)/info->vxlsize[1];
    ftmp_l[6] = - hMatGetVal(&(info->M_vxl2lph), 0, 0)/info->vxlsize[0];
    ftmp_l[7] = - hMatGetVal(&(info->M_vxl2lph), 1, 0)/info->vxlsize[0];
    ftmp_l[8] = hMatGetVal(&(info->M_vxl2lph), 2, 0)/info->vxlsize[0];
    ftmp_l[9] = - hMatGetVal(&(info->M_vxl2lph), 0, 2)/info->vxlsize[2];
    ftmp_l[10] = - hMatGetVal(&(info->M_vxl2lph), 1, 2)/info->vxlsize[2];
    ftmp_l[11] = hMatGetVal(&(info->M_vxl2lph), 2, 2)/info->vxlsize[2];
    hvec cent_vxl, cent_lph;
    hVecInit(&cent_vxl, (info->dim[0]-1.)/2., (info->dim[1]-1.)/2., (info->dim[2]-1.)/2., 1.);
    hMatMultiplyhVec(&(info->M_vxl2lph), &cent_vxl, &cent_lph);
    ftmp_l[12] = -hVecGetVal(&cent_lph, 0);
    ftmp_l[13] = -hVecGetVal(&cent_lph, 1);
    ftmp_l[14] = hVecGetVal(&cent_lph, 2);
    ByteSwap((char*)ftmp_l, (char*)ftmp_b, 15, FLOAT);
    fwrite((char *)ftmp_b, sizeof(float), 15, fid);
    char unused[194];
    fwrite((char *)unused, sizeof(char), 194, fid); 
}


void
VolsInfoFromMGH(volinfo *info, FILE *fid) { // ADDED BY DOMINIC, modeled on VolInfoFromMGH().
    int tmp_b[7], tmp_l[7];
    fread((char *)tmp_b, sizeof(int), 7, fid);
    ByteSwap((char*)tmp_b, (char*)tmp_l, 7, INT);
    int version = tmp_l[0];
    int width = tmp_l[1];
    int height = tmp_l[2];
    int depth = tmp_l[3];
    int nframes = tmp_l[4];
    int dim[4];
    dim[0]=height; // This and
    dim[1]=width;  // this are backwards!!!
    dim[2]=depth;
    dim[3]=nframes;
    int dtype  = tmp_l[5];
    int dof = tmp_l[6];
    short ras_good;
    fread((char *)(&ras_good), sizeof(short), 1, fid);
    float ftmp_b[15], ftmp_l[15];
    fread((char *)ftmp_b, sizeof(float), 15, fid);
    ByteSwap((char*)ftmp_b, (char*)ftmp_l, 15, FLOAT);
    float xsize=ftmp_l[0];
    float ysize=ftmp_l[1];
    float zsize=ftmp_l[2];
    float x_r=ftmp_l[3];
    float x_a=ftmp_l[4];
    float x_s=ftmp_l[5];
    float y_r=ftmp_l[6];
    float y_a=ftmp_l[7];
    float y_s=ftmp_l[8];
    float z_r=ftmp_l[9];
    float z_a=ftmp_l[10];
    float z_s=ftmp_l[11];
    float c_r=ftmp_l[12];
    float c_a=ftmp_l[13];
    float c_s=ftmp_l[14];
    char unused[194];
    fread((char *)unused, sizeof(char), 194, fid);
    // Calculate M_vxl2lph; 
    hmat Mvxl2lph, Mpat2grad;
    hMatEye(&Mpat2grad);
    hMatEye(&Mvxl2lph);
    hvec rvec, cvec, dvec;
    hVecInit(&rvec, -ysize*y_r, -ysize*y_a, ysize*y_s, 0);
    hMatSetCol(&Mvxl2lph, &rvec, 0);
    hVecInit(&cvec, -xsize*x_r, -xsize*x_a, xsize*x_s, 0);
    hMatSetCol(&Mvxl2lph, &cvec, 1);
    hVecInit(&dvec, -zsize*z_r, -zsize*z_a, zsize*z_s, 0);
    hMatSetCol(&Mvxl2lph, &dvec, 2);
    hvec cent_vxl, cent_lph;
    hVecInit(&cent_vxl, -(height-1.)/2., -(width-1.)/2., -(depth-1.)/2., 1.);
    hVecInit(&cent_lph, -c_r,-c_a, c_s, 0);
    hvec T1, T;
    hMatMultiplyhVec(&Mvxl2lph, &cent_vxl, &T1);
    hVecAdd(&T1, &cent_lph, &T);
    hMatSetCol(&Mvxl2lph, &T, 3);
    
    DATATYPE datatype;
    switch(dtype) {
    case 0: datatype = UINT8;
      break;
    case 1: datatype = INT;
      break;
    case 2: datatype = LONG;
      break;
    case 3: datatype = FLOAT;
      break;
    case 4: datatype = UINT16;
    }
    
    //VolInfoFill(info, 3, dim, &Mvxl2lph, &Mpat2grad, dtype, 0);
    VolInfoFill(info, 4, dim, &Mvxl2lph, &Mpat2grad, datatype, 0);
}


//static int (*myclose)(FILE *);  // declare function pointer

int fVols2MGH(fvols* vols, char *fn, DATATYPE dtype) { // ADDED BY DOMINIC
  
  // static int (*myclose)(FILE *);  // declare function pointer  
  FILE *fid = 0;
  char *ext;
  int  gzipped = 0;
  int  frame;
  
  ext = strrchr(fn, '.') ;
  if(ext) {
    char command[512];
    ++ext;
    if( !strcmp(ext, "mgz") || strstr(fn, "mgh.gz") ) { // if mgz, then it is compressed;
                                                        // route stdout to a file
      gzipped = 1;
      myclose = pclose; // assign function pointer for closing
      // pipe writeto "gzip" open
      // pipe can executed under shell and thus understands >
      strcpy(command, "gzip -f -c > ");
      strcat(command, fn);
      errno = 0;
      fid = popen(command, "w");
      if(!fid) {
	errno = 0;
	printf("\n %s: cannot open file %s for writing\n", __FILE__, fn);
	fflush(stdout);
	return EOF;
      }
      if(errno) {
	pclose(fid);
	errno = 0;
	printf("\n %s: gzip had error writing file %s\n", __FILE__, fn);
	fflush(stdout);
	return EOF;
      }
    }
    else if( !strcmp(ext, "mgh") ) {
      fid = fopen(fn, "wb") ;
      myclose = fclose; // assign function pointer for closing
      if(!fid) {
	errno = 0;
	printf("\n %s: cannot open file %s for writing\n", __FILE__, fn);
	fflush(stdout);
	return EOF;
      }
    }
  }
  
  
  size_t sizeT_Written;
  //FILE *fid = fopen(fn, "wb");
  if(fid) {
    VolsInfo2MGH(&(vols->info), fid, dtype);
    
    // Write  data
    int nb = GetNumberOfBytes(dtype);
    int tsize = vols->info.volvxlcount;
    char *data_b, *data_l;
    data_b = (char *) malloc(tsize*nb);
    data_l = (char *) malloc(tsize*nb);
    
    for(frame = 0; frame < vols->info.dim[3]; frame++) {
      CastfData(tsize, vols->data[frame], data_l, dtype, 0);
      ByteSwap(data_l, data_b, tsize, dtype);
      sizeT_Written = fwrite((char *)data_b, nb, tsize, fid);
      //printf("%s:%d  sizeT_Written = %d\n",__FILE__, __LINE__, sizeT_Written);
      //fflush(stdout);
    }
    
    free(data_b);
    free(data_l);
  }
  else {
    printf("\n cannot open file: %s\n", fn);
    fflush(stdout);
  }
  //return fclose(fid);
  return myclose(fid);
}


fvols* fVolsFromMGH(char *fn, int bLoadData) { // ADDED BY DOMINIC
  //static int (*myclose)(FILE *);  // declare function pointer
  FILE *fid = 0;
  char *ext;
  int  gzipped = 0;
  char command[512];
  int  frame;
  
  ext = strrchr(fn, '.');
  if (ext) {
    ++ext;
    
    if( !strcmp(ext, "mgz") || strstr(fn, "mgh.gz") ) {
      gzipped = 1;
      myclose = pclose;          // assign function pointer for closing
      strcpy(command, "zcat ");  // OR strcpy(command, "gunzip -c ");
      strcat(command, fn);
      
      errno = 0;
      fid = popen(command, "r"); // popen(): read and write to a unix pipe. E.g., fid = popen("ls -l", "r");
      if(!fid) {
	errno = 0;
        printf("\n%s command: %s\n", __FILE__, command);
      }
      if(errno) {
	pclose(fid);
	errno = 0;
        printf("\n%s command: %s\n", __FILE__, command);
      }
    }
    else if(!strcmp(ext, "mgh")) {
      myclose = fclose;          // assign function pointer for closing
      fid = fopen(fn, "rb");
      if(!fid) {
	errno = 0;
      }
    }
  }
  
  
  //FILE *fid = fopen(fn, "rb");
  fvols* vols = NULL;
  if(fid) {
    volinfo info;
    VolsInfoFromMGH(&info, fid);
    // Allocate memory and read data
    vols = fVolsNew(info.dim, &(info.M_vxl2lph), &(info.M_pat2grad), NULL);
    if(vols && bLoadData) {
      
      int nb = GetNumberOfBytes(info.dtype);
      int tsize = vols->info.volvxlcount;
      char *data_b, *data_l;
      data_b = (char *) malloc(tsize*nb);
      data_l = (char *) malloc(tsize*nb);
      
      for(frame = 0; frame < vols->info.dim[3]; frame++) {
	//vols->data[i] = img+(i*vols->info.volvxlcount);
	
	fread((char *)data_b, nb, tsize, fid);
	ByteSwap(data_b, data_l, tsize, info.dtype);
	CastfData(tsize, vols->data[frame], data_l, info.dtype, 1);
      }
      
      free(data_b);
      free(data_l);
      
    }
    //fclose(fid);
    myclose(fid);
  }
  else
    printf("\n cannot open file: %s", fn);
  return vols;
}



/*
 * functions for fVols
*/
fvols * fVolsNew(int* dim, hmat * mvxl2lph, hmat * mpat2grad, float* img)
{
    fvols *vols = (fvols *) malloc(sizeof(fvols)) ;
    if (vols)
    {
        VolInfoFill(&(vols->info), 4, dim, mvxl2lph, mpat2grad, FLOAT, 0);
        
        vols->data = (float **)malloc(vols->info.dim[3]*sizeof(float *));
        int i;
        if (vols->data)
        {
            if (img == NULL) 
            {
                vols->info.ballocated = 1;
                for (i=0;i<vols->info.dim[3];i++)
                    vols->data[i] = (float *)malloc(vols->info.volbytecount);
            }
            else
            {
                for (i=0;i<vols->info.dim[3];i++)
                    vols->data[i] = img+(i*vols->info.volvxlcount);
                vols->info.ballocated = 0;
            }
        }
    }
    return vols;
}

fvols * fVolsNewVectorField(fvol* dL, fvol* dP, fvol* dH, int bcopydata)
{
    int i;
    fvols *vols = (fvols *) malloc(sizeof(fvols)) ;
    if (vols)
    {
        int dim[4];
        for (i=0;i<3;i++)
            dim[i] = dL->info.dim[i];
        dim[3] =3;
        VolInfoFill(&(vols->info), 4, dim, &(dL->info.M_vxl2lph), &(dL->info.M_pat2grad), FLOAT, 0);
        
        vols->data = (float **)malloc(vols->info.dim[3]*sizeof(float *));
       
        if (vols->data)
        {
            if (bcopydata)
            {
                vols->info.ballocated = 1;
                for (i=0;i<vols->info.dim[3];i++)
                    vols->data[i] = (float *)malloc(vols->info.volbytecount);
                 memcpy(vols->data[0], dL->data, vols->info.volbytecount);
                 memcpy(vols->data[1], dP->data, vols->info.volbytecount);
                 memcpy(vols->data[2], dH->data, vols->info.volbytecount);
            }
            else
            {
                vols->data[0] = dL->data;
                vols->data[1] = dP->data;
                vols->data[2] = dH->data;
                vols->info.ballocated = 0;
            }  
        }
    }
    return vols;
}

fvols * fVolsIntenistyAtl(fvol* meanvol, fvol* stdvol, int bcopydata)
{
    int i;
    fvols *vols = (fvols *) malloc(sizeof(fvols)) ;
    if (vols)
    {
        int dim[4];
        for (i=0;i<3;i++)
            dim[i] = meanvol->info.dim[i];
        dim[3] =2;
        VolInfoFill(&(vols->info), 4, dim, &(meanvol->info.M_vxl2lph), &(meanvol->info.M_vxl2lph), FLOAT, 0);
        
        vols->data = (float **)malloc(vols->info.dim[3]*sizeof(float *));
       
        if (vols->data)
        {
            if (bcopydata)
            {
                vols->info.ballocated = 1;
                for (i=0;i<vols->info.dim[3];i++)
                    vols->data[i] = (float *)malloc(vols->info.volbytecount);
                 memcpy(vols->data[0], meanvol->data, vols->info.volbytecount);
                 memcpy(vols->data[1], stdvol->data, vols->info.volbytecount);
            }
            else
            {
                vols->data[0] = meanvol->data;
                vols->data[1] = stdvol->data;
                vols->info.ballocated = 0;
            }  
        }
    }
    return vols;
}

void fVolsDelete(fvols* vols)
{
    if (vols != NULL)
    {
        if (vols->info.ballocated == 1)
        {
            for (int i=0;i<vols->info.dim[3];i++)
                free(vols->data[i]);
        }
        free(vols->data);
        free(vols);
        vols = NULL;
    }
}

fvols* fVolsCopy(fvols* vols, int bcopydata)
{
    fvols *volsc = (fvols *) malloc(sizeof(fvols)) ;
    if (vols && volsc)
    {
        VolInfoCopy(&(vols->info), &(volsc->info));
        
        volsc->data = (float **)malloc(volsc->info.dim[3]*sizeof(float *));
        int i,j;
        if (volsc->data)
        {
            volsc->info.ballocated = 1;
            for (i=0;i<volsc->info.dim[3];i++)
                   volsc->data[i] = (float *)malloc(volsc->info.volbytecount);
            if (bcopydata == 1)
            {
                for (i=0;i<volsc->info.dim[3];i++)
                     memcpy(volsc->data[i], vols->data[i], volsc->info.volbytecount);
            }
        }
    }
    return volsc;
}


/**
 * Create an fvol struct form a specified frame of an fvols struct.
 * If bcopydata==1, it  will copy voxel values to the new volume.
 * ADDED BY DOMINIC
 */
fvol* fVolCopyFrameFromfVols(fvols* vols, int frame, int bcopydata) {
  
  fvol* volc = (fvol *)malloc(sizeof(fvol));
  
  VolInfoCopy(&(vols->info), &(volc->info));
  
  // Now, reset some values in volc->info...
  int ndim = vols->info.ndim - 1;   // ndim = 3;
  volc->info.ndim = ndim;
  volc->info.dim[ndim] = 1;
  
  int i,j;
  int vxlcount = 1;
  
  for(i=0;i<ndim;i++)
    vxlcount *= volc->info.dim[i];
  
  volc->info.vxlcount  = vxlcount;
  volc->info.bytecount = vxlcount * volc->info.nbytes;
  
  // Allocate memory
  volc->info.ballocated = 1;
  volc->data = (float *)malloc(volc->info.bytecount);
  if(volc->data == NULL)
    return NULL;
  
  if(bcopydata == 1)
    memcpy(volc->data, vols->data[frame], volc->info.volbytecount);
  
  return volc;
}


/**
 * Copy an fvols into another already existing and allocated fvols of
 * the same voxel dimensions. ADDED BY DOMINIC
 */
char
fVols2volsCopy(fvols* vols, fvols* volsc) {
  
  if(vols->info.dim[1] != volsc->info.dim[1] ||
     vols->info.dim[0] != volsc->info.dim[0] ||
     vols->info.dim[2] != volsc->info.dim[2] ||
     vols->info.dim[3] != volsc->info.dim[3] ) {
    printf("\n %s: volumes are different sizes!\n",__FILE__);
    return 1; // ==> error. Caller must check error code.
  }
  
  VolInfoCopy(&(vols->info), &(volsc->info));
  int i;
  
  for (i = 0; i < volsc->info.dim[3]; i++)
    memcpy(volsc->data[i], vols->data[i], volsc->info.volbytecount);
  
#if 0
  for(frame = 0; frame < vols->info.dim[3]; frame++)
#pragma omp parallel for
    for (i=0;i<vols->info.vxlcount;i++)
      volsc->data[frame][i] = vols->data[frame][i];
#endif
  
  return 0;   // no error.
  
}



void fVolsPrint(fvols* vols)
{
    VolInfoPrint(&(vols->info));
    printf("\n Data address: %u \n",vols->data);
}

void fVolsSetVal(fvols *vols, int i, int j, int k, int l, float val)
{
    vols->data[l][sub2ind3(&(vols->info), i, j, k)] = val;   
}

float fVolsGetVal(fvols *vols, int i, int j, int k, int l)
{
    //float *fp = (float *)(vols->data[l]);
    return vols->data[l][sub2ind3(&(vols->info), i, j, k)];
}

void fVolsSetVals(fvols *vols, int i, int j, int k, float* val)
{
    int ind = sub2ind3(&(vols->info), i, j, k);
    for (int l=0;l<vols->info.dim[3];l++)
         vols->data[l][ind] = val[l]; 
        
}
void fVolsGetVals(fvols *vols, int i, int j, int k, float* val)
{
    int ind = sub2ind3(&(vols->info), i, j, k);
    for (int l=0;l<vols->info.dim[3];l++)
         val[l] = vols->data[l][ind] ;
}

int fVols2CXV(fvols *vols, char *fn)
{
    FILE *fid = fopen(fn, "wb");
    if (fid)
    {
        char h[6]={'C', 'T', 'X', 'V', 'O', 'L'};
        fwrite((char *)h, sizeof(char), 6, fid);
        VolInfo2CXV(&(vols->info), fid);
        for (int i=0;i<vols->info.dim[3];i++)
            fwrite((char *)(vols->data[i]), vols->info.nbytes, vols->info.volvxlcount, fid);
    }
    else
        printf("\n cannot open file: %s", fn);
    return fclose(fid);
}

fvols* fVolsFromCXV(char *fn, int bLoadData)
{
    FILE *fid = fopen(fn, "rb");
    fvols *vols = NULL;
    if (fid)
    {
        char H[6]={'C', 'T', 'X', 'V', 'O', 'L'};
        char h[6];
        fread((char *)h, sizeof(char), 6, fid);
        if (strncmp(h, H, 6)!=0)
        {
            printf("\n %s is not a cxv file", fn);
            fclose(fid);
            return vols;
        }
        
        volinfo info;
        VolInfoFromCXV(&info, fid);
        if (info.ndim != 4) 
        {
            printf("\n %s is not a 4-D volume. Number of dimension is:%d", fn, info.ndim);
            fclose(fid);
            return vols;
        }
        
        if (info.dtype != FLOAT) 
        {
            printf("\n %s's datatype is not FLOAT. it is type:%d", fn, info.dtype);
            fclose(fid);
            return vols;
        }
        
        // Allocate memory and read data
        vols = fVolsNew(info.dim, &(info.M_vxl2lph), &(info.M_pat2grad), NULL);
        if (vols && bLoadData)
        {
            for (int i=0;i<vols->info.dim[3];i++)
                fread((char *)(vols->data[i]), vols->info.nbytes, vols->info.volvxlcount, fid);
        }
    }
    else
        printf("\n cannot open file: %s", fn);
    fclose(fid);
    
    return vols;
}

fvols* fVolsNewAtlVols(int nvols)
{
    int dim[4];
    dim[0]=256;dim[1]=256;dim[2]=256;dim[3]=nvols;
    hvec dv;
    hmat mvxl2lph, mpat2grad;
    hMatEye(&mpat2grad);
    hVecInit(&dv, 0., 0., -1., 0.);
    hMatSetCol(&mvxl2lph, &dv, 0);
    hVecInit(&dv, 1., 0., 0., 0.);
    hMatSetCol(&mvxl2lph, &dv, 1);
    hVecInit(&dv, 0., -1., 0., 0.);
    hMatSetCol(&mvxl2lph, &dv, 2);
    hVecInit(&dv, -127.5, 127.5, 127.5, 1.);
    hMatSetCol(&mvxl2lph, &dv, 3);
    return fVolsNew(dim, &mvxl2lph, &mpat2grad, 0);
}

void fVolVxlScale(fvol *vol, double s, double a)
{
	int i;
	#pragma omp parallel for
	for (i=0;i<vol->info.vxlcount;i++)
		vol->data[i] = vol->data[i]*s+a;
	
}
