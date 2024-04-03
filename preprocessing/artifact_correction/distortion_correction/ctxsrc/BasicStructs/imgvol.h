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
 *  imgvol.h                                         
 *  Created on: 13-Dec-2006 18:22:57                      
 *  Implementation of the Class imgvol       
 ****************************************************/

#ifndef IMGVOL_H_
#define IMGVOL_H_

#include <stdio.h> // For FILE
#include "hvecmat.h"


/* Deal with Windows*/

#ifdef _WIN32
#ifndef DLLTYPE
#define DLLTYPE  __declspec(dllimport)
#endif
#else
#define DLLTYPE
#endif


#ifdef __cplusplus
extern "C" {
#endif



typedef enum 
{
    UINT8,
    INT,
    LONG,
    FLOAT,
    UINT16
} DATATYPE;

/**
 * Image Volume Info 
 */
typedef struct 
{
    /**
     * 1: memory allocated by fVol operation
     * 0: memory alloacted by outside.
     */
    int ballocated;
    /*
     * nbytes: Number of bytes for this datatype
     * totalbytes: Total Bytes
     * volbytecount: byte count per volume
     */
    DATATYPE dtype;
    int nbytes;
    int bytecount;
    int volbytecount; 
    /**
     * ndim: Number of dimension
     * vxlcount: Total number of voxles
     * volvxlcount: voxel count per volume
     * dimension of image volume
     * dim[0]: Number of Rows (Column i direction)
     * dim[1]: Number of Columns (Row j direction)
     * dim[2]: Number of Slices (Slice k direction)
     * dim[3]: Number of Volumes (Slice k direction) 
     *      /k  Slice Direction
     *     /         
     *    -----> j Row direction 
     *    |
     *    | Column direction
     *   \|
     *   i
     */
    int ndim;
    int dim[4];
    int vxlcount;
    int volvxlcount;
    /**
     * voxel size in mm
     * vxlsize[0]: vxlsize in Column i Direction
     * vxlsize[1]: vxlsize in Row j Direction
     * vxlsize[2]: vxlsize in Slcie k Direction
     */
    float vxlsize[3];
    /**
     * Mapping from Patient LPH coordinate to Gradient LPH coordinates
     */
    hmat M_pat2grad;
    /**
     * Mapping from Patient's Voxel coordinate to Patient LPH coordinate
     * (l,p,h)-->(i,j,k)
     */
    hmat M_vxl2lph;

} volinfo;

                                    
DLLTYPE void VolInfoFill(volinfo *info, int ndim, int* dim, hmat * mvxl2lph, hmat * mpat2grad, DATATYPE dtype, int ballocated);
DLLTYPE void VolInfoCopy(volinfo *in, volinfo *inc);
DLLTYPE void VolInfoPrint(volinfo *info);
DLLTYPE void VolInfo2MGH(volinfo *info, FILE *fid, DATATYPE dtype);
DLLTYPE void VolsInfo2MGH(volinfo *info, FILE *fid, DATATYPE dtype); // ADDED BY DOMINIC
DLLTYPE void VolInfoFromMGH(volinfo *info, FILE *fid);
DLLTYPE void VolsInfoFromMGH(volinfo *info, FILE *fid); // ADDED BY DOMINIC
DLLTYPE void VolInfo2CXV(volinfo *info, FILE *fid);
DLLTYPE void VolInfoFromCXV(volinfo *info, FILE *fid);
DLLTYPE void ind2sub(volinfo *info, int ind, int *sub);
DLLTYPE int sub2ind3(volinfo *info, int i, int j, int k);
DLLTYPE int ind4(const int *ps, int i,int j, int k, int l);
DLLTYPE int ind2(const int *ps, int i,int j);
DLLTYPE int ind3(const int *ps, int i,int j, int k);


/**
 * Image Volume with float data
 */
typedef struct 
{
	volinfo info;
	float * data;
} fvol ;

DLLTYPE fvol * fVolNew(int* dim, hmat * mvxl2lph, hmat * mpat2grad, float* img);
DLLTYPE void fVolDelete(fvol* vol);
DLLTYPE fvol* fVolCopy(fvol* vol, int bcopydata);
DLLTYPE char fVol2volCopy(fvol* vol, fvol* volc); // ADDED BY DOMINIC
DLLTYPE void fVolPrint(fvol* vol);
DLLTYPE void fVolFill(fvol *vol, float);
DLLTYPE void fVolRescaleI(fvol *vol, float val);  // ADDED BY DOMINIC
DLLTYPE void fVolSetVal(fvol *vol, int i, int j, int k, float val);
DLLTYPE float fVolGetVal(fvol *vol, int i, int j, int k);
DLLTYPE int fVol2MGH(fvol *vol, char *fn, DATATYPE dtype);
DLLTYPE fvol* fVolFromMGH(char *fn, int bLoadData);
DLLTYPE int fVol2CXV(fvol *vol, char *fn);
DLLTYPE fvol* fVolFromCXV(char *fn, int bLoadData);
DLLTYPE fvol* fVolNewAtlVol();
DLLTYPE int fVolGetVxlCountTH(fvol *vol, int th, int bgreater, unsigned int *ind);
DLLTYPE void fVolVxlScale(fvol *vol, double s, double a);

/**
 * Image Volume with unsigned char data (0-255). Good for segmentation and mask
 */
typedef struct 
{
    volinfo info;
    unsigned char * data;
} ucvol;

DLLTYPE ucvol * ucVolNew(int* dim, hmat * mvxl2lph, hmat * mpat2grad, unsigned char* img);
DLLTYPE void ucVolDelete(ucvol* vol);
DLLTYPE ucvol* ucVolCopy(ucvol* vol, int bcopydata);
DLLTYPE void ucVolPrint(ucvol* vol);
DLLTYPE void ucVolFill(ucvol *vol, unsigned char val);
DLLTYPE void ucVolSetVal(ucvol *vol, int i, int j, int k, unsigned char val);
DLLTYPE unsigned char ucVolGetVal(ucvol *vol, int i, int j, int k);
DLLTYPE int ucVol2MGH(ucvol *vol, char *fn);
DLLTYPE ucvol* ucVolFromMGH(char *fn, int bLoadData);
DLLTYPE ucvol* ucVolNewAtlVol();
DLLTYPE int ucVol2CXV(ucvol *vol, char *fn);
DLLTYPE ucvol* ucVolFromCXV(char *fn, int bLoadData);
DLLTYPE int ucVolGetVxlCountTH(ucvol *vol, int th, int bgreater, unsigned int *ind);

/**
 * Image Volumes with float data. Good for Vector field or multi volume in the same space and shape.
 */
typedef struct 
{
    volinfo info;
    float ** data;
} fvols;

DLLTYPE fvols * fVolsNew(int* dim, hmat * mvxl2lph, hmat * mpat2grad, float* img);
DLLTYPE fvols * fVolsNewVectorField(fvol* dL, fvol* dP, fvol* dH, int bcopydata);
DLLTYPE fvols * fVolsIntenistyAtl(fvol* meanvol, fvol* stdvol, int bcopydata);
DLLTYPE void fVolsDelete(fvols* vols);
DLLTYPE fvols* fVolsCopy(fvols* vols, int bcopydata);
DLLTYPE char fVols2volsCopy(fvols* vols, fvols* volsc);       // ADDED BY DOMINIC
DLLTYPE int fVols2MGH(fvols* vols, char *fn, DATATYPE dtype); // ADDED BY DOMINIC
DLLTYPE fvols* fVolsFromMGH(char *fn, int bLoadData);         // ADDED BY DOMINIC
DLLTYPE void fVolsPrint(fvols* vols);
DLLTYPE void fVolsSetVal(fvols *vols, int i, int j, int k, int l, float val);
DLLTYPE float fVolsGetVal(fvols *vols, int i, int j, int k, int l);
DLLTYPE void fVolsSetVals(fvols *vols, int i, int j, int k, float* val);
DLLTYPE void fVolsGetVals(fvols *vols, int i, int j, int k, float* val);
DLLTYPE int fVols2CXV(fvols *vols, char *fn);
DLLTYPE fvols* fVolsFromCXV(char *fn, int bLoadData);
DLLTYPE fvols* fVolsNewAtlVols(int nvols);

DLLTYPE fvol* fVolCopyFrameFromfVols(fvols* vols, int frame, int bcopydata); // ADDED BY DOMINIC


#ifdef __cplusplus
} /* closing brace for extern "C" */
#endif
#endif 
