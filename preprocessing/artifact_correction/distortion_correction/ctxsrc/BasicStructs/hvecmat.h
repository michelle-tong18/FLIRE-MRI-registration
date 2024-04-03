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
 *  hvecmat.h                                         
 *  Created on: 13-Dec-2006 08:32:08                      
 *  Implementation of the Class hvecmat       
 *  Original author: Gennan Chen                     
 ****************************************************/

#ifndef HVECMAT_H_
#define HVECMAT_H_
#include <xmmintrin.h>

#ifdef __cplusplus
extern "C" {
#endif


/* Deal with Windows*/
#ifdef _WIN32
#ifndef DLLTYPE
#define DLLTYPE  __declspec(dllimport)
#endif
#else
#define DLLTYPE
#endif





/**
 * Vector and Matrix for Homogeneous coodinates by using SIMD instructions
 */

/**
 * 4 x 1 vector
 */
typedef struct 
{
	__m128 fVect4;
} hvec;

/**
 * 4 X 4 matrix 
 */
typedef struct 
{
	__m128 rows[4];
} hmat;



/* hVec functions*/
DLLTYPE void hVecInit(hvec *v1,float val1, float val2, float val3, float val4); 
DLLTYPE void hVec2fPtr(hvec *v1,float * fptr);
DLLTYPE void hVecSetVal(hvec *v1, int i, float val);
DLLTYPE float hVecGetVal(hvec *v1, int i);
DLLTYPE void hVecFromfPtr(float * fptr, hvec* v1);
DLLTYPE float hVecDot(hvec *v1, hvec *v2);
DLLTYPE float hVecDot3(hvec *v1, hvec *v2);
DLLTYPE void hVecCrossProduct3(hvec *v1, hvec *v2, hvec *vout);
DLLTYPE double hVecComputeAngle(hvec *v1, hvec *v2);
DLLTYPE void hVecPrint(hvec *v);
DLLTYPE float hVecNorm(hvec *v);
DLLTYPE float* hVecGetfPtr(hvec *v);
DLLTYPE void hVecMultiplyfVal(hvec *v, float val, hvec* vo);
DLLTYPE void hVecMultiplyhVec(hvec *v1, hvec *v2, hvec* vo);
DLLTYPE float hVecSum(hvec *v);
DLLTYPE void hVecAdd(hvec *v1, hvec *v2, hvec* vo);
DLLTYPE void hVecSubtracts(hvec *v1, hvec *v2, hvec* vo);

/* hMat functions*/
DLLTYPE void hMat2fPtr(hmat *mat,float * fptr);
DLLTYPE void hMatFromfPtr(float * fptr, hmat *mat);
DLLTYPE void hMatEye(hmat *mat);
DLLTYPE void hMatSetVal(hmat *mat, int i, int j, float val);
DLLTYPE float hMatGetVal(hmat *mat, int i, int j);
DLLTYPE void hMatMultiplyhVec(hmat *mat, hvec *vin, hvec *vout);
DLLTYPE void hMatMultiplyhMat(hmat *matL, hmat *matR, hmat *matout);
DLLTYPE void hMatInverse(hmat *matin, hmat *matout);
DLLTYPE void hMatTranspose(hmat *matin, hmat *matout);
DLLTYPE void hMatRotI(hmat *mat, double radian);
DLLTYPE void hMatRotJ(hmat *mat, double radian);
DLLTYPE void hMatRotK(hmat *mat, double radian);
DLLTYPE void hMatTranslate(hmat *mat, double ti, double tj, double tk);
DLLTYPE void hMatPrint(hmat *mat);
DLLTYPE void hMatSetRow(hmat *mat, hvec *v, int i);
DLLTYPE void hMatSetCol(hmat *mat, hvec *v, int j);
DLLTYPE void hMatGetRow(hmat *mat, hvec *v, int i);
DLLTYPE void hMatGetCol(hmat *mat, hvec *v, int j);
DLLTYPE void hMatGetNormCol(hmat *mat, hvec *v, int j);
DLLTYPE float * hMatGetRowfPtr(hmat *mat, int i);
DLLTYPE int hMat2File(hmat *mat, char * fn);
DLLTYPE int hMatFromFile(char *fn, hmat *mat);
DLLTYPE void hMatCopy(hmat *mi, hmat *mc);

#ifdef __cplusplus
} /* closing brace for extern "C" */
#endif

#endif 
