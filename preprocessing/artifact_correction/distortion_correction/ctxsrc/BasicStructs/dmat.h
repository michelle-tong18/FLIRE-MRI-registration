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

#ifndef DMAT_H_
#define DMAT_H_


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


typedef struct 
{
	int rows;
	int cols;
	double *data;
	int ballocated;;
} dmat;


/* mat functions*/
DLLTYPE dmat* dMatNew(int rows, int cols, double *data);
DLLTYPE void dMatDel(dmat* mat);
DLLTYPE dmat* dMatCopy(dmat* mat, int bcopydata);
DLLTYPE dmat* dMatEye(int n);
DLLTYPE void dMatPrint(dmat *mat, char *fn);
DLLTYPE dmat* dMatfromFile(int rows, int cols, char *fn);
DLLTYPE int dMatSameSize(dmat * m1, dmat *m2);

DLLTYPE void dMatSetVal(dmat * mat, int i, int j, double val);
DLLTYPE double dMatGetVal(dmat * mat, int i, int j);
DLLTYPE dmat* dMatGetRow(dmat * mat, int i);
DLLTYPE dmat* dMatGetCol(dmat * mat, int i);

// blas 1
DLLTYPE dmat * dMatAdd(dmat *m1, dmat *m2, double a);
DLLTYPE dmat * dMatScal(dmat *m,  double a);
DLLTYPE double dMatNorm2(dmat *m, int n);
DLLTYPE double dMatDot(dmat *m1, dmat *m2,  int n);

// blas 2


// blas 3
DLLTYPE dmat* dMatMultiply(dmat * matL, dmat *matR, double alpha, double beta);

// lapack
DLLTYPE dmat* dMatInv(dmat * mat);
DLLTYPE dmat* dMatInvhMat(hmat *m, hmat *mi);

//other
DLLTYPE dmat* dMatCross3(dmat *m1, dmat *m2);
DLLTYPE double dMatComputeAngle(dmat *m1, dmat *m2, int n);
DLLTYPE dmat* dMatRotI(double radian);
DLLTYPE dmat* dMatRotJ(double radian);
DLLTYPE dmat* dMatRotK(double radian);
DLLTYPE dmat* dMatTranslate(double ti, double tj, double tk);
#ifdef __cplusplus
} /* closing brace for extern "C" */
#endif

#endif
