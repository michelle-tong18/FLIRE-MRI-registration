// ==================================================
// Copyright (c) 2010 Dominic Holland and Anders Dale
// All rights reserved
// ==================================================

#include <cmath>                // for fabs().
#include "interpolation.h"      // CTX stuff. fvol

#define SQR(A) ((A)*(A))
#define WRITE_JACOBIAN 0


void
applyDeformationField(fvols* vols, fvol* dVol, const int direction, const bool jic) {
  
  const INTERPMETHOD interpMethod ( CUBIC );   // NEAREST  LINEAR   CUBIC
  const INTERPEDGE   interpEdge   ( ZEROPAD ); // MIRROR   ZEROPAD  THROWAWAY
  const int bPositive ( 0 );
  
  int width   ( vols->info.dim[1] );
  int height  ( vols->info.dim[0] );
  int depth   ( vols->info.dim[2] );
  int nframes ( vols->info.dim[3] );
  
  int bcopydata ( 0 );
  fvols* volsC ( fVolsCopy(vols, bcopydata) );
  
#if WRITE_JACOBIAN
  fvol* jacobianVol  ( fVolCopy(dVol, bcopydata) );
  fvol* jacobianVoly ( fVolCopy(dVol, bcopydata) );
  fvol* jacobianVolx ( fVolCopy(dVol, bcopydata) );
  fvol* jacobianVolz ( fVolCopy(dVol, bcopydata) );
#endif
  
  for(int x = 0; x < width; ++x)
    for(int y = 0; y < height; ++y)
      for(int z = 0; z < depth; ++z) {
	
	float dy  ( fVolGetVal(dVol, y, x, z) );
	
	int   ym  ( (y == 0 ? 0 : y-1) );
	float dym ( fVolGetVal(dVol, ym, x, z) );
	
	int   yp  ( (y == height-1 ? height-1 : y+1) );
	float dyp ( fVolGetVal(dVol, yp, x, z) );

	int   xm  ( (x == 0 ? 0 : x-1) );
	float dxm ( fVolGetVal(dVol, y, xm, z) );
	
	int   xp  ( (x == width-1 ? width-1 : x+1) );
	float dxp ( fVolGetVal(dVol, y, xp, z) );
	
	int   zm  ( (z == 0 ? 0 : z-1) );
	float dzm ( fVolGetVal(dVol, y, x, zm) );
	
	int   zp  ( (z == depth-1 ? depth-1 : z+1) );
	float dzp ( fVolGetVal(dVol, y, x, zp) );
	
	const float widtha ( 2.0 ); // Original width around current voxel, ya.
	
	float widthb;
	float yShift;
	
	float widthbx ( 2.0 );
	float widthbz ( 2.0 );
	
	if(direction == 1) {
	  widthb = widtha+dyp-dym;  // New width around current voxel in forward.
	  yShift = y+dy;
	  widthbx = widtha+dxp-dxm;
	  widthbz = widtha+dzp-dzm;
	}
	else {
	  widthb = widtha-dyp+dym;  // New width around current voxel in reverse.
	  yShift = y-dy;
	  widthbx = widtha-dxp+dxm;
	  widthbz = widtha-dzp+dzm;
	}
	
	hvec voxShift;
	hVecInit(&voxShift, yShift, x, z, 1.0);
	
	//float val ( 0.0 );
	//GetVxlVal(&voxShift, &val, uncorrectedMri, interpMethod, interpEdge, bPositive);
	
	float vals[nframes];        // Initialize to all zero
	for(int i=0; i<nframes; ++i)
	  vals[i] = 0.0;
	
	GetVxlVals(&voxShift, vals, vols, interpMethod, interpEdge, bPositive);
        
        if( jic ) { // Jacobian intensity correction.
	  float iScale ( fabs(widthb/widtha) );
	  for(int i=0; i<nframes; ++i)
	    vals[i] *= iScale;
	}
	
	fVolsSetVals(volsC, y, x, z, vals);
	
#if WRITE_JACOBIAN
	float iScalex ( fabs(widthbx/widtha) );
	float iScalez ( fabs(widthbz/widtha) );
	float iScale ( fabs(widthb/widtha) );
	fVolSetVal(jacobianVoly, y, x, z, iScale-1.0);
	fVolSetVal(jacobianVolx, y, x, z, iScalex-1.0);
	fVolSetVal(jacobianVolz, y, x, z, iScalez-1.0);
	fVolSetVal(jacobianVol,  y, x, z, sqrt(SQR(iScale-1.0) + SQR(iScalex-1.0) + SQR(iScalez-1.0)));
#endif
      }
  
  fVols2volsCopy(volsC, vols);
#if WRITE_JACOBIAN
  fVol2MGH(jacobianVoly, "/space/md10/2/holland/b0unwarp/lucinehG_3T_20090318.161256_1/unwarpSEaxial/jyRev.mgz", FLOAT);
  fVol2MGH(jacobianVolx, "/space/md10/2/holland/b0unwarp/lucinehG_3T_20090318.161256_1/unwarpSEaxial/jxRev.mgz", FLOAT);
  fVol2MGH(jacobianVolz, "/space/md10/2/holland/b0unwarp/lucinehG_3T_20090318.161256_1/unwarpSEaxial/jzRev.mgz", FLOAT);
  fVol2MGH(jacobianVol,  "/space/md10/2/holland/b0unwarp/lucinehG_3T_20090318.161256_1/unwarpSEaxial/jRev.mgz",  FLOAT);
  fVolDelete(jacobianVol);
  fVolDelete(jacobianVoly);
  fVolDelete(jacobianVolx);
  fVolDelete(jacobianVolz);
#endif  
  fVolsDelete(volsC);
  
}

#undef SQR
