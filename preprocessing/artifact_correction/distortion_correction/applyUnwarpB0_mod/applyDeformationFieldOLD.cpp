#include <cmath>                // for fabs().
#include "interpolation.h"      // CTX stuff. fvol


void
applyDeformationField(fvols* vols, fvol* dVol, const int direction) {
  
  const INTERPMETHOD interpMethod ( CUBIC );   // NEAREST  LINEAR   CUBIC
  const INTERPEDGE   interpEdge   ( ZEROPAD ); // MIRROR   ZEROPAD  THROWAWAY
  const int bPositiveYES ( 1 );
  
  int width   ( vols->info.dim[1] );
  int height  ( vols->info.dim[0] );
  int depth   ( vols->info.dim[2] );
  int nframes ( vols->info.dim[3] );
  
  int dim[3] = {height, width, depth};
  
  float* correctedVol ( new float[width*height*depth] );
  
  for(int frame = 0; frame < nframes; ++frame) {    
    
    for(int x = 0; x < width; ++x)
      for(int y = 0; y < height; ++y)
	for(int z = 0; z < depth; ++z) {
	  
	  float dy  ( fVolGetVal(dVol, y, x, z) );
	  
	  int   ym  ( (y == 0 ? 0 : y-1) );
	  float dym ( fVolGetVal(dVol, ym, x, z) );
	  
	  int   yp  ( (y == height-1 ? height-1 : y+1) );
	  float dyp ( fVolGetVal(dVol, yp, x, z) );
	  
	  const float widtha ( 2.0 ); // Original width around current voxel, ya.
	  
	  float widthb;
	  float yShift;
	  
	  if(direction == 1) {
	    widthb = widtha+dyp-dym;  // New width around current voxel in forward.
	    yShift = y+dy;
	  }
	  else {
	    widthb = widtha-dyp+dym;  // New width around current voxel in reverse.
	    yShift = y-dy;
	  }
	  
	  hvec voxShift;
	  hVecInit(&voxShift, yShift, x, z, 1.0);
	  
	  float val ( 0.0 );
	  GetVxlVal(&voxShift, &val, uncorrectedMri, interpMethod, interpEdge, bPositiveYES);
	  /**
	   * Get Vxl Value for a single point from Volumns
	   */
	  //GetVxlVals(hvec * vxl, float* val, fvols *vols, INTERPMETHOD method, INTERPEDGE edge, int bPositive)	  


	  if(val < 0.0)
	    val = 0.0;
	  
	  float iScale ( fabs(widthb/widtha) );
	  val *= iScale;
	  
	  int voxel ( x*height*depth + y*depth + z );
	  correctedVol[voxel] = val;
	  
	}
    
    // Reset frame segment of incoming vol of frames
    for(int x = 0; x < width; ++x)
      for(int y = 0; y < height; ++y)
	for(int z = 0; z < depth; ++z)  {
	  
	  float val ( correctedVol[voxel] );
	  fVolsSetVals(vols, y, x, z, frame, val);
	}
  }

  delete[] correctedVol;
  
}
