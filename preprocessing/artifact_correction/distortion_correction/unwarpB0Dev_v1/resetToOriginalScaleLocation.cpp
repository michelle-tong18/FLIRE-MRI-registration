#include <string>
#include <vector>
#include <iostream>
#include "createRescaledMri.h"
#include "interpolation.h" // CTX stuff. fvol

using std::string;
using std::vector;
using std::clog;
using std::endl;


void
resetToOriginalScaleLocation(fvol*& subjectMri,
			     const vector<int>& origSize,
			     const vector<int>& targetShiftValues,
			     const float scaleFac,
			     const bool cubicInterp,
			     const bool rescale) {
  
  
  if(!subjectMri)
    return; // Nothing to do.
  
  const INTERPMETHOD interpMethod ( cubicInterp ? CUBIC : LINEAR ); // NEAREST  LINEAR   CUBIC
  const INTERPEDGE   interpEdge   ( ZEROPAD );                      // MIRROR   ZEROPAD  THROWAWAY
  const int bPositive ( 0 );
  
  if(interpMethod == CUBIC && scaleFac != 1)
    clog << __FILE__ << ":\n Using cubic interpolation when rescaling image...\n" << endl;
  
  const int width  ( subjectMri->info.dim[1] );
  const int height ( subjectMri->info.dim[0] );
  const int depth  ( subjectMri->info.dim[2] );
  //int dim[3] = {height, width, depth};
  
  const int wTMin ( targetShiftValues[0] );
  const int hTMin ( targetShiftValues[1] );
  const int dTMin ( targetShiftValues[2] );
  
  const int wTMax ( targetShiftValues[3] );
  const int hTMax ( targetShiftValues[4] );
  const int dTMax ( targetShiftValues[5] );
  
  const float scaleFacInv ( 1.0/scaleFac );
  const int bcopydataNO  ( 0 );
  
  if( origSize[0] == width  &&  origSize[1] == height  &&  origSize[2] == depth ) {
    // Squeeze images by scale factor; scaleFac >= 1.0
    fvol* subjectMriAux1 ( fVolCopy(subjectMri, bcopydataNO) );
    
    if( scaleFac != 1 ) // Generally, scaleFac >= 1
      for(int x = 0; x < width; ++x)
	for(int y = 0; y < height; ++y)
	  for(int z = 0; z < depth; ++z) {
	    
	    float valS ( 0.0 );
	    
	    hvec vox;
	    hVecInit(&vox, y*scaleFac, x*scaleFac, z*scaleFac, 1.0);
	    
	    GetVxlVal(&vox, &valS, subjectMri, interpMethod, interpEdge, bPositive);
	    
	    if(rescale)
	      valS *= scaleFacInv;
	    
	    fVolSetVal(subjectMriAux1, y, x, z, valS);
	  }
    
    // Translate subject image by shift values
    fvol* subjectMriAux2 ( fVolCopy(subjectMri, bcopydataNO) );
    
    for(int x = wTMin; x < wTMax; ++x)
      for(int y = hTMin; y < hTMax; ++y)
	for(int z = dTMin; z < dTMax; ++z) {
	  
	  float val ( 0.0 );
	  if( scaleFac != 1 ) // Generally, scaleFac >= 1
	    val = fVolGetVal(subjectMriAux1, y-hTMin, x-wTMin, z-dTMin);
	  else
	    val = fVolGetVal(subjectMri, y-hTMin, x-wTMin, z-dTMin);
	  
	  fVolSetVal(subjectMriAux2, y, x, z, val);
	}
    
    
    // Reset incoming MRI
    fVol2volCopy(subjectMriAux2, subjectMri);
    
    fVolDelete(subjectMriAux1);
    fVolDelete(subjectMriAux2);
    
  }
  else { // origSize[0] != width  and/or  origSize[1] != height  and/or  origSize[2] 1= depth
    // Squeeze images by scale factor; scaleFac >= 1.0
    fvol* subjectMriAux1 ( fVolCopy(subjectMri, bcopydataNO) );
    
    if( scaleFac != 1 )
      for(int x = 0; x < width; ++x)
	for(int y = 0; y < height; ++y)
	  for(int z = 0; z < depth; ++z) {
	    
	    float valS ( 0.0 );
	    
	    hvec vox;
	    hVecInit(&vox, y*scaleFac, x*scaleFac, z*scaleFac, 1.0);
	    
	    GetVxlVal(&vox, &valS, subjectMri, interpMethod, interpEdge, bPositive);
	    
	    if(rescale)
	      valS *= scaleFacInv;
	    
	    fVolSetVal(subjectMriAux1, y, x, z, valS);
	  }
    
    // Plonk subject image where it should be in the full original volume.
    
    // Use of createRescaledMri() will result in the LPH center being the same for the new Aux and old non-Aux volumes.
    const bool fillDataNO ( false );
    const float xsize ( scaleFac * subjectMri->info.vxlsize[1] );
    const float ysize ( scaleFac * subjectMri->info.vxlsize[0] );
    const float zsize ( scaleFac * subjectMri->info.vxlsize[2] );
    
    fvol* subjectMriAux2 ( createRescaledMri(origSize[0], origSize[1], origSize[2], xsize, ysize, zsize, subjectMri, fillDataNO) );
    
    //Use of fVolNew() will result in the LPH center being different for the new Aux and old non-Aux volumes.
    //int origDim[3] = {origSize[1], origSize[0], origSize[2]};
    //fvol* subjectMriAux2 ( fVolNew(origDim, &(subjectMri->info.M_vxl2lph), &(subjectMri->info.M_pat2grad), NULL) );
    
    for(int x = wTMin; x < wTMax; ++x)
      for(int y = hTMin; y < hTMax; ++y)
	for(int z = dTMin; z < dTMax; ++z) {
	  
	  float val ( 0.0 );
	  if( scaleFac != 1 ) // Generally, scaleFac >= 1
	    val = fVolGetVal(subjectMriAux1, y-hTMin, x-wTMin, z-dTMin);
	  else
	    val = fVolGetVal(subjectMri, y-hTMin, x-wTMin, z-dTMin);
	  
	  fVolSetVal(subjectMriAux2, y, x, z, val);
	}
    
    fVolDelete(subjectMriAux1);
    // Reset incoming MRI
    fVolDelete(subjectMri);
    subjectMri = subjectMriAux2;
  }
  
}
