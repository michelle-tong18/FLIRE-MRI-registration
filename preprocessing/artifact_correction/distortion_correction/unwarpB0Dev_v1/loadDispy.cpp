#include <iostream>
#include <string>
#include <cstdlib>        // for exit()
#include <cmath>          // for round()
#include "mriClass.h"
#include "imgvol.h"       // CTX stuff. fvol
#include "inputParams.h"
#include "rescaleMri.h"

using std::string;


void
loadDispy(const string& dyFile,
	  const int& numVox,
	  const int& voxStep,
	  const int nvoxNewZbdry,
	  const InputParams& p,
	  double* dispy) {
  
  MriClass dyMriC ( dyFile, nvoxNewZbdry );
  fvol*    dyMri  ( const_cast<fvol*>(dyMriC.getMri()) );
  
  if( p.getResample() ) {  // See also ~/mri/dev/register/src/nonlinearReg.cpp
    //int bcopydataYES  ( 1 );
    //fvol* dyMriOrig ( fVolCopy(dyMri, bcopydataYES) ); // Free in this block.
    
    const float xsOrig (dyMri->info.vxlsize[1]);
    const float ysOrig (dyMri->info.vxlsize[0]);
    const float zsOrig (dyMri->info.vxlsize[2]);
    
    const int   wOrig  ( static_cast<int>(dyMri->info.dim[1]) );
    const int   hOrig  ( static_cast<int>(dyMri->info.dim[0]) );
    const int   dOrig  ( static_cast<int>(dyMri->info.dim[2]) ); // THIS if probably 4 greater than what you think it is!
                                                                      //    SEE mriClass.cpp  p.getNvoxNewZbdry()   mriClass.cpp and solveHessian.cpp
    const int   ddTmp  ( 2*p.getNvoxNewZbdry() );                     // 4  SEE mriClass.cpp  p.getNvoxNewZbdry()   mriClass.cpp and solveHessian.cpp
    
    const float wsOrig ( xsOrig*wOrig );
    const float hsOrig ( ysOrig*hOrig );
    const float dsOrig ( zsOrig*(dOrig-ddTmp) );
    
    const float xs ( p.getXs() );      // Transform input images to this (x,y,z) voxel size (mm^3).
    const float ys ( p.getYs() );      // Usually input images are not isotropic (not cubic voxels). 
    const float zs ( p.getZs() );      // Transform to cubic, estimate defpormation field, and finally transform everything back to orignal voxel size.
    
    const int   w  ( static_cast<int>(round(wsOrig/xs)) );         // Casting float to int rounds down; want the nearest int, not the floor.
    const int   h  ( static_cast<int>(round(hsOrig/ys)) );         // w, h, d are the new numbers of voxels in each dimension.
    const int   d  ( static_cast<int>(round(dsOrig/zs + ddTmp)) );
    
    bool  sincInterp ( false );  // Use trilinear for displacement field interpolation
    rescaleMri(w, h, d, xs, ys, zs, dyMri, sincInterp);       // This routine should have been called resampleMri. Duh.
    
    fVolRescaleI(dyMri, (ysOrig/ys)); // Scale, in preparation for the new voxel y-size
  }
  
  int width  ( dyMri->info.dim[1] );
  int height ( dyMri->info.dim[0] );
  int depth  ( dyMri->info.dim[2] );
  
  int numVoxX =  width/voxStep + ( width%voxStep ? 1 : 0);
  int numVoxY = height/voxStep + (height%voxStep ? 1 : 0);
  int numVoxZ =  depth/voxStep + ( depth%voxStep ? 1 : 0);
  int numVoxp = numVoxX*numVoxY*numVoxZ;
  
  if( numVoxp != numVox ) {
    std::cerr << "numVoxp = " << numVoxp << ",   numVox = " << numVox << std::endl;
    exit(1);
  }
  
  for(int x = 0; x < numVoxX; ++x)
    for(int y = 0; y < numVoxY; ++y)
      for(int z = 0; z < numVoxZ; ++z)  {
	int voxel = x*numVoxY*numVoxZ + y*numVoxZ + z;
	int xa = x*voxStep;
	int ya = y*voxStep;
	int za = z*voxStep;
	dispy[voxel] = fVolGetVal(dyMri, ya, xa, za);
      }
  
}
