// ==================================================
// Copyright (c) 2010 Dominic Holland and Anders Dale
// All rights reserved
// ==================================================

#include <fstream>           // to open files for reading & writing
#include <cmath>             // for round()
#include <ctime>             // for clock()
#include "inputParams.h"
#include "mriClass.h"
#include "intensityScaleImages.h"
#include "mriToFloatArray.h" // FOR TESTING LIBRARY ONLY!!!!!!!
#include "imageInfo.h"       // FOR TESTING LIBRARY ONLY!!!!!!!
#include "computeDeformationField.h"
#include "writeOutput.h"
#include "imgvol.h"          // CTX stuff. fvol

#include <iostream>
#include <string>            // For memset() and strings.
#include <vector>
#include "rescaleMri.h"
#include "resetToOriginalScaleLocation.h"
#include "blurImages.h"

#include "writeImage.h"

using std::string;
using std::vector;

using std::clog;
using std::endl;


void
processSubject(const InputParams& p) {
  
  MriClass forwardMriC ( p.getForwardImageInFileName(), p.getNvoxNewZbdry() );
  fvol*    forwardMri  ( const_cast<fvol*>(forwardMriC.getMri()) );
  
  MriClass reverseMriC ( p.getReverseImageInFileName(), p.getNvoxNewZbdry() );
  fvol*    reverseMri  ( const_cast<fvol*>(reverseMriC.getMri()) );
  
  if(p.getScaleImages())
    intensityScaleImages(forwardMri, reverseMri, p.getImageMax());
    
  if( p.getResample() ) {  // See also ~/mri/dev/register/src/nonlinearReg.cpp
    int bcopydataYES  ( 1 );
    fvol* forwardMriOrig ( fVolCopy(forwardMri, bcopydataYES) ); // Free in this block.
    fvol* reverseMriOrig ( fVolCopy(reverseMri, bcopydataYES) ); // Free in this block.
    
    const float xsOrig (forwardMri->info.vxlsize[1]);
    const float ysOrig (forwardMri->info.vxlsize[0]);
    const float zsOrig (forwardMri->info.vxlsize[2]);
    
    const int   wOrig  ( static_cast<int>(forwardMri->info.dim[1]) );
    const int   hOrig  ( static_cast<int>(forwardMri->info.dim[0]) );
    const int   dOrig  ( static_cast<int>(forwardMri->info.dim[2]) ); // THIS if probably 4 greater than what you think it is!
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
    
    bool  sincInterp ( true );   // <--- Was using this!
    //bool  sincInterp ( false );
    rescaleMri(w, h, d, xs, ys, zs, forwardMri, sincInterp);       // This routine should have been called resampleMri. Duh.
    rescaleMri(w, h, d, xs, ys, zs, reverseMri, sincInterp);
    
    // Use the p that was already built.
    clock_t start ( clock() );
    std::clog << "p.getKernelWidthMax() = " << p.getKernelWidthMax() <<std::endl;
    
    float* dispy ( computeDeformationField(forwardMri, reverseMri, p) );
    
    clock_t secondsUsed ( (clock() - start)/(CLOCKS_PER_SEC) );
    
    std::clog << __FILE__ << ": Duration: " << secondsUsed << " seconds" << std::endl;
    
    //writeOutput(forwardMri, reverseMri, dispy, p);  // Write output in the isotropic resolution used when doing the calculations.
    
    std::clog << "Rescaling displacement field by ys/ysOrig = " << (ys/ysOrig) <<std::endl;    
    int bcopydataNO  ( 0 );
    fvol* dispyMri ( fVolCopy(forwardMri, bcopydataNO) ); // Free in this block.
    for(int x = 0; x < w; ++x)
      for(int y = 0; y < h; ++y)
	for(int z = 0; z < d; ++z) {
	  int   voxel   ( x*(h * d) + y*d + z );
	  float dy      ( dispy[voxel] );
	  //dy *= (ys/ysOrig); // Scale, in preparation for the original voxel y-size
          if (dy>10.0 | dy <-10) dy = 0; // Do nothing.
	  fVolSetVal(dispyMri, y, x, z, dy);
	}
    
    const int kernelWidth ( p.getSmoothDfieldKernelWidth() );  // 3  Odd positive definite number    
    if (kernelWidth > 1) { // Smoothe displacement field
      std::clog << "Smoothing displacement field with kernelWidth = " << kernelWidth <<std::endl;
      //if (kernelWidth%2) kernelWidth += 1;  // Keep it odd.
      int bcopydataNO  ( 0 );
      fvol* dVolBlur ( fVolCopy(dispyMri, bcopydataNO) );        // Free in this block.
      blurImage(dispyMri, dVolBlur, kernelWidth);
      fVol2volCopy(dVolBlur, dispyMri);
    //fVol2MGH(dispyMri, "/space/md10/2/holland/ana/Exam159_Brain/dVolBlur.mgz", FLOAT);
      fVolDelete(dVolBlur);
    }
    fVolRescaleI(dispyMri, (ys/ysOrig)); // Scale, in preparation for the original voxel y-size
    
    std::clog << "Resampling displacement to original resolution" <<std::endl;
    rescaleMri(wOrig, hOrig, dOrig, xsOrig, ysOrig, zsOrig, dispyMri, false);
    delete[] dispy;
    dispy = new float[wOrig*hOrig*dOrig];  // Reuse dispy
    memset(dispy, 0, wOrig*hOrig*dOrig*sizeof(float));
    // This is all a bit of round-the-house way of doing things, but it works.
    for(int x = 0; x < wOrig; ++x)
      for(int y = 0; y < hOrig; ++y)
	for(int z = 0; z < dOrig; ++z) {
	  int   voxel   ( x*(hOrig * dOrig) + y*dOrig + z );
	  dispy[voxel] = fVolGetVal(dispyMri, y, x, z);
	}
    writeOutput(forwardMriOrig, reverseMriOrig, dispy, p);  // Write output in the original resolution of the inout images.
    fVolDelete(forwardMriOrig);
    fVolDelete(reverseMriOrig);
    fVolDelete(dispyMri);
  } else { // No resampling to isotropic: find the deformation field in the native resolution of the inout images directly.
    // Use the p that was already built.
    clock_t start ( clock() );
    std::clog << "p.getKernelWidthMax() = " << p.getKernelWidthMax() <<std::endl;
    
    float* dispy ( computeDeformationField(forwardMri, reverseMri, p) );
    
    clock_t secondsUsed ( (clock() - start)/(CLOCKS_PER_SEC) );
    
    std::clog << __FILE__ << ": Duration: " << secondsUsed << " seconds" << std::endl;
    
    const int kernelWidth ( p.getSmoothDfieldKernelWidth() );  // 3  Odd positive definite number    
    if (kernelWidth > 1) { // Smoothe displacement field
      const int   w  ( static_cast<int>(forwardMri->info.dim[1]) );
      const int   h  ( static_cast<int>(forwardMri->info.dim[0]) );
      const int   d  ( static_cast<int>(forwardMri->info.dim[2]) ); // THIS if probably 4 greater than what you think it is!
      
      int bcopydataNO  ( 0 );
      fvol* dispyMri ( fVolCopy(forwardMri, bcopydataNO) ); // Free in this block.
      for(int x = 0; x < w; ++x)
	for(int y = 0; y < h; ++y)
	  for(int z = 0; z < d; ++z) {
	    int   voxel   ( x*(h * d) + y*d + z );
	    float dy      ( dispy[voxel] );
	    if (dy>10.0 | dy <-10) dy = 0; // Do nothing.
	    fVolSetVal(dispyMri, y, x, z, dy);
	  }
      
      std::clog << "Smoothing displacement field with kernelWidth = " << kernelWidth <<std::endl;
      //if (kernelWidth%2) kernelWidth += 1;  // Keep it odd.
      
      fvol* dVolBlur ( fVolCopy(dispyMri, bcopydataNO) );        // Free in this block.
      blurImage(dispyMri, dVolBlur, kernelWidth);
      fVol2volCopy(dVolBlur, dispyMri);
      //fVol2MGH(dispyMri, "/space/md10/2/holland/ana/Exam159_Brain/dVolBlur.mgz", FLOAT);
      fVolDelete(dVolBlur);
      
      memset(dispy, 0, w*h*d*sizeof(float));      
      // This is all a bit of round-the-house way of doing things, but it works.
      for(int x = 0; x < w; ++x)
	for(int y = 0; y < h; ++y)
	  for(int z = 0; z < d; ++z) {
	    int   voxel   ( x*(h * d) + y*d + z );
	    dispy[voxel] = fVolGetVal(dispyMri, y, x, z);
	  }
      fVolDelete(dispyMri);
    }
    
    writeOutput(forwardMri, reverseMri, dispy, p);
    delete[] dispy;
  }
  
  fVolDelete(forwardMri);
  fVolDelete(reverseMri);
}
