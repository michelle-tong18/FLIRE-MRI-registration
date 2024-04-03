#include <string>
#include <iostream>
#include <cmath>
#include "rescaleMri.h"

#include "interpolation.h" // CTX stuff. fvol

using std::string;
using std::clog;
using std::endl;


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// This first one makes volS the same size as volT, i.e., same number
// of voxels and same voxel sizes as in volT, but keeps its own
// direction cosines x_l, x_p, x_h, ditto y & z, and keeps its own LPH
// center c_l, c_p, and c_h.
// 
// The second, overloaded, version below just gives incoming mri the
// specified voxel slzes and number (and again keeps its own
// direction cosines and LPH center).
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


bool
rescaleMri(fvol* volT, fvol*& volS, const bool useCubic) {
  
  if(!volS || !volT)
    return false; // Nothing to do.
  
  int widthT   ( volT->info.dim[1] );
  int heightT  ( volT->info.dim[0] );
  int depthT   ( volT->info.dim[2] );
  
  float xsizeT ( volT->info.vxlsize[1] );
  float ysizeT ( volT->info.vxlsize[0] );
  float zsizeT ( volT->info.vxlsize[2] );
  
  return ( rescaleMri(widthT, heightT, depthT, xsizeT, ysizeT, zsizeT, volS, useCubic) );
}


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////// O V E R L O A D I N G ///////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////


bool
rescaleMri(const int   widthC, const int   heightC, const int   depthC,
	   const float xsizeC, const float ysizeC,  const float zsizeC, fvol*& volS, const bool useCubic) {
  
  if(!volS)
    return false; // Nothing to do.
  
  int widthS  ( volS->info.dim[1] );
  int heightS ( volS->info.dim[0] );
  int depthS  ( volS->info.dim[2] );
  int frames  ( volS->info.dim[3] );
  
  float xsizeS ( volS->info.vxlsize[1] );
  float ysizeS ( volS->info.vxlsize[0] );
  float zsizeS ( volS->info.vxlsize[2] );
  
  if( widthS  == widthC  &&
      heightS == heightC &&
      depthS  == depthC  &&
      xsizeS  == xsizeC  &&
      ysizeS  == ysizeC  &&
      zsizeS  == zsizeC )
    return false; // Nothing to do.
  
  clog<<__FILE__<<":\n"
      <<"widthS  = "<<widthS <<"\t"<<"widthC  = "<<widthC <<"\n"
      <<"heightS = "<<heightS<<"\t"<<"heightC = "<<heightC<<"\n"
      <<"depthS  = "<<depthS <<"\t"<<"depthC  = "<<depthC <<"\n\n"
      <<"xsizeS  = "<<xsizeS <<"\t"<<"xsizeC  = "<<xsizeC <<"\n"
      <<"ysizeS  = "<<ysizeS <<"\t"<<"ysizeC  = "<<ysizeC <<"\n"
      <<"zsizeS  = "<<zsizeS <<"\t"<<"zsizeC  = "<<zsizeC <<endl;
  
  if(widthC  != widthS)  clog<<"widthC-widthS   = "<<widthC-widthS  <<endl;
  if(heightC != heightS) clog<<"heightC-heightS = "<<heightC-heightS<<endl;
  if(depthC  != depthS)  clog<<"depthC-depthS   = "<<depthC-depthS  <<endl;
  
  if(xsizeC != xsizeS)   clog<<"xsizeC-xsizeS   = "<<xsizeC-xsizeS<<endl;
  if(ysizeC != ysizeS)   clog<<"ysizeC-ysizeS   = "<<ysizeC-ysizeS<<endl;
  if(zsizeC != zsizeS)   clog<<"zsizeC-zsizeS   = "<<zsizeC-zsizeS<<endl;
  
  const float eps ( 1.0e-2);
  if( widthC == widthS && heightC == heightS && depthC == depthS )
    if( fabs(xsizeC-xsizeS) < eps && fabs(ysizeC-ysizeS) < eps && fabs(zsizeC-zsizeS) < eps )
      return false; // Nothing to do.
  
  clog<<__FILE__<<"\nRescaling image...\n"<<endl;
  if(useCubic)
    clog<<__FILE__<<"\nUsing cubic interpolation...\n"<<endl;
  
  ///////////////////////////////////////////////
  // Create a new volA to hold the rescaled volS.
  
  float ysizeA ( ysizeC ); // NEW voxel size
  float xsizeA ( xsizeC );
  float zsizeA ( zsizeC );
  
  float y_lA ( hMatGetVal(&(volS->info.M_vxl2lph), 0, 0)/ysizeS ); // OLD direction cosines
  float y_pA ( hMatGetVal(&(volS->info.M_vxl2lph), 1, 0)/ysizeS );
  float y_hA ( hMatGetVal(&(volS->info.M_vxl2lph), 2, 0)/ysizeS );
  
  float x_lA ( hMatGetVal(&(volS->info.M_vxl2lph), 0, 1)/xsizeS ); // OLD
  float x_pA ( hMatGetVal(&(volS->info.M_vxl2lph), 1, 1)/xsizeS );
  float x_hA ( hMatGetVal(&(volS->info.M_vxl2lph), 2, 1)/xsizeS );
  
  float z_lA ( hMatGetVal(&(volS->info.M_vxl2lph), 0, 2)/zsizeS ); // OLD
  float z_pA ( hMatGetVal(&(volS->info.M_vxl2lph), 1, 2)/zsizeS );
  float z_hA ( hMatGetVal(&(volS->info.M_vxl2lph), 2, 2)/zsizeS );
  
  // Build new vox --> lph matrix, Mvxl2lphA, for volA using the
  // NEW voxel sizes and the OLD direction cosines.
  
  hmat Mvxl2lphA;
  hMatEye(&Mvxl2lphA);
  
  hmat Mpat2grad;
  hMatEye(&Mpat2grad);
  
  hvec rvec, cvec, dvec;
  
  hVecInit(&rvec, ysizeA*y_lA, ysizeA*y_pA, ysizeA*y_hA, 0);
  hMatSetCol(&Mvxl2lphA, &rvec, 0);
  
  hVecInit(&cvec, xsizeA*x_lA, xsizeA*x_pA, xsizeA*x_hA, 0);
  hMatSetCol(&Mvxl2lphA, &cvec, 1);
  
  hVecInit(&dvec, zsizeA*z_lA, zsizeA*z_pA, zsizeA*z_hA, 0);
  hMatSetCol(&Mvxl2lphA, &dvec, 2);
  
  int heightA  ( heightC );  // NEW
  int widthA   ( widthC  );
  int depthA   ( depthC  );
  
  int dimA[4] = {heightA, widthA, depthA, frames};
  
  hvec cent_vxlS; // OLD
  hVecInit(&cent_vxlS, (heightS-1.0)/2.0, (widthS-1.0)/2.0, (depthS-1.0)/2.0, 1.0);
  
  hvec cent_lph;  // OLD; volA will have same c_l, c_p, c_h (or c_r, c_a, c_s) as volS.
  hMatMultiplyhVec(&(volS->info.M_vxl2lph), &cent_vxlS, &cent_lph);
  
  hVecSetVal(&cent_lph, 3, 0.0); // Change 4th component of cent_lph from 1.0 to 0.0
  
  hvec cent_vxlANeg;
  hVecInit(&cent_vxlANeg, -(heightA-1.0)/2.0, -(widthA-1.0)/2.0, -(depthA-1.0)/2.0, 1.0);
  
  hvec T1, T;
  hMatMultiplyhVec(&Mvxl2lphA, &cent_vxlANeg, &T1);
  hVecAdd(&T1, &cent_lph, &T);
  hMatSetCol(&Mvxl2lphA, &T, 3);
  
  fvol* volA( fVolNew(dimA, &Mvxl2lphA, &Mpat2grad, NULL) ); // The new image as a tabla raza.
  
  // Now, load up volA with volS rescaled to fit the sepcified dimensions.
  
  hvec lph, voxS, voxA;
  hmat lph2VoxS;
  hMatInverse(&(volS->info.M_vxl2lph), &lph2VoxS); // volS->info.M_vxl2lph is vox --> lph matrix for volS.
  
  const INTERPMETHOD interpMethod ( useCubic ? CUBIC : LINEAR ); // NEAREST  LINEAR   CUBIC
  const INTERPEDGE   interpEdge   ( THROWAWAY );                 // MIRROR   ZEROPAD  THROWAWAY
  int bPositive ( 0 );
  
  for(int i=0; i<heightA; ++i)
    for(int j=0; j<widthA; ++j)
      for(int k=0; k<depthA; ++k) {
	
	hVecInit(&voxA, i, j, k, 1.0);
	hMatMultiplyhVec(&Mvxl2lphA, &voxA, &lph);
	hMatMultiplyhVec(&lph2VoxS, &lph, &voxS);
	
	float valS;
	
	GetVxlVal(&voxS, &valS, volS, interpMethod, interpEdge, bPositive);
	
	fVolSetVal(volA, i, j, k, valS);
      }
  
  fVolDelete(volS);
  volS = volA;
  
  return true; // input volS rescaled.
}



/**********************************************************************************************************************
 **********************************************************************************************************************
 **********************************************************************************************************************
#define MRIgetVoxelToRasXform   extract_i_to_r
#define MRIgetRasToVoxelXform   extract_r_to_i

MATRIX *extract_r_to_i(MRI *mri)
{
  MATRIX *m_ras_to_voxel, *m_voxel_to_ras ;

  m_voxel_to_ras = extract_i_to_r(mri) ;
  m_ras_to_voxel = MatrixInverse(m_voxel_to_ras, NULL) ;
  MatrixFree(&m_voxel_to_ras) ;
  return(m_ras_to_voxel) ;
}
  

MATRIX *extract_i_to_r(MRI *mri)
{
  MATRIX *m;
  float m11, m12, m13, m14;
  float m21, m22, m23, m24;
  float m31, m32, m33, m34;
  float ci, cj, ck;

  m = MatrixAlloc(4, 4, MATRIX_REAL);
  if(m == NULL)
  {
    ErrorReturn(NULL, (ERROR_BADPARM, "extract_i_to_r(): error allocating matrix"));
  }

  m11 = mri->info.vxlsize[1] * mri->x_r;  m12 = mri->info.vxlsize[0] * mri->y_r;  m13 = mri->info.vxlsize[2] * mri->z_r;
  m21 = mri->info.vxlsize[1] * mri->x_a;  m22 = mri->info.vxlsize[0] * mri->y_a;  m23 = mri->info.vxlsize[2] * mri->z_a;
  m31 = mri->info.vxlsize[1] * mri->x_s;  m32 = mri->info.vxlsize[0] * mri->y_s;  m33 = mri->info.vxlsize[2] * mri->z_s;

  ci = (mri->info.dim[1]) / 2.0;
  cj = (mri->info.dim[0]) / 2.0;
  ck = (mri->info.dim[2]) / 2.0;
  
  m14 = mri->c_r - (m11 * ci + m12 * cj + m13 * ck);
  m24 = mri->c_a - (m21 * ci + m22 * cj + m23 * ck);
  m34 = mri->c_s - (m31 * ci + m32 * cj + m33 * ck);

  stuff_four_by_four(m, m11, m12, m13, m14, 
		        m21, m22, m23, m24, 
		        m31, m32, m33, m34, 
		        0.0, 0.0, 0.0, 1.0);

  return(m);

}
***********************************************************************************************************************
***********************************************************************************************************************
***********************************************************************************************************************/
