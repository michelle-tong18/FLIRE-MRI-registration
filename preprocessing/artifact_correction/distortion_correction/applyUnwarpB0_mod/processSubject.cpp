// ==================================================
// Copyright (c) 2010 Dominic Holland and Anders Dale
// All rights reserved
// ==================================================

#include <ctime>             // for clock()
#include <string>
#include "inputParams.h"
#include "applyDeformationField.h"
#include "imgvol.h"          // CTX stuff. fvol
#include "loadFvols.h"
#include "loadFvol.h"
#include "writeImage.h"

using std::string;


void
processSubject(const InputParams& p) {
  
  fvols* forwardVols ( loadFvols( p.getForwardImageInFileName() ) );
  fvols* reverseVols ( loadFvols( p.getReverseImageInFileName() ) );
  fvol*  dVol        ( loadFvol ( p.getDisplacementFieldFileName() ) );
  
  string outDir               ( p.getOutDir() + "/" );
  string forwardCorrectedFile ( outDir + p.getForwardImageOutFileName() );
  string reverseCorrectedFile ( outDir + p.getReverseImageOutFileName() );
  
  if(forwardVols) {
    applyDeformationField(forwardVols, dVol, 1, p.getJic());
    writeImage(forwardCorrectedFile, forwardVols);
    fVolsDelete(forwardVols);
  }
  
  if(reverseVols) {
    applyDeformationField(reverseVols, dVol, -1, p.getJic());
    writeImage(reverseCorrectedFile, reverseVols);
    fVolsDelete(reverseVols);
  }
  
  fVolDelete(dVol);
  
}
