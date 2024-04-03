#ifndef GUARD_applyDeformationField
#define GUARD_applyDeformationField

#include "interpolation.h"      // CTX stuff. fvol

void
applyDeformationField(fvols* vols, fvol* dVol, const int direction, const bool jic);

#endif
