#ifndef GUARD_resetToOriginalScaleLocation
#define GUARD_resetToOriginalScaleLocation

#include <vector>
#include "imgvol.h"              // CTX stuff. fvol

void
resetToOriginalScaleLocation(fvol*& subjectMri,
			     const std::vector<int>& origSize,
			     const std::vector<int>& targetShiftValues,
			     const float scaleFac,
			     const bool sincInterp = false,
			     const bool rescale = false);

#endif
