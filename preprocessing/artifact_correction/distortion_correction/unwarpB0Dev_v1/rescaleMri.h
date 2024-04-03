#ifndef GUARD_rescaleMri
#define GUARD_rescaleMri

#include "imgvol.h" // CTX stuff. fvol

bool
rescaleMri(fvol* atlasMri, fvol*& subjectMri, const bool useCubic = false);

bool
rescaleMri(const int   awidth, const int   aheight, const int   adepth,
	   const float axsize, const float aysize,  const float azsize, fvol*& subjectMri, const bool useCubic = false);

#endif
