#ifndef GUARD_loadDispy
#define GUARD_loadDispy

#include <string>

void
loadDispy(const std::string& dyFile,
	  const int& numVox,
	  const int& voxStep,
	  const int nvoxNewZbdry,
	  const InputParams& p,
	  double* dispy);

#endif
