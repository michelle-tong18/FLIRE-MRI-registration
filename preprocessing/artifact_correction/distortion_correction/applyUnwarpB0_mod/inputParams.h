// ==================================================
// Copyright (c) 2010 Dominic Holland and Anders Dale
// All rights reserved
// ==================================================

#ifndef GUARD_inputParams
#define GUARD_inputParams

#include <iostream>
#include <string>
#include <vector>


class InputParams {
  
 public:
  InputParams():
    forwardImageOutFileName ("fB0uw.mgz"),
    reverseImageOutFileName ("rB0uw.mgz"),
    displacementFieldFileName ("d.mgz"),
    jic(false), // applyJacobianIntensityCorrection
    outDir ("./") {};
    
  InputParams(int argc, char** argv):
    forwardImageOutFileName ("fB0uw.mgz"),
    reverseImageOutFileName ("rB0uw.mgz"),
    displacementFieldFileName ("d.mgz"),
    jic(false), // applyJacobianIntensityCorrection
    outDir ("./") { readParameters(argc, argv); }
  
  void readParameters(int argc, char**);
  std::ostream& print(std::ostream&) const;
  
  const std::string& getForwardImageInFileName()       const { return forwardImageInFileName; }
  const std::string& getReverseImageInFileName()       const { return reverseImageInFileName; }
  const std::string& getDisplacementFieldFileName()    const { return displacementFieldFileName; }
  
  const std::string& getOutDir()                       const { return outDir; }
  const std::string& getForwardImageOutFileName()      const { return forwardImageOutFileName; }
  const std::string& getReverseImageOutFileName()      const { return reverseImageOutFileName; }
  
  const bool&        getJic()                          const { return jic; } // applyJacobianIntensityCorrection
  
 private:
  std::string      forwardImageInFileName;
  std::string      reverseImageInFileName;
  std::string      displacementFieldFileName;
  
  std::string      outDir;
  std::string      forwardImageOutFileName;
  std::string      reverseImageOutFileName;
  
  bool jic;
  
};

// Need this declaration to handle expressions like "cout << p;" in main().
std::ostream& operator<<(std::ostream& os, const InputParams& p);

void errorReadFile(const char*, const char*);
void errorReadFile(const char*, const std::string);
void errorReadFile(std::ifstream&, const char*, int) throw(std::string);
void errorReadFile(const std::string, const std::string);

#endif
