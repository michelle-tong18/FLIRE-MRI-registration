#include <iostream>


void
usage() {
  
  std::cout <<
    "\n\nUsage:\n\n"
    "This program applies a deformation field to unwarp B0-distorted images, e.g., EPI and MPRAGE.\n"
    "A forward- or reverse-distorted image must be supplied (or both together), specified with -f or -r, respectively.\n"
    "The displacement field must also be supplied, specified with -d.\n"
    "These may be supplied as part of an optional input-parameters file,\n"
    "possibly called \"inputParams.txt\", specified with -ip.\n"
    "Options may be supplied in any order.\n"
    "Command line values supersede any values supplied in the optional parameters file.\n"
    "Thus, for example,\n\n"
    "./applyUnwarpB0 -ip inputParams.txt\n"
    "or use defaults\n"
    "./applyUnwarpB0 -f forwardImageIn.mgz -d d.mgz\n"
    "or\n"
    "./applyUnwarpB0 -r reverseImageIn.mgz -d d.mgz\n\n"
    "Full list of options:\n"
    "[-ip] <inputParamsFile.txt> \n"
    "[-f]  <forwardImageIn.mgz> \n"
    "[-r]  <reverseImageIn.mgz> \n"
    "[-fo] <forwardImageOut.mgz> \n"
    "[-ro] <reverseImageOut.mgz> \n"
    "[-d]  <displacementFieldFileName.mgz> \n"
    "[-od] <outDir> (A '/' must precede the text name of the directory.)\n"
    "[-jic]         To apply Jacobian Intenisty Correction\n"
	    << std::endl;
}
