#include <iostream>


void
usage() {
  
  std::cout <<
    "Usage:\n\n"
    "This program unwarps B0-distorted images, e.g., EPI and MPRAGE.\n"
    "Forward- and reverse-distorted images must be supplied, specified with -f and -r.\n"
    "These may be supplied as part of an optional input-parameters file,\n"
    "possibly called \"inputParams.txt\", specified with -ip.\n"
    "Options may be supplied in any order.\n"
    "Command line values supersede any values supplied in the optional parameters file.\n"
    "Thus, for example,\n"
    "./unwarpB0 -ip inputParams.txt\n"
    "or use defaults\n"
    "./unwarpB0 -f forwardImageIn.mgz -r reverseImageIn.mgz\n"
    "or\n"
    "./unwarpB0 -ip inputParams.txt -f forwardImageIn.mgz -r reverseImageIn.mgz\n\n"
    "Full list of options:\n"
    "[-defaults]        Print defalut parameter values and exit \n"
    "[-ip]              <inputParamsFile.txt> \n"
    "[-f]               <forwardImageIn.mgz> \n"
    "[-r]               <reverseImageIn.mgz> \n"
    "[-fo]              <forwardImageOut.mgz> \n"
    "[-ro]              <reverseImageOut.mgz> \n"
    "[-di]              <displacementFieldInFileName.mgz> \n"
    "[-do]              <displacementFieldOutFileName.mgz> \n"
    "[-od]              <outDirFullPath> \n"
    "[-restart]         To restart using a previously caclulated displacement field\n"
    "[-voxStep]         <int n> Take only every nth voxel in each dimension; n >= 1\n"
    "[-nchunksZ]        <int> To divide system up into this many overlaping segments\n"
    "[-resample]        To reset the voxel size to xs, ys, zs (default is 4mm,4mm,4mm)\n"
    "[-xs]              <float> voxel x-dimension in which to perform estimation of deformation field\n"
    "[-ys]              <float> voxel y-dimension in which to perform estimation of deformation field\n"
    "[-zs]              <float> voxel z-dimension in which to perform estimation of deformation field\n"
    "[-smoothDfieldKernelWidth] <int> To smooth displacement field before final application and writing out [odd and positive; 1 ==> do nothing\n"
    "[-scaleImages]     Flag to scale images so that max value is given after -imageMax\n"
    "[-imageMax]        <float>\n"
    "[-kernelWidthMax]  <int> Pixel width of Gaussian kernel for blurring; an odd number >= 3, or 0\n"
    "[-kernelWidthStep] <even int>  An even number!\n"
    "[-nvoxNewZbdry]    <int> Pad z-boundaries with this many identical layers; 2 (default) is good\n"
    "[-lambda1]         <float> Coefficiennt of displacement norm (size) term in cost fn. > 0 \n"
    "[-lambda2]         <float> Coefficiennt of gradient in cost fn. > 0\n"
    "[-lambda2P]        <float> Coefficiennt of alternative gradient in cost fn. > 0\n"
	    << std::endl;
}
