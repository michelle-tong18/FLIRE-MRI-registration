#!/bin/csh
#set x
cd BasicStructs/
#setenv LD_LIBRARY_PATH  "$GSL_DIR/lib":"$LD_LIBRARY_PATH"
setenv LD_LIBRARY_PATH   ${PWD}:${LD_LIBRARY_PATH}
echo ${LD_LIBRARY_PATH}
setenv DYLD_LIBRARY_PATH ${PWD}:${DYLD_LIBRARY_PATH}
cd ..

cd Interpolation/
setenv LD_LIBRARY_PATH   ${PWD}:${LD_LIBRARY_PATH}
echo ${LD_LIBRARY_PATH}
setenv DYLD_LIBRARY_PATH ${PWD}:${DYLD_LIBRARY_PATH}
cd ..

#set MKLPATH=/home/holland/intel/mkl/10.0.1.014
#set MKLPATH=/opt/intel/cmkl/10.0.2.018
set MKLPATH=/usr/pubsw/packages/opt/intel/mkl/10.0.2.018
set MKLLIBDIR=${MKLPATH}/lib/em64t
setenv LD_LIBRARY_PATH ${MKLLIBDIR}:${LD_LIBRARY_PATH}
echo ${LD_LIBRARY_PATH}
setenv DYLD_LIBRARY_PATH ${MKLLIBDIR}:${DYLD_LIBRARY_PATH}
echo ${DYLD_LIBRARY_PATH}
