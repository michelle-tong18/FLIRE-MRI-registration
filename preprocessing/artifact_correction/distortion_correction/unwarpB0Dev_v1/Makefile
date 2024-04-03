include MakeMacros
include MakeImplicit
include ObsList

#EXE = unwarpB0Dev
EXE = epic # EPI Correction

##MKLDIR  = /opt/intel/mkl/10.0.2.018
##MKLDIR  = /usr/pubsw/packages/opt/intel/mkl/10.0.2.018
##MKLLIBS = -L${MKLDIR}/lib/em64t -lmkl -lmkl_lapack -lguide -lpthread
##MKLLIBS = -L/opt/intel/cmkl/10.0.2.018/lib/em64t -lmkl -lmkl_lapack -lguide -lpthread
#MKLLIBS = -lmkl -lmkl_lapack -lguide -lpthread
MKLLIBS = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

# ip113 ip114
GSLLIBS = # NOT NEEDED! /usr/pubsw/packages/lib/libgsl.a # /usr/local/lib/libgsl.a # -lgsl #-lgslcblas
BLASLIB = -L/usr/lib/gcc/x86_64-redhat-linux/3.4.6 -lg2c #-lblas -lg2c

# ip18 ip15 clusters
#GSLLIBS = -lgsl # /usr/pubsw/packages/lib/libgsl.a # /usr/local/lib/libgsl.a # -lgsl #-lgslcblas
#BLASLIB = -lg2c # -L/usr/lib/gcc/x86_64-redhat-linux/3.4.6 -lg2c #-lblas -lg2c

# fftw3.h in /usr/local/include
# whereis libfftw3f.a
# libfftw3f: /usr/local/lib/libfftw3f.la /usr/local/lib/libfftw3f.a
INCS = -I/usr/local/include $(CTXINCBS) $(CTXINCINTERP)
LIBS = -lfftw3f $(CTXLIBBS) $(CTXLIBINTERP) $(MKLLIBS) $(BLASLIB) $(GSLLIBS) -lm
include MakeCTX

include $(CTXDIR)/oscpu.def  # See ~/mri/dev/unwarpB0Dev/MakeCTX_SAVE
#include /home/holland/mri/ctxsrcDev/oscpu.def

$(EXE): $(OBS) *.h
	$(LDRC) $(INCS) $(LDFLAGS) $(OBS) -o $@ $(LIBS)

install: $(EXE)
	cp $(EXE) /home/ccconlin/Code/Distortion_correction/bin/

clean:
	rm -f *.o *~ $(EXE)
