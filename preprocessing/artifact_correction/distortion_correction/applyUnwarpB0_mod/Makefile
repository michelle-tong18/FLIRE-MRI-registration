include MakeMacros
include MakeImplicit
include ObsList

EXE = applyUnwarpB0

#MKLLIBS = -L/opt/intel/cmkl/10.0.2.018/lib/em64t -lmkl -lmkl_lapack -lguide -lpthread
#MKLDIR  = /opt/intel/mkl/10.0.2.018
#MKLDIR  = /usr/pubsw/packages/opt/intel/mkl/10.0.2.018
#MKLLIBS = -L${MKLDIR}/lib/em64t -lmkl -lmkl_lapack -lguide -lpthread
MKLLIBS = -lmkl -lmkl_lapack -lguide -lpthread

# ip113 ip114
GSLLIBS = # NOT NEEDED! /usr/pubsw/packages/lib/libgsl.a # /usr/local/lib/libgsl.a # -lgsl #-lgslcblas
BLASLIB = -L/usr/lib/gcc/x86_64-redhat-linux/3.4.6 -lg2c #-lblas -lg2c

# ip18 ip15 clusters
# GSLLIBS = # NOT NEEDED! -lgsl # /usr/pubsw/packages/lib/libgsl.a # /usr/local/lib/libgsl.a # -lgsl #-lgslcblas
# BLASLIB = -lg2c # -L/usr/lib/gcc/x86_64-redhat-linux/3.4.6 -lg2c #-lblas -lg2c


INCS = -I/usr/local/include $(CTXINCBS) $(CTXINCINTERP)
LIBS = -lfftw3f $(CTXLIBBS) $(CTXLIBINTERP) $(MKLLIBS) $(BLASLIB) $(GSLLIBS) -lm
include MakeCTX

include $(CTXDIR)/oscpu.def

$(EXE): $(OBS) *.h
	$(LDRC) $(INCS) $(LDFLAGS) $(OBS) -o $@ $(LIBS)

install: $(EXE)
	cp $(EXE) /home/ccconlin/Code/Distortion_correction/bin

clean:
	rm -f *.o *~ $(EXE)
