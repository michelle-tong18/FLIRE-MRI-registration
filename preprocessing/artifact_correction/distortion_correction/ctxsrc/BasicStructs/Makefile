# Makefile for ctxlab
#
# VF7K-R4677DX3
SHELL = /bin/sh

#
# Compilation Flags (-g == debug; change to optimize for release)
#
LIBCXX = icc -m32 #g++ #gcc g++ icc
EXECXX = icc -m32 #g++ #g++ gcc

LIBSRCS = hvecmat.c imgvol.c surface.c dmat.c
LIBLINK = corstruct
EXESRCS = test.c
EXEOUTS = test
include ../oscpu.def
LIBOUT = libcorstruct.$(LIBEXT)


test:lib	
	$(EXECXX) $(EXECFLAGS)  $(EXELINKFLAGS) $(EXESRCS) -o $(EXEOUTS)
lib:
	$(LIBCXX) $(LIBCFLAGS) $(LIBSRCS) 
	$(LIBCXX) -o $(LIBOUT) *.o $(LIBLINKFLAGS)

install:
	cp $(LIBOUT) ../../lib/
	cp *.h ../../include

clean:
	rm -f *.o $(LIBOUT) ../../lib/$(LIBOUT) $(EXEOUTS)
