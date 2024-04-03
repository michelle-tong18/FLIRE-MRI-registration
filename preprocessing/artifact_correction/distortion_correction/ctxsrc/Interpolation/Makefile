# Makefile for ctxlab
#
SHELL = /bin/sh

#
# Compilation Flags (-g == debug; change to optimize for release)
#
LIBCXX = icc -m32 #g++ #gcc g++ icc
EXECXX = icc -m32 #g++ #gcc g++ gcc

INCBS = -I../BasicStructs
LIBBS = -L../BasicStructs -lcorstruct
LIBSRCS = interpolation.c
LIBLINK = corinterp
EXERESAMPLE = volresample
EXEMORPH = volmorph
EXEOUTS = $(EXERESAMPLE) $(EXEMORPH)
include ../oscpu.def
LIBOUT = libcorinterp.$(LIBEXT)


all: $(EXEOUTS)

$(EXEMORPH): lib
	$(EXECXX) $(EXECFLAGS) $(INCBS) $(LIBBS) $(EXELINKFLAGS) $(EXEMORPH).c -o $(EXEMORPH)

$(EXERESAMPLE): lib
	$(EXECXX) $(EXECFLAGS) $(INCBS) $(LIBBS) $(EXELINKFLAGS) $(EXERESAMPLE).c -o $(EXERESAMPLE)

lib:
	$(LIBCXX) $(LIBCFLAGS) $(INCBS) $(LIBSRCS) 
	$(LIBCXX) -o $(LIBOUT) *.o $(LIBLINKFLAGS) $(LIBBS)

install:
	cp $(LIBOUT) ../../lib/
	cp *.h ../../include
	cp $(EXEOUTS) ../../bin/

clean:
	rm -f *.o $(LIBOUT) ../../lib/$(LIBOUT) $(EXEOUTS)
