# makefile for fft library

# be sure to define env variables:
#	CC = gcc or another ansi compiler
#	ARCH = whatever is appropriate

# choose an FFT package, CFLAGS - may be  overridden by master Makefile
# FFT=NR
# FFTLIB=
# FFT=SGI
# FFTLIB=-L../sgi -lfft
# FFT=FFTPACK
# FFTLIB=-L../fftpack -lfftpack_c
FFT=FFTW
FFTLIB=

# CFLAGS=-DFFTPACK_FFT
CFLAGS=-DFFTW_FFT


OBJECTS = fft_$(FFT).o myfft.o rlft3.o fourn.o shmalloc.o
INCLUDES = myfft.h ../utils/error.h 
HEADERS = myfft.h fft_NR.h fft_SGI.h fft_FFTPACK.h
DOCDIR = ../../doc/fftlib

SHELL = /bin/sh

all : libfft.a

libfft.a : $(OBJECTS)
	ar r libfft.a $?

ranlib :
	ranlib libfft.a

install : libfft.a $(RANLIB)
	@echo +++++++++++++ installing libfft.a in ../../lib/$(ARCH)/ +++++++++++++
	@ if test -d ../../lib/$(ARCH) ; then \
		echo  +++++++++++++ dest-dir exists +++++++++++++; else mkdir ../../lib/$(ARCH) ; fi
	cp libfft.a ../../lib/$(ARCH)

docs :
	@echo making html Man Pages for libfft
	@if test ! -d $(DOCDIR) ; then mkdir $(DOCDIR) ; fi
	@cp docs_head.html docs.html
	@for h in $(HEADERS) ; do \
		echo '<h2><a name="'$$h'">'$$h'</a></h2><pre>' >> docs.html ; \
		cat $$h >> docs.html ; \
		echo '</pre><p>back to <a href="#INDEX">index</a><hr>' >> docs.html ; \
	done
	@cat docs_tail.html >> docs.html
	@cp docs.html $(DOCDIR)/docs.html
	makecatmanpage fftlib $(HEADERS) > $(MANDIR)/cat3/fftlib.3

clean:
	rm -f *.o libfft.a

myfft.o : myfft.c $(INCLUDES)

fft_$(FFT).o : fft_$(FFT).c $(INCLUDES)

