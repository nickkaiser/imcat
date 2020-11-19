# Makefile to make registration tools

TOOLS = acfregister lintransfit fitstack warpcat fit2cats fitdistortion
SCRIPTS = registercats scalerottrans mergestack \
	registercats2 getstacktrans warpstack merge2cats


SHELL = /bin/sh
BINDIR = ../../../bin/$(ARCH)
SCRIPTDIR = ../../../bin/scripts
FFTLIB = -L$(IMCATDIR)/lib/$(ARCH) -lfftpack_c
LIBS = -L../../../lib/$(ARCH) -lcat -lim -lfft $(FFTLIB) -lutils $(FLIBS) -lm
DOCDIR = ../../../doc/imcattools/registration
CFLAGS = -I$(IMCATDIR)/src -DFFTPACK_FFT

all : $(TOOLS)


clean :
	rm -f *.o core
	for i in $(TOOLS); do rm -f $$i ; done

install : $(TOOLS)
	@echo installing $(SCRIPTS) in $(SCRIPTDIR)
	for i in $(SCRIPTS) ; do \
		if test -h $(SCRIPTDIR)/$$i ; then rm $(SCRIPTDIR)/$$i ; fi ;\
		ln -s ../../src/imcattools/registration/$$i $(SCRIPTDIR)/$$i ;\
	done
	@echo installing $(TOOLS) in $(BINDIR) ....
	@ if test ! -d $(BINDIR) ; then  mkdir $(BINDIR) ; fi
	for i in $(TOOLS) ; do cp $$i $(BINDIR)/$$i ; done

docs :
	@echo making html Man Pages for imcattools_registration
	@if test ! -d $(DOCDIR) ; then mkdir -p $(DOCDIR) ; fi
	makemanpages imcattools_registration $(MANDIR) $(DOCDIR) $(TOOLS) $(SCRIPTS)

acfregister : acfregister.o
	$(CC) acfregister.o $(LIBS) -o acfregister

lintransfit : lintransfit.o
	$(CC) lintransfit.o $(LIBS) -o lintransfit

fitstack : fitstack.o fitstack_read.o
	$(CC) fitstack.o fitstack_read.o $(LIBS) -o fitstack

fit2cats : fit2cats.o
	$(CC) fit2cats.o $(LIBS) -o fit2cats

warpcat : warpcat.o
	$(CC) warpcat.o $(LIBS) -o warpcat

fitdistortion : fitdistortion.o
	$(CC) fitdistortion.o $(LIBS) -o fitdistortion
