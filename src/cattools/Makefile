# Makefile for cattools

TOOLS=pastecats catstats pair gridavg maskcat fit2Dpolymodel gen2Dpolymodel \
	diff2Dpolymodel int2Dpolymodel rebin gen2Dpolymodelimage
SCRIPTS=catcats
all : $(TOOLS)

INCLUDE=../catlib/cat.h
SHELL=/bin/sh
BINDIR = ../../bin/$(ARCH)
SCRIPTDIR = ../../bin/scripts
LIBS = -L../../lib/$(ARCH) -lcat -lutils -lm
DOCDIR = ../../doc/cattools


# stuff for smselect: now superseded
# SUNLIBS = -lsuntool -lsunwindow -lpixrect
# SUNLIBS =
# XLIBS =  -lX11
# XLIBS = -lgl -lX11
# libraries
# SMLIBS = -L../../lib/sm -lsmplotsub -lsmdevices -lsmutils $(XLIBS) $(SUNLIBS)
         

clean :
	rm -f *.o core
	for i in $(TOOLS); do rm -f $$i ; done

install : $(TOOLS)
	@echo installing scripts in $(SCRIPTDIR) ....
	@ if test ! -d $(SCRIPTDIR) ; then mkdir $(SCRIPTDIR) ; fi
	for i in $(SCRIPTS); do \
		if test -h $(SCRIPTDIR)/$$i ; then rm -f $(SCRIPTDIR)/$$i ; fi ; \
		ln -s ../../src/cattools/$$i $(SCRIPTDIR)/$$i ; \
	done
	@echo installing cattools in $(BINDIR) ....
	@ if test ! -d $(BINDIR) ; then mkdir $(BINDIR) ; fi
	for i in $(TOOLS); do cp $$i $(BINDIR)/$$i ; done

docs :
	@echo making html Man Pages for cattools
	@if test ! -d $(DOCDIR) ; then mkdir $(DOCDIR) ; fi
	makemanpages cattools $(MANDIR) $(DOCDIR) $(TOOLS) $(SCRIPTS)

smselect : smselect.o
	$(CC) smselect.o $(LIBS) $(SMLIBS) -lm -o smselect

pastecats : pastecats.o
	$(CC) pastecats.o $(LIBS) -o pastecats

catstats : catstats.o
	$(CC) catstats.o $(LIBS) -o catstats

pair : pair.o
	$(CC) pair.o $(LIBS) -o pair

gridavg : gridavg.o
	$(CC) gridavg.o $(LIBS) -o gridavg

maskcat : maskcat.o
	$(CC) maskcat.o $(LIBS) -o maskcat

fit2Dpolymodel : fit2Dpolymodel.o
	$(CC) fit2Dpolymodel.o $(LIBS) -o fit2Dpolymodel

gen2Dpolymodel : gen2Dpolymodel.o
	$(CC) gen2Dpolymodel.o $(LIBS) -o gen2Dpolymodel

gen2Dpolymodelimage : gen2Dpolymodelimage.o
	$(CC) gen2Dpolymodelimage.o $(LIBS) -lim -o gen2Dpolymodelimage

diff2Dpolymodel : diff2Dpolymodel.o
	$(CC) diff2Dpolymodel.o $(LIBS) -o diff2Dpolymodel

int2Dpolymodel : int2Dpolymodel.o
	$(CC) int2Dpolymodel.o $(LIBS) -o int2Dpolymodel

rebin : rebin.o
	$(CC) rebin.o $(LIBS) -o rebin
