# makefile for object finding tools

TOOLS = mosaicfit fitgeometry fitmagshifts fitgeometry2 fitgeometry3 fitpolymodelmagshift fitmagshiftsplus
SCRIPTS = makemosaiccats mergemosaiccats1 transformmosaiccat \
	mergemosaiccats2 makemosaicstack makemosaicstacks makewarpedalbum printimages \
	combinestacks fitshifts getfieldgeomfromheaders makestack

# retired  mosaicmap mosaicmap2 


#libraries
LIBS = -L../../../lib/$(ARCH) -lim -lfft -lcat -lutils -lm


BINDIR = ../../../bin/$(ARCH)
SCRIPTDIR = ../../../bin/scripts
DOCDIR = ../../../doc/imcattools/mosaics

SHELL = /bin/sh

all : $(TOOLS)

install : $(TOOLS)
	@echo installing mosaic tools in $(BINDIR) ....
	for i in $(SCRIPTS); do \
		if test -h $(SCRIPTDIR)/$$i ; then rm -f $(SCRIPTDIR)/$$i ; fi ; \
		ln -s ../../src/imcattools/mosaics/$$i $(SCRIPTDIR)/$$i ; \
	done
	@ if test -d $(BINDIR) ; then \
		echo $(BINDIR) exists ; else mkdir $(BINDIR) ; fi
	for i in $(TOOLS); do \
		cp $$i $(BINDIR)/$$i ; \
	done

clean:
	rm -f core *.o
	for i in $(TOOLS); do rm -f $$i ; done

docs :
	@echo making html Man Pages for imcattools_mosaics
	@if test ! -d $(DOCDIR) ; then mkdir -p $(DOCDIR) ; fi
	makemanpages imcattools_mosaics $(MANDIR) $(DOCDIR) $(TOOLS) $(SCRIPTS)

##
## the tools:
##

mosaicfit : mosaicfit.o
	$(CC) mosaicfit.o $(LIBS) -o mosaicfit

mosaicmap : mosaicmap.o
	$(CC) mosaicmap.o $(LIBS) -o mosaicmap

mosaicmap2 : mosaicmap2.o
	$(CC) mosaicmap2.o $(LIBS) -o mosaicmap2

fitgeometry : fitgeometry.o
	$(CC) fitgeometry.o $(LIBS) -o fitgeometry

fitmagshifts : fitmagshifts.o
	$(CC) fitmagshifts.o $(LIBS) -o fitmagshifts

fitmagshiftsplus : fitmagshiftsplus.o
	$(CC) fitmagshiftsplus.o $(LIBS) -o fitmagshiftsplus

fitpolymodelmagshift : fitpolymodelmagshift.o
	$(CC) fitpolymodelmagshift.o $(LIBS) -o fitpolymodelmagshift

fitgeometry2 : fitgeometry2.o
	$(CC) fitgeometry2.o $(LIBS) -o fitgeometry2

fitgeometry2.o : fitgeometry2.c
	$(CC) -c fitgeometry2.c -I$(IMCATDIR)/src -o fitgeometry2.o

fitgeometry3 : fitgeometry3.o
	$(CC) fitgeometry3.o $(LIBS) -o fitgeometry3

fitgeometry3.o : fitgeometry3.c
	$(CC) -c fitgeometry3.c -I$(IMCATDIR)/src -o fitgeometry3.o















