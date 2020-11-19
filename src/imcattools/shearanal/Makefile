# makefile for shear analysis tools

TOOLS = etprofile massmap massmap_mp massmap_ft makeshearimage maketestshear \
	sheartogradkappa apmass
SCRIPTS =

#libraries
LIBS = -L../../../lib/$(ARCH) -lim -lfft $(FFTLIB) -lutils $(FLIBS) -lm


BINDIR = ../../../bin/$(ARCH)
SCRIPTDIR = ../../../bin/scripts
DOCDIR = ../../../doc/imcattools/shearanal

SHELL = /bin/sh

all : $(TOOLS)

install : $(TOOLS)
	@echo installing shear analysis tools in $(BINDIR) ....
#	for i in $(SCRIPTS); do \
#		if test -h $(SCRIPTDIR)/$$i ; then rm -f $(SCRIPTDIR)/$$i ; fi ; \
#		ln -s ../../src/imcattools/psfcorrection/$$i $(SCRIPTDIR)/$$i ; \
#	done
	@ if test -d $(BINDIR) ; then \
		echo $(BINDIR) exists ; else mkdir $(BINDIR) ; fi
	for i in $(TOOLS); do \
		cp $$i $(BINDIR)/$$i ; \
	done

clean:
	rm -f core *.o
	for i in $(TOOLS); do rm -f $$i ; done

docs :
	@echo making html Man Pages for imcattools_shearanal
	@if test ! -d $(DOCDIR) ; then mkdir -p $(DOCDIR) ; fi
	makemanpages imcattools_shearanal $(MANDIR) $(DOCDIR) $(TOOLS) $(SCRIPTS)

##
## the tools:
##

etprofile : etprofile.o
	$(CC) etprofile.o $(LIBS) -o etprofile


massmap : massmap.o
	$(CC) massmap.o $(LIBS) -o massmap

massmap_mp : massmap_mp.o
	$(CC) massmap_mp.o $(LIBS) -o massmap_mp

massmap_ft : massmap_ft.o
	$(CC) massmap_ft.o $(LIBS) -o massmap_ft

makeshearimage : makeshearimage.o
	$(CC) makeshearimage.o $(LIBS) -o makeshearimage

maketestshear : maketestshear.o
	$(CC) maketestshear.o $(LIBS) -o maketestshear

sheartogradkappa : sheartogradkappa.o
	$(CC) sheartogradkappa.o $(LIBS) -o sheartogradkappa

apmass : apmass.o kernels.o
	$(CC) apmass.o kernels.o $(LIBS) -o apmass
















