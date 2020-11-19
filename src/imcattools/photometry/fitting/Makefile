# makefile for photometry_fitting tools

TOOLS = fit1object gaussfit betafit
SCRIPTS = testscript

#libraries
LIBS = -L$(IMCATDIR)/lib/$(ARCH) -lim -lfft $(FFTLIB) -lcat -lutils -lm


BINDIR = $(IMCATDIR)/bin/$(ARCH)
SCRIPTDIR = $(IMCATDIR)/bin/scripts
DOCDIR = $(IMCATDIR)/doc/imcattools/photometry/fitting

OBJS = 
SHELL = /bin/sh

all : $(TOOLS)

install : $(TOOLS)
	@echo installing imcattools/photometry/fitting in $(BINDIR) ....
	for i in $(SCRIPTS); do \
		if test -h $(SCRIPTDIR)/$$i ; then rm -f $(SCRIPTDIR)/$$i ; fi ; \
		ln -s $(IMCATDIR)/src/imcattools/photometry/fitting/$$i $(SCRIPTDIR)/$$i ; \
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
	@echo making html Man Pages for imcattools_photometry_fitting
	@if test ! -d $(DOCDIR) ; then mkdir -p $(DOCDIR) ; fi
	makemanpages imcattools_photometry_fitting $(MANDIR) $(DOCDIR) $(TOOLS) $(SCRIPTS)

##
## the tools:
##


fit1object : fit1object.o gengaussfit.o
	$(CC) fit1object.o gengaussfit.o $(LIBS) -o fit1object

gaussfit : gaussfitmain.o gaussfitn.o
	$(CC) gaussfitmain.o gaussfitn.o $(LIBS) -o gaussfit

betafit : betafitmain.o betafit.o
	$(CC) betafitmain.o betafit.o $(LIBS) -o betafit













