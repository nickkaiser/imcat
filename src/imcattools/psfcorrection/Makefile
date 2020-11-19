# makefile for psf-correction tools

TOOLS = efit ecorrect
SCRIPTS = ecorrect2

#libraries
LIBS = -L../../../lib/$(ARCH) -lcat -lutils -lm


BINDIR = ../../../bin/$(ARCH)
SCRIPTDIR = ../../../bin/scripts
DOCDIR = ../../../doc/imcattools/psfcorrection

SHELL = /bin/sh

all : $(TOOLS)

install : $(TOOLS)
	@echo installing psf-correction tools in $(BINDIR) ....
	for i in $(SCRIPTS); do \
		if test -h $(SCRIPTDIR)/$$i ; then rm -f $(SCRIPTDIR)/$$i ; fi ; \
		ln -s ../../src/imcattools/psfcorrection/$$i $(SCRIPTDIR)/$$i ; \
	done
	@ if test -d $(BINDIR) ; then \
		echo $(BINDIR) exists ; else mkdir $(BINDIR) ; fi
	for i in $(TOOLS); do \
		mv $$i $(BINDIR)/$$i ; \
	done

clean:
	rm -f core *.o
	for i in $(TOOLS); do rm -f $$i ; done

docs :
	@echo making html Man Pages for imcattools_psfcorrection
	@if test ! -d $(DOCDIR) ; then mkdir -p $(DOCDIR) ; fi
	makemanpages imcattools_psfcorrection $(MANDIR) $(DOCDIR) $(TOOLS) $(SCRIPTS)

##
## the tools:
##

efit : efit.o efit_stuff.o
	$(CC) efit.o efit_stuff.o $(LIBS) -o efit


ecorrect : ecorrect.o efit_stuff.o
	$(CC) ecorrect.o efit_stuff.o $(LIBS) -o ecorrect


















