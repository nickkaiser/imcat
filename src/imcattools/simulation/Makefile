# Makefile for simulation tools

TOOLS = makecosmocat makemockimage circimfromcat
SCRIPTS = 

LIBS = -L$(IMCATDIR)/lib/$(ARCH) -lim -lcat -lutils -lm
INCLUDES = -I$(IMCATDIR)/src

BINDIR = $(IMCATDIR)/bin/$(ARCH)
SCRIPTDIR = $(IMCATDIR)/bin/scripts
DOC = $(IMCATDIR)/doc/simulation.txt
DOCDIR = $(IMCATDIR)/doc/imcattools/simulation


SHELL = /bin/sh

all : $(TOOLS)

install : $(TOOLS)
	@echo installing simulation tools in $(BINDIR) ....
#	for i in $(SCRIPTS); do \
#		if test -h $(SCRIPTDIR)/$$i ; then rm -f $(SCRIPTDIR)/$$i ; fi ; \
#		ln -s ../../src/imcattools/simulation/$$i $(SCRIPTDIR)/$$i ; \
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
	@echo making html Man Pages for imcattools_simulation
	@if test ! -d $(DOCDIR) ; then mkdir -p $(DOCDIR) ; fi
	makemanpages imcattools_simulation $(MANDIR) $(DOCDIR) $(TOOLS) $(SCRIPTS)

##
## the tools:
##

makecosmocat : makecosmocat.o
	$(CC) makecosmocat.o $(LIBS) -o makecosmocat

makemockimage : makemockimage.o
	$(CC) makemockimage.o $(LIBS) -o makemockimage

circimfromcat : circimfromcat.o
	$(CC) circimfromcat.o $(LIBS) -o circimfromcat




##
## the objects
##


makecosmocat.o : makecosmocat.c
	$(CC) -c makecosmocat.c $(INCLUDES) -o makecosmocat.o

makemockimage.o : makemockimage.c
	$(CC) -c makemockimage.c $(INCLUDES) -o makemockimage.o
	
circimfromcat.o : circimfromcat.c
	$(CC) -c circimfromcat.c $(INCLUDES) -o circimfromcat.o

	
