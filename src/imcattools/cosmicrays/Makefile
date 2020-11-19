# Makefile for cosmic ray removal tools

TOOLS = removecrs  gethotpix
SCRIPTS = 

LIBS = -L$(IMCATDIR)/lib/$(ARCH) -lim -lcat -lutils -lm
INCLUDES = -I$(IMCATDIR)/src

BINDIR = $(IMCATDIR)/bin/$(ARCH)
SCRIPTDIR = $(IMCATDIR)/bin/scripts
DOCDIR = $(IMCATDIR)/doc/imcattools/cosmicrays

SHELL = /bin/sh

all : $(TOOLS)

install : $(TOOLS)
	@echo installing cosmic ray removal tools in $(BINDIR) ....
#	for i in $(SCRIPTS); do \
#		if test -h $(SCRIPTDIR)/$$i ; then rm -f $(SCRIPTDIR)/$$i ; fi ; \
#		ln -s ../../src/imcattools/cosmicrays/$$i $(SCRIPTDIR)/$$i ; \
#	done
#	@ if test -d $(BINDIR) ; then \
#		echo $(BINDIR) exists ; else mkdir $(BINDIR) ; fi
	for i in $(TOOLS); do \
		cp $$i $(BINDIR)/$$i ; \
	done

clean:
	rm -f core *.o
	for i in $(TOOLS); do rm -f $$i ; done

docs :
	@echo making html Man Pages for imcattools_cosmicrays
	@if test ! -d $(DOCDIR) ; then mkdir -p $(DOCDIR) ; fi
	makemanpages imcattools_cosmicrays $(MANDIR) $(DOCDIR) $(TOOLS) $(SCRIPTS)

##
## the tools:
##

removecrs : removecrs.o
	$(CC) removecrs.o $(LIBS) -o removecrs

gethotpix : gethotpix.o
	$(CC) gethotpix.o $(LIBS) -o gethotpix





##
## the objects
##


removecrs.o :	removecrs.c
	$(CC) -c removecrs.c $(INCLUDES) -o removecrs.o

gethotpix.o :	gethotpix.c
	$(CC) -c gethotpix.c $(INCLUDES) -o gethotpix.o

