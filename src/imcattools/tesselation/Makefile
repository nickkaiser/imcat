# Makefile for tesselation tools

TOOLS = triangulatecat
SCRIPTS = 

LIBS = -L$(IMCATDIR)/lib/$(ARCH) -lcat -lim -lutils $(FLIBS) -lm
INCLUDES = -I$(IMCATDIR)/src

BINDIR = $(IMCATDIR)/bin/$(ARCH)
SCRIPTDIR = $(IMCATDIR)/bin/scripts
DOCDIR = $(IMCATDIR)/doc/imcattools/tesselation

SHELL = /bin/sh

all : $(TOOLS)

install : $(TOOLS)
	@echo installing tesselation removal tools in $(BINDIR) ....
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
	@echo making html Man Pages for imcattools_tesselation
	@if test ! -d $(DOCDIR) ; then mkdir -p $(DOCDIR) ; fi
	makemanpages imcattools_tesselation $(MANDIR) $(DOCDIR) $(TOOLS) $(SCRIPTS)

##
## the tools:
##

triangulatecat : triangulatecat.o triangle.o average.o painttriangle.o fixedges.o
	$(CC) triangulatecat.o triangle.o average.o painttriangle.o fixedges.o $(LIBS) -o triangulatecat

triangulatecat.o : triangulatecat.c
	$(CC) -c triangulatecat.c $(INCLUDES) -o triangulatecat.o

triangle.o : triangle.c
	$(CC) -c triangle.c -DTRILIBRARY -o triangle.o

fixedges.o : fixedges.c
	$(CC) -c fixedges.c $(INCLUDES) -o fixedges.o
