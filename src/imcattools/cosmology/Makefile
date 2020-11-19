# Makefile for cosmology tools

TOOLS = frw nfw
SCRIPTS = 

LIBS = -L$(IMCATDIR)/lib/$(ARCH) -lim -lutils -lm
INCLUDES = -I$(IMCATDIR)/src

BINDIR = $(IMCATDIR)/bin/$(ARCH)
SCRIPTDIR = $(IMCATDIR)/bin/scripts
DOCDIR = $(IMCATDIR)/doc/imcattools/cosmology

SHELL = /bin/sh

all : $(TOOLS)

install : $(TOOLS)
	@echo installing cosmology tools in $(BINDIR) ....
#	for i in $(SCRIPTS); do \
#		if test -h $(SCRIPTDIR)/$$i ; then rm -f $(SCRIPTDIR)/$$i ; fi ; \
#		ln -s ../../src/imcattools/cosmology/$$i $(SCRIPTDIR)/$$i ; \
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
	@echo making Man Pages for imcattools_cosmology
	@if test ! -d $(DOCDIR) ; then mkdir -p $(DOCDIR) ; fi
	makemanpages imcattools_cosmology $(MANDIR) $(DOCDIR) $(TOOLS) $(SCRIPTS)


##
## the tools:
##

frw : frw.o
	$(CC) frw.o -lm -o frw

frw.o : frw.c

nfw : nfw.o
	$(CC) nfw.o $(LIBS) -lm -o nfw

nfw.o : nfw.c
	$(CC) -c nfw.c $(INCLUDES) -o nfw.o
