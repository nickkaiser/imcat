# makefile for pgplot stuff

TOOLS = contour plotcat
SCRIPTS = epsfcompose

# defaults (may be overridden by top level Makefile)
PGLIBDIR = $(IMCATDIR)/lib/$(ARCH)/pgplot/
XLIBS = -L/usr/X11R6/lib -lXt -lSM -lICE -lXext -lX11
# XLIBS = -lX11
FLIBS=/home/kaiser/research/imcat/lib/Linux/dummy_MAIN.o -lf2c -lm
# FLIBS = -lUfor -lfor -lFutil -lm -lots

# IMCATLIBS = -L../../../lib/$(ARCH) -lim -lcat -lutils
IMCATLIBS = -L$(IMCATDIR)/lib/$(ARCH) -lim -lcat $(IMCATDIR)/lib/$(ARCH)/libutils.a
PGLIBS = -L$(PGLIBDIR) -lcpgplot -lpgplot $(PGEXTRAS)
CLIBS = -lm

BINDIR = ../../../bin/$(ARCH)
SCRIPTDIR = ../../../bin/scripts
DOCDIR = ../../../doc/imcattools/pgplot

SHELL = /bin/sh

all : $(TOOLS)

install : $(TOOLS)
	@echo installing imcattools_pgplot in $(BINDIR) ....
	for i in $(SCRIPTS); do \
		if test -h $(SCRIPTDIR)/$$i ; then rm -f $(SCRIPTDIR)/$$i ; fi ; \
		ln -s ../../src/imcattools/pgplot/$$i $(SCRIPTDIR)/$$i ; \
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
	@echo making html Man Pages for imcattools_pgplot
	@if test ! -d $(DOCDIR) ; then mkdir -p $(DOCDIR) ; fi
	makemanpages imcattools_pgplot  $(MANDIR) $(DOCDIR) $(TOOLS) $(SCRIPTS)

##
## the tools:
##


contour : contour.o dummy_MAIN.o
	$(FC) -o contour contour.o $(PGLIBS) $(XLIBS) $(IMCATLIBS) $(FLIBS) $(CLIBS)

contour.o : contour.c
	$(CC) -c contour.c $(CFLAGS) -I$(PGLIBDIR) -o contour.o

plotcat : plotcatmain.o getrange.o cpgellipses.o rectangles.o lcfilter.o dummy_MAIN.o
	$(FC) -o plotcat plotcatmain.o getrange.o cpgellipses.o rectangles.o \
		lcfilter.o $(IMCATLIBS) $(PGLIBS) $(XLIBS) $(FLIBS) $(CLIBS)

plotcatmain.o : plotcatmain.c
	$(CC) -c plotcatmain.c $(CFLAGS) -I$(PGLIBDIR) -o plotcatmain.o 

getrange.o : getrange.c
	$(CC) -c getrange.c -o getrange.o

cpgellipses.o : cpgellipses.c
	$(CC) -c cpgellipses.c $(CFLAGS) -I$(PGLIBDIR) -o cpgellipses.o

rectangles.o : rectangles.c
	$(CC) -c rectangles.c -o rectangles.o

lcfilter.o : lcfilter.c
	$(CC) -c lcfilter.c -o lcfilter.o


