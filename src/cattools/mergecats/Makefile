# Makefile to make mergecats

TOOLS = mergecats unmergecat getclosepairs getisolatedobjects cleancat \
	findtrailobjects nuketrailobjects
OBJECTS = mergecatsmain.o grid.o getclosepairs.o cleancat.o


SHELL = /bin/sh
BINDIR = ../../../bin/$(ARCH)
LIBS = -L../../../lib/$(ARCH) -lcat -lim -lutils -lm
DOCDIR = ../../../doc/cattools/mergecats

all : $(TOOLS)


clean :
	rm -f *.o core
	for i in $(TOOLS); do rm -f $$i ; done

install : $(TOOLS)
	@echo installing $(TOOLS) in $(BINDIR) ....
	@ if test ! -d $(BINDIR) ; then mkdir $(BINDIR) ; fi
	for i in $(TOOLS); do cp $$i $(BINDIR)/$$i ; done

docs :
	@echo making html Man Pages for cattools_mergecats
	@if test ! -d $(DOCDIR) ; then mkdir $(DOCDIR) ; fi
	makemanpages cattools_mergecats $(MANDIR) $(DOCDIR) $(TOOLS) $(SCRIPTS)

mergecats : mergecatsmain.o grid.o
	$(CC) mergecatsmain.o grid.o $(LIBS) -o mergecats

unmergecat : unmergecat.o
	$(CC) unmergecat.o $(LIBS) -o unmergecat

getclosepairs : getclosepairs.o grid.o
	$(CC) getclosepairs.o grid.o $(LIBS) -o getclosepairs

getisolatedobjects : getisolatedobjects.o grid.o
	$(CC) getisolatedobjects.o grid.o $(LIBS) -o getisolatedobjects

cleancat : cleancat.o grid.o
	$(CC) cleancat.o grid.o $(LIBS) -o cleancat

findtrailobjects : findtrailobjects.o grid.o
	$(CC) findtrailobjects.o grid.o $(LIBS) -o findtrailobjects

nuketrailobjects : nuketrailobjects.o
	$(CC) nuketrailobjects.o $(LIBS) -o nuketrailobjects
