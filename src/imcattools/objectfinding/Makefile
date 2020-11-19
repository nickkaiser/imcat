# makefile for object finding tools

TOOLS = findpeaks getbadpix hfindpeaks hcat2cat makechart\
	makestamps minisky surfdens
SCRIPTS = chunkyfp

#libraries
LIBS = -L../../../lib/$(ARCH) -lim -lfft $(FFTLIB) -lcat -lutils $(FLIBS) -lm


BINDIR = ../../../bin/$(ARCH)
SCRIPTDIR = ../../../bin/scripts
DOC = ../../../doc/objectfinding.txt
DOCDIR = ../../../doc/imcattools/objectfinding

SHELL = /bin/sh

all : $(TOOLS)

install : $(TOOLS)
	@echo installing imcattools in $(BINDIR) ....
	for i in $(SCRIPTS); do \
		if test -h $(SCRIPTDIR)/$$i ; then rm -f $(SCRIPTDIR)/$$i ; fi ; \
		ln -s ../../src/imcattools/objectfinding/$$i $(SCRIPTDIR)/$$i ; \
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
	@echo making html Man Pages for imcattools_objectfinding
	@if test ! -d $(DOCDIR) ; then mkdir -p $(DOCDIR) ; fi
	makemanpages imcattools_objectfinding $(MANDIR) $(DOCDIR) $(TOOLS) $(SCRIPTS)


##
## the tools:
##


findpeaks :  findpeaks.o
	$(CC) findpeaks.o $(LIBS) -o findpeaks

hfindpeaks :  hfindpeaks.o
	$(CC) hfindpeaks.o $(LIBS) -o hfindpeaks

hcat2cat : hcat2cat.o
	$(CC) hcat2cat.o $(LIBS) -o hcat2cat

makechart : makechart.o
	$(CC) makechart.o $(LIBS) -o makechart

surfdens : surfdens.o
	$(CC) surfdens.o $(LIBS) -o surfdens

minisky : minisky.o
	$(CC) minisky.o $(LIBS) -o minisky

makestamps : makestamps.o
	$(CC) makestamps.o $(LIBS) -o makestamps

combinestamps : combinestamps.o
	$(CC) combinestamps.o $(LIBS) -o combinestamps

getbadpix : getbadpix.o
	$(CC) getbadpix.o $(LIBS) -o getbadpix

fp : fp.o
	$(CC) fp.o $(LIBS) -o fp



















