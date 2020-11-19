# Makefile for optics tools

TOOLS = phasetopsf trackpeak stackpsfs makeotf makevonkarmanS make2Dotf otftopsf
SCRIPTS = 

LIBS = -L$(IMCATDIR)/lib/$(ARCH) -lcat -lim -lfft $(FFTLIB) -lutils $(FLIBS) -lm
INCLUDES = -I$(IMCATDIR)/src

BINDIR = $(IMCATDIR)/bin/$(ARCH)
SCRIPTDIR = $(IMCATDIR)/bin/scripts
DOCDIR = $(IMCATDIR)/doc/imcattools/optics

SHELL = /bin/sh

all : $(TOOLS)

install : $(TOOLS)
	@echo installing cosmic ray removal tools in $(BINDIR) ....
#	for i in $(SCRIPTS); do \
#		if test -h $(SCRIPTDIR)/$$i ; then rm -f $(SCRIPTDIR)/$$i ; fi ; \
#		ln -s ../../src/imcattools/optics/$$i $(SCRIPTDIR)/$$i ; \
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
	@echo making Man Pages for imcattools_optics
	@if test ! -d $(DOCDIR) ; then mkdir -p $(DOCDIR) ; fi
	makemanpages imcattools_optics $(MANDIR) $(DOCDIR) $(TOOLS) $(SCRIPTS)




##
## the tools:
##

trackpeak.o : trackpeak.c
	$(CC) -c trackpeak.c $(INCLUDES) -o trackpeak.o

trackpeak : trackpeak.o poidev.o gammln.o findpeak.o
	$(CC) trackpeak.o poidev.o gammln.o findpeak.o $(LIBS) -lm -o trackpeak

phasetopsf.o : phasetopsf.c
	$(CC) -c phasetopsf.c $(INCLUDES) -DFFTPACK_FFT -o phasetopsf.o

phasetopsf : phasetopsf.o
	$(CC) phasetopsf.o $(FFTLIB) $(LIBS) -lm -o phasetopsf

stackpsfs.o : stackpsfs.c
	$(CC) -c stackpsfs.c $(INCLUDES) -o stackpsfs.o

stackpsfs: stackpsfs.o
	$(CC) stackpsfs.o $(LIBS) -lm -o stackpsfs

makeotf : makeotf.o
	$(CC) makeotf.o $(FFTLIB) $(LIBS) -lm -o makeotf

makeotf.o : makeotf.c
	$(CC) -c makeotf.c $(CFLAGS) $(INCLUDES) -o makeotf.o

makevonkarmanS : makevonkarmanS.o
	$(CC) makevonkarmanS.o $(LIBS) -lm -o makevonkarmanS

makevonkarmanS.o : makevonkarmanS.c
	$(CC) -c makevonkarmanS.c $(INCLUDES) -o makevonkarmanS.o

make2Dotf : make2Dotf.o
	$(CC) make2Dotf.o $(FFTLIB) $(LIBS) -lm -o make2Dotf

make2Dotf.o : make2Dotf.c
	$(CC) -c make2Dotf.c $(INCLUDES) -o make2Dotf.o

otftopsf : otftopsf.o
	$(CC) otftopsf.o $(LIBS) -o otftopsf

otftopsf.o : otftopsf.c
	$(CC) -c otftopsf.c $(INCLUDES) -o otftopsf.o

