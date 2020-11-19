# Makefile for waves tools

TOOLS = edw evolvescalar evolvemaxwell2D evolverealscalar evolvecomplexscalar
SCRIPTS = generate_dw

LIBS = -L$(IMCATDIR)/lib/$(ARCH) -lcat -lim -lfft $(FFTLIB) -lutils $(FLIBS) -lm
INCLUDES = -I$(IMCATDIR)/src

BINDIR = $(IMCATDIR)/bin/$(ARCH)
SCRIPTDIR = $(IMCATDIR)/bin/scripts
DOCDIR = $(IMCATDIR)/doc/imcattools/waves

SHELL = /bin/sh

all : $(TOOLS)

install : $(TOOLS)
	@echo installing waves tools in $(BINDIR) ....
	for i in $(SCRIPTS); do \
		if test -h $(SCRIPTDIR)/$$i ; then rm -f $(SCRIPTDIR)/$$i ; fi ; \
		ln -s ../../src/imcattools/waves/$$i $(SCRIPTDIR)/$$i ; \
	done
#	@ if test -d $(BINDIR) ; then \
#		echo $(BINDIR) exists ; else mkdir $(BINDIR) ; fi
	for i in $(TOOLS); do \
		cp $$i $(BINDIR)/$$i ; \
	done

clean:
	rm -f core *.o
	for i in $(TOOLS); do rm -f $$i ; done

docs :
	@echo making Man Pages for imcattools_waves
	@if test ! -d $(DOCDIR) ; then mkdir -p $(DOCDIR) ; fi
	makemanpages imcattools_waves $(MANDIR) $(DOCDIR) $(TOOLS) $(SCRIPTS)


##
## the tools:
##

evolverealscalar : evolverealscalar.o
	$(CC) evolverealscalar.o  -L/opt/local/lib -lgsl $(LIBS) -o evolverealscalar

evolverealscalar.o : evolverealscalar.c
	$(CC) evolverealscalar.c -c $(INCLUDES) -I/opt/local/include -o evolverealscalar.o

evolvemaxwell2D : evolvemaxwell2D.o
	$(CC) evolvemaxwell2D.o $(LIBS) -o evolvemaxwell2D

evolvemaxwell2D.o : evolvemaxwell2D.c
	$(CC) evolvemaxwell2D.c -c $(INCLUDES) -o evolvemaxwell2D.o

evolvescalar : evolvescalar.o
	$(CC) evolvescalar.o $(LIBS) -o evolvescalar

evolvescalar.o : evolvescalar.c
	$(CC) evolvescalar.c -c $(INCLUDES) -o evolvescalar.o

evolvecomplexscalar : evolvecomplexscalar.o
	$(CC) evolvecomplexscalar.o $(LIBS) -o evolvecomplexscalar

evolvecomplexscalar.o : evolvecomplexscalar.c
	$(CC) evolvecomplexscalar.c -c $(INCLUDES) -o evolvecomplexscalar.o

edw : edw.o
	$(CC) edw.o $(FFTLIB) $(LIBS) -o edw

edw.o : edw.c
	$(CC) edw.c -c $(CFLAGS) $(INCLUDES) -o edw.o
