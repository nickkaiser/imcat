# makefile for photometry tools

TOOLS = getsky apphot getshapes getfitsval makekernel makegausskernel \
	getshapes2 getshapes3 getshapes1 modelpsf profile deproject \
	makedensity makepeff makepeff2 getdonutstats
SCRIPTS =

#libraries
LIBS = -L../../../lib/$(ARCH) -lim -lfft $(FFTLIB) -lcat -lutils $(FLIBS) -lm


BINDIR = ../../../bin/$(ARCH)
SCRIPTDIR = ../../../bin/scripts
DOCDIR = ../../../doc/imcattools/photometry

SHELL = /bin/sh

all : $(TOOLS)

install : $(TOOLS)
	@echo installing imcattools in $(BINDIR) ....
#	for i in $(SCRIPTS); do \
#		if test -h $(SCRIPTDIR)/$$i ; then rm -f $(SCRIPTDIR)/$$i ; fi ; \
#		ln -s ../../src/imcattools/objectfinding/$$i $(SCRIPTDIR)/$$i ; \
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
	@echo making html Man Pages for imcattools_photometry
	@if test ! -d $(DOCDIR) ; then mkdir -p $(DOCDIR) ; fi
	makemanpages  imcattools_photometry $(MANDIR) $(DOCDIR) $(TOOLS) $(SCRIPTS)

##
## the tools:
##


getsky : getskymain.o getsky.o zap.o
	$(CC) getskymain.o getsky.o zap.o $(LIBS) -o getsky

apphot : apphotmain.o apphot.o zap.o
	$(CC) apphotmain.o apphot.o zap.o $(LIBS) -o apphot

getshapes : getshapesmain.o getshape.o zap.o
	$(CC) getshapesmain.o getshape.o zap.o $(LIBS) -o getshapes

getfitsval : getfitsval.o
	$(CC) getfitsval.o $(LIBS) -o getfitsval

makekernel : makekernel.o
	$(CC) makekernel.o $(LIBS) -o makekernel

makegausskernel : makegausskernel.o
	$(CC) makegausskernel.o $(LIBS) -o makegausskernel

getshapes2 : getshapes2.o getshape2.o zap.o
	$(CC) getshapes2.o getshape2.o zap.o $(LIBS) -o getshapes2

getshapes3 : getshapes3.o getshape3.o zap.o
	$(CC) getshapes3.o getshape3.o zap.o $(LIBS) -o getshapes3

getshapes1 : getshapes1main.o getshape1.o zap.o
	$(CC) getshapes1main.o getshape1.o zap.o $(LIBS) -o getshapes1

modelpsf : modelpsf.o modelpsf.o
	$(CC) modelpsf.o $(LIBS) -o modelpsf

profile : profile.o
	$(CC) -o profile profile.o $(LIBS)

profile.o : profile.c
	$(CC) -c profile.c -o profile.o

deproject : deproject.o
	$(CC) -o deproject deproject.o $(LIBS)

deproject.o : deproject.c
	$(CC) -c deproject.c -o deproject.o

makepeff : makepeff.o
	$(CC) makepeff.o $(LIBS) -o makepeff

makepeff2 : makepeff2.o
	$(CC) makepeff2.o $(LIBS) -o makepeff2

makedensity : makedensity.o
	$(CC) makedensity.o $(LIBS) -o makedensity

getdonutstats : getdonutstats.o
	$(CC) getdonutstats.o $(LIBS) -o getdonutstats













