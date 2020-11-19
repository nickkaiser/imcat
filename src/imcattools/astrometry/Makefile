# Makefile for astrometry tools

TOOLS = readusnoacat readusnobcat getxsfromradec hmstodecimal decimaltohms \
dmstodecimal decimaltodms getworldcoords convertcoords getwcsinfo setwcsinfo
SCRIPTS = getUSNOAcatfromESOhtml

LIBS = -L$(IMCATDIR)/lib/$(ARCH) -lim -lcat -lutils -lm
INCLUDES = -I$(IMCATDIR)/src

BINDIR = $(IMCATDIR)/bin/$(ARCH)
SCRIPTDIR = $(IMCATDIR)/bin/scripts
DOC = $(IMCATDIR)/doc/astrometry.txt
DOCDIR = $(IMCATDIR)/doc/imcattools/astrometry

WCSOBJECTS = wcs.o platepos.o hget.o worldpos.o wcscon.o hput.o slasubs.o

SHELL = /bin/sh

all : $(TOOLS)

install : $(TOOLS)
	@echo installing astrometry tools in $(BINDIR) ....
	for i in $(SCRIPTS); do \
		if test -h $(SCRIPTDIR)/$$i ; then rm -f $(SCRIPTDIR)/$$i ; fi ; \
		ln -s ../../src/imcattools/astrometry/$$i $(SCRIPTDIR)/$$i ; \
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
	@echo making html Man Pages for imcattools_astrometry
	@if test ! -d $(DOCDIR) ; then mkdir -p $(DOCDIR) ; fi
	makemanpages imcattools_astrometry $(MANDIR) $(DOCDIR) $(TOOLS) $(SCRIPTS)


##
## the tools:
##


readusnoacat : readusnoacat.o getxfromradecfunc.o radecio.o
	$(CC) readusnoacat.o getxfromradecfunc.o radecio.o $(LIBS) -o readusnoacat

readusnobcat : readusnobcat.o getxfromradecfunc.o radecio.o
	$(CC) readusnobcat.o getxfromradecfunc.o radecio.o $(LIBS) -o readusnobcat

getxsfromradec : getxsfromradec.o getxfromradecfunc.o radecio.o
	$(CC) getxsfromradec.o getxfromradecfunc.o radecio.o $(LIBS) -o getxsfromradec

hmstodecimal : hmstodecimal.o radecio.o
	$(CC) hmstodecimal.o radecio.o $(LIBS) -o hmstodecimal

decimaltohms : decimaltohms.o radecio.o
	$(CC) decimaltohms.o radecio.o $(LIBS) -o decimaltohms

dmstodecimal : dmstodecimal.o radecio.o
	$(CC) dmstodecimal.o radecio.o $(LIBS) -o dmstodecimal

decimaltodms : decimaltodms.o radecio.o
	$(CC) decimaltodms.o radecio.o $(LIBS) -o decimaltodms

getworldcoords : getworldcoords.o getxfromradecfunc.o $(WCSOBJECTS)
	$(CC) getworldcoords.o getxfromradecfunc.o $(WCSOBJECTS) $(LIBS) -o getworldcoords

convertcoords : convertcoords.o radecio.o $(WCSOBJECTS)
	$(CC) convertcoords.o radecio.o $(WCSOBJECTS) $(LIBS) -o convertcoords

getwcsinfo : getwcsinfo.o $(WCSOBJECTS)
	$(CC) getwcsinfo.o $(WCSOBJECTS) $(LIBS) -o getwcsinfo

setwcsinfo : setwcsinfo.o $(WCSOBJECTS)
	$(CC) setwcsinfo.o $(WCSOBJECTS) $(LIBS) -o setwcsinfo


##
## the objects
##

readusnoacat.o : readusnoacat.c radecio.h
	$(CC) -c readusnoacat.c $(INCLUDES) -o readusnoacat.o

getxfromradecfunc.o : getxfromradecfunc.c
	$(CC) -c getxfromradecfunc.c $(INCLUDES) -o getxfromradecfunc.o

getxsfromradec.o : getxsfromradec.c radecio.h
	$(CC) -c getxsfromradec.c $(INCLUDES) -o getxsfromradec.o

radecio.o : radecio.c 
	$(CC) -c radecio.c $(INCLUDES) -o radecio.o

hmstodecimal.o : hmstodecimal.c 
	$(CC) -c hmstodecimal.c $(INCLUDES) -o hmstodecimal.o

decimaltohms.o : decimaltohms.c 
	$(CC) -c decimaltohms.c $(INCLUDES) -o decimaltohms.o

dmstodecimal.o : dmstodecimal.c 
	$(CC) -c dmstodecimal.c $(INCLUDES) -o dmstodecimal.o

decimaltodms.o : decimaltodms.c 
	$(CC) -c decimaltodms.c $(INCLUDES) -o decimaltodms.o

getworldcoords.o : getworldcoords.c 
	$(CC) -c getworldcoords.c $(INCLUDES) -o getworldcoords.o

convertcoords.o : convertcoords.c radecio.h
	$(CC) -c convertcoords.c $(INCLUDES) -o convertcoords.o

getwcsinfo.o : getwcsinfo.c 
	$(CC) -c getwcsinfo.c $(INCLUDES) -o getwcsinfo.o

setwcsinfo.o : setwcsinfo.c 
	$(CC) -c setwcsinfo.c $(INCLUDES) -o setwcsinfo.o
