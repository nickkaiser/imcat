# Makefile for fitstops

all: fitstops

BINDIR = ../../../bin/$(ARCH)
LIB = -L../../../lib/$(ARCH) -lim -lutils -lm
SHELL = /bin/sh
DOCDIR = ../../../doc/imtools/fitstops

OBJECTS = fitstops.o makec.o psutils.o printimage.o contourplot.o surfaceplot.o tonry3d.o

install : fitstops
	@ echo ++++++++++++ installing fitstops in $(BINDIR) .... ++++++++++++
	@ if test -d $(BINDIR) ; then \
		echo ++++++ $(BINDIR) exists ++++++; else mkdir $(BINDIR) ; fi
	cp fitstops $(BINDIR)
	rm -f $(BINDIR)/print_image
	ln -s fitstops $(BINDIR)/print_image

fitstops : $(OBJECTS)
	$(CC) $(OBJECTS) $(LIB) -o fitstops


clean:
	rm -f core *.o fitstops

docs :
	@echo making html Man Pages for fitstops
	@if test ! -d $(DOCDIR) ; then mkdir $(DOCDIR) ; fi
	makemanpages imtools_fitstops $(MANDIR) $(DOCDIR) fitstops

