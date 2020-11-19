# Makefile for iis

all: iis

BINDIR = ../../../bin/$(ARCH)
LIB = -L../../../lib/$(ARCH) -lcdl -lim -lcat -lutils -lm $(SOCKET)
SHELL = /bin/sh
DOCDIR = ../../../doc/imtools/iis

OBJECTS = iis.o

install : iis
	@ echo ++++++++++++ installing iis in $(BINDIR) .... ++++++++++++
	@ if test -d $(BINDIR) ; then \
		echo ++++++ $(BINDIR) exists ++++++; else mkdir $(BINDIR) ; fi
	cp iis $(BINDIR)

iis : $(OBJECTS)
	$(CC) $(OBJECTS) $(LIB) -o iis

iis.o : iis.c

clean:
	rm -f core *.o iis

docs :
	@echo making html Man Pages for iis
	@if test ! -d $(DOCDIR) ; then mkdir $(DOCDIR) ; fi
	makemanpages imtools_iis $(MANDIR) $(DOCDIR) $(TOOLS) $(SCRIPTS)

