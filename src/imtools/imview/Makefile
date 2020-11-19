# makefile for xim

# general libraries
LIBS = -lm

# X-libraries
XMLIBS = -lXm -lXt -lX11 -lPW 

# imcat libraries
IMCATLIBDIR = ../../../lib/$(ARCH)
IMCATLIBS = -L$(IMCATLIBDIR) -lfft -lim -lutils

# includes
IMH = ../../imlib/fits.h ../../utils/error.h 

BINDIR = ../../../bin/$(ARCH)

SHELL = /bin/sh

all : imview

install : imview
	@echo installing imview in $(BINDIR) ....
	@ if test -d $(BINDIR) ; then \
		echo dest-dir exists ; else mkdir $(BINDIR) ; fi
	cp imview $(BINDIR)

clean:
	rm -f *.o imview

imview : imview.o imview_stuff.o
	$(CC) imview.o imview_stuff.o $(IMCATLIBS) $(XMLIBS) $(LIBS) -o imview

