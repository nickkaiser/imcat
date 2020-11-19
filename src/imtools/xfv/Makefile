# makefile for xfv

# general libraries
LIBS = -lm

# the following library definitions will usually be overridden by top level Makefile
# X libs
XLIBS = -L/usr/X11R6/lib -lXt -lSM -lICE -lXext -lX11
# default athena widget libraries
XAWLIBS = -L/usr/X11R6/lib -lXaw -lXmu -lXext -lXt -lX11
# alternative athena widget libraries for aat solaris systems
# XAWLIBS = -L/usr/local/X11/lib -lXaw -lXmu -lXext -lXt -lX11 -lsocket
# alternative athena widget libraries for stego.ifa
# XAWLIBS = -L/usr/local/X11R6/lib -lXaw -lXmu -lXext -lXt -lX11

# imlib libraries
IMLIB = -ljpeg -lpng -ltiff -lz -lgif -lImlib
# IMLIB = -ljpeg -lpng -ltiff -lz -lgif -lImlib2

# imcat libraries
IMCATLIBDIR = ../../../lib/$(ARCH)
IMCATLIBS = -L$(IMCATLIBDIR) -lfft -lim -lutils

# includes
IMH = ../../imlib/fits.h ../../utils/error.h ../../utils/stats_stuff.h
XIMH = front.h back.h


BINDIR = ../../../bin/$(ARCH)
DOCDIR = ../../../doc/imtools/xfv

SHELL = /bin/sh

all : xfv

install : xfv
	@echo ++++++++++++ installing xfv in $(BINDIR) ....++++++++++++
	@ if test -d $(BINDIR) ; then \
		echo ++++++++++++ dest-dir exists ++++++++++++ ; else mkdir $(BINDIR) ; fi
	cp xfv $(BINDIR)

clean:
	rm -f *.o xfv

docs :
	@echo making html Man Pages for xfv
	@if test ! -d $(DOCDIR) ; then mkdir $(DOCDIR) ; fi
	makemanpages imtools_xfv $(MANDIR) $(DOCDIR) xfv

# straight X version
xfv : 	front.o back.o
	$(CC) front.o back.o $(IMCATLIBS) $(XAWLIBS) $(LIBS) -o xfv 

# Imlib version
# xfv : xfv.o
#	$(CC) xfv.o $(XLIBS) $(IMCATLIBS) $(IMLIB) $(LIBS) -o xfv

xfv.o : xfv.c
	$(CC) -c xfv.c $(CFLAGS) -o xfv.o -I$(IMCATDIR)/src
 
front.o :  front.c $(XIMH) $(IMH)

back.o : back.c $(XIMH) $(IMH)

color.o : color.c $(XIMH) $(IMH)

picclick.o : picclick.c $(XIMH) $(IMH)

