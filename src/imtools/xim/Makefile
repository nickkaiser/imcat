# makefile for xim

# general libraries
LIBS = -lm

# default athena widget libraries
# XAWLIBS = -L/usr/local/X11/lib -lXaw -lXmu -lXext -lXt -lX11
# alternative athena widget libraries for aat solaris systems
# XAWLIBS = -L/usr/local/X11/lib -lXaw -lXmu -lXext -lXt -lX11 -lsocket
# alternative athena widget libraries for stego.ifa
# XAWLIBS = -L/usr/local/X11R6/lib -lXaw -lXmu -lXext -lXt -lX11

# imcat libraries
IMCATLIBDIR = ../../../lib/$(ARCH)
IMCATLIBS = -L$(IMCATLIBDIR) -lfft -lim -lutils

# includes
IMH = ../../imlib/fits.h ../../utils/error.h ../../utils/stats_stuff.h
XIMH = ximfront.h ximback.h ximbuttons.h ximcolor.h


BINDIR = ../../../bin/$(ARCH)
DOCDIR = ../../../doc/imtools/xim

SHELL = /bin/sh

all : xim

install : xim
	@echo ++++++++++++ installing xim in $(BINDIR) ....++++++++++++
	@ if test -d $(BINDIR) ; then \
		echo ++++++++++++ dest-dir exists ++++++++++++ ; else mkdir $(BINDIR) ; fi
	cp xim $(BINDIR)

clean:
	rm -f *.o xim

docs :
	@echo making html Man Pages for xim
	@if test ! -d $(DOCDIR) ; then mkdir $(DOCDIR) ; fi
	makemanpages imtools_xim $(MANDIR) $(DOCDIR) xim

xim : 	ximfront.o ximback.o ximbuttons.o ximcolor.o print_stuff.o
	$(CC) ximfront.o ximback.o ximbuttons.o ximcolor.o print_stuff.o \
	$(IMCATLIBS) $(XAWLIBS) $(LIBS) -o xim 


ximfront.o :  ximfront.c $(XIMH) $(IMH)

ximback.o : ximback.c $(XIMH) $(IMH)

ximbuttons.o : ximbuttons.c $(XIMH) $(IMH)

ximcolor.o :  ximcolor.c $(XIMH) $(IMH)

print_stuff.o : print_stuff.c
