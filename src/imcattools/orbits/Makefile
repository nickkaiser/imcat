# Makefile for orbits tools

TOOLS = airmass convertorbits evolveorbits fitorbit getdet getmoid laplace3 \
 makeobs_circ makeobs_inertial orbanim r2n orbs2obs planetephem
SCRIPTS = 

LIBS = -L$(IMCATDIR)/lib/$(ARCH) -lim -lutils -lm
INCLUDES = -I$(IMCATDIR)/src

# MYORBITLIBS = kepler_heasley/libkepler.a orbitutils/liborbitutils.a
MYORBITLIBS = orbitutils/liborbitutils.a
# ORBITLIBS = $(MYORBITLIBS) $(FLIBS)
ORBITLIBS = $(MYORBITLIBS) -L$(ORBFITLIBDIR) -lsuit -lm
# ORBITINCLUDES = -Ikepler_heasley -Iorbitutils
ORBITINCLUDES = -Iorbitutils

BINDIR = $(IMCATDIR)/bin/$(ARCH)
SCRIPTDIR = $(IMCATDIR)/bin/scripts
DOCDIR = $(IMCATDIR)/doc/imcattools/orbits

SHELL = /bin/sh

all : $(TOOLS)

install : $(TOOLS)
	@echo installing orbits tools in $(BINDIR) ....
#	for i in $(SCRIPTS); do \
#		if test -h $(SCRIPTDIR)/$$i ; then rm -f $(SCRIPTDIR)/$$i ; fi ; \
#		ln -s ../../src/imcattools/orbits/$$i $(SCRIPTDIR)/$$i ; \
#	done
#	@ if test -d $(BINDIR) ; then \
#		echo $(BINDIR) exists ; else mkdir $(BINDIR) ; fi
	for i in $(TOOLS); do \
		cp $$i $(BINDIR)/$$i ; \
	done

clean:
	rm -f core *.o
	for i in $(TOOLS); do rm -f $$i ; done
	cd kepler_heasley ; make clean ; make empty
	cd orbitutils ; make clean ; make empty

docs :
	@echo making Man Pages for imcattools_orbits
	@if test ! -d $(DOCDIR) ; then mkdir -p $(DOCDIR) ; fi
	makemanpages imcattools_orbits $(MANDIR) $(DOCDIR) $(TOOLS) $(SCRIPTS)

# the libs
# kepler_heasley/libkepler.a :
#	cd kepler_heasley ; make libkepler.a

orbitutils/liborbitutils.a :
	cd orbitutils ; make liborbitutils.a

##
## the tools:
##

airmass : airmass.o $(MYORBITLIBS)
	gcc airmass.o $(ORBITLIBS) -o airmass

airmass.o : airmass.c
	gcc -c airmass.c $(ORBITINCLUDES) -o airmass.o

makeobs_inertial : makeobs_inertial.o
	gcc makeobs_inertial.o $(ORBITLIBS) -o makeobs_inertial

makeobs_inertial.o : makeobs_inertial.c
	gcc -c makeobs_inertial.c $(ORBITINCLUDES) -o makeobs_inertial.o

makeobs_circ : makeobs_circ.o
	gcc makeobs_circ.o $(ORBITLIBS) -o makeobs_circ

makeobs_circ.o : makeobs_circ.c
	gcc -c makeobs_circ.c $(ORBITINCLUDES) -o makeobs_circ.o

laplace3 : laplace3.o
	gcc laplace3.o $(ORBITLIBS) -o laplace3

laplace3.o : laplace3.c
	gcc -c laplace3.c $(ORBITINCLUDES) -o laplace3.o

tcl_evolve : tcl_evolve.o
	gcc tcl_evolve.o $(ORBITLIBS) -o tcl_evolve

tcl_evolve.o : tcl_evolve.c
	gcc -c tcl_evolve.c $(ORBITINCLUDES) -o tcl_evolve.o

r2n : r2n.o
	gcc r2n.o $(ORBITLIBS) -o r2n

r2n.o : r2n.c
	gcc -c r2n.c $(ORBITINCLUDES) -o r2n.o

fitorbit : fitorbit.o
	g77 -g -C fitorbit.o $(ORBITLIBS) -o fitorbit

fitorbit.o : fitorbit.c 
	gcc -c fitorbit.c $(ORBITINCLUDES) -o fitorbit.o

tcl_evolveN : tcl_evolveN.o
	gcc tcl_evolveN.o $(LIBS) $(ORBITLIBS) -o tcl_evolveN

tcl_evolveN.o : tcl_evolveN.c
	gcc -c tcl_evolveN.c $(INCLUDES) $(ORBITINCLUDES) -o tcl_evolveN.o

evolveorbits : evolveorbits.o
	g77 -g -C evolveorbits.o $(LIBS) $(ORBITLIBS) -o evolveorbits

evolveorbits.o : evolveorbits.c
	gcc -c evolveorbits.c $(INCLUDES) $(ORBITINCLUDES) -o evolveorbits.o

orbanim : orbanim.o
	gcc orbanim.o $(LIBS) $(ORBITLIBS) -o orbanim

orbanim.o : orbanim.c
	gcc -c orbanim.c $(INCLUDES) $(ORBITINCLUDES) -o orbanim.o

orbs2obs : orbs2obs.o
	gcc orbs2obs.o $(ORBITLIBS) -o orbs2obs

orbs2obs.o : orbs2obs.c
	gcc -c orbs2obs.c $(ORBITINCLUDES) -o orbs2obs.o

getdet : getdet.o
	gcc getdet.o $(ORBITLIBS) -o getdet

getdet.o : getdet.c
	gcc -c getdet.c $(ORBITINCLUDES) -o getdet.o

convertorbits : convertorbits.o
	g77 -g -C convertorbits.o $(ORBITLIBS) -o convertorbits

convertorbits.o : convertorbits.c
	gcc -c convertorbits.c $(ORBITINCLUDES) -o convertorbits.o

getmoid : getmoid.o
	g77 -g -C getmoid.o $(ORBITLIBS) -o getmoid

getmoid.o : getmoid.c
	gcc -c getmoid.c $(ORBITINCLUDES) -o getmoid.o

planetephem : planetephem.o orbitutils/jpleph.o
	g++ planetephem.o orbitutils/jpleph.o -lm -o planetephem

planetephem.o : planetephem.c
	$(CC) -c planetephem.c $(ORBITINCLUDES) -o planetephem.o


