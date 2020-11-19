# master Makefile for the imcat stuff

# be sure to define the environment variables $(IMCATDIR), $(CC) and $(ARCH)
# IMCATDIR is where you install imcat
# CC should be an ansi compatible compiler (e.g. gcc)
# ARCH should identify the architecture

# I have defined various default macros for library locations etc
# in this file.  You can override these by modifying this makefile

# typical build:
#	> make libs			# make the libraries
#	> make libs MODE=install	# install them
#	> make tools			# make the tools
#	> make tools MODE=install	# install them
#	> make docs			# make the html man pages
# you will need to install imcat/src/doc/imcat_search.pl in your
# local cgi-bin directory by hand if you want to use it.
# if all goes well, to clean up, do
#	> make MODE=clean
# to nuke the libraries and executables for current architecture, do
#	> make empty
# to nuke all the libraries, executables and scripts
#	> make empty_all

# default X athena widget libraries
XAWLIBS = "-L/usr/local/X11/lib -lXaw -lXmu -lXext -lXt -lX11"
# alternatives:
# for aao solaris machines I used
# XAWLIBS = -L/usr/local/X11/lib -lXaw -lXmu -lXext -lXt -lX11 -lsocket
# for stego.ifa I used
# XAWLIBS = -L/usr/local/X11R6/lib -lXaw -lXmu -lXext -lXt -lX11

# default X motif libraries (only used for imview)
XMLIBS = "-lXm -lXt -lX11 -lPW"

# default FFT is numerical recipes (very slow + inefficient)
# IMCATFFT = NR
# IMCATFFTLIB=
# on IRIX I use
IMCATFFT = SGI
IMCATFFTLIB = "-L$(IMCATDIR)/lib/sgi_fft -lsgi_fft"



# defaults for imcattools/pgplot are a bit trickier

# first we specify the pgplot libraries
PGLIBDIR = $(IMCATDIR)/lib/$(ARCH)/pgplot/

# we also need to tell the c-compiler where to find the fortran
# libraries.  To get these do 'f77 -v temp.f' and peruse the output!
# for OSF1 this is
# FLIBS = "-lUfor -lfor -lFutil -lm -lots"
# for IRIX this is
FLIBS = "-lF77 -lm -lU77 -lI77 -lisam -lmpc"
# for aat sun4-solaris I use gcc but linked with f77 so have to specify the library
# CLIBS = "-L/usr/local/lib/gcc-lib/sparc-sun-solaris2.4/2.7.2 -lgcc"

# finally we need to specify the X-libs
# at CITA we just use
XLIBS = -lX11
# though we used to have to prepend '-lgl' for sgi machines
# and for aat solaris machines I used
# XLIBS = "-L/opt/X11R5/lib -lX11 -lsocket"



# default ranlib
RANLIB = ranlib
# if you don;t have ranlib then use
# RANLIB =


# you should not need to change anything below

# directory for html format man pages (mine is actually a symbolic link to my
# www space)
DOCDIR = ../doc

SHELL = /bin/sh

all : libs tools

libs : utils_ fftlib_ imlib_ catlib_

tools : cattools_ cattools_lc cattools_mergecats\
	dev_ nk_ nk_tablestuff \
	imtools_ imtools_xim imtools_iis imtools_ic \
	imcattools_objectfinding imcattools_mosaics \
	imcattools_psfcorrection imcattools_registration \
	imcattools_shearanal imcattools_pgplot

empty :
	rm -rf ../bin/$(ARCH)/*
	rm -rf ../lib/$(ARCH)/*

empty_all :
	rm -rf ../bin/*
	rm -rf ../lib/*

docs :
	@echo copying doc/mainindex.html to $(DOCDIR)
	@if test ! -d $(DOCDIR) ; then mkdir $(DOCDIR) ; fi
	cp doc/mainindex.html $(DOCDIR)/mainindex.html
	cp doc/imcat_search.html $(DOCDIR)/imcat_search.html
	make tools MODE=docs
	make libs MODE=docs
	@echo don't forget to install imcat_search.pl!

# the libraries:
utils_ :
	cd utils ; make $(MODE) RANLIB=$(RANLIB)

fftlib_ :
	cd fftlib ; make $(MODE) RANLIB=$(RANLIB) FFT=$(IMCATFFT)

catlib_ :
	cd catlib ; make $(MODE) RANLIB=$(RANLIB)

imlib_ :
	cd imlib ; make $(MODE) RANLIB=$(RANLIB)



# the tools:
cattools_ :
	cd cattools ; make $(MODE)

cattools_lc :
	cd cattools/lc ; make $(MODE)

cattools_tk :
	cd cattools/tk ; make $(MODE)

cattools_mergecats :
	cd cattools/mergecats ; make $(MODE)



dev_ :
	cd dev ; make $(MODE)

nk_ :
	cd nk ; make $(MODE)

nk_tablestuff :
	cd nk/tablestuff ; make $(MODE)



imtools_ :
	cd imtools ; make $(MODE) FFTLIB=$(IMCATFFTLIB)

imtools_xim :
	cd imtools/xim ; make $(MODE) XAWLIBS=$(XAWLIBS)

imtools_iis :
	cd imtools/iis ; make $(MODE)

imtools_imview :
	cd imtools/imview ; make $(MODE) XMLIBS=$(XMLIBS)

imtools_ic :
	cd imtools/ic ; make $(MODE)



imcattools_objectfinding :
	cd imcattools/objectfinding ; make $(MODE) FFTLIB=$(IMCATFFTLIB)

imcattools_mosaics :
	cd imcattools/mosaics ; make $(MODE)

imcattools_psfcorrection :
	cd imcattools/psfcorrection ; make $(MODE)

imcattools_registration :
	cd imcattools/registration ; make $(MODE) FFTLIB=$(IMCATFFTLIB)

imcattools_shearanal :
	cd imcattools/shearanal ; make $(MODE) FFTLIB=$(IMCATFFTLIB)

imcattools_pgplot :
	cd imcattools/pgplot ; make $(MODE) XLIBS=$(XLIBS) FLIBS=$(FLIBS) \
		PGLIBDIR=$(PGLIBDIR)



