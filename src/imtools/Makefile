# makefile for imtools

TOOLS = acf album asciitofits boxavg bscale changesex\
	colorbox colorize combinestamps combineimages complex convolve cycleimage\
	fft fillcores fitstoascii fitstocat fitstomovie flatten flip getplane \
	getplanes grad growmagic helicalscan ic_stream imhead invgrad kappa2stuff linesmode\
	lensmap makecolourimage make_image makesubimage\
	mapbynumericdef mask maskfits overscancorr pixdist projectfits \
	scrunch simulate sliceimage smooth spinflip stackplanes stats transformimage\
	translate_fft unscrunch warpimage transposebits unpackextensions getimagemoments

# discontinued: read4shooter rescale
SCRIPTS = logscaleimage statstocat


# general libraries
# SGI has fast malloc
# LIBS = $(FFTLIB) -lm -lmalloc
# DEC
# LIBS = $(FFTLIB) -lm
# workaround for suns
LIBS = -L../../lib/$(ARCH) -lim -lfft $(FFTLIB) -lutils $(FLIBS) -lm

OBJ = 

# groups of files for dependencies
IMH = ../utils/error.h ../utils/arrays.h ../imlib/filters.h \
../utils/stats_stuff.h ../fftlib/myfft.h


# sm libraries:
# sm libraries for CITA SGIs
# SMLIBS = -L/usr/local/lib -lSMplotsub -lSMdevices -lSMutils -lgl -lX11 
# sm libraries for hare
# SMLIBS = -L/usr/local/lib -lSMplotsub -lSMdevices -lSMutils -L/usr/local/X11/lib -lX11 
# sm libraries for bobcat
# SMLIBS = -L/usr/local/lib -lplotsub -ldevices -lutils -L/usr/local/X11/lib -lX11 
# sm libraries for CITA SUNs
# SMLIBS = -L/usr/local/lib -lSMplotsub -lSMdevices -lSMutils  -lX11 \
# 	 -lsuntool -lsunwindow -lpixrect
# sm libraries for CFHT
#SMLIBS = -L/usr/local/lib/sm -lplotsub -ldevices -lutils  -lX11 \
#        -lsuntool -lsunwindow -lpixrect
# sm libraries for stego at ifa
# SMLIBS = libplotsub.a libdevices.a libutils.a  -lX11 \
#        -lsuntool -lsunwindow -lpixrect
# sm libraries for AAO SUNs
SMLIBS = -L/usr/local/lib -lplotsub -ldevices -lutils  -lX11 \
	 -lsuntool -lsunwindow -lpixrect

BINDIR = ../../bin/$(ARCH)
SCRIPTDIR = ../../bin/scripts
DOCDIR = ../../doc/imtools

SHELL = /bin/sh

all : $(TOOLS)

install : $(TOOLS)
	@echo installing imtool scripts in $(SCRIPTDIR) ....
	for i in $(SCRIPTS); do \
		if test -h $(SCRIPTDIR)/$$i ; then rm -f $(SCRIPTDIR)/$$i ; fi ; \
		ln -s ../../src/imtools/$$i $(SCRIPTDIR)/$$i ; \
	done
	@echo installing imtools in $(BINDIR) ....
	@ if test -d $(BINDIR) ; then \
		echo $(BINDIR) exists ; else mkdir $(BINDIR) ; fi
	for i in $(TOOLS); do \
		cp $$i $(BINDIR)/$$i ; \
	done

clean:
	rm -f *.o core
	for i in $(TOOLS); do rm -f $$i ; done

docs :
	@echo making Man Pages for Imtools
	@if test ! -d $(DOCDIR) ; then mkdir $(DOCDIR) ; fi
	makemanpages imtools $(MANDIR) $(DOCDIR) $(TOOLS) $(SCRIPTS)


##
## the tools:
##

getimagemoments :  $(OBJ) getimagemoments.o
	$(CC) $(OBJ) getimagemoments.o $(LIBS) -o getimagemoments

make_image :  $(OBJ) make_image.o
	$(CC) $(OBJ) make_image.o $(LIBS) -o make_image

simulate :  $(OBJ) simulate.o
	$(CC) $(OBJ) simulate.o $(LIBS) -o simulate

print_image : $(OBJ) print_image.o
	$(CC) print_image.o $(OBJ) $(LIBS) -o print_image

smooth :  $(OBJ) smooth.o
	$(CC) $(OBJ) smooth.o $(LIBS) -o smooth

stats :  stats.o $(OBJ)
	$(CC) $(OBJ) stats.o $(LIBS) -o stats

scrunch :	scrunch.o $(OBJ)
	$(CC) scrunch.o $(OBJ) $(LIBS) -o scrunch 

imhead :	imhead.o $(OBJ)
	$(CC) imhead.o $(OBJ) $(LIBS) -o imhead

mask :	mask.o $(OBJ) 
	$(CC) mask.o $(OBJ) $(LIBS) -o mask 

addfitscomment : addfitscomment.o $(OBJ) 
	$(CC) addfitscomment.o $(OBJ) $(LIBS) -o addfitscomment

combineimages : combineimages.o avgpixstack.o $(OBJ)
	$(CC) combineimages.o avgpixstack.o $(OBJ) $(LIBS) -o combineimages 

combinestamps : combinestamps.o combine_stuff.o $(OBJ)
	$(CC) combinestamps.o combine_stuff.o $(OBJ) $(LIBS) -o combinestamps 

makecolourimage :	makecolourimage.o $(OBJ)
	$(CC) makecolourimage.o $(OBJ) $(LIBS) -o makecolourimage 

colorbox :	colorbox.o $(OBJ)
	$(CC) colorbox.o $(OBJ) $(LIBS) -o colorbox 

colorize :      colorize.o $(OBJ)
	$(CC) colorize.o $(OBJ) $(LIBS) -o colorize

complex :      complex.o $(OBJ)
	$(CC) complex.o $(OBJ) $(LIBS) -o complex

makesubimage :	makesubimage.o $(OBJ)
	$(CC) makesubimage.o $(OBJ) $(LIBS) -o makesubimage

unscrunch :	unscrunch.o $(OBJ) 
	$(CC) unscrunch.o $(OBJ) $(LIBS) -o unscrunch 

maskfits :	maskfits.o $(OBJ) 
	$(CC) maskfits.o $(OBJ) $(LIBS) -o maskfits 

contour :	contour.o $(OBJ) 
	$(CC) contour.o $(OBJ) $(LIBS) -o contour 

profile :	profile.o $(OBJ) 
	$(CC) profile.o $(OBJ) $(LIBS) -o profile 

lensmap :	lensmap.o $(OBJ) 
	$(CC) lensmap.o  $(OBJ) $(LIBS) -o lensmap 

spinflip :	spinflip.o $(OBJ) 
	$(CC) spinflip.o $(OBJ) $(LIBS) -o spinflip 

fitstoascii :	fitstoascii.o $(OBJ) 
	$(CC) fitstoascii.o $(OBJ) $(LIBS) -o fitstoascii 

fitstocat :	fitstocat.o $(OBJ) 
	$(CC) fitstocat.o $(OBJ) $(LIBS) -o fitstocat 

asciitofits :	asciitofits.o $(OBJ) 
	$(CC) asciitofits.o $(OBJ) $(LIBS) -o asciitofits 

read4shooter :	read4shooter.o $(OBJ) 
	$(CC) read4shooter.o $(OBJ) $(LIBS) -o read4shooter 

album :	album.o $(OBJ) 
	$(CC) album.o $(OBJ) $(LIBS) -o album 

acf :	acf.o $(OBJ) 
	$(CC) acf.o $(OBJ) $(LIBS) -o acf

pixdist :   pixdist.o $(OBJ)
	$(CC) pixdist.o $(OBJ) $(LIBS) -o pixdist 

transformimage :   transformimage.o $(OBJ)
	$(CC) transformimage.o $(OBJ) $(LIBS) -o transformimage 

changesex :   changesex.o $(OBJ)
	$(CC) changesex.o $(OBJ) $(LIBS) -o changesex 

bscale :   bscale.o $(OBJ)
	$(CC) bscale.o $(OBJ) $(LIBS) -o bscale 

kappa2stuff :   kappa2stuff.o $(OBJ)
	$(CC) kappa2stuff.o $(OBJ) $(LIBS) -o kappa2stuff

overscancorr : $(OBJ) overscancorr.o
	$(CC) $(OBJ) overscancorr.o $(LIBS) -o overscancorr

fft : fft.o
	$(CC) fft.o $(LIBS) -o fft

convolve : convolve.o
	$(CC) convolve.o $(LIBS) -o convolve

cycleimage : cycleimage.o
	$(CC) cycleimage.o $(LIBS) -o cycleimage

linesmode : linesmode.o
	$(CC) linesmode.o $(LIBS) -o linesmode


grad : grad.o
	$(CC) grad.o $(LIBS) -o grad

invgrad : invgrad.o
	$(CC) invgrad.o $(LIBS) -o invgrad


translate_fft : translate_fft.o
	$(CC) translate_fft.o $(LIBS) -o translate_fft

mapbynumericdef : mapbynumericdef.o
	$(CC) mapbynumericdef.o $(LIBS) -o mapbynumericdef

warpimage : warpimage.o
	$(CC) warpimage.o $(LIBS) -o warpimage

fft_task : fft_task.o
	$(CC) fft_task.o $(LIBS) -o fft_task

fillcores : $(OBJ) fillcores.o
	$(CC) $(OBJ) fillcores.o $(LIBS) -o fillcores

flip :	flip.o $(OBJ) 
	$(CC) flip.o $(OBJ) $(LIBS) -o flip 

flatten :	flatten.o $(OBJ) 
	$(CC) flatten.o $(OBJ) $(LIBS) -o flatten

growmagic :	growmagic.o $(OBJ) 
	$(CC) growmagic.o $(OBJ) $(LIBS) -o growmagic

getplane :      getplane.o $(OBJ)
	$(CC) getplane.o $(OBJ) $(LIBS) -o getplane

getplanes :      getplanes.o $(OBJ)
	$(CC) getplanes.o $(OBJ) $(LIBS) -o getplanes

stackplanes : stackplanes.o $(OBJ)
	$(CC) stackplanes.o $(OBJ) $(LIBS) -o stackplanes

transposebits : transposebits.o $(OBJ)
	$(CC) transposebits.o $(OBJ) $(LIBS) -o transposebits

boxavg : boxavg.o $(OBJ)
	$(CC) boxavg.o $(OBJ) $(LIBS) -o boxavg

projectfits : projectfits.o $(OBJ)
	$(CC) projectfits.o $(OBJ) $(LIBS) -o projectfits

unpackextensions : unpackextensions.o $(OBJ)
	$(CC) unpackextensions.o $(OBJ) $(LIBS) -o unpackextensions

helicalscan : helicalscan.o
	$(CC) helicalscan.o $(LIBS) -o helicalscan

fitstomovie : fitstomovie.o
	$(CC) fitstomovie.o $(LIBS) -o fitstomovie

sliceimage : sliceimage.o
	$(CC) sliceimage.o $(LIBS) -o sliceimage

ic_stream : ic_stream.o
	$(CC) ic_stream.o $(LIBS) -o ic_stream


##
## imtool objects
##

make_image.o : make_image.c make_image.h $(IMH)

simulate.o : simulate.c simulate.h $(IMH)

print_image.o : print_image.c $(IMH)

smooth.o : smooth.c $(IMH)

stats.o : stats.c $(IMH)

scrunch.o :	scrunch.c ../imlib/scrunch_stuff.h $(IMH)

imhead.o :	imhead.c $(IMH)

mask.o :	mask.c $(IMH)

pastiche.o :	pastiche.c $(IMH)

makecolourimage.o :	makecolourimage.c $(IMH)

makesubimage.o :	makesubimage.c $(IMH)

unscrunch.o :	unscrunch.c $(IMH)

maskfits.o :	maskfits.c $(IMH)

contour.o :	contour.c ../imlib/sm_stuff.h $(IMH)

profile.o :	profile.c ../imlib/sm_stuff.h $(IMH)

lensmap.o : lensmap.c lensmap.h ../imlib/map.h ../imlib/map_private.h $(IMH)

spinflip.o : spinflip.c $(IMH)

fitstoascii.o : fitstoascii.c $(IMH)

fitstocat.o : fitstocat.c $(IMH)

asciitofits.o : asciitofits.c $(IMH)

read4shooter.o : read4shooter.c $(IMH)

album.o : album.c $(IMH)

acf.o : acf.c $(IMH)

pixdist.o : pixdist.c $(IMH)

transformimage.o : transformimage.c ../imlib/map.h $(IMH)

changesex.o : changesex.c $(IMH)

bscale.o : bscale.c $(IMH)

kappa2stuff.o : kappa2stuff.c $(IMH)

mapbynumericdef.o : mapbynumericdef.c $(IMH)

warpimage.o : warpimage.c $(IMH)

flip.o : flip.c $(IMH)

flatten.o : flatten.c $(IMH)

growmagic.o : growmagic.c $(IMH)

getplane.o : getplane.c $(IMH)

getplanes.o : getplanes.c $(IMH)

stackplanes.o : stackplanes.c $(IMH)

transposebits.o : transposebits.c $(IMH)

complex.o : complex.c $(IMH)

unpackextensions.o : unpackextensions.c $(IMH)

helicalscan.o : helicalscan.c $(IMH)

fitstomovie.o : fitstomovie.c mpegtemplate.h $(IMH)

sliceimage.o : sliceimage.c $(IMH)

ic_stream.o : ic_stream.c $(IMH)






