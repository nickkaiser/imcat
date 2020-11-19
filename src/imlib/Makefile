# Makefile for libim.a

# be sure to have environment variables:
#	CC = gcc
#	ARCH = whatever is appropriate

SHELL = /bin/sh

OBJECTS = fits.o convertarray.o filters.o map.o scrunch_stuff.o sm_stuff.o \
	fastmap.o
INCLUDES = ../utils/error.h ../fftlib/myfft.h ../utils/arrays.h fits.h

HEADERS = fits.h convertarray.h filters.h map.h map_private.h sm_stuff.h scrunch_stuff.h
DOCDIR = ../../doc/imlib

all : libim.a

libim.a : $(OBJECTS)
	ar r libim.a $?

ranlib : libim.a
	ranlib libim.a

install : libim.a $(RANLIB)
	@ echo +++++++++++ installing libim.a in ../../lib/$(ARCH) +++++++++++
	@ if test -d ../../lib/$(ARCH) ; then \
		echo  +++++++++++ dest-dir exists  +++++++++++ ; else mkdir ../../lib/$(ARCH) ; fi
	cp libim.a ../../lib/$(ARCH)

docs :
	@echo making html Man Pages for libim
	@if test ! -d $(DOCDIR) ; then mkdir $(DOCDIR) ; fi
	@cp docs_head.html docs.html
	@for h in $(HEADERS) ; do \
		echo '<h2><a name="'$$h'">'$$h'</a></h2><pre>' >> docs.html ; \
		cat $$h >> docs.html ; \
		echo '</pre><p>back to <a href="#INDEX">index</a><hr>' >> docs.html ; \
	done
	@cat docs_tail.html >> docs.html
	@cp docs.html $(DOCDIR)/docs.html
	makecatmanpage imlib $(HEADERS) > $(MANDIR)/cat3/imlib.3

clean :
	rm -f *.o libim.a	


fastmap.o : fastmap.c map.h

map.o : map.c map.h map_private.h $(INCLUDES)

fits.o : fits.c $(INCLUDES)

filters.o : filters.c filters.h $(INCLUDES)

scrunch_stuff.o : scrunch_stuff.c scrunch_stuff.h $(INCLUDES)

sm_stuff.o : $(INCLUDES)



LIBS = -L../../lib/$(ARCH) -lim -lutils -lm

maptest : maptest.o
	$(CC) -o maptest maptest.o $(LIBS)

maptest.o : maptest.c $(INCLUDES)

