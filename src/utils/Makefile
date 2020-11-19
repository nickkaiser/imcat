# makefile for utils

# be sure to have env variables
#	CC = some ansi c compiler
#	ARCH = whatever is appropriate

SHELL = /bin/sh

# object files:
OBJECTS = lu.o error.o arrays.o stats_stuff.o ran1.o nrutil.o \
fmedian.o smpopup.o linmodel.o brent.o f1dim.o frprmn.o linmin.o mnbrak.o \
powell.o modefunc.o fitstatus.o ipbuff.o args.o colormaps.o \
lmodel.o lmodelio.o lmodelcalculus.o gaussdev.o iostream.o

# headers
HEADERS = arrays.h error.h fmedian.h linmodel.h lu.h nrutil.h \
ran1.h smpopup.h stats_stuff.h brent.h f1dim.h frprmn.h linmin.h mnbrak.h \
powell.h modefunc.h fitstatus.h ipbuff.h args.h colormaps.h gaussdev.h \
lmodel.h iostream.h

DOCDIR = ../../doc/utils

all : libutils.a dummy_MAIN.o

libutils.a : $(OBJECTS)
	ar r libutils.a $(OBJECTS)

dummy_MAIN.o : dummy_MAIN.c

install : libutils.a dummy_MAIN.o $(RANLIB)
	@ echo +++++++++++ installing libutils.a in ../../lib/$(ARCH) +++++++++++
	@ if test -d ../../lib/$(ARCH) ; then \
		echo +++++++++++ dest-dir exists +++++++++++; else mkdir ../../lib/$(ARCH) ; fi
	cp libutils.a dummy_MAIN.o ../../lib/$(ARCH)

ranlib : libutils.a
	ranlib libutils.a


clean :
	rm -f *.o libutils.a

docs :
	@echo making html Man Pages for utils
	@if test ! -d $(DOCDIR) ; then mkdir $(DOCDIR) ; fi
	@cp docs_head.html docs.html
	@for h in $(HEADERS) ; do \
		echo $$h ; \
		echo '<h2><a name="'$$h'">'$$h'</a></h2><pre>' >> docs.html ; \
		cat $$h >> docs.html ; \
		echo '</pre><p>back to <a href="#INDEX">index</a><hr>' >> docs.html ; \
	done
	@cat docs_tail.html >> docs.html
	@cp docs.html $(DOCDIR)/docs.html
	makecatmanpage utils $(HEADERS) > $(MANDIR)/cat3/utils.3

lu.o : lu.c lu.h

error.o : error.c error.h

arrays.o :  arrays.c $(IMH)

fmedian.o : fmedian.c fmedian.h error.h

stats_stuff.o : stats_stuff.c $(IMH)

nrutil.o : nrutil.c

smpopup.o : smpopup.c smpopup.h

ran1.o : ran1.c

linmodel.o : linmodel.c linmodel.h lu.h

modefunc.o : modefunc.c modefunc.h

fitstatus.o : fitstatus.c fitstatus.h

args.o : args.c args.h

lmodel.o : lmodel.c lmodel.h

lmodelio.o : lmodelio.c lmodel.h

lmodelcalculus.o : lmodelcalculus.c lmodel.h

iostream.o : iostream.c iostream.h

iostreamtest : iostreamtest.o iostream.o
	gcc iostreamtest.o iostream.o error.o -o iostreamtest

