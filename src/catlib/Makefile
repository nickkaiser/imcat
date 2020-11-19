# Makefile for libcat.a

# be sure to have environment variables:
#	CC = gcc
#	ARCH = whatever is appropriate

SHELL = /bin/sh
HEADERS = cat.h error.h
DOCDIR = ../../doc/catlib

all : libcat.a

libcat.a : cat.o
	ar r libcat.a cat.o

install : libcat.a $(RANLIB)
	@ echo ++++++++++++ installing libcat.a in ../../lib/$(ARCH) ++++++++++++
	@ if test -d ../../lib/$(ARCH) ; then \
		echo ++++++++++++ dest-dir ++++++++++++ exists ; else mkdir ../../lib/$(ARCH) ; fi
	cp libcat.a ../../lib/$(ARCH)

ranlib : libcat.a
	ranlib libcat.a


docs :
	@echo making html Man Pages for catlib
	@if test ! -d $(DOCDIR) ; then mkdir $(DOCDIR) ; fi
	@cp docs_head.html docs.html
	@for h in $(HEADERS) ; do \
		echo '<h2><a name="'$$h'">'$$h'</a></h2><pre>' >> docs.html ; \
		cat $$h >> docs.html ; \
		echo '</pre><p>back to <a href="#INDEX">index</a><hr>' >> docs.html ; \
	done
	@cat docs_tail.html >> docs.html
	@cp docs.html $(DOCDIR)/docs.html
	makecatmanpage catlib $(HEADERS) > $(MANDIR)/cat3/catlib.3

clean :
	rm -f cat.o libcat.a	
