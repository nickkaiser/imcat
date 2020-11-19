TOOLS = makelist transformlist pairs triplets matchtriplets matchpairs findcluster\
	mergelists
# SCRIPTS = scalerottrans
BINDIR = ../../bin/$(ARCH)
SCRIPTDIR = ../../bin/scripts
LIBS = -L../../lib/$(ARCH) -lim -lfft -lutils -lm
SHELL = /bin/sh
DOCDIR = ../../doc/dev

all : $(TOOLS)

install : $(TOOLS)
#	@ if test ! -d $(SCRIPTDIR) ; then  mkdir $(SCRIPT) ; fi
#	for i in $(SCRIPTS); do \
#		if test -f $(SCRIPTDIR)/$$i ; then rm -f $(SCRIPTDIR)/$$i ; fi ; \
#		ln -s ../../src/dev/$$i $(SCRIPTDIR)/$$i ; \
#	done
	@echo installing dev-tools in $(BINDIR)
	@ if test ! -d $(BINDIR) ; then  mkdir $(BINDIR) ; fi
	for i in $(TOOLS); do \
		cp $$i $(BINDIR)/$$i ; \
	done

clean:
	rm -f *.o core
	for i in $(TOOLS); do rm -f $$i ; done

docs :
	@echo making html Man Pages for dev
	@if test ! -d $(DOCDIR) ; then mkdir $(DOCDIR) ; fi
	makemanpages dev $(MANDIR) $(DOCDIR) $(TOOLS) $(SCRIPTS)

makelist : makelist.o
	$(CC) makelist.o -lm -o makelist

transformlist : transformlist.o
	$(CC) transformlist.o -lm -o transformlist

pairs : pairs.o
	$(CC) pairs.o -lm -o pairs

triplets : triplets.o
	$(CC) triplets.o -lm -o triplets

matchtriplets : matchtriplets.o
	$(CC) matchtriplets.o -lm -o matchtriplets

matchpairs : matchpairs.o
	$(CC) matchpairs.o -lm -o matchpairs

findcluster : findcluster.o
	$(CC) findcluster.o -lm -o findcluster

lintransformfit : lintransformfit.o
	$(CC) lintransformfit.o $(LIBS) -o lintransformfit

mosaicfit : mosaicfit.o
	$(CC) mosaicfit.o $(LIBS) -lm -o mosaicfit

mosaicmap : mosaicmap.o
	$(CC) mosaicmap.o $(LIBS) -lm -o mosaicmap

register : register.o
	$(CC) register.o $(LIBS) -lm -o register

mergelists : mergelists.o
	$(CC) mergelists.o -lm -o mergelists
