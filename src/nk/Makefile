# makefile for nk directory scripts

SCRIPTS = notabs reverse translateps bundle overwrite rotateps tfilter \
	travesty doto pscat sendbundle thist truepath echoeval remote stree Top tsort \
	 fixpgplotps tproc thist xevtocat clickstocat \
	makeimcatman makemanpages makecatmanpage makegroffmanpage makehtmlmanpage maketexmanpage

SCRIPTSWITHMANPAGES = makeimcatman makemanpages makecatmanpage makegroffmanpage \
	makehtmlmanpage maketexmanpage doto xevtocat clickstocat


SHELL = /bin/sh
DOCDIR = ../../doc/nk

all:
	@ echo +++++++++++++ this directory only contains scripts +++++++++++++

BINDIR = ../../bin/scripts

install :
	if test -d $(BINDIR) ; then echo ++++++ $(BINDIR) exists ++++ ; else mkdir $(BINDIR) ; fi ;
	for i in $(SCRIPTS); do \
		if test -h $(BINDIR)/$$i ; then rm $(BINDIR)/$$i ; fi ; \
		ln -s ../../src/nk/$$i $(BINDIR)/$$i ; \
	done


docs :
	@if test ! -d $(DOCDIR) ; then mkdir $(DOCDIR) ; fi
	@cp aaaindex.html $(DOCDIR)/aaaindex.html
	makemanpages nk $(MANDIR) $(DOCDIR) $(SCRIPTSWITHMANPAGES)

clean :
	@ echo ++++++++++++++ nothing to clean up here... ++++++++++++++++++

