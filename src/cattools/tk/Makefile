# makefile for tk directory scripts

SCRIPTS = xyplot eplot
SHELL = /bin/sh

all:
	@ echo +++++++++++++ this directory only contains scripts +++++++++++++

BINDIR = ../../../bin/scripts

install :
	if test ! -d $(BINDIR) ; then mkdir $(BINDIR) ; fi ;
	for i in $(SCRIPTS); do \
		if test -h $(BINDIR)/$$i ; then rm $(BINDIR)/$$i ; fi ; \
		ln -s ../../src/cattools/tk/$$i $(BINDIR)/$$i ; \
	done


clean :
	@ echo ++++++++++++++ nothing to clean up here... ++++++++++++++++++

