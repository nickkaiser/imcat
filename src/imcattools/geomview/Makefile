# Makefile for geomview tools

TOOLS = fits3Dviewer
SCRIPTS = geomviewsetup

LIBS = -L$(IMCATDIR)/lib/$(ARCH) -lim -lutils -lm
INCLUDES = -I$(IMCATDIR)/src

BINDIR = $(IMCATDIR)/bin/$(ARCH)
SCRIPTDIR = $(IMCATDIR)/bin/scripts
DOCDIR = $(IMCATDIR)/doc/imcattools/geomview

# XFORMSLIBS = -L$(GEOMVIEWLIBDIR)/xforms/FORMS -lforms
# GEOMVIEWINCLUDES = -I$(GEOMVIEWLIBDIR)/xforms/FORMS -I$(GEOMVIEWLIBDIR)/include
XFORMSLIBS = -L$(GEOMVIEWLIBDIR)/lib -lforms
GEOMVIEWINCLUDES = -I$(GEOMVIEWLIBDIR)/xforms/FORMS -I$(GEOMVIEWLIBDIR)/include -I$(GEOMVIEWLIBDIR)/include/geomview

OBJECTS = main.o control_panel.o callbacks.o mesh.o

SHELL = /bin/sh

all : $(TOOLS)

install : $(TOOLS)
	@echo installing geomview tools in $(BINDIR) ....
	for i in $(SCRIPTS); do \
		if test -h $(SCRIPTDIR)/$$i ; then rm -f $(SCRIPTDIR)/$$i ; fi ; \
		ln -s ../../src/imcattools/geomview/$$i $(SCRIPTDIR)/$$i ; \
	done
#	@ if test -d $(BINDIR) ; then \
#		echo $(BINDIR) exists ; else mkdir $(BINDIR) ; fi
	for i in $(TOOLS); do \
		cp $$i $(BINDIR)/$$i ; \
	done

clean:
	rm -f core *.o
	for i in $(TOOLS); do rm -f $$i ; done

docs :
	@echo making Man Pages for imcattools_geomview
	@if test ! -d $(DOCDIR) ; then mkdir -p $(DOCDIR) ; fi
	makemanpages imcattools_geomview $(MANDIR) $(DOCDIR) $(TOOLS) $(SCRIPTS)


##
## the tools:
##

fits3Dviewer : $(OBJECTS)
	$(CC) $(OBJECTS) $(IMCATLIBS) $(XFORMSLIBS) $(XLIBS) $(LIBS) -o fits3Dviewer

.c.o:
	cc $(GEOMVIEWINCLUDES) $(INCLUDES) -c $<

