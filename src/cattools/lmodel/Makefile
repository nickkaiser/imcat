# Makefile for cattools

TOOLS = fitlmodel generatelmodel generatelmodelimage difflmodel
SCRIPTS =
# all : $(TOOLS) test
all : $(TOOLS)

SHELL=/bin/sh
BINDIR = $(IMCATDIR)/bin/$(ARCH)
SCRIPTDIR = $(IMCATDIR)/bin/scripts
# LIBS = -L$(IMCATDIR)/lib/$(ARCH) -lcat -lutils -lm
LIBS = -L$(IMCATDIR)/lib/$(ARCH) -lutils -lcat -lim -lutils -lm
DOCDIR = $(IMCATDIR)/doc/cattools/lmodel
OBJS = 

clean :
	rm -f *.o core
	for i in $(TOOLS); do rm -f $$i ; done

install : $(TOOLS)
#	@echo installing scripts in $(SCRIPTDIR) ....
#	@ if test ! -d $(SCRIPTDIR) ; then mkdir $(SCRIPTDIR) ; fi
#	for i in $(SCRIPTS); do \
#		if test -h $(SCRIPTDIR)/$$i ; then rm -f $(SCRIPTDIR)/$$i ; fi ; \
#		ln -s ../../src/cattools/lmodel/$$i $(SCRIPTDIR)/$$i ; \
#	done
	@echo installing cattools in $(BINDIR) ....
	@ if test ! -d $(BINDIR) ; then mkdir $(BINDIR) ; fi
	for i in $(TOOLS); do cp $$i $(BINDIR)/$$i ; done

docs :
	@echo making html Man Pages for cattools/lmodel
	@if test ! -d $(DOCDIR) ; then mkdir $(DOCDIR) ; fi
	makemanpages cattools_lmodel  $(MANDIR) $(DOCDIR) $(TOOLS) $(SCRIPTS)

fitlmodel : fitlmodel.o $(OBJS)
	$(CC) fitlmodel.o $(OBJS) $(LIBS) -o fitlmodel

fitlmodel.o : fitlmodel.c
	$(CC) -c -I$(IMCATDIR)/src fitlmodel.c -o fitlmodel.o

generatelmodel : generatelmodel.o $(OBJS)
	$(CC) generatelmodel.o $(OBJS) $(LIBS) -o generatelmodel

generatelmodel.o : generatelmodel.c
	$(CC) -c -I$(IMCATDIR)/src generatelmodel.c -o generatelmodel.o

generatelmodelimage : generatelmodelimage.o $(OBJS)
	$(CC) generatelmodelimage.o $(OBJS) $(LIBS) -o generatelmodelimage

generatelmodelimage.o : generatelmodelimage.c
	$(CC) -c -I$(IMCATDIR)/src generatelmodelimage.c -o generatelmodelimage.o

difflmodel : difflmodel.o $(OBJS)
	$(CC) difflmodel.o $(OBJS) $(LIBS) -o difflmodel

difflmodel.o : difflmodel.c
	$(CC) -c -I$(IMCATDIR)/src difflmodel.c -o difflmodel.o

lmodel.o : lmodel.c lmodel.h
	$(CC) -c -I$(IMCATDIR)/src lmodel.c -o lmodel.o

lmodelio.o : lmodelio.c lmodel.h
	$(CC) -c -I$(IMCATDIR)/src lmodelio.c -o lmodelio.o

lmodelcalculus.o : lmodelcalculus.c lmodel.h
	$(CC) -c -I$(IMCATDIR)/src lmodelcalculus.c -o lmodelcalculus.o

test : test.o $(OBJS)
	$(CC) test.o $(OBJS) $(LIBS) -o test

test.o : test.c lmodel.h
	$(CC) -c -I$(IMCATDIR)/src test.c -o test.o

