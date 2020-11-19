# Makefile for ic

all: ic

BINDIR = ../../../bin/$(ARCH)
LIB = -L../../../lib/$(ARCH) -lim -lutils -lm
SHELL = /bin/sh
DOCDIR = ../../../doc/imtools/ic

OBJECTS = ic.o getop.o

install : ic
	@echo installing imtools in $(BINDIR) ....
	@ if test -d $(BINDIR) ; then \
		echo $(BINDIR) exists ; else mkdir $(BINDIR) ; fi
	cp ic $(BINDIR)

ic : $(OBJECTS)
	$(CC) $(OBJECTS) $(LIB) -o ic

clean:
	rm -f core *.o ic

docs :
	@echo making html Man Pages for ic
	@if test ! -d $(DOCDIR) ; then mkdir $(DOCDIR) ; fi
	makemanpages imtools_ic $(MANDIR) $(DOCDIR) ic
