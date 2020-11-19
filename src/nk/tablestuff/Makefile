# Makefile to make tablestuff directory

TOOLS = tablearith tablefilter tablesort mova tablecumulate tableop tablesum

all : $(TOOLS)

SHELL=/bin/sh
BINDIR = ../../../bin/$(ARCH)
DOCDIR = ../../../doc/nk/tablestuff

# suffix rule (adds the -lm library invocation)
.c :
	$(CC) $(CFLAGS) $(LDFLAGS) $<  -o $@ -lm

clean :
	rm -f *.o core
	for i in $(TOOLS); do rm -f $$i ; done

install : $(TOOLS) 
	@echo installing tablestuff in $(BINDIR) ....
	@ if test -d $(BINDIR) ; then \
		echo $(BINDIR) exists ; else mkdir $(BINDIR) ; fi 
	for i in $(TOOLS) ; do \
		cp $$i $(BINDIR)/$$i ; \
	done

docs :
	@echo making html Man Pages for nk_tablestuff
	@if test ! -d $(DOCDIR) ; then mkdir -p $(DOCDIR) ; fi
	makemanpages nk_tablestuff $(MANDIR) $(DOCDIR) $(TOOLS) $(SCRIPTS)
