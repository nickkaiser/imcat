# Makefile to make lc

TOOLS=lc makerandcat makegridcat
all : $(TOOLS)

CFLAGS=-g
INCLUDE=../../catlib/cat.h stack.h getop.h error.h operators.h rpn.h lc.h
OP_INCLUDE=op_math.h op_logic.h op_deref.h op_vector.h op_dot.h op_rand.h op_misc.h \
op_lintrans.h op_strings.h
OP_OBJECTS=op_math.o op_logic.o op_deref.o op_vector.o op_dot.o op_rand.o op_misc.o \
op_lintrans.o op_strings.o
SHELL=/bin/sh
BINDIR = ../../../bin/$(ARCH)
LIBS = -L../../../lib/$(ARCH) -lcat -lim -lutils -lm
DOCDIR = ../../../doc/cattools/lc

clean :
	rm -f *.o core testrpn test $(TOOLS)

install : $(TOOLS)
	@echo installing imtools in $(BINDIR) ....
	@ if test -d $(BINDIR) ; then \
		echo $(BINDIR) exists ; else mkdir $(BINDIR) ; fi
	for i in $(TOOLS); do cp $$i $(BINDIR)/$$i ; done


docs :
	@echo making html Man Pages for cattools_lc
	@if test ! -d $(DOCDIR) ; then mkdir $(DOCDIR) ; fi
	makemanpages  cattools_lc $(MANDIR) $(DOCDIR) $(TOOLS) $(SCRIPTS)

lc : lc.o rpn.o stack.o getop.o operators.o $(OP_OBJECTS)
	$(CC) lc.o rpn.o stack.o getop.o operators.o $(OP_OBJECTS) $(LIBS) -o lc

lc.o : lc.c $(INCLUDE) lc.man
	$(CC) -c lc.c -o lc.o

test : test.o
	$(CC) test.o $(LIBS) -o test

test.o : test.c $(INCLUDE)

cat.o : cat.c cat.h

testrpn : testrpn.o rpn.o cat.o stack.o getop.o operators.o
	$(CC) testrpn.o rpn.o cat.o stack.o getop.o operators.o $(LIBS) -o testrpn

testrpn.o : testrpn.c $(INCLUDE) 

rpn.o : rpn.c $(INCLUDE) 

stack.o : stack.c $(INCLUDE) 

getop.o : getop.c $(INCLUDE) 

operators.o : operators.c $(INCLUDE) $(OP_INCLUDE)

maketestcat : maketestcat.c
	$(CC) maketestcat.c -lm -o maketestcat

makerandcat : makerandcat.c
	$(CC) makerandcat.c $(LIBS) -o makerandcat

makegridcat : makegridcat.c
	$(CC) makegridcat.c $(LIBS) -o makegridcat
