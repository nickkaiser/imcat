# Makefile for ic

all: ic

BINDIR = ../../../bin/$(ARCH)
LIB = -L../../../lib/$(ARCH) -lim -lutils -lm
SHELL = /bin/sh

OBJECTS = ic.o getop.o operators.o op_math.o op_logic.o

install : ic
	@echo installing imtools in $(BINDIR) ....
	@ if test -d $(BINDIR) ; then \
		echo $(BINDIR) exists ; else mkdir $(BINDIR) ; fi
	mv ic $(BINDIR)

ic : $(OBJECTS)
	$(CC) $(OBJECTS) $(LIB) -o ic

clean:
	rm -f core *.o
