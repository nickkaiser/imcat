# Makefile to test NR amoeba programme

LIBS = -L../../../lib/$(ARCH) -lcat -lim -lfft $(FFTLIB) -lutils -lm
NROBJS = powell_dbl.o linmin_dbl.o mnbrak_dbl.o brent_dbl.o f1dim_dbl.o

all : fitstack testpowell

fitstack: fitstack.o fitstack_read.o fitstack_modefunc.o $(NROBJS)
	$(CC) fitstack.o fitstack_read.o fitstack_modefunc.o $(NROBJS) $(LIBS) -o fitstack


testpowell: testpowell.o $(NROBJS)
	$(CC) testpowell.o $(NROBJS) $(LIBS) -o testpowell
