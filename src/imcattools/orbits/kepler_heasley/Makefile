all : libkepler.a

OBJECTS = kepler.o kepler_jh.o dummy_MAIN.o

LIBS = -lf2c libkepler.a -lm

libkepler.a : $(OBJECTS)
	ar r libkepler.a $(OBJECTS)
	ranlib libkepler.a

kepler_jh.c : kepler_jh.f
	f2c -AC kepler_jh.f

kepler_jh.o : kepler_jh.c
	gcc -c kepler_jh.c -o kepler_jh.o

testkepler : testkepler.o libkepler.a vectors.o
	gcc testkepler.o vectors.o $(LIBS) -o testkepler

clean :
	rm -f *.o

empty :
	rm -f libkepler.a testkepler


