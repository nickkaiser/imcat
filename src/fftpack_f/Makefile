# Makefile for the fftpack library

OBJECTS= cfftb.o  cosqb1.o  ezfftb.o  passb5.o  radb3.o   radf5.o   rffti1.o \
cfftb1.o cosqf.o   ezfftf.o  passf.o   radb4.o   radfg.o   sinqb.o \
cfftf.o  cosqf1.o  ezffti.o  passf2.o  radb5.o   rfftb.o   sinqf.o \
cfftf1.o cosqi.o   passb.o   passf3.o  radbg.o   rfftb1.o  sinqi.o \
cffti.o  cost.o    passb2.o  passf4.o  radf2.o   rfftf.o   sint.o \
cffti1.o costi.o   passb3.o  passf5.o  radf3.o   rfftf1.o  sint1.o \
cosqb.o   ezfft1.o  passb4.o  radb2.o   radf4.o   rffti.o   sinti.o

SHELL = /bin/sh
FFLAGS = -fast

all : libfftpack_f.a

libfftpack_f.a : $(OBJECTS)
	ar r libfftpack_f.a $?

ranlib :
	ranlib libfftpack_f.a

install : libfftpack_f.a $(RANLIB)
	@echo +++++++++++++ installing libfftpack_f.a in ../../lib/$(ARCH)/ +++++++++++++
	@ if test -d ../../lib/$(ARCH) ; then \
		echo  +++++++++++++ dest-dir exists +++++++++++++; else mkdir ../../lib/$(ARCH) ; fi
	cp libfftpack_f.a ../../lib/$(ARCH)

docs :

clean:
	rm -f *.o libfftpack_f.a


