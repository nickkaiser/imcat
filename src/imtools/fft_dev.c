#define	usage "\n\n\n\
NAME\n\
	fft --- take the fast fourier transform of a fits image\n\
\n\
SYNOPSIS\n\
	fft [options...]\n\
\n\
DESCRIPTION\n\
	'fft' creates a three dimensional image with dimensions\n\
	N1 = Nx, N2 = Ny, N3 = 2 whose zeroth and first planes\n\
	contain the real/imaginary parts of the discrete fft of the\n\
	input image. Zero spatial fequency lives at Nx/2, Ny/2.\n\
	Options are:\n\
		-s		# output 16 bit image\n\
		-i		# output 32 bit int format image\n\
		-n		# no MAGIC substitution\n\
		-m magicval	# replacement for MAGIC pixels (0)\n\
		-c		# don't cycle input image.\n\
		-I		# perform inverse transform.\n\
		-C		# complex fft   \n\
	Default output image format is BITPIX = -32 (i.e. 4-byte real).\n\
\n\
	With the numerical-recipes fft package this only works for\n\
	image dimensions 2^N.\n\
\n\
	This may be useful for `psf-surgery':\n\
	You could 'fft' a composite psf image; modify the fft e.g.\n\
	by dividing it into some desired circular psf, and then apply\n\
	the fft as a filter to your original image using 'smooth -F'\n\
	and --- hey presto! --- your psf is now circular.\n\
\n\
	By default, and when doing a forward transform\n\
	the input image will be cycled by N1/2, N2/2, but\n\
	you can override this with -c option.\n\
\n\
	There is redundancy in the fft image created from real input since\n\
	f(-k) = f^*(-k).\n\
\n\
	Use -C option to take tranform of a complex function.\n\
\n\
AUTHOR\n\
	Nick Kaiser:  kaiser@cita.utoronto.ca\n\
\n\n\n"		

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "../imlib/fits.h"
#include "../utils/error.h"
#include "../utils/arrays.h"
#include "../fftlib/myfft.h"

main(int argc, char *argv[])	
{
	float 		**fin, **fout, magicval;
	int		arg = 1, x, y, N1, N2, iplane, nplanes, fixmagic, docycle, doinverse, pixtype, complexfft;
	fitsheader	*fits;
	fft_type	fk, fk2;

	
	/* defaults */
	magicval = 0.0;
	fixmagic = 1;
	docycle = 1;
	doinverse = 0;
	pixtype = FLOAT_PIXTYPE;
	complexfft = 0;
	nplanes = 1;
	
	while (arg < argc) {
		switch (argv[arg++][1]) {
			case 's':
				pixtype = SHORT_PIXTYPE;
				break;
			case 'i':
				pixtype = INT_PIXTYPE;
				break;
			case 'n':
				fixmagic = 0;
				break;
			case 'm':
				sscanf(argv[arg++], "%f", &magicval);
				break;
			case 'c':
				docycle = 0;
				break;
			case 'I':
				doinverse = 1;
				break;
			case 'C':
				complexfft = 1;
				break;
			default:
				error_exit(usage);
				break;
		}
	}


        fits = readfitsheader(stdin);
        N1 = fits->n[0];
        N2 = fits->n[1];
	if (doinverse || complexfft) {
		allocFloatArray(&fin, N1, 2 * N2);
		for (y = 0; y < 2 * N2; y++) {
			readfitsline(fin[y], fits);
		}
		if (!complexfft) {
			fits->ndim = 2;
		}
	} else {
		allocFloatArray(&fin, N1, N2);
		allocFloatArray(&fout, N1, 2 * N2);
		if (fits->ndim = 3) {
			/* mutiplane input image */
			fits->ndim = 4;
			nplanes = fits->n[3] = fits->n[2];
		} else {
			fits->ndim = 3;
		}
		fits->n[2] = 2;
	}
	fits->extpixtype = pixtype;
	add_comment(argc, argv, fits);

	if (!doinverse) {
	    writefitsheader(fits);
	    for (iplane = 0; iplane < nplanes; iplane++) {
		readfitsplane((void **)fin, fits);
		if (docycle) {
			if (complexfft) {
				cycleimage(fin, N1, N2, N1 / 2, N2 / 2);
				cycleimage(fin + N2, N1, N2, N1 / 2, N2 / 2);
			} else {
				cycleimage(fin, N1, N2, N1 / 2, N2 / 2);
			}
		}
		if (fixmagic) {
			if (complexfft) {
				substitute(fin, N1, 2 * N2, magicval);
			} else {
				substitute(fin, N1, N2, magicval);
			}
		}
		if (complexfft) {
			forward_cfft(fin, N1, N2, fout);
		} else {
			alloc_fft(&fk, N1, N2);
			forward_fft(fin, N1, N2, fk);
			freeFloatArray(fin, N1, N2);
			get_fft(fk, N1, N2, fout);
		}
		for (y = 0; y < 2 * N2; y++) {
			writefitsline(fout[y], fits);
		}
	    }
	    writefitstail(fits);
	} else {
		if (complexfft) {
			allocFloatArray(&fout, N1, 2 * N2);
			inverse_cfft(fin, N1, N2, fout);
		} else {
			alloc_fft(&fk, N1, N2);
			allocFloatArray(&fout, N1, N2);
			set_fft(fk, N1, N2, fin);
			freeFloatArray(fin, N1, N2 * 2);
			inverse_fft(fk, N1, N2, fout);
		}
		if (docycle) {
			if (complexfft) {
				cycleimage(fout, N1, N2, -N1 / 2, -N2 / 2);
				cycleimage(fout + N2, N1, N2, -N1 / 2, -N2 / 2);
			} else {
				cycleimage(fout, N1, N2, -N1 / 2, -N2 / 2);
			}
		}
		if (complexfft) {
			writefitsheader(fits);
			for (y = 0; y < 2 * N2; y++) {
				writefitsline(fout[y], fits);
			}
			writefitstail(fits);
		} else {
			write2Dfloatimage(fout, fits);
		}
	}

	exit(0);
}











