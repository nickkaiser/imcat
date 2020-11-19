#define	usage "\n\n\n\
NAME\n\
	acf --- calculates 2-D autocorrelation function of a fits image\n\
\n\
SYNOPSIS\n\
	acf [options...]\n\
\n\
DESCRIPTION\n\
	'acf' reads a fits file from stdin and writes the autocorrelation\n\
	function to stdout. If the input image is fin(r) then we compute\n\
		fout(r) = sum_r' fin(r') fin(r' + r) / (N1 * N2)\n\
	and resulting image is wrapped so that zero lag is at pixel (N1/2, N2/2).\n\
\n\
	Options are:\n\
		-p		# calculate power spectrum instead\n\
		-c f1 f2	# cross correlate named fits files\n\
		-P pixtype	# output pixtype (FLOAT_PIXTYPE)\n\
		-n		# no MAGIC substitution\n\
		-m magicval	# replacement for MAGIC pixels (0)   \n\
\n\
	Power is defined so that white noise with variance sigma2\n\
	will produce P = sigma2. Power is translated so that zero frequency\n\
	lies at N1/2, N2/2.\n\
\n\
	With '-c fits0 fits1' we compute\n\
\n\
		c(r) = sum_r' f0(r') f1(r' + r) / (N1 * N2)\n\
\n\
	to compute instead\n\
\n\
		c(r) = sum_r' f0(r') f1(r - r') / (N1 * N2)\n\
\n\
	you should first rotate the first image by 180 degrees.\n\
\n\
	Supply file name '-' to read from stdin.\n\
\n\
	Default output image format is BITPIX = -32 (i.e. 4-byte float).\n\
\n\
	If source image is 3D fits then result is a 3D image where each\n\
	plane is acf, or power spectrum, of corresponding input frame.  Does\n\
	not work in -c mode.\n\
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
#include "../utils/args.h"
#include "../utils/arrays.h"
#include "../fftlib/myfft.h"

main(int argc, char *argv[])	
{
	float 		**theresult, **f1, **f2, mul, magicval;
	int		arg = 1, nplanes = 1, iplane, N1, N2, M1, M2, dopower, doccf, i, j, fixmagic, oppixtype;
	fitsheader	*fitsh1, *fitsh2, *fitsho;
	fft_type	fk1, fk2;
	FILE		*fits1, *fits2;
	char		*flag, *fitsname1, *fitsname2;
	
	/* defaults */
	dopower = 0;
	doccf = 0;
	magicval = 0.0;
	fixmagic = 1;
	oppixtype = FLOAT_PIXTYPE;

	/* parse arguments */
	argsinit(argc, argv, usage);
	while (flag = getflag()) {
		switch (flag[0]) {
			case 'p':
				if (doccf)
					error_exit("acf: error cannot do ccf and P(k)\n");
				dopower = 1;
				break;
			case 'c':
				if (dopower)
					error_exit("acf: error cannot do ccf and P(k)\n");
				doccf = 1;
				fitsname1 = getargs();
				fitsname2 = getargs();
				if (strcmp(fitsname1, "-")) {
					fits1 = fopen(fitsname1, "r");
					if (!fits1) {
						error_exit("acf: failed to open fits file 1\n");
					}
				} else {
					fits1 = stdin;
				}
				if (strcmp(fitsname2, "-")) {
					fits2 = fopen(fitsname2, "r");
					if (!fits2) {
						error_exit("acf: failed to open fits file 2\n");
					}
				} else {
					fits2 = stdin;
				}
				break;
			case 'P':
				oppixtype = getargi();
				break;
			case 'n':
				fixmagic = 0;
				break;
			case 'm':
				magicval = getargf();
				break;
			default:
				error_exit(usage);
				break;
		}
	}


	if (doccf) {
		fitsh1 = readfitsheader(fits1);
		N1 = fitsh1->n[0];
		N2 = fitsh1->n[1];
		fitsh2 = readfitsheader(fits2);
		M1 = fitsh2->n[0];
		M2 = fitsh2->n[1];
		if (M1 != N1 || M2 != N2) {
			error_exit("acf: fits files must be same size\n");
		}
	} else {
		fitsh1 = readfitsheader(stdin);
		N1 = fitsh1->n[0];
		N2 = fitsh1->n[1];
		if (fitsh1->ndim == 3) {
			nplanes = fitsh1->n[2];
		}
	}

	fitsho = copyfitsheader(fitsh1);
	setextpixtype(fitsho, oppixtype);
	add_comment(argc, argv, fitsho);
	writefitsheader(fitsho);

	allocFloatArray(&f1, N1, N2);
	alloc_fft(&fk1, N1, N2);
	if (doccf) {
		allocFloatArray(&f2, N1, N2);
		alloc_fft(&fk2, N1, N2);
	}
	allocFloatArray(&theresult, N1, N2);

	for (iplane = 0; iplane < nplanes; iplane++) {
		readfitsplane((void **) f1, fitsh1);
		if (fixmagic) {
			substitute(f1, N1, N2, magicval);
		}
		forward_fft(f1, N1, N2, fk1);
		if (dopower) {
			power(fk1, N1, N2, theresult, N1/2, N2/2);
		} else {
			if (doccf) {
				readfitsplane((void **) f2, fitsh2);
				if (fixmagic) {
					substitute(f2, N1, N2, magicval);
				}
				forward_fft(f2, N1, N2, fk2);
			} else {
				fk2 = fk1;
			}
			ccf(fk1, fk2, N1, N2, theresult, N1/2, N2/2);
		}
		writefitsplane(theresult, fitsho);
	}
	writefitstail(fitsho);
	exit(0);
}











