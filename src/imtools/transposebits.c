#define usage "\n\n\n\
NAME\n\
	transposebits - swap order of bits and numbers in a FITS file\n\
\n\
SYNOPSIS\n\
	transposebits [-i]\n\
\n\
DESCRIPTION\n\
	\"transposebits\" reads a FITS file a line at a time and transposes the\n\
	order of bits and numbers such that if the input bit\n\
	array is b[x][i], where 0 <= x < N1, and 0 <= i < |BITPIX|\n\
	then the output line is the bit array b[i][x]. Thus the\n\
	N1 lowest order bits will be output first, followed by the N1\n\
	second lowest and so on.  This may aid compression.\n\
\n\
	Use -i option to perform the inverse transpose.\n\
\n\
AUTHOR\n\
	Nick Kaiser:  kaiser@hawaii.edu\n\
\n\n\n"


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>


#include "../imlib/fits.h"
#include "../utils/error.h"

int		main(int argc, char *argv[])	
{
	int		N, N1, N2, b, nlines, dim, doinverse, arg = 1;
	fitsheader	*fits;
	void		*fin, *fout;
	unsigned short	*fsin, *fsout;
	unsigned int	*fiin, *fiout;
	int		i, x, it, xt, iin, xin, iout, xout;
	FILE		*theipf, *theopf;

	/* defaults */
	theipf = stdin;
	theopf = stdout;
	doinverse = 0;

	/* parse args */
	if (argc > 2) {
		error_exit(usage);
	}
	while (arg < argc) {
		if (argv[arg][0] != '-') {
			error_exit(usage);
		}
		switch (argv[arg++][1]) {
			case 'i':
				doinverse = 1;
				break;
			default:
				error_exit(usage);
		}
	}

	/* read the fits header, set opstream, and count total number of lines */
	fits = readfitsheader(theipf);
	fits->opstream = theopf;
	N1 = fits->n[0];
	nlines = 1;
	for (dim = 1; dim < fits->ndim; dim++) {
		nlines *= fits->n[dim];
	}
	add_comment(argc, argv, fits);
	writefitsheader(fits);

	/* allocate input and output arrays */
	N = (fits->extpixtype > 0 ? fits->extpixtype : -(fits->extpixtype));
	fin  = calloc(N1, sizeof(N / 8));
	fout = calloc(N1, sizeof(N / 8));
	fsin  = (unsigned short *) fin;
	fsout = (unsigned short *) fout;
	fiin  = (unsigned int *) fin;
	fiout = (unsigned int *) fout;

	while (nlines--) {
		fread(fin, N / 8, N1, fits->ipstream);
		/* zero the output line */
		memset(fout, 0, N1 * N / 8);
		for (x = 0; x < N1; x++) {
			for (i = 0; i < N; i++) {
				b = i * N1 + x;
				xt = b / N;
				it = b - N * xt;
				if (doinverse) {
					xin = xt;
					iin = it;
					xout = x;
					iout = i;
				} else {
					xin = x;
					iin = i;
					xout = xt;
					iout = it;
				}
				switch (N) {
					case 16:
						fsout[xout] = fsout[xout] | (((fsin[xin] >> iin) & (unsigned short) 1) << iout);
						break;
					case 32:
						fiout[xout] = fiout[xout] | (((fiin[xin] >> iin) & (unsigned int) 1) << iout);
						break;
					default:
						error_exit("transposebits: I can only handle 16 or 32 bit images\n");
				}
			}
		}
		fwrite(fout, N / 8, N1, fits->opstream);
	}

	writefitstail(fits);
	exit(0);
}



