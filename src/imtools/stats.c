/**
 **
 ** calculate simple statistics for an image
 **
 **/

#define	usage "\n\n\n\
NAME\n\
	stats -- calculate simple statistics for an image\n\
\n\
SYNOPSIS\n\
	stats 	[option...]\n\
		-m n	# ignore outer n pixel margin (def = 0)\n\
		-v stat	# just output value for 'stat', which\n\
			can be one of:\n\
			N1 N2 pixtype min max mean mode median\n\
			lquart uquart sigma\n\
		-s stat bzero	# subtract statistic and output image\n\
		-d stat bzero	# divide by statistic and output image\n\
\n\
DESCRIPTION\n\
	\"stats\" reads a 2-D fits file from stdin\n\
	and writes the descriptive statistics listed above to stdout.\n\
	Mode, sigma are computed as described in the 'catstats' man\n\
	page.\n\
\n\
	If stats is given an image of dimensionality 3 then\n\
	it will generate a lc-format output giving the stats\n\
	for the NAXIS3 planes (each of size NAXIS2 x NAXIS1)\n\
	and with an index 'i = 0 ... NAXIS3 - 1' giving the plane number.\n\
	With an image of dimensionality 4 it generates statistics\n\
	for the NAXIS4 x NAXIS3 planes, and the index i becomes\n\
	a 2-vector with i[0] = 0...NAXIS3-1, i[1] = 0...NAXIS4-1,\n\
	and similarly for higher dimensions.\n\
\n\
SEE ALSO\n\
	catstats(1)\n\
\n\
AUTHOR\n\
	Nick Kaiser:  kaiser@cita.utoronto.ca\n\
\n\n\n"		

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "../utils/stats_stuff.h"
#include "../utils/error.h"
#include "../utils/arrays.h"
#include "../imlib/fits.h"

int main(int argc, char *argv[])	
{
	int		i, ii, j, nplanes, arg = 1, N1, N2, margin = 0, statarg = 0, ind[MAX_FITS_DIM];
	float		**f;
	fstatsrec 	srec;
	fitsheader	*fits;
	FILE		*opf;
	char		lcstring[1024];
	int		subtractstat = 0, dividebystat = 0, ix, iy;
	float		bzero, fsub;

	opf = stdout;

	while (arg < argc) {
		if (*argv[arg] != '-') {
			error_exit(usage);
		}
		switch (*(argv[arg++]+1)) {
			case 'm':
				if (1 != sscanf(argv[arg++], "%d", &margin))
					error_exit(usage);
				break;
			case 'u':
				error_exit(usage);
				break;
			case 'v':
				statarg = arg++;
				break;
			case 's':
				subtractstat = 1;
				statarg = arg++;
				if (1 != sscanf(argv[arg++], "%f", &bzero))
					error_exit(usage);
				break;
			case 'd':
				dividebystat = 1;
				statarg = arg++;
				if (1 != sscanf(argv[arg++], "%f", &bzero))
					error_exit(usage);
				break;
			default:
				return (0);
				break;
		}
	}
	
	fits = readfitsheader(stdin);
	N1 = fits->n[0];
	N2 = fits->n[1];
	allocFloatArray(&f, N1, N2);

	nplanes = 1;
	for (i = 2; i < fits->ndim; i++) {
		nplanes *= fits->n[i];
		ind[i] = 0;
	}

	if (subtractstat || dividebystat) {
		fits->bzero = bzero;
		if (!(fits->bscaling)) {
			fits->bscaling = 1;
			fits->bscale = 1;
		}
		writefitsheader(fits);
	} else {
	if (nplanes > 1) {		/* we output vi lc */
		sprintf(lcstring, "lc -C -x -d -a 'history: stats' -H 'N1 = %d' -H 'N2 = %d' -H 'bitpix = %d' -N '1 %d i' ", 
			N1, N2, fits->extpixtype, fits->ndim - 2);
		if (statarg) {
			strcat(lcstring, "-n ");
			strcat(lcstring, argv[statarg]);
			strcat(lcstring, " ");
		} else {
			strcat(lcstring, "-n min -n max -n mean -n mode -n median -n lquart -n uquart -n sigma -n goodpix -n badpix");
		}
		opf = popen(lcstring, "w");
		if (!opf) {
			error_exit("stats: failed to open lc-command pipe for output\n");
		}
	}
	}


	if (statarg) {
		if (!strcmp(argv[statarg], "N1")) {
			fprintf(opf, "%d\n", N1);
			exit(0);
		}
		if (!strcmp(argv[statarg], "N2")) {
			fprintf(opf, "%d\n", N2);
			exit(0);
		}
		if (!strcmp(argv[statarg], "pixtype")) {
			fprintf(opf, "%d\n", fits->extpixtype);
			exit(0);
		}
	}

	for (i = 0; i < nplanes; i++) {
		readfitsplane((void **)f, fits);
		fdo_stats(f, N1, N2, margin, &srec);
		if (subtractstat || dividebystat) {
			fsub = 0.0;
			if (!strcmp(argv[statarg], "min"))
				fsub = srec.fmin;
			if (!strcmp(argv[statarg], "max")) 
				fsub = srec.fmax;
			if (!strcmp(argv[statarg], "mean"))
				fsub = srec.fmean;
			if (!strcmp(argv[statarg], "mode"))
				fsub = srec.fmode;
                        if (!strcmp(argv[statarg], "median"))
				fsub = srec.fmedian;
			if (!strcmp(argv[statarg], "lquart"))
				fsub = srec.flowerquartile;
			if (!strcmp(argv[statarg], "uquart"))
				fsub = srec.fupperquartile;
			for (iy = 0; iy < N2; iy++) {
				for (ix = 0; ix < N1; ix++) {
					if (f[iy][ix] != FLOAT_MAGIC) {
						f[iy][ix] = (subtractstat ? f[iy][ix] - fsub : f[iy][ix] / fsub);
					}
				}
			}
			writefitsplane((void **)f, fits);
		} else {
		if (nplanes > 1) {
			ii = i;
			for (j = 2; j < fits->ndim; j++) {
				ind[j] = ii % fits->n[j];
				ii /= fits->n[j];
				fprintf(opf, "%d ", ind[j]);
			}
		}
		if (statarg) {
			if (!strcmp(argv[statarg], "min")) {
				fprintf(opf, "%g\n", srec.fmin);
			}
			if (!strcmp(argv[statarg], "max")) {
				fprintf(opf, "%g\n", srec.fmax);
			}
			if (!strcmp(argv[statarg], "mean")) {
				fprintf(opf, "%g\n", srec.fmean);
			}
			if (!strcmp(argv[statarg], "mode")) {
				fprintf(opf, "%g\n", srec.fmode);
			}
			if (!strcmp(argv[statarg], "median")) {
				fprintf(opf, "%g\n", srec.fmedian);
			}
			if (!strcmp(argv[statarg], "lquart")) {
				fprintf(opf, "%g\n", srec.flowerquartile);
			}
			if (!strcmp(argv[statarg], "uquart")) {
				fprintf(opf, "%g\n", srec.fupperquartile);
			}
			if (!strcmp(argv[statarg], "sigma")) {
				fprintf(opf, "%g\n", srec.sigma);
			}
		} else {
			if (nplanes == 1) {
				fprintf(opf, "%4d x %4d image\n", N1, N2);
				fprintf(opf, "BITPIX = %4d\n", fits->extpixtype);
				fprintf(opf, "min = %g; max = %g\n", srec.fmin, 
						srec.fmax);
				fprintf(opf, "mean       = %g\n", srec.fmean);
				fprintf(opf, "mode       = %g\n", srec.fmode);
				fprintf(opf, "median     = %g\n", srec.fmedian);
				fprintf(opf, "quartiles  = %g %g\n", srec.flowerquartile,
									srec.fupperquartile);
				fprintf(opf, "sigma      = %g\n",  srec.sigma);
				fprintf(opf, "goodpix    = %ld\n", srec.goodpix);
				fprintf(opf, "samplesize = %ld\n", srec.samplesize);
				fprintf(opf, "badpix     = %ld\n", srec.badpix);
			} else {
				fprintf(opf, "%14.8g %14.8g %14.8g %14.8g %14.8g %14.8g %14.8g %14.8g %d %d\n",
					srec.fmin, srec.fmax, srec.fmean, srec.fmode, srec.fmedian,
					srec.flowerquartile, srec.fupperquartile, srec.sigma,
					srec.goodpix, srec.badpix);
			}
		}
		}
	}
	if (nplanes > 1) {
		if (subtractstat || dividebystat) {
			writefitstail(fits);
		} else {
			pclose(opf);
		}
	}
	exit (0);
}





