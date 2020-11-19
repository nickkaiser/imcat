#define	usage "\n\n\n\
NAME\n\
	modelpsf --- generate a 2D polynomial model for psf\n\
\n\
SYNOPSIS\n\
	modelpsf starcat fitsimage [options...]\n\
		-n N		# postage stamp image size (32)\n\
		-l lmax		# max order for fit (1)\n\
		-g		# generate psf model\n\
		-s		# read stamps image\n\
		-w f_bg		# use background in weighting\n\
\n\
DESCRIPTION\n\
	'modelpsf' first reads a catalogue of stars from 'starcat'\n\
	and fits for the following model for the psf:\n\
\n\
		fpsf(x; xobj) = sum_l fpsf_l(x) f_l(xobj)\n\
\n\
	where f_l(xobj) are polynomial mode functions and the fpsf_l(x)\n\
	are the image valued 'mode amplitudes'.\n\
\n\
	The stars need to be normalised in flux, so it is necessary\n\
	that the input catalogues contain an entry\n\
	'flux' in addition to the position vector 'x[2]'.\n\
\n\
	By default it generates and reads a postage stamp album image\n\
	of these stars from 'fitsimage'. We then\n\
	do a least squares fit to obtain a final image model psf\n\
	amplitudes fpsf_l(x).\n\
\n\
	With -s option 'fitsimage' is supplied as a ready made stamp\n\
	album.\n\
\n\
	By default, objects recieve weight proportional to their\n\
	flux, which is optimum in the limit of negligible sky background.\n\
	Use -w option to supply a background value (input image is assumed to\n\
	have had background subtracted however), and then weight\n\
	objects in proportion to flux^2 / (f + f_bg).\n\
\n\
	With the -g option, the fits image is instead interpreted\n\
	as the nmodes x N x N model psf image and we generate\n\
	an nstars x N x N image containing the synthesised\n\
	models.  In this mode the -n flag is ignored and\n\
	the input catalogue need contain only a position vector\n\
	x[2].\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@hawaii.edu\n\
\n\n\n"		


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "../../imlib/fits.h"
#include "../../utils/error.h"
#include "../../utils/args.h"
#include "../../utils/ipbuff.h"
#include "../../utils/modefunc.h"
#include "../../utils/lu.h"

#define HEADER_SIZE	1024

main(int argc, char *argv[])	
{
	int		arg = 3;
	char		*starcatname, *fitsimagename, string[128], header[HEADER_SIZE], *flag;
	int		N, nstars, istar, x, y;
	int		l, m, mode0, mode1, lmin, lmax, nmodes, *ll, *mm, *indx;
	fitsheader	*fits;
	FILE		*fitsf, *catf, *lcpipe;
	float		***fpsf, ***fmodel;
	double		**ipbuff, *X, flux, fbg, weight;
	double		**A, *B, fmode0, fmode1, det;
	int		generatemodel, makestamps, dobgweighting;


	/* defaults */
	N = 32;
	lmin = 0;
	lmax = 1;
	generatemodel = 0;
	makestamps = 1;
	dobgweighting = 0;
	
	/* parse args */
	argsinit(argc, argv, usage);
	if (argc < 3) {
		error_exit(usage);
	}
	starcatname 	= getargs();
	fitsimagename 	= getargs();
	while (flag = getflag()) {
		switch (flag[0]) {
			case 'n':
				N = getargi();
				break;
			case 'l':
				lmax = getargi();
				break;
			case 'g':
				generatemodel = 1;
				break;
			case 's':
				makestamps = 0;
				break;
			case 'w':
				fbg = getargd();
				dobgweighting = 1;
				break;	
			default:
				error_exit(usage);
				break;
			}
	}

	/* read the catalogue */
	if (generatemodel) {
		sprintf(string, "lc -b -o x < %s", starcatname);
	} else {
		sprintf(string, "lc -b -o x flux < %s", starcatname);
	}
	catf = popen(string, "r");
	if (!catf) {
		error_exit("modelpsf: unable to open star catalogue\n");
	}
	if (generatemodel) {
		ipbuff = readdoublebuff(2, catf, &nstars);
	} else {
		ipbuff = readdoublebuff(3, catf, &nstars);
	}
	pclose(catf);

	if (generatemodel) {
		fitsf = fopen(fitsimagename, "r");
		if (!fitsf) {
			error_exit("modelpsf: unable to open fits image for input\n");
		}
		fits = readfitsheader(fitsf);
		if (fits->ndim != 3) {
			error_exit("modelpsf: fits image must be 3-dimensional psf model\n");
		}
		nmodes = fits->n[2];
		N = fits->n[0];
	} else {
		/* compute number of modes */
		nmodes = 0;
		for (l = 0; l <= lmax; l++) {
			nmodes += (l + 1);
		}
	}

	/* set up mode indices */
	ll = (int *) calloc(nmodes, sizeof(int));
	mm = (int *) calloc(nmodes, sizeof(int));
	mode0 = 0;
	for (l = lmin; l <= lmax; l++) {
		for (m = 0; m <= l; m++) {
			ll[mode0] = l;
			mm[mode0++] = m;
		}			
	}

	/* allocate space for the model amplitudes */
	fmodel = (float ***) calloc(nmodes, sizeof(float **));
	for (mode0 = 0; mode0 < nmodes; mode0++) {
		allocFloatArray(&(fmodel[mode0]), N, N);
	}

	/* allocate space for fpsf */
	fpsf = (float ***) calloc(nstars, sizeof(float **));
	for (istar = 0; istar < nstars; istar++) {
		allocFloatArray(&(fpsf[istar]), N, N);
	}

	if (generatemodel) {
		/* read the model amplitude images */
		for (mode0 = 0; mode0 < nmodes; mode0++) {
			for (y = 0; y < N; y++) {
				readfitsline(fmodel[mode0][y], fits);
			}
		}
		/* compute the models */
		for (istar = 0; istar < nstars; istar++) {
			X = ipbuff[istar];
			for (mode0 = 0; mode0 < nmodes; mode0++) {
				fmode0 = f(ll[mode0], mm[mode0], X);
				for (y = 0; y < N; y++) {
					for (x = 0; x < N; x++) {
						fpsf[istar][y][x] += fmode0 * fmodel[mode0][y][x];
					}
				}
			}
		}
		fits->n[2] = nstars;
		N = fits->n[0];
		add_comment(argc, argv, fits);
		writefitsheader(fits);
		for (istar = 0; istar < nstars; istar++) {
			for (y = 0; y < N; y++) {
				writefitsline(fpsf[istar][y], fits);
			}
		}
		writefitstail(fits);
		exit(0);
	}

	/* allocate space for linear algebra */
	B = (double *) calloc(nmodes, sizeof(double));
	A = (double **) calloc(nmodes, sizeof(double *));
	for (mode0 = 0; mode0 < nmodes; mode0++) {
                 A[mode0] = (double *) calloc(nmodes, sizeof(double));
	}
        indx = (int *) calloc(nmodes, sizeof(int));

	/* generate, if necessary, and read the stamp album */
	if (makestamps) {
		sprintf(string, "makestamps -N flux -b %d -f %s < %s", N, fitsimagename, starcatname);
		fitsf = popen(string, "r");
	} else {
		fitsf = fopen(fitsimagename, "r");
	}
	if (!fitsf) {
		error_exit("modelpsf: unable to open fits image\n");
	}	

	/* read the stamp album */
	fits = readfitsheader(fitsf);
	if (nstars != fits->n[2]) {
		error_exit("modelpsf: mismatch between nstars and image dimensions - bailing\n");
	}
	for (istar = 0; istar < nstars; istar++) {
		for (y = 0; y < N; y++) {
			readfitsline(fpsf[istar][y], fits);
		}
	}
	pclose(fitsf);

	/* now we loop over pixels */
	for (y = 0; y < N; y++) {
		for (x = 0; x < N; x++) {
			/* zero arrays */
			for (mode0 = 0; mode0 < nmodes; mode0++) {
				B[mode0] = 0.0;
				for (mode1 = 0; mode1 < nmodes; mode1++) {
					A[mode0][mode1] = 0.0;
				}
			}
			/* accumulate arrays */
			for (istar = 0; istar < nstars; istar++) {
				X = ipbuff[istar];
				flux = ipbuff[istar][2];
				if (dobgweighting) {
					weight = flux / (fabs(fpsf[istar][y][x]) + fbg / flux);
				} else {
					weight = flux;
				}
				if (fpsf[istar][y][x] != FLOAT_MAGIC) {
					for (mode0 = 0; mode0 < nmodes; mode0++) {
						fmode0 = f(ll[mode0], mm[mode0], X);
						B[mode0] += fmode0 * fpsf[istar][y][x] * weight;
						for (mode1 = 0; mode1 < nmodes; mode1++) {
							fmode1 = f(ll[mode1], mm[mode1], X);
							A[mode0][mode1] += fmode0 * fmode1 * weight;
						}
					}
				}
			}
			/* solve the linear equations */
			myludcmp(A, nmodes, indx, &det);
			mylubksb(A, nmodes, indx, B);
			/* extract the mode amplitudes */
			for (mode0 = 0; mode0 < nmodes; mode0++) {
				fmodel[mode0][y][x] = B[mode0];
			}
		}
	}

	/* now we output the mode amplitude images */
	fits->n[2] = nmodes;
	add_comment(argc, argv, fits);
	writefitsheader(fits);
	for (mode0 = 0; mode0 < nmodes; mode0++) {
		for (y = 0; y < N; y++) {
			writefitsline(fmodel[mode0][y], fits);
		}
	}
	writefitstail(fits);

	/* all done */
	exit(0);
}
