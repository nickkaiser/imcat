/*
 * fitgeometry.c
 */


#define usage "\n\
NAME\n\
	fitgeometry --- fit layout of a set of images\n\
\n\
SYNOPSIS\n\
	fitgeometry \n\
		-c nimages	# number of images (8)\n\
		-l lmax		# max order for distortion polynomials (1)\n\
		-d parfilebase	# basename for transformation parameter files\n\
		-w weightfac	# relative weight for pairs involving reference cat (1)\n\
\n\
DESCRIPTION\n\
	'fitgeometry' reads from stdin the result of merging a set of\n\
	overlapping catalogues, which must contain at least entries for\n\
	a pair of spatial coordinate vectors 'x[2][2]' and image numbers c[2].\n\
	It then fits a model in which the coordinates of the object on\n\
	the zeroth image (the 'reference image') are related to those on\n\
	the c'th image by\n\
			x_0 = x_c + sum_m a_cm f_m(x_c)\n\
	where mode function f_m are polynomials up to order lmax in x.\n\
\n\
	The solution is obtained by minimising squared residuals in x_0 space.\n\
\n\
	'fitgeometry' can be used in various ways.  One application is to\n\
	generate accurately registered (but generally somewhat distorted) images from\n\
	a set of dithered images from a mosaic camera. To do this one must\n\
	first generate a set of overlapping catalogues, one for each chip,\n\
	by 'growing' the coordinate system for some reference exposure. This is done\n\
	by finding a low order polynomial transformation which maps successive\n\
	exposures onto the reference exposure.  Once this is done the overlapping\n\
	catalogues can be merged in pairs and the concatenation of all the merges\n\
	fed to 'fitgeometry' which will then find a solution for the layout of\n\
	the chips on some idealised 'detector plane' (whose coordinates coincide\n\
	with pixel coordinates on one of the chips --- the zeroth chip ---  which may\n\
	be chosen arbitrarily).  This does not take out telescope distortion, and\n\
	(especially with small offsets between the dithered exposures) the process is\n\
	liable to introduce additional distortion.\n\
\n\
	A second application is to register a set of CCD images to some external\n\
	'reference' catalogue --- such as a catalogue derived from the digital\n\
	sky survey or from the USNOA catalogue.  In this case it may be useful to\n\
	use the '-w' option with a small argument to downweight the contribution\n\
	to the 'chi-squared' from pairs with one element in the reference\n\
	catalogue to reflect their relatively poor precision.\n\
\n\
	By default, the transformation parameters a_{lm} are written to stdout\n\
	as a concatenation of '.par' format files (which is not particularly\n\
	useful), but with the -d option you can specify a basename\n\
	 and fitgeometry will create a set of files\n\
	'parfilebase'c.par for c = 0, nimages - 1, where the first of these\n\
	contains mode coefficients which are all zero. You will likely\n\
	want to make 'parfilebase' a directory, in which case make\n\
	sure you terminate it with '/' and remember to 'mkdir' it first.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@ifa.hawaii.edu\n\
\n"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../../utils/modefunc.h"
#include "../../utils/lu.h"

#define SCALE 1.0

static int	nm;
int	I(int c, int m);

main(int argc, char *argv[])
{
	/* number of pairs, exposures, etc */
	static	int nM, lmax, *ll, *mm, *indx;
	/* model parameters */
	static	double	*a[2];
	/* matrix stuff */
	static	double	**A[2], *B[2], det;
	/* parameter file */
	FILE	*lcpipe, *tempf;
	/* arg counter */
	int	arg = 1;
	/* miscellaneous indices */
	int	i, j, l, m, n, M, c[2], nchips, chip;
	/* input record */
	double *inputrecord, *x[2], *cd;
	/* output stuff */
	char	*parfilebase, *parfilename;
	int	makeparfiles;
	char	*vardef[1];
	/* weight and weight factor */
	double	w, weightfac;

	/* defaults */
	lmax = 1;
	nchips = 8;
	makeparfiles = 0;
	parfilename = NULL;
	vardef[0] = "1 2 a";
	weightfac = 1.0;

	/* parse args */
	while (arg < argc) {
		if (argv[arg][0] != '-') {
			error_exit(usage);
		}
		switch (argv[arg++][1]) {
			case 'l':
				if (argc - arg < 1) {
					error_exit(usage);
				}
				if ((1 != sscanf(argv[arg++], "%d", &lmax))) {
						error_exit(usage);
				}
				break;
			case 'c':
				if (argc - arg < 1) {
					error_exit(usage);
				}
				if ((1 != sscanf(argv[arg++], "%d", &nchips))) {
						error_exit(usage);
				}
				break;
			case 'd':
				if (argc - arg < 1) {
					error_exit(usage);
				}
				parfilebase = argv[arg++];
				makeparfiles = 1;
				break;
			case 'w':
				if (argc - arg < 1) {
					error_exit(usage);
				}
				if ((1 != sscanf(argv[arg++], "%lf", &weightfac))) {
						error_exit(usage);
				}
				break;
			case 'u':
			default:
				error_exit(usage);
		}
	}

	/* copy args to modefunc argstring */
	modefunc_addargcomment(argc, argv);

	/* compute number of coefficients per chip*/
	nm = 0;
	for (l = 0; l <= lmax; l++) {
		nm += (l + 1);
	}
	/* and total number of modes */
	nM = (nchips - 1) * nm;

	/* allocate arrays for linear algebra */
	for (i = 0; i < 2; i++) {
        	B[i] = (double *) calloc(nM, sizeof(double));
        	A[i] = (double **) calloc(nM, sizeof(double *));
        	for (M = 0; M < nM; M++) {
                	A[i][M] = (double *) calloc(nM, sizeof(double));
       		}
	}
        indx = (int *) calloc(nM, sizeof(int));

	/* allocate space for the model parameters (with space for nm zeros)*/
	for (i = 0; i < 2; i++) {
		a[i] = (double *) calloc(nM + nm, sizeof(double));
	}

	/* set up arrays of l, m values */
	ll = (int *) calloc(nm, sizeof(int));
	mm = (int *) calloc(nm, sizeof(int));
	M = 0;
	for (l = 0; l <= lmax; l++) {
		for (m = 0; m <= l; m++) {
			ll[M] = l;
			mm[M++] = m;
		}			
	}

	/* allocate space for input record */
	inputrecord = (double *) calloc(6, sizeof(double));
	x[0] 	= inputrecord + 0;
	x[1] 	= inputrecord + 2;
	cd 	= inputrecord + 4;

	/* open lc pipe for input */
        if (!(lcpipe = popen("lc -b -o x c", "r"))) {
                fprintf(stderr, "fit2cats: readmergedcat: unable to open lc-pipe for input\n");
                exit(-1);
        }
		

	/* accumulate arrays */
	while (fread(inputrecord, sizeof(double), 6, lcpipe)) {
		c[0] = (int) floor(cd[0]);
		c[1] = (int) floor(cd[1]);
		if (!c[0] || !c[1]) {
			w = weightfac;
		} else {
			w = 1.0;
		}
		if ((c[0] >= nchips) || (c[0] < 0) || (c[1] >= nchips) || (c[1] < 0)) {
			error_exit("fitgeometry: chip number out of range\n");
		}
		for (i = 0; i < 2; i++) {
			for (m = 0; m < nm; m++) {
				if (c[0]) {
					B[i][I(c[0], m)] -= w * (x[0][i] - x[1][i]) * f(ll[m], mm[m], x[0]);
				}
				if (c[1]) {
					B[i][I(c[1], m)] -= w * (x[1][i] - x[0][i]) * f(ll[m], mm[m], x[1]);
				}
				for (n = 0; n < nm; n++) {
					if (c[0]) {
						A[i][I(c[0], m)][I(c[0], n)] += w * f(ll[m], mm[m], x[0]) * f(ll[n], mm[n], x[0]);
					}
					if (c[0] && c[1]) {
						A[i][I(c[1], m)][I(c[0], n)] -= w * f(ll[m], mm[m], x[1]) * f(ll[n], mm[n], x[0]);
						A[i][I(c[0], m)][I(c[1], n)] -= w * f(ll[m], mm[m], x[0]) * f(ll[n], mm[n], x[1]);
					}
					if (c[1]) {
						A[i][I(c[1], m)][I(c[1], n)] += w * f(ll[m], mm[m], x[1]) * f(ll[n], mm[n], x[1]);
					}
				}
			}
		}
	}

	for (i = 0; i < 2; i++) {
		/* solve the linear equations */
		myludcmp(A[i], nM, indx, &det);
		mylubksb(A[i], nM, indx, B[i]);
		/* extract the model coefficients */
		for (M = 0; M < nM; M++) {
			a[i][M] = B[i][M];
		}
	}

	/* write the 'null' parameter file */
	a[0] += nM;
	a[1] += nM;
	if (makeparfiles) {
		parfilename = (char *) calloc(1024, sizeof(char));
		sprintf(parfilename, "%s%d%s", parfilebase, 0, ".par");
	}
	write2Dpolymodel(parfilename, nm, ll, mm, 2, a, 1, vardef, "x");
	a[0] -= nM;
	a[1] -= nM;
	/* write the real parameter files */
	for (chip = 1; chip < nchips; chip++) {
		if (makeparfiles) {
			sprintf(parfilename, "%s%d%s", parfilebase, chip, ".par");
		}
		write2Dpolymodel(parfilename, nm, ll, mm, 2, a, 1, vardef, "x");
		a[0] += nm;
		a[1] += nm;
	}
	exit(0);
}


int	I(int c, int m)
{
	return (nm * (c - 1) + m);
}



