/*
 * fitdistortion.c
 */


#define usage "\n\
NAME\n\
	fitdistortion --- fit for field distortion\n\
\n\
SYNOPSIS\n\
	fitdistortion ne [options...] \n\
		-l lmax		# max order for field distortion polynomials (3)\n\
		-d parfilebase	# basename for transformation param files\n\
		-o x0 y0	# spatial origin for mode functions\n\
\n\
DESCRIPTION\n\
	'fitdistortion' reads from stdin a set of pairs of coords x[2][2]\n\
	and exposure numbers e[2] (as prodced by merging a set of catalogues\n\
	pair by pair, and with 0 <= e < ne) and fits a model in which sky\n\
	coord 'r' is related to chip coordinate 'x' by\n\
\n\
		r = x + sum_M a_M f_M(x) + sum_N a_{eN} f_N(x)\n\
\n\
	where the f_M's are polynomials of order l = 2 to l_max and\n\
	describe the field distortion and the f_N's are polynomials\n\
	of order 0, 1 (i.e. linear transformations) which describe the\n\
	telescope pointing as well as possible scale changes, or shear\n\
	introduced by atmospheric refraction etc.  With the -o option\n\
	these become polynomials position relative to given spatial\n\
	origin. We let exposure e = 0 define the orientation and scale\n\
	of the image (so a_{0N} = 0).\n\
\n\
	The model is only approximate; really we should have:\n\
\n\
		r = x + sum_M a_M f_M(x) + sum_N a_{eN} f_N(x + sum_M a_M f_M(x))\n\
\n\
	which is equivalent to the simpler form above in the limit\n\
	of small mode coefficient amplitudes (i.e small rotations between\n\
	fields etc.\n\
	By default, the transformation parameters a_{lm} are written to stdout\n\
	as a concatenation of '.par' format files (which is not particularly\n\
	useful), but with the -d option you can specify a basename (which\n\
	might be a preexisting directory) and fitgeometry will create a set of files\n\
	'parfilebase''e'.par with e = 0... ne-1 containing the coefficients\n\
	a_{eN} describing the pointings (a 'null' parameter file for the\n\
	e = 0 exposure is provided for convenience) and parfilebase'dist.par\n\
	containing the distortion parameters a_M.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@ifa.hawaii.edu\n\
\n"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utils/modefunc.h"
#include "utils/lu.h"

main(int argc, char *argv[])
{
	/* number of modes, exposures, etc */
	int nM, nm, ne, lmax, *ll, *mm, *indx, e0, e1, *ee;
	/* model parameters */
	static	double	*a[2];
	/* matrix stuff */
	static	double	**A[2], *B[2], det;
	/* input pipe, temp par file */
	FILE	*lcpipe, *tempf;
	/* arg counter */
	int	arg = 1;
	/* miscellaneous indices */
	int	i, e, l, m, M, N;
	/* input record */
	double *inputrecord, *x[2], *ed;
	/* output stuff */
	char	*parfilebase, *parfilename;
	int	makeparfiles;
	/* mode func values */
	double	fM0, fM1, fN0, fN1;
	/* spatial origin */
	double	x0[2];
	char	*vardef[1];

	/* defaults */
	lmax = 3;
	makeparfiles = 0;
	parfilename = NULL;
	vardef[0] = "1 2 a";

	/* get the number of exposures */
	if ((1 != sscanf(argv[arg++], "%d", &ne))) {
			error_exit(usage);
	}

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
			case 'd':
				if (argc - arg < 1) {
					error_exit(usage);
				}
				parfilebase = argv[arg++];
				makeparfiles = 1;
				break;
			case 'o':
				if (argc - arg < 1) {
					error_exit(usage);
				}
				if ((1 != sscanf(argv[arg++], "%lf", &(x0[0]))) ||
					(1 != sscanf(argv[arg++], "%lf", &(x0[1])))) {
						error_exit(usage);
				}
				setorigin(x0);
				break;
			case 'u':
			default:
				error_exit(usage);
		}
	}

	/* first check that we can write to the parfiledir if neccessary */
	if (makeparfiles) {
		parfilename = (char *) calloc(128, sizeof(char));
		sprintf(parfilename, "%s%s", parfilebase, "dist.par");
		if (!(tempf = fopen(parfilename, "w"))) {
			error_exit("fitdistortion: unable to open par file for output\n");
		}
		/* and immediately close and remove it */
		close(tempf);
		remove(parfilename);
	}

	/* copy args to modefunc argstring */
	modefunc_addargcomment(argc, argv);


	/* compute number of modes*/
	nm = 0;
	for (l = 2; l <= lmax; l++) {
		nm += (l + 1);
	}
	nM = nm + 3 * (ne - 1);

	/* allocate arrays for linear algebra */
	for (i = 0; i < 2; i++) {
        	B[i] = (double *) calloc(nM, sizeof(double));
        	A[i] = (double **) calloc(nM, sizeof(double *));
        	for (M = 0; M < nM; M++) {
                	A[i][M] = (double *) calloc(nM, sizeof(double));
       		}
	}
        indx = (int *) calloc(nM, sizeof(int));

	/* allocate space for the model parameters */
	for (i = 0; i < 2; i++) {
		a[i] = (double *) calloc(nM + 3, sizeof(double));
	}

	/* set up arrays of l, m, e values */
	ll = (int *) calloc(nM, sizeof(int));
	mm = (int *) calloc(nM, sizeof(int));
	ee = (int *) calloc(nM, sizeof(int));
	M = 0;
	for (l = 2; l <= lmax; l++) {
		for (m = 0; m <= l; m++) {
			ll[M] = l;
			mm[M++] = m;
		}			
	}
	for (e = 1; e < ne; e++) {
		for (l = 0; l <= 1; l++) {
			for (m = 0; m <= l; m++) {
				ee[M] = e;
				ll[M] = l;
				mm[M++] = m;
			}
		}
	}
	

	/* allocate space for input record */
	inputrecord = (double *) calloc(6, sizeof(double));
	x[0] 	= inputrecord + 0;
	x[1] 	= inputrecord + 2;
	ed 	= inputrecord + 4;

	/* open lc pipe for input */
        if (!(lcpipe = popen("lc -b -o x e", "r"))) {
                fprintf(stderr, "fit2cats: readmergedcat: unable to open lc-pipe for input\n");
                exit(-1);
        }
		

	/* accumulate arrays */
	while (fread(inputrecord, sizeof(double), 6, lcpipe)) {
		e0 = (int) floor(ed[0]);
		e1 = (int) floor(ed[1]);
		if ((e0 >= ne) || (e0 < 0) || (e1 >= ne) || (e1 < 0)) {
			error_exit("fitgeometry: exp number out of range\n");
		}
		for (i = 0; i < 2; i++) {
			for (M = 0; M < nM; M++) {
				if (M < nm) {
					fM0 = f(ll[M], mm[M], x[0]);
					fM1 = f(ll[M], mm[M], x[1]);
				} else {
					fM0 = (e0 == ee[M] ? f(ll[M], mm[M], x[0]) : 0.0);
					fM1 = (e1 == ee[M] ? f(ll[M], mm[M], x[1]) : 0.0);
				}
				B[i][M] -= (x[0][i] - x[1][i]) * (fM0 - fM1);
				for (N = 0; N < nM; N++) {
					if (N < nm) {
						fN0 = f(ll[N], mm[N], x[0]);
						fN1 = f(ll[N], mm[N], x[1]);
					} else {
						fN0 = (e0 == ee[N] ? f(ll[N], mm[N], x[0]) : 0.0);
						fN1 = (e1 == ee[N] ? f(ll[N], mm[N], x[1]) : 0.0);
					}
					A[i][M][N] += (fM0 - fM1) * (fN0 - fN1);
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

	if (makeparfiles) {
		sprintf(parfilename, "%s%s", parfilebase, "dist.par");
	}
	write2Dpolymodel(parfilename, nm, ll, mm, 2, a, 1, vardef, "x");
	ll += nm;
	mm += nm;
	for (e = 0; e < ne; e++) {
		if (makeparfiles) {
			sprintf(parfilename, "%s%d%s", parfilebase, e, ".par");
		}
		if (e == 0) {
			a[0] += nM;
			a[1] += nM;
			write2Dpolymodel(parfilename, 3, ll, mm, 2, a, 1, vardef, "x");
			a[0] -= (nM - nm);
			a[1] -= (nM - nm);			
		} else {
			write2Dpolymodel(parfilename, 3, ll, mm, 2, a, 1, vardef, "x");
			ll += 3;
			mm += 3;
			a[0] += 3;
			a[1] += 3;			
		}
	}

	exit(0);
}

