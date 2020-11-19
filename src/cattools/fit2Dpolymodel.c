/*
 * fit2Dpolymodel.c
 */


#define usage "\n\
NAME\n\
	fit2Dpolymodel --- fit a 2-dimensional spatial polynomial model to a cat\n\
\n\
SYNOPSIS\n\
	fit2Dpolymodel xname lmin lmax var [options...] \n\
		-o x0 y0	# spatial origin for mode functions\n\
\n\
DESCRIPTION\n\
	'fit2Dpolymodel' reads from stdin a catalogue which must contain\n\
	at least a 2-vector 'xname[2]' and a numerical variable 'var'\n\
	(which may be a scalar, vector or tensor of arbitrary rank)\n\
	and determines the coefficients of a 2-D spatial poly model:\n\
\n\
		var_model = sum_l sum_m a_lm f_lm(x)\n\
\n\
	where the f_lm's are 2-D polynomials:\n\
\n\
		f_lm(x) = x[0]^(l - m) x[1]^m\n\
\n\
	with lmin <= l <= lmax and 0 <= m <= l.\n\
\n\
	The model coefficients are output as an lc-format catalogue\n\
	containing the indices l, m, followed by the mode coefficients\n\
	which have the same names as the 'var' variable.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@ifa.hawaii.edu\n\
\n"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../utils/modefunc.h"
#include "../utils/lu.h"

main(int argc, char *argv[])
{
	/* number of modes  etc */
	int nmodes, lmin, lmax, *ll, *mm, *indx;
	/* model parameters */
	static	double	**a;
	/* matrix stuff */
	static	double	***A, **B, det;
	/* input pipe */
	FILE	*lcpipe;
	/* arg counter */
	int	arg = 1;
	/* miscellaneous indices */
	int	i, l, m, N, M, asize, nvar;
	/* input record */
	double *inputrecord, *x, *var, fN, fM;
	/* spatial origin */
	double	x0[2];
	char	*vardef[MODEFUNC_MAX_VARS], *xname, *varname, lccommand[1024];

	if (argc < 5) {
		error_exit(usage);
	}

	/* get required arguments */
	xname = argv[arg++];
	if ((1 != sscanf(argv[arg++], "%d", &lmin))) {
			error_exit(usage);
	}
	if ((1 != sscanf(argv[arg++], "%d", &lmax))) {
			error_exit(usage);
	}
	varname = argv[arg++];

	/* parse args */
	while (arg < argc) {
		if (argv[arg][0] != '-') {
			error_exit(usage);
		}
		switch (argv[arg++][1]) {
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

	/* copy args to modefunc argstring */
	modefunc_addargcomment(argc, argv);

	/* compute number of modes*/
	nmodes = 0;
	for (l = lmin; l <= lmax; l++) {
		nmodes += (l + 1);
	}

	/* open lc-pipe for input */
	sprintf(lccommand, "lc -b %s %s", xname, varname);
	lcpipe = popen(lccommand, "r");
	if (!lcpipe) {
		error_exit("fit2Dpolymodel: failed to open lc-pipe for input\n");
	}
	getvars(lcpipe, &nvar, vardef, &asize);
	asize -= 2;
	if (strncmp(vardef[0], "1 2", 3)) {
		error_exit("fit2Dpolymodel: 'xname' should be a 2-vector!\n");
	}

	/* allocate arrays for linear algebra */
	B = (double **) calloc(asize, sizeof(double *));
	A = (double ***) calloc(asize, sizeof(double **));
	for (i = 0; i < asize; i++) {
        	B[i] = (double *) calloc(nmodes, sizeof(double));
        	A[i] = (double **) calloc(nmodes, sizeof(double *));
        	for (M = 0; M < nmodes; M++) {
                	A[i][M] = (double *) calloc(nmodes, sizeof(double));
       		}
	}
        indx = (int *) calloc(nmodes, sizeof(int));

	/* allocate space for the model parameters */
	a = (double **) calloc(asize, sizeof(double *));
	for (i = 0; i < asize; i++) {
		a[i] = (double *) calloc(nmodes, sizeof(double));
	}

	/* set up arrays of l, m values */
	ll = (int *) calloc(nmodes, sizeof(int));
	mm = (int *) calloc(nmodes, sizeof(int));
	M = 0;
	for (l = lmin; l <= lmax; l++) {
		for (m = 0; m <= l; m++) {
			ll[M] = l;
			mm[M++] = m;
		}			
	}
	

	/* allocate space for input record */
	inputrecord = (double *) calloc(2 + asize, sizeof(double));
	x = inputrecord;
	var = inputrecord + 2;



	/* accumulate arrays */
	while (fread(inputrecord, sizeof(double), 2 + asize, lcpipe)) {
		for (i = 0; i < asize; i++) {
			for (M = 0; M < nmodes; M++) {
				fM = f(ll[M], mm[M], x);
				B[i][M] += var[i] * fM;
				for (N = 0; N < nmodes; N++) {
					fN = f(ll[N], mm[N], x);
					A[i][M][N] += fM * fN;
				}

			}
		}
	}

	for (i = 0; i < asize; i++) {
		/* solve the linear equations */
		myludcmp(A[i], nmodes, indx, &det);
		mylubksb(A[i], nmodes, indx, B[i]);
		/* extract the model coefficients */
		for (M = 0; M < nmodes; M++) {
			a[i][M] = B[i][M];
		}
	}

	write2Dpolymodel(NULL, nmodes, ll, mm, asize, a, 1, vardef + 1, xname);

	exit(0);
}

