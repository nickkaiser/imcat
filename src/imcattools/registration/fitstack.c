/*
 * fitstack.c
 */


#define usage "\n\
NAME\n\
	fitstack --- fit for transformation coefficients for a set of exposures\n\
\n\
SYNOPSIS\n\
	fitstack nexp [options....]\n\
		-l lmin lmax		# min and max order for distortion polynomials\n\
		-i niter		# number of iterations (1)\n\
		-c outcat		# output catalogue containing magc, r-values\n\
		-o x0 y0		# origin for spatial coordinates (0, 0)\n\
		-d distparfilename	# file for distortion model parameters\n\
		-t transparfilename	# file for linear transformation params\n\
\n\
DESCRIPTION\n\
	'fitstack' reads from stdin a catalogue containing the result of\n\
	merging all pairs of cats for a stack of 'nexp' images (as created\n\
	by 'mergestack') and which must contain entries for\n\
	spatial coords 'x[2][2]', magnitude 'mag[2]' and exposure number 'exp[2]'.\n\
	It then fits a model in which sky coords (in frame defined by exposure-0) are\n\
		r = r_e + dphi_e r_e + d_e\n\
	where the 2x2 matrix dphi allows for rotations between\n\
	exposures and possibly atmospheric refraction, and we set\n\
	dphi = d = 0 for the 0th exposure.\n\
	It will also optionally then fit for distortion of telescope\n\
	using a model in which sky coords r are related to detector coords x by\n\
		r_e = x_e + sum a_m f_m(x_e)\n\
	where m labels the modes, and where each mode coefficient a_m\n\
	is a 2-vector and the modes are polynomials\n\
	of order 'lmin' through 'lmax'.\n\
	For lmin = 2, lmax = 4 say, the modes are:\n\
		x^2, xy, y^2, x^3, x^2 y, x y^2, y^3,\n\
		x^4, x^3 y, x^2 y^2, x y^3, y^4.\n\
	With -o option, these become polynomials in position relative\n\
	to specified spatial origin.\n\
	We also read magnitudes, which we model as:\n\
		m_e = m + M_e\n\
	where m is the true magnitude and M_e is the magnitude\n\
	offset the e'th exposure (relative to exp-0).\n\
	Use the '-c' option to generate a merged catalogue which contains,\n\
	in addition to the source catalogue values, the sky coordinates 'r'\n\
	and also corrected magnitudes 'magc', which can then be filtered\n\
	to remove bad pairs and then fed back to 'fitstack' to improve\n\
	the solution.\n\
	See also fitstack.tex.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@cita.utoronto.ca\n\
\n"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fitstack.h"
#include "utils/modefunc.h"
#include "fitstack_read.h"

#define SCALE 1.0

/* global nprs-length vectors of data */
double 		***x, **mag, ***re = NULL, ***r = NULL, ***z = NULL;
int		**e;

/* number of pairs, exposures, etc */
static	int	nprs, nexp, nmodes, lmin, lmax, *ll, *mm;
/* translation model parameters */
static	double	*d[2], *phi[2][2];
/* distortion model parameters */
static	double	*a[2];
/* magnitude offset parameters */
static	double	*dm;
/* matrix stuff 't' = translation; 'd' = distortion; 'm' = magnitudes */
static	int	*indxt, *indxd, *indxm, ntcoefft, ndcoefft, nmcoefft, dodistfit;
static	double	**At, *Bt, **Ad, *Bd, **Am, *Bm, det, *Ct[2], *Cd[2], *Cm;
	
main(int argc, char *argv[])
{
	/* arg counter, indices */
	int	arg = 2, i, j, l, m, n, M, niter;
	double	chi2, x0[2] = {0.0, 0.0};
	char	*outputcatfilename;
	int	dooutputrcat;
	/* parameter file */
	char	*distparfile, *transparfile;
	char	*vardef[1];

	/* defaults */
	dodistfit = 0;
	niter = 1;
	dooutputrcat = 0;
	distparfile = NULL;
	transparfile = NULL;
	vardef[0] = "1 1 a";

	/* parse args */
	if (argc < 2 || argv[1][0] == '-')
		error_exit(usage);
	if (1 != sscanf(argv[1], "%d", &nexp))
		error_exit(usage);
	while (arg < argc) {
		if (argv[arg][0] != '-') {
			error_exit(usage);
		}
		switch (argv[arg++][1]) {
			case 'l':
				dodistfit = 1;
				if (argc - arg < 2) {
					error_exit(usage);
				}
				if ((1 != sscanf(argv[arg++], "%d", &lmin)) ||
					(1 != sscanf(argv[arg++], "%d", &lmax))) {
						error_exit(usage);
				}
				break;
			case 'i':
				if (argc - arg < 1) {
					error_exit(usage);
				}
				if (1 != sscanf(argv[arg++], "%d", &niter)) {
						error_exit(usage);
				}
				break;
			case 'c':
				if (argc - arg < 1) {
					error_exit(usage);
				}
				outputcatfilename = argv[arg++];
				dooutputrcat = 1;
				break;
			case 'o':
				if (argc - arg < 2) {
					error_exit(usage);
				}
				if (1 != sscanf(argv[arg++], "%lf", &(x0[0])) ||
					1 != sscanf(argv[arg++], "%lf", &(x0[1]))) {
						error_exit(usage);
				}
				break;
			case 'd':
				if (argc - arg < 1) {
					error_exit(usage);
				}
				distparfile = argv[arg++];
				break;
			case 't':
				if (argc - arg < 1) {
					error_exit(usage);
				}
				transparfile = argv[arg++];
				break;
			case 'u':
			default:
				error_exit(usage);
		}
	}	
	setorigin(x0);

	modefunc_addargcomment(argc, argv);

	/* read the merged catalogue */
	nprs = readmergedcat(nexp);
	fprintf(stderr, "# read %d pairs\n", nprs);

	/* compute number of coefficients */
	ntcoefft = 6 * (nexp - 1);
	nmcoefft = nexp - 1;
	if (dodistfit) {
		nmodes = 0;
		for (l = lmin; l <= lmax; l++) {
			nmodes += (l + 1);
		}
		ndcoefft = 2 * nmodes;
	} else {
		nmodes = 0;
	}

	allocatearrays();

	/* set up arrays of l, m values */
	if (dodistfit) {
		ll = (int *) calloc(nmodes, sizeof(int));
		mm = (int *) calloc(nmodes, sizeof(int));
		M = 0;
/*		ll[M] = 1;
		mm[M++] = 1;*/
		for (l = lmin; l <= lmax; l++) {
			for (m = 0; m <= l; m++) {
				ll[M] = l;
				mm[M++] = m;
			}			
		}
	}


	make_r();
	chi2 = chisquared();
	fprintf(stderr, "# chisquared = %13.8lg\n", chi2);

	fitextinctions();
	for (i = 0; i < niter; i++) {
		make_re();
		fittranslations();
		make_z();
		fitdistortions();
		make_r();
		chi2 = chisquared();
		fprintf(stderr, "# chisquared = %13.8lg\n", chi2);
	}
	printtranslations(transparfile);
	write2Dpolymodel(distparfile, nmodes, ll, mm, 2, a, 1, vardef, "x");
	if (dooutputrcat) {
		outputrcat(outputcatfilename);
	}
	exit(0);
}



int	make_r(void)
{
	int	i, j, n, M, ipr;
	double	xx[2];

	/* allocate space if necessary */
	if (!r) {
		r = allocpositionvector();
	}

	for (ipr = 0; ipr < nprs; ipr++) {
		for (n = 0; n < 2; n++) {
			xx[0] = x[n][0][ipr];
			xx[1] = x[n][1][ipr];
			for (i = 0; i < 2; i++) {
				r[n][i][ipr] = x[n][i][ipr] + d[i][e[n][ipr]];
				for (j = 0; j < 2; j++) {
					r[n][i][ipr] += phi[i][j][e[n][ipr]] * x[n][j][ipr];
					for (M = 0; M < nmodes; M++) {
						if (i == j) {
							r[n][i][ipr] += a[i][M] * f(ll[M], mm[M], xx);
						}
						r[n][i][ipr] += phi[i][j][e[n][ipr]] * a[j][M] * f(ll[M], mm[M], xx);
					}
				}
			}
		}
	}
}


int	make_re(void)
{
	int	i, n, M, ipr;
	double	xx[2];

	/* allocate space if necessary */
	if (!re) {
		re = allocpositionvector();
	}

	for (ipr = 0; ipr < nprs; ipr++) {
		for (n = 0; n < 2; n++) {
			xx[0] = x[n][0][ipr];
			xx[1] = x[n][1][ipr];
			for (i = 0; i < 2; i++) {
				re[n][i][ipr] = x[n][i][ipr];
				for (M = 0; M < nmodes; M++) {
					re[n][i][ipr] += a[i][M] * f(ll[M], mm[M], xx);
				}
			}
		}
	}
}


int	make_z(void)
{
	int	i, j, n, ipr;

	/* allocate space if necessary */
	if (!z) {
		z = allocpositionvector();
	}

	for (ipr = 0; ipr < nprs; ipr++) {
		for (n = 0; n < 2; n++) {
			for (i = 0; i < 2; i++) {
				z[n][i][ipr] = x[n][i][ipr] + d[i][e[n][ipr]];
				for (j = 0; j < 2; j++) {
					z[n][i][ipr] += phi[i][j][e[n][ipr]] * x[n][j][ipr];
				}
			}
		}
	}
}


int	fitextinctions(void)
{
	int	i, j, ipr, n, m;

	/* zero the matrices and vectors */
	for (i = 0; i < nmcoefft; i++) {
		Bm[i] = 0.0;
		for (j = 0; j < nmcoefft; j++) {
			Am[i][j] = 0.0;
		}
	}

	for (ipr = 0; ipr < nprs; ipr++) {
		for (i = 0; i < nmcoefft; i++) {
			Cm[i] = 0.0;
		}
		for (m = 1; m < nexp; m++) {
			if (e[0][ipr] == m) {
				Cm[m - 1] += 1.0;
			}
			if (e[1][ipr] == m) {
				Cm[m - 1] -= 1.0;
			}
		}
		/* accumulate A matrix, B-vector */
		for (m = 0; m < nmcoefft; m++) {
			Bm[m] += (mag[0][ipr] - mag[1][ipr]) * Cm[m];
			for (n = 0; n < nmcoefft; n++) {
				Am[m][n] += Cm[m] * Cm[n];
			}
		}
	}
	/* solve the linear equations */
	myludcmp(Am, nmcoefft, indxm, &det);
	mylubksb(Am, nmcoefft, indxm, Bm);
	/* extract the model coefficients */
	for (m = 1; m < nexp; m++) {
		dm[m] = Bm[m - 1];
	}		
}



int	fittranslations(void)
{
	int	i, j, n, m, ipr, ne;

	ne = nexp - 1;

	/* zero the matrices and vectors */
	for (i = 0; i < ntcoefft; i++) {
		Bt[i] = 0.0;
		for (j = 0; j < ntcoefft; j++) {
			At[i][j] = 0.0;
		}
	}

	for (ipr = 0; ipr < nprs; ipr++) {
		for (i = 0; i < ntcoefft; i++) {
			Ct[0][i] = 0.0;
			Ct[1][i] = 0.0;
		}
		if (m = e[0][ipr]) {
			Ct[0][0 * ne + m - 1] += 1.0;
			Ct[1][1 * ne + m - 1] += 1.0;
			Ct[0][2 * ne + m - 1] += re[0][0][ipr];
			Ct[0][3 * ne + m - 1] += re[0][1][ipr];
			Ct[1][4 * ne + m - 1] += re[0][0][ipr];
			Ct[1][5 * ne + m - 1] += re[0][1][ipr];
		}
		if (m = e[1][ipr]) {
			Ct[0][0 * ne + m - 1] -= 1.0;
			Ct[1][1 * ne + m - 1] -= 1.0;
			Ct[0][2 * ne + m - 1] -= re[1][0][ipr];
			Ct[0][3 * ne + m - 1] -= re[1][1][ipr];
			Ct[1][4 * ne + m - 1] -= re[1][0][ipr];
			Ct[1][5 * ne + m - 1] -= re[1][1][ipr];
		}
		/* accumulate A matrix, B-vector */
		for (i = 0; i < 2; i++) {
			for (m = 0; m < ntcoefft; m++) {
				Bt[m] -= (re[0][i][ipr] - re[1][i][ipr]) * Ct[i][m];
				for (n = 0; n < ntcoefft; n++) {
					At[m][n] += Ct[i][m] * Ct[i][n];
				}
			}
		}
	}
	/* solve the linear equations */
	myludcmp(At, ntcoefft, indxt, &det);
	mylubksb(At, ntcoefft, indxt, Bt);
	/* extract the model coefficients */
	for (m = 1; m < nexp; m++) {
		d[0][m] = Bt[0 * ne + m - 1];
		d[1][m] = Bt[1 * ne + m - 1];
		phi[0][0][m] = Bt[2 * ne + m - 1];
		phi[0][1][m] = Bt[3 * ne + m - 1];
		phi[1][0][m] = Bt[4 * ne + m - 1];
		phi[1][1][m] = Bt[5 * ne + m - 1];
	}		
}



int	fitdistortions(void)
{
	int	i, j, n, m, M, ipr;
	double	x0[2], x1[2];

	/* zero the matrices and vectors */
	for (i = 0; i < ndcoefft; i++) {
		Bd[i] = 0.0;
		for (j = 0; j < ndcoefft; j++) {
			Ad[i][j] = 0.0;
		}
	}

	for (ipr = 0; ipr < nprs; ipr++) {
		for (M = 0; M < nmodes; M++) {
			x0[0] = x[0][0][ipr];
			x0[1] = x[0][1][ipr];
			x1[0] = x[1][0][ipr];
			x1[1] = x[1][1][ipr];
			Cd[0][0 * nmodes + M] = (1 + phi[0][0][e[0][ipr]]) * f(ll[M], mm[M], x0) -
						(1 + phi[0][0][e[1][ipr]]) * f(ll[M], mm[M], x1);
			Cd[0][1 * nmodes + M] = phi[0][1][e[0][ipr]] * f(ll[M], mm[M], x0) -
						phi[0][1][e[1][ipr]] * f(ll[M], mm[M], x1);
			Cd[1][0 * nmodes + M] = phi[1][0][e[0][ipr]] * f(ll[M], mm[M], x0) -
						phi[1][0][e[1][ipr]] * f(ll[M], mm[M], x1);
			Cd[1][1 * nmodes + M] = (1 + phi[1][1][e[0][ipr]]) * f(ll[M], mm[M], x0) -
						(1 + phi[1][1][e[1][ipr]]) * f(ll[M], mm[M], x1);
		}
		/* accumulate A matrix, B-vector */
		for (i = 0; i < 2; i++) {
			for (m = 0; m < ndcoefft; m++) {
				Bd[m] -= (z[0][i][ipr] - z[1][i][ipr]) * Cd[i][m];
				for (n = 0; n < ndcoefft; n++) {
					Ad[m][n] += Cd[i][m] * Cd[i][n];
				}
			}
		}
	}
	/* solve the linear equations */
	myludcmp(Ad, ndcoefft, indxd, &det);
	mylubksb(Ad, ndcoefft, indxd, Bd);
	/* extract the model coefficients */
	for (M = 0; M < nmodes; M++) {
		a[0][M] = Bd[0 * nmodes + M];
		a[1][M] = Bd[1 * nmodes + M];
	}		
}



int	allocatearrays(void)
{
	int	i, j, m;

	/* allocate arrays for linear algebra */
        Bt = (double *) calloc(ntcoefft, sizeof(double));
        Ct[0] = (double *) calloc(ntcoefft, sizeof(double));
        Ct[1] = (double *) calloc(ntcoefft, sizeof(double));
        At = (double **) calloc(ntcoefft, sizeof(double *));
        for (m = 0; m < ntcoefft; m++) {
                At[m] = (double *) calloc(ntcoefft, sizeof(double));
        }
        indxt = (int *) calloc(ntcoefft, sizeof(int));
        Bm = (double *) calloc(nmcoefft, sizeof(double));
        Cm = (double *) calloc(nmcoefft, sizeof(double));
        Am = (double **) calloc(nmcoefft, sizeof(double *));
        for (m = 0; m < nmcoefft; m++) {
                Am[m] = (double *) calloc(nmcoefft, sizeof(double));
        }
        indxm = (int *) calloc(nmcoefft, sizeof(int));
	if (dodistfit) {
        	Bd = (double *) calloc(ndcoefft, sizeof(double));
        	Cd[0] = (double *) calloc(ndcoefft, sizeof(double));
        	Cd[1] = (double *) calloc(ndcoefft, sizeof(double));
        	Ad = (double **) calloc(ndcoefft, sizeof(double *));
        	for (m = 0; m < ndcoefft; m++) {
               		Ad[m] = (double *) calloc(ndcoefft, sizeof(double));
       		}
        	indxd = (int *) calloc(ndcoefft, sizeof(int));
	}

	/* allocate space for the model parameters */
	for (i = 0; i < 2; i++) {
		d[i] = (double *) calloc(nexp, sizeof(double));
		if (dodistfit) {
			a[i] = (double *) calloc(nmodes, sizeof(double));
		}
		for (j = 0; j < 2; j++) {
			phi[i][j] = (double *) calloc(nexp, sizeof(double));
		}
	}
	dm = (double *) calloc(nexp, sizeof(double));
}



double	***allocpositionvector(void)
{
	double	***rr;
	int	i, j;

	rr = (double ***) calloc(2, sizeof(double **));
	if (!rr)
		error_exit("fitstack: allocpositionvector: memory allocation failure\n");
	for (i = 0; i < 2; i++) {
		rr[i] = (double **) calloc(2, sizeof(double *));
		if (!(rr[i]))
			error_exit("fitstack: allocpositionvector: memory allocation failure\n");
		for (j = 0; j < 2; j++) {
			rr[i][j] = (double *) calloc(nprs, sizeof(double));
			if (!(rr[i][j]))
				error_exit("fitstack: allocpositionvector: memory allocation failure\n");
		}
	}
	return(rr);
}

int	printtranslations(char	*parfile)
{
	int	i;
	FILE	*parf;

	if (!(parf = fopen(parfile, "w"))) {
		error_exit("printtranslations: unable to open parameter file\n");
	}
	fprintf(parf, "%d exposures\n", nexp);
	fprintf(parf, "#        d[0]          d[1]     phi[0][0]     phi[0][1]     phi[1][0]     phi[1][1]            dm\n");
	for (i = 0; i < nexp; i++) {
		fprintf(parf, "%13.8lg %13.8lg %13.8lg %13.8lg %13.8lg %13.8lg %13.8lg\n",
			d[0][i], d[1][i], 1 + phi[0][0][i], phi[0][1][i], phi[1][0][i], 1 + phi[1][1][i], dm[i]);  
	}
}




double	chisquared(void)
{
	int	i, ipr;
	double 	chi2;

	chi2 = 0.0;
	for (ipr = 0; ipr < nprs; ipr++) {
		for (i = 0; i < 2; i++) {
			chi2 += (r[0][i][ipr] - r[1][i][ipr]) * (r[0][i][ipr] - r[1][i][ipr]);
		}
	}

	return(chi2);
}


int	outputrcat(char *outputcatfilename)
{
	FILE	*lcpipe;
	int	ipr;
	char	outpipe[128];

	sprintf(outpipe, "lc -C -N '2 2 2 x' -N '2 2 2 r' -N '1 2 mag' -N '1 2 magc' -N '1 2 exp' > %s", outputcatfilename);
	lcpipe = popen(outpipe, "w");
	if (!lcpipe) {
		fprintf(stderr, "fitstack: failed to open %s\n", outpipe);
		exit(1);
	}
	for (ipr = 0; ipr < nprs; ipr++) {
		fprintf(lcpipe, "%13.8lg %13.8lg %13.8lg %13.8lg %13.8lg %13.8lg %13.8lg %13.8lg %13.8lg %13.8lg %13.8lg %13.8lg %13.8lg %13.8lg\n", 
			x[0][0][ipr], x[0][1][ipr], x[1][0][ipr], x[1][1][ipr],
			r[0][0][ipr], r[0][1][ipr], r[1][0][ipr], r[1][1][ipr],
			mag[0][ipr], mag[1][ipr], mag[0][ipr] - dm[e[0][ipr]], mag[1][ipr] - dm[e[1][ipr]],
			(double) e[0][ipr], (double) e[1][ipr]);
	}
	pclose(outpipe);
}
