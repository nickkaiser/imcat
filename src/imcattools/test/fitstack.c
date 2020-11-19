/*
 * fitstack.c
 */


#define usage "\n\
NAME\n\
	fitstack --- fit for transformation coefficients for a set of exposures\n\
\n\
SYNOPSIS\n\
	fitstack nexp [options....]\n\
		-l lmin lmax	# min and max order for distortion polynomials\n\
		-t ftol		# fractional tolerance (1.e-10)\n\
\n\
DESCRIPTION\n\
	'fitstack' reads from stdin a catalogue containing the result of\n\
	merging all pairs of cats for a stack of 'nexp' images (as created\n\
	by 'mergestacks1' or 'mergestacks2') and which must contain entries for\n\
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
	of order 'lmin' through 'lmax'. 
	For lmin = 2, lmax = 4 say, the modes are:\n\
		x^2, xy, y^2, x^3, x^2 y, x y^2, y^3,\n\
		x^4, x^3 y, x^2 y^2, x y^3, y^4.\n\
	We also read magnitudes, which we model as:\n\
		m_e = m + M_e\n\
	where m is the true magnitude and M_e is the magnitude\n\
	offset the e'th exposure (relative to exp-0).\n\
	See also fitstack.tex.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@cita.utoronto.ca\n\
\n"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fitstack.h"
#include "fitstack_modefunc.h"
#include "fitstack_read.h"
#include "amoeba_dbl.h"
#include "amotry_dbl.h"

#define SCALE 1.0

/* global nprs-length vectors of data */
double 		***x, **mag, ***r = NULL, ***re;
int		**e;

/* number of pairs, exposures, etc */
static	int	nprs, nexp, nmodes, ncoefft, lmin, lmax, *ll, *mm, dodistfit;
/* translation model parameters */
static	double	*d[2], *phi[2][2];
/* distortion model parameters */
static	double	*a[2];

main(int argc, char *argv[])
{
	/* arg counter, indices */
	int	arg = 2, i, j, l, m, n, M, ne;
	double	chi2;
	/* stuff for powell - v = vector of coefficients */
	double	**xi, ftol, *v, fret;
	int	iter;

	/* defaults */
	dodistfit = 0;
	ftol = 1.e-10;

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
			case 't':
				dodistfit = 1;
				if (argc - arg < 1) {
					error_exit(usage);
				}
				if ((1 != sscanf(argv[arg++], "%lf", &ftol))) {
						error_exit(usage);
				}
				break;
			case 'u':
			default:
				error_exit(usage);
		}
	}	

	/* read the merged catalogue */
	nprs = readmergedcat(nexp);

	re = x;

	/* compute number of coefficients */
	ncoefft = 6 * (nexp - 1);
	nmodes = 0;
	if (dodistfit) {
		for (l = lmin; l <= lmax; l++) {
			nmodes += (l + 1);
		}
	}
	ncoefft += 2 * nmodes;
	allocatearrays();


	/* set up arrays of l, m values */
	if (dodistfit) {
		ll = (int *) calloc(nmodes, sizeof(int));
		mm = (int *) calloc(nmodes, sizeof(int));
		M = 0;
		for (l = lmin; l <= lmax; l++) {
			for (m = 0; m <= l; m++) {
				ll[M] = l;
				mm[M++] = m;
			}			
		}
	}


	v = (double *) calloc(ncoefft, sizeof(double));

	make_r();
	chi2 = chisquared(v);
	fprintf(stderr, "chisquared = %13.8lg\n", chi2);

	fittranslations();
	printtranslations();

	ne = nexp - 1;
	for (n = 0; n < ne; n++) {
		v[0 * ne + n] = d[0][n+1];
		v[1 * ne + n] = d[1][n+1];
		v[2 * ne + n] = phi[0][0][n+1];
		v[3 * ne + n] = phi[0][1][n+1];
		v[4 * ne + n] = phi[1][0][n+1];
		v[5 * ne + n] = phi[1][1][n+1];
	}

	make_r();
	chi2 = chisquared(v);
	fprintf(stderr, "chisquared = %13.8lg\n", chi2);

	xi = (double **) calloc(ncoefft, sizeof(double *));
	for (i = 0; i < ncoefft; i++) {
		xi[i] = (double *) calloc(ncoefft, sizeof(double));
		xi[i][i] = 1.0;
	}
	powell(v, xi, ncoefft, ftol, &iter, &fret, chisquared);	
	fprintf(stderr, "%d iterations\n", iter);
 
	chi2 = chisquared(v);
	fprintf(stderr, "chisquared = %13.8lg\n", chi2);
	printtranslations();
	printdistortions();
	outputrcat();
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





int	allocatearrays(void)
{
	int	i, j, m;

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

int	printtranslations(void)
{
	int	i;

	for (i = 0; i < nexp; i++) {
		fprintf(stderr, "%13.8lg %13.8lg %13.8lg %13.8lg %13.8lg %13.8lg\n",
			d[0][i], d[1][i], phi[0][0][i], phi[0][1][i], phi[1][0][i], phi[1][1][i]);  
	}
}


int	printdistortions(void)
{
	int	M;

	for (M = 0; M < nmodes; M++) {
		fprintf(stderr, "%3d %3d %13.8lg %13.8lg\n", ll[M], mm[M], a[0][M], a[1][M]);  
	}
}


double	chisquared(double *v)
{
	int	ne, n, m, i, j, ipr;
	double 	chi2;

	ne = nexp - 1;
	for (n = 0; n < ne; n++) {
		d[0][n+1] = v[0 * ne + n];
		d[1][n+1] = v[1 * ne + n];
		phi[0][0][n+1] = v[2 * ne + n];
		phi[0][1][n+1] = v[3 * ne + n];
		phi[1][0][n+1] = v[4 * ne + n];
		phi[1][1][n+1] = v[5 * ne + n];
	}
	for (m = 0; m < nmodes; m++) {
		a[0][m] = v[6 * ne + m];
		a[1][m] = v[6 * ne + nmodes + m];
	}

	make_r();
	chi2 = 0.0;
	for (ipr = 0; ipr < nprs; ipr++) {
		for (i = 0; i < 2; i++) {
			chi2 += (r[0][i][ipr] - r[1][i][ipr]) * (r[0][i][ipr] - r[1][i][ipr]);
		}
	}

	return(chi2);
}


int	outputrcat(void)
{
	FILE	*lcpipe;
	int	ipr;


	lcpipe = popen("lc -C -N '2 2 2 x' -N '1 2 exp' ", "w");
	for (ipr = 0; ipr < nprs; ipr++) {
		fprintf(lcpipe, "%13.8lg %13.8lg %13.8lg %13.8lg %13.8lg %13.8lg\n", 
			r[0][0][ipr], r[0][1][ipr], r[1][0][ipr], r[1][1][ipr],
			(double) e[0][ipr], (double) e[1][ipr]);
	}
}



int	fittranslations(void)
{
	int	i, j, n, m, ipr, ne;
	int	*indxt, ntcoefft;
	double	**At, *Bt, det, *Ct[2];

	ne = nexp - 1;
	ntcoefft = 6 * ne;

        Bt = (double *) calloc(ntcoefft, sizeof(double));
        Ct[0] = (double *) calloc(ntcoefft, sizeof(double));
        Ct[1] = (double *) calloc(ntcoefft, sizeof(double));
        At = (double **) calloc(ntcoefft, sizeof(double *));
        for (m = 0; m < ntcoefft; m++) {
                At[m] = (double *) calloc(ntcoefft, sizeof(double));
        }
        indxt = (int *) calloc(ntcoefft, sizeof(int));	

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
		for (m = 1; m < nexp; m++) {
			if (e[0][ipr] == m) {
				Ct[0][0 * ne + m - 1] += 1.0;
				Ct[1][1 * ne + m - 1] += 1.0;
				Ct[0][2 * ne + m - 1] += re[0][0][ipr];
				Ct[0][3 * ne + m - 1] += re[0][1][ipr];
				Ct[1][4 * ne + m - 1] += re[0][0][ipr];
				Ct[1][5 * ne + m - 1] += re[0][1][ipr];
			}
			if (e[1][ipr] == m) {
				Ct[0][0 * ne + m - 1] -= 1.0;
				Ct[1][1 * ne + m - 1] -= 1.0;
				Ct[0][2 * ne + m - 1] -= re[1][0][ipr];
				Ct[0][3 * ne + m - 1] -= re[1][1][ipr];
				Ct[1][4 * ne + m - 1] -= re[1][0][ipr];
				Ct[1][5 * ne + m - 1] -= re[1][1][ipr];
			}
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
