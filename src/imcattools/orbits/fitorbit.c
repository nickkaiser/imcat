#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vectors.h"
#include "Ffunc.h"
#include "gaussdev.h"
#include "tcl.h"
#include "kepler.h"

#define usage "\nNAME\n\
	fitorbit - find least squares orbit solution\n\
\n\
SYNOPSIS\n\
	fitorbit initsolcat np [-tcl dtmax]\n\
\n\
DESCRIPTION\n\
	fitorbit finds a minimum chi-squared solution for a set of observations\n\
	given an initial solution.\n\
\n\
	fitorbit first reads an lc format catalog initsolcat containing the\n\
	a preliminary solution ra0[3], va0[3], as well as the position re0[3] and\n\
	velocity ve0[3] of the Earth at t0 the initial time t0 (though these just\n\
	get passed through.  It also contains the number of observations nt.\n\
\n\
	It then reads a set of nt observations re[3], rho[3], nobs[3], sigma\n\
	and t from obscat, and finds a minimum chi^2 by propagating from the\n\
	t0 solution and computing the residuals.  By default it uses Kepler elements\n\
	but with -tcl flag will use time-centered leapfrog with maximum timestep\n\
	dtmax.  Minimum is found using Powell's method.\n\
\n\
	Here re[] is the geocenter, rho[] is the position of the observatory\n\
	realtive to the geocenter and nobs is a unit vector in the direction\n\
	ra[] - re[].\n\
\n\
	fitorbit outputs a series of np+2 sets of phase-space coordinates\n\
	r[], v[] where the first set is the observer, the second is the\n\
	asteroid solution and subsequent np sets of coordinates are solutions\n\
	obtained by randomly perturbing the observations with Gaussian\n\
	distributed errors with rms value sigma.  Also output is the minimum\n\
	value of chi-squared for the solution, and a particle number p (-1 for\n\
	the observer, 0 for the actual solution and 1 ... np for the virtual\n\
	asteroids.\n\
\n\
SEE ALSO\n\
	makeobs_inertial makeobs_circ laplace3 tcl_evolve r2n tcl_evolveN orbs2obs\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@hawaii.edu\n\n"

/* globals */
double	**re, **rho, **nobs, **nin, *n, *t, t0, *sigma, *ra, *va, arcsec, dtmax, *iota[2], *n0;
int	nt, dotcl;
/* powell stuff */
double	**powell_xi, powell_ftol, powell_fret;
/* amoeba stuff */
double	**amoeba_p, *amoeba_y, amoeba_ftol;
/* the orbit structure */
keplerorbit	*theorbit;
/* derivatives */
double	*chi_i, **chi_ij, **gx;

double	chisquared(double *x);
void 	powell(double p[], double **xi, int n, double ftol, int *iter, double *fret,
        double (*func)(double []));
void 	amoeba(double **p, double y[], int ndim, double ftol, double (*funk)(double []), int *nfunk);
void	findmin_powell(double *r0, double *v0, double *x);
void	findmin_amoeba(double *r0, double *v0, double *x, double dx);
void	getbasis(double *n, double *i0, double *i1);
void	getderivs(double *r0, double *v0, double dx);
void	getchiplanes(double *x, int N, double dx);

main (int argc, char *argv[])
{
	double	*ra0, *va0, *re0, *ve0, dtmp, chi2min, *x;
	int	i, j, it, ip, np;
	/* I/O stuff */
	FILE	*ipf, *initsolcatf, *opf;
	char	*initsolcatfilename, lccom[1024];

	/* defaults */
	dotcl = 0;

	/* parse args */
	if ((argc != 3) && (argc != 5)) {
		fprintf(stderr, usage);
		exit(-1);
	}
	initsolcatfilename = argv[1];
	if (1 != sscanf(argv[2], "%d", &np)){
		fprintf(stderr, usage);
		exit(-1);
	}
	if (argc == 5) {
		dotcl = 1;
		if (1 != sscanf(argv[4], "%lf", &dtmax)){
			fprintf(stderr, usage);
			exit(-1);
		}
	}		

	/* arc-second in radians */
	arcsec = M_PI / (180.0 * 3600.0);

	/* allocate the initial solution */
	ra0 = (double *) calloc(3, sizeof(double));
	va0 = (double *) calloc(3, sizeof(double));
	re0 = (double *) calloc(3, sizeof(double));
	ve0 = (double *) calloc(3, sizeof(double));
	/* and the internal arrays */
	ra = (double *) calloc(3, sizeof(double));
	va = (double *) calloc(3, sizeof(double));
	n  = (double *) calloc(3, sizeof(double));
	/* the phase space coordinates and powell_xi */
	x = (double *) calloc(6, sizeof(double));
	powell_xi = (double **) calloc(6, sizeof(double *));
	for (i = 0; i < 6; i++) {
		powell_xi[i] = (double *) calloc(6, sizeof(double));
	}
	/* for amoeba minimization */
	amoeba_p = (double **) calloc(7, sizeof(double *));
	for (i = 0; i < 7; i++) {
		amoeba_p[i] = (double *) calloc(6, sizeof(double));
	}
	amoeba_y = (double *) calloc(7, sizeof(double));
	/* basis vectors */
	n0  = (double *) calloc(3, sizeof(double));
	iota[0]  = (double *) calloc(3, sizeof(double));
	iota[1]  = (double *) calloc(3, sizeof(double));
	/* orbit elements */
	theorbit = (keplerorbit *) calloc(1, sizeof(keplerorbit));
	/* gradient and second derivatives */
	gx = (double **) calloc(3, sizeof(double *)) + 1;
	for (i = -1; i <= 1; i++) {
		gx[i] = (double *) calloc(6, sizeof(double));
	}
	chi_i = (double *) calloc(6, sizeof(double));
	chi_ij = (double **) calloc(6, sizeof(double *));
	for (i = 0; i < 6; i++) {
		chi_ij[i] = (double *) calloc(6, sizeof(double));
	}

	/* open the pipe for the initial solution */
	sprintf(lccom, "lc -b -o ra0 va0 re0 ve0 t0 nt < %s", initsolcatfilename);
	ipf = popen(lccom, "r");
	if (!ipf) {
		fprintf(stderr, "fitorbit : failed to open lc-pipe for input\n");
		exit(-1);
	}

	/* read the initial solution */
	fread(ra0, sizeof(double), 3, ipf);
	fread(va0, sizeof(double), 3, ipf);
	fread(re0, sizeof(double), 3, ipf);
	fread(ve0, sizeof(double), 3, ipf);
	fread(&t0, sizeof(double), 1, ipf);
	fread(&dtmp, sizeof(double), 1, ipf);
	nt = (int) dtmp;
	pclose(ipf);

	/* allocate observation data */
	t = (double *) calloc(nt, sizeof(double));
	sigma = (double *) calloc(nt, sizeof(double));
	re = (double **) calloc(nt, sizeof(double *));
	rho = (double **) calloc(nt, sizeof(double *));
	nobs = (double **) calloc(nt, sizeof(double *));
	nin = (double **) calloc(nt, sizeof(double *));
	for (it = 0; it < nt; it++) {
		re[it] = (double *) calloc(3, sizeof(double));
		rho[it] = (double *) calloc(3, sizeof(double));
		nobs[it] = (double *) calloc(3, sizeof(double));
		nin[it] = (double *) calloc(3, sizeof(double));
	}

	/* open the pipe for the observations */
	ipf = popen("lc -b -o t sigma re rho nobs", "r");
	if (!ipf) {
		fprintf(stderr, "fitorbit : failed to open lc-pipe for obscat\n");
		exit(-1);
	}

	/* read the observation data */
	for (it = 0; it < nt; it++) {
		fread(t + it, sizeof(double), 1, ipf);
		fread(sigma + it, sizeof(double), 1, ipf);
		fread(re[it], sizeof(double), 3, ipf);
		fread(rho[it], sizeof(double), 3, ipf);
		fread(nin[it], sizeof(double), 3, ipf);
	}
	pclose(ipf);

	/* compute basis vectors */
	copy(n0, nin[0]);
	getbasis(n0, iota[0], iota[1]);

	/* generate the header for output */
	system("lc -C -b -N '1 3 r' -N '1 3 v' -n chisquared -n p < /dev/null");

	/* output the Earth's coordinates */
	chi2min = 0.0;
	fwrite(re0, sizeof(double), 3, stdout);
	fwrite(ve0, sizeof(double), 3, stdout);
	fwrite(&chi2min, sizeof(double), 1, stdout);
	dtmp = -1.0;
	fwrite(&dtmp, sizeof(double), 1, stdout);

	for (ip = 0; ip <= np; ip++) {
		/* add the noise to nin to make nobs */
		for (it = 0; it < nt; it++) {
			for (i = 0; i < 3; i++) {
				nobs[it][i] = nin[it][i];
				if (ip) {
					nobs[it][i] += sigma[it] * gaussdev();
				}
			}
			scale(nobs[it], 1.0 / length(nobs[it]));
		}
		for (i = 0; i < 3; i++) {
			x[i] = ra0[i];
			x[3+i] = va0[i];
		}
		/* chi2min = chisquared(x); */
		/* compute the gradient and second derivatives
		getderivs(ra0, va0, 0.1 * sigma[0]); */
		/* find the minimum of chi^2 using powell's method */
		findmin_powell(ra0, va0, x);
		/* getderivs(x, x + 3, 1.e-10 * sigma[0]);
		if (ip == 1) {
			getchiplanes(x, 64, 2.e-2 * sigma[0]);
			exit(0);
		} */
		/* find the minimum of chi^2 using downhill simplex method
		findmin_amoeba(ra0, va0, x, 0.1 * sigma[0]); */
		chi2min = chisquared(x);
		/* output the result */
		fwrite(x, sizeof(double), 6, stdout);
		fwrite(&chi2min, sizeof(double), 1, stdout);
		dtmp = (double) ip;
		fwrite(&dtmp, sizeof(double), 1, stdout);
	}
	exit(0);
}

void	getbasis(double *n, double *i0, double *i1)
{
	double	theta, phi;

	theta = acos(n[2]);
	phi = atan2(n[1], n[0]);
        assign(i0, - sin(phi), cos(phi), 0.0);
        assign(i1, cos(theta) * cos(phi), cos(theta) * sin(phi), - sin(theta));
}

double	chisquared(double *x)
{
	int	i, it, ns;
	double	chi2, dn;

	chi2 = 0.0;
	for (it = 0; it < nt; it++) {
		/* copy the current phase-space coords */
		for (i = 0; i < 3; i++) {
			ra[i] = x[i];
			va[i] = x[3 + i];
		}
		if (dotcl) {
			ns = (int) ceil(fabs(t[it] - t0) / dtmax);
			if (ns) {
				/* evolve the phase-space coordinates */
				tcl((t[it] - t0) / ns, ns, ra, va);
			}
		} else {
			cartesiantokepler(ra, va, theorbit);
			theorbit->M += (t[it] - t0) / pow(theorbit->a, 1.5);
			keplertocartesian(theorbit, ra, va);
		}
		for (i = 0; i < 3; i++) {
			n[i] = ra[i] - re[it][i] - rho[it][i];
		}
		scale(n, 1.0 / length(n));
		for (i = 0; i < 3; i++) {
			dn = nobs[it][i] - n[i];
			chi2 += dn * dn / (sigma[it] * sigma[it]);
		}
	}
	return(chi2);
}

void	findmin_powell(double *r0, double *v0, double *x)
{
	int	powell_iter, i, j;

	/* create the directions and initial position for powell */
	for (i = 0; i < 6; i++) {
		for (j = 0; j < 6; j++) {
			powell_xi[i][j] = 0.0;
		}
	}
	for (i = 0; i < 3; i++) {
		powell_xi[0][i] = n0[i];
		powell_xi[1][i] = iota[0][i];
		powell_xi[2][i] = iota[1][i];
		powell_xi[3][3 + i] = n0[i];
		powell_xi[4][3 + i] = iota[0][i];
		powell_xi[5][3 + i] = iota[1][i];
		x[i] = r0[i];
		x[3 + i] = v0[i];
	}
	
	/* run powell to find the minimum chi^2 */
	powell_ftol = 1.e-15;
	powell(x, powell_xi, 6, powell_ftol, &powell_iter, &powell_fret, chisquared);
}

void	findmin_amoeba(double *r0, double *v0, double *x, double dx)
{
	int	nfunccalls, i, j;

	/* create the simplex */
	for (j = 0; j < 7; j++) {
		for (i = 0; i < 3; i++) {
			amoeba_p[j][i] = r0[i];
			amoeba_p[j][i+3] = v0[i];
		}
		if (j) {
			amoeba_p[j][j-1] += dx;
		}
		amoeba_y[j] = chisquared(amoeba_p[j]);
	}

	amoeba_ftol = 1.e-8;
	amoeba(amoeba_p, amoeba_y, 6, amoeba_ftol, chisquared, &nfunccalls);
	fprintf(stderr, "# number of calls: %d\n", nfunccalls);
	for (i = 0; i < 6; i++) {
		x[i] = amoeba_p[0][i];
	}
}

void	getderivs(double *r0, double *v0, double dx)
{
	int	i, j, k;

	for (i = 0; i < 3; i++) {
		gx[0][i] = r0[i];
		gx[0][3+i] = v0[i];
	}
	for (j = 0; j < 6; j++) {
		for (k = 0; k < 6; k++) {
			gx[-1][k] = gx[1][k] = gx[0][k];
		}
		gx[-1][j] -= dx;
		gx[+1][j] += dx;
		fprintf(stderr, "# %14.8le %14.8le %14.8le\n", 
			chisquared(gx[0]), chisquared(gx[-1]) - chisquared(gx[0]), chisquared(gx[+1]) - chisquared(gx[0]));
	}
}

void	getchiplanes(double *x, int N, double dx)
{
	int 	i, j, k, I, J, N2, nplanes;
	float	*f;
	FILE	*opf;
	char	com[1024];
	double	*xtmp;

	xtmp = (double *) calloc(6, sizeof(double));
	f = (float *) calloc(N * N, sizeof(float));

	sprintf(com, "imhead -g -32 4 %d %d 5 6 > tmp.fits", N, N);
	system(com);
	opf = fopen("tmp.fits", "a");

	N2 = N / 2;
	for (I = 0; I < 6; I++) {
		for (J = 0; J < 6; J++) {
			if (J != I) {
				for (k = 0; k < 6; k++) {
					xtmp[k] = x[k];
				}
				for (i = 0; i < N; i++) {
					for (j = 0; j < N; j++) {
						xtmp[I] = x[I] + (i - N2) * dx;
						xtmp[J] = x[J] + (j - N2) * dx;
						f[i * N + j] = (float) chisquared(xtmp);
					}
				}
				fwrite(f, sizeof(float), N * N, opf);
			}
		}
	}
	fclose(opf);

	free(f);
}
