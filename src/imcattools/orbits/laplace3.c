#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vectors.h"
#include "Ffunc.h"
#include "gaussdev.h"
#include "makedddotr.h"

#define usage "\nNAME\n\
	laplace3 - compute an approximate orbit from 3 observations\n\
\n\
SYNOPSIS\n\
	laplace3\n\
\n\
DESCRIPTION\n\
	laplace3 reads an lc format catalog (probably generated by makeobs)\n\
	from stdin.  This calatog must contain 3 lines, each containing\n\
	at least:\n\
		t		time\n\
		sigma	astrometry precision in radians\n\
		re[]	earth's position\n\
		ve[]	earth's velocity\n\
		rho[]	position of observatory vs earth\n\
		ra[]	asteroid's position\n\
		va[]	asteroid's velocity\n\
\n\
	It then solves for the phase space coordinates of the asteroid\n\
	at t=0.  Currently it does this graphically.\n\
\n\
	Times are given in units of earths dynamical time (approx\n\
	57.21 days in this model).\n\
\n\
	Distances are expressed in AU.\n\
\n\
SEE ALSO\n\
	maketestpscoords.pl makeobs\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@hawaii.edu\n\n"

#define	DTINY	1.0e-20
#define FTOL	1.0e-8
#define NITERMAX	100

int	newtonraphson(double *d, double Fobs, double mu);

main (int argc, char *argv[])
{
	/* input values */
	double	*t, **re, **ve, **rho, **nobs, **ncor, **ra_in, **va_in, *sigma;
	double	*n, *dotn, *ddotn, *iota[2], *ra, *va, *dr;
	double	dt, ddotn0, ddotn1, re0, re1, Fobs, d, mu, rre, r, ndot, vpara;
	double	dmax, dmin, dlogd, ndotrho, lastF, newF, lastd, newd, dFdd, sigmad;
	int	i, it, nt, niter, firstd, signF;
	FILE	*ipf, *opf;

	/* parse args */
	if (argc != 1) {
		fprintf(stderr, usage);
		exit(-1);
	}

	/* number of observations */
	nt = 3;

	/* allocate input data */
	t = (double *) calloc(nt, sizeof(double));
	sigma = (double *) calloc(nt, sizeof(double));
	re = (double **) calloc(nt, sizeof(double *));
	ve = (double **) calloc(nt, sizeof(double *));
	rho = (double **) calloc(nt, sizeof(double *));
	nobs = (double **) calloc(nt, sizeof(double *));
	ra_in = (double **) calloc(nt, sizeof(double *));
	va_in = (double **) calloc(nt, sizeof(double *));
	for (it = 0; it < nt; it++) {
		re[it] = (double *) calloc(nt, sizeof(double));
		ve[it] = (double *) calloc(nt, sizeof(double));
		rho[it] = (double *) calloc(nt, sizeof(double));
		nobs[it] = (double *) calloc(nt, sizeof(double));
		ra_in[it] = (double *) calloc(nt, sizeof(double));
		va_in[it] = (double *) calloc(nt, sizeof(double));
	}
	/* allocate internal arrays */
	ncor = (double **) calloc(nt, sizeof(double *));
	for (it = 0; it < nt; it++) {
		ncor[it] = (double *) calloc(nt, sizeof(double));
	}
	n = (double *) calloc(3, sizeof(double));
	dotn = (double *) calloc(3, sizeof(double));
	ddotn = (double *) calloc(3, sizeof(double));
	iota[0] = (double *) calloc(3, sizeof(double));
	iota[1] = (double *) calloc(3, sizeof(double));
	ra = (double *) calloc(3, sizeof(double));
	va = (double *) calloc(3, sizeof(double));

	/* open the pipe for the observations */
	ipf = popen("lc -b -o t sigma re ve ra va rho n", "r");
	if (!ipf) {
		fprintf(stderr, "laplace3 : failed to open lc-pipe for input\n");
		exit(-1);
	}

	/* open lc-pipe for output */
	opf = popen("lc -C -n d -n mu -n F -n dFdd -n re1 -n sigmad -N '1 3 ra' -N '1 3 va' -N '1 3 ra_true' -N '1 3 va_true' -N '1 3 re' -N '1 3 ve'", "w");
	if (!opf) {
		fprintf(stderr, "laplace3 : failed to open lc-pipe for output\n");
		exit(-1);
	}

	while (1) {
		/* read the observations */
		for (it = 0; it < nt; it++) {
			if (!fread(t + it, sizeof(double), 1, ipf)) {
				pclose(ipf);
				pclose(opf);
				exit(0);
			}
			fread(sigma + it, sizeof(double), 1, ipf);
			fread(re[it], sizeof(double), 3, ipf);
			fread(ve[it], sizeof(double), 3, ipf);
			fread(ra_in[it], sizeof(double), 3, ipf);	
			fread(va_in[it], sizeof(double), 3, ipf);
			fread(rho[it], sizeof(double), 3, ipf);
			fread(nobs[it], sizeof(double), 3, ipf);
		}

		/* this is |re[]| */
		rre = length(re[1]);

		/* this version solves graphically by simply stepping down in d-values */
		dmax = 10.0;
		dmin = 0.1;
		dlogd = 0.02;
		firstd = 1;

		for (d = dmax; d > dmin; d *= (1 - dlogd)) {
			/* the physical distance */
			r = rre * d;
	
			/* compute corrected directions */
			for (it = 0; it < nt; it++) {
				ndotrho = dot(nobs[it], rho[it]);
				for (i = 0; i < 3; i++) {
					ncor[it][i] = nobs[it][i] + (rho[it][i] - ndotrho * nobs[it][i]) / r;
				}
			}
	
			/* get estimates of n[], dotn[], ddotn[] */	
			/* we should solve A N = B here, but for now we use the solution for nt=3 equal interval case */
			dt = t[2] - t[1];
			copy(n, ncor[1]);
			diff(dotn, ncor[2], ncor[0]);
			scale(dotn, 0.5 / dt);
			/* compute ddotn */
			for (i = 0; i < 3; i++) {
				ddotn[i] = (ncor[2][i] - 2 * ncor[1][i] + ncor[0][i]) / (dt * dt);
			}
	
			/* compute mu */
			mu = dot(re[1], n) / length(re[1]);

			/* create the basis vectors */
			copy(iota[0], dotn);
			scale(iota[0], 1.0 / length(iota[0]));
			iota[1][0] = n[1] * iota[0][2] - n[2] * iota[0][1];
			iota[1][1] = n[2] * iota[0][0] - n[0] * iota[0][2];
			iota[1][2] = n[0] * iota[0][1] - n[1] * iota[0][0];

			/* compute ddotn[01], re[01] */
			ddotn0 = dot(ddotn, iota[0]);
			ddotn1 = dot(ddotn, iota[1]);
			re0 = dot(re[1], iota[0]);
			re1 = dot(re[1], iota[1]);

			/* compute Fobs */
			if (fabs(re1) > DTINY) {
				Fobs = ddotn1 * pow(rre, 4.0) / re1;
			} else {
				fprintf(stderr, "laplace3 : re1 is too small - bailing\n");
				exit(-1);
			}
fprintf(stderr, "Fobs = %lf\n", Fobs);
			if (firstd) {
				signF = (F(d, mu) > Fobs ? 1 : -1);
				firstd = 0;
			} else {
				if (signF * (F(d, mu) - Fobs) < 0.0) {		/* we have passed a solution */
					signF *= -1;
					newd = d;
					newF = F(d, mu) - Fobs;
					dFdd = (newF - lastF) / (newd - lastd);
					if (fabs(dFdd) > DTINY) {
						d = lastd - lastF / dFdd;
					}
					/* compute the uncertainty in distance */
					if (re1 > DTINY) {
						sigmad = sigma[1] * sqrt(6.0) / (fabs(dFdd) * re1 * dt * dt);
					} else {
						sigmad = 0.0;
					}
					/* compute the position and velocity */
					r = d * rre;
					ndot = length(dotn);
					vpara = (0.5 * r / ndot) * (re0 * ddotn1 / re1 - ddotn0);
					for (i = 0; i < 3; i++) {
						ra[i] = re[1][i] + r * n[i];
						va[i] = ve[1][i] + r * dotn[i] + vpara * n[i];
					}
					fprintf(opf, "%14.8lg %14.8lg %14.8lg %14.8lg %14.8lg %14.8lg\n", d, mu, F(d, mu), dFdd, re1, sigmad);
					fprintvec(ra, opf);
					fprintvec(va, opf);
					fprintvec(ra_in[1], opf);
					fprintvec(va_in[1], opf);
					fprintvec(re[1], opf);
					fprintvec(ve[1], opf);
				}
			}
			lastF = F(d, mu) - Fobs;
			lastd = d;
		}
	}
}

int	newtonraphson(double *d, double Fobs, double mu)
{
	int	niter;
	double	dl;

	/* solve for d */
	fprintf(stderr, "# solving for d...\n");
	dl = 3.0;	/* initial d */
	niter = 0;
	while (fabs(Fobs - F(dl, mu)) > FTOL) { 
		dl += (Fobs - F(dl, mu)) / dFdd(dl, mu);
		niter++;
		if (niter > NITERMAX) {
			fprintf(stderr, "orbit : niter > NITERMAX\n");
			exit(-1);
		}
	}
	fprintf(stderr, "# d=%14.8lf\n# F=%14.8lf\n# niter = %d\n", dl, F(dl, mu), niter);
	*d = dl;
}

