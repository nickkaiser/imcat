/*
 * efit.c
 */

#define	usage "\n\n\n\
NAME\n\
	efit --- generate model for stellar ellipticities\n\
SYNOPSIS\n\
	efit [option...]\n\
		-m emax		# reject stars with |e| > emax (0.4)\n\
		-x numax	# don't use stars with |eres| > numax * sigma (3.0)\n\
		-n N		# image size (2048)\n\
		-o order	# order of Taylor series (1)\n\
		-e ename	# name for ellipticity (e)\n\
		-i niter	# number of times to iterate (2)\n\
\n\
DESCRIPTION\n\
	\"efit\" reads a catalogue of stars from stdin and writes a set\n\
	of coefficients for Taylor series expansion of the e / psm values\n\
	to stdout.\n\
	We first reject stars with |e| > emax, make the fit, and calculate\n\
	the rms deviation from the fit. We then reject outliers\n\
	with |eres| > numax * sigma and then we refit to get refined\n\
	coefficients.  Use '-i' option to iterate rejection-fitting\n\
	cycle. Maximum order is 6.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@cita.utoronto.ca\n\
\n\n\n"		


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "../../utils/error.h"
#include "efit_stuff.h"
#include "../../utils/lu.h"
#include "../../imlib/fits.h"

#define MAX_STARS 10000


double		x[2][MAX_STARS];
double		e[2][MAX_STARS], eres[2][MAX_STARS], ereslimit;
int		nobj, nmodes, noutliers;
int		maxorder = 6, Nmodes[7] = {1, 3, 6, 10, 15, 21, 28};
double		**A, *B[2], *a[2];

void		fillarrays(void);
void		makeresiduals(double *sigmae);


main(int argc, char *argv[])	
{
	int		arg = 1; 
	int		comc, N, iobj, pol, mode, l, m, *indx, order, niter;
	int		ntotal = 0, nrejected = 0;
	char		lcstring[128], *ename, defename[2] = "e", line[1024];
	double		emax, eesum, nsum, sigmae, numax;
	double		d;
	double		E0, E1, x0, x1, obje[2], psm[2][2], objx[2];
	FILE		*lcpipe;
	
	/* defaults */
	emax = 0.4;
	N = 2048;
	numax = 3.0;
	order = 1;
	ename = defename;
	niter = 2;
	
	while (arg < argc) {
		if (*argv[arg] == '-') {
			switch (*(argv[arg++]+1)) {
				case 'm':
					if (1 != sscanf(argv[arg++], "%lf", &emax))
						error_exit(usage);
					break;
				case 'n':
					if (1 != sscanf(argv[arg++], "%d", &N))
						error_exit(usage);
					break;
				case 'o':
					if (1 != sscanf(argv[arg++], "%d", &order))
						error_exit(usage);
					if (order > maxorder)
						error_exit("efit: order too high!\n");
					break;
				case 'x':
					if (1 != sscanf(argv[arg++], "%lf", &numax))
						error_exit(usage);
					break;
				case 'e':
					ename = argv[arg++];
					break;
				case 'i':
					if (1 != sscanf(argv[arg++], "%d", &niter))
						error_exit(usage);
					break;
				default:
					error_exit(usage);
					break;
			}
		} else {
			error_exit(usage);
		}
	}
	
	setframesize(N);
	
	/* read em */
	eesum = nsum = 0.0;
	iobj = 0;
		sprintf(lcstring, "lc -o x %s psm", ename);
	if (!(lcpipe = popen(lcstring, "r")))
		error_exit("efit: failed to open lc-pipe for input\n");
	while(fgets(line, 1024, lcpipe)) {
		sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf", &(objx[0]),  &(objx[1]), 
				&(obje[0]),  &(obje[1]), &(psm[0][0]), &(psm[0][1]), &(psm[1][0]), &(psm[1][1]));
		if ((obje[0] == 0.0 && obje[1] == 0.0) || obje[0] * obje[0] + obje[1] * obje[1] > emax * emax) {
			nrejected++;
			continue;
		}
		x[0][iobj] = objx[0];
		x[1][iobj] = objx[1];
		e[0][iobj] = 2 * obje[0] / (psm[0][0] + psm[1][1]);
		e[1][iobj] = 2 * obje[1] / (psm[0][0] + psm[1][1]);
		eesum += e[0][iobj] * e[0][iobj] + e[1][iobj] * e[1][iobj];
		nsum += 2.0;
		iobj++;
	}
	pclose(lcpipe);
	nobj = iobj;
	fprintf(stderr, "# %d objects accepted\n# %d objects rejected\n", nobj, nrejected);
	fprintf(stderr, "#      raw rms e / P = %13.8lg\n", sigmae = sqrt(eesum / nsum));
	
	/* allocate arrays */
	nmodes = Nmodes[order];
	indx = (int *) calloc(nmodes, sizeof(int));
	a[0] = (double *) calloc(nmodes, sizeof(double));
	a[1] = (double *) calloc(nmodes, sizeof(double));
	B[0] = (double *) calloc(nmodes, sizeof(double));
	B[1] = (double *) calloc(nmodes, sizeof(double));
	A = (double **) calloc(nmodes, sizeof(double *));
	for (mode = 0; mode < nmodes; mode++)
		A[mode] = (double *) calloc(nmodes, sizeof(double));
	
	/* do fit for all objects */
	fprintf(stderr, "# fitting %d order model (%d coeffs) to all objects...\n", order, nmodes);
	ereslimit = 100.0;
	fillarrays();
	myludcmp(A, nmodes, indx, &d);
	mylubksb(A, nmodes, indx, B[0]);
	for (l = 0; l < nmodes; l++)
		a[0][l] = B[0][l];
	mylubksb(A, nmodes, indx, B[1]);
	for (l = 0; l < nmodes; l++)
		a[1][l] = B[1][l];

	while (niter > 0) {
		makeresiduals(&sigmae);
		fprintf(stderr, "# residual rms e / P (per component) = %13.8lg\n", sigmae);
		ereslimit = numax * sigmae;
		fprintf(stderr, "# refitting rejecting outliers with e_res > %lg\n", ereslimit);
		fillarrays();
		fprintf(stderr, "# %d outliers ignored\n", noutliers);
		myludcmp(A, nmodes, indx, &d);
		mylubksb(A, nmodes, indx, B[0]);
		for (l = 0; l < nmodes; l++)
			a[0][l] = B[0][l];
		mylubksb(A, nmodes, indx, B[1]);
		for (l = 0; l < nmodes; l++)
			a[1][l] = B[1][l];
		niter--;
	}
	makeresiduals(&sigmae);
	fprintf(stderr, "# final residual rms e / P (per component) = %13.8lg\n", sigmae);
	writeamplitudes(nmodes, order, N, a);
	exit(0);
}



void		fillarrays(void)
{
	int	pol, im, jm, iobj;
	double	ee;
	
	for (im = 0; im < nmodes; im++) {
		for (pol = 0; pol < 2; pol++)
			B[pol][im] = 0.0;
		for (jm = 0; jm < nmodes; jm++)
			A[im][jm] = 0.0;
	}
	
	noutliers = 0;
	for (iobj = 0; iobj < nobj; iobj++) {
		if (fabs(eres[0][iobj]) > ereslimit || fabs(eres[1][iobj]) > ereslimit) {
			noutliers++;
			continue;
		}
		for (im = 0; im < nmodes; im++) {
			for (pol = 0; pol < 2; pol++)
				B[pol][im] += e[pol][iobj] * g(im, x[0][iobj], x[1][iobj]);
			for (jm = 0; jm < nmodes; jm++) {
				A[im][jm] += g(im, x[0][iobj], x[1][iobj]) * g(jm, x[0][iobj], x[1][iobj]);
			}
		}
	}
}



void	makeresiduals(double *sigmae)
{
	double	eesum, nsum;
	int	iobj, pol, m;

	eesum = nsum = 0.0;
	for (iobj = 0; iobj < nobj; iobj++) {
		for (pol = 0; pol < 2; pol++) {
			eres[pol][iobj] = e[pol][iobj];
			for (m = 0; m < nmodes; m++)
				eres[pol][iobj] -= a[pol][m] * g(m, x[0][iobj], x[1][iobj]);
			if (fabs(eres[pol][iobj]) < ereslimit) {
				eesum += eres[pol][iobj] * eres[pol][iobj];
				nsum += 1.0;
			}
		}
	}
	*sigmae = sqrt(eesum / nsum);
}









