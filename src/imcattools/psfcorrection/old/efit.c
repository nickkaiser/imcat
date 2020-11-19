/*
 * efit.c
 */

#define	usage "\n\n\n\
NAME\n\
		efit --- generate model for stellar ellipticities\n\
SYNOPSIS\n\
		efit	[option...] a.cat b.cat ....\n\
			-m emax		# reject stars with |e| > emax (0.4)\n\
			-x numax	# don't use stars with |eres| > numax * sigma (3.0)\n\
			-n N		# image size (2048)\n\
			-o order	# output this order model to stdout (1)\n\
			-r		# output e's and residuals\n\
			-e ename	# name for ellipticity (e)\n\
\n\
DESCRIPTION\n\
		\"efit\" fits the p = psm^{-1} e values of stars from a list of cat files\n\
		to a model of constant offset for each frame plus a global 1st or\n\
		2nd order Taylor expansion\n\
		Does least squares 1st order fit, then refits ignoring outliers > numax * sigma\n\
		finally refits 2nd order Taylor expansion\n\
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


double		**x[2];
double		**e[2], **eres[2], ereslimit;
int		nframes = 0, *nobj, nmodes, noutliers;
double		**A, *B[2], *a[2];
void		fillarrays(void);
void		makeresiduals(double *sigmae);
double		modelrms(int mode1, int mode2);
void		writeamplitudes(int order, int framesize);
void		printoffsets(void);
void		printgrads(void);
void		printgradgrads(void);


main(int argc, char *argv[])	
{
	int		arg = 1; 
	FILE		*catf;
	int		comc, N, iobj, frame, pol, mode, l, m, *indx, order;
	int		ntotal = 0, nrejected = 0;
	char		*comv[MAX_COMMENTS], lcstring[128], *ename, defename[2] = "e", line[1024];
	double		emax, eesum, nsum, sigmae, numax;
	double		d;
	int		outputresiduals = 0;
	double		E0, E1, x0, x1, obje[2], psm[2][2], objx[2];
	
	/* defaults */
	emax = 0.4;
	N = 2048;
	numax = 3.0;
	order = 1;
	ename = defename;
	
	while (arg < argc) {
		if (*argv[arg] == '-') {
			switch (*(argv[arg++]+1)) {
				case 'm':
					if (EOF == sscanf(argv[arg++], "%lf", &emax))
						error_exit(usage);
					break;
				case 'n':
					if (EOF == sscanf(argv[arg++], "%d", &N))
						error_exit(usage);
					break;
				case 'o':
					if (EOF == sscanf(argv[arg++], "%d", &order))
						error_exit(usage);
					break;
				case 'x':
					if (EOF == sscanf(argv[arg++], "%lf", &numax))
						error_exit(usage);
					break;
				case 'r':
					outputresiduals = 1;
					break;
				case 'e':
					ename = argv[arg++];
					break;
				default:
					error_exit(usage);
					break;
			}
		} else {
			nframes++;
			arg++;
		}
	}
	
	nobj = (int *) calloc(nframes, sizeof(int));
	x[0] = (double **) calloc(nframes, sizeof(double *));
	x[1] = (double **) calloc(nframes, sizeof(double *));
	e[0] = (double **) calloc(nframes, sizeof(double *));
	e[1] = (double **) calloc(nframes, sizeof(double *));
	eres[0] = (double **) calloc(nframes, sizeof(double *));
	eres[1] = (double **) calloc(nframes, sizeof(double *));
	setnframes(nframes);
	setframesize(N);
	
	/* count em and allocate space for x, y, e[0], e[1] */
	for (frame = 0; frame < nframes; frame++) {
		arg = argc - nframes + frame;
		sprintf(lcstring, "lc -o x %s psm < %s", ename, argv[arg]);
		if (!(catf = popen(lcstring, "r")))
			error_exit("efit: failed to open lc-pipe for input\n");
		while(fgets(line, 1024, catf)) {
			sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf", &(objx[0]),  &(objx[1]), 
				&(obje[0]),  &(obje[1]), &(psm[0][0]), &(psm[0][1]), &(psm[1][0]), &(psm[1][1]));
			if ((obje[0] == 0.0 && obje[1] == 0.0) || obje[0] * obje[0] + obje[1] * obje[1] > emax * emax) {
				nrejected++;
				continue;
			}
			nobj[frame]++;
			ntotal++;
		}
		pclose(catf);
		x[0][frame] = (double *) calloc(nobj[frame], sizeof(double));
		x[1][frame] = (double *) calloc(nobj[frame], sizeof(double));
		e[0][frame] = (double *) calloc(nobj[frame], sizeof(double));
		e[1][frame] = (double *) calloc(nobj[frame], sizeof(double));
		eres[0][frame] = (double *) calloc(nobj[frame], sizeof(double));
		eres[1][frame] = (double *) calloc(nobj[frame], sizeof(double));
	}
	fprintf(stderr, "# %d objects rejected\n", nrejected);
	fprintf(stderr, "# %d objects read\n", ntotal);
	/* read em */
	eesum = nsum = 0.0;
	for (frame = 0; frame < nframes; frame++) {
		iobj = 0;
		arg = argc - nframes + frame;
		sprintf(lcstring, "lc -o x %s psm < %s", ename, argv[arg]);
		if (!(catf = popen(lcstring, "r")))
			error_exit("efit: failed to open lc-pipe for input\n");
		while(fgets(line, 1024, catf)) {
			sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf", &(objx[0]),  &(objx[1]), 
				&(obje[0]),  &(obje[1]), &(psm[0][0]), &(psm[0][1]), &(psm[1][0]), &(psm[1][1]));
			if ((obje[0] == 0.0 && obje[1] == 0.0) || obje[0] * obje[0] + obje[1] * obje[1] > emax * emax) {
				continue;
			}
			x[0][frame][iobj] = objx[0];
			x[1][frame][iobj] = objx[1];
			e[0][frame][iobj] = 2 * obje[0] / (psm[0][0] + psm[1][1]);
			e[1][frame][iobj] = 2 * obje[1] / (psm[0][0] + psm[1][1]);
			eesum += e[0][frame][iobj] * e[0][frame][iobj] + e[1][frame][iobj] * e[1][frame][iobj];
			nsum += 2.0;
			iobj++;
		}
		pclose(catf);
	}
	fprintf(stderr, "#      raw rms e / P = %.3f\n", sigmae = sqrt(eesum / nsum));
	
	/* allocate arrays */
	nmodes = nframes + 5;
	indx = (int *) calloc(nmodes, sizeof(int));
	a[0] = (double *) calloc(nmodes, sizeof(double));
	a[1] = (double *) calloc(nmodes, sizeof(double));
	B[0] = (double *) calloc(nmodes, sizeof(double));
	B[1] = (double *) calloc(nmodes, sizeof(double));
	A = (double **) calloc(nmodes, sizeof(double *));
	for (mode = 0; mode < nmodes; mode++)
		A[mode] = (double *) calloc(nmodes, sizeof(double));
	
	/* do 1st order fit for all objects */
	fprintf(stderr, "# fitting 1st order model\n");
	nmodes = nframes + 2;
	ereslimit = 100.0;
	fillarrays();
	myludcmp(A, nmodes, indx, &d);
	mylubksb(A, nmodes, indx, B[0]);
	for (l = 0; l < nmodes; l++)
		a[0][l] = B[0][l];
	mylubksb(A, nmodes, indx, B[1]);
	for (l = 0; l < nmodes; l++)
		a[1][l] = B[1][l];
	makeresiduals(&sigmae);
	fprintf(stderr, "# residual rms e / P (per component) = %.4f\n", sigmae);
	
	/* now redo the fit rejecting outliers */
	fprintf(stderr, "#\n#\n# refitting rejecting outliers with e_res > %.1f x sigma....\n", numax);
	ereslimit = numax * sigmae;
	fillarrays();
	fprintf(stderr, "# %d outliers ignored\n", noutliers);
	myludcmp(A, nmodes, indx, &d);
	mylubksb(A, nmodes, indx, B[0]);
	for (l = 0; l < nmodes; l++)
		a[0][l] = B[0][l];
	mylubksb(A, nmodes, indx, B[1]);
	for (l = 0; l < nmodes; l++)
		a[1][l] = B[1][l];
	if (order == 1 && !outputresiduals) {
		writeamplitudes(order, N);
	}
	makeresiduals(&sigmae);
	fprintf(stderr, "# residual rms e / P (per component) = %.4f\n", sigmae);
	fprintf(stderr, "# rms model contributions:\n");
	fprintf(stderr, "# total model rms = %.4f\n", modelrms(0, nmodes - 1));
	fprintf(stderr, "#         offsets = %.4f\n", modelrms(0, nframes - 1));
	fprintf(stderr, "#          linear = %.4f\n", modelrms(nframes, nframes + 1));
	
	/* now redo the fit for 2nd order taylor expansion */
	nmodes = nframes + 5;
	fprintf(stderr, "#\n#\n# refitting with 2nd derivatives\n");
	ereslimit = numax * sigmae;
	fillarrays();
	fprintf(stderr, "# %d outliers ignored\n", noutliers);
	myludcmp(A, nmodes, indx, &d);
	mylubksb(A, nmodes, indx, B[0]);
	for (l = 0; l < nmodes; l++)
		a[0][l] = B[0][l];
	mylubksb(A, nmodes, indx, B[1]);
	for (l = 0; l < nmodes; l++)
		a[1][l] = B[1][l];
	if (order == 2 && !outputresiduals) {
		writeamplitudes(order, N);
	}
	makeresiduals(&sigmae);
	fprintf(stderr, "# residual rms ellipticity (per component) = %.4f\n", sigmae);
	
	/* now we calculate the model rms averaged over galaxies */
	fprintf(stderr, "# rms model contributions:\n");
	fprintf(stderr, "# total model rms = %.4f\n", modelrms(0, nmodes - 1));
	fprintf(stderr, "#         offsets = %.4f\n", modelrms(0, nframes - 1));
	fprintf(stderr, "#          linear = %.4f\n", modelrms(nframes, nframes + 1));
	fprintf(stderr, "#       quadratic = %.4f\n", modelrms(nframes + 2, nframes + 4));


	/* the idea here is to output residual after subtracting the frame offsets */
	/* with position along gradient direction in order to show gradient clearly */	
	if (outputresiduals) {
		fprintf(stdout, "# the idea here is to output residual after subtracting the frame offsets\n");
		fprintf(stdout, "# with position along gradient direction in order to show gradient clearly\n");
		fprintf(stdout, "#       x0         E0         x1         E1      eres0      eres5\n");
		for (frame = 0; frame < nframes; frame++) {
			for (iobj = 0; iobj < nobj[frame]; iobj++) {
				E0 = e[0][frame][iobj] - a[0][frame];
				E1 = e[1][frame][iobj] - a[1][frame];
				x0 = (x[0][frame][iobj] - N / 2) * a[0][nframes] + (x[1][frame][iobj] - N / 2) * a[0][nframes + 1];
				x0 /= sqrt(a[0][nframes] * a[0][nframes] + a[0][nframes + 1] * a[0][nframes + 1]);
				x1 = (x[0][frame][iobj] - N / 2) * a[1][nframes] + (x[1][frame][iobj] - N / 2) * a[1][nframes + 1];
				x1 /= sqrt(a[1][nframes] * a[1][nframes] + a[1][nframes + 1] * a[1][nframes + 1]);
				fprintf(stdout, "%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n",
					x0, E0, x1, E1, eres[0][frame][iobj], eres[1][frame][iobj]);
			}
		}
	}
	exit(0);
}



void		fillarrays(void)
{
	int	pol, im, jm, frame, iobj;
	double	ee;
	
	for (im = 0; im < nmodes; im++) {
		for (pol = 0; pol < 2; pol++)
			B[pol][im] = 0.0;
		for (jm = 0; jm < nmodes; jm++)
			A[im][jm] = 0.0;
	}
	
	noutliers = 0;
	for (frame = 0; frame < nframes; frame++) {
		for (iobj = 0; iobj < nobj[frame]; iobj++) {
			if (fabs(e[0][frame][iobj]) > ereslimit || fabs(e[1][frame][iobj]) > ereslimit) {
				noutliers++;
				continue;
			}
			for (im = 0; im < nmodes; im++) {
				for (pol = 0; pol < 2; pol++)
					B[pol][im] += e[pol][frame][iobj] * 
						g(im, frame, x[0][frame][iobj], x[1][frame][iobj]);
				for (jm = 0; jm < nmodes; jm++) {
					A[im][jm] += g(im, frame, x[0][frame][iobj], x[1][frame][iobj]) *
						g(jm, frame, x[0][frame][iobj], x[1][frame][iobj]);
				}
			}
		}
	}
}



void	makeresiduals(double *sigmae)
{
	double	eesum, nsum;
	int	frame, iobj, pol, m;

	eesum = nsum = 0.0;
	for (frame = 0; frame < nframes; frame++) {
		for (iobj = 0; iobj < nobj[frame]; iobj++) {
			for (pol = 0; pol < 2; pol++) {
				eres[pol][frame][iobj] = e[pol][frame][iobj];
				for (m = 0; m < nmodes; m++)
					eres[pol][frame][iobj] -= a[pol][m] * g(m, frame, x[0][frame][iobj], x[1][frame][iobj]);
				if (fabs(eres[pol][frame][iobj]) < ereslimit) {
					eesum += eres[pol][frame][iobj] * eres[pol][frame][iobj];
					nsum += 1.0;
				}
			}
		}
	}
	*sigmae = sqrt(eesum / nsum);
}



double		modelrms(int mode1, int mode2)
/* calculate rms for modes mode1 to mode2 inclisive */
{
	double	emodel[2], eesum = 0.0, nsum = 0.0;
	int	frame, pol, mode, iobj;
	
	for (frame = 0; frame < nframes; frame++) {
		for (iobj = 0; iobj < nobj[frame]; iobj++) {
			for (pol = 0; pol < 2; pol++) {
				emodel[pol] = 0.0;
				for (mode = mode1; mode <= mode2; mode++) {
					emodel[pol] += a[pol][mode] * g(mode, frame, 
						x[0][frame][iobj], x[1][frame][iobj]);
				}
				eesum += emodel[pol] * emodel[pol];
				nsum += 1.0;
			}
		}
	}
	return (sqrt(eesum / nsum));
}


void		writeamplitudes(int order, int framesize)
{
		fprintf(stdout, "%d modes\n", nmodes);
		fprintf(stdout, "%d order model\n", order);
		fprintf(stdout, "%d frames\n", nframes);
		fprintf(stdout, "%d pixels wide\n", framesize);
		switch (order) {
			case 1:
				printoffsets();
				printgrads();
				break;
			case 2:
				printoffsets();
				printgrads();
				printgradgrads();
				break;
			default:
				error_exit("writeamplitudes: bad order\n");
				break;
		}
}


void	printoffsets(void)
{
	int	l;
	
	fprintf(stdout, "# offsets:\n");
	fprintf(stdout, "# mode       e[0]       e[1]\n");
	for (l = 0; l < nframes; l++)
		fprintf(stdout, "%6d %10.3e %10.3e\n", l, a[0][l], a[1][l]);
}



void	printgrads(void)
{
	int	l;
	
	fprintf(stdout, "# gradients:\n");
	fprintf(stdout, "# mode       g[0]       g[1]\n");
	for (l = nframes; l < nframes + 2; l++)
		fprintf(stdout, "%6d %10.3e %10.3e\n", l, a[0][l], a[1][l]);
}	



void	printgradgrads(void)
{
	int	l;
	
	fprintf(stdout, "# 2nd derivatives:\n");
	fprintf(stdout, "# mode       g[0]       g[1]\n");
	for (l = nframes + 2; l < nframes + 5; l++)
		fprintf(stdout, "%6d %10.3e %10.3e\n", l, a[0][l], a[1][l]);
}	



