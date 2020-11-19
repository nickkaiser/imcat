#define usage "\n\
NAME\n\
	fitstack --- fit for transformation coefficients for a set of exposures\n\
\n\
SYNOPSIS\n\
	fitstack nexp lmax\n\
\n\
DESCRIPTION\n\
	'fitstack' reads from stdin a catalogue containing the result of\n\
	merging all pairs of cats for a stack of 'nexp' images (as created\n\
	by 'mergestacks1' or 'mergestacks2') and which must contain entries for\n\
	spatial coords 'x', magnitude 'mag' and exposure number 'exp'.\n\
	It then fits a model in which sky coords (in frame defined by exposure-0) are\n\
		r = r_e + dphi_e r_e + d_e\n\
	where the 2x2 matrix dphi allows for rotations between\n\
	exposures and possibly atmospheric refraction, and we set\n\
	dphi = d = 0 for the 0th exposure.\n\
	Sky coords are related to detector (chip) coords by\n\
		r_e = x + sum a_m f_m(x)\n\
	where m labels the modes, and where each mode coefficient a_m\n\
	is a 2-vector, f_1 = y, and the other modes are polynomials\n\
	of order 2 through 'lmax', and describe distortions of the\n\
	field. For ord = 4 say, the modes are:\n\
		y, x^2, xy, y^2, x^3, x^2 y, x y^2, y^3,\n\
		x^4, x^3 y, x^2 y^2, x y^3, y^4.\n\
	Model is linearised - and so only valid for small\n\
	dhi_e (d_e can be large though) - and we\n\
	solve for model coefficients by minimising squared residuals.\n\
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
#include <math.h>
#include "fitstack_modefunc.h"

#define SCALE 1.0

main(int argc, char *argv[])
{
	int	arg = 3, i, j, m, n, *indx, ncoefft, mcoefft, nobjects, e, ep;
	double	**A, *B, det, *C[2], x[2], xp[2], *d[2], *phi[2][2], *a[2];
	char	line[1024];
	int	nexp, nM, M, l, lmax, dbase[2], phibase[2][2], abase[2], *ll, *mm;
	FILE	*lcpipe;
	double	mag, magp, **Am, *Bm, *Cm, *dm;
	double	inputbuff[8];

	/* parse args */
	if (argc < 3 || argv[1][0] == '-')
		error_exit(usage);
	if (1 != sscanf(argv[1], "%d", &nexp))
		error_exit(usage);	
	if (1 != sscanf(argv[2], "%d", &lmax))
		error_exit(usage);
/*
	while (arg < argc) {
		if (argv[arg][0] != '-') {
			error_exit(usage);
		}
		switch (argv[arg++][1]) {
			case 'a':
				parfilename = argv[arg++];
				applytrans = 1;
				break;
			default:
				error_exit(usage);
				break;
		}
	}
*/

	/* figure the number of coefficients and the `base' indices where they live */
	dbase[0] = 0;
	dbase[1] = nexp - 1;
	phibase[0][0] = 2 * (nexp - 1);
	phibase[0][1] = 3 * (nexp - 1);
	phibase[1][0] = 4 * (nexp - 1);
	phibase[1][1] = 5 * (nexp - 1);
	nM = 0;
	for (l = 2; l <= lmax; l++) {
		nM += (l + 1);
	}
	ncoefft = 6 * (nexp - 1) + 2 * nM;
	abase[0] = 6 * (nexp - 1);
	abase[1] = abase[0] + nM;

	/* set up arrays of l, m values */
	ll = (int *) calloc(nM, sizeof(int));
	mm = (int *) calloc(nM, sizeof(int));
	M = 0;
	ll[M] = 1;
	mm[M++] = 1;
	for (l = 2; l <= lmax; l++) {
		for (m = 0; m <= l; m++) {
			ll[M] = l;
			mm[M++] = m;			
		}
	}
	
	/* allocate the big arrays */
        B = (double *) calloc(ncoefft, sizeof(double));
        C[0] = (double *) calloc(ncoefft, sizeof(double));
        C[1] = (double *) calloc(ncoefft, sizeof(double));
        A = (double **) calloc(ncoefft, sizeof(double *));
        for (m = 0; m < ncoefft; m++) {
                A[m] = (double *) calloc(ncoefft, sizeof(double));
        }
        indx = (int *) calloc(ncoefft, sizeof(int));

	/* similarly for the magnitude fitting stuff */
	mcoefft = nexp - 1;
	Bm = (double *) calloc(mcoefft, sizeof(double));
	Cm = (double *) calloc(mcoefft, sizeof(double));
	Am = (double **) calloc(mcoefft, sizeof(double));
	for (m = 0; m < mcoefft; m++) {
		Am[m] =  (double *) calloc(mcoefft, sizeof(double));
	}		
	
	/* allocate space for the model parameters */
	for (i = 0; i < 2; i++) {
		d[i] = (double *) calloc(nexp, sizeof(double));
		a[i] = (double *) calloc(nM, sizeof(double));
		for (j = 0; j < 2; j++) {
			phi[i][j] = (double *) calloc(nexp, sizeof(double));
		}
	}
	dm = (double *) calloc(nexp, sizeof(double));
	
	if (!(lcpipe = popen("lc -b -o x exp mag", "r"))) {
		fprintf(stderr, "mosaicfit: unable to open lc-pipe for input\n");
		exit(-1);
	}
	nobjects = 0;
	while (fread(inputbuff, sizeof(double), 8, lcpipe)) {
		x[0] 	= inputbuff[0];
		x[1] 	= inputbuff[1];
		xp[0] 	= inputbuff[2];
		xp[1] 	= inputbuff[3];
		e 	= (int) (inputbuff[4]);
		ep 	= (int) (inputbuff[5]);
		mag 	= inputbuff[6];
		magp 	= inputbuff[7];
		for (i = 0; i < 2; i++) {
			x[i] /= SCALE;
			xp[i] /= SCALE;
		}
		if (e < 0 || e > (nexp - 1) || ep < 0 || ep > (nexp - 1)) {
			fprintf(stderr, "mosaicfit: exposure number out of allowed range\n");
			exit(-1);
		}
		nobjects++;
		/* calculate C-vectors */
		/* first the exposure dependent terms */
		for (m = 1; m < nexp; m++) {
			for (n = 0; n < 2; n++) {
				for (i = 0; i < 2; i++) {
					C[n][dbase[i] + m - 1] = 0.0;
					for (j = 0; j < 2; j++) {
						C[n][phibase[i][j] + m - 1] = 0.0;
					}
				}
			}
			if (e == m) {
				C[0][dbase[0] + m - 1] += 1.0;
				C[1][dbase[1] + m - 1] += 1.0;
				C[0][phibase[0][0] + m - 1] += x[0];
				C[0][phibase[0][1] + m - 1] += x[1];
				C[1][phibase[1][0] + m - 1] += x[0];
				C[1][phibase[1][1] + m - 1] += x[1];
			}
			if (ep == m) {
				C[0][dbase[0] + m - 1] -= 1.0;
				C[1][dbase[1] + m - 1] -= 1.0;
				C[0][phibase[0][0] + m - 1] -= xp[0];
				C[0][phibase[0][1] + m - 1] -= xp[1];
				C[1][phibase[1][0] + m - 1] -= xp[0];
				C[1][phibase[1][1] + m - 1] -= xp[1];
			}
		}
		/* now the distortion terms */
		for (M = 0; M < nM; M++) {
			C[0][abase[0] + M] = C[1][abase[1] + M] = f(ll[M], mm[M], x) - f(ll[M], mm[M], xp);
		}
		/* accumulate A matrix, B-vector */
		for (i = 0; i < 2; i++) {
			for (m = 0; m < ncoefft; m++) {
				B[m] -= (x[i] - xp[i]) * C[i][m];
				for (n = 0; n < ncoefft; n++) {
					A[m][n] += C[i][m] * C[i][n];
				}
			}
		}
		/* now we do the magnitude terms */
		/* first the exposure terms...*/
		for (m = 1; m < nexp; m++) {
			Cm[m - 1] = 0.0;
			if (e == m) {
				Cm[m - 1] += 1.0;
			}
			if (ep == m) {
				Cm[m - 1] -= 1.0;
			}
		}
		/* accumulate Am matrix, Bm-vector */
		for (m = 0; m < mcoefft; m++) {
			Bm[m] += (mag - magp) * Cm[m];
			for (n = 0; n < mcoefft; n++) {
				Am[m][n] += Cm[m] * Cm[n];
			}
		}
	}


	if (nobjects <= ncoefft) {
		fprintf(stderr, "mosaicfit: too few objects\n");
		exit(-1);
	}

	/* solve the linear equations */
	myludcmp(A, ncoefft, indx, &det);
	mylubksb(A, ncoefft, indx, B);

	/* extract the model coefficients */
	for (e = 1; e < nexp; e++) {
		for (i = 0; i < 2; i++) {
			d[i][e] = B[dbase[i] + e - 1];
			for (j = 0; j < 2; j++) {
				phi[i][j][e] = B[phibase[i][j] + e - 1];
			}
		}
	}
	for (M = 0; M < nM; M++) {
		for (i = 0; i < 2; i++) {
			a[i][M] = B[abase[i] + M];
		}
	}
	
	/* solve equations for the magnitude offsets */
	myludcmp(Am, mcoefft, indx, &det);
	mylubksb(Am, mcoefft, indx, Bm);

	for (e = 1; e < nexp; e++) {
		dm[e] = Bm[e - 1];
	}

	/* output the result */
	fprintf(stdout, "# fitstack output\n%s\t# number of exposures\n%s\t# max order\n", argv[1], argv[2]);
	fprintf(stdout, "#        d[0]          d[1]     phi[0][0]     phi[0][1]     phi[1][0]     phi[1][1]          dmag\n");
	for (e = 0; e < nexp; e++) {
		for (i = 0; i < 2; i++) {
			fprintf(stdout, "%13.8lg ", d[i][e]);
		}
		for (i = 0; i < 2; i++) {
			for (j = 0; j < 2; j++) {
				fprintf(stdout, "%13.8lg ", phi[i][j][e]);
			}
		}
		fprintf(stdout, "%13.8lg \n", dm[e]);
	}
	fprintf(stdout, "#    l      m          a[0]          a[1]\n");
	for (M = 0; M < nM; M++) {
		fprintf(stdout, "%6d %6d %13.8lg %13.8lg\n", ll[M], mm[M], a[0][M], a[1][M]);
	}

	exit(0);
}



