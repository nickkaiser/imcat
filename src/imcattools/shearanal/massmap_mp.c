#define	usage "\n\n\n\
NAME\n\
	massmap_mp - maximum probability reconstruction\n\
\n\
SYNOPSIS\n\
		massmap_mp [option...] <input.cat >output.image\n\
			-g ng			# output grid size (64)\n\
			-n N			# input image size (2048)\n\
			-k kmax			# maximum k-value (5)\n\
			-a alpha		# power parameter (0.01)\n\
			-e ename		# name for 2-vector ellipticity (e)\n\
			-x xname		# name for 2-vector spatial coordinate (x)\n\
\n\
DESCRIPTION\n\
		\"massmap_mp\" reads x[2], e[2] from a catalogue and calculates foreground\n\
		surface mass density field from background galaxy ellipticities a la KS.\n\
		Uses maximum probability technique diagonal stripe controlled by alpha.\n\
		Use alpha << 1 for weak regularization.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@cita.utoronto.ca\n\
\n\n\n"


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "../../utils/arrays.h"
#include "../../utils/error.h"
#include "../../imlib/fits.h"
#include "../../utils/lu.h"

double		modefunc(int m, double i, double j);
double		*k[2];
#define		MAX_GALS	20000

typedef struct object {
	double	i, j, e[1];
} object;

#define	PI	M_PI

main(int argc, char *argv[])	
{
	int		arg = 1;
	object		obj;
	FILE		*catf;
	char		errorstring[1024];
	int		m, n, ng, N, ig, jg, i, j, kmax, ki, kj, nmodes, pol, *indx;
	double		scale;			/* output grid points = pixels x scale */
	double		x[2], e[2];
	int		ngal, igal;	
	double		dk, kk, *chi[2], **A, *B, *a, d, aa, alpha, **modef, chichi;
	float		**sigma;
	FILE		*tempf;
	char		defename[32] = "e", defxname[32] = "x", *ename, *xname, lcstring[128];
	char		line[256], tempfile[128], systemstring[128];
	FILE		*lcpipe;
	fitsheader	*fits;

	/* defaults */
	ng = 64;
	N = 2048;
	kmax = 5;
	alpha = 0.01;
	xname = defxname;
	ename = defename;
	
	/* get the optional arguments */
	while (arg < argc) {
		if (argv[arg][0] == '-') {			/* an option argument */
			switch (argv[arg++][1]) {
				case 'g':
					if (1 != sscanf(argv[arg++], "%d", &ng))
						error_exit(usage);
					break;					
				case 'n':
					if (1 != sscanf(argv[arg++], "%d", &N))
						error_exit(usage);
					break;					
				case 'k':
					if (1 != sscanf(argv[arg++], "%d", &kmax))
						error_exit(usage);
					break;					
				case 'a':
					if (1 != sscanf(argv[arg++], "%lf", &alpha))
						error_exit(usage);
					break;					
				case 'e':
					ename = argv[arg++];
					break;
				case 'x':
					xname = argv[arg++];
					break;
				default:
					error_exit(usage);
					break;
			}
		} else {
			error_exit(usage);
		}
	}

	/* set option dependent variables */
	scale = (double) N / (double) ng;
	dk = 2 * PI / (double) N;

	/* set up the output image */
	fits = new2Dfitsheader(ng, ng, FLOAT_PIXTYPE);
	allocFloatArray(&sigma, ng, ng);

	/* count the number of k-modes (times two) */
	nmodes = 0;
	for (ki = -kmax; ki <= kmax; ki++) {
		for (kj = 0; kj <= kmax; kj++) {
			if ((kj == 0 && ki <= 0) || (ki * ki + kj * kj > kmax * kmax))
				continue;
			nmodes += 2;
		}
	}
	fprintf(stderr, "# nmodes = %d\n", nmodes);

	/* allocate space for k, A, B matrices, chi array etc */
	k[0] = (double *) calloc(nmodes, sizeof(double));
	k[1] = (double *) calloc(nmodes, sizeof(double));
	chi[0] = (double *) calloc(nmodes, sizeof(double));
	chi[1] = (double *) calloc(nmodes, sizeof(double));
	a = (double *) calloc(nmodes, sizeof(double));
	B = (double *) calloc(nmodes, sizeof(double));
	A = (double **) calloc(nmodes, sizeof(double *));
	for (m = 0; m < nmodes; m++) {
		A[m] = (double *) calloc(nmodes, sizeof(double));
	}
	indx = (int *) calloc(nmodes, sizeof(int));
	modef = (double **) calloc(nmodes, sizeof(double *));

	/* read catalogue, count objects and write them to tempfile */
	sprintf(tempfile, "%d.tmp", getpid());
	if (!(tempf = fopen(tempfile, "w")))
		error_exit("massmap_mp: unable to open temp file\n");
	sprintf(lcstring, "lc -o %s %s", xname, ename);
	if (!(lcpipe = popen(lcstring, "r"))) {
		error_exit("etprofile: failed to open lc-pipe for input\n");
	}
	ngal = 0;
	while (fgets(line, 255, lcpipe)) {
		fprintf(tempf, "%s", line);
		ngal++;
	}
	fclose(tempf);

	/* allocate modef arrray */
	for (m = 0; m < nmodes; m++)
		modef[m] = (double *) calloc(ngal, sizeof(double));		

	/* set up k, chi arrays */
	m = 0;
	for (ki = -kmax; ki <= kmax; ki++) {
		for (kj = 0; kj <= kmax; kj++) {
			kk = ki * ki + kj * kj;
			if ((kj == 0 && ki <= 0) || (kk > kmax * kmax))
				continue;
			chi[0][m] = chi[0][m + 1] = (ki * ki - kj * kj) / kk;			
			chi[1][m] = chi[1][m + 1] = 2 * ki * kj / kk;
			k[0][m] = k[0][m + 1] = ki * dk;			
			k[1][m] = k[1][m + 1] = kj * dk;
			m += 2;			
		}
	}

	/* read catalogue and accumulate A, B arrays */
	if (!(tempf = fopen(tempfile, "r")))
		error_exit("massmap_mp: unable to open temp file for input\n");
	ngal = 0;
	while (fgets(line, 255, tempf)) {
		if (4 != sscanf(line, "%lf %lf %lf %lf", &(x[0]), &(x[1]), &(e[0]), &(e[1])))
			error_exit("massmap: input format error\n");
		/* convert to old style i,j convention */
		obj.i = x[1];
		obj.j = x[0];
		obj.e[0] = -e[0];
		obj.e[1] = e[1];
		for (m = 0; m < nmodes; m++) {
			modef[m][ngal] = modefunc(m, obj.i, obj.j);
			B[m] += (obj.e[0] * chi[0][m]  + obj.e[1] * chi[1][m]) *  modef[m][ngal];
		}
		ngal++;
	}
	fclose(tempf);
	sprintf(systemstring, "rm %s", tempfile);
	system(systemstring);

	for (m = 0; m < nmodes; m++) {
		for (n = 0; n <= m; n++) {
			chichi = (chi[0][m] * chi[0][n] + chi[1][m] * chi[1][n]);
			for (igal = 0; igal < ngal; igal++) {
				A[m][n] += chichi * modef[m][igal] * modef[n][igal];
			}
		}
	}
	for (m = 0; m < nmodes; m++) {
		for (n = m+1; n < nmodes; n++) {
			A[m][n] = A[n][m];
		}
	}

	/* add regularising stuff */
	for (m = 0; m < nmodes; m++)
		A[m][m] += alpha * ngal;

	/* solve A * a = b */
	for (m = 0; m < nmodes; m++)
		a[m] = B[m];
	myludcmp(A, nmodes, indx, &d);
	mylubksb(A, nmodes, indx, a);

	/* make the reconstruction */
	for (ig = 0; ig < ng; ig++) {
		i = floor(0.5 + scale * ig);
		for (jg = 0; jg < ng; jg++) {
			j = floor(0.5 + scale * jg);
			for (m = 0; m < nmodes; m++) {
				sigma[ig][jg] += a[m] * modefunc(m, i, j);
			}
		}	
	}	

	/* output */
	add_comment(argc, argv, fits);
	write2Dfloatimage(sigma, fits);
	exit(0);
}



double	modefunc(int m, double i, double j)
{
	if (2 * (m / 2) == m)
		return(cos(k[0][m] * i + k[1][m] * j));
	else
		return(sin(k[0][m] * i + k[1][m] * j));
}

#undef PI
