#define	usage "\n\n\n\
SYNOPSIS\n\
		composite	[option...] a.cat b.cat ....\n\
			-r	min max	# range of half-light radii (default 2-3)\n\
			-l  	min max	# range of log luminosity (default 4-8)\n\
			-m	M	# array size (32)\n\
			-s		# scale by r * r / l\n\
			-i 		# generate an image: text -> stderr\n\
\n\
DESCRIPTION\n\
		\"composite()\" forms a MODAL galaxy surface brightness profile\n\
		from a collection of cat files created by findpeaks().\n\
		Outputs modal SB + ellpticity alpha parameter to stdout\n\
\n\n\n"		
	

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "fits.h"
#include "error.h"
#include "stats_stuff.h"
#include "object_stuff.h"
#include "magic.h"
#include "catio.h"
#include "arrays.h"
#include "composite.h"

#define MAX_VALS	20000
#define N_ALPHA		6

float	lbar, rbar;
int		scale;

main(int argc, char *argv[])	
{
	short 		**f, **fimage;
	float		***fval;
	int		**index;
	int		arg = 1; 
	object 		obj;
	cathead		thecat;
	FILE		*fitsf, *catf, *textf, *fitsoutf;
	int			comc;
	char		*comv[MAX_COMMENTS], argstring[512];
	float		rmin, rmax, lmin, lmax;
	int		ngals = 0, r, N;
	float		sumrn4df[N_ALPHA], sumrn3f[N_ALPHA], alpha[N_ALPHA];
	float		n[N_ALPHA] = {0.0, -0.5, -1.0, -1.5, -2.0, -2.5};
	int		i, j, M, image, ialpha;
	float		**fmode;
	float		*avfmode, *nn, *df;
	float		lquart, uquart, median, sigma;
	
	/* defaults */
	rmin = 2;
	rmax = 3;
	lmin = exp(4.0);
	lmax = exp(8.0);
	M = 32;
	image = 0;
	scale = 0;
	
	while (arg < argc) {
		if (*argv[arg] == '-') {
			switch (*(argv[arg++]+1)) {
				case 'r':
					if (EOF == sscanf(argv[arg++], "%f", &rmin))
						error_exit(usage);
					if (EOF == sscanf(argv[arg++], "%f", &rmax))
						error_exit(usage);
					break;
				case 'l':
					if (EOF == sscanf(argv[arg++], "%f", &lmin))
						error_exit(usage);
					if (EOF == sscanf(argv[arg++], "%f", &lmax))
						error_exit(usage);
					lmin = exp(lmin);
					lmax = exp(lmax);
					break;
				case 'i':
					image = 1;
					break;
				case 's':
					scale = 1;
					break;
				case 'm':
					if (EOF == sscanf(argv[arg++], "%d", &M))
						error_exit(usage);
					break;
				default:
					error_exit(usage);
					break;
			}
		} else {
			break;			/* have read all the args */
		}
	}
	rbar = 0.5 * (rmin + rmax);
	lbar = sqrt(lmin * lmax);

	/* allocate space for big arrays */
	fval = (float ***) calloc(M, sizeof(float **)) + M / 2;
	index = (int **) calloc(M, sizeof(int *)) + M / 2;
	fmode = (float **) calloc(M, sizeof(float *)) + M / 2;
	for (i = - M / 2; i < M / 2; i++) {
		fval[i] = (float **) calloc(M, sizeof(float *)) + M / 2;
		index[i] = (int *) calloc(M, sizeof(int)) + M / 2;
		fmode[i] = (float *) calloc(M, sizeof(float)) + M / 2;
		for (j = - M / 2; j < M / 2; j++) {
			fval[i][j] = (float *) calloc(MAX_VALS, sizeof(float));
		}
	}

	/* loop over files and fill arrays */
	for ( ; arg < argc; arg++) {
		catf = fopen(argv[arg], "r");
		if (!catf)
			error_exit(usage);
		set_cat_ipf(catf);
		read_cat_head(&thecat);
		fitsf = fopen(thecat.image, "r");
		if (!fitsf)
			error_exit("composite: unable to open fits file for input");
		set_fits_ipf(fitsf);
		read_fits(&f, &N, &N, &comc, comv);
		while (read_object(&obj)) {
			if (obj.rh > rmin && obj.rh < rmax && obj.l > lmin && obj.l < lmax) {
				ngals++;
				accumulate(obj, fval, index, M, f, N);
			}
		}
		fclose(fitsf);
		fclose(catf);
		freeShortArray(f, N, N);
	}

	/* now find the modes */
	for (i = - M /  2; i < M / 2; i++) {
		for (j = - M / 2; j < M / 2; j++) {
			liststats(fval[i][j], index[i][j], &(fmode[i][j]), &median, 
						&lquart, &uquart, &sigma);
		}
	}


	if (image)
		textf = stderr;
	else
		textf = stdout;
	
	avfmode = (float *) calloc(M / 2, sizeof(float));
	nn = (float *) calloc(M / 2, sizeof(float));
	df = (float *) calloc(M / 2, sizeof(float));
	for (i = - M /  2; i < M / 2; i++) {
		for (j = - M / 2; j < M / 2; j++) {
			if (i != 0 || j != 0)
				r = floor(0.5 + sqrt(i * i + j * j));
			else
				r = 0;
			if (r >= M / 2)
				continue;
			nn[r] += 1.0;
			avfmode[r] += fmode[i][j];
		}
	}

	for (r = 0; r < M / 2; r++) {
		if (nn[r] > 0.0)
			avfmode[r] /= nn[r];
	}
	for (r = 0; r < M / 2 - 1; r++) {
		if (r == 0)
			df[r] = avfmode[0] - avfmode[1];
		else
			df[r] = 0.5 * (avfmode[r-1] - avfmode[r+1]);
	}
	df[M / 2 - 1] = 0;
	for (ialpha = 0; ialpha < N_ALPHA; ialpha++)
		sumrn4df[ialpha] = sumrn3f[ialpha] = 0;
	fprintf(textf, "# composite profile from %d objects\n", ngals);
	fprintf(textf, "# %8.3f < r < %8.3f\n", rmin, rmax);
	fprintf(textf, "# %8.3f < ln(l) < %8.3f\n# alpha's for n = ", log(lmin), log(lmax));
	for (ialpha = 0; ialpha < N_ALPHA; ialpha++) {
		fprintf(textf, " %6.3f", n[ialpha]);
	}
	fprintf(textf, "\n#  r      <f>    df/dr alpha[1] alpha[2] ...\n");
	for (r = 0; r < M / 2; r++) {
		fprintf(textf, "%4d %8.3lf %8.3lf", 
			r, avfmode[r], df[r]);
		for (ialpha = 0; ialpha < N_ALPHA; ialpha++) {
			sumrn4df[ialpha] += pow((double) r, 4 + n[ialpha]) * df[r];
			sumrn3f[ialpha] += pow((double) r, 3 + n[ialpha]) * avfmode[r];
			if (sumrn3f[ialpha] > 0.0)
				alpha[ialpha] = sumrn4df[ialpha] / 
					((4 + n[ialpha]) * sumrn3f[ialpha]);
			else
				alpha[ialpha] = 0.0;
			fprintf(textf, " %8.3lf", alpha[ialpha]);
		}
		fprintf(textf, "\n");
	}
	if (image) {
		allocShortArray(&fimage, M, M);
		add_comment(argc, argv, &comc, comv);
		for (i = - M / 2; i < M / 2; i++) {
			for (j = - M / 2; j < M / 2; j++) {
				fimage[i + M / 2][j + M / 2] = 100 * fmode[i][j];
			}
		}
		write_fits(fimage, M, M, comc, comv);
	}
	exit(0);
}


void		accumulate(object obj, float ***fval, int **index, int M, short **f, int N)
{
	int 	i, j, di, dj;
	
	for (di = - M / 2; di < M / 2; di++) {
		for (dj = - M / 2; dj <  M / 2; dj++) {
			i = obj.i + di;
			j = obj.j + dj;
			if (i < 0 || i >= N || j < 0 || j >= N)
				continue;
			if (f[i][j] == MAGIC)
				continue;
			if (scale) 
				fval[di][dj][ index[di][dj] ] = f[i][j] * 
					lbar * obj.rh * obj.rh / (obj.l * rbar * rbar);
			else
				fval[di][dj][ index[di][dj] ] = f[i][j];
			index[di][dj]++;
			if (index[di][dj] >= MAX_VALS)
				error_exit("accumulate: bomb\n");
		}
	}
}





