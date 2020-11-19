#define usage "\n\n\n\
NAME\n\
	makeotf --- make fast guiding otf\n\
SYNOPSIS\n\
	makeotf [-N N] [-R r_outer] [-r r0] [-D D] [-o opfile] [-z zmax] [-i imname]\n\
\n\
DESCRIPTION\n\
	'makeotf' computes the OTF gk(z) for perfect fast guiding.\n\
\n\
	Options are\n\
		-N N		# image size in pixels (256)\n\
		-R r_outer	# outer scale in m (infinite)\n\
		-r r0		# Fried length in m (0.4)\n\
		-D D		# telescope diameter (1.6)\n\
		-e e		# obscuration = D_2 / D (0.0)\n\
		-o opfile	# output file\n\
		-z zmax		# upper limit for integerized z (N/2)\n\
		-d dz		# step in integer z (1)\n\
		-i imname	# output image 'imname' and exit\n\
		-g xg yg	# location of the guide star\n\
\n\
	Output is a table (in lc or human readable format) containing:\n\
		z		# argument of transfer function (in m)\n\
		S(z)		# atmospheric phase structure function\n\
		gknatl		# atmospheric OTF = exp(-S/2)\n\
		gkdiff		# diffraction limited OTF\n\
		gktilt		# fast guiding OTF\n\
		gkfried		# Fried approximation\n\
	except with -i option in which case it will output FITS image\n\
	'imname' which can be one of\n\
		T, S, Sg, TT, TTxS, WxS, TxT\n\
\n\
	Currently, -g option only works with pure Kolmogorov (infinite outer\n\
	scale) spectrum.\n\
\n\
AUTHOR\n\
	Nick Kaiser:  kaiser@hawaii.edu\n\
\n\n\n"


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "imlib/fits.h"
#include "utils/error.h"
#include "utils/args.h"
#include "fftlib/myfft.h"

void	gradx(float **f, int N1, int N2, float **fx);
void	mult(float **f, int N1, int N2, float fac);
void	zero(float **f, int N1, int N2);


int		main(int argc, char *argv[])	
{
	char		*flag, *opfilename, opstring[256], argstring[256], vonkarmanstring[256], *imname, tmpstring[256];
	char		lcstrbase[256], lcstr0[256], lcstr1[256];
	int		N, N2, NN, finiteouterscale, x, xp, y, xx, yy, z, dz, zmax, M1, M2, doimage, offsetguidestar;
	double		e, area, Tval, TTval, r_outer, r0, D, rap, pixsize, WWxSint, S0;
	double		gkdiff, gknatl, gktilt, gkmove, gkfried, xg, yg;
	float		**r, **T, **S, **Sg, **TT, **TTxS, **WxS, **tmp0, **TxT;
	fitsheader	*fits;
	fft_type	TTk, Sk, Tk;
	FILE		*opf = NULL, *spipe;

	/* defaults */
	finiteouterscale = 0;
	r_outer = 0.0;
	N 	= 256;
	r0	= 0.4;
	D	= 1.6;
	opfilename = NULL;
	zmax	= 0;
	doimage = 0;
	offsetguidestar = 0;
	xg = yg = 0.0;
	dz = 1;
	e = 0.0;

	/* parse args */
	argsinit(argc, argv, usage);
	while (flag = getflag()) {
		switch (flag[0]) {
			case 'N':
				N = getargi();
				break;
			case 'R':
				r_outer = getargd();
				finiteouterscale = 1;
				break;
			case 'r':
				r0 = getargd();
				break;
			case 'D':
				D = getargd();
				break;
			case 'e':
				e = getargd();
				break;
			case 'o':
				opfilename = getargs();
				break;
			case 'z':
				zmax = getargi();
				break;
			case 'i':
				doimage = 1;
				imname = getargs();
				break;
			case 'g':
				offsetguidestar = 1;
				xg = getargd();
				yg = getargd();
				break;
			case 'd':
				dz = getargi();
				break;
			case 'u':
			default:
				error_exit(usage);
		}
	}

	if (offsetguidestar && finiteouterscale) {
		error_exit("makeotf: illegal combination of -R and -g options\n");
	}	

	/* aperture size in pixels for critical sampling */
	rap = 0.25 * N;
	/* compute pixel size */
	pixsize = D / (2 * rap);
	/* area in pixels */
	area = M_PI * rap * rap * (1 - e * e);
	TTval = 1.0 / area;
	Tval = sqrt(TTval);
	/* half and squared box size */
	N2 = N / 2;
	NN = N * N;
	/* convert Fried length to pixels */
	r0 /= pixsize;
	/* convert guide star offset to pixels */
	xg /= pixsize;
	yg /= pixsize;
	/* set zmax */
	if (!zmax) {
		zmax = N2;
	}
	if (zmax > N2) {
		error_exit("makeotf: zmax must be less than or equal to N/2");
	}

	/* generate the aperture function T, TT = T * T and structure function S */
	allocFloatArray(&T, N, N);
	allocFloatArray(&r, N, N);
	allocFloatArray(&TT, N, N);
	allocFloatArray(&TxT, N, N);
	for (y = 0; y < N; y++) {
		yy = (y - N2) * (y - N2);
		for (x = 0; x < N; x++) {
			xx = (x - N2) * (x - N2);
			r[y][x] = sqrt(xx + yy);
			if (r[y][x] < rap && r[y][x] > e * rap) {
				T[y][x]  = Tval;
				TT[y][x] = TTval;
			}
		}
	}

	/* make structure function S */
	if (finiteouterscale) {
		sprintf(vonkarmanstring, "makevonkarmanS %lf %lf %lf %d | circimfromcat -F S -n %d", 
			r_outer, r0 * pixsize, pixsize, (int) ceil(N2 * sqrt(2.0)), N);
		spipe = popen(vonkarmanstring, "r");
		if (!spipe) {
			error_exit("makeotf: failed to open process to generate von Karman S(z)\n");
		}
		read2Dfloatimage(&S, &M1, &M2, &fits, spipe);
		if (M1 != N || M2 != N) {
			error_exit("makeotf: problem generating von Karman S\n");
		}
		Sg = S;
	} else {
		allocFloatArray(&S, N, N);
		if (offsetguidestar) {
			allocFloatArray(&Sg, N, N);
		} else {
			Sg = S;
		}
		for (y = 0; y < N; y++) {
			for (x = 0; x < N; x++) {
				S[y][x] = 6.88 * pow(r[y][x] / r0, 5.0 / 3.0);
				if (offsetguidestar) {
					Sg[y][x] = 6.88 * pow(((x - N2 - xg) * (x - N2 - xg) +
						(y - N2 - yg) * (y - N2 - yg))	/ (r0 * r0), 5.0 / 6.0);
				}
			}
		}
	}

	/* generate TxT */
	alloc_fft(&Tk, N, N);
	forward_fft(T, N, N, Tk);
	ccf(Tk, Tk, N, N, TxT, N2, N2);

	/* generate WxS for unshifted S() */
	allocFloatArray(&TTxS, N, N);	
	allocFloatArray(&WxS, N, N);
	alloc_fft(&TTk, N, N);
	alloc_fft(&Sk, N, N);
	forward_fft(TT, N, N, TTk);
	forward_fft(S, N, N, Sk);
	ccf(TTk, Sk, N, N, TTxS, N2, N2);
	mult(TTxS, N, N, (float) NN);
	gradx(TTxS, N, N, WxS);

	/* compute WWxS integral */
	WWxSint = 0.0;
	allocFloatArray(&tmp0, N, N);
	gradx(WxS, N, N, tmp0);
	for (y = 0; y < N; y++) {
		for (x = 0; x < N; x++) {
			WWxSint += tmp0[y][x] * TT[y][x];
		}
	}

	/* now compute WxS with S = Sg */
	forward_fft(Sg, N, N, Sk);
	forward_fft(TT, N, N, TTk);
	ccf(TTk, Sk, N, N, TTxS, N2, N2);
	mult(TTxS, N, N, (float) NN);
	gradx(TTxS, N, N, WxS);
	
	/* for debugging */
	if (doimage) {
		fits = new2Dfitsheader(N, N, FLOAT_PIXTYPE);
		if (!strcmp(imname, "T")) {
			write2Dfloatimage(T, fits);
		}
		if (!strcmp(imname, "S")) {
			write2Dfloatimage(S, fits);
		}
		if (!strcmp(imname, "Sg")) {
			write2Dfloatimage(Sg, fits);
		}
		if (!strcmp(imname, "TT")) {
			write2Dfloatimage(TT, fits);
		}
		if (!strcmp(imname, "TxT")) {
			write2Dfloatimage(TxT, fits);
		}
		if (!strcmp(imname, "TTxS")) {
			write2Dfloatimage(TTxS, fits);
		}
		if (!strcmp(imname, "WxS")) {
			write2Dfloatimage(WxS, fits);
		}
		exit(0);
	}

	argsToString(argc, argv, argstring);
	sprintf(tmpstring, "-H 'r0 = %14.8lg' -H 'D = %14.8lg' -H 'r_outer = %14.8lg' -H 'dz = %14.8lg' -H 'sigma_tilt = %14.8lg'", 
		r0 * pixsize, D, r_outer, dz * pixsize, 0.5 * WWxSint / (pixsize * pixsize)); 
	sprintf(lcstrbase, "lc -C -x -a 'history: %s' -n z -n S -n gknatl -n gkdiff -n gktilt -n gkfried %s", 
		argstring, tmpstring);
	sprintf(lcstr0, "%s < /dev/null", lcstrbase);
	system(lcstr0);
	if (opfilename) {
		sprintf(lcstr1, "%s > %s", lcstrbase, opfilename);
		opf = popen(lcstr1, "w");
	}

	for (z = 0; z < zmax; z += dz) {
		S0 = S[N2][N2 + z];
		gktilt = 0.0;
		gknatl = exp(-0.5 * S0);
		gkdiff = NN * TxT[N2][N2 + z];
		gkmove = exp(-0.25 * z * z * WWxSint);
		S0 += 0.5 * z * z * WWxSint;
		for (y = 0; y < N; y++) {
			for (x = 0; x < N - z; x++) {
				xp = x + z;
				if (T[y][x] > 0.0 && T[y][xp] > 0.0) {
					gktilt += TTval * exp(-0.5 * (S0 + z * (WxS[y][x] - WxS[y][xp])));
				}
			}
		}
		gknatl = (gknatl > 1.e-30 ? gknatl : 0.0);
		gkmove = (gkmove > 1.e-30 ? gkmove : 1.e-30);
		gkfried = gknatl * gkdiff / gkmove;
		sprintf(opstring, "  %14.8lg %14.8lg %14.8lg %14.8lg %14.8lg %14.8lg\n", z * pixsize, S[N2][N2 + z], gknatl, gkdiff, gktilt, gkfried);
		fprintf(stdout, opstring);
		if (opf) {
			fprintf(opf, opstring);
		}
	}	

	exit(0);
}


void	gradx(float **f, int N1, int N2, float **fx)
{
	int x, y;

	for (y = 0; y < N2; y++) {
		for (x = 0; x < N1 - 1; x++) {
			fx[y][x] = f[y][x+1] - f[y][x];
		}
	}
}

void	mult(float **f, int N1, int N2, float fac)
{
	int x, y;

	for (y = 0; y < N2; y++) {
		for (x = 0; x < N1; x++) {
			f[y][x] *= fac;
		}
	}
}


void	zero(float **f, int N1, int N2)
{
	int x, y;

	for (y = 0; y < N2; y++) {
		for (x = 0; x < N1; x++) {
			f[y][x] = 0.0;
		}
	}
}
