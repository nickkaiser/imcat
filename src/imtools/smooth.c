#define	usage "\n\
NAME\n\
	smooth --- spatially filter a fits image\n\
\n\
SYNOPSIS\n\
	smooth [option...]\n\
\n\
DESCRIPTION\n\
	'smooth' reads a fits image from standard input and writes a\n\
	smoothed version to standard output.\n\
\n\
	Various types of smoothing are provided:\n\
		-k	m rf		# smooth with m x m gaussian kernel\n\
		-b  	m		# m x m box filter\n\
		-t 			# tukey-style running median\n\
		-f 	a s1 s2		# fft filter (1+k^2 s1^2)^-a/2 exp(-0.5 k^2 s2^2)\n\
		-g	a b phi		# gaussian: major/minor = a/b, pos angle phi [deg]\n\
		-K	r		# fft kolmogorov turb: exp(-0.5 (k r)^(5/3))\n\
		-p	alpha		# power law: transfer function = k^(alpha)\n\
		-e	r gamma		# generalized exponential: exp(-(k r)^gamma)\n\
\n\
AUTHOR\n\
	Nick Kaiser:  kaiser@hawaii.edu\n\
\n\n\n"		

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "../imlib/fits.h"
#include "../utils/error.h"
#include "../imlib/filters.h"
#include "../utils/arrays.h"
#include "../fftlib/myfft.h"
#ifdef DUAL_PROC
#include "../fftlib/shmalloc.h"
#endif

#define	SCHECHTERFILTER		0
#define	GAUSSFILTER		1
#define KOLMOGOROVFILTER	2
#define POWERLAWFILTER		3
#define EXPONENTIALFILTER	4

#define MAGIC FLOAT_MAGIC

main(int argc, char *argv[])	
{
	float 	**f, **fs, magicsubstitute = 0.0, **temp;
	int	arg = 1, N1, N2, N3 = 1, iplane, mk, mb, i, j;
	int	dofft = 0, dokernel = 0, dotukey = 0, dobox = 0, ffttype = 0; 
	float	**ff, a, sigma1, sigma2, A, B, phi, rfk, rkolmogorov, powerlawalpha;
	float	exponentialsigma, exponentialgamma;
	fitsheader	*fits;
#ifdef DUAL_PROC
	int		fr_shmid;
#endif
	
	while (arg < argc) {
		if (*argv[arg] != '-')
			error_exit(usage);
		switch (*(argv[arg++]+1)) {
			case 'k':
				if (1 != sscanf(argv[arg++], "%d", &mk))
					error_exit(usage);
				if (1 != sscanf(argv[arg++], "%f", &rfk))
					error_exit(usage);
				dokernel = 1;
				break;
			case 'b':
				if (1 != sscanf(argv[arg++], "%d", &mb))
					error_exit(usage);
				dobox = 1;
				break;
			case 't':
				dotukey = 1;
				break;
			case 'f':
				if (1 != sscanf(argv[arg++], "%f", &a))
					error_exit(usage);
				if (1 != sscanf(argv[arg++], "%f", &sigma1))
					error_exit(usage);
				if (1 != sscanf(argv[arg++], "%f", &sigma2))
					error_exit(usage);
				dofft = 1;
				ffttype = SCHECHTERFILTER;
				break;
			case 'g':
				if (1 != sscanf(argv[arg++], "%f", &A))
					error_exit(usage);
				if (1 != sscanf(argv[arg++], "%f", &B))
					error_exit(usage);
				if (1 != sscanf(argv[arg++], "%f", &phi))
					error_exit(usage);
				dofft = 1;
				ffttype = GAUSSFILTER;
				break;
			case 'K':
				if (1 != sscanf(argv[arg++], "%f", &rkolmogorov))
					error_exit(usage);
				dofft = 1;
				ffttype = KOLMOGOROVFILTER;
				break;
			case 'p':
				if (1 != sscanf(argv[arg++], "%f", &powerlawalpha))
					error_exit(usage);
				dofft = 1;
				ffttype = POWERLAWFILTER;
				break;
			case 'e':
				if (1 != sscanf(argv[arg++], "%f", &exponentialsigma))
					error_exit(usage);
				if (1 != sscanf(argv[arg++], "%f", &exponentialgamma))
					error_exit(usage);
				dofft = 1;
				ffttype = EXPONENTIALFILTER;
				break;
			default:
				error_exit(usage);
				break;
		}
	}
	
#ifdef DUAL_PROC
	fr_shmid = read2Dfloatimage_shm(&f, &N1, &N2, &fits, stdin);
	set_frshmid(fr_shmid);
#else
	fits = readfitsheader(stdin);
	/* read2Dfloatimage(&f, &N1, &N2, &fits, stdin); */
	N1 = fits->n[0];
	N2 = fits->n[1];
	if (fits->ndim == 3) {
		N3 = fits->n[2];
	}
	add_comment(argc, argv, fits);
	writefitsheader(fits);
	allocFloatArray(&fs, N1, N2);
	allocFloatArray(&f, N1, N2);
#endif

	for (iplane = 0; iplane < N3; iplane++) {
	readfitsplane((void **) f, fits);
		
	/* if (dokernel || dobox) 
		allocFloatArray(&fs, N1, N2); */
		
	if (dotukey)
		tukey(f, N1, N2);
	if (dobox) {	
		block_filter(f, fs, N1, N2, mb);
		temp = f; f = fs; fs = temp;
	}
	if (dokernel) {
		gaussian_kernel_filter(f, fs, N1, N2, mk, rfk);
		temp = f; f = fs; fs = temp;
	}
	if (dofft) {
#ifdef DUAL_PROC
		fs = f;
#else
		allocFloatArray(&fs, N1, N2);
#endif
		switch (ffttype) {
			case SCHECHTERFILTER:
				schecterfilter(f, N1, N2, fs, sigma1, sigma2, a, magicsubstitute);
				break;
			case GAUSSFILTER:
				gaussfilter(f, N1, N2, fs, A, B, phi, magicsubstitute);
				break;
			case KOLMOGOROVFILTER:
				kolmogorovfilter(f, N1, N2, fs, rkolmogorov, magicsubstitute);
				break;
			case POWERLAWFILTER:
				powerlawfilter(f, N1, N2, fs, powerlawalpha, magicsubstitute);
				break;
			case EXPONENTIALFILTER:
				exponentialfilter(f, N1, N2, fs, exponentialsigma, exponentialgamma, magicsubstitute);
				break;
			default:
				error_exit("smooth: bad fft filter type\n");
		}
#ifndef DUAL_PROC
		for (i = 0; i < N2; i++)
			for (j = 0; j < N1; j++)
				if (f[i][j] == MAGIC)
					fs[i][j] = MAGIC;
#endif
		f = fs;
	}

	writefitsplane((void **) f, fits);
	}
writefitstail(fits);

#ifdef DUAL_PROC
	shmfree(fr_shmid);
#endif
	exit(0);
}








