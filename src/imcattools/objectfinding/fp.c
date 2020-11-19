#define	usage "\n\n\
NAME\n\
        fp --- hierarchical object finder\n\
\n\
SYNOPSIS\n\
        fp [-u] [-r r1 r2] [-n nplanes]\n\
\n\
DESCRIPTION\n\
        fp .....\n\
\n\
AUTHOR\n\
        Nick Kaiser --- kaiser@hawaii.edu\n\
\n\n"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "../../imlib/fits.h"
#include "../../utils/error.h"
#include "../../utils/args.h"
#include "../../imlib/filters.h"
#include "../../utils/arrays.h"

#define MAGIC 		FLOAT_MAGIC
#define PI		M_PI
#define TINY		1.e-20

main(int argc, char *argv[])	
{
	float 		**f0, **f1, **f2, magicsubstitute;
	fitsheader	*fitsin, *fitsout;
	double		r1, r2;
	int		nplanes, N1, N2;
	char		*flag;
	
	/* defaults */
	r1 = 2.0;
	r2 = 4.0;
	nplanes = 16;
	magicsubstitute = 0.0;

	argsinit(argc, argv, usage);
	while (flag = getflag()) {
		switch (flag[0]) {
			case 'r':
				r1 = getargd();
				r2 = getargd();
				break;
			case 'n':
				nplanes = getargi();
				break;
			default:
				error_exit(usage);
		}
	}

	read2Dfloatimage(&f0, &N1, &N2, &fitsin, stdin);
	
	allocFloatArray(&f1, N1, N2);
	allocFloatArray(&f2, N1, N2);

	mexicanfilter(f0, N1, N2, f1, r1, 2 * r1, magicsubstitute);
	mexicanfilter(f0, N1, N2, f2, r2, 2 * r2, magicsubstitute);

	fitsout = copyfitsheader(fitsin);
	add_comment(argc, argv, fitsout);
	fitsout->ndim = 3;
	fitsout->n[2] = 2;
	writefitsheader(fitsout);
	writefitsplane((void **) f1, fitsout);	
	writefitsplane((void **) f2, fitsout);	
	writefitstail(fitsout);
	exit(0);
}




#define ZMAX	10

/*
 * mexican hat filter function 
 */
float	thefilterfunction(float ki, float kj)
{
	float	kk, z1, z2;
		
	kk = ki * ki + kj * kj;
	z1 = kk * sigma1;
	if (z1 > ZMAX)
		return (0.0);
	z2 = kk * sigma2;
	if (z2 > ZMAX) {
		return (exp(-z1));
	} else {
		return (exp(-z1) - exp(-z2));
	}
}

#undef ZMAX






