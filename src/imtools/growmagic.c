#define usage "\n\n\n\
NAME\n\
	growmagic --- expand MAGIC regions of an image\n\
SYNOPSIS\n\
	growmagic [-n n]\n\
\n\
DESCRIPTION\n\
	\"growmagic\" reads a source image from stdin. It makes\n\
	a copy and then sets to MAGIC any pixel for which any\n\
	of its eight neighbours have the MAGIC value.\n\
\n\
	With '-n n' option it performs this operation n times.\n\
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
#include "../utils/args.h"

#define max(a,b)    ((a)<(b)?(b):(a)) 
#define min(a,b)    ((a)<(b)?(a):(b))

#define MAGIC FLOAT_MAGIC

int		main(int argc, char *argv[])	
{
	int		ix, iy, N1, N2, jx, jy, nop;
	fitsheader	*fits;
	float		**fin, **fout, **ftemp;
	char		*flag;
	
	/* defaults */
	nop = 1;

	argsinit(argc, argv, usage);
	while (flag = getflag()) {
		switch (flag[0]) {
			case 'n':
				nop = getargi();
				break;
			default:
				error_exit(usage);
		}
	}


	read2Dfloatimage(&fin, &N1, &N2, &fits, stdin);
	allocFloatArray(&fout, N1, N2);

	while (nop) {
		/* make a copy */	
		for (iy = 0; iy < N2; iy++) {
			for (ix = 0; ix < N1; ix++) {
				fout[iy][ix] = fin[iy][ix];
			}
		}
		/* nuke pixels around magic ones */
		for (iy = 0; iy < N2; iy++) {
			for (ix = 0; ix < N1; ix++) {
				if (fin[iy][ix] == MAGIC) {
					for (jy = max(0, iy-1); jy <= min(N2-1, iy+1); jy++) {
						for (jx = max(0, ix-1); jx <= min(N1-1, ix+1); jx++) {
							fout[jy][jx] = MAGIC;
						}
					}
				}
			}
		}
		/* swap images */
		ftemp = fout;
		fout = fin;
		fin = ftemp;
		nop--;
	}

	add_comment(argc, argv, fits);
	write2Dfloatimage(fin, fits);
	exit(0);
}



