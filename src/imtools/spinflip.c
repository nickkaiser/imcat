#define usage "\n\n\n\
NAME\n\
	spinflip - apply 90 degree rotations or parity flips to FITS image\n\
\n\
SYNOPSIS\n\
	spinflip mii mij mji mjj\n\
\n\
DESCRIPTION\n\
	\"spinflip\" reflects images or rotates them through multiples\n\
	of 90 degrees according to\n\
		di' = Mii * di + Mij * dj\n\
		dj' = Mji * di + Mjj * dj\n\
	with matrix elements M = -1, 0, 1 and with the constraint\n\
		Mii^2 + Mij^2 = Mji^2 + Mjj^2  = 1\n\
\n\
AUTHOR\n\
	Nick Kaiser:  kaiser@cita.utoronto.ca\n\
\n\n\n"


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "../imlib/fits.h"
#include "../utils/error.h"

int		main(int argc, char *argv[])	
{
	int		i, j, ip, jp, N, N1, N2, mii, mij, mji, mjj;
	fitsheader	*fits;
	float		**fin, **fout;
	
	if (argc != 5)
			error_exit(usage);
	if (1 != sscanf(argv[1], "%d", &mii) ||
		1 != sscanf(argv[2], "%d", &mij) ||
		1 != sscanf(argv[3], "%d", &mji) ||
		1 != sscanf(argv[4], "%d", &mjj)) error_exit(usage);
		
	if ((mii * mii + mij * mij) != 1 || (mji * mji + mjj * mjj) != 1)
		error_exit(usage);
	
	read2Dfloatimage(&fin, &N1, &N2, &fits, stdin);

	if (N1 != N2)
		error_exit("spinflip: I can only handle square images\n");
	N = N1;
	allocFloatArray(&fout, N, N);
	
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			ip = mii * i + mij * j;
			if (mii < 0 || mij < 0)
				ip = N - 1 + ip;
			jp = mji * i + mjj * j;
			if (mji < 0 || mjj < 0)
				jp = N - 1 + jp;
			fout[ip][jp] = fin[i][j];
		}
	}
	
	add_comment(argc, argv, fits);
	write2Dfloatimage(fout, fits);
	exit(0);
}



