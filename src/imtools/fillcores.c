#define	usage "\n\n\n\
NAME\n\
	fillcores --- patch up holes in stellar cores\n\
\n\
SYNOPSIS\n\
	fillcores threshold replaceval [corevalue]\n\
\n\
DESCRIPTION\n\
	fillcores is a crude kludge to patch up holes in stellar\n\
	cores as sometimes get produced by clumsy image arithmetic.\n\
	By default, any pixel with MAGIC value which has a neighbour with\n\
	value exceeding 'threshold' will be replaced by 'replaceval'.\n\
	If the optional third argument is given then any pixel\n\
	with value 'corevalue' will be replaced by 'replaceval'.\n\
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
#include "../utils/stats_stuff.h"
#include "../utils/arrays.h"

main(int argc, char *argv[])	
{
	float 	**f, **fout, coreval, threshold, replaceval;
	int	arg = 1, N1, N2, i, j, in, jn, replaced;
	fitsheader	*fits;
	
	if (argc != 3 && argc != 4) {
		error_exit(usage);
	}
	if ((1 != sscanf(argv[1], "%f", &threshold)) ||
		(1 != sscanf(argv[2], "%f", &replaceval))) {
			error_exit(usage);
	}
	if (argc == 4) {
		if (1 != sscanf(argv[3], "%f", &coreval)) {
			error_exit(usage);
		}
	} else {
		coreval = FLOAT_MAGIC;
	}
	
	
	read2Dfloatimage(&f, &N1, &N2, &fits, stdin);		
	allocFloatArray(&fout, N1, N2);

	for (i = 0; i < N2; i++) {
		for (j = 0; j < N1; j++) {
			fout[i][j] = f[i][j];
			if (f[i][j] != coreval) {
				continue;
			}
			replaced = 0;
			for (in = i - 1; in <= i + 1; in++) {
				if (in < 0 || in >= N2)
					continue;
				for (jn = j - 1; jn <= j + 1; jn++) {
					if (jn < 0 || jn >= N1)
						continue;
					if (jn == j && in == i)
						continue;
					if (f[in][jn] > threshold) { 
						fout[i][j] = replaceval;
						replaced = 1;
						break;
					}
				}
				if (replaced) {
					break;
				}
			}
		}
	}

	add_comment(argc, argv, fits);
	write2Dfloatimage(fout, fits);
	exit(0);
}







