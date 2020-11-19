#define	usage "\n\n\n\
NAME\n\
	removecrs - replace cosmic rays with MAGIC\n\
\n\
SYNOPSIS\n\
	removecrs	[option...]\n\
		-c	nuc	# threshold for centre pixel(default = 20)\n\
		-s  	nus	# threshold for surrounding pixels (default 5)\n\
		-n	nn	# required number of neighbours (default = 4)\n\n\
\n\
DESCRIPTION\n\
	'removecrs' replaces cosmic rays with MAGIC.\n\n\
	A pixel is replaced provided f > mode + nuc * sigma and\n\
	at least nn immediate neighbours have f < mode + nus * sigma\n\
	Reads from stdin and writes to stdout\n\
\n\
AUTHOR\n\
	Nick Kaiser:  kaiser@cita.utoronto.ca\n\
\n\n\n"		

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "imlib/fits.h"
#include "utils/error.h"
#include "utils/stats_stuff.h"
#include "utils/arrays.h"

#define MAGIC FLOAT_MAGIC

main(int argc, char *argv[])	
{
	float 		**f, **fout;
	int		arg = 1, N1, N2, i, j, in, jn, margin = 0;
	float		nuc, nus;
	fitsheader	*fits;
	fstatsrec	srec;
	int		nn, nngood;
	
	/* defaults */
	nuc = 20;
	nus = 5;
	nn = 4;
	
	while (arg < argc) {
		if (*argv[arg] != '-')
			error_exit(usage);
		switch (*(argv[arg++]+1)) {
			case 'c':
				if (1 != sscanf(argv[arg++], "%f", &nuc))
					error_exit(usage);
				break;
			case 's':
				if (1 != sscanf(argv[arg++], "%f", &nus))
					error_exit(usage);
				break;
			case 'n':
				if (1 != sscanf(argv[arg++], "%d", &nn))
					error_exit(usage);
				break;
			default:
				error_exit(usage);
				break;
		}
	}
	
	read2Dfloatimage(&f, &N1, &N2, &fits, stdin);		
	fdo_stats(f, N1, N2, margin, &srec);

	allocFloatArray(&fout, N1, N2);

	for (i = 0; i < N2; i++) {
		for (j = 0; j < N1; j++) {
			if (f[i][j] == MAGIC || f[i][j] < srec.fmode + nuc * srec.sigma) {
				fout[i][j] = f[i][j];
				continue;
			}
			nngood = 0;
			for (in = i - 1; in <= i + 1; in++) {
				if (in < 0 || in >= N2)
					continue;
				for (jn = j - 1; jn <= j + 1; jn++) {
					if (jn < 0 || jn >= N1)
						continue;
					if (jn == j && in == i)
						continue;
					if (f[in][jn] == MAGIC)
						continue;
					if (f[in][jn] < srec.fmode + nus * srec.sigma) 
						nngood++;
				}
			}
			if (nngood >= nn)
				fout[i][j] = MAGIC;
			else
				fout[i][j] = f[i][j];
		}
	}

	add_comment(argc, argv, fits);
	write2Dfloatimage(fout, fits);
	exit(0);
}







