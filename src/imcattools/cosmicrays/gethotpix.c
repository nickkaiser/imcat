#define	usage "\n\n\n\
NAME\n\
	gethotpix - get candidate cosmic ray pixels\n\
SYNOPSIS\n\
	gethotpix fmin rmax\n\
\n\
DESCRIPTION\n\
	'gethotpix' reads a 2-dimensional FITS image from stdin\n\
	and sends to stdout a catalogue containing the locations\n\
	of pixels with f > fmin and an immediate N,S,E or W\n\
	neighbour with r = f_neighbour / f < rmax.\n\
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

#define MAGIC FLOAT_MAGIC

#define BUFFSIZE 4

void	writehotpixel(int x, int y, float f, float fn);

main(int argc, char *argv[])	
{
	float 		**f, **fout;
	int		arg = 1, N1, N2, x, y;
	float		fmin, rmax;
	fitsheader	*fits;
	FILE		*op;
	
	/* parse args */
	if (argc != 3) {
		error_exit(usage);
	}	
	if (1 != sscanf(argv[1], "%f", &fmin))
		error_exit(usage);
	if (1 != sscanf(argv[2], "%f", &rmax))
		error_exit(usage);
	
	/* read the image */
	read2Dfloatimage(&f, &N1, &N2, &fits, stdin);

	/* generate the header */
	system("lc -C -b -N '1 2 x' -n f -n fn -a gethotpix -x < /dev/null");		

	for (y = 0; y < N2; y++) {
		for (x = 0; x < N1; x++) {
			if (f[y][x] > fmin) {
				if (x < (N1 - 1)) {
					if (f[y][x + 1] < rmax * f[y][x]) {
						writehotpixel(x, y, f[y][x], f[y][x + 1]);
					}
				}
				if (x > 0) {
					if (f[y][x - 1] < rmax * f[y][x]) {
						writehotpixel(x, y, f[y][x], f[y][x - 1]);
					}
				}
				if (y < (N2 - 1)) {
					if (f[y + 1][x] < rmax * f[y][x]) {
						writehotpixel(x, y, f[y][x], f[y + 1][x]);
					}
				}
				if (y > 0) {
					if (f[y - 1][x] < rmax * f[y][x]) {
						writehotpixel(x, y, f[y][x], f[y - 1][x]);
					}
				}
			}
		}
	}
	exit(0);
}


void	writehotpixel(int x, int y, float f, float fn)
{
	static double	buff[BUFFSIZE];


	buff[0] = (double) x;
	buff[1] = (double) y;
	buff[2] = (double) f;
	buff[3] = (double) fn;
	fwrite(buff, sizeof(double), BUFFSIZE, stdout);
}




